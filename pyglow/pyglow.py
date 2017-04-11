"""
Global variable indicating if IRI 2016 has been initialized with
the contents of the ionosphere global index (ig_rz.dat) and Ap/F10.7
index (apf107.dat) files. IRI 2016 initialization is required only
once per session.
"""
__INIT_IRI16 = False

def _igrf_tracefield(dn, lat, lon, alt, target_ht, step):
    import numpy as np

    # Go North:
    lla_north = _igrf_tracefield_hemis(dn, lat, lon, alt,\
            target_ht, step)

    # Go South:
    lla_south = _igrf_tracefield_hemis(dn, lat, lon, alt,\
            target_ht, -step)

    # Stack them together:
    lla = np.vstack( [ np.flipud(lla_north)[:-1,:], lla_south ])

    return lla


def _igrf_tracefield_hemis(dn, lat, lon, alt, target_ht, step):
    import numpy as np
    import coord
    from numpy import array as arr

    lat = float(lat)
    lon = float(lon)
    alt = float(alt)
    target_ht = float(target_ht)
    step = float(step)

    target_ht = target_ht*1e3
    step = step*1e3

    lla = np.array([lat, lon, alt*1e3])

    lla_field = lla

    """ Step 1: trace the field along a given direction """
    TOLERANCE = 10 # [m]
    i = 0
    while (lla[2] > target_ht):
        # convert to ECEF:
        ecef = coord.lla2ecef(lla)

        # Grab field line information:
        p = Point(dn, lla[0], lla[1], lla[2]/1e3)
        p.run_igrf()

        # coordinates follow pyglow's convention:
        # (x -> east, y -> north, z -> up)

        N = p.By # North
        E = p.Bx # East
        D = -p.Bz # Down
        A = p.B  # Total

        # Step along the field line
        ecef_new =  ecef + coord.ven2ecef(lla,[(-D / A * step), (E / A * step),
                                               (N / A * step)])

        # Convert to lla coordinates:
        lla = coord.ecef2lla(ecef_new)

        # add the field line to our collection:
        lla_field = np.vstack( [lla_field, lla] )
        i = i + 1

    """ Step 2: Make the last point close to target_ht """
    while (abs(lla[2]-target_ht) > TOLERANCE):
        my_error = lla[2]-target_ht
        old_step = step
        # find out how much we need to step by:
        step = -np.sign(step)*abs(lla[2]-target_ht)

        ecef = coord.lla2ecef(lla)

        p = Point(dn, lla[0], lla[1], lla[2]/1e3)
        p.run_igrf()

        N = p.Bx # North
        E = p.By # East
        D = p.Bz # Down
        A = p.B  # Total

        # trace the field, but use the modified step:
        ecef_new = ecef + coord.ven2ecef(lla, np.array([(-D/A), (E/A), N/A])
                                         * step/(-D/A))
        # TODO : I changed this, is this correct?
        lla = coord.ecef2lla(ecef_new)

    # replace last entry with the point close to target_ht:
    lla_field[-1,:] = lla

    return lla_field

def Line(dn, lat, lon, alt, target_ht=90., step=15.):
    '''
    pts = Line(dn, lat, lon, alt, target_ht=90., step=15.)

    Return a list of instances of Point by
    tracing along the geomagnetic field line.

    target_ht : altitude to quit tracing at
                for N and S hemispheres [km]
    step : approximate step to take between points [km]
    '''
    llas = _igrf_tracefield(dn, lat, lon, alt, target_ht, step)
    pts = []
    for lla in llas:
        pts.append( Point(dn, lla[0], lla[1], lla[2]/1e3) )
    return pts

def update_kpap(years=None):
    '''
    Update the Kp and Ap indices used in pyglow.
    The files will be downloaded from noaa to your pyglow
    installation directory.

    update_kpap(years=None)

    Inputs:
    ------
    years : (optional) a list of years to download.
            If this input is not provided, the full
            range of years starting from 1932 to the
            current year will be downloaded.
    '''
    from datetime import date,timedelta
    import shutil
    import urllib2
    from contextlib import closing
    import pyglow

    # Load all data up until today
    if years is None: years=range(1932, date.today().year + 1)

    pyglow_dir =\
            '/'.join(pyglow.__file__.split("/")[:-1]) + "/kpap/"

    for year in years:
        src = 'ftp://ftp.ngdc.noaa.gov/'\
                + 'STP/GEOMAGNETIC_DATA/INDICES/KP_AP/%4i'\
                % (year,)
        des = pyglow_dir + "%4i" % (year,)
        print "\nDownloading"
        print src
        print "to"
        print des
        try:
            with closing(urllib2.urlopen(src)) as r:
                with open(des, 'wb') as f:
                    shutil.copyfileobj(r, f)
        except IOError as e:
            print 'Failed downloading data for year %i. File does not exist' \
                % year


def update_dst(years=None):
    '''
    Update the Dst index files used in pyglow.
    The files will be downloaded from WDC Kyoto
    to your pyglow installation directory.

    update_dst(years=None)

    Inputs:
    ------

    years : (optional) a list of years to download.
            If this input is not provided, the full
            range of years starting from 2005 to the
            current year will be downloaded. Pre-2005
            files are shipped with pyglow.
    '''
    from datetime import date, timedelta
    import urllib2
    from contextlib import closing
    import pyglow


    def download_dst(year, month, des):
        '''
        Helper function to earch for the appropriate location
        and download the DST index file from WDC Kyoto for the
        given month and year. Save it to the specified file "des".
        Return True if successful, False if not.
        '''
        # There are three possible sources of data. Search for
        # them in the following order:
        # 1) Final
        # 2) Provisional
        # 3) Realtime
        year_month = '%i%02i' % (year, month)
        wgdc_fn = 'dst%s%02i.for.request' % (str(year)[2:],month)
        src_final = 'http://wdc.kugi.kyoto-u.ac.jp/dst_final/%s/%s' % \
                          (year_month, wgdc_fn)
        src_provisional = 'http://wdc.kugi.kyoto-u.ac.jp/dst_provisional/%s/%s'\
                          % (year_month, wgdc_fn)
        src_realtime = 'http://wdc.kugi.kyoto-u.ac.jp/dst_realtime/%s/%s' \
                       % (year_month, wgdc_fn)

        success = False
        for src in [src_final, src_provisional, src_realtime]:
            try:
                with closing(urllib2.urlopen(src)) as r:
                    contents = r.read()
                    # If that succeeded, then the file exists
                    print "\nDownloading"
                    print src
                    print "to"
                    print des
                    with open(des,'w') as f:
                        f.write(contents)
                    success = True
                    break
            except urllib2.HTTPError:
                pass
        return success

    # Read files from 2005 until today. Pre-2005
    # files are shipped with pyglow.
    if years is None:
        years = range(2005,date.today().year + 1)
    pyglow_dir =\
            '/'.join(pyglow.__file__.split("/")[:-1]) + "/dst/"

    for year in years:
        for month in range(1,13):
            des = '%s%i%02i' % (pyglow_dir,year,month)
            download_dst(year, month, des)

def update_ae(years = None):
    '''
    Update the AE index files used in pyglow.
    The files will be downloaded from WDC Kyoto
    to your pyglow installation directory.

    update_ae(years=None)

    Inputs:
    ------

    years : (optional) a list of years to download.
            If this input is not provided, the full
            range of years starting from 2005 to the
            current year will be downloaded. Pre-2005
            files are shipped with pyglow.
    '''
    from datetime import date, timedelta
    import urllib2
    from contextlib import closing
    import pyglow
    import glob

    def download_ae(year, month, des):
        '''
        Helper function to earch for the appropriate location
        and download the AE index file from WDC Kyoto for the
        given month and year. Save it to the specified file "des".
        Return True if successful, False if not.
        '''
        # There are three possible sources of data. Search for
        # them in the following order:
        # 1) Final
        # 2) Provisional
        # 3) Realtime
        year_month = '%i%02i' % (year, month)
        wgdc_fn = 'ae%s%02i.for.request' % (str(year)[2:],month)
        #src_final = 'http://wdc.kugi.kyoto-u.ac.jp/ae_final/%s/%s' \
        #            % (year_month, wgdc_fn)
        src_provisional = 'http://wdc.kugi.kyoto-u.ac.jp/ae_provisional/%s/%s'\
                          % (year_month, wgdc_fn)
        src_realtime = 'http://wdc.kugi.kyoto-u.ac.jp/ae_realtime/%s/%s' % \
                       (year_month, wgdc_fn)

        success = False
        for src in [src_provisional, src_realtime]:
            try:
                with closing(urllib2.urlopen(src)) as r:
                    contents = r.readlines()
                    # If that succeeded, then the file exists
                    print "\nDownloading"
                    print src
                    print "to"
                    print des
                    with open(des,'w') as f:
                        # this shrinks the filesize to hourly
                        for c in contents:
                            f.write("%s%s%s\n"%(c[12:18],c[19:21],c[394:400]))
                    success = True
                    break
            except urllib2.HTTPError:
                pass
        return success

    # Read files from 2000 until today. Pre-2000
    # files are shipped with pyglow.
    if years is None:
        years = range(2000,date.today().year + 1)
    pyglow_dir =\
            '/'.join(pyglow.__file__.split("/")[:-1]) + "/ae/"

    for year in years:
        for month in range(1,13):
            des = '%s%i%02i' % (pyglow_dir,year,month)
            download_ae(year, month, des)


def update_indices(years = None):
    '''
    Update all geophysical indices (e.g., kp, dst).

    update_indices(years=None)

    Inputs:
    ------

    years : (optional) a list of years to download.
            If this input is not provided, default
            values will be used.
    '''
    update_kpap(years=years)
    update_dst(years=years)
    update_ae(years=years)


class Profile(object):
    '''Object containing geomagnetic indices and model values for an altitude
    profile

    Parameters
    ------------
    dn : (datetime.datetime)
        Datetime object denoting universal time for profile
    lat : (float)
        Geographic latitude of profile in degrees
    lon : (float)
        Geographic longitude of profile in degrees
    alt_min : (float)
        Lower limit of altitude profile in km (default=90.0)
    alt_max : (float)
        Upper limit of altitude profile in km (default=2000.0)
    alt_step : (float)
        Altitude steps for profile in km (default=10.0)
    user_ind : (bool)
        Allow user specified indices (default=False)

    Returns
    ---------
    self : Profile class object with attributes
    '''

    def __init__(self, dn, lat, lon, alt_min=90.0, alt_max=2000.0,
                 alt_step=10.0, user_ind=False):
        '''Object that contains an altitude profile
        '''
        from get_apmsis import get_apmsis
        import numpy as np

        # record input:
        self.dn = dn
        self.lat = lat
        self.lon = lon
        self.alt_step = alt_step
        self.alt = np.arange(alt_min, alt_max + alt_step, alt_step)

        # Warn if date is too early
        if self.dn.year < 1932:
            raise ValueError('Date cannot be before 1932!')

        self.doy = self.dn.timetuple().tm_yday
        self.utc_sec = self.dn.hour*3600. + self.dn.minute*60.
        self.utc_hour = self.dn.hour
        self.slt_hour = np.mod(self.utc_sec/3600. + self.lon/15., 24)
        self.iyd = np.mod(self.dn.year,100)*1000 + self.doy

        # initialize variables:
        # ---------------------

        # for kp, ap function
        self.kp       = np.nan
        self.ap       = np.nan
        self.f107     = np.nan
        self.f107a    = np.nan
        self.kp_daily = np.nan
        self.ap_daily = np.nan
        self.apmsis   = [np.nan,]*7
        self.dst      = np.nan
        self.ae       = np.nan

        # for iri:
        self.iri_version = np.nan
        self.ne = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        ions = ['O+', 'H+', 'HE+', 'O2+', 'NO+']
        self.ni={ion:np.ones(shape=self.alt.shape, dtype=float) * np.nan
                 for ion in ions}

        self.Ti = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.Te = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.Tn_iri = np.ones(shape=self.alt.shape, dtype=float) * np.nan

        self.NmF2 = np.nan
        self.hmF2 = np.nan
        self.NmF1 = np.nan
        self.hmF1 = np.nan
        self.NmE = np.nan
        self.hmE = np.nan
        self.NmD = np.nan
        self.hmD = np.nan

        # for msis:
        self.msis_version = np.nan
        self.Tn_msis = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        neuts = ['HE','O','N2','O2','AR','H','N','O_anomalous']
        self.nn = {neu:np.ones(shape=self.alt.shape, dtype=float) * np.nan
                   for neu in neuts}
        self.rho = np.ones(shape=self.alt.shape, dtype=float) * np.nan

        # for hwm 93/07:
        self.u = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.v = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.hwm_version = np.nan

        # for igrf:
        self.igrf_version = np.nan
        self.Bx = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.By = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.Bz = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.B = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.dip = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.dec = np.ones(shape=self.alt.shape, dtype=float) * np.nan

        # for run_airglow:
        self.ag_stat = False
        self.ag6300 = np.ones(shape=self.alt.shape, dtype=float) * np.nan
        self.ag7774 = np.ones(shape=self.alt.shape, dtype=float) * np.nan

        if not user_ind:
            # call the indice models:
            self.get_indices()
            self.apmsis = get_apmsis(self.dn)

    def __repr__(self):
        import numpy as np
        # Print time and location information
        out = "{:>36s}{:}\n".format("date and time = ", self.dn)
        loc_line = "lat, lon, alt min, max, step (km) = "
        out = "{:s}{:>36s}".format(out, loc_line)
        out = "{:s}{:.2f}, {:.2f}, ".format(out, self.lat, self.lon)
        out = "{:s}{:.2f}, {:.2f}, ".format(out, self.alt[0], self.alt[-1])
        out = "{:s}{:.2f}\n\n".format(out, self.alt_step)

        # Print geophysical indices:
        out = "{:s}Geophysical Indices:\n---------------------\n".format(out)
        out = "{:s}{:>36s}{:4.2f}\n".format(out, "kp = ", self.kp)
        out = "{:s}{:>36s}{:4.2f}\n".format(out, "ap = ", self.ap)
        out = "{:s}{:>36s}{:4.2f}\n".format(out, "f107 = ", self.f107)
        out = "{:s}{:>36s}{:4.2f}\n".format(out, "f107a = ", self.f107a)
        out = "{:s}{:>36s}{:3.2f}\n".format(out, "Dst = ", self.dst)
        out = "{:s}{:>36s}{:3.2f}\n".format(out, "ae = ", self.ae)

        # Specify which models have been loaded
        out = "{:s}\nModel Data Available For:\n-----------\n".format(out)
        out = "{:s}{:>36s}{:}\n".format(out, "IGRF = ", self.igrf_version
                                        if not np.isnan(self.igrf_version)
                                        else False)
        out = "{:s}{:>36s}{:}\n".format(out, "HWM = ", self.hwm_version
                                        if not np.isnan(self.hwm_version)
                                        else False)
        out = "{:s}{:>36s}{:}\n".format(out, "IRI = ", self.iri_version
                                        if not np.isnan(self.iri_version)
                                        else False)
        out = "{:s}{:>36s}{:}\n".format(out, "MSIS = ", self.msis_version
                                        if not np.isnan(self.msis_version)
                                        else False)
        out = "{:s}{:>36s}{:}".format(out, "Airglow = ", self.ag_stat)

        return out

    def get_indices(self):
        from get_kpap import get_kpap

        (self.kp, self.ap, self.f107, self.f107a, self.kp_daily, self.ap_daily,
         self.dst, self.ae) = get_kpap(self.dn)
        return self

    @staticmethod
    def init_iri16():
        """ If required (depending on the global variable *__INIT_IRI16*),
        initialize IRI 2016. Return True if the model was initialized and False
        otherwise.
        """
        if not globals()['__INIT_IRI16']:
            from iri16py import read_ig_rz, readapf107
            read_ig_rz()
            readapf107()
            globals()['__INIT_IRI16'] = True
            return True
        return False

    def run_iri(self, NmF2=None, hmF2=None, version=2016, compute_Ne=True,
                compute_Te_Ti=True, compute_Ni=True, keep_workdir=True,
                debug=False):
        """ Run IRI model at point time/location and update the object state
        accordingly.

        Parameters
        -------------
        NmF2 : (NoneType or float)
           Density of the F2 layer in per cubic cm, see documentation in
           IRI_SUB (default=None)
        hmF2 : (NoneType or float)
           Height of the F2 layer in km, see documantaiton in IRI_SUB
           (default=None)
        version : (int)
            IRI Model version.  Currently 2016 and 2012 are valid.
            (default=2016)
        compute_Ne :  (bool)
            Calculate electron density (default=True)
        compute_Te_Ti : (bool)
            Calculate ion and electron temperatures (default=True)
        compute_Ni : (bool)
            Calculate ion population densities (default=True)
        keep_workdir : (bool)
            IRI must be run in the IRI directory.  This keywork specifies
            whether to return to the current working directory (True) or stay
            in theIRI directory after running (False). (default=True)
        debug : (bool)
            Output debugging information (default=False)

        Returns
        ---------
        Updates profile class object parameters:
        ne = electron density (per cubic cm)
        ni = ion densities (per cubic cm) for ['O+', 'H+', 'HE+', 'O2+', 'NO+']
        Ti = ion temperature (K)
        Te = electron temperature (K)
        Tn_iri = neutral temperature (K)
        NmF2 = ion density at the F2 peak (per cubic cm)
        hmF2 = height of the F2 peak (km)
        NmF1 = ion density at the F1 peak, if present (per cubic cm)
        hmF1 = height of the F1 peak, if present (km)
        NmE = ion density at the E peak (per cubic cm)
        hmE = height of the E peak (km)
        NmD = ion density at the D peak (per cubic cm)
        hmD = height of the D peak (km)
        """
        from iri12py import iri_sub as iri12
        from iri16py import iri_sub as iri16
        import os
        import sys
        import pyglow
        import numpy as np

        if version==2016:
            iri_data_stub = '/iri16_data/'
            iri = iri16
            init_iri = Profile.init_iri16
        elif version==2012:
            iri_data_stub = '/iri12_data/'
            iri = iri12
            init_iri = lambda: False
        else:
            raise ValueError('Invalid version of \'%i\' for IRI.' % (version) +\
                    '\nEither 2016 (default) or 2012 is valid.')

        if debug:
            print "version = {:d}".format(version)

        jf = np.ones((50,)) # JF switches
        # Standard IRI model flags
        #             | FORTRAN Index
        #             |
        #             V
        jf[3]  = 0 #  4 B0,B1 other model-31
        jf[4]  = 0 #  5  foF2 - URSI
        jf[5]  = 0 #  6  Ni - RBV-10 & TTS-03
        jf[20] = 0 # 21 ion drift not computed
        jf[22] = 0 # 23 Te_topside (TBT-2011)
        jf[27] = 0 # 28 spreadF prob not computed
        jf[28] = 0 # 29 (29,30) => NeQuick
        jf[29] = 0 # 30
        jf[32] = 0 # 33 Auroral boundary model off
                   #    (Brian found a case that stalled IRI when on)
        jf[34] = 0 # 35 no foE storm update

        # Not standard, but outputs same as values as standard so not an issue
        jf[21] = 0 # 22 ion densities in m^-3 (not %)
        jf[33] = 0 # 34 turn messages off

        if not compute_Ne:
            jf[0] = 0

        if not compute_Te_Ti:
            jf[1] = 0

        if not compute_Ni:
            jf[2] = 0

        oarr = np.zeros((100,))

        if NmF2 is not None:
            # use specified F2 peak density
            jf[7] = 0
            oarr[0] = NmF2 * 100.**3  # IRI expects [m^{-3}]

        if hmF2 is not None:
            # use specified F2 peak height
            jf[8] = 0
            oarr[1] = hmF2

        # IRI must be run in it's own directory
        my_pwd = os.getcwd()
        iri_data_path = '/'.join(pyglow.__file__.split("/")[:-1]) \
                        + iri_data_stub
        if debug:
            print "changing directory to \n", iri_data_path

        if my_pwd != iri_data_path:
            os.chdir(iri_data_path)

        # Run IRI
        init_iri()
        outf = iri(jf, 0, self.lat, self.lon, self.dn.year, -self.doy,
                   (self.utc_sec/3600.+25.), self.alt[0], self.alt[-1],
                   self.alt_step, oarr)
        
        if keep_workdir:
            if debug: print "returning to working directory \n", my_pwd
            os.chdir(my_pwd)

        # Save desired output
        if compute_Te_Ti:
            self.Te = outf[3][0:self.alt.shape[0]] # electron temperature (K)
            self.Ti = outf[2][0:self.alt.shape[0]] # ion temperature (K)

        self.Tn_iri = outf[1][0:self.alt.shape[0]] # neutral temperature (K)

        # Electron and ion densities, converting from /cubic m to /cubic cm
        self.ne = outf[0][0:self.alt.shape[0]] * 1.0e-6
        self.ni['O+'] = outf[4][0:self.alt.shape[0]] * 1.0e-6
        self.ni['H+'] = outf[5][0:self.alt.shape[0]] * 1.0e-6
        self.ni['HE+'] = outf[6][0:self.alt.shape[0]] * 1.0e-6
        self.ni['O2+'] = outf[7][0:self.alt.shape[0]] * 1.0e-6
        self.ni['NO+'] = outf[8][0:self.alt.shape[0]] * 1.0e-6

        self.NmF2 = oarr[0] * 1.0e-6
        self.hmF2 = oarr[1]
        # F1 peak is not always present
        if oarr[2] != -1.0 and oarr[3] != -1.0:
            self.NmF1 = oarr[2] * 1.0e-6
            self.hmF1 = oarr[3]
        self.NmE = oarr[4] * 1.0e-6
        self.hmE = oarr[5]
        self.NmD = oarr[6] * 1.0e-6
        self.hmD = oarr[7]

        if not np.isnan(self.hmF2):
            self.iri_version = version

        return self

    def run_msis(self, version=2000):
        from msis00py import gtd7 as msis00
        import numpy as np

        if version==2000:
            msis = msis00
        else:
            raise ValueError('Invalid version of \'%i\' for MSIS.' % (version)+\
                    '\n2000 (default) is valid.')

        msis_lon = np.mod(self.lon,360)
        for ia, aa in enumerate(self.alt):
            [d,t] = msis(self.doy, self.utc_sec, aa, self.lat, msis_lon,
                         self.slt_hour, self.f107a, self.f107, self.apmsis, 48)

            # Save output
            self.Tn_msis[ia] = t[1] # neutral temperature from MSIS (K)
            self.nn['HE'][ia] = d[0] # [items/cm^3]
            self.nn['O'][ia] = d[1] # [items/cm^3]
            self.nn['N2'][ia] = d[2] # [items/cm^3]
            self.nn['O2'][ia] = d[3] # [items/cm^3]
            self.nn['AR'][ia] = d[4] # [items/cm^3]
            # [5] is below
            self.nn['H'][ia] = d[6] # [items/cm^3]
            self.nn['N'][ia] = d[7] # [items/cm^3]
            self.nn['O_anomalous'][ia] = d[8] # [items/cm^3]

            self.rho[ia] = d[5] # total mass density [grams/cm^3]

        if not np.isnan(self.rho.all()):
            self.msis_version = version
            
        return self

    def run_hwm(self, version=2014):
        if version==2014:
            self.run_hwm14()
        elif version==2007:
            self.run_hwm07()
        elif version==1993:
            self.run_hwm93()
        else:
            raise ValueError('Invalid version of \'%i\' for HWM.' % (version) +\
                    '\nEither 2014 (default), 2007, or 1993 is valid.')
        return self

    def run_hwm93(self):
        from hwm93py import gws5 as hwm93
        import numpy as np

        self.hwm_version = 1993
        hwm_lon = np.mod(self.lon,360)
        for ia,aa in enumerate(self.alt):
            w = hwm93(self.iyd, self.utc_sec, aa, self.lat, hwm_lon,
                      self.slt_hour, self.f107a, self.f107, self.ap_daily)
            self.v[ia] = w[0]
            self.u[ia] = w[1]
        return self

    def run_hwm07(self):
        from hwm07py import hwmqt as hwm07
        import pyglow
        import os
        import numpy as np

        self.hwm_version = 2007
        my_pwd = os.getcwd()
        hwm07_data_path = '/'.join(pyglow.__file__.split("/")[:-1]) + \
                          "/hwm07_data/"
        os.chdir(hwm07_data_path)
        
        hwm_lon = np.mod(self.lon,360)
        for ia,aa in enumerate(self.alt):
            w = hwm07(self.iyd, self.utc_sec, aa, self.lat, hwm_lon,
                      self.slt_hour, self.f107a, self.f107, [np.nan, self.ap])
            self.v[ia] = w[0]
            self.u[ia] = w[1]
            
        os.chdir(my_pwd)
        return self

    def run_hwm14(self):
        from hwm14py import hwm14
        import pyglow
        import os
        import numpy as np

        my_pwd = os.getcwd()
        hwm14_data_path = '/'.join(pyglow.__file__.split("/")[:-1]) + \
                          "/hwm14_data/"
        os.chdir(hwm14_data_path)

        hwm_lon = np.mod(self.lon,360)
        for ia,aa in enumerate(self.alt):
            self.v[ia], self.u[ia] = hwm14(self.iyd, self.utc_sec, aa, self.lat,
                                           hwm_lon, np.nan, np.nan, np.nan,
                                           [np.nan, self.ap])
            
        os.chdir("%s" % my_pwd)
        self.hwm_version = 2014
        return self

    def run_igrf(self, version=12):
        from igrf11py import igrf11syn as igrf11
        from igrf12py import igrf12syn as igrf12
        import warnings
        import numpy as np

        if version==12:
            igrf=igrf12
        elif version==11:
            igrf=igrf11
        else:
            raise ValueError('Invalid version of \'%i\' for IGRF.' % (version)+\
                             '\nVersion 12 (default) and 11 are valid.')

        igrf_lat = 90.0 - self.lat
        igrf_lon = np.mod(self.lon, 360)
        warnings.warn("Caution: IGRF coordinates have been recently changed" + \
                      " to\n Bx -> positive eastward\n By -> positive " + \
                      "northward\n Bz -> positive upward\n")
        for ia,aa in enumerate(self.alt):
            x, y, z, f = igrf(0, self.dn.year, 1, aa, igrf_lat, igrf_lon)

            # Note that the changes here match
            # coordinate convention with other models
            # (i.e., HWM), that is:
            # (x -> east, y -> north, z -> up)
            #
            # IGRF gives (x -> north, y -> east, z -> down)
            # B in [T] (note x/y switch here and negation of z)
            self.Bx[ia] = y * 1.0e-9
            self.By[ia] = x * 1.0e-9
            self.Bz[ia] = -z * 1.0e-9
            self.B[ia] =  f * 1e-9
            h = np.sqrt(x**2 + y**2)
            self.dip[ia] = np.degrees(np.arctan2(z, h))
            self.dec[ia] = np.degrees(np.arctan2(y, x))

        if not np.isnan(self.B.all()):
            self.igrf_version = version
            
        return self

    def run_airglow(self):
        ''' Computes airglow intensities

        After running, self.ag6300 has the 630.0-nm volume emission
        rate (ph/cm^3/s) and self.ag7774 has the 777.4-nm volume
        emission rate.

        History
        ------
        9/9/13 -- implemented into pyglow
                  based on Jonathan J. Makela's MATLAB
                  and subsequent python code
        9/13/16 -- added 7774 calculation
        '''
        import numpy as np
        
        # let's see if IRI and MSIS have been executed
        # if not, run the appropriate models:
        if np.isnan(self.ne).all():
            self.run_iri()
        if np.isnan(self.nn['O2']).all():
            self.run_msis()

        # Perform 630.0-nm calculation
        te = self.Te / 300.0
        ti = self.Ti / 300.0

        # These coefs are from Link and Cogger, JGR 93(A9), 988309892, 1988
        K1_6300 = 3.23e-12 * np.exp(3.72 / ti - 1.87 / ti**2)
        K2_6300 = 2.78e-13 * np.exp(2.07 / ti - 0.61 / ti**2)
        K3_6300 = 2.0e-11 * np.exp(111.8 / self.Tn_msis)
        K4_6300 = 2.9e-11 * np.exp(67.5 / self.Tn_msis)
        K5_6300 = 1.6e-12 * self.Te**0.91
        b6300 = 1.1

        # Corrected values from Link and Cogger, JGR, 94(A2), 1989
        a1D = 7.45e-3
        a6300 = 5.63e-3

        # Calculate O+ assuming mixture of ions (Link and Cogger, 1989)
        a1 = 1.95e-7 * te**-0.7
        a2 = 4.00e-7 * te**-0.9
        Oplus = self.ne / (1.0 + K2_6300 * self.nn['N2'] / a2 / self.ne
                           + K1_6300 * self.nn['O2'] / a1 / self.ne)

        AGNumerator = a6300 / a1D * b6300 * K1_6300 * Oplus * self.nn['O2']
        AGDenominator = 1.0 + (K3_6300 * self.nn['N2'] + K4_6300 * self.nn['O2']
                               + K5_6300 * self.ne) / a1D
        self.ag6300 = AGNumerator / AGDenominator

        # Perform the 777.4-nm calculation

        # These coefs are from a number of sources (see Makela's
        # dissertation Table 2.3 or Makela et al., "Ionospheric
        # topography maps using multiple-wavelength all-sky images",
        # JGR, 106, 29161--29174, 2001.)
        alpha1_7774 = 7.8e-13
        beta_7774 = 0.42
        K1_7774 = 1.3e-15
        K2_7774 = 1.5e-7
        K3_7774 = 1.4e-10

        # Makela shows that Oplus may be replaced by Ne below in his
        # dissertation (see p. 25) with only a 0.5% impact. Since we
        # have Oplus at hand, we use it in the calculation.
        V7774_rr = alpha1_7774 * Oplus * self.ne

        V7774_ii_num = beta_7774 * K1_7774 * K2_7774 * self.nn['O'] * Oplus \
                       * self.ne
        V7774_ii_den = K2_7774 * Oplus + K3_7774 * self.nn['O']

        self.ag7774 = V7774_rr + V7774_ii_num / V7774_ii_den

        if not np.isnan(self.ag7774.all()):
            self.ag_stat = True
            
        return self

# ---------

class Point(Profile):
    def __init__(self, dn, lat, lon, alt, user_ind=False):
        Profile.__init__(self, dn, lat, lon, alt_min=alt, alt_max=alt,
                         user_ind=user_ind)

        # record input that differs from Profile:
        self.alt_float = self.alt[0]

    def flatten(self):
        '''Change variable dimention to reflect single point in profile
        '''
        import numpy as np
        
        # for iri:
        if not np.isnan(self.iri_version):
            try:
                self.ne = self.ne[0]
                self.ni = {ion : self.ni[ion][0] for ion in self.ni.keys()}
                self.Ti = self.Ti[0]
                self.Te = self.Te[0]
                self.Tn_iri = self.Tn_iri[0]
            except: pass

        # for msis:
        if not np.isnan(self.msis_version):
            try:
                self.Tn_msis = self.Tn_msis[0]
                self.nn = {neu : self.nn[neu][0] for neu in self.nn.keys()}
                self.rho = self.rho[0]
            except: pass

        # for hwm 93/07:
        if not np.isnan(self.hwm_version):
            try:
                self.u = self.u[0]
                self.v = self.v[0]
            except: pass

        # for igrf:
        if not np.isnan(self.igrf_version):
            try:
                self.Bx = self.Bx[0]
                self.By = self.By[0]
                self.Bz = self.Bz[0]
                self.B = self.B[0]
                self.dip = self.dip[0]
                self.dec = self.dec[0]
            except: pass

        # for run_airglow:
        if self.ag_stat:
            try:
                self.ag6300 = self.ag6300[0]
                self.ag7774 = self.ag7774[0]
            except: pass

# ---------


if __name__=="__main__":
    from datetime import datetime, timedelta
    import numpy as np

    dn = datetime(2002, 3, 25, 12, 0, 0)
    lat = 42.5
    lon = 0.
    alt = 250.
    o = Point(dn, lat, lon, alt)


    o.run_iri()
    o.run_msis()
    o.run_hwm93()
    o.run_hwm07()
    o.run_igrf()



    #l = Line(dn, lat, lon, alt)


    '''
    dn_start = datetime(2000,1,1)
    dn_end   = datetime(2005,1,1)
    dns = [dn_start + timedelta(days=kk) for kk in range((dn_end-dn_start).days)]
    f107a = []
    for dn in dns:
        p = Point(dn, 0, 0, 300)
        f107a.append( p.f107a )

    from matplotlib.pyplot import *
    figure(1); clf()
    plot(dns, f107a)
    grid()
    draw(); show();
    kmin = np.argmin(f107a)
    kmax = np.argmax(f107a)
    print dns[kmin], f107a[kmin]
    print dns[kmax], f107a[kmax]
    '''

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import contextlib
from datetime import date, timedelta
import glob
import numpy as np
import os
import shutil
import sys
import warnings
import urllib.request, urllib.error, urllib.parse

from . import coord
from .get_kpap import get_kpap
from .get_apmsis import get_apmsis
from hwm93py import gws5 as hwm93
from hwm07py import hwmqt as hwm07
from hwm14py import hwm14
from igrf11py import igrf11syn as igrf11
from igrf12py import igrf12syn as igrf12
from iri12py import iri_sub as iri12
from iri16py import iri_sub as iri16
from iri16py import read_ig_rz, readapf107
from msis00py import gtd7 as msis00

"""
Global variable indicating if IRI 2016 has been initialized with
the contents of the ionosphere global index (ig_rz.dat) and Ap/F10.7
index (apf107.dat) files. IRI 2016 initialization is required only
once per session.
"""
__INIT_IRI16 = False
__version__ = '1.2'

DIR_FILE = os.path.dirname(__file__)


class Point(object):

    def __init__(
        self,
        dn,
        lat,
        lon,
        alt,
        user_ind=False,
    ):
        """
        An instance of Point is the fundamental data object
        for running each climatological model.

        Instantation of a Point initializes member variables,
        and also grabs the corresponding geophysical indices.

        :param dn: datetime.datetime object
        :param lat: Latitude [degrees]
        :param lon: Longitude [degrees]
        :param alt: Altitude [km]
        :param user_ind: (optional) Boolean switch to calculate
                         geophysical indices. If True, then it
                         is up to the user to assign geophysical
                         indices to the Point
        """

        nan = float('nan')

        # Record input:
        self.dn = dn
        self.lat = lat
        self.lon = lon
        self.alt = alt

        # Error if date is too early
        if self.dn.year < 1932:
            raise ValueError('Date cannot be before 1932!')

        # Time variables:
        self.doy = self.dn.timetuple().tm_yday
        self.utc_sec = self.dn.hour*3600. + self.dn.minute*60.
        self.utc_hour = self.dn.hour
        self.slt_hour = np.mod(self.utc_sec/3600. + self.lon/15., 24)
        self.iyd = np.mod(self.dn.year,100)*1000 + self.doy

        # For kp, ap function
        self.kp       = nan
        self.ap       = nan
        self.f107     = nan
        self.f107a    = nan
        self.f107p    = nan # previous day's F10.7
        self.kp_daily = nan
        self.ap_daily = nan
        self.apmsis   = [nan,]*7
        self.dst      = nan
        self.ae       = nan

        # For iri:
        self.ne = nan
        ions = ['O+', 'H+', 'HE+', 'O2+', 'NO+']
        self.ni={}
        for ion in ions:
            self.ni[ion] = nan

        self.Ti = nan
        self.Te = nan
        self.Tn_iri = nan

        self.NmF2 = nan
        self.hmF2 = nan

        # For msis:
        self.Tn_msis = nan
        self.nn = {}
        for neutral in ['HE','O','N2','O2','AR','H','N','O_anomalous']:
            self.nn[neutral] = nan
        self.rho = nan

        # For hwm 93/07:
        self.u = nan
        self.v = nan
        self.hwm_version = nan

        # For igrf:
        self.Bx  = nan
        self.By  = nan
        self.Bz  = nan
        self.B   = nan
        self.dip = nan
        self.dec = nan

        # For run_airglow:
        self.ag6300 = nan
        self.ag7774 = nan

        # Flag for user indices:
        self.user_ind = user_ind

        if not self.user_ind:
            # Call the indice models:
            self.get_indices()
            self.apmsis = get_apmsis(self.dn)


    def get_indices(self):
        """
        Retreives geophysical indices.
        """
        self.kp, self.ap, self.f107, self.f107a, self.f107p, \
                self.kp_daily, self.ap_daily, self.dst, self.ae  \
                        = get_kpap(self.dn)
        return self


    @staticmethod
    def init_iri16():
        """
        If required (depending on the global variable *__INIT_IRI16*),
        initialize IRI 2016. Return `True` if the model was
        initialized and `False` otherwise.
        """
        if not globals()['__INIT_IRI16']:
            read_ig_rz()
            readapf107()
            globals()['__INIT_IRI16'] = True
            return True
        else:
            return False


    def run_iri(
        self,
        NmF2=None,
        hmF2=None,
        version=2016,
        compute_Ne=True,
        compute_Te_Ti=True,
        compute_Ni=True,
        debug=False,
    ):
        """
        Run IRI model at point time/location and update the object state
        accordingly. If *NmF2* (in [cm^{-3}}]) or *hmF2* (in [km]) are
        specified, input them to the model (see documentation for
        IRI_SUB)). Override the model with *version* --- valid options
        are currently 2016 or 2012. Output debugging information if
        *debug* is true. The toggles *compute_Ne*, *compute_Te_Ti*,
        and *compute_Ni* control, respectively, whether electron
        density, electron and ion temperatures, and ion density are
        computed (restricting the model to only what is required can
        reduce run time) or set to `NaN`.
        """

        if version==2016:
            iri_data_stub = 'iri16_data/'
            iri = iri16
            init_iri = Point.init_iri16
        elif version==2012:
            iri_data_stub = 'iri12_data/'
            iri = iri12
            init_iri = lambda: False
        else:
            raise ValueError(
                "Invalid version of {} for IRI.\n".format(version) + \
                        "Either 2016 (default) or 2012 is valid."
            )

        if debug:
            print("Version = {}".format(version))

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

        if self.user_ind:
            # Set jf(25) switch to false (in Fortran)
            #   which is jf[24] in Python
            jf[24] = 0

            # Set jf(32) switch to false (in Fortran)
            #   which is jf[31] in Python
            jf[31] = 0

            # Store user indice for F10.7 in oarr:
            oarr[40] = self.f107

            # Store user index for F10.7 81 day average in oarr:
            oarr[45] = self.f107a

            # Reference:
            # https://github.com/timduly4/pyglow/issues/34#issuecomment-340645358

        # Get current directory:
        my_pwd = os.getcwd()

        # IRI data path.  We need to change directories
        # into where the IRI data are located in order
        # to run IRI:
        iri_data_path = os.path.join(
            DIR_FILE,
            iri_data_stub,
        )
        if debug:
            print("Changing directory to {}".format(iri_data_path))

        os.chdir(iri_data_path)
        init_iri()
        outf = iri(
            jf,
            0,
            self.lat,
            self.lon,
            int(self.dn.year),
            -self.doy,
            (self.utc_sec/3600.+25.),
            self.alt,
            self.alt+1,
            1,
            oarr,
        )
        os.chdir(my_pwd)

        if compute_Te_Ti:
            self.Te = outf[3,0] # Electron temperature from IRI (K)
            self.Ti = outf[2,0] # Ion temperature from IRI (K)
        else:
            self.Te = float('NaN')
            self.Ti = float('NaN')

        self.Tn_iri = outf[1,0] # Neutral temperature from IRI (K)

        self.ne        = outf[0,0] # Electron density (m^-3)
        self.ni['O+']  = outf[4,0] # O+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['H+']  = outf[5,0] # H+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['HE+'] = outf[6,0] # HE+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['O2+'] = outf[7,0] # O2+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['NO+'] = outf[8,0] # NO+ Density (%, or m^-3 with JF(22) = 0)

        self.NmF2 = oarr[0]
        self.hmF2 = oarr[1]

        if compute_Ne:
            self.ne    = outf[0,0] # Electron density (m^-3)
        else:
            self.ne    = float('NaN')

        # Densities are now in cm^-3:
        self.ne        = self.ne        / 100.**3 # [items/cm^3]
        self.ni['O+']  = self.ni['O+']  / 100.**3 # [items/cm^3]
        self.ni['H+']  = self.ni['H+']  / 100.**3 # [items/cm^3]
        self.ni['HE+'] = self.ni['HE+'] / 100.**3 # [items/cm^3]
        self.ni['O2+'] = self.ni['O2+'] / 100.**3 # [items/cm^3]
        self.ni['NO+'] = self.ni['NO+'] / 100.**3 # [items/cm^3]
        self.NmF2      = self.NmF2      / 100.**3 # [items/cm^3]

        return self


    def run_msis(self, version=2000):
        """
        Method to call MSIS model
        """

        if version==2000:
            msis = msis00
        else:
            raise ValueError(
                "Invalid version of '{}' for MSIS.\n".format(version) + \
                        "2000 (default) is valid."
            )

        [d,t] = msis(
            self.doy,
            self.utc_sec,
            self.alt,
            self.lat,
            np.mod(self.lon,360),
            self.slt_hour,
            self.f107a,
            self.f107p,
            self.apmsis,
            48,
        )
        self.Tn_msis = t[1] # neutral temperature from MSIS (K)

        self.nn = {}
        self.nn['HE'] = d[0] # [items/cm^3]
        self.nn['O']  = d[1] # [items/cm^3]
        self.nn['N2'] = d[2] # [items/cm^3]
        self.nn['O2'] = d[3] # [items/cm^3]
        self.nn['AR'] = d[4] # [items/cm^3]
        # [5] is below
        self.nn['H']  = d[6] # [items/cm^3]
        self.nn['N']  = d[7] # [items/cm^3]
        self.nn['O_anomalous'] = d[8] # [items/cm^3]

        self.rho = d[5] # total mass density [grams/cm^3]

        return self


    def run_hwm(self, version=2014):
        """
        Wrapper to call various HWM models
        """
        if version==2014:
            self._run_hwm14()
        elif version==2007:
            self._run_hwm07()
        elif version==1993:
            self._run_hwm93()
        else:
            raise ValueError(
                "Invalid version of {} for HWM.\n".format(version) +\
                        "Either 2014 (default), 2007, or 1993 is valid."
            )

        return self


    def _run_hwm93(self):
        """
        HWM 1993 Climatological model.

        """

        w = hwm93(
            self.iyd,
            self.utc_sec,
            self.alt,
            self.lat,
            np.mod(self.lon, 360),
            self.slt_hour,
            self.f107a,
            self.f107,
            self.ap_daily,
        )
        self.v = w[0]
        self.u = w[1]
        self.hwm_version = '93'

        return self


    def _run_hwm07(self):
        """
        HWM 2007 Climatological model.

        """

        my_pwd = os.getcwd()

        hwm07_data_path = os.path.join(DIR_FILE, "hwm07_data/")

        os.chdir(hwm07_data_path)
        aphwm07 = [float('NaN'), self.ap]
        w = hwm07(
            self.iyd,
            self.utc_sec,
            self.alt,
            self.lat,
            np.mod(self.lon,360),
            self.slt_hour,
            self.f107a,
            self.f107,
            aphwm07,
        )

        # Change back to original directory:
        os.chdir(my_pwd)
        self.v = w[0]
        self.u = w[1]
        self.hwm_version = '07'

        return self


    def _run_hwm14(self):
        """
        HWM 2014 Climatological model.

        """

        my_pwd = os.getcwd()

        hwm14_data_path = os.path.join(DIR_FILE, "hwm14_data/")

        os.chdir(hwm14_data_path)

        v, u = hwm14(
            self.iyd,
            self.utc_sec,
            self.alt,
            self.lat,
            np.mod(self.lon, 360),
            np.nan,
            np.nan,
            np.nan,
            [np.nan,self.ap],
        )

        # Change back to original directory:
        os.chdir(my_pwd)
        self.v = v
        self.u = u
        self.hwm_version = '14'

        return self


    def run_igrf(self, version=12):
        """
        Run the IGRF climatological model
        """

        if version==12:
            igrf=igrf12
        elif version==11:
            igrf=igrf11
        else:
            raise ValueError(
                "Invalid version of {} for IGRF.\n".format(version) + \
                        "Version 12 (default) and 11 are valid."
            )

        x, y, z, f = igrf(
            0,
            self.dn.year,
            1,
            self.alt,
            90.-self.lat,
            np.mod(self.lon, 360),
        )

        h = np.sqrt(x**2 + y**2)
        dip = 180./np.pi * np.arctan2(z,h)
        dec = 180./np.pi * np.arctan2(y,x)

        # Note that the changes here match
        # coordinate convention with other models
        # (i.e., HWM), that is:
        # (x -> east, y -> north, z -> up)
        #
        # IGRF gives (x -> north, y -> east, z -> down)
        warnings.warn(
            "Caution: IGRF coordinates have been recently changed to\n" + \
                    "Bx -> positive eastward\n" + \
                    "By -> positive northward\n" + \
                    " Bz -> positive upward\n"
        )
        self.Bx  =  y/1e9 # [T] (positive eastward) (note x/y switch here)
        self.By  =  x/1e9 # [T] (positive northward) (note x/y switch here)
        self.Bz  = -z/1e9 # [T] (positive upward) (note negation here)
        self.B   =  f/1e9 # [T]

        self.dip = dip
        self.dec = dec

        return self


    def run_airglow(self):
        """
        Computes airglow intensities


        After running, self.ag6300 has the 630.0-nm volume emission
        rate (ph/cm^3/s) and self.ag7774 has the 777.4-nm volume
        emission rate.


        History
        ------
        9/9/13 -- implemented into pyglow
                  based on Jonathan J. Makela's MATLAB
                  and subsequent python code
        9/13/16 -- added 7774 calculation
        """

        # Let's see if IRI and MSIS have been executed.
        # If not, run the appropriate models:
        if np.isnan(self.ne):
            self.run_iri()
        if np.isnan(self.nn['O2']):
            self.run_msis()

        # Perform 630.0-nm calculation
        Ne = self.ne        # electron density [cm^-3]
        Tn = self.Tn_msis   # neutral temperature [K]
        Ti = self.Ti        # ion temperature [K]
        Te = self.Te        # electron temperature [K]
        O2 = self.nn['O2']  # O2 density [cm^-3]
        N2 = self.nn['N2']  # N2 density [cm^-3]
        O  = self.nn['O']   # O density [cm^-3]

        te = Te/300.
        ti = Ti/300.

        # These coefs are from Link and Cogger, JGR 93(A9), 988309892, 1988
        K1_6300 = 3.23e-12*np.exp(3.72/ti - 1.87/ti**2)
        K2_6300 = 2.78e-13*np.exp(2.07/ti - 0.61/ti**2)
        K3_6300 = 2.0e-11*np.exp(111.8/Tn)
        K4_6300 = 2.9e-11*np.exp(67.5/Tn)
        K5_6300 = 1.6e-12*Te**0.91
        b6300 = 1.1
        a1D = 7.45e-3    # Corrected value form Link and Cogger, JGR, 94(A2), 1989
        a6300 = 5.63e-3  # Corrected value form Link and Cogger, JGR, 94(A2), 1989

        # Calculate O+ assuming mixture of ions (also from Link and Cogger, 1988)
        a1 = 1.95e-7*te**-0.7
        a2 = 4.00e-7*te**-0.9
        Oplus = Ne/(1.+K2_6300*N2/a2/Ne + K1_6300*O2/a1/Ne)

        AGNumerator = a6300/a1D*b6300*K1_6300*Oplus*O2
        AGDenominator = 1.+(K3_6300*N2+K4_6300*O2+K5_6300*Ne)/a1D
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
        V7774_rr = alpha1_7774 * Oplus * Ne

        V7774_ii_num = beta_7774 * K1_7774 * K2_7774 * O * Oplus * Ne
        V7774_ii_den = K2_7774 * Oplus + K3_7774 * O

        self.ag7774 = V7774_rr + V7774_ii_num / float(V7774_ii_den)

        return self

# ---------

def _igrf_tracefield(dn, lat, lon, alt, target_ht, step):
    """
    Helper function to trace along a magnetic field line using IGRF


    :param dn: datetime.datetime object of requested trace
    :param lat: Latitude [degrees]
    :param lon: Longitude [degrees]
    :param alt: Altitude [km]
    :param target_ht: Altitude to stop trace [km]
    :param step: Step size of trace [km]

    :return lla: (latitude, longitude, altitude) data structure of trace
                 with dimentions [Nsteps x 3]

    """

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
    """
    Helper function to trace along a magnetic field line using IGRF
    for only one hemisphere

    :param dn: datetime.datetime object of requested trace
    :param lat: Latitude [degrees]
    :param lon: Longitude [degrees]
    :param alt: Altitude [km]
    :param target_ht: Altitude to stop trace [km]
    :param step: Step size of trace [km]

    :return lla: (latitude, longitude, altitude) data structure of trace
                 with dimentions [Nsteps x 3]

    """

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
        ecef_new =  ecef + coord.ven2ecef(
            lla,
            [-D/A*step, E/A*step, N/A*step],
        )

        # Convert to lla coordinates:
        lla = coord.ecef2lla(ecef_new)

        # add the field line to our collection:
        lla_field = np.vstack([lla_field, lla])
        i = i + 1

    """ Step 2: Make the last point close to target_ht """
    while (abs(lla[2]-target_ht) > TOLERANCE):

        my_error = lla[2]-target_ht
        old_step = step
        # Find out how much we need to step by:
        step = -np.sign(step)*abs(lla[2]-target_ht)

        ecef = coord.lla2ecef(lla)

        p = Point(dn, lla[0], lla[1], lla[2]/1e3)
        p.run_igrf()

        N = p.Bx # North
        E = p.By # East
        D = p.Bz # Down
        A = p.B  # Total

        # Trace the field, but use the modified step:
        ecef_new = ecef + coord.ven2ecef(
            lla,
            np.array([-D/A, E/A, N/A]) * step/(-D/A),
        )

        # TODO : I changed this, is this correct?
        lla = coord.ecef2lla(ecef_new)

    # replace last entry with the point close to target_ht:
    lla_field[-1,:] = lla

    return lla_field

def Line(dn, lat, lon, alt, target_ht=90., step=15.):
    """
    Return a list of instances of Point by
    tracing along the geomagnetic field line.

    pts = Line(dn, lat, lon, alt, target_ht=90., step=15.)

    :param dn: datetime.datetime object of requested trace
    :param lat: Latitude [degrees]
    :param lon: Longitude [degrees]
    :param alt: Altitude [km]
    :param target_ht: Altitude to stop trace [km]
    :param step: Step size of trace [km]

    """
    llas = _igrf_tracefield(dn, lat, lon, alt, target_ht, step)
    pts = []

    for lla in llas:
        pts.append(
            Point(dn, lla[0], lla[1], lla[2]/1e3)
        )

    return pts


def update_kpap(years=None):
    '''
    Update the Kp and Ap indices used in pyglow.
    The files will be downloaded from noaa to your pyglow
    installation directory.

    update_kpap(years=None)

    :param years: (optional) a list of years to download.
            If this input is not provided, the full
            range of years starting from 1932 to the
            current year will be downloaded.
    '''

    # Load all data up until today
    if years is None:
        years = range(1932, date.today().year + 1)[::-1] # reverse

    pyglow_dir = os.path.join(DIR_FILE, "kpap/")

    for year in years:
        src = 'ftp://ftp.ngdc.noaa.gov/'\
                + 'STP/GEOMAGNETIC_DATA/INDICES/KP_AP/%4i'\
                % (year,)
        des = pyglow_dir + "%4i" % (year,)
        print("\nDownloading\n{src}\nto\n{des}".format(src=src, des=des))
        try:
            with contextlib.closing(urllib.request.urlopen(src)) as r:
                with open(des, 'wb') as f:
                    shutil.copyfileobj(r, f)
        except IOError as e:
            print(
                'Failed downloading data for year {}. File does not exist'.format(year)
            )


def update_dst(years=None):
    """
    Update the Dst index files used in pyglow.
    The files will be downloaded from WDC Kyoto
    to your pyglow installation directory.

    update_dst(years=None)

    :param years: (optional) a list of years to download.
            If this input is not provided, the full
            range of years starting from 2005 to the
            current year will be downloaded. Pre-2005
            files are shipped with pyglow.
    """


    def download_dst(year, month, des):
        """
        Helper function to earch for the appropriate location
        and download the DST index file from WDC Kyoto for the
        given month and year. Save it to the specified file "des".
        Return True if successful, False if not.
        """
        # There are three possible sources of data. Search for
        # them in the following order:
        # 1) Final
        # 2) Provisional
        # 3) Realtime
        year_month = '%i%02i' % (year, month)
        wgdc_fn = 'dst%s%02i.for.request' % (str(year)[2:],month)
        src_final       = 'http://wdc.kugi.kyoto-u.ac.jp/dst_final/%s/%s' % (year_month, wgdc_fn)
        src_provisional = 'http://wdc.kugi.kyoto-u.ac.jp/dst_provisional/%s/%s' % (year_month, wgdc_fn)
        src_realtime    = 'http://wdc.kugi.kyoto-u.ac.jp/dst_realtime/%s/%s' % (year_month, wgdc_fn)

        success = False
        for src in [src_final, src_provisional, src_realtime]:
            try:
                with contextlib.closing(urllib.request.urlopen(src)) as r:
                    contents = r.read().decode('utf8')
                    # If that succeeded, then the file exists
                    print("\nDownloading\n{src}\nto\n{des}".format(src=src, des=des))
                    with open(des, 'w') as f:
                        f.write(contents)
                    success = True
                    break
            except urllib.error.HTTPError:
                pass
        return success

    # Read files from 2005 until today. Pre-2005
    # files are shipped with pyglow.
    if years is None:
        years = list(range(2005, date.today().year + 1))[::-1] # reversed
    pyglow_dir = os.path.join(DIR_FILE, "dst/")

    for year in years:
        for month in range(1, 13):
            des = '%s%i%02i' % (pyglow_dir,year,month)
            download_dst(year, month, des)

    return


def update_ae(years = None):
    '''
    Update the AE index files used in pyglow.
    The files will be downloaded from WDC Kyoto
    to your pyglow installation directory.

    update_ae(years=None)

    :param years: (optional) a list of years to download.
            If this input is not provided, the full
            range of years starting from 2005 to the
            current year will be downloaded. Pre-2005
            files are shipped with pyglow.
    '''

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
        #src_final       = 'http://wdc.kugi.kyoto-u.ac.jp/ae_final/%s/%s' % (year_month, wgdc_fn)
        src_provisional = 'http://wdc.kugi.kyoto-u.ac.jp/ae_provisional/%s/%s' %\
                (year_month, wgdc_fn)
        src_realtime    = 'http://wdc.kugi.kyoto-u.ac.jp/ae_realtime/%s/%s' %\
                (year_month, wgdc_fn)

        success = False
        for src in [src_provisional, src_realtime]:
            try:
                with contextlib.closing(urllib.request.urlopen(src)) as r:
                    contents = r.readlines()
                    # If that succeeded, then the file exists
                    print("\nDownloading\n{src}\nto\n{des}".format(src=src, des=des))
                    with open(des,'w') as f:
                        # this shrinks the filesize to hourly
                        for c in contents:
                            c = c.decode('utf8')
                            f.write("%s%s%s\n"%(c[12:18],c[19:21],c[394:400]))
                    success = True
                    break
            except urllib.error.HTTPError:
                pass
        return success

    # Read files from 2000 until today. Pre-2000
    # files are shipped with pyglow.
    if years is None:
        years = list(range(2000, date.today().year + 1))[::-1] # reversed
    pyglow_dir = os.path.join(DIR_FILE, "ae/")

    for year in years:
        for month in range(1,13):
            des = '%s%i%02i' % (pyglow_dir,year,month)
            download_ae(year, month, des)


def update_indices(years=None):
    '''
    Update all geophysical indices (e.g., KP, DST, AE).

    update_indices(years=None)

    :param years: (optional) a list of years to download.
            If this input is not provided, default
            values will be used.
    '''

    update_kpap(years=years)
    update_dst(years=years)
    update_ae(years=years)

    return

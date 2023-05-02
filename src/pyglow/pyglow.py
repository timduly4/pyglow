from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import str  # noqa E402
from builtins import range  # noqa E402
from builtins import object  # noqa E402
from past.utils import old_div  # noqa E402
import numpy as np  # noqa E402
import warnings  # noqa E402

from . import coord  # noqa E402
from .location_time import LocationTime  # noqa E402
from . import constants  # noqa E402
from .iri import IRI  # noqa E402
from .hwm import HWM  # noqa E402
from .msis import MSIS  # noqa E402
from .igrf import IGRF  # noqa E402
from .geophysical_indices import Indice  # noqa E402

# Code version:
__version__ = constants.VERSION


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

        # Error if date is too early:
        if self.dn.year < 1932:
            raise ValueError('Date cannot be before 1932!')

        # Form location/time data structure:
        self.location_time = LocationTime(dn, lat, lon, alt)

        # Initialize IRI:
        self.iri = IRI()
        self.ne = self.iri.ne
        self.ni = self.iri.ni
        self.Ti = self.iri.Ti
        self.Te = self.iri.Te
        self.Tn_iri = self.iri.Tn
        self.NmF2 = self.iri.NmF2
        self.hmF2 = self.iri.hmF2

        # Flag for user indices:
        self.user_ind = user_ind

        # Geophysical indices:
        self.indice = Indice(dn)
        if not self.user_ind:
            # Call the indice models:
            self.indice.run()

            # Post warning if all the indices are NaN:
            if self.indice.all_nan():
                warnings.warn(
                    "All geophysical indices are NaN. "
                    "Consider updating the indices."
                )

        # Assign indice data:
        self.kp = self.indice.kp
        self.ap = self.indice.ap
        self.f107 = self.indice.f107
        self.f107a = self.indice.f107a
        self.f107p = self.indice.f107p
        self.kp_daily = self.indice.kp_daily
        self.ap_daily = self.indice.ap_daily
        self.dst = self.indice.dst
        self.ae = self.indice.ae
        self.ap1 = self.indice.ap1
        self.apmsis = self.indice.apmsis

        # For msis:
        self.msis = MSIS()
        self.Tn_msis = self.msis.Tn
        self.nn = self.msis.nn
        self.rho = self.msis.rho

        # For hwm 93/07/14:
        self.hwm = HWM()
        self.u = self.hwm.u
        self.v = self.hwm.v
        self.hwm_version = self.hwm.hwm_version
        self.hwm_dwm = self.hwm.hwm_dwm

        # For igrf:
        self.igrf = IGRF()
        self.Bx = self.igrf.Bx
        self.By = self.igrf.By
        self.Bz = self.igrf.Bz
        self.B = self.igrf.B
        self.dip = self.igrf.dip
        self.dec = self.igrf.dec

        # For run_airglow:
        self.ag6300 = nan
        self.ag7774 = nan

        return

    def __str__(self):
        """ String representation of pyglow class """

        pyglow_str = "pyglow.Point: dn = {dn}, lat = {lat:3.2f} [deg], "\
                     "lon = {lon:3.2f} [deg], alt = {alt:3.2f} [km]".format(
                        dn=self.dn.strftime("%Y-%m-%d %H:%M:%S"),
                        lat=self.lat,
                        lon=self.lon,
                        alt=self.alt,
                     )

        return pyglow_str

    def __repr__(self):
        """ Representation value of pyglow class """

        pyglow_repr = "pyglow.Point({dn}, {lat}, {lon}, {alt})".format(
                        dn=self.dn.__repr__(),
                        lat=self.lat,
                        lon=self.lon,
                        alt=self.alt,
                     )

        return pyglow_repr

    def run_iri(
        self,
        version=2016,
        NmF2=None,
        hmF2=None,
        compute_Ne=True,
        compute_Te_Ti=True,
        compute_Ni=True,
        f1_layer=True,
        bil2000=False,
    ):
        """
        Executes IRI and assigns results to instance.

        :param version: Version of IRI to run
        :param NmF2: User-specified NmF2 [cm^-3]
        :param hmF2: User-specified hmF2 [km]
        :param compute_Ne: Switch to compute Ne
        :param compute_Te_Ti: Switch to compute Te and Ti
        :param compute_Ni: Switch to compute Ni
        :param f1_layer: If True (default) include F1-layer (JF switches 19,20 = True)
        :param bil2000: If True, use Bil-2000 model for bottomside (JF switch 4). Default IRI is False
        """

        # Check if user supplies indices:
        if self.user_ind:
            f107 = self.f107
            f107a = self.f107a
            if np.isnan(f107) or np.isnan(f107a):
                raise ValueError(
                    "Cannot assign f107 or f017a to NaN when executing IRI"
                )
        else:
            f107 = None
            f107a = None

        # Execute IRI:
        self.iri.run(
            self.location_time,
            version=version,
            NmF2=NmF2,
            hmF2=hmF2,
            compute_Ne=compute_Ne,
            compute_Te_Ti=compute_Te_Ti,
            compute_Ni=compute_Ni,
            f107=f107,
            f107a=f107a,
            f1_layer=f1_layer,
            bil2000=bil2000,
        )

        # Assign output of IRI:
        self.ne = self.iri.ne
        self.ni = self.iri.ni
        self.Ti = self.iri.Ti
        self.Te = self.iri.Te
        self.Tn_iri = self.iri.Tn
        self.NmF2 = self.iri.NmF2
        self.hmF2 = self.iri.hmF2

        return self

    def run_hwm(self, version=2014, dwm = 'on'):
        """
        Executes HWM and assigns results to instance.

        Call examples:
        - run_hwm(version=2014, dwm='on')  # Default setting: Standard DWM
                                           # this can cause discontinuities
                                           # in the wind due to the usage of
                                           # the 3-hour ap to drive DWM.

        - run_hwm(version=2014, dwm='off') # Default, with ap=-1



        - run_hwm(version=2014, dwm='smooth')   # This is our modified DWM:
                                                # It interpolates the 3-hour
                                                # ap index into a smooth curve
                                                # using a cubic interpolation.

        :param version: Version of HWM to run
        """

        # Run HWM:
        self.hwm.run(
            self.location_time,
            version,
            dwm,
            f107=self.f107,
            f107a=self.f107a,
            ap=self.ap,
            ap_daily=self.ap_daily,
            ap1=self.ap1
        )

        # Assign output:
        self.u = self.hwm.u
        self.v = self.hwm.v
        self.hwm_version = self.hwm.hwm_version
        self.hwm_dwm = self.hwm.hwm_dwm

        return self

    def run_msis(self, version=2000):
        """
        Run the MSIS climatological model

        :param version: Version of MSIS to run
        """

        # Run MSIS:
        self.msis.run(
            self.location_time,
            self.f107a,
            self.f107p,
            self.apmsis,
            version=version,
            )

        # Assign output:
        self.Tn_msis = self.msis.Tn
        self.nn = self.msis.nn
        self.rho = self.msis.rho

        return self

    def run_igrf(self, version=12):
        """
        Run the IGRF climatological model

        :param version: Version of IGRF to run
        """

        # Run IGRF:
        self.igrf.run(self.location_time, version)

        # Assign output:
        self.Bx = self.igrf.Bx
        self.By = self.igrf.By
        self.Bz = self.igrf.Bz
        self.B = self.igrf.B
        self.dip = self.igrf.dip
        self.dec = self.igrf.dec

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
                  and subsequent Python code
        9/13/16 -- added 7774 calculation
        """

        # Let's see if IRI and MSIS have been executed. If not, run the
        # appropriate models:
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
        O_density = self.nn['O']   # O density [cm^-3]

        te = Te/300.
        ti = Ti/300.

        # These coefs are from Link and Cogger, JGR 93(A9), 988309892, 1988
        K1_6300 = 3.23e-12*np.exp(3.72/ti - 1.87/ti**2)
        K2_6300 = 2.78e-13*np.exp(2.07/ti - 0.61/ti**2)
        K3_6300 = 2.0e-11*np.exp(111.8/Tn)
        K4_6300 = 2.9e-11*np.exp(67.5/Tn)
        K5_6300 = 1.6e-12*Te**0.91
        b6300 = 1.1
        a1D = 7.45e-3  # Corr. value from Link and Cogger, JGR, 94(A2), 1989
        a6300 = 5.63e-3  # Corr. value form Link and Cogger, JGR, 94(A2), 1989

        # Calculate O+ assuming mixture of ions
        # (also from Link and Cogger, 1988):
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

        V7774_ii_num = beta_7774 * K1_7774 * K2_7774 * O_density * Oplus * Ne
        V7774_ii_den = K2_7774 * Oplus + K3_7774 * O_density

        self.ag7774 = V7774_rr + V7774_ii_num / float(V7774_ii_den)

        return self


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
    lla_north = _igrf_tracefield_hemis(
        dn,
        lat,
        lon,
        alt,
        target_ht,
        step,
    )

    # Go South:
    lla_south = _igrf_tracefield_hemis(
        dn,
        lat,
        lon,
        alt,
        target_ht,
        -step,
    )

    # Stack them together:
    lla = np.vstack(
        [
            np.flipud(lla_north)[:-1, :],
            lla_south,
        ]
    )

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
    TOLERANCE = 10  # [m]
    i = 0
    while (lla[2] > target_ht):
        # convert to ECEF:
        ecef = coord.lla2ecef(lla)

        # Grab field line information:
        p = Point(dn, lla[0], lla[1], lla[2]/1e3)
        p.run_igrf()

        # coordinates follow pyglow's convention:
        # (x -> east, y -> north, z -> up)

        N = p.By  # North
        E = p.Bx  # East
        D = -p.Bz  # Down
        A = p.B   # Total

        # Step along the field line
        ecef_new = ecef + coord.ven2ecef(
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

        # Find out how much we need to step by:
        step = -np.sign(step)*abs(lla[2]-target_ht)

        ecef = coord.lla2ecef(lla)

        p = Point(dn, lla[0], lla[1], lla[2]/1e3)
        p.run_igrf()

        N = p.Bx  # North
        E = p.By  # East
        D = p.Bz  # Down
        A = p.B   # Total

        # Trace the field, but use the modified step:
        ecef_new = ecef + coord.ven2ecef(
            lla,
            np.array([-D/A, E/A, N/A]) * step/(-D/A),
        )

        # TODO : I changed this, is this correct?
        lla = coord.ecef2lla(ecef_new)

    # replace last entry with the point close to target_ht:
    lla_field[-1, :] = lla

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
    :param target_ht: (optional) Altitude to stop trace [km]
    :param step: (optional) Step size of trace [km]

    """
    llas = _igrf_tracefield(dn, lat, lon, alt, target_ht, step)
    pts = []

    for lla in llas:
        pts.append(
            Point(dn, lla[0], lla[1], lla[2]/1e3)
        )

    return pts

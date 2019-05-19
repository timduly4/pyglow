import numpy as np

from msis00py import gtd7 as msis00  # noqa E402
from .constants import nan

CONSTITUENTS = ['HE', 'O', 'N2', 'O2', 'AR', 'H', 'N', 'O_anomalous']


class MSIS(object):

    def __init__(self):

        # Initialize member variables:
        self.nn = {constituent: nan for constituent in CONSTITUENTS}
        self.Tn = nan
        self.rho = nan

        pass

    def run(self, location_time, f107a, f107p, apmsis, version=2000):
        """
        Method to call MSIS model

        :param location_time: Instance of LocationTime
        :param f107: f107 indice
        :param f107a: f107a indice
        :param apmsis: ap indice array for MSIS
        """

        if version == 2000:
            msis = msis00
        else:
            raise ValueError(
                "Invalid version of '{}' for MSIS.\n".format(version) +
                "2000 (default) is valid."
            )

        # Call MSIS:
        [d, t] = msis(
            location_time.doy,
            location_time.utc_sec,
            location_time.alt,
            location_time.lat,
            np.mod(location_time.lon, 360),
            location_time.slt_hour,
            f107a,
            f107p,
            apmsis,
            48,
        )

        # Neutral temperature:
        self.Tn = t[1]  # neutral temperature from MSIS (K)

        # Constituent densities:
        self.nn['HE'] = d[0]  # [items/cm^3]
        self.nn['O'] = d[1]  # [items/cm^3]
        self.nn['N2'] = d[2]  # [items/cm^3]
        self.nn['O2'] = d[3]  # [items/cm^3]
        self.nn['AR'] = d[4]  # [items/cm^3]
        # (FYI, [5] is below)
        self.nn['H'] = d[6]  # [items/cm^3]
        self.nn['N'] = d[7]  # [items/cm^3]
        self.nn['O_anomalous'] = d[8]  # [items/cm^3]

        # Total mass density:
        self.rho = d[5]  # total mass density [grams/cm^3]

        return self

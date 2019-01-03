
import numpy as np


class LocationTime(object):

    def __init__(self, dn, lat, lon, alt):

        self.dn = dn
        self.lat = lat
        self.lon = lon
        self.alt = alt

        # Year:
        self.year = dn.year

        # Day of year:
        self.doy = dn.timetuple().tm_yday

        # UTC seconds:
        self.utc_sec = dn.hour*3600. + dn.minute*60.

        # UTC Hour:
        self.utc_hour = dn.hour

        # Solar local time hour:
        self.slt_hour = np.mod(self.utc_sec/3600. + self.lon/15., 24)

        # Integer year and doy, e.g. 2018039
        self.iyd = np.mod(dn.year, 100)*1000 + self.doy

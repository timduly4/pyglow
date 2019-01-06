import numpy as np

from igrf11py import igrf11syn as igrf11
from igrf12py import igrf12syn as igrf12
from .constants import nan


class IGRF(object):

    def __init__(self):

        # Member variables:
        self.Bx = nan
        self.By = nan
        self.Bz = nan
        self.B = nan
        self.dip = nan
        self.dec = nan

    def run(self, location_time, version):
        """
        Run the IGRF climatological model

        :param location_time: Instance of LocationTime
        :param version: Version of IGRF to run
        """

        if version == 12:
            igrf = igrf12
        elif version == 11:
            igrf = igrf11
        else:
            raise ValueError(
                "Invalid version of {} for IGRF.\n".format(version) +
                "Version 12 (default) and 11 are valid."
            )

        # Run the IGRF model:
        x, y, z, f = igrf(
            0,
            location_time.dn.year,
            1,
            location_time.alt,
            90.-location_time.lat,
            np.mod(location_time.lon, 360),
        )

        # Compute dip and declination angles:
        h = np.sqrt(x**2 + y**2)
        dip = 180./np.pi * np.arctan2(z, h)
        dec = 180./np.pi * np.arctan2(y, x)

        # Note that the changes here match coordinate convention with other
        # models (i.e., HWM), that is:
        #
        #   (x -> east, y -> north, z -> up)
        #
        # IGRF gives (x -> north, y -> east, z -> down)
        #

        # Assign output:
        self.Bx = y/1e9  # [T] (positive eastward) (note x/y switch here)
        self.By = x/1e9  # [T] (positive northward) (note x/y switch here)
        self.Bz = -z/1e9  # [T] (positive upward) (note negation here)
        self.B = f/1e9  # [T]

        self.dip = dip
        self.dec = dec

        return self

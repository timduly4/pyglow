from datetime import datetime
import math
import unittest

from src.pyglow import IGRF
from src.pyglow import LocationTime
from src.pyglow.constants import DIR_FILE as pyglow_file
print("pyglow file: {}".format(pyglow_file))


class TestIgrf(unittest.TestCase):

    def setUp(self):

        self.igrf = IGRF()

        # Set up LocationTime instance:
        dn = datetime(2010, 3, 23, 15, 30)
        lat = 30
        lon = -80
        alt = 250
        self.location_time = LocationTime(dn, lat, lon, alt)

    def tearDown(self):

        pass

    def test_igrf_run(self):
        """ Simple IGRF run """

        for version in [11, 12]:
            # Run IGRF:
            self.igrf.run(self.location_time, version)

            print(self.igrf.Bx)

            # Make sure we have a IGRF result:
            self.assert_igrf_result(self.igrf)

    def assert_igrf_result(self, igrf):
        """ Ensures that we have an IGRF result """

        self.assertFalse(math.isnan(igrf.Bx))
        self.assertFalse(math.isnan(igrf.By))
        self.assertFalse(math.isnan(igrf.Bz))
        self.assertFalse(math.isnan(igrf.B))
        self.assertFalse(math.isnan(igrf.dip))
        self.assertFalse(math.isnan(igrf.dec))

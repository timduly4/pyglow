from datetime import datetime
import math
import unittest

from src.pyglow import HWM
from src.pyglow import LocationTime
from src.pyglow.constants import DIR_FILE as pyglow_file
print("pyglow file: {}".format(pyglow_file))


class TestHwm(unittest.TestCase):

    def setUp(self):

        # Instantiate:
        self.hwm = HWM()

        # Set up LocationTime instance:
        dn = datetime(2010, 3, 23, 15, 30)
        lat = 30
        lon = -80
        alt = 250
        self.location_time = LocationTime(dn, lat, lon, alt)

        # Geophysical indices:
        self.ap = 2.0
        self.ap_daily = 1.0
        self.f107 = 80
        self.f107a = 80

    def tearDown(self):
        pass

    def test_run_all_versions(self):
        """ Simply HWM 93 run """

        for version in [1993, 2007, 2014]:

            # Run HWM:
            self.hwm.run(self.location_time, version,
                         f107=self.f107, f107a=self.f107a,
                         ap=self.ap, ap_daily=self.ap_daily)

            # Make sure we have a result:
            self.assert_hwm_result(self.hwm)

    def assert_hwm_result(self, hwm):
        """ Make sure we have a HWM result """

        self.assertFalse(math.isnan(hwm.u))
        self.assertFalse(math.isnan(hwm.v))
        self.assertTrue(hwm.hwm_version)

from datetime import datetime
import math
import unittest

from src.pyglow import HWM
from src.pyglow import LocationTime
from src.pyglow import Indice
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

        self.indice = Indice(dn)
        # For simplicity, we can set the indices manually:
        self.indice.ap = 2.0
        self.indice.ap_daily = 1.0
        self.indice.f107 = 80
        self.indice.f107a = 80

    def tearDown(self):
        pass

    def test_run_93(self):
        """ Simply HWM 93 run """

        self.hwm.run(self.location_time, self.indice, 1993)

        # Make sure we have a result:
        self.assert_hwm_result(self.hwm)

    def test_run_07(self):
        """ Simply HWM 07 run """

        self.hwm.run(self.location_time, self.indice, 2007)

        # Make sure we have a result:
        self.assert_hwm_result(self.hwm)

    def test_run_14(self):
        """ Simply HWM 14 run """

        self.hwm.run(self.location_time, self.indice, 2014)

        # Make sure we have a result:
        self.assert_hwm_result(self.hwm)

    def assert_hwm_result(self, hwm):
        """ Make sure we have a HWM result """

        self.assertFalse(math.isnan(hwm.u))
        self.assertFalse(math.isnan(hwm.v))
        self.assertTrue(hwm.hwm_version)

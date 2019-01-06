from datetime import datetime
import math
import unittest

from src.pyglow import MSIS
from src.pyglow import LocationTime
from src.pyglow import Indice
from src.pyglow.constants import DIR_FILE as pyglow_file
print("pyglow file: {}".format(pyglow_file))


class TestMSIS(unittest.TestCase):

    def setUp(self):

        self.msis = MSIS()

        # Set up LocationTime instance:
        dn = datetime(2010, 3, 23, 15, 30)
        lat = 30
        lon = -80
        alt = 250
        self.location_time = LocationTime(dn, lat, lon, alt)

        # Helper for indices:
        indice = Indice(dn)
        indice.run()

        # Indices used for MSIS call:
        self.f107a = indice.f107a
        self.f107p = indice.f107p
        self.apmsis = indice.apmsis

    def tearDown(self):

        pass

    def test_msis_run(self):

        # Run the MSIS climatological model:
        self.msis.run(
            self.location_time,
            self.f107a,
            self.f107p,
            self.apmsis,
        )

        # Make sure we have a MSIS result:
        self.assert_msis_result(self.msis)

    def assert_msis_result(self, msis):

        self.assertFalse(math.isnan(msis.Tn))
        for constituent in msis.nn:
            self.assertFalse(math.isnan(msis.nn[constituent]))
        self.assertFalse(math.isnan(msis.rho))

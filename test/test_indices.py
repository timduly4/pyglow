
from datetime import datetime
import math
import unittest

from src.pyglow import Indice
from src.pyglow.constants import DIR_FILE as pyglow_file
print("pyglow file: {}".format(pyglow_file))


class TestIndice(unittest.TestCase):

    def setUp(self):

        # Instantiate indice data structure:
        dn = datetime(2010, 3, 23, 15, 30)
        self.indice = Indice(dn)

    def tearDown(self):

        pass

    def test_run(self):
        """ Retrieval of geophysical indices """

        self.indice.run()
        self.assertFalse(math.isnan(self.indice.f107))

    def test_all_nan(self):

        # Nominal case:
        self.indice.run()
        self.assertFalse(self.indice.all_nan())

    def test_all_nan_not_run(self):

        # Not running indices:
        self.assertTrue(self.indice.all_nan())

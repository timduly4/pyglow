""" Unit tests for pyglow """

from datetime import datetime
import unittest

from ipdb import set_trace as db
from pprint import pprint

import pyglow


class TestPyglow(unittest.TestCase):

    def setUp(self):

        dn = datetime(2017, 3, 23, 12)
        lat = 40
        lon = -88
        alt = 250
        self.pt = pyglow.Point(dn, lat, lon, alt)

    def tearDown(self):
        pass

    def test_version(self):

        # Make sure we have a version associated with pyglow:
        version = pyglow.__version__
        print("pyglow: v{}".format(version))

        self.assertTrue(version)

    def test_string_representation(self):

        print("'str' representation = {}".format(self.pt))
        print("'repr' representation = {}".format(self.pt.__repr__()))

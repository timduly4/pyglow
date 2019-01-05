""" Unit tests for pyglow """

from datetime import datetime
import math
import unittest

from ipdb import set_trace as db
from pprint import pprint

from src import pyglow
print("pyglow file: {}".format(pyglow.__file__))


class TestPoint(unittest.TestCase):

    def setUp(self):

        dn = datetime(2007, 3, 23, 12)
        lat = 40
        lon = -88
        alt = 250
        self.pt = pyglow.Point(dn, lat, lon, alt)

        self.pt_user_ind = pyglow.Point(dn, lat, lon, alt, user_ind=True)
        self.pt_user_ind.f107 = 100
        self.pt_user_ind.f107a = 100

    def tearDown(self):
        pass

    def test_version(self):
        """ Version associated with pyglow """

        # Make sure we have a version associated with pyglow:
        version = pyglow.__version__
        print("pyglow: v{}".format(version))

        self.assertTrue(version)

    def test_string_representation(self):
        """ String representations of the Point class """

        print("'str' representation = {}".format(self.pt))
        print("'repr' representation = {}".format(self.pt.__repr__()))

    def test_run_iri(self):
        """ Interface to IRI """

        # Nominal run:

        # Execute IRI
        self.pt.run_iri()

        # Make sure we have an IRI result:
        self.assert_iri_result(self.pt)

        # User specified indices:
        self.pt_user_ind.run_iri()

        # Make sure we have an IRI result:
        self.assert_iri_result(self.pt_user_ind)

        # Make sure the results from the default and user-defined indices
        # differ:
        self.assertNotEqual(
            self.pt.ne,
            self.pt_user_ind.ne,
        )

    def assert_iri_result(self, pt):
        """ Ensures that we have an IRI result """

        self.assertFalse(math.isnan(pt.ne))
        ions = pt.ni.keys()
        for ion in ions:
            self.assertFalse(math.isnan(pt.ni[ion]))
        self.assertFalse(math.isnan(pt.Ti))
        self.assertFalse(math.isnan(pt.Te))
        self.assertFalse(math.isnan(pt.Tn_iri))
        self.assertFalse(math.isnan(pt.NmF2))
        self.assertFalse(math.isnan(pt.hmF2))

    def test_run_hwm(self):
        """ Interface to HWM """

        # Nominal run:

        # Execute IRI
        for version in [1993, 2007, 2014]:
            self.pt.run_hwm(version=version)

            # Make sure we have an IRI result:
            self.assert_hwm_result(self.pt)

    def assert_hwm_result(self, pt):
        """ Ensures that we have an hwm result """

        self.assertFalse(math.isnan(pt.u))
        self.assertFalse(math.isnan(pt.v))
        self.assertTrue(pt.hwm_version)

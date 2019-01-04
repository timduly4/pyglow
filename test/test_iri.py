from datetime import datetime
import math
import unittest

from src.pyglow import IRI
from src.pyglow import LocationTime
from src.pyglow.constants import DIR_FILE as pyglow_file
print("pyglow file: {}".format(pyglow_file))


class TestIri(unittest.TestCase):

    def setUp(self):

        self.iri = IRI()

        # Set up LocationTime instance:
        dn = datetime(2010, 3, 23, 15, 30)
        lat = 30
        lon = -80
        alt = 250
        self.location_time = LocationTime(dn, lat, lon, alt)

    def tearDown(self):

        pass

    def test_iri_run(self):
        """ Simple IRI run """

        # Run IRI:
        self.iri.run(self.location_time)

        # Make sure we have a IRI result:
        self.assert_iri_result(self.iri)

    def test_iri_versions(self):
        """ Versions of IRI """

        for version in [2012, 2016]:

            # Run IRI:
            self.iri.run(self.location_time, version=version)

            # Make sure we have a IRI result:
            self.assert_iri_result(self.iri)

    def test_iri_nmf2(self):
        """ Input NmF2 of IRI """

        # Run IRI:
        self.iri.run(self.location_time, NmF2=467145.0)

        # Make sure we have a IRI result:
        self.assert_iri_result(self.iri)

    def test_iri_hmf2(self):
        """ Input hmF2 of IRI """

        # Run IRI:
        self.iri.run(self.location_time, hmF2=300.0)

        # Make sure we have a IRI result:
        self.assert_iri_result(self.iri)

    def test_iri_compute_ne(self):
        """ compute_Ne switch """

        # Run IRI:
        self.iri.run(self.location_time, compute_Ne=False)

        # Ne should not have a result:
        self.assertTrue(math.isnan(self.iri.ne))

    def test_iri_compute_te_ti(self):
        """ compute_Te_Ti switch """

        # Run IRI:
        self.iri.run(self.location_time, compute_Te_Ti=False)

        # Te/Ti should not have a result:
        self.assertTrue(math.isnan(self.iri.Te))
        self.assertTrue(math.isnan(self.iri.Ti))

    def test_iri_compute_ni(self):
        """ compute_Ni switch """

        # Run IRI:
        self.iri.run(self.location_time, compute_Ni=False)

        # Ni should not have a result:
        ions = self.iri.ni.keys()
        for ion in ions:
            self.assertTrue(math.isnan(self.iri.ni[ion]))

    def assert_iri_result(self, iri):
        """ Ensures that we have an IRI result """

        self.assertFalse(math.isnan(iri.ne))
        ions = iri.ni.keys()
        for ion in ions:
            self.assertFalse(math.isnan(iri.ni[ion]))
        self.assertFalse(math.isnan(iri.Ti))
        self.assertFalse(math.isnan(iri.Te))
        self.assertFalse(math.isnan(iri.Tn))
        self.assertFalse(math.isnan(iri.NmF2))
        self.assertFalse(math.isnan(iri.hmF2))

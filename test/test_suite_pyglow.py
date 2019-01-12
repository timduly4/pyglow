import unittest

from .test_iri import TestIri
from .test_hwm import TestHwm
from .test_msis import TestMSIS
from .test_point import TestPoint
from .test_indices import TestIndice

suite = unittest.TestSuite()

# Tests from each module:
suite.addTest(unittest.makeSuite(TestIri))
suite.addTest(unittest.makeSuite(TestHwm))
suite.addTest(unittest.makeSuite(TestMSIS))
suite.addTest(unittest.makeSuite(TestPoint))
suite.addTest(unittest.makeSuite(TestIndice))

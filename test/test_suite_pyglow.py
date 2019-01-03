import unittest

from . import test_iri
from . import test_point

loader = unittest.TestLoader()
suite = unittest.TestSuite()

# Tests from each module:
suite.addTests(loader.loadTestsFromModule(test_iri))
suite.addTests(loader.loadTestsFromModule(test_point))

runner = unittest.TextTestRunner(verbosity=0)
result = runner.run(suite)

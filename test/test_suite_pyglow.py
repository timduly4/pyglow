import unittest

from . import test_iri
from . import test_point
from . import test_indices

loader = unittest.TestLoader()
suite = unittest.TestSuite()

# Tests from each module:
suite.addTests(loader.loadTestsFromModule(test_iri))
suite.addTests(loader.loadTestsFromModule(test_point))
suite.addTests(loader.loadTestsFromModule(test_indices))

runner = unittest.TextTestRunner(verbosity=0)
result = runner.run(suite)

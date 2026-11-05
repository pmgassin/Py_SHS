#!/usr/bin/env python3
"""Run the whole PySHS test-suite (standard library only):

    python tests/run_tests.py            # everything (~2 minutes)
    python tests/run_tests.py -f         # fast subset (skips the slow convergence tests)
"""
import os
import sys
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

if __name__ == "__main__":
    fast = "-f" in sys.argv
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for mod in ("test_core", "test_models", "test_units_angles", "test_golden", "test_tools", "test_housekeeping") + (() if fast else ("test_grids",)):
        suite.addTests(loader.loadTestsFromName(mod))
    t0 = time.time()
    res = unittest.TextTestRunner(verbosity=1).run(suite)
    print("%d tests in %.0f s%s" % (res.testsRun, time.time() - t0, " (fast subset)" if fast else ""))
    sys.exit(0 if res.wasSuccessful() else 1)

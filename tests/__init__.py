# -*- coding: utf-8 -*-

"""
Python function for MDStudio_SMARTCyp module, run as:
::
    test = module_test_suite()
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(test)
"""

import os
import sys
import unittest
import logging

# Init basic logging
logging.basicConfig(level=logging.INFO)

# Add modules in package to path so we can import them
modulepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
sys.path.insert(0, modulepath)


def module_test_suite():
    """
    Run MDStudio_SMARTCyp module unit tests
    """

    testpath = os.path.join(os.path.dirname(__file__), 'module')
    loader = unittest.TestLoader()
    suite = loader.discover(testpath, pattern='module_*.py')
    return suite

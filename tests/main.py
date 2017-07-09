# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 19:23:27 2017

@author: edgar
"""

import unittest
import sys

loader = unittest.TestLoader()
tests = loader.discover('.')
testRunner = unittest.runner.TextTestRunner()
test_results = testRunner.run(tests)

sys.exit(0)
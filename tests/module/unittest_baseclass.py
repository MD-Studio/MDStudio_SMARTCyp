# -*- coding: utf-8 -*-

"""
file: unittest_baseclass.py

Python 2/3 unittest compatibility class
"""

import unittest
import sys

version = sys.version_info
MAJOR_PY_VERSION = sys.version_info.major
PY_VERSION = '{0}.{1}'.format(version.major, version.minor)

# Unicode test
UNICODE_TYPE = str
if MAJOR_PY_VERSION == 2:
    UNICODE_TYPE = unicode


class UnittestPythonCompatibility(unittest.TestCase):

    def assertItemsEqual(self, expected_seq, actual_seq, msg=None):
        """
        Universal assertItemsEqual method.

        Python 2.x has assertItemsEqual but it is assertCountEqual in 3.x.
        """

        if MAJOR_PY_VERSION == 2:
            return super(UnittestPythonCompatibility, self).assertItemsEqual(expected_seq, actual_seq, msg=msg)
        return super(UnittestPythonCompatibility, self).assertCountEqual(expected_seq, actual_seq, msg=msg)

    def assertNestedObjects(self, ref, tar):
        """
        Test equality of nested objects of different type(s)

        :param ref: reference object to validate against
        :param tar: target object to validate
        """

        self.assertEqual(type(ref), type(tar))

        if isinstance(ref, (list, tuple)):
            self.assertItemsEqual(ref, tar)

        elif isinstance(ref, dict):

            ref_key_val = [(k, v) for k, v in ref.items() if not isinstance(v, (dict, list, tuple))]
            tar_key_val = [(k, v) for k, v in tar.items() if not isinstance(v, (dict, list, tuple))]

            self.assertItemsEqual(ref_key_val, tar_key_val)

            for key, val in ref.items():
                if isinstance(val, (dict, list, tuple)):
                    self.assertNestedObjects(ref[key], tar[key])

    def assertFileEqual(self, ref, tar):
        """
        Test equality between files by comparing all non-empty lines

        :param ref: reference file content to validate against
        :param tar: target file content to validate
        """

        reference_lines = [l.strip() for l in ref.readlines() if len(l.strip())]
        target_lines = [l for l in tar.split('\n') if l]

        self.assertListEqual(reference_lines, target_lines)

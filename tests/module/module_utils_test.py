# -*- coding: utf-8 -*-

"""
file: module_utils_test.py

Unit tests MDStudio_SMARTCyp utility functions
"""

import os
import shutil
import unittest
import platform

from mdstudio_smartcyp import __package_path__
from mdstudio_smartcyp.utils import prepare_work_dir

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
PLANTS_EXEC = os.path.join(__package_path__, 'bin/plants_{0}'.format(platform.system().lower()))


class UtilsTest(unittest.TestCase):

    tempdirs = []

    def tearDown(self):
        """
        tearDown method called after each unittest to cleanup
        the temporary directories
        """

        for tmpdir in self.tempdirs:
            if os.path.exists(tmpdir):
                shutil.rmtree(tmpdir)

    def test_prepare_work_dir_def(self):
        """
        Test default behaviour, create unique tmp dir in system temp.
        """

        path = prepare_work_dir()
        self.tempdirs.append(path)

        self.assertTrue(os.path.isdir(path))
        self.assertTrue(os.access(path, os.W_OK))

    def test_prepare_work_dir_prefix(self):
        """
        Test creation unique temp dir with 'prefix' in basename.
        """

        path = prepare_work_dir(prefix='docking-')
        self.tempdirs.append(path)

        self.assertTrue(os.path.isdir(path))
        self.assertTrue(os.path.basename(path).startswith('docking-'))
        self.assertTrue(os.access(path, os.W_OK))

    def test_prepare_work_dir_suffix(self):
        """
        Test creation unique temp dir with 'suffix' in basename.
        """

        path = prepare_work_dir(suffix='-user')
        self.tempdirs.append(path)

        self.assertTrue(os.path.isdir(path))
        self.assertTrue(os.path.basename(path).endswith('-user'))
        self.assertTrue(os.access(path, os.W_OK))

    def test_prepare_work_dir_preffix_suffix(self):
        """
        Test creation unique temp dir with 'prefix' and 'suffix' in basename.
        """

        path = prepare_work_dir(prefix='docking-', suffix='-user')
        self.tempdirs.append(path)

        self.assertTrue(os.path.isdir(path))
        self.assertTrue(os.path.basename(path).startswith('docking-'))
        self.assertTrue(os.path.basename(path).endswith('-user'))
        self.assertTrue(os.access(path, os.W_OK))

    def test_prepare_work_dir_create(self):
        """
        Test creation of unique temp dir in custom base dir
        """

        to_create = os.path.join(FILEPATH, 'tmp_user_dir')

        path = prepare_work_dir(path=to_create)
        self.tempdirs.append(to_create)

        self.assertEqual(os.path.dirname(path), to_create)
        self.assertTrue(os.path.isdir(path))
        self.assertTrue(os.access(path, os.W_OK))

    def test_prepare_work_dir_nocreate(self):
        """
        Test creation of unique temp dir in custom base dir, raise IOError
        when base dir not present
        """

        to_create = os.path.join(FILEPATH, 'tmp_user_dir')

        self.assertRaises(IOError, prepare_work_dir, path=to_create, create=False)

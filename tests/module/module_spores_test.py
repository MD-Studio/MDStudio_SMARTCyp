# -*- coding: utf-8 -*-

"""
file: module_plants_test.py

Unit tests for PlantsDocking methods
"""

import os
import glob
import shutil
import unittest
import platform

from mdstudio_smartcyp import __package_path__
from mdstudio_smartcyp.spores_run import SporesRunner

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
PLANTS_EXEC = os.path.join(__package_path__, 'bin/spores_{0}'.format(platform.system().lower()))


class SporesTest(unittest.TestCase):

    workdir = None
    ligand_file = os.path.join(FILEPATH, 'ligand.mol2')
    protein_file = os.path.join(FILEPATH, 'protein.mol2')

    @classmethod
    def setUpClass(cls):
        """
        SporesTest class setup

        Read structure files for docking
        """

        with open(cls.protein_file, 'r') as pfile:
            cls.protein = pfile.read()

        with open(cls.ligand_file, 'r') as lfile:
            cls.ligand = lfile.read()

    @classmethod
    def tearDownClass(cls):
        """
        Remove all created result files
        """

        for sporesdir in glob.glob(os.path.join(FILEPATH, 'spores-*/')):
            if os.path.isdir(sporesdir):
                shutil.rmtree(sporesdir)

    def test_spores_faultyexec(self):
        """
        Docking is unable to start if the PLANTS executable is not found
        """

        spores = SporesRunner(base_work_dir=FILEPATH,
                              exec_path='/Users/_dummy_user/smartcyp/tests/spores')

        self.assertFalse(spores.run(self.protein, self.ligand))
        self.assertIsNone(spores.workdir)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_spores(self):

        pass

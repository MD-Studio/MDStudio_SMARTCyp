# -*- coding: utf-8 -*-

"""
file: module_combined_test.py

Unit tests for combined SMARTCyp and PLANTS Cyp SOM prediction
"""

import os
import glob
import shutil
import unittest
import platform

from mdstudio_smartcyp import __package_path__
from mdstudio_smartcyp.combined_prediction import CombinedPrediction
from mdstudio_smartcyp.plants_run import PlantsDocking
from tests.module.unittest_baseclass import UnittestPythonCompatibility

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
PLANTS_EXEC = os.path.join(__package_path__, 'bin/plants_{0}'.format(platform.system().lower()))


class CombinedPredictionTest(UnittestPythonCompatibility):

    workdir = None
    prediction = None
    ligand_file = os.path.join(FILEPATH, 'ligand.mol2')

    @classmethod
    def setUpClass(cls):
        """
        CombinedPredictionTest class setup

        Read structure files for docking
        """

        with open(cls.ligand_file, 'r') as lfile:
            cls.ligand = lfile.read()

    @classmethod
    def tearDownClass(cls):
        """
        Remove all created docking result dirs
        """

        for dockdir in glob.glob(os.path.join(FILEPATH, 'docking-*/')):
            if os.path.isdir(dockdir):
                shutil.rmtree(dockdir)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_a_som_prediction(self):
        """
        A working plants docking, the number of generated docking poses should
        equal the 'cluster_structures' setting and the results dictionary
        """

        som = CombinedPrediction(base_work_dir=FILEPATH)
        prediction = som.run(self.ligand)

        self.assertTrue(all(['Docking' in i for i in  prediction.values()]))
        self.assertTrue(all(['SMARTCyp' in i for i in  prediction.values()]))

        self.__class__.prediction = prediction

    def test_b_get_structures(self):

        if self.__class__.prediction:

            poses = []
            for val in self.__class__.prediction.values():
                if val['Docking'] == 1.0:
                    poses = [i for i in val.keys() if i not in ('Docking', 'SMARTCyp')]
                    break

            self.assertTrue(len(poses) > 0)
            docking = PlantsDocking(base_work_dir=FILEPATH)
            with open('OUT.pdb', 'w') as outpdb:
                outpdb.write(docking.get_structures(poses, output_format='pdb', include_protein=True))

    def test_c_get_results(self):

        if self.__class__.prediction:

            poses = []
            for val in self.__class__.prediction.values():
                if val['Docking'] == 1.0:
                    poses = [i for i in val.keys() if i not in ('Docking', 'SMARTCyp')]
                    break

            docking = PlantsDocking(base_work_dir=FILEPATH)
            stats = docking.get_results(poses)

            clusters = [i['CLUSTER'] for i in stats.values()]
            self.assertTrue(len(set(clusters)) >= 1)

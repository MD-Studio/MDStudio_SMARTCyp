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
from mdstudio_smartcyp.plants_run import PlantsDocking
from mdstudio_smartcyp.plants_run import MDStudioException
from tests.module.unittest_baseclass import UnittestPythonCompatibility

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
PLANTS_EXEC = os.path.join(__package_path__, 'bin/plants_{0}'.format(platform.system().lower()))


class PlantsDockingTest(UnittestPythonCompatibility):

    workdir = None
    ligand_file = os.path.join(FILEPATH, 'ligand.mol2')
    protein_file = os.path.join(FILEPATH, 'protein.mol2')

    @classmethod
    def setUpClass(cls):
        """
        PlantsDockingTest class setup

        Read structure files for docking
        """

        with open(cls.protein_file, 'r') as pfile:
            cls.protein = pfile.read()

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

    def run_plants(self):

        plants = PlantsDocking(base_work_dir=FILEPATH, bindingsite_center=[-0.989, 3.261, 0.826])

        did_run_successfully = True
        if self.__class__.workdir is None:
            did_run_successfully = plants.run(self.protein, self.ligand)
            self.__class__.workdir = plants.workdir
        else:
            plants.workdir = self.__class__.workdir

        return did_run_successfully, plants

    def test_plants_faultyexec(self):
        """
        Docking is unable to start if the PLANTS executable is not found
        """

        plants = PlantsDocking(workdir=FILEPATH,
                               exec_path='/Users/_dummy_user/smartcyp/tests/plants',
                               bindingsite_center=[-0.989, 3.261, 0.826])

        self.assertRaises(MDStudioException, plants.run, self.ligand, self.protein)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_plants_docking(self):
        """
        A working plants docking, the number of generated docking poses should
        equal the 'cluster_structures' setting and the results dictionary
        """

        did_run_successfully, plants = self.run_plants()
        self.assertTrue(did_run_successfully)

        outputfiles = glob.glob('{0}/_entry_00001_conf_*.mol2'.format(plants.workdir))
        self.assertEqual(len(outputfiles), plants.config['cluster_structures'])
        self.assertEqual(len(outputfiles), len(plants.get_results()))

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_plants_docking_get_structures(self):
        """
        Get specific poses based on results file
        """

        did_run_successfully, plants = self.run_plants()
        self.assertTrue(did_run_successfully)

        results = plants.get_results()
        paths = [v['PATH'] for v in list(results.values())[0:5]]
        mol = plants.get_structures(paths)

        self.assertIsNotNone(mol)

        molcount = 0
        for line in mol.split('\n'):
            if line.startswith('@<TRIPOS>MOLECULE'):
                molcount += 1
        self.assertEqual(molcount, 5)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_plants_docking_get_structures_not_exist(self):
        """
        Get specific poses for docking run that no longer exists
        """

        plants = PlantsDocking(base_work_dir=FILEPATH, bindingsite_center=[-0.989, 3.261, 0.826])

        paths = ['not_exist/_entry_00001_conf_{0}.mol2'.format(i) for i in range(10)]
        self.assertRaises(MDStudioException, plants.get_structures, paths)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_plants_docking_get_results_not_exist(self):
        """
        Get specific results for docking run that no longer exists
        """

        plants = PlantsDocking(base_work_dir=FILEPATH, bindingsite_center=[-0.989, 3.261, 0.826])

        paths = ['not_exist/_entry_00001_conf_{0}.mol2'.format(i) for i in range(10)]
        self.assertRaises(MDStudioException, plants.get_results, paths)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_plants_docking_get_results(self):
        """
        Get results for subset of structures, will redo clustering
        """

        did_run_successfully, plants = self.run_plants()
        self.assertTrue(did_run_successfully)

        # Get all results
        results = plants.get_results()
        self.assertEqual(len(results.values()), 50)

        # Get results subset
        paths = [v['PATH'] for v in list(results.values())[0:40]]

        results = plants.get_results(paths)
        self.assertEqual(len(results), 40)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_plants_docking_redo_clustering(self):
        """
        Redo clustering
        """

        did_run_successfully, plants = self.run_plants()
        self.assertTrue(did_run_successfully)

        results = plants.get_results()
        clust_orig = [b['CLUSTER'] for b in results.values()]

        plants.update(threshold=4.0)
        results = plants.get_results()
        clust_new = [b['CLUSTER'] for b in results.values()]

        self.assertNotEqual(clust_orig, clust_new)

    @unittest.skipIf(not os.path.exists(PLANTS_EXEC), 'This test requires proprietary software')
    def test_plants_docking_wrong_structures(self):
        """
        Faulty docking with a ligand uploaded as protein
        """

        plants = PlantsDocking(base_work_dir=FILEPATH, bindingsite_center=[-0.989, 3.261, 0.826])
        self.assertRaises(MDStudioException, plants.run, self.ligand, self.protein)

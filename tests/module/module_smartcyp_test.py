# -*- coding: utf-8 -*-

"""
file: module_smartcyp_test.py

Unit tests for SmartCypRunner methods
"""

import os
import json
import unittest

from mdstudio_smartcyp.smartcyp_run import SmartCypRunner

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))


class SmartCypRunnerTests(unittest.TestCase):

    def setUp(self):
        """
        Init a SmartCypRunner class
        """

        self.scr = SmartCypRunner()
        self.assertTrue(os.path.exists(self.scr.tempdir))

    def tearDown(self):
        """
        Check if temporary dir is cleaned
        """

        self.assertFalse(os.path.exists(self.scr.tempdir))

    def test_smartcyp_smiles(self):
        """
        Test default run of SMARTCyp on SMILES string
        """

        json_result = self.scr.run(mol='O=Cc1ccc(s1)c2cccnc2', is_smiles=True)

        reference = open(os.path.join(FILEPATH, 'result.json'), 'r')
        json_reference = json.load(reference)

        self.assertDictEqual(json_result['result'], json_reference)

    def test_smartcyp_mol2(self):
        """
        Test default run of SMARTCyp on a mol2 file
        """

        with open(os.path.join(FILEPATH, 'ligand.mol2'), 'r') as infile:
            molfile = infile.read()

        json_result = self.scr.run(mol=molfile)

        reference = open(os.path.join(FILEPATH, 'result.json'), 'r')
        json_reference = json.load(reference)

        self.assertDictEqual(json_result['result'], json_reference)

    def test_smartcyp_sdf(self):
        """
        Test default run of SMARTCyp on a sdf file
        """

        with open(os.path.join(FILEPATH, 'ligand.sdf'), 'r') as infile:
            molfile = infile.read()

        json_result = self.scr.run(mol=molfile)

        reference = open(os.path.join(FILEPATH, 'result.json'), 'r')
        json_reference = json.load(reference)

        self.assertDictEqual(json_result['result'], json_reference)

    def test_smartcyp_noempcorr(self):
        """
        Test default run of SMARTCyp on SMILES string not use the empirical
        N-oxidation correction
        """

        json_result = self.scr.run(mol='O=Cc1ccc(s1)c2cccnc2', is_smiles=True, noempcorr=True)

        reference = open(os.path.join(FILEPATH, 'result.json'), 'r')
        json_reference = json.load(reference)

        self.assertDictEqual(json_result['result'], json_reference)

    def test_smartcyp_csv(self):
        """
        Test SMARTCyp CSV export
        """

        csv_result = self.scr.run(mol='O=Cc1ccc(s1)c2cccnc2', is_smiles=True, output_format='csv')

        with open(os.path.join(FILEPATH, 'result.csv'), 'r') as reference:
            self.assertEqual(reference.read(), csv_result['result'])

    def test_smartcyp_html(self):
        """
        Test SMARTCyp HTML export
        """

        html_result = self.scr.run(mol='O=Cc1ccc(s1)c2cccnc2', is_smiles=True, output_format='html')

    def test_smartcyp_png(self):
        """
        Test SMARTCyp PNG image export
        """

        png_result = self.scr.run(mol='O=Cc1ccc(s1)c2cccnc2', is_smiles=True, output_png=True)

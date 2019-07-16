# -*- coding: utf-8 -*-

"""
file: module_spores_rest_test.py

Unit tests for the REST interface to SporesRunner methods
"""

import os
import unittest
import requests

from mdstudio_smartcyp.spores_run import spores_version_info

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
URL = 'http://localhost:8081'


def import_mol2(mol):

    lines = []
    read = False
    for line in mol.split('\n'):

        if line.startswith('@<TRIPOS>ATOM'):
            read = True
            continue

        if line.startswith('@<TRIPOS>BOND'):
            read = False
            break

        if read:
            line = line.strip()
            lines.append(line.split())

    return lines


def test_localhost_connection():

    try:
        requests.get(URL)
    except requests.exceptions.ConnectionError:
        return False

    return True


class SporesRestTest(unittest.TestCase):

    ligand_sdf = open(os.path.join(FILEPATH, 'ligand.sdf'), 'rb')
    ligand_mol2 = open(os.path.join(FILEPATH, 'ligand.mol2'), 'rb')

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_spores_info_get(self):
        """
        Test default spores_info get response
        """

        response = requests.get('{0}/spores_info'.format(URL))

        rest_response = response.json()
        func_response = spores_version_info()

        self.assertDictEqual(rest_response, func_response)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_spores_info_post(self):
        """
        Post to get endpoint not allowed
        """

        response = requests.post('{0}/spores_info'.format(URL))

        rest_response = response.json()
        self.assertEqual(rest_response.get('status', 0), 405)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_spores_sdf(self):
        """
        Test default smartcyp post response, JSON results
        """

        data = {'input_format': 'mol2', 'spores_mode': 'complete'}
        files = {'mol': self.ligand_mol2}
        response = requests.post('{0}/spores'.format(URL), files=files, data=data)

        rest_response = response.json()
        reffile = open(os.path.join(FILEPATH, 'ligand.mol2'))
        self.assertListEqual(import_mol2(reffile.read()), import_mol2(rest_response['content']))

    # @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    # def test_smartcyp_html(self):
    #     """
    #     Test default smartcyp post response, HTML result
    #     """
    #
    #     data = {'smiles': 'O=Cc1ccc(s1)c2cccnc2', 'output_format': 'html'}
    #     response = requests.post('{0}/smartcyp'.format(URL), data=data)
    #
    #     rest_response = response.json()
    #     self.assertTrue('result' in rest_response)
    #
    #     parser = HTMLParser()
    #     self.assertIsNone(parser.feed(rest_response['result']))
    #
    # @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    # def test_smartcyp_csv(self):
    #     """
    #     Test default smartcyp post response, CSV result
    #     """
    #
    #     data = {'smiles': 'O=Cc1ccc(s1)c2cccnc2', 'output_format': 'csv'}
    #     response = requests.post('{0}/smartcyp'.format(URL), data=data)
    #
    #     rest_response = response.json()
    #     with open(os.path.join(FILEPATH, 'result.csv'), 'r') as reference:
    #         self.assertEqual(reference.read(), rest_response['result'])
    #
    # @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    # def test_smartcyp_fileupload(self):
    #     """
    #     Test default smartcyp post response, file upload
    #     """
    #
    #     files = {'mol': self.ligand_file}
    #     response = requests.post('{0}/smartcyp'.format(URL), files=files)
    #
    #     rest_response = response.json()
    #     self.assertTrue('result' in rest_response)
    #     self.assertEqual(len(rest_response['result']), 13)
    #
    # @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    # def test_smartcyp_missing_required(self):
    #     """
    #     Test default smartcyp post response, file upload
    #     """
    #
    #     response = requests.post('{0}/smartcyp'.format(URL))
    #
    #     rest_response = response.json()
    #     self.assertEqual(response.status_code, 401)
    #     self.assertFalse('result' in rest_response)

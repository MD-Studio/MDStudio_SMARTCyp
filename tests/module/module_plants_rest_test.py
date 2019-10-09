# -*- coding: utf-8 -*-

"""
file: module_plants_rest_test.py

Unit tests for the REST interface to the PlantsDocking methods
"""

import os
import unittest
import requests
import random

from mdstudio_smartcyp.plants_run import plants_version_info
from tests.module.unittest_baseclass import UnittestPythonCompatibility

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
URL = 'http://localhost:8081'
DOCKRESULTS = {u'ACC', u'PLPtotal', u'SCORE_NORM_WEIGHT', u'DON', u'SCORE_RB_PEN', u'TOTAL_SCORE',
               u'PLPpartsteric', u'SCORE_RB_PEN_NORM_CRT_HEVATOMS', u'PLPpartmetal', u'PATH', u'LIG_NUM_CLASH',
               u'LIG_NUM_NO_CONTACT', u'SCORE_NORM_CRT_WEIGHT', u'PLPparthbond', u'CLUSTER', u'PLPpartburpolar',
               u'TRIPOS_TORS', u'SCORE_NORM_HEVATOMS', u'SCORE_NORM_CONTACT', u'UNUSED_ACC', u'MEAN', u'CHEMpartmetal',
               u'CHEMparthbondCHO', u'CHEMparthbond', u'ATOMS_OUTSIDE_BINDINGSITE', u'UNUSED_DON', u'PLPpartrepulsive',
               u'CHEMPLP_CLASH2', u'SCORE_NORM_CRT_HEVATOMS', u'LIG_NUM_CONTACT'}


def test_localhost_connection():

    try:
        requests.get(URL)
    except requests.exceptions.ConnectionError:
        return False

    return True


def count_poses_mulimol(mol):

    pose_count = 0
    for line in mol.split('\n'):
        if line.startswith('@<TRIPOS>ATOM'):
            pose_count += 1

    return pose_count


class PlantsRestTest(UnittestPythonCompatibility):

    ligand_file = open(os.path.join(FILEPATH, 'ligand.mol2'), 'rb')
    protein_file = open(os.path.join(FILEPATH, 'protein.mol2'), 'rb')

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_info_get(self):
        """
        Test default docking_info get response.
        Remove 'citation' as it requires a change in encoding
        """

        response = requests.get('{0}/plants_docking_info'.format(URL))

        rest_response = response.json()
        if 'citation' in rest_response:
            del rest_response['citation']

        func_response = plants_version_info()
        if 'citation' in func_response:
            del func_response['citation']

        self.assertDictEqual(rest_response, func_response)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_info_post(self):
        """
        Post to get endpoint not allowed
        """

        response = requests.post('{0}/plants_docking_info'.format(URL))

        rest_response = response.json()
        self.assertEqual(rest_response.get('status', 0), 405)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_info_get_data(self):
        """
        Docking_info does not accept arguments but will still return results
        if any given. Remove 'citation' as it requires a change in encoding
        """

        response = requests.get('{0}/plants_docking_info'.format(URL), params={'input': 10})

        rest_response = response.json()
        if 'citation' in rest_response:
            del rest_response['citation']

        func_response = plants_version_info()
        if 'citation' in func_response:
            del func_response['citation']

        self.assertDictEqual(rest_response, func_response)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking(self):
        """
        Test default docking post response
        """

        files = {'ligand_file': self.ligand_file, 'protein_file': self.protein_file}
        data = {'bindingsite_center': [-0.989, 3.261, 0.826]}
        response = requests.post('{0}/plants_docking'.format(URL), files=files, data=data)

        rest_response = response.json()
        self.assertEqual(len(rest_response), 50)
        self.assertTrue(all([set(v.keys()) == DOCKRESULTS for v in rest_response.values()]))

        selection = random.sample(rest_response.keys(), 25)

        # Get docking results for selection
        response = requests.post('{0}/plants_docking_statistics'.format(URL),
                                 data={'paths': [rest_response[n]['PATH'] for n in selection]})

        rest_stats_response = response.json()

        self.assertListEqual(sorted(rest_stats_response.keys()), sorted(selection))

        # Get docking solutions for selection
        response = requests.post('{0}/plants_docking_structures'.format(URL),
                                 data={'paths': [rest_response[n]['PATH'] for n in selection]})

        rest_poses_response = response.text
        self.assertEqual(count_poses_mulimol(rest_poses_response), len(selection))

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_get_statistics_noresults(self):
        """
        Test call to docking_get_statistics while results no (longer) available
        """

        dummy_selection = ['docking-noexists/_entry_00001_conf_{0}.mol2'.format(n) for n in range(10)]

        # Get docking results for selection
        response = requests.post('{0}/plants_docking_statistics'.format(URL), data={'paths': dummy_selection})

        self.assertEqual(response.status_code, 401)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_wrong_structures(self):
        """
        Faulty docking with a ligand uploaded as protein
        """

        files = {'ligand_file': self.ligand_file, 'protein_file': self.ligand_file}
        data = {'bindingsite_center': [-0.989, 3.261, 0.826]}
        response = requests.post('{0}/plants_docking'.format(URL), files=files, data=data)

        self.assertEqual(response.status_code, 401)

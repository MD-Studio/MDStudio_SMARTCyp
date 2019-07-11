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

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
URL = 'http://localhost:8081'
DOCKRESULTS = ['ACC', 'ATOMS_OUTSIDE_BINDINGSITE', 'CHEMPLP_CLASH2', 'CHEMparthbond', 'CHEMparthbondCHO',
               'CHEMpartmetal', 'CLUSTER', 'DON', 'LIG_NUM_CLASH', 'LIG_NUM_CONTACT', 'LIG_NUM_NO_CONTACT', 'MEAN',
               'PATH', 'PLPpartburpolar', 'PLPparthbond', 'PLPpartmetal', 'PLPpartrepulsive', 'PLPpartsteric',
               'PLPtotal', 'SCORE_NORM_CONTACT', 'SCORE_NORM_CRT_HEVATOMS', 'SCORE_NORM_CRT_WEIGHT',
               'SCORE_NORM_HEVATOMS', 'SCORE_NORM_WEIGHT', 'SCORE_RB_PEN', 'SCORE_RB_PEN_NORM_CRT_HEVATOMS',
               'TOTAL_SCORE', 'TRIPOS_TORS', 'UNUSED_ACC', 'UNUSED_DON']


def test_localhost_connection():

    try:
        result = requests.get(URL)
    except requests.exceptions.ConnectionError:
        return False

    return True


def count_poses_mulimol(mol):

    pose_count = 0
    for line in mol.split('\n'):
        if line.startswith('@<TRIPOS>ATOM'):
            pose_count += 1

    return pose_count


class PlantsRestTest(unittest.TestCase):

    ligand_file = open(os.path.join(FILEPATH, 'ligand.mol2'), 'rb')
    protein_file = open(os.path.join(FILEPATH, 'protein.mol2'), 'rb')

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_info_get(self):
        """
        Test default docking_info get response
        """

        response = requests.get('{0}/docking_info'.format(URL))

        rest_response = response.json()
        func_response = plants_version_info()

        self.assertDictEqual(rest_response, func_response)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_info_post(self):
        """
        Post to get endpoint not allowed
        """

        response = requests.post('{0}/docking_info'.format(URL))

        rest_response = response.json()
        self.assertEqual(rest_response.get('status', 0), 405)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_info_get_data(self):
        """
        Docking_info does not accept arguments but will still return results
        if any given.
        """

        response = requests.get('{0}/docking_info'.format(URL), params={'input': 10})

        rest_response = response.json()
        func_response = plants_version_info()

        self.assertDictEqual(rest_response, func_response)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking(self):
        """
        Test default docking post response
        """

        files = {'ligand_file': self.ligand_file, 'protein_file': self.protein_file}
        data = {'bindingsite_center': [7.79934, 9.49666, 3.39229]}
        response = requests.post('{0}/docking'.format(URL), files=files, data=data)

        rest_response = response.json()
        self.assertEqual(len(rest_response), 50)
        self.assertTrue(all([list(v.keys()) == DOCKRESULTS for v in rest_response.values()]))

        selection = random.sample(rest_response.keys(), 25)

        # Get docking results for selection
        response = requests.post('{0}/docking_get_statistics'.format(URL),
                                 data={'paths': [rest_response[n]['PATH'] for n in selection]})

        rest_stats_response = response.json()

        self.assertListEqual(sorted(rest_stats_response.keys()), sorted(selection))

        # Get docking solutions for selection
        response = requests.post('{0}/docking_get_structures'.format(URL),
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
        response = requests.post('{0}/docking_get_statistics'.format(URL), data={'paths': dummy_selection})

        self.assertEqual(response.status_code, 401)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_docking_wrong_structures(self):
        """
        Faulty docking with a ligand uploaded as protein
        """

        files = {'ligand_file': self.ligand_file, 'protein_file': self.ligand_file}
        data = {'bindingsite_center': [7.79934, 9.49666, 3.39229]}
        response = requests.post('{0}/docking'.format(URL), files=files, data=data)

        self.assertEqual(response.status_code, 401)



# -*- coding: utf-8 -*-

"""
file: module_smartcyo_rest_test.py

Unit tests for the REST interface to the SmartCypRunner methods
"""

import os
import unittest
import requests
import json

from mdstudio_smartcyp.smartcyp_run import smartcyp_version_info
from . import MAJOR_PY_VERSION

if MAJOR_PY_VERSION == 2:
    from HTMLParser import HTMLParser
else:
    from html.parser import HTMLParser

FILEPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../files/'))
URL = 'http://localhost:8081'


def test_localhost_connection():

    try:
        result = requests.get(URL)
    except requests.exceptions.ConnectionError:
        return False

    return True


class SmartcypRestTest(unittest.TestCase):

    ligand_file = open(os.path.join(FILEPATH, 'ligand.sdf'), 'rb')

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_smartcyp_info_get(self):
        """
        Test default smartcyp_info get response
        """

        response = requests.get('{0}/smartcyp_info'.format(URL))

        rest_response = response.json()
        func_response = smartcyp_version_info()

        self.assertDictEqual(rest_response, func_response)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_smartcyp_info_post(self):
        """
        Post to get endpoint not allowed
        """

        response = requests.post('{0}/smartcyp_info'.format(URL))

        rest_response = response.json()
        self.assertEqual(rest_response.get('status', 0), 405)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_smartcyp_json(self):
        """
        Test default smartcyp post response, JSON results
        """

        data = {'smiles': 'O=Cc1ccc(s1)c2cccnc2'}
        response = requests.post('{0}/smartcyp'.format(URL), data=data)

        reference = open(os.path.join(FILEPATH, 'result.json'), 'r')
        json_reference = json.load(reference)

        rest_response = response.json()
        self.assertDictEqual(rest_response['result'], json_reference)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_smartcyp_html(self):
        """
        Test default smartcyp post response, HTML result
        """

        data = {'smiles': 'O=Cc1ccc(s1)c2cccnc2', 'output_format': 'html'}
        response = requests.post('{0}/smartcyp'.format(URL), data=data)

        rest_response = response.json()
        self.assertTrue('result' in rest_response)

        parser = HTMLParser()
        self.assertIsNone(parser.feed(rest_response['result']))

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_smartcyp_csv(self):
        """
        Test default smartcyp post response, CSV result
        """

        data = {'smiles': 'O=Cc1ccc(s1)c2cccnc2', 'output_format': 'csv'}
        response = requests.post('{0}/smartcyp'.format(URL), data=data)

        rest_response = response.json()
        with open(os.path.join(FILEPATH, 'result.csv'), 'r') as reference:
            self.assertEqual(reference.read(), rest_response['result'])

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_smartcyp_fileupload(self):
        """
        Test default smartcyp post response, file upload
        """

        files = {'mol': self.ligand_file}
        response = requests.post('{0}/smartcyp'.format(URL), files=files)

        rest_response = response.json()
        self.assertTrue('result' in rest_response)
        self.assertEqual(len(rest_response['result']), 13)

    @unittest.skipIf(not test_localhost_connection(), 'MDStudio_SMARTCyp REST service not running on: {0}'.format(URL))
    def test_smartcyp_missing_required(self):
        """
        Test default smartcyp post response, file upload
        """

        response = requests.post('{0}/smartcyp'.format(URL))

        rest_response = response.json()
        self.assertEqual(response.status_code, 401)
        self.assertFalse('result' in rest_response)

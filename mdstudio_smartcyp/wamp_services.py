# -*- coding: utf-8 -*-

"""
file: wamp_services.py

WAMP service methods the module exposes.
"""

import os
import re

from autobahn.wamp import RegisterOptions

from mdstudio.api.endpoint import endpoint
from mdstudio.component.session import ComponentSession

from mdstudio_smartcyp import __smartcyp_version__, __smartcyp_citation__, __supported_models__
from mdstudio_smartcyp.smartcyp_run import SmartCypRunner

smiles_regex = re.compile('^([^J][A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+)$')


def mol_validate_file_object(path_file):
    """
    Validate a MDStudio path_file object

    - Check if 'content' is a InChI or SMILES string and set extension
      (mol_format)
    - If no 'content' check if path exists

    :param path_file: path_file object
    :type path_file:  :py:dict

    :return:          validated path_file object
    :rtype:           :py:dict
    """

    content = path_file['content']
    if content is not None:

        # SMILES and InChI are single line strings
        if len(content.split('\n')) == 1:

            # Test for InChI type
            if content.startswith('InChI='):
                path_file['extension'] = 'inchi'

            # Test for SMILES
            if smiles_regex.match(path_file['content']):
                path_file['extension'] = 'smi'

    elif path_file['path'] is not None and os.path.exists(path_file['path']):

        with open(path_file['path']) as pf:
            path_file['content'] = pf.read()

    return path_file


class SmartCypWampApi(ComponentSession):
    """
    WAMP SMARTCyp endpoint methods
    """

    def authorize_request(self, uri, claims):
        return True

    @endpoint('info', 'info_request', 'info_response',
              options=RegisterOptions(invoke=u'roundrobin'))
    def smartcyp_info(self, request, claims):
        """
        Returns an informative summary of the supported SMARTCyp version
        and configuration.
        """

        info_dict = {'version': __smartcyp_version__,
                     'models': __supported_models__,
                     'citation': __smartcyp_citation__}

        return info_dict

    @endpoint('predict', 'predict_request', 'predict_response',
              options=RegisterOptions(invoke=u'roundrobin'))
    def smartcyp_prediction(self, request, claims):
        """
        Run a SMARTCyp prediction for a molecule

        :param request:
        :param claims:
        :return:
        """

        mol = mol_validate_file_object(request['mol'])

        smartcyp = SmartCypRunner(log=self.log)
        result_dict = smartcyp.run(mol['content'],
                                   is_smiles=mol['extension'] == 'smi',
                                   output_format=request['output_format'],
                                   noempcorr=request['noempcorr'])

        return result_dict

# -*- coding: utf-8 -*-

"""
file: wamp_services.py

WAMP service methods the module exposes.
"""

from autobahn.wamp import RegisterOptions

from mdstudio.api.endpoint import endpoint
from mdstudio.component.session import ComponentSession

from mdstudio_smartcyp.smartcyp_run import SmartCypRunner, smartcyp_version_info, mol_validate_file_object


class SmartCypWampApi(ComponentSession):
    """
    WAMP SMARTCyp endpoint methods
    """

    def authorize_request(self, uri, claims):

        return True

    @endpoint('info', 'info_request', 'info_response', options=RegisterOptions(invoke=u'roundrobin'))
    def smartcyp_info(self, request, claims):
        """
        Returns an informative summary of the supported SMARTCyp version
        and configuration.
        """

        return smartcyp_version_info()

    @endpoint('predict', 'predict_request', 'predict_response', options=RegisterOptions(invoke=u'roundrobin'))
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

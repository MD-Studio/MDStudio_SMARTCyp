# -*- coding: utf-8 -*-

"""
Unittests for the MDStudio SMARTCyp WAMP API

.. note::

   You are required to first start the core MDStudio infrastructure and the
   MDStudio_SMARTCyp package as WAMP microservice in order to run these tests
"""

import os

from mdstudio.deferred.chainable import chainable
from mdstudio.component.session import ComponentSession
from mdstudio.runner import main

root = os.path.abspath(os.path.join(os.path.dirname(__file__), 'files/'))


def create_path_file_obj(path):
    """
    Encode the input files
    """
    extension = os.path.splitext(path)[1]
    with open(path, 'r') as f:
        content = f.read()

    return {
        u'path': path, u'content': content,
        u'extension': extension}


class MDStudioSMARTCypWAMPTests(ComponentSession):

    def authorize_request(self, uri, claims):
        return True

    def assertEqual(self, first, second):

        if first == second:
            print('ok')
        else:
            raise AssertionError('Instances not the same')

    def assertIsInstance(self, first, second):

        if isinstance(first, second):
            print('ok')
        else:
            raise AssertionError('Not instance of {0} ({1})'.format(repr(second), type(first)))

    def assertItemsEqual(self, first, second):

        if set(first) == set(second):
            print('ok')
        else:
            raise AssertionError('Items not equal')

    @chainable
    def on_run(self):

        print('Test SMARTCyp get info endpoint')
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.smartcyp_info', {})
        self.assertEqual(result['status'], 'completed')
        self.assertIsInstance(result['info'], dict)
        self.assertItemsEqual(result['info'].keys(), (u'version', u'models', u'citation'))

        print('Test SMARTCyp prediction from SMILES')
        request = {'ligand_file': {'content': 'O=Cc1ccc(s1)c2cccnc2', 'extension': 'smi', 'path': None}}
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.smartcyp', request)
        self.assertEqual(result['status'], 'completed')
        self.assertIsInstance(result['result'], dict)
        self.assertItemsEqual(result['result'].keys(), (u'C.6', u'N.12', u'C.4', u'C.5', u'C.2', u'C.3', u'C.8', u'C.9',
                                             u'O.1', u'S.7', u'C.10', u'C.11', u'C.13'))

        print('Test SMARTCyp prediction for ligand file')
        request = {'ligand_file': create_path_file_obj(os.path.join(root, "ligand.mol2"))}
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.smartcyp', request)
        self.assertEqual(result['status'], 'completed')
        self.assertIsInstance(result['result'], dict)
        self.assertItemsEqual(result['result'].keys(), (u'C.6', u'N.12', u'C.4', u'C.5', u'C.2', u'C.3', u'C.8', u'C.9',
                                                        u'O.1', u'S.7', u'C.10', u'C.11', u'C.13'))

        print('Test PLANTS docking info endpoint')
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.docking_info', {})
        self.assertEqual(result['status'], 'completed')
        self.assertIsInstance(result['info'], dict)
        self.assertItemsEqual(result['info'].keys(), (u'version', u'default_settings', u'citation'))

        print('Test PLANTS docking with unsupported SDF file format')
        request = {'bindingsite_center': [-0.989, 3.261, 0.826],
                   'ligand_file': create_path_file_obj(os.path.join(root, "ligand.sdf")),
                   'protein_file': create_path_file_obj(os.path.join(root, "protein.mol2"))}
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.docking', request)
        self.assertEqual(result['status'], 'failed')

        print('Test PLANTS docking defaults')
        request = {'bindingsite_center': [-0.989, 3.261, 0.826],
                   'ligand_file': create_path_file_obj(os.path.join(root, "ligand.mol2")),
                   'protein_file': create_path_file_obj(os.path.join(root, "protein.mol2"))}
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.docking', request)
        self.assertEqual(result['status'], 'completed')
        self.assertIsInstance(result['result'], dict)
        self.assertEqual(len(result['result']), 50)

        print('Test PLANTS structure retrieval')
        means = [n[u'PATH'] for n in result['result'].values() if n['MEAN'] == True]
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.docking_structures',
                                 {'paths': means})
        self.assertEqual(len(result['result']), 1)
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.docking_structures',
                                 {'paths': means, 'create_ensemble': False})
        self.assertEqual(len(result['result']), len(means))

        print('Test SPORES get info endpoint')
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.spores_info', {})
        self.assertEqual(result['status'], 'completed')
        self.assertIsInstance(result['info'], dict)
        self.assertItemsEqual(result['info'].keys(), (u'version', u'default_settings', u'citation'))

        print('Test SPORES run')
        request = {'mol': create_path_file_obj(os.path.join(root, "ligand.mol2")), 'input_format': 'mol2'}
        result = yield self.call(u'mdgroup.mdstudio_smartcyp.endpoint.spores', request)
        self.assertEqual(result['status'], 'completed')
        self.assertIsInstance(result['result'], dict)
        self.assertEqual(result['result']['content'].startswith('@<TRIPOS>'), True)


if __name__ == "__main__":

    main(MDStudioSMARTCypWAMPTests)

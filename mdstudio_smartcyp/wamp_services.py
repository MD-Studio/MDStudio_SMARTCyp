# -*- coding: utf-8 -*-

"""
file: wamp_services.py

WAMP service methods the module exposes.
"""

import os

from autobahn.wamp import RegisterOptions
from mdstudio.api.endpoint import endpoint
from mdstudio.component.session import ComponentSession

from mdstudio_smartcyp.smartcyp_run import SmartCypRunner, smartcyp_version_info
from mdstudio_smartcyp.plants_run import PlantsDocking, plants_version_info
from mdstudio_smartcyp.spores_run import SporesRunner
from mdstudio_smartcyp.utils import mol_validate_file_object


def encoder(file_path):
    """
    Encode the content of `file_path` into a simple dict.
    """
    extension = os.path.splitext(file_path)[1]

    return {"path": file_path, "extension": extension.lstrip('.'),
            "content": None, "encoding": "utf8"}


def encode_file(d):
    """ Serialize the file containing the molecular representation"""
    d['path'] = encoder(d['path'])
    return d


class SmartCypWampApi(ComponentSession):
    """
    WAMP SMARTCyp endpoint methods
    """

    def authorize_request(self, uri, claims):

        return True

    @endpoint('smartcyp_info', 'info_request', 'info_response', options=RegisterOptions(invoke=u'roundrobin'))
    def smartcyp_info(self, request, claims):
        """
        Returns an informative summary of the supported SMARTCyp version
        and configuration.
        """

        return smartcyp_version_info()

    @endpoint('smartcyp', 'smartcyp_request', 'smartcyp_response', options=RegisterOptions(invoke=u'roundrobin'))
    def smartcyp_prediction(self, request, claims):
        """
        Run a SMARTCyp prediction for a molecule

        :param request:
        :param claims:
        :return:
        """

        # Validate input path_file object for mol
        mol = mol_validate_file_object(request['mol'])

        # Run smartcyp
        smartcyp = SmartCypRunner(log=self.log, workdir=request.get('workdir'))
        result_dict = smartcyp.run(mol['content'],
                                   is_smiles=mol['extension'] == 'smi',
                                   output_format=request['output_format'],
                                   noempcorr=request['noempcorr'])

        # Remove temporary working directory
        smartcyp.delete()

        return result_dict

    @endpoint('docking_info', 'info_request', 'info_response', options=RegisterOptions(invoke=u'roundrobin'))
    def plants_info(self, request, claims):
        """
        Returns an informative summary of the supported PLANTS version
        and configuration.
        """

        return plants_version_info()

    @endpoint('docking', 'docking_request', 'docking_response', options=RegisterOptions(invoke='roundrobin'))
    def plants_docking(self, request, claims):
        """
        Perform a PLANTS (Protein-Ligand ANT System) molecular docking.
        For a detail description of the input see the file:
        schemas/endpoints/docking-request.v1.json
        """

        # Validate input path_file object for protein and ligand file
        protein_file = mol_validate_file_object(request['protein_file'])
        ligand_file = mol_validate_file_object(request['ligand_file'])

        # Run docking
        docking = PlantsDocking(workdir=os.environ.get('BASE_WORK_DIR', request.get('base_workdir')), **request)
        success = docking.run(protein_file['content'], ligand_file['content'])

        if success:
            status = 'completed'
            results = docking.get_results()
            output = {key: encode_file(value) for key, value in results.items()}

        else:
            self.log.error('PLANTS docking failed')
            return {'status': 'failed'}

        return {'status': status, 'output': output}

    @endpoint('docking_statistics', 'docking_statistics_request', 'docking_statistics_response',
              options=RegisterOptions(invoke='roundrobin'))
    def plants_docking_statistics(self, request, claims):
        """
        Return PLANTS docking statistics for particular docking solutions run previously.
        Clustering will also be redone and optionally adjusted.
        """

        base_dir = os.environ.get('BASE_WORK_DIR', request.get('base_workdir'))
        base_path = list(set([os.path.dirname(p) for p in request['paths']]))
        if len(base_path) > 1:
            self.log.error('Unable to combine results for more then one docking run')
            return {'status': 'failed'}

        if not os.path.exists(os.path.join(base_dir or '', base_path[0])):
            self.log.error('Docking results (no longer) exist: {0}'.format(os.path.basename(base_path[0])))
            return {'status': 'failed'}

        docking = PlantsDocking(base_work_dir=base_dir)
        docking.workdir = base_path[0]
        docking.update(**request)

        results = docking.get_results(structures=request['paths'])
        if results:
            return {'status': 'completed', 'output': results}

        self.log.error('PLANTS docking failed')
        return {'status': 'failed'}

    @endpoint('docking_structures', 'docking_structures_request', 'docking_structures_response',
              options=RegisterOptions(invoke='roundrobin'))
    def plants_docking_structures(self, request, claims):
        """
        Return PLANTS docking statistics for particular docking solutions run previously.
        Clustering will also be redone and optionally adjusted.
        """

        base_dir = os.environ.get('BASE_WORK_DIR', request.get('base_workdir'))
        base_path = list(set([os.path.dirname(p) for p in request['paths']]))
        if len(base_path) > 1:
            self.log.error('Unable to combine results for more then one docking run')
            return {'status': 'failed'}

        if not os.path.exists(os.path.join(base_dir or '', base_path[0])):
            self.log.error('Docking results (no longer) exist: {0}'.format(os.path.basename(base_path[0])))
            return {'status': 'failed'}

        docking = PlantsDocking(base_work_dir=base_dir)
        docking.workdir = base_path[0]
        results = docking.get_structures(**request)

        return results

    @endpoint('spores', 'spores_request', 'spores_response', options=RegisterOptions(invoke='roundrobin'))
    def spores_run(self, request, claims):
        """
        Perform a SPORES (Structure PrOtonation and REcognition System) structure preparation.
        For a detail description of the input see the file:
        schemas/endpoints/spores-request.v1.json
        """

        # Validate input path_file object for mol
        mol = mol_validate_file_object(request['mol'])

        spores = SporesRunner(log=self.log, base_work_dir=request.get('workdir'))
        result_dict = spores.run(mol['content'], mode=request['spores_mode'], input_format=request['input_format'])
        spores.delete()

        return result_dict

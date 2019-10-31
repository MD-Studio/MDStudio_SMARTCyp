# -*- coding: utf-8 -*-

"""
file: wamp_services.py

WAMP service methods the module exposes.
"""

import os

from autobahn.wamp import RegisterOptions
from mdstudio.api.endpoint import endpoint
from mdstudio.component.session import ComponentSession

from mdstudio_smartcyp.combined_prediction import CombinedPrediction
from mdstudio_smartcyp.smartcyp_run import SmartCypRunner, smartcyp_version_info
from mdstudio_smartcyp.plants_run import PlantsDocking, plants_version_info
from mdstudio_smartcyp.spores_run import spores_version_info, SporesRunner
from mdstudio_smartcyp.utils import mol_validate_file_object, MDStudioException


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

    @endpoint('smartcyp_info', 'smartcyp_info_request', 'smartcyp_info_response',
              options=RegisterOptions(invoke=u'roundrobin'))
    def smartcyp_info(self, request, claims):
        """
        Returns an informative summary of the supported SMARTCyp version
        and configuration.
        """

        info = smartcyp_version_info()
        if info:
            return {'status': 'completed', 'info': info}

        return {'status': 'failed'}

    @endpoint('smartcyp', 'smartcyp_request', 'smartcyp_response', options=RegisterOptions(invoke=u'roundrobin'))
    def smartcyp_prediction(self, request, claims):
        """
        Run a SMARTCyp prediction for a molecule

        :param request:
        :param claims:
        :return:
        """

        # Validate input path_file object for mol
        ligand_file = mol_validate_file_object(request['ligand_file'])

        # Run smartcyp
        base_dir = os.environ.get('BASE_WORK_DIR', request.get('base_work_dir'))
        smartcyp = SmartCypRunner(log=self.log, base_work_dir=base_dir)
        result_dict = smartcyp.run(ligand_file['content'],
                                   is_smiles=ligand_file['extension'] == 'smi',
                                   output_format=request['output_format'],
                                   noempcorr=request['noempcorr'])

        # Remove temporary working directory
        smartcyp.delete()

        result_dict['status'] = 'completed' if result_dict['result'] is not None else 'failed'
        return result_dict

    @endpoint('docking_info', 'docking_info_request', 'docking_info_response',
              options=RegisterOptions(invoke=u'roundrobin'))
    def plants_info(self, request, claims):
        """
        Returns an informative summary of the supported PLANTS version
        and configuration.
        """

        info = plants_version_info()
        if info:
            return {'status': 'completed', 'info': info}

        return {'status': 'failed'}

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
        base_dir = os.environ.get('BASE_WORK_DIR', request.get('base_work_dir'))

        for drop_key in ('protein_file', 'ligand_file', 'base_work_dir'):
            if drop_key in request:
                del request[drop_key]

        docking = PlantsDocking(log=self.log, base_work_dir=base_dir, **request)

        if docking.run(protein_file['content'], ligand_file['content']):
            return {'status': 'completed', 'result': docking.get_results()}

        self.log.error('PLANTS docking failed')
        return {'status': 'failed'}

    @endpoint('docking_statistics', 'docking_statistics_request', 'docking_statistics_response',
              options=RegisterOptions(invoke='roundrobin'))
    def plants_docking_statistics(self, request, claims):
        """
        Return PLANTS docking statistics for particular docking solutions run previously.
        Clustering will also be redone and optionally adjusted.
        """

        base_dir = os.environ.get('BASE_WORK_DIR', request.get('base_work_dir'))
        if 'base_work_dir' in request:
            del request['base_work_dir']

        docking = PlantsDocking(log=self.log, base_work_dir=base_dir, **request)

        try:
            results = docking.get_results(structures=request.get('paths'))
        except MDStudioException as error:
            self.log.error(repr(error))
            return {'status': 'failed'}

        if results:
            return {'status': 'completed', 'result': results}

        self.log.error('PLANTS docking failed')
        return {'status': 'failed'}

    @endpoint('docking_structures', 'docking_structures_request', 'docking_structures_response',
              options=RegisterOptions(invoke='roundrobin'))
    def plants_docking_structures(self, request, claims):
        """
        Return PLANTS docked structures
        """

        base_dir = os.environ.get('BASE_WORK_DIR', request.get('base_work_dir'))
        docking = PlantsDocking(log=self.log, base_work_dir=base_dir)

        try:
            results = docking.get_structures(structures=request.get('paths'),
                                             output_format=request.get('output_format', 'mol2'),
                                             include_protein=request.get('include_protein', False),
                                             create_ensemble=request.get('create_ensemble', True))
        except MDStudioException as error:
            self.log.error(repr(error))
            return {'status': 'failed'}

        # Encode files
        if not isinstance(results, list):
            results = [results]

        output = [{"path": None, "extension": request.get('output_format', 'mol2'), "content": mol, "encoding": "utf8"}
                  for mol in results]

        return {'status': 'completed', 'result': output}

    @endpoint('spores_info', 'spores_info_request', 'spores_info_response',
              options=RegisterOptions(invoke=u'roundrobin'))
    def spores_info(self, request, claims):
        """
        Returns an informative summary of the supported SPORES version
        and configuration.
        """

        info = spores_version_info()
        if info:
            return {'status': 'completed', 'info': info}

        return {'status': 'failed'}

    @endpoint('spores', 'spores_request', 'spores_response', options=RegisterOptions(invoke='roundrobin'))
    def spores_run(self, request, claims):
        """
        Perform a SPORES (Structure PrOtonation and REcognition System) structure preparation.
        For a detail description of the input see the file:
        schemas/endpoints/spores-request.v1.json
        """

        # Validate input path_file object for mol
        mol = mol_validate_file_object(request['mol'])

        spores = SporesRunner(log=self.log, base_work_dir=os.environ.get('BASE_WORK_DIR', request.get('base_work_dir')))
        try:
            result_dict = spores.run(mol['content'], mode=request['spores_mode'], input_format=request['input_format'])
        except Exception as e:
            self.log.error(str(e))
            return {'status': 'failed'}
        finally:
            spores.delete()

        return {'status': 'completed', 'result': result_dict}

    @endpoint('som_prediction', 'som_prediction_request', 'som_prediction_response',
              options=RegisterOptions(invoke='roundrobin'))
    def som_prediction(self, request, claims):
        """
        Run a REST based SOM prediction run
        """

        # Validate input path_file object for protein and ligand file
        ligand_file = mol_validate_file_object(request['ligand_file'])

        # Run combined structure/reactivity prediction
        base_dir = os.environ.get('BASE_WORK_DIR', request.get('base_work_dir'))
        sompred = CombinedPrediction(base_work_dir=base_dir, **request)
        prediction = sompred.run(ligand_file, filter_clusters=request['filter_clusters'])

        if prediction:
            return {'status': 'completed', 'output': prediction}

        self.log.error('SOM prediction failed')
        return {'status': 'failed'}

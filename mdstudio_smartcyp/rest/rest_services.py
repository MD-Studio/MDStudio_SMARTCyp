# -*- coding: utf-8 -*-

"""
file: rest_services.py

REST service methods the module exposes.
"""

import os

from werkzeug import FileStorage

from mdstudio_smartcyp.smartcyp_run import SmartCypRunner
from mdstudio_smartcyp.plants_run import PlantsDocking
from mdstudio_smartcyp.spores_run import SporesRunner
from mdstudio_smartcyp.utils import mol_validate_file_object, MDStudioException
from mdstudio_smartcyp.combined_prediction import CombinedPrediction


def som_prediction(ligand_file, base_work_dir=None, cyp='3A4', filter_clusters=True, explicit_oxygen=False,
                   smartcyp_score_label=None, **kwargs):
    """
    Run a REST based SOM prediction run

    :param ligand_file:          ligand structure MOL2 file
    :type ligand_file:           :py:str
    :param base_work_dir:        optional work directory to (temporary) store
                                 PLANTS docking results.
    :type base_work_dir:         :py:str
    :param cyp:                  CYP isoform to make prediction for
    :type cyp:                   :py:str
    :param filter_clusters:      make prediction for clustered docking results
                                 only
    :type filter_clusters:       :py:bool
    :param explicit_oxygen:      Use protein structure with explicit oxygen on
                                 the heme
    :type explicit_oxygen:       :py:bool
    :param smartcyp_score_label: SMARTCyp output 'score' values to use for
                                 prediction
    :type smartcyp_score_label:  :py:str
    :param kwargs:               additional docking configuration parameters
    :type kwargs:                :py:dict

    :return:                     PLANTS docking statistics (content of features.csv)
                                 file.
    :rtype:                      :py:dict
    """

    if isinstance(ligand_file, FileStorage):
        ligand_file = ligand_file.read().decode('utf-8')
    else:
        return 'Unsupported protein file structure: {0}'.format(type(ligand_file)), 401

    # Run combined structure/reactivity prediction
    sompred = CombinedPrediction(base_work_dir=os.environ.get('BASE_WORK_DIR', base_work_dir), cyp=cyp,
                                 explicit_oxygen=explicit_oxygen, smartcyp_score_label=smartcyp_score_label, **kwargs)
    prediction = sompred.run(ligand_file, filter_clusters=filter_clusters)

    if prediction:
        return prediction

    return 'SOM prediction failed', 401


def plants_docking(protein_file, ligand_file, base_work_dir=None, **kwargs):
    """
    Run a REST based PLANTS docking run

    :param protein_file:    protein structure MOL2 file
    :type protein_file:     :py:str
    :param ligand_file:     ligand structure MOL2 file
    :type ligand_file:      :py:str
    :param base_work_dir:   optional work directory to (temporary) store PLANTS
                            docking results.
    :type base_work_dir:    :py:str

    :return:                PLANTS docking statistics (content of features.csv)
                            file.
    :rtype:                 :py:dict
    """

    if isinstance(protein_file, FileStorage):
        protein_file = protein_file.read().decode('utf-8')
    else:
        return 'Unsupported protein file structure: {0}'.format(type(protein_file)), 401

    if isinstance(ligand_file, FileStorage):
        ligand_file = ligand_file.read().decode('utf-8')
    else:
        return 'Unsupported protein file structure: {0}'.format(type(ligand_file)), 401

    # Run docking
    docking = PlantsDocking(base_work_dir=os.environ.get('BASE_WORK_DIR', base_work_dir), **kwargs)

    results = None
    try:
        if docking.run(protein_file, ligand_file):
            results = docking.get_results()
    except MDStudioException as error:
        return repr(error), 401

    if results:
        return results

    return 'PLANTS docking failed', 401


def plants_docking_statistics(paths=None, **kwargs):
    """
    Return PLANTS docking statistics for particular docking solutions run previously.
    Clustering will also be redone and optionally adjusted.

    :param paths: list of docking solution paths
    :type paths:  :py:list

    :return:      PLANTS docking statistics (content of features.csv) file.
    :rtype:       :py:dict
    """

    docking = PlantsDocking(base_work_dir=os.environ.get('BASE_WORK_DIR'), **kwargs)

    try:
        results = docking.get_results(structures=paths)
    except MDStudioException as error:
        return repr(error), 401

    if results:
        return results

    return 'PLANTS docking failed', 401


def plants_docking_structures(paths=None, output_format='mol2', include_protein=False, **kwargs):
    """
    Return PLANTS docking statistics for particular docking solutions run previously.
    Clustering will also be redone and optionally adjusted.

    :param paths:           list of docking solution paths
    :type paths:            :py:list
    :param output_format:   return ensemble in MOL2 or PDB format
    :type output_format:    :py:str
    :param include_protein: return ensemble including the protein
                            structure or only the ligand structures
    :type include_protein:  :py:bool

    :return:                docking structures as combined Tripos MOL2 file
    :rtype:                 :py:str
    """

    docking = PlantsDocking(base_work_dir=os.environ.get('BASE_WORK_DIR'), **kwargs)

    try:
        results = docking.get_structures(paths, output_format=output_format, include_protein=include_protein)
    except MDStudioException as error:
        return repr(error), 401

    if results:
        return results

    return 'PLANTS docking failed', 401


def smartcyp_prediction(ligand_file=None, smiles=None, output_format='json', noempcorr=False, output_png=False):
    """
    Run a REST based SMARTCyp prediction for a molecule

    :param ligand_file:   molecule to make prediction for
    :type ligand_file:    :py:str
    :param smiles:        molecule as SMILES or InChI string
    :type smiles:         :py:str
    :param output_format: output format as CSV, JSON or HTML
    :type output_format:  :py:str
    :param noempcorr:     do not use the empirical N-oxidation correction
                          (smartcyp >= v2.3)
    :type noempcorr:      :py:bool
    :param output_png:    export PNG image files for the prediction

    :return:              SMARTCyp results as JSON, CSV or HTML
    :rtype:               :py:str
    """

    # Input molecule is either a SMILES/InChI string or file object
    path_file_object = {'content': smiles, 'extension': 'smi' if smiles else None, 'path': None}
    if isinstance(ligand_file, FileStorage):
        path_file_object['content'] = ligand_file.read().decode('utf-8')
        if '.' in os.path.basename(ligand_file.filename):
            path_file_object['extension'] = ligand_file.filename.split('.')[-1]
        path_file_object['path'] = ligand_file.filename

    if path_file_object['content'] is None:
        return "Ligand input from file 'mol' or as 'smiles' string required", 401

    # Validate path_file_object: valid SMILES or InChI strings
    path_file_object = mol_validate_file_object(path_file_object)

    # Run SMARTCyp
    smartcyp = SmartCypRunner()
    try:
        result_dict = smartcyp.run(path_file_object['content'],
                                   is_smiles=path_file_object['extension'] == 'smi',
                                   input_format=path_file_object.get('extension', 'mol2'),
                                   output_format=output_format,
                                   noempcorr=noempcorr,
                                   output_png=output_png)
    except Exception as e:
        return str(e), 500
    finally:
        smartcyp.delete()

    return result_dict


def spores_run(mol, spores_mode='complete', input_format='mol2', base_work_dir=None):
    """
    Perform a SPORES (Structure PrOtonation and REcognition System) structure preparation.
    For a detail description of the input see the file:
    schemas/endpoints/spores-request.v1.json

    :param mol:           SPORES input structure
    :type mol:            :py:str
    :param spores_mode:   SPORES running mode
    :type spores_mode:    :py:str
    :param input_format:  Input structure format
    :type input_format:   :py:str
    :param base_work_dir: optional work directory to (temporary) store SPORES
                          results.
    :type base_work_dir:  :py:str

    :return:             SPORES processed structure
    :rtype:              :py:str
    """

    path_file_object = {'content': None, 'extension': input_format, 'path': None}
    if isinstance(mol, FileStorage):
        path_file_object['content'] = mol.read().decode('utf-8')
        if '.' in os.path.basename(mol.filename):
            path_file_object['extension'] = mol.filename.split('.')[-1]
        path_file_object['path'] = mol.filename

    if path_file_object['content'] is None:
        return "Structure input from file 'mol' required", 401

    spores = SporesRunner(base_work_dir=os.environ.get('BASE_WORK_DIR', base_work_dir))
    try:
        result_dict = spores.run(path_file_object['content'], mode=spores_mode, input_format=input_format)
    except Exception as e:
        return str(e), 500
    finally:
        spores.delete()

    if result_dict is None or 'content' not in result_dict:
        return 'SPORES processing of file {0} failed'.format(path_file_object['path']), 401

    return result_dict['content']

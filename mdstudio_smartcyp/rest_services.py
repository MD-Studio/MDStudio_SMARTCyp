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
from mdstudio_smartcyp.utils import mol_validate_file_object


def plants_docking(protein_file, ligand_file, bindingsite_center, workdir=None, **kwargs):
    """
    Run a REST based PLANTS docking run
    :return:
    """

    if isinstance(protein_file, FileStorage):
        protein_file = protein_file.read().decode('utf-8')
    else:
        return 'Unsuported protein file structure: {0}'.format(type(protein_file)), 401

    if isinstance(ligand_file, FileStorage):
        ligand_file = ligand_file.read().decode('utf-8')
    else:
        return 'Unsuported protein file structure: {0}'.format(type(ligand_file)), 401

    # Format the binding site center
    try:
        kwargs['bindingsite_center'] = [float(n) for n in bindingsite_center.split(',')]
    except BaseException:
        return 'Malformed bindingsite_center: {0}'.format(bindingsite_center)

    # Run docking
    docking = PlantsDocking(workdir=workdir, **kwargs)
    success = docking.run(protein_file, ligand_file)

    if success:
        results = docking.results()
        docking.delete()
        return results

    docking.delete()
    return 'PLANTS docking failed', 401


def smartcyp_prediction(mol=None, smiles=None, output_format='json', noempcorr=False, output_png=False):
    """
    Run a REST based SMARTCyp prediction for a molecule

    :param mol:           molecule to make prediction for
    :type mol:            :py:str
    :param smiles:        molecule as SMILES or InChI string
    :type smiles:         :py:str
    :param is_smiles:     is the molecule a SMILES string
    :type is_smiles:      :py:bool
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
    if isinstance(mol, FileStorage):
        path_file_object['content'] = mol.read().decode('utf-8')
        if '.' in os.path.basename(mol.filename):
            path_file_object['extension'] = mol.filename.split('.')[-1]
        path_file_object['path'] = mol.filename

    if path_file_object['content'] is None:
        return "Ligand input from file 'mol' or as 'smiles' string required", 401

    # Validate path_file_object: valid SMILES or InChI strings
    path_file_object = mol_validate_file_object(path_file_object)

    # Run SMARTCyp
    smartcyp = SmartCypRunner()
    try:
        result_dict = smartcyp.run(path_file_object['content'],
                                   is_smiles=path_file_object['extension'] == 'smi',
                                   output_format=output_format,
                                   noempcorr=noempcorr,
                                   output_png=output_png)
    except Exception as e:
        return str(e), 500
    finally:
        smartcyp.delete()

    return result_dict


def spores_run(mol, spores_mode, input_format='mol2'):
    """
    Perform a SPORES (Structure PrOtonation and REcognition System) structure preparation.
    For a detail description of the input see the file:
    schemas/endpoints/spores-request.v1.json

    :param mol:          SPORES input structure
    :type mol:           :py:str
    :param spores_mode:  SPORES running mode
    :type spores_mode:   :py:str
    :param input_format: Input structure format
    :type input_format:  :py:str

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

    spores = SporesRunner()
    try:
        result_dict = spores.run(path_file_object['content'], spores_mode, input_format=input_format)
    except Exception as e:
        return str(e), 500
    finally:
        spores.delete()

    if result_dict is None:
        return 'SPORES processing of file {0} failed'.format(path_file_object['path']), 401

    return result_dict

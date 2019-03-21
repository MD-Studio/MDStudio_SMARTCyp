# -*- coding: utf-8 -*-

"""
file: rest_services.py

REST service methods the module exposes.
"""

import os

from werkzeug import FileStorage

from mdstudio_smartcyp.smartcyp_run import SmartCypRunner, mol_validate_file_object


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
        path_file_object['filepath'] = mol.filename

    if path_file_object['content'] is None:
        return "Ligand input from file 'mol' or as 'smiles' string required", 401

    # Validate path_file_object: valid SMILES or InChI strings
    path_file_object = mol_validate_file_object(path_file_object)

    # Run SMARTCyp
    smartcyp = SmartCypRunner()
    result_dict = smartcyp.run(path_file_object['content'],
                               is_smiles=path_file_object['extension'] == 'smi',
                               output_format=output_format,
                               noempcorr=noempcorr,
                               output_png=output_png)

    return result_dict

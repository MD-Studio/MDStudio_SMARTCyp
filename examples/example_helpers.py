# -*- coding: utf-8 -*-

import os
import logging

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

no_wrap_div = '<div style="white-space: nowrap">{}{}</div>'


def show_predictions(mol, prob, cutoff=0.75):
    """
    Build 2D depiction of rdkit molecule
    Label atoms with index number and SOM prediction probabilities if
    above cutoff.
    """

    prob_dict = dict([(int(a.split('.')[1]), b) for a, b in prob.items()])

    Chem.rdDepictor.Compute2DCoords(mol)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx() + 1
        atom.SetProp('molAtomMapNumber', '{0}-{1:.2f}'.format(idx, prob_dict.get(idx, 0.0)))

    drawer = rdMolDraw2D.MolDraw2DSVG(400, 250)
    drawer.DrawMolecule(mol, highlightAtoms=[i - 1 for i in prob_dict if prob_dict[i] >= cutoff])
    drawer.FinishDrawing()

    return drawer.GetDrawingText().replace('svg:', '')


def get_dataset():

    data = {'aminopyrine': {'som': [12, 13], 'cyp': '2D6', 'smiles': 'CN(C1=C[NH+](N(C1=O)c1ccccc1)C)C'},
            'estradiol': {'som': [8, 9, 11, 15], 'cyp': '2D6', 'smiles': 'Oc1ccc2c(c1)CC[C@@H]1[C@@H]2CC[C@]2([C@H]1CC[C@@H]2O)C'},
            'methadone': {'som': [20, 21], 'cyp': '2D6', 'smiles': 'CCC(=O)C(c1ccccc1)(c1ccccc1)C[C@H]([NH+](C)C)C'}}

    for case in data:
        mol2 = '{0}.mol2'.format(case)
        if os.path.exists(mol2):
            data[case]['file'] = open(mol2, 'rb').read()
        else:
            logging.warn('Test case {0} does not exist'.format(mol2))
            del data[case]

    return data

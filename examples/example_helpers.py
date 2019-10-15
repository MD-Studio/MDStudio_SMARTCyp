# -*- coding: utf-8 -*-

import os
import logging
import pandas

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

no_wrap_div = '<div style="white-space: nowrap">{}{}{}</div>'


def show_predictions(mol, prob, cutoff=0.75, show_prob_label=False):
    """
    Build 2D depiction of rdkit molecule
    Label atoms with index number and SOM prediction probabilities if
    above cutoff.
    """

    prob_dict = dict([(int(a.split('.')[1]), b) for a, b in prob.items()])

    Chem.rdDepictor.Compute2DCoords(mol)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx() + 1
        if show_prob_label:
            atom.SetProp('molAtomMapNumber', '{0}-{1:.2f}'.format(idx, prob_dict.get(idx, 0.0)))
        else:
            atom.SetProp('molAtomMapNumber', '{0}'.format(idx))

    drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
    drawer.DrawMolecule(mol, highlightAtoms=[i - 1 for i in prob_dict if prob_dict[i] >= cutoff])
    drawer.FinishDrawing()

    return drawer.GetDrawingText().replace('svg:', '')


def get_dataset():

    data = {'aminopyrine': {'som': [12, 13], 'cyp': '2D6', 'smiles': 'CN(C1=C[NH+](N(C1=O)c1ccccc1)C)C'},
            'estradiol': {'som': [8, 9, 11, 15], 'cyp': '2D6', 'smiles': 'Oc1ccc2c(c1)CC[C@@H]1[C@@H]2CC[C@]2([C@H]1CC[C@@H]2O)C'},
            'methadone': {'som': [20, 21], 'cyp': '2D6', 'smiles': 'CCC(=O)C(c1ccccc1)(c1ccccc1)C[C@H]([NH+](C)C)C'}}

    for case in data:
        mol2 = '{0}.mol2'.format(case)
        sdf = '{0}.sdf'.format(case)
        if os.path.exists(mol2) and os.path.exists(sdf):
            data[case]['mol2'] = open(mol2, 'rb').read()
            data[case]['sdf'] = open(sdf, 'rb').read()
        else:
            logging.warn('Test case {0} does not exist'.format(case))
            del data[case]

    return data


def process_fame_results(fame_out, df):
    """
    Add Fame results to pandas DataFrame and determine cutoff propensity
    """
    atom_labels = dict([(int(a.split('.')[1]), a) for a in df.index])

    fame_pred = {}
    cutoff = 1.0
    if 'predictions' in fame_out:

        for pred in fame_out['predictions']:
            fame_pred = dict([(atom_labels[a['atomID']], a['probability']) for a in pred['atomPredictions']])

            fame_positives = [a['probability'] for a in pred['atomPredictions'] if a['decision']]
            if (fame_positives):
                cutoff = min(fame_positives)

            break

    df['Fame'] = pandas.Series(fame_pred)

    return df, cutoff


def style_dataframe(df, smartcyp_cutoff=0.75, docking_cutoff=0.75, fame_cutoff=1.0):

    def highlight_som(df, color='#f08783'):
        df1 = pandas.DataFrame('', index=df.index, columns=df.columns)
        df1.loc[(df['SMARTCyp'] >= smartcyp_cutoff), 'SMARTCyp'] = 'background-color: {}'.format(color)
        df1.loc[(df['Docking'] >= docking_cutoff), 'Docking'] = 'background-color: {}'.format(color)
        df1.loc[(df['Fame'] >= fame_cutoff), 'Fame'] = 'background-color: {}'.format(color)
        return df1

    return df.style.apply(highlight_som, axis=None).set_precision(3)

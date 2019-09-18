# -*- coding: utf-8 -*-

"""
Combining SMARTCyp SOM prediction with PLANTS docking based SOM predictions
"""

import os
import logging
import pandas

from fnmatch import fnmatch
from interact import System
from interact.interactions.charged import eval_heme_coordination

from mdstudio_smartcyp import __module__, __package_path__
from mdstudio_smartcyp.smartcyp_run import SmartCypRunner
from mdstudio_smartcyp.plants_run import PlantsDocking
from mdstudio_smartcyp.utils import (parse_tripos_atom, merge_protein_ligand_mol2,
                                     hydrophobic_atom_count, molecular_weight)

logger = logging.getLogger(__module__)

# Cyp conformation decision tree
cyp_conf = {'3A4': [{'min_mw': 0, 'max_mw': 354, 'conf': '3UA1_apo_5901.mol2'},
                    {'min_mw': 354, 'max_mw': 500, 'conf': '3UA1_apo_6091.mol2'},
                    {'min_mw': 354, 'max_mw': 9999, 'conf': '3UA1_BCP_3521.mol2'}],
            '1A2': [{'conf': '1A2_nathan.mol2'}],
            '2D6': [{'min_mw': 0, 'max_mw': 280, 'conf': '2D6_PPD_70_216.mol2'},
                    {'min_mw': 280, 'max_mw': 9999, 'max_hydrophob': 18, 'conf': '2D6_CHZ_170_79.mol2'},
                    {'min_mw': 280, 'max_mw': 9999, 'conf': '2D6_TMF_70_3.mol2'}]
            }


class CombinedPrediction(object):
    """
    Class for CYP based SOM predictions combining reactivity based SMARTCyp
    predictions with PLANTS structure based docking predictions.

    :param log:                  Python logger instance
    :type log:                   :py:logging
    :param base_work_dir:        optional work directory to (temporary) store
                                 PLANTS docking results.
    :type base_work_dir:         :py:str
    :param cyp:                  CYP isoform to make prediction for
    :type cyp:                   :py:str
    :param filter_clusters:      make prediction for clustered docking results
                                 only
    :type filter_clusters:       :py:bool
    :param smartcyp_score_label: SMARTCyp output 'score' values to use for
                                 prediction
    :type smartcyp_score_label:  :py:str
    :param kwargs:               additional docking configuration parameters
    :type kwargs:                :py:dict
    """

    def __init__(self, log=logger, base_work_dir=None, cyp='3A4', smartcyp_score_label='Score',  **kwargs):

        self.log = log
        self.base_work_dir = base_work_dir
        self.docking_config = kwargs
        self.smartcyp_score_label = smartcyp_score_label
        self._workdir = None

        self.cyp = cyp.upper()
        self.smartcyp_results = None
        self.docking_results = None
        self.combined = None

    def cyp_decision_tree(self, lig_mol2_atoms):
        """
        Determine best CYP isoform structure conformation to use based on the
        decision trees from previous research that use ligand molecular weight
        and optionally the number of hydrophobic atoms to decide on the best
        conformation. Decision trees used:

          CYP   Mol Weight / nr Hfob    Conformation                    Reference
        * 3A4   < 354                   3UA1_apo_5901                   Bon
        * 3A4   354 - 500               3UA1_apo_6091                   Bon
        * 3A4   > 500                   3UA1_BCP_3521   3UA1_BCP_5291   Bon
        * 1A2   overall                 1A2_nathan                      Nathan
        * 2D6   < 280                   2D6_PPD_70_216                  Hritz
        * 2D6   > 280, < 18 n_hydroph   2D6_CHZ_170_79                  Hritz
        * 2D6   > 280, > 18 n_hydroph   2D6_TMF_70_3                    Hritz

        :param lig_mol2_atoms:  Tripos MOL2 atom records as returned by `parse_tripos_atom`
        :type lig_mol2_atoms:   :py:dict

        :return:                Path to protein structure
        :rtype:                 :py:str
        """

        molw = molecular_weight(lig_mol2_atoms)
        n_hydrophob = hydrophobic_atom_count(lig_mol2_atoms)

        choice = None
        for conf in cyp_conf[self.cyp]:
            if conf['min_mw'] < molw < conf['max_mw']:

                if n_hydrophob > conf.get('max_hydrophob', 9999):
                    continue

                choice = conf['conf']

        self.log.info('Use {0} conformation {1}, molecular weight: {2:.3f} and hydrophobic atom count: {3}'.format(
            self.cyp, choice, molw, n_hydrophob))

        protein = open(os.path.join(__package_path__, 'data/{0}'.format(choice)), 'r').read()
        return protein

    def combine_docking_smartcyp(self, hemecoor, pose_count):
        """
        Combine SMARTCyp reactivity based SOM prediction with docking based
        prediction.

        This function aims to provide a normalized score between 0 and 1
        representing the probability of a atom being a SOM for SMARTCyp and
        docking based prediction strategies by:

        * PLANTS docking: count possible SOMs identified by MDInteract heme
          coordination for every docking pose (all poses, clustered poses or
          custom selection). Dividing the count sum by the number of poses
          for every atom. Normalize by the atom with the maximum value from
          the previous output.
        *

        Count the number of times ligand atoms participate in possible heme
        coordination according to the MDInteract results.
        Normalize the counts between 0 and 1 by dividing by the number of
        evaluated poses (`pose_count`). The latter is either all generated
        poses or the ones filtered for non-clustered poses.

        Normalize the SMARTCyp score used for prediction
        (`smartcyp_score_label`) between 0 and 1 for easy comparison to docking
        results.

        :param hemecoor:    MDInteract heme coordination evaluation for docking
                            poses
        :type hemecoor:     :pandas:DataFrame
        :param pose_count:  number of evaluated docking poses (population)
        :type pose_count:   :py:int

        :return:            accumulated SOM prediction
        :rtype:             :pandas:DataFrame
        """

        ec = ['Docking', 'SMARTCyp']
        poses = list(hemecoor['pose'].unique())
        docmat = pandas.DataFrame(0, index=self.smartcyp_results['Atom_id'], columns=poses + ec)

        # Add MDInteract predicted SOMS for each docking pose to docmat
        for pose in poses:
            pred = hemecoor[hemecoor['pose'] == pose]['source', 'serial']
            if not pred.empty:
                docmat.loc[pred, pose] = 1

        # Normalize som counts
        norm = docmat[poses].sum(axis=1) / pose_count
        docmat['Docking'] = norm / max(norm)

        # Add SMARTCyp prediction. Normalize defined score column
        norm_smartcyp_score = 1.0 - (self.smartcyp_results[self.smartcyp_score_label] /
                                     max(self.smartcyp_results[self.smartcyp_score_label]))
        docmat['SMARTCyp'] = norm_smartcyp_score.values

        # Change index and columns
        docmat.columns = list(self.docking_results[self.docking_results['POSE'].isin(poses)]['PATH']) + ec
        docmat.index = self.smartcyp_results.index

        return docmat

    def run(self, ligand, filter_clusters=True):
        """
        Run combined SOM prediction

        * Determine Cyp isoform using decision tree
        * Run SMARTCyp
        * Run PLANTS docking
        * Run MDInteract heme coordination on docking poses
        * Combine and return prediction results

        :param ligand:           ligand in Tripos MOL2 format
        :type ligand:            :py:str
        :param filter_clusters:  filter docking psoes on clusters
        :type filter_clusters:   :py:bool

        :return:                 combined prediction results
        :rtype:                  :py:dict
        """

        # Get CYP isoform to use
        isoforms = [c for c in cyp_conf.keys() if fnmatch(self.cyp, c)]
        if not len(isoforms):
            self.log.error('Unsupported CYP isoform: {0}'.format(self.cyp))
            return
        elif len(isoforms) > 1:
            isoforms = [c for c in isoforms if '*' not in c]

        self.cyp = isoforms[0]

        # Determine protein conformation to use
        lig_mol2_atoms = parse_tripos_atom(ligand)
        protein = self.cyp_decision_tree(lig_mol2_atoms)

        # Run SMARTCyp
        smartcyp = SmartCypRunner()
        smartcyp_results = smartcyp.run(ligand, is_smiles=False)
        if smartcyp_results['result'] is None:
            self.log.error('Error running SMARTCyp')
            return

        # Import SMARTCyp results as Pandas DataFrame
        self.smartcyp_results = pandas.DataFrame.from_dict(smartcyp_results['result'], orient='index')

        # Perform PLANTS docking
        docking = PlantsDocking(base_work_dir=self.base_work_dir, bindingsite_center=[-0.989, 3.261, 0.826],
                                **self.docking_config)
        success = docking.run(protein, ligand)

        if success:
            self.docking_results = pandas.DataFrame.from_dict(docking.get_results(), orient='index')
            self.docking_results['POSE'] = [int(f.split('_')[-1]) for f in self.docking_results.index]
        else:
            self.log.error('Failed to run PLANTS docking')
            docking.delete()
            return

        # Store results data in docking results dir
        self.smartcyp_results.to_csv(os.path.join(docking.workdir, 'smartcyp.csv'))
        self.docking_results.to_csv(os.path.join(docking.workdir, 'docking.csv'))

        # Prepare PDB ensemble and system MOL2
        poses = list(self.docking_results['PATH'])
        ensemble_pdb = os.path.join(docking.workdir, 'ensemble.pdb')
        with open(ensemble_pdb, 'w') as epdb:
            epdb.write(docking.get_structures(poses, output_format='pdb', include_protein=True))

        system_mol2 = os.path.join(docking.workdir, 'system.mol2')
        with open(system_mol2, 'w') as smol:
            ligp = os.path.join(docking.base_work_dir, poses[0])
            smol.write(merge_protein_ligand_mol2(protein, open(ligp).read()))

        # Run heme-coordination detection
        molsys = System(ensemble_pdb, mol2file=system_mol2)
        lig_resname = set([n['subst_name'][0:3] for n in lig_mol2_atoms.values()])
        rings = None
        combined = []
        for frame, nr in molsys.iter_frames(auto_chunk=False):

            logging.info('Evaluate heme-coordination on docking pose: {0}'.format(nr + 1))

            frame.distances()
            ls = frame[frame['resName'].isin(lig_resname)]

            if rings is None:
                rings = ls.find_rings()

            cf = ls.contacts(ls.neighbours())
            cf = eval_heme_coordination(cf, molsys.topology, rings=rings)
            cf['pose'] = nr + 1

            f = cf[cf['contact'] != 'nd']
            if not f.empty:
                combined.append(f)

        hemecoor = pandas.concat(combined)
        hemecoor.to_csv(os.path.join(docking.workdir, 'hemecoor.csv'))

        pose_count = len(self.docking_results)
        if filter_clusters:
            clustered_poses = list(self.docking_results[self.docking_results['CLUSTER'] != 0]['POSE'])
            pose_count = len(clustered_poses)
            hemecoor = hemecoor[hemecoor['pose'].isin(clustered_poses)]

        self.combined = self.combine_docking_smartcyp(hemecoor, pose_count)
        self.combined.to_csv(os.path.join(docking.workdir, 'prediction.csv'))

        return self.combined.to_dict(orient='index')

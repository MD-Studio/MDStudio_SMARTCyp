# -*- coding: utf-8 -*-

"""
Classes for running the SPORES software:

  SPORES: Structure PrOtonation and REcognition System.
  "Influence of Protonation, Tautomeric, and Stereoisomeric States on Protein-Ligand Docking Results"
  J.Chem. Inf. Model. 49: 1535-1546 (2009). T. ten Brink and T.E Exner
"""

import logging
import os

from mdstudio_smartcyp import __module__, __spores_path__, __spores_version__, __spores_citation__
from mdstudio_smartcyp.utils import prepare_work_dir, RunnerBaseClass

logger = logging.getLogger(__module__)


def spores_version_info():
    """
    :return: information on the packaged SPORES version, supported models
             and software citation reference.
    :rtype:  :py:dict
    """

    info_dict = {'version': __spores_version__,
                 'citation': __spores_citation__,
                 'default_settings': {'input_format': 'mol2', 'spores_mode': 'complete'}}

    return info_dict


class SporesRunner(RunnerBaseClass):

    def __init__(self, log=logger, workdir=None):

        self.log = log
        self.workdir = prepare_work_dir(path=workdir, prefix='spores-')
        self.exec_path = __spores_path__

    def run(self, mol, mode='complete', input_format='mol2'):
        """
        Run SMARTCyp predictions

        Runs a SMARTCyp prediction for molecule `mol` in a system temporary
        directory and returns the content of the prediction .csv file as
        a dictionary.

        :param mol:           molecule to run SPORES on
        :type mol:            :py:str
        :param mode:          SPORES execution mode
        :type mode:           :py:str
        :param input_format:  Input structure format
        :type input_format:   :py:str

        :return:              SPORES processed structure
        :rtype:               :py:dict
        """

        # Copy files to working directory
        input_file = 'structure.{0}'.format(input_format)
        output_file = 'structure_{0}.mol2'.format(mode)
        with open(os.path.join(self.workdir, input_file), 'w') as protein_file:
            protein_file.write(mol)

        # Build and run SPORES CMD
        self.cmd_runner([self.exec_path, '--mode', mode, input_file, output_file])

        output_file_path = os.path.join(self.workdir, output_file)
        if os.path.isfile(output_file_path):

            result = {'extension': output_file_path.split('.')[-1],
                      'path': output_file_path,
                      'encoding': 'utf8'}

            with open(output_file_path, 'r') as spores_output:
                result['content'] = spores_output.read()

            return result
        else:
            self.log.error('SPORES failed to create output file {0}'.format(output_file))
            return None

# -*- coding: utf-8 -*-

import subprocess
import tempfile
import glob
import csv
import os
import shutil
import logging

from mdstudio_smartcyp import __smartcyp_path__, __module__

logger = logging.getLogger(__module__)
result_types = {'Ranking': int, '2DSASA': float, 'N+Dist': int, '2D6ranking': int, 'Energy': float,
                'Molecule': int, 'Span2End': int, 'Score': float, 'Relative Span': float,
                '2Cscore': float, 'Atom': str, '2Cranking': int, 'COODist': int, '2D6score': float}


class SmartCypRunner(object):

    def __init__(self, log=logger):

        self.log = log

    def _parse_csv(self, csvfile):
        """
        Parse SMARTCyp results .csv file

        Returns a JSON style dictionary with CSV rows keyed by atom identifier.
        Rows are keyed by the .csv header name. Values are coursed to their
        appropriate type with 'null' values to python None.

        :param csvfile: CSV file name
        :type csvfile:  :py:dict

        :return:        content of .csv file as dictionary
        :rtype:         :py:dict
        """

        csv_rows = {}
        with open(csvfile) as csvfile:
            reader = csv.DictReader(csvfile)
            title = reader.fieldnames
            for row in reader:
                drow = {title[i]: row[title[i]] for i in range(len(title))}

                # Type convert
                for key, value in drow.items():
                    if not key in result_types:
                        self.log.warning('Unknow results parameter: {0}'.format(key))

                    convert_to = result_types[key]
                    try:
                        # 'null' to None
                        if value == 'null':
                            drow[key] = None
                        else:
                            drow[key] = convert_to(value)
                    except ValueError:
                        self.log.warning('Unable to convert {0}: {1} to {2}'.format(key, value, convert_to))

                csv_rows[drow['Atom']] = drow

        return csv_rows

    def run(self, mol, is_smiles=False):
        """
        Run SMARTCyp predictions

        Runs a SMARTCyp prediction for molecule `mol` in a system temporary
        directory and returns the content of the prediction .csv file as
        a dictionary.

        :param mol:         molecule to make prediction for
        :type mol:          :py:str
        :param is_smiles:   is the molecule a SMILES string
        :type is_smiles:    :py:bool

        :return:            SMARTCyp prediction results
        :rtype:             :py:dict
        """

        # Make a temporary directory
        tempdir = tempfile.mkdtemp()
        self.log.debug('Make temporary directory at: {0}'.format(tempdir))

        # Build CMD
        cmd = ['java', '-jar', __smartcyp_path__, '-printall']
        if is_smiles:
            cmd.extend(['-smiles', mol])
            self.log.info('SMARTCyp prediction for SMILES string: {0}'.format(mol))

        else:
            with open('{0}/ligand.mol2'.format(tempdir), 'w') as ligfile:
                ligfile.write(mol)

            cmd.append(os.path.join(tempdir, 'ligand.mol2'))

        # Run SMARTCyp
        cmd.extend(['-outputdir', tempdir])
        self.log.debug('SMARTCyp CMD: {0}'.format(' '.join(cmd)))
        out = subprocess.call(cmd)

        # Get output CSV
        csvfile = glob.glob('{0}/*.csv'.format(tempdir))
        result = {}
        if len(csvfile):
            result = self._parse_csv(csvfile[0])
        else:
            self.log.error('SMARTCyp did not create a results .csv file')

        # Remove the temporary directory again
        shutil.rmtree(tempdir)

        return result

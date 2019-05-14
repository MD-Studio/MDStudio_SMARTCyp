# -*- coding: utf-8 -*-

"""
IO and subprocess related utility functions.
"""

import os
import sys
import subprocess
import logging
import re
import tempfile
import shutil
import csv

# Library and function compatibility
if sys.version_info[0] < 3:
    from cStringIO import StringIO
else:
    from io import StringIO

logger = logging.getLogger(__name__)
smiles_regex = re.compile('^([^J][A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+)$')


class RunnerBaseClass(object):

    def delete(self):
        """
        Remove the temporary working directory
        """

        # Remove the temporary directory again
        self.log.debug('Remove working directory: {0}'.format(self.workdir))

        shutil.rmtree(self.workdir)
        self.workdir = None

    def cmd_runner(self, cmd):
        """
        Common Command Line Interface runner

        :param cmd:      CLI commands
        :type cmd:       :py:list
        :param workdir:  working directory to run command in
        :type workdir:   :py:str

        :return:         subprocess stdout and stderr pipes
        """

        # Run cli command
        self.log.info('Execute cli process: {0}'.format(' '.join(cmd)))
        try:
            process = subprocess.Popen(cmd, cwd=self.workdir,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as err:
            self.log.error('Process failed:', err)
        else:
            self.log.info('Process returncode: {0}'.format(process.returncode))
            output, errors = process.communicate()

            if errors:
                self.log.error(errors.decode('utf-8'))
            if output:
                self.log.info(output.decode('utf-8'))


def _schema_to_data(schema, data=None, defdict=None):

    default_data = defdict or {}

    properties = schema.get('properties', {})

    for key, value in properties.items():
        if 'default' in value:
            if 'properties' in value:
                default_data[key] = _schema_to_data(value)
            else:
                default_data[key] = value.get('default')

    # Update with existing data
    if data:
        default_data.update(data)

    return default_data


def mol_validate_file_object(path_file):
    """
    Validate a MDStudio/REST path_file object

    - Check if 'content' is a InChI or SMILES string and set extension
    - If no 'content' check if path exists and import file

    :param path_file: path_file object
    :type path_file:  :py:dict

    :return:          validated path_file object
    :rtype:           :py:dict
    """

    content = path_file['content']
    if content is not None:

        # SMILES and InChI are single line strings
        if len(content.split('\n')) == 1:

            # Test for InChI type
            if content.startswith('InChI='):
                path_file['extension'] = 'inchi'

            # Test for SMILES
            if smiles_regex.match(path_file['content']):
                path_file['extension'] = 'smi'

    elif path_file['path'] is not None and os.path.exists(path_file['path']):

        with open(path_file['path']) as pf:
            path_file['content'] = pf.read()

    return path_file


def prepare_work_dir(path=None, prefix='', suffix='', create=True):
    """
    Prepare a unique workdir directory of the form:

        /dir/prefix-<unique ID>-suffix

    The directory basename is created using the `tempfile` package. If path
    is defined an attempt will be made to create the directory at that
    location otherwise the system `tmp` dir will be used.

    Using Pythons tempfile.mkdtemp a temporary directory is created in the
    most secure manner possible. There are no race conditions in the
    directory’s creation. The directory is readable, writable, and searchable
    only by the creating user ID.

    :param path:   target path to prepare the working directory in
    :type path:    :py:str
    :param prefix: prefix for directory basename
    :type prefix:  :py:str
    :param suffix: suffix for directory basename
    :type suffix:  :py:str
    :param create: create path if it does not exist
    :type create:  bool

    :return:       path to working directory
    :rtype:        :py:str
    :raise:        IOError, no writable directory or unable to create it
    """

    # Create temporary directory if path not defined
    if not path:
        path = tempfile.mkdtemp(prefix=prefix, suffix=suffix)

    # Check if target path exists
    else:
        path = os.path.abspath(path)
        if not os.path.exists(path):
            if create:
                logger.info('Path does not exists try creating it: {0}'.format(path))

                try:
                    os.makedirs(path)
                except OSError:
                    raise IOError('Unable to create path: {0}'.format(path))
            else:
                raise IOError('Path does not exist: {0}'.format(path))

        # Create temporary directory in path
        path = tempfile.mkdtemp(prefix=prefix, suffix=suffix, dir=path)

    # Is target directory writable
    if not os.access(path, os.W_OK):
        raise IOError('Working directory not writable: {0}'.format(path))

    logger.debug('Create working directory: {0}'.format(path))

    return path


def create_multi_mol2(mol2_file_paths):
    """
    Create a multi-molecule MOL2 file by concatenating
    single MOL2 files.

    :param mol2_file_paths: single MOL2 file paths
    :type mol2_file_paths:  :py:list

    :return:                multi-mol2 file
    :rtype:                 :py:str
    """

    multimol2 = StringIO()
    for path in mol2_file_paths:
        if os.path.exists(path):

            with open(path, 'r') as singlemol2:
                multimol2.write(singlemol2.read())

    multimol2.seek(0)
    return multimol2.read()


def import_plants_csv(result_dir, structures=None, files=('features.csv', 'ranking.csv')):
    """
    Import PLANTS results csv files

    :param result_dir: PLANTS results directory
    :type result_dir:  :py:str
    :param structures: only import files in structure selection, import all
                       by default
    :type structures:  :py:list
    :param files:      CSV files to import
    :type files:       :py:list, py:tuple

    :return:           docking results
    :rtype:            :py:dict
    """

    results = {}
    for resultcsv in files:
        resultcsv = os.path.join(result_dir, resultcsv)
        if os.path.isfile(resultcsv):

            with open(resultcsv, 'r') as csvfile:
                reader = csv.DictReader(csvfile)
                for i, row in enumerate(reader):
                    mol2 = row.get('TOTAL_SCORE', i)
                    path = os.path.join(result_dir, '{0}.mol2'.format(mol2))

                    # Only import structure selection if needed
                    if structures is not None and not path in structures:
                        continue

                    row = dict(row)
                    if None in row:
                        del row[None]

                    results[mol2] = row
                    results[mol2]['path'] = path
            break

    return results

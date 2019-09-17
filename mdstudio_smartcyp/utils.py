# -*- coding: utf-8 -*-

"""
MDStudio_SMARTCyp utility functions
"""

import os
import sys
import subprocess
import logging
import re
import tempfile
import shutil
import glob
import time

from threading import Event, Thread

# Library and function compatibility
if sys.version_info[0] < 3:
    from cStringIO import StringIO
else:
    from io import StringIO

logger = logging.getLogger(__name__)
smiles_regex = re.compile('^([^J][A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+)$')
molmass = {'Ru': 101.072, 'Re': 186.2071, 'Rf': 267.0, 'Rg': 282.0, 'Ra': 226.0, 'Rb': 85.46783, 'Rn': 222.0,
           'Rh': 102.905502, 'Be': 9.01218315, 'Ba': 137.3277, 'Bh': 270.0, 'Bi': 208.980401, 'Bk': 247.0,
           'Br': 79.904, 'Og': 294.0, 'H': 1.008, 'P': 30.9737619985, 'Os': 190.233, 'Es': 252.0, 'Hg': 200.5923,
           'Ge': 72.6308, 'Gd': 157.253, 'Ga': 69.7231, 'Pr': 140.907662, 'Pt': 195.0849, 'Pu': 244.0, 'C': 12.011,
           'Pb': 207.21, 'Pa': 231.035882, 'Pd': 106.421, 'Cd': 112.4144, 'Po': 209.0, 'Pm': 145.0, 'Hs': 269.0,
           'Ho': 164.930332, 'Uue': 315.0, 'Hf': 178.492, 'K': 39.09831, 'He': 4.0026022, 'Md': 258.0, 'Mg': 24.305,
           'Mc': 289.0, 'Mo': 95.951, 'Mn': 54.9380443, 'O': 15.999, 'Mt': 278.0, 'S': 32.06, 'W': 183.841,
           'Zn': 65.382, 'Eu': 151.9641, 'Zr': 91.2242, 'Er': 167.2593, 'Nh': 286.0, 'Ni': 58.69344, 'No': 259.0,
           'Na': 22.989769282, 'Nb': 92.906372, 'Nd': 144.2423, 'Ne': 20.17976, 'Np': 237.0, 'Fr': 223.0,
           'Fe': 55.8452, 'Fl': 289.0, 'Fm': 257.0, 'B': 10.81, 'F': 18.9984031636, 'Sr': 87.621, 'N': 14.007,
           'Kr': 83.7982, 'Si': 28.085, 'Sn': 118.7107, 'Sm': 150.362, 'V': 50.94151, 'Sc': 44.9559085, 'Sb': 121.7601,
           'Sg': 269.0, 'Se': 78.9718, 'Co': 58.9331944, 'Cn': 285.0, 'Cm': 247.0, 'Cl': 35.45, 'Ca': 40.0784,
           'Cf': 251.0, 'Ce': 140.1161, 'Xe': 131.2936, 'Lu': 174.96681, 'Cs': 132.905451966, 'Cr': 51.99616,
           'Cu': 63.5463, 'La': 138.905477, 'Ts': 294.0, 'Li': 6.94, 'Lv': 293.0, 'Tl': 204.38, 'Tm': 168.934222,
           'Lr': 266.0, 'Th': 232.03774, 'Ti': 47.8671, 'Te': 127.603, 'Tb': 158.925352, 'Tc': 98.0, 'Ta': 180.947882,
           'Yb': 173.0451, 'Db': 268.0, 'Dy': 162.5001, 'Ds': 281.0, 'I': 126.904473, 'U': 238.028913, 'Y': 88.905842,
           'Ac': 227.0, 'Ag': 107.86822, 'Ir': 192.2173, 'Am': 243.0, 'Al': 26.98153857, 'As': 74.9215956,
           'Ar': 39.9481, 'Au': 196.9665695, 'At': 210.0, 'In': 114.8181}


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
        # TODO: add timeout to prevent infinite jobs
        was_successfull = True
        self.log.info('Execute cli process: {0}'.format(' '.join(cmd)))
        try:
            process = subprocess.Popen(cmd, cwd=self.workdir,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as err:
            self.log.error('Process failed:', err)
            was_successfull = False
        else:
            self.log.info('Process returncode: {0}'.format(process.returncode))
            output, errors = process.communicate()

            if errors:
                self.log.error(errors.decode('utf-8'))
            if output:
                self.log.info(output.decode('utf-8'))

        return was_successfull


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


def hydrophobic_atom_count(mol2_dict, hphob_types=('C.1', 'C.2', 'C.3', 'C.ar', 'S.3')):
    """
    Calculate the number of hydrophobic atoms in a Tripos MOL2 file based on
    SYBYL atom types in `hphob_types`.

    :param mol2_dict:   Tripos MOL2 atom records as returned by `parse_tripos_atom`
    :type mol2_dict:    :py:dict
    :param hphob_types: hydrophic atoms as SYBYL atom type
    :type hphob_types:  :py:tuple, :py:list

    :return:            hydrophobic atom count
    :rtype:             :py:int
    """

    return sum([1 for atoms in mol2_dict.values() if atoms['atom_type'] in hphob_types])


def molecular_weight(mol2_dict):
    """
    Calculate the molecular weight of a ligand based on the atom elements
    contained in the Tripos MOL2 SYBYL atom types.

    :param mol2_dict:   Tripos MOL2 atom records as returned by `parse_tripos_atom`
    :type mol2_dict:    :py:dict

    :return:            molecular weight
    :rtype:             :py:float
    """

    mol_weight = 0.0
    for atom in mol2_dict.values():
        element = atom['atom_type'].split('.')[0]
        if element in molmass:
            mol_weight += molmass[element]
        else:
            logging.warning('Unknown element: {0}'.format(element))

    return mol_weight


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


def parse_tripos_atom(mol2):
    """
    Parse Tripos MOL2 ATOM records to a dictionary

    :param mol2:    Tripos MOL2 file as string
    :type mol2:     :py:str

    :return:        Tripos atom records
    :rtype:         :py:dict
    """

    read = False
    atom_dict = {}
    headers = ('atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge')
    for line in mol2.split('\n'):

        if line.startswith('@<TRIPOS>ATOM'):
            read = True
            continue

        if line.startswith('@<TRIPOS>BOND'):
            read = False
            break

        if read:
            line = line.split()
            if line:
                atom_dict[int(line[0])] = dict(zip(headers, line[1:]))

    for value in atom_dict.values():
        for key in ('charge', 'x', 'y', 'z'):
            value[key] = float(value[key])
        value['subst_id'] = int(value['subst_id'])

    return atom_dict


def parse_tripos_bond(mol2):
    """
    Parse Tripos MOL2 BOND records to a dictionary

    :param mol2:    Tripos MOL2 file as string
    :type mol2:     :py:str

    :return:        Tripos bond records
    :rtype:         :py:dict
    """

    read = False
    bond_dict = {}
    headers = ('b_start', 'b_end', 'b_type')
    for line in mol2.split('\n'):

        if line.startswith('@<TRIPOS>BOND'):
            read = True
            continue

        if line.startswith('@<TRIPOS>') and read:
            read = False
            break

        if read:
            line = line.split()
            if line:
                bond_dict[int(line[0])] = dict(zip(headers, line[1:]))

    for value in bond_dict.values():
        for key in ('b_start', 'b_end'):
            value[key] = int(value[key])

    return bond_dict


def mol2_to_pdb(mol2, record='ATOM', chain='A'):

    pdb = []
    for atom_id, values in mol2.items():
        atom = [record]
        atom.append(atom_id)
        atom.append(values['atom_name'])
        atom.append(values['subst_name'][0:3])
        atom.append(chain)
        atom.append(values['subst_id'])
        atom.append(values['x'])
        atom.append(values['y'])
        atom.append(values['z'])
        atom.append('1.00')

        pdb.append(atom)

    return pdb


def create_multi_mol2(mol2_file_paths, protein=None):
    """
    Create a multi-molecule MOL2 file by concatenating
    single MOL2 files.

    :param mol2_file_paths: single MOL2 file paths
    :type mol2_file_paths:  :py:list

    :return:                multi-mol2 file
    :rtype:                 :py:str
    """

    multimol2 = StringIO()
    if protein:

        with open(protein, 'r') as singleprotein:
            multimol2.write(singleprotein.read())

    for path in mol2_file_paths:
        if os.path.exists(path):

            with open(path, 'r') as singlemol2:
                multimol2.write(singlemol2.read())

    multimol2.seek(0)
    return multimol2.read()


def create_multi_pdb(mol2_file_paths, protein=None):
    """
    Create a multi-molecule PDB file by converting independent MOL2 files
    to PDB format and concatenating them in a single multi MODEL PDB file

    :param mol2_file_paths: single MOL2 file paths
    :type mol2_file_paths:  :py:list

    :return:                multi MODEL PDB file
    :rtype:                 :py:str
    """

    prot_pdb = None
    if protein:

        with open(protein, 'r') as singleprotein:
            prot_mol2 = parse_tripos_atom(singleprotein.read())
            prot_pdb = mol2_to_pdb(prot_mol2)

    multi_pdb = StringIO()
    for model, path in enumerate(mol2_file_paths, start=1):

        multi_pdb.write('MODEL {0}\n'.format(model))
        if prot_pdb:
            for line in prot_pdb:
                multi_pdb.write('{0:6}{1:>5} {2:^5}{3:>3} {4}{5:>4}    {6:8.3f}{7:8.3f}{8:8.3f}  {9:6}\n'.format(*line))
            multi_pdb.write('TER\n')

        with open(path, 'r') as singlemol2:
            lig_mol2 = parse_tripos_atom(singlemol2.read())
            lig_pdb = mol2_to_pdb(lig_mol2, record='HETATM', chain='B')

            for line in lig_pdb:
                multi_pdb.write('{0:6}{1:>5} {2:^5}{3:>3} {4}{5:>4}    {6:8.3f}{7:8.3f}{8:8.3f}  {9:6}\n'.format(*line))

        multi_pdb.write('ENDMDL\n')
    multi_pdb.write('END\n')

    multi_pdb.seek(0)
    return multi_pdb.read()


def merge_protein_ligand_mol2(protein, ligand, name='system'):
    """
    Merge a protein and ligand structure in Tripos MOL2 format together
    as one structure system.

    :param protein:    Tripos MOL2 file of the protein as string
    :type protein:     :py:str
    :param protein:    Tripos MOL2 file of the ligand as string
    :type protein:     :py:str
    :param name:       new system name
    :type name:        :py:str

    :return:           combined mol2 file
    :rtype:            :py:str
    """

    prot_atom = parse_tripos_atom(protein)
    prot_bond = parse_tripos_bond(protein)
    lig_atom = parse_tripos_atom(ligand)
    lig_bond = parse_tripos_bond(ligand)

    # Stats
    total_atoms = len(prot_atom) + len(lig_atom)
    total_bonds = len(prot_bond) + len(lig_bond)
    id_trans_dict = {}

    merged_mol = StringIO()
    merged_mol.write('@<TRIPOS>MOLECULE\n{0}\n'.format(name))
    merged_mol.write('{0} {1} 1\n'.format(total_atoms, total_bonds))
    merged_mol.write('SMALL\nGASTEIGER\n\n')

    merged_mol.write('@<TRIPOS>ATOM\n')
    for i, a in enumerate(sorted(prot_atom.keys()), start=1):
        merged_mol.write('{0:>7}  {1:8}{2:9.4f} {3:9.4f} {4:9.4f} {5:<5}{6:>4}  {7:8} {8:9.4f}\n'.format(i,
            prot_atom[a]['atom_name'], prot_atom[a]['x'], prot_atom[a]['y'], prot_atom[a]['z'],
            prot_atom[a]['atom_type'], prot_atom[a]['subst_id'], prot_atom[a]['subst_name'], prot_atom[a]['charge']))

    for i, a in enumerate(sorted(lig_atom.keys()), start=i+1):
        id_trans_dict[a] = i
        merged_mol.write('{0:>7}  {1:8}{2:9.4f} {3:9.4f} {4:9.4f} {5:<5}{6:>4}  {7:8} {8:9.4f}\n'.format(i,
            lig_atom[a]['atom_name'], lig_atom[a]['x'], lig_atom[a]['y'], lig_atom[a]['z'], lig_atom[a]['atom_type'],
            lig_atom[a]['subst_id'], lig_atom[a]['subst_name'], lig_atom[a]['charge']))

    merged_mol.write('@<TRIPOS>BOND\n')
    for i, a in enumerate(sorted(prot_bond.keys()), start=1):
        merged_mol.write('{0:>6}{1:>6}{2:>6} {3:>4}\n'.format(i, prot_bond[a]['b_start'], prot_bond[a]['b_end'],
                                                              prot_bond[a]['b_type']))

    for i, a in enumerate(sorted(lig_bond.keys()), start=i+1):
        merged_mol.write('{0:>6}{1:>6}{2:>6} {3:>4}\n'.format(i, id_trans_dict[lig_bond[a]['b_start']],
                                                              id_trans_dict[lig_bond[a]['b_end']],
                                                              lig_bond[a]['b_type']))

    merged_mol.seek(0)
    return merged_mol.read()


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
    docking_dir_name = os.path.basename(result_dir)
    for resultcsv in files:
        resultcsv = os.path.join(result_dir, resultcsv)
        if os.path.isfile(resultcsv):

            header = []
            with open(resultcsv, 'r') as csv_file:
                for line in csv_file.readlines():

                    line = line.strip().split(',')
                    if not header:
                        header = line
                        continue

                    # Only import structure selection if needed
                    mol2 = line[0]
                    path = os.path.join(result_dir, '{0}.mol2'.format(mol2))
                    if structures is not None and path not in structures:
                        continue

                    row = {}
                    for i, val in enumerate(line[1:]):

                        if not len(val):
                            row[header[i]] = None

                        elif '.' in val:
                            row[header[i]] = float(val)

                        else:
                            row[header[i]] = int(val)

                    row['PATH'] = os.path.join(docking_dir_name, '{0}.mol2'.format(mol2))
                    results[mol2] = row
            break

    return results


class PeriodicCleanup(object):
    """
    Asynchronous periodic cleanup class checking a base working directory for
    result directories starting with 'docking-' that are more then
    `result_storage_time` hours old every `period` seconds and removing them.
    """

    def __init__(self, base_work_dir, result_storage_time, period=60):
        """

        :param base_work_dir:        base directory containing 'docking-' dirs.
        :type base_work_dir:         :py:str
        :param result_storage_time:  maximum time in hours results are stored
                                     on disc.
        :type result_storage_time:   :py:int
        :param period:               event interval time in seconds
        :type period:                :py:int
        """

        self.base_work_dir = base_work_dir
        self.result_storage_time = result_storage_time * 3600
        self.period = period

        logging.info('Start periodic cleanup of "{0}" every {1} sec. removing results > {2} hours old'.format(
            base_work_dir, period, result_storage_time))

        self.clean_event = Event()

    def start(self):
        """
        Start the Event thread
        """

        self.clean_event.clear()
        self.proc = Thread(target=self.cleanup)
        self.proc.start()

    def stop(self):
        """
        Nearly immediately kills the Periodic function
        """

        self.clean_event.set()
        self.proc.join()

    def cleanup(self):

        while True:
            self.clean_event.wait(self.period)
            if self.clean_event.is_set():
                break

            for dockdir in glob.glob('{0}/docking-*/'.format(self.base_work_dir)):

                if os.path.isdir(dockdir):
                    if int(time.time()) - int(os.path.getmtime(dockdir)) >= self.result_storage_time:

                        logging.info('Periodic cleanup, remove: {0}'.format(dockdir))
                        shutil.rmtree(dockdir)

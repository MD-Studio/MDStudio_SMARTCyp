# -*- coding: utf-8 -*-

"""
Classes for protein-ligand docking using the PLANTS software:

  PLANTS: Protein-Ligand ANT System.
  "An ant colony optimization approach to flexible protein-ligand docking" Swarm Intell. 1, 115-134 (2007).
  O. Korb, T. Stützle, T.E. Exner
  URL: http://www.tcd.uni-konstanz.de/research/plants.php
"""

import logging
import os
import json
import copy
import glob

from mdstudio_smartcyp import __module__, __package_path__, __plants_path__, __plants_version__, __plants_citation__
from mdstudio_smartcyp.plants_conf import PLANTS_CONF_FILE_TEMPLATE
from mdstudio_smartcyp.utils import (_schema_to_data, RunnerBaseClass, prepare_work_dir, create_multi_mol2,
                                     create_multi_pdb, import_plants_csv, atom_count, MDStudioException)
from mdstudio_smartcyp.clustering import coords_from_mol2, ClusterStructures

logger = logging.getLogger(__module__)

PLANTS_DOCKING_SCHEMA = os.path.join(__package_path__, 'schemas/endpoints/docking_request.v1.json')
settings = _schema_to_data(json.load(open(PLANTS_DOCKING_SCHEMA)))


def plants_version_info():
    """
    :return: information on the packaged PLANTS version, supported models
             default configuration and software citation reference.
    :rtype:  :py:dict
    """

    default_settings = copy.deepcopy(settings)
    for param in ('base_work_dir', 'bindingsite_center', 'output_dir'):
        del default_settings[param]

    info_dict = {'version': __plants_version__,
                 'citation': __plants_citation__,
                 'default_settings': default_settings}

    return info_dict


class PlantsDocking(RunnerBaseClass):
    """
    Class for running a protein-ligand docking using the PLANTS docking
    software:

      PLANTS: Protein-Ligand ANT System.
      "An ant colony optimization approach to flexible protein-ligand
      docking" Swarm Intell. 1, 115-134 (2007). O. Korb, T. Stützle, T.E. Exner
      URL: http://www.tcd.uni-konstanz.de/research/plants.php

    This class is compatible with PLANTS versions 1.1 and 1.2.
    If not otherwise defined, the PLANTS executable files are available
    in the bin directory of the lie_plants_docking package suffixed by the
    OS identifier as returned by `sys.platform`.
    Support is available for all of PLANTS default configuration options
    described in sections 1.0 of the PLANTS manual.

    Run a PLANTS docking as:
    ::
        docking = PlantsDocking(plants_config_dict)
        docking.run(protein, ligand)
        results_json = docking.results()

    :param log:           Python logger instance
    :type log:            :py:logging
    :param base_work_dir: base directory for unique docking results dirs.
    :type base_work_dir:  :py:str
    :param kwargs:        additional keyword arguments are considered as
                          PLANTS configuration options
    :type kwargs:         :py:dict
    """

    def __init__(self, log=logger, base_work_dir=None, **kwargs):

        self.log = log
        self.base_work_dir = base_work_dir
        self._workdir = None

        self.config = copy.deepcopy(settings)
        self.config['exec_path'] = __plants_path__
        self.config.update(kwargs)

    def __setitem__(self, key, value):
        """
        __setitem__ overload.

        Set values using dictionary style access, fallback to
        default __setattr__

        :param key:   attribute name
        :type key:    str
        :param value: attribute value
        """

        if key in self.config or key in settings:
            self.config[key] = value
        else:
            dict.__setattr__(self, key, value)

    def _absolute_paths(self, rel_paths=None):
        """
        Construct absolute paths based relative 'docking workdir/pose name'
        paths

        :param rel_paths:   docking pose structure path IDs for which to
                            return results. Defaults to all poses.
        :type rel_paths:    :py:list

        :return:            absolute paths
        :rtype:             :py:list
        """

        if rel_paths is None:
            if self.workdir:
                rel_paths = [path for path in glob.glob(os.path.join(self.workdir, '*_entry_*_conf_*.mol2'))]

        if isinstance(rel_paths, str):
            rel_paths = [rel_paths]

        if isinstance(rel_paths, (list, tuple)):
            rel_paths = [os.path.join(self.base_work_dir, path) if not os.path.exists(path) else path for
                         path in rel_paths]

            # Cannot combine results from different docking runs
            if len(set([os.path.dirname(path) for path in rel_paths])) > 1:
                raise MDStudioException('Unable to combine results for more then one docking run')

            # Do the results still exist
            if any([not os.path.exists(path) for path in rel_paths]):
                raise MDStudioException('Docking results (no longer) exist: {0}'.format(os.path.basename(rel_paths[0])))

        if rel_paths:

            # Set the docking working directory
            self.workdir = os.path.dirname(rel_paths[0])

        return rel_paths

    @property
    def workdir(self):
        """
        Get the validated PLANTS working directory as absolute path
        """

        return self._workdir

    @workdir.setter
    def workdir(self, wdir):
        """
        Set the PLANTS working directory. This may be a new empty directory
        ready for a new PLANTS docking run or an existing directory for an
        analysis. In the latter case, the existence of the directory is
        checked as it may have been deleted already by the temporary cleanup
        routines of the service.

        :param wdir: working directory path
        :type wdir:  :py:str
        """

        if isinstance(wdir, str):

            if self.base_work_dir and self.base_work_dir not in wdir:
                wdir = os.path.join(self.base_work_dir, wdir)

            wdir = os.path.abspath(wdir)
            if not os.path.isdir(wdir):
                raise MDStudioException('Docking results (no longer) exist: {0}'.format(wdir))
            self.config['workdir'] = wdir
            self._workdir = wdir
        else:
            self._workdir = None

    def update(self, config=None, **kwargs):
        """
        Update the configuration settings for the docking method from a
        dictionary of custom settings.

        Configuration options (keys) are validated against the set of
        allowed (known) option names for the docking method.

        :param config: custom settings
        :type config:  :py:class:`dict`
        """

        if not config:
            config = kwargs

        for key, value in config.items():
            if key in settings:
                self.config[key] = value
            else:
                self.log.warn('PLANTS configuration file has no setting named: {0}'.format(key))
        self.log.info('Override PLANTS configuration for options: {0}'.format(', '.join(config.keys())))

    def get_results(self, structures=None, do_cluster=True):
        """
        Return PLANTS results

        PLANTS general docking results are stored in the features.csv
        and ranking.csv Comma Separated Value files.
        Ranking is a subset of features. The latter contains additional
        scoring function specific data and is parsed in favour of the
        first.

        Results are parsed into a dictionary with the docking pose
        identifier as key. These identifiers are already sorted by
        PLANTS docking score.

        :param structures:      docking pose structure path IDs for which to
                                return results. Defaults to all poses.
        :type structures:       :py:list

        :return:                general PLANTS docking results
        :rtype:                 :py:dict
        """

        # Structure selection to return results for
        structures = self._absolute_paths(structures)

        # Read docking results: first try features.csv, else ranking.csv
        results = import_plants_csv(self.workdir, structures)

        if do_cluster:

            # Run a clustering
            xyz = coords_from_mol2(structures)
            try:
                c = ClusterStructures(xyz, labels=list(results.keys()))
            except AssertionError as e:
                logging.error(e)
                return None

            clusters = c.cluster(threshold=self.config.get('threshold', 8.0),
                                 criterion=self.config.get('criterion', 'maxclust'),
                                 min_cluster_count=self.config.get('min_cluster_size', 2))

            for structure, res in clusters.items():
                results[structure].update(res)

            # Plot cluster results
            clusterplot = os.path.join(self.workdir, 'cluster_dendrogram.pdf')
            c.plot(to_file=clusterplot)

        return results

    def get_structures(self, structures, output_format='mol2', include_protein=False, create_ensemble=True):
        """
        Return docked poses in MOL2 or PDB file format. If `create_ensemble`
        all poses will be concatenated as single multi-molecule (ensemble)
        structure file.

        The structures to return are identified by their unique path ID
        listed in the docking results as the 'PATH' variable. It consists of
        the '<docking results directory name>/<pose name>' that retained on
        the server for the duration defined by 'result_storage_time'.

        :param structures:      docking pose structure path IDs.
        :type structures:       :py:list
        :param output_format:   structure file format as MOL2 or PDB
        :type output_format:    :py:str
        :param include_protein: return structures including the protein
                                or only the ligand structures
        :type include_protein:  :py:bool
        :param create_ensemble: concatenate multiple docked poses as a single
                                multi-molecule ensemble file.

        :return:                list of multiple molecular structures or one
                                single (ensemble) structure in MOL2 or PDB
                                format.
        :rtype:                 :py:str, :py:list
        """

        structures = self._absolute_paths(structures)

        protein = None
        if include_protein:
            protein = os.path.join(self.workdir, 'protein.mol2')

        self.log.debug('Return {0} structures for {1}'.format(len(structures), self.workdir))

        # Return structure ensemble in PDB or MOL2 format
        if output_format == 'mol2':
            if create_ensemble:
                return create_multi_mol2(structures, protein=protein)
            return [create_multi_mol2([mol], protein=protein) for mol in structures]
        elif output_format == 'pdb':
            if create_ensemble:
                return create_multi_pdb(structures, protein=protein)
            return [create_multi_pdb([mol], protein=protein) for mol in structures]

    def run(self, protein, ligand, mode='screen'):
        """
        Run a PLANTS docking for a given protein and ligand in mol2
        format in either 'screen' or 'rescore' mode.

        A docking run requires the following PLANTS configuration arguments
        to be defined:
        * exec_path: path to the PLANTS executable
        * workdir: a working directory to write docking results to
        * bindingsite_center: target ligand binding site in the protein defined
          as a 3D coordinate
        The `run` function will exit if any of these requirements are not
        resolved.

        The PLANTS_CONF_FILE_TEMPLATE serves as a template where
        option values are replaced by format placeholders with the
        same name as the keys in the configuration dictionary.

        .. note:: the PLANTS write_multi_mol2 parameter is turned off by
                  default to enable seperate clustering and result retrieval
                  by the user.

        :param protein: protein 3D structure in mol2 format
        :type protein:  str
        :param ligand:  ligand 3D structure in mol2 format
        :type ligand:   str
        :param mode:    PLANTS execution mode as either virtual
                        screening 'screen' or rescoring 'rescore'
        :type mode:     str

        :return:        boolean to indicate successful docking
        :rtype:         bool
        """

        # Check required PLANTS configuration arguments
        exec_path = self.config.get('exec_path')
        if not os.path.exists(exec_path):
            raise MDStudioException('Plants executable not available at: {0}'.format(exec_path))

        if not os.access(exec_path, os.X_OK):
            raise MDStudioException('Plants executable {0} does not have exacutable permissions'.format(exec_path))

        if sum(self.config['bindingsite_center']) == 0 or len(self.config['bindingsite_center']) != 3:
            raise MDStudioException('Malformed binding site center definition: {0}'.format(
                self.config['bindingsite_center']))

        # Create a working directory
        self.workdir = prepare_work_dir(path=self.base_work_dir, prefix='docking-')
        self.log.info('Created docking directory {0}'.format(self.workdir))

        # Copy files to working directory
        if os.path.isfile(protein):
            self.config['protein_file'] = protein
        else:
            with open(os.path.join(self.workdir, 'protein.mol2'), 'w') as protein_file:
                protein_file.write(protein)
                self.config['protein_file'] = 'protein.mol2'

        if os.path.isfile(ligand):
            self.config['ligand_file'] = ligand
        else:
            with open(os.path.join(self.workdir, 'ligand.mol2'), 'w') as ligand_file:
                ligand_file.write(ligand)
                self.config['ligand_file'] = 'ligand.mol2'

        lig_atm_count = atom_count(os.path.join(self.workdir, self.config['ligand_file']))
        prot_atm_count = atom_count(os.path.join(self.workdir, self.config['protein_file']))
        if lig_atm_count > prot_atm_count:
            self.delete()
            raise MDStudioException('Ligand structure contains more atoms than protein structure. Swapped input?')

        # Write PLANTS configuration file
        conf_file = os.path.join(self.workdir, 'plants.config')
        with open(conf_file, 'w') as conf:
            conf.write(PLANTS_CONF_FILE_TEMPLATE.format(**self.config))

        success = self.cmd_runner([exec_path, '--mode', mode, 'plants.config'])
        if not success or not len(glob.glob(os.path.join(self.workdir, '*_entry_*_conf_*.mol2'))):
            success = False
            self.delete()

        return success

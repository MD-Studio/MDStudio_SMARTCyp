# -*- coding: utf-8 -*-

import itertools
import logging
import numpy
import matplotlib

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

# Init maplotlib
from matplotlib import style
matplotlib.use('Agg')  # Use Agg for non-interactive plotting
style.use('ggplot')  # Because of AttributeError: Unknown property color_cycle bug in Pandas 1.7.1 with Matplotlib 1.5.0

import matplotlib.pyplot as plt


def coords_from_mol2(mol2_files):
    """
    Extract XYZ coordinates from a mol2 file

    :param mol2_files: mol2 file paths to import
    :type mol2_files:  list
    :return:           coordinates as list of numpy arrays
    """

    sets = []
    for mol2 in mol2_files:

        with open(mol2, 'r') as structure_file:
            read = False
            coords = []
            for line in structure_file.readlines():
                line = line.strip()
                if line.endswith('ATOM'):
                    read = True
                elif line.endswith('BOND'):
                    break
                elif read:
                    line = line.split()
                    coords.append(list(map(float, line[2:5])))
            if coords:
                sets.append(coords)

    return numpy.array(sets)


def _rotate(c1, c2):
    """
    Rotate matrix c1 onto matrix c2 using Kabsch algorithm
    """
    rotated = _kabsch(c1, c2)

    # Rotate c1
    return numpy.dot(c1, rotated)


def _centroid(x):
    """
    Calculate the centroid from a vectorset x
    """

    return sum(x) / len(x)


def _kabsch(pmat, qmat):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    pmat unto matrix qmat so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point pmat and qmat,
    centered around the center-of-mass. Each vector set is represented as
    an NxD matrix, where D is the the dimension of the space.

    The algorithm works in three steps:
        * A translation of pmat and qmat
        * The computation of a covariance matrix C
        * Computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    :param pmat: (N, number of points)x(D, dimension) matrix
    :param qmat: (N, number of points)x(D, dimension) matrix

    :return:     the rotation matrix
    :rtype:      object
    """

    # Computation of the covariance matrix
    covmat = numpy.dot(numpy.transpose(pmat), qmat)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    v, s, w = numpy.linalg.svd(covmat)
    d = (numpy.linalg.det(v) * numpy.linalg.det(w)) < 0.0

    if d:
        s[-1] = -s[-1]
        v[:, -1] = -v[:, -1]

    # Create Rotation matrix U
    return numpy.dot(v, w)


def rmsd(c1, c2):

    delta = c1 - c2
    delta_sq = delta * delta
    return numpy.sqrt(delta_sq.sum() / len(c1))


def kabsch(c1, c2):
    """
    Rotate matrix c1 onto c2 and calculate the RMSD
    """

    c1 -= _centroid(c1)
    c2 -= _centroid(c2)

    c1 = _rotate(c1, c2)
    return rmsd(c1, c2)


class ClusterStructures(object):
    """
    Cluster analysis on sets of structures that are identical in atom count,
    atom type and atom order such as docking poses.

    The class works on plain sets of coordinates and thereby assumes that
    the above mentioned atom constraints are satisfied.
    A calculated pairwise distance matrix (**pdist**) is used as similarity measure
    in a hierarchical clustering analysis.

    Supported metrics of the **pdist** are:

        * **rmsd**:
          plane Root Mean Square Deviation between two coordinate
          sets without fitting.
        * **kabsch**:
          first fit the two coordinates set before calculating the
          RMSD value. Minimize the Root Mean Square Deviation
          between the two coordinate sets by calculating the
          optimal rotation matrix.

          This uses the Kabsch algorithm:

          * https://en.wikipedia.org/wiki/Kabsch_algorithm

    The **pdist** is calculated once at class construction and is reused to
    calculate different cluster flavors using the class `cluster` method.

    :param xyz:    structure xyz coordinate sets
    :type xyz:     list of numpy.ndarray's
    :param metric: one of the supported distance metrics
    :type metric:  str
    :param labels: structure identifiers such as docking pose ID's
                   corresponding to the entries in the coordinate set.
    :type labels:  list
    """

    logger = logging.getLogger(__name__)

    def __init__(self, xyz, metric='rmsd', labels=None):

        self.xyz = xyz
        self.metric = metric
        self.labels = labels or range(len(labels))

        # All xyz coordinate sets need to be of type numpy.ndarray
        if not all([isinstance(coords, numpy.ndarray) for coords in self.xyz]):
            raise AssertionError('Structure coordinates need to be of type numpy.ndarray')

        # Equality in number and order of atoms for all coordinate sets is assumed
        if not len(set([coords.size for coords in self.xyz])) == 1:
            raise AssertionError('Structure coordinates have an unequal number of atoms')

        # Optional list of labels (e.a. structure id's) need to match coordinate set in length
        if not len(self.labels) == len(self.xyz):
            raise AssertionError('Number of labels is not matching number of coordinate sets')

        self._condensed_distance_matrix = self._build_pdist()
        self._clusters = []
        self._clusters_filtered = {}

    def __len__(self):
        """
        Return the number of clustered structures
        """

        len(self._clusters)

    def __str__(self):
        """
        Return a nicely formatted summary of the clustering
        """

        summary = [
            'Clustering {0} structures'.format(len(self.xyz)),
            'Metric for pairwise distance matrix: {0}'.format(self.metric),
            'Metric for hierarchical clustering: {0}'.format(self.method),
            'Cluster selection criterion: {0} with parameter {1}\n'.format(self.criterion, self.threshold),
            'Clusters: {0}, coverage: {1:.2f}%'.format(self.cluster_count, self.coverage * 100),
        ]

        return '\n'.join(summary)

    def _build_pdist(self):
        """
        Construct condensed pairwise distance matrix using the `metric`
        defined in the class constructor

        :return:       condensed distance matrix
        :rtype:        list
        """

        _metric_func = globals().get(self.metric)
        if not _metric_func:
            raise LookupError(
                '{0} class does not know about "{1}" pdist metric'.format(type(self).__name__, self.metric))

        dst = []
        for pair in itertools.combinations(self.xyz, 2):
            dst.append(_metric_func(*pair))

        return dst

    @property
    def clustered_structures(self):
        """
        For all structures that have been assigned to a cluster,
        return the cluster ID's

        :rtype: list
        """

        return [sid for sid in self._clusters_filtered if self._clusters_filtered[sid].get('CLUSTER', 0) != 0]

    @property
    def cluster_medians(self):
        """
        Return ID's for structures that are representative for the
        cluster center determined as the structure for which the
        similarity distance metric is closest to the average of the
        given cluster.

        :return: dictionary of cluster ID/structure ID key/value pairs
        :rtype:  :py:class:`dict`
        """

        return {self._clusters_filtered[sid]['CLUSTER']: sid for sid
                in self._clusters_filtered
                if self._clusters_filtered[sid].get('MEAN', False)}

    @property
    def cluster_count(self):

        return len(self.cluster_medians)

    @property
    def coverage(self):
        """
        Cluster coverage as the fraction of structures assigned to a cluster

        :rtype: float
        """

        return len(self.clustered_structures) / float(len(self.xyz))

    def plot(self, to_file='cluster_dendrogram.pdf'):
        """
        Plot the cluster dendrogram.
        """

        fig = plt.figure()

        ddata = dendrogram(self._linkage,
                           labels=self.labels,
                           color_threshold=self._linkage[-(max(self._clusters) - 1), 2])
        annotate_above = 1
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        max_d = None
        if max_d:
            fig.axhline(y=max_d, c='k')

        if not to_file:
            return fig

        fig.savefig(to_file)

    def cluster(self, threshold=5, method='single', criterion='maxclust', min_cluster_count=1):
        """
        Cluster the structures using hierarchical clustering methods on
        the calculated pairwise distance matrix.

        The 'maxclust' criterion is used by default to determined the optimal number
        of clusters between 1 and 'threshold' that will assign a maximum number of structures
        to each cluster

        :param threshold:         pairwise distance cutoff used to populate returned
                                  clusters. Same as in scipy.cluster.hierarchy.fcluster
        :type threshold:          :py:float
        :param method:            hierarchical clustering methods as defined in
                                  the scipy.cluster.hierarchy.linkage method.
                                  Options are: single, complete, average, weighted,
                                  centroid, median and ward.
        :type method:             str
        :param criterion:         method to use for flattening clusters from the
                                  hierarchical clustering as defined in the
                                  scipy.cluster.hierarchy.fcluster method.
                                  Options are: inconsistent, distance, maxclust
                                  monocrit and maxclust_monocrit. 'Threshold'
                                  is used as threshold for each of these methods
        :type criterion:          str
        :param min_cluster_count: minimal number of structures in a cluster
        :type min_cluster_count:  int

        :return:                  clustering results
        :rtype:                   :py:class:`dict`
        """

        self.method = method
        self.criterion = criterion
        self.threshold = threshold

        self._linkage = linkage(self._condensed_distance_matrix, method=method)
        self._clusters = fcluster(self._linkage, threshold, criterion=criterion)

        self._clusters_filtered = {}
        sqmatr = squareform(self._condensed_distance_matrix)
        for n in range(1, max(self._clusters) + 1):
            cl = numpy.where(self._clusters == n)
            if len(cl[0]) >= min_cluster_count:

                # Get one structure as representative of the cluster
                if len(cl[0]) == 1:
                    lc = [0]
                else:
                    a = sqmatr[cl[0]]
                    a = a[:, cl[0]]
                    m = numpy.mean(numpy.hstack([a[i][:i] for i in range(a.shape[0])]))
                    x = (numpy.abs(a - m)).argmin()
                    lc = numpy.where(a == a.flat[x])[0]

                meanpose = self.labels[cl[0][lc[0]]]
                for idx in cl[0]:
                    self._clusters_filtered[self.labels[idx]] = {'CLUSTER': n, 'MEAN': self.labels[idx] == meanpose}
            else:
                self.logger.debug('Cluster {0} contains less that {1} structures ({2}). Dropping'.format(n,
                                                                                        min_cluster_count, len(cl[0])))
                for idx in cl[0]:
                    self._clusters_filtered[self.labels[idx]] = {'CLUSTER': 0, 'MEAN': 0}

        self.logger.info('Cluster {0} structures. pdist method: {1}, cluster method: {2}, criterion: {3}, tolerance: '
                         '{4}, minimum cluster size: {5}'.format(len(self.xyz), self.metric, self.method,
                                                                 self.criterion, self.threshold, min_cluster_count))
        self.logger.info('Resolved {0} clusters, coverage of {1}%'.format(self.cluster_count, self.coverage * 100))

        return self._clusters_filtered

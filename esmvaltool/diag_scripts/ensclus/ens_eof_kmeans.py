"""Find the most representative ensemble member for each cluster."""

import collections
import datetime
import math
import os

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

# User-defined libraries
from eof_tool import eof_computation
from read_netcdf import read_n_2d_fields


def ens_eof_kmeans(dir_output, name_outputs, numens, numpcs, perc, numclus):
    """Find the most representative ensemble member for each cluster.

    METHODS:
    - Empirical Orthogonal Function (EOF) analysis of the input file
    - K-means cluster analysis applied to the retained
      Principal Components (PCs)
    OUTPUT:
    Frequency
    """
    print('The name of the output files will be <variable>_{0}.txt'
          .format(name_outputs))
    print('Number of ensemble members: {0}'.format(numens))

    outfiles = []
    model = name_outputs.split("_")[1]
    print('Model: {0}'.format(model))
    # Either perc (cluster analysis is applied on a number of PCs
    # such as they explain 'perc' of total variance) or numpcs
    # (number of PCs to retain) is set.
    # numpcs has priority over perc, ignored if it is set to 0
    if numpcs:
        numpcs = int(numpcs)
        print('Number of principal components: {0}'.format(numpcs))
    else:
        print('Percentage of variance explained: {0}'.format(perc))

    print('Number of clusters: {0}'.format(numclus))

    # Reading the netCDF file of N 2Dfields of anomalies, saved by ens_anom.py
    ifile = os.path.join(dir_output, 'ens_anomalies_{0}.nc'
                         .format(name_outputs))
    var, _, lat, _ = read_n_2d_fields(ifile)
    print('var dim: (numens x lat x lon)={0}'.format(var.shape))

    # Compute EOFs (Empirical Orthogonal Functions)
    # and PCs (Principal Components) with respect to ensemble memeber
    print('_________________________________________________________')
    print('EOF analysis:')
    # --------------------------------------------------------------------
    _, _, _, pcs_unscal0, eofs_unscal0, varfrac = eof_computation(var, lat)

    acc = np.cumsum(varfrac * 100)
    if numpcs:
        exctperc = acc[numpcs - 1]
    else:
        # Find how many PCs explain a certain percentage of variance
        # (find the mode relative to the percentage closest to perc,
        #  but bigger than perc)
        numpcs = min(enumerate(acc), key=lambda x: x[1] <= perc)[0] + 1
        print('\nThe number of PCs that explain the percentage closest '
              'to {0}% of variance (but grater than {0}%) is {1}'
              .format(perc, numpcs))
        exctperc = min(enumerate(acc), key=lambda x: x[1] <= perc)[1]
    print('(the first {0} PCs explain exactly the {1}% of variance)'
          .format(numpcs, "%.2f" % exctperc))

    # ____________Compute k-means analysis using a subset of PCs
    print('_________________________________________________________')
    print('k-means analysis using a subset of PCs')

    pcs = pcs_unscal0[:, :numpcs]

    clus = KMeans(n_clusters=numclus, n_init=2000,
                  init='k-means++', tol=1e-4,
                  max_iter=1000, random_state=42)
    start = datetime.datetime.now()
    clus.fit(pcs)
    end = datetime.datetime.now()
    print('k-means algorithm took me %s seconds' % (end - start))

    centroids = clus.cluster_centers_          # shape---> (numclus,numpcs)
    labels = clus.labels_                      # shape---> (numens,)

    print('\nClusters are identified for {0} PCs (explained variance {1}%)'
          .format(numpcs, "%.2f" % exctperc))
    print('PCs dim: (number of ensemble members, number of PCs)={0}, '
          'EOF dim: (number of ensemble members, lat, lon)={1}'
          .format(pcs_unscal0[:, :numpcs].shape, eofs_unscal0[:numpcs].shape))
    print('Centroid coordinates dim: (number of clusters, number of PCs)={0}, '
          'labels dim: (number of ensemble members,)={1}\n'
          .format(centroids.shape, labels.shape))

    # ____________Save labels
    namef = os.path.join(dir_output, 'labels_{0}.txt'.format(name_outputs))
    np.savetxt(namef, labels, fmt='%d')
    outfiles.append(namef)

    # ____________Compute cluster frequencies
    clusters = []
    for nclus in range(numclus):
        clu = list(np.where(labels == nclus)[0])
        frq = len(clu) * 100 / len(labels)
        clusters.append([nclus, frq, clu])
    print('Cluster labels:')
    print([clusters[ncl][0] for ncl in range(numclus)])
    print('Cluster frequencies (%):')
    print([round(clusters[ncl][1], 3) for ncl in range(numclus)])
    print('Cluster members:')
    print([clusters[ncl][2] for ncl in range(numclus)])

    # ____________Find the most representative ensemble member for each cluster
    print('_________________________________________________________')
    print('In order to find the most representative ensemble member for each '
          'cluster\n(which is the closest member to the cluster centroid)')
    print('the Euclidean distance between cluster centroids and each ensemble '
          'member is computed in the PC space')
    print('_________________________________________________________')
    print('Check: cluster #1 centroid coordinates vector dim {0} should be '
          'the same as the member #1 PC vector dim {1}\n'
          .format(centroids[1, :].shape, pcs[1, :].shape))

    print('_________________________________________________________')
    print('In order to study the spread of each cluster,')
    print('the standard deviation of the distances between each member '
          'in a cluster and the cluster centroid is computed in the PC space')
    stat_output = []
    repres = []
    for nclus in range(numclus):
        members = clusters[nclus][2]
        norm = np.empty([numclus, len(members)])
        for mem, ens in enumerate(members):
            normens = centroids[nclus, :] - pcs[ens, :]
            norm[nclus, mem] = math.sqrt(sum(normens**2))
        print('the distances between centroid of cluster {0} and its '
              'belonging members {1} are:\n{2}'
              .format(nclus, members, np.round(norm[nclus], 3)))
        print('MINIMUM DISTANCE WITHIN CLUSTER {0} IS {1} --> member #{2}'
              .format(nclus, round(norm[nclus].min(), 3),
                      members[np.where(norm[nclus] ==
                                       norm[nclus].min())[0][0]]))
        repres.append(members[np.where(norm[nclus] ==
                                       norm[nclus].min())[0][0]])
        print('MAXIMUM DISTANCE WITHIN CLUSTER {0} IS {1} --> member #{2}'
              .format(nclus, round(norm[nclus].max(), 3),
                      members[np.where(norm[nclus] ==
                                       norm[nclus].max())[0][0]]))
        print('INTRA-CLUSTER STANDARD DEVIATION FOR CLUSTER {0} IS {1}\n'
              .format(nclus, norm[nclus].std()))

        d_stat = collections.OrderedDict()
        d_stat['cluster'] = nclus
        d_stat['member'] = members
        d_stat['d_to_centroid'] = np.round(norm[nclus], 3)
        d_stat['intra-clus_std'] = norm[nclus].std()
        d_stat['d_min'] = round(norm[nclus].min(), 3)
        d_stat['d_max'] = round(norm[nclus].max(), 3)
        d_stat['freq(%)'] = round(clusters[nclus][1], 3)
        stat = pd.DataFrame(d_stat)
        stat_output.append(stat)

    # ____________Save the most representative ensemble members
    namef = os.path.join(dir_output, 'repr_ens_{0}.txt'.format(name_outputs))
    np.savetxt(namef, repres, fmt='%i')
    outfiles.append(namef)

    stat_output = pd.concat(stat_output, axis=0)
    # ____________Save statistics of cluster analysis
    namef = os.path.join(dir_output, 'statistics_clustering_{0}.txt'
                         .format(name_outputs))
    outfiles.append(namef)
    with open(namef, 'w') as text_file:
        text_file.write(stat_output.__repr__())

    return outfiles

#*********************************
#        ens_eof_kmeans          *
#*********************************

# Standard packages
import numpy as np
import sys
import os
from sklearn.cluster import KMeans
import datetime
import math
import pandas as pd
import collections

def ens_eof_kmeans(dir_OUTPUT,name_outputs,numens,numpcs,perc,numclus):
    '''
    \nGOAL:
    Find the most representative ensemble member for each cluster.
    METHODS:
    - Empirical Orthogonal Function (EOF) analysis of the input file
    - K-means cluster analysis applied to the retained Principal Components (PCs)
    OUTPUT:
    Frequency 
    '''

    # User-defined libraries
    from read_netcdf import read_N_2Dfields
    from eof_tool import eof_computation
    
    print('***********************************OUTPUT***********************************')
    print('The name of the output files will be <variable>_{0}.txt'.format(name_outputs))
    print('Number of ensemble members: {0}'.format(numens))
    
    model=name_outputs.split("_")[1]
    print('Model: {0}'.format(model))
    # Either perc (cluster analysis is applied on a number of PCs such as they explain
    # 'perc' of total variance) or numpcs (number of PCs to retain) is set:
    if numpcs!='no':
        numpcs=int(numpcs)
        print('Number of principal components: {0}'.format(numpcs))
    
    if perc!='no':
        perc=int(perc)
        print('Percentage of explained variance: {0}%'.format(perc))
    
    if (perc=='no' and numpcs=='no') or (perc!='no' and numpcs!='no'):
        raise ValueError('You have to specify either "perc" or "numpcs".')
    
    print('Number of clusters: {0}'.format(numclus))

    #____________Reading the netCDF file of N 2Dfields of anomalies, saved by ens_anom.py
    ifile=os.path.join(dir_OUTPUT,'ens_anomalies_{0}.nc'.format(name_outputs))
    var, varunits, lat, lon = read_N_2Dfields(ifile)
    print('var dim: (numens x lat x lon)={0}'.format(var.shape))
    
    
    #____________Compute EOFs (Empirical Orthogonal Functions)
    #____________and PCs (Principal Components) with respect to ensemble memeber
    print('____________________________________________________________________________________________________________________')
    print('EOF analysis')
    #----------------------------------------------------------------------------------------
    solver, pcs_scal1, eofs_scal2, pcs_unscal0, eofs_unscal0, varfrac = eof_computation(var,varunits,lat,lon)
        
    acc=np.cumsum(varfrac*100)
    if perc!='no':
        # Find how many PCs explain a certain percentage of variance
        # (find the mode relative to the percentage closest to perc, but bigger than perc)
        numpcs=min(enumerate(acc), key=lambda x: x[1]<=perc)[0]+1
        print('\nThe number of PCs that explain the percentage closest to {0}% of variance (but grater than {0}%) is {1}'.format(perc,numpcs))
        exctperc=min(enumerate(acc), key=lambda x: x[1]<=perc)[1]
    if numpcs!='no':
        exctperc=acc[numpcs-1]
    print('(the first {0} PCs explain exactly the {1}% of variance)'.format(numpcs,"%.2f" %exctperc))
    

    #____________Compute k-means analysis using a subset of PCs
    print('____________________________________________________________________________________________________________________')
    print('k-means analysis using a subset of PCs')
    print('____________________________________________________________________________________________________________________')
    #----------------------------------------------------------------------------------------
    PCs=pcs_unscal0[:,:numpcs]

    clus=KMeans(n_clusters=numclus, n_init=600, max_iter=1000)
    
    start = datetime.datetime.now()
    clus.fit(PCs)
    end = datetime.datetime.now()
    print('k-means algorithm took me %s seconds' %(end-start))
    
    centroids=clus.cluster_centers_          # shape---> (numclus,numpcs)
    labels=clus.labels_                      # shape---> (numens,)
    
    print('\nClusters are identified for {0} PCs (explained variance {1}%)'.format(numpcs, "%.2f" %exctperc))
    print('PCs dim: (number of ensemble members, number of PCs)={0}, EOF dim: (number of ensemble members, lat, lon)={1}'.format(pcs_unscal0[:,:numpcs].shape,eofs_unscal0[:numpcs].shape))
    print('Centroid coordinates dim: (number of clusters, number of PCs)={0}, labels dim: (number of ensemble members,)={1}\n'.format(centroids.shape,labels.shape))
    
    #____________Save labels
    namef=os.path.join(dir_OUTPUT,'labels_{0}.txt'.format(name_outputs))
    np.savetxt(namef,labels,fmt='%d')
    
    #____________Compute cluster frequencies
    L=[]
    for nclus in range(numclus):
        cl=list(np.where(labels==nclus)[0])
        fr=len(cl)*100/len(labels)
        L.append([nclus,fr,cl])
    print('Cluster labels:')
    print([L[ncl][0] for ncl in range(numclus)])
    print('Cluster frequencies (%):')
    print([round(L[ncl][1],3) for ncl in range(numclus)]) 
    print('Cluster members:')
    print([L[ncl][2] for ncl in range(numclus)])
    
    #____________Find the most representative ensemble member for each cluster
    print('____________________________________________________________________________________________________________________')
    print('In order to find the most representative ensemble member for each cluster\n(which is the closest member to the cluster centroid)')
    print('the Euclidean distance between cluster centroids and each ensemble member is computed in the PC space')
    print('____________________________________________________________________________________________________________________')
    # 1)
    print('Check: cluster #1 centroid coordinates vector dim {0} should be the same as the member #1 PC vector dim {1}\n'.format(centroids[1,:].shape,PCs[1,:].shape))
    #print('\nIn the PC space, the distance between:')
    norm=np.empty([numclus,numens])
    finalOUTPUT=[]
    repres=[]
    for nclus in range(numclus):
        for ens in range(numens):
            normens=centroids[nclus,:]-PCs[ens,:]
            norm[nclus,ens]=math.sqrt(sum(normens**2))
            #print('The distance between centroid of cluster {0} and member {1} is {2}'.format(nclus,ens,round(norm[nclus,ens],3)))
        print('The distances between centroid of cluster {0} and member #0 to #{1} are:\n{2}'.format(nclus,numens-1,np.round(norm[nclus],3)))
        print('MINIMUM DISTANCE FOR CLUSTER {0} IS {1} --> member #{2}'.format(nclus,round(norm[nclus].min(),3),list(np.where(norm[nclus] == norm[nclus].min())[0])))
        repres.append(np.where(norm[nclus] == norm[nclus].min())[0][0])
        print('MAXIMUM DISTANCE FOR CLUSTER {0} IS {1} --> member #{2}\n'.format(nclus,round(norm[nclus].max(),3),list(np.where(norm[nclus] == norm[nclus].max())[0])))

        txt='Closest ensemble member/members to centroid of cluster {0} is/are {1}\n'.format(nclus,list(np.where(norm[nclus] == norm[nclus].min())[0]))
        finalOUTPUT.append(txt)
    with open(os.path.join(dir_OUTPUT,'RepresentativeEnsembleMembers_{0}.txt'.format(name_outputs)), "w") as text_file:
        text_file.write(''.join(str(e) for e in finalOUTPUT))

    #____________Save the most representative ensemble members
    namef=os.path.join(dir_OUTPUT,'repr_ens_{0}.txt'.format(name_outputs))
    np.savetxt(namef,repres,fmt='%i')
        
    print('____________________________________________________________________________________________________________________')
    print('In order to study the spread of each cluster,')
    print('the standard deviation of the distances between each member in a cluster and the cluster centroid is computed in the PC space')
    print('____________________________________________________________________________________________________________________')
    print('\nIn the PC space:')
    statOUTPUT=[]
    for nclus in range(numclus):
        members=L[nclus][2]
        norm=np.empty([numclus,len(members)])
        for mem in range(len(members)):
            #print('mem=',mem)
            ens=members[mem]
            #print('ens',ens)
            normens=centroids[nclus,:]-PCs[ens,:]
            norm[nclus,mem]=math.sqrt(sum(normens**2))
            #print('norm=',norm[nclus],norm.dtype)
        print('the distances between centroid of cluster {0} and its belonging members {1} are:\n{2}'.format(nclus,members,np.round(norm[nclus],3)))
        print('MINIMUM DISTANCE WITHIN CLUSTER {0} IS {1} --> member #{2}'.format(nclus,round(norm[nclus].min(),3),members[np.where(norm[nclus] == norm[nclus].min())[0][0]]))
        print('MAXIMUM DISTANCE WITHIN CLUSTER {0} IS {1} --> member #{2}'.format(nclus,round(norm[nclus].max(),3),members[np.where(norm[nclus] == norm[nclus].max())[0][0]]))
        print('INTRA-CLUSTER STANDARD DEVIATION FOR CLUSTER {0} IS {1}\n'.format(nclus,norm[nclus].std()))
        
        d_stat=collections.OrderedDict()
        d_stat['cluster']=nclus
        d_stat['member']=members
        d_stat['d_to_centroid']=np.round(norm[nclus],3)
        d_stat['intra-clus_std']=norm[nclus].std()
        d_stat['d_min']=round(norm[nclus].min(),3)
        d_stat['d_max']=round(norm[nclus].max(),3)
        d_stat['freq(%)']=round(L[nclus][1],3)
        stat=pd.DataFrame(d_stat)
        statOUTPUT.append(stat)
    statOUTPUT = pd.concat(statOUTPUT, axis=0)
    #____________Save statistics of cluster analysis
    namef=os.path.join(dir_OUTPUT,'statistics_clutering_{0}.txt'.format(name_outputs))
    with open(namef, 'w') as text_file:
        text_file.write(statOUTPUT.__repr__())

    return


#========================================================

if __name__ == '__main__':
    print('This program is being run by itself')
    
    print('**************************************************************')
    print('Running {0}'.format(sys.argv[0]))
    print('**************************************************************')
    dir_OUTPUT    = sys.argv[1]  # OUTPUT DIRECTORY
    name_outputs  = sys.argv[2]  # name of the outputs
    numens        = int(sys.argv[3])  # number of ensemble members
    numpcs        = sys.argv[4]  # number of retained PCs
    perc          = sys.argv[5]  # percentage of explained variance by PCs
    numclus       = int(sys.argv[6])  # number of clusters
    
    ens_eof_kmeans(dir_OUTPUT,name_outputs,numens,numpcs,perc,numclus)
    
else:
    print('ens_eof_kmeans is being imported from another module')

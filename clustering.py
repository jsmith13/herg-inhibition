### This script contains functions for clustering compounds using Morgan connectivity fingerprints and the Butina clustering algorithm.

## Import required libraries
# chemistry functions
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

# other libraries
import numpy as np
from matplotlib import pyplot as plt

## declare a function cluster_compounds
# SMILES: a list of SMILES strings
# threshhold: the distance threshold to use for assigning clusters
# neighbors: the number of neighboring atoms to fingerprint
# raw_clusters: a boolean flag, returns an unsorted list of clusters when True
# returns a list containing cluster assignments for each element of SMILES
def cluster_compounds(SMILES, threshold, neighbors = 2, raw_clusters = False):
    # convert SMILES strings first into RDKit format, then into Morgan fingerprints
    fingerprints = list(map(lambda x: GetMorganFingerprintAsBitVect(x, neighbors), map(Chem.MolFromSmiles, SMILES)))
    
    # calculate Tanimoto distance matrix between molecular fingerprints
    # declare an empty list to hold the distances
    distances = []
    
    # calculate Tanimoto distances between molecular fingerprints
    # the clustering algorithm requires a list containing only unique distances
    for i in range(len(fingerprints)):
        for j in range(i+1, len(fingerprints)):
            distances.append(1 - DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j]))
            
    # cluster the compounds using the Butina method
    clusters = Butina.ClusterData(distances, nPts = len(SMILES), distThresh = threshold, isDistData = True)
    
    # declare a list to hold cluster assignments
    assignments = [np.nan for i in range(len(SMILES))]
    
    # map cluster assignments to elements in SMILES
    for cluster in range(len(clusters)):
        for compound in clusters[cluster]:
            assignments[compound] = cluster
    
    # return cluster assignments sorted or not
    if raw_clusters == True:
        return(clusters)
    else:
        return(assignments)


## declare a function tune_clustering_threshold
# SMILES: a list of SMILES strings
# thresholds: a list of distances to test as thresholds
# neighbors: the number of neighboring atoms to fingerprint
# plot: a boolean flag indicating whether to return a rough plot of the results
# returns the number of clusters and median cluster size
def tune_clustering_threshold(SMILES, thresholds, neighbors = 2, plot = False):
    # generate clusters for each distance in thresholds
    clusters = list(map(lambda x: cluster_compounds(SMILES, x, neighbors, raw_clusters = True), thresholds))
    
    # calculate number of clusters and median cluster size, excluding singletons
    number_of_clusters = list(map(lambda x: len([len(cluster) for cluster in x if len(cluster) > 1]), clusters))
    median_cluster_size = list(map(lambda x: np.median([len(cluster) for cluster in x if len(cluster) > 1]), clusters))
    
    # plot the results if requested
    if plot == True:
        # create an 1 x 2 grid of axes
        fig, axs = plt.subplots(1, 2, tight_layout = True)
        
        # add scatterplots
        axs[0].plot(thresholds, number_of_clusters)
        axs[1].plot(thresholds, median_cluster_size)
        
        # label the axis
        axs[0].set_xlabel("distance")
        axs[0].set_ylabel("number of clusters")
        axs[1].set_xlabel("distance")
        axs[1].set_ylabel("median cluster size")
        
        # display the figure
        fig.show()
    
    # return the number of clusters and median cluster size as a tuple of lists
    return(number_of_clusters, median_cluster_size)
    

### This script contains a function for clustering compounds using Morgan connectivity fingerprints and the Butina clustering algorithm.

## Import required libraries
# chemistry functions
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

# other libraries
import numpy as np


## declare a function cluster_compounds
## SMILES: a list of SMILES strings
## neighbors: the number of neighboring atoms to fingerprint
def cluster_compounds(SMILES, threshold, neighbors = 2):
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
            assignments = cluster
    
    # return cluster assignments
    return(clusters)


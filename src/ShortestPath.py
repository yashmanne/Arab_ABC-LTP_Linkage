# Import necessary Libraries
import os
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import collections
import pickle # used to store cleaned data
import numpy.ma as ma # masking for numpy Arrays
boolVal = False
attedDataName = 'Overlapping_GenesMR.p'
bioGRIDfileName = 'BioGRID_with_ATTEDMR.p'
attedData = pickle.load(open(attedDataName, 'rb'))
genesGRID = pickle.load(open(bioGRIDfileName, 'rb'))
# Initializes Networkx Graph
genesG = nx.Graph()
# Adds edges for each interaction
for i in range(genesGRID.shape[0]):
    genesG.add_edge(genesGRID['Entrez Gene Interactor A'].loc[i], genesGRID['Entrez Gene Interactor B'].loc[i], weight = genesGRID['MR'].loc[i])

def nonConvertible(orig,trans): 
    '''
    Checks if all orginal values are in new list of values 
    and returns the ones that aren't
    Inputs: List of original values, list of new values
    Output: List of all original values that are not in the new value list
    ''' 
    nc = []
    for i in range(len(orig)):
        if orig[i] not in trans:
            nc.append(orig[i])
        else:
            continue
    return nc

def existing(lst, reference): 
    '''
    Checks if items of lst are in reference, returns a list of existing values
    Inputs: List to check, list to cross-reference with
    Output: List of all values that also exist in the reference list.
    ''' 
    exists = [] # stores the existing values
    NE = [] # stores the non-existing values
    for x in lst:
        if x in reference:
            exists.append(x)
        else:
            NE.append(x)
    return exists, NE

def ShortestDistances(Graph,grpA,grpB,dataFrame = True): #using function caused Kernel to crash due to memory
    '''
    Finds the shortest distances in a graph between all members of group 1 
    and all members of group 2 and stores into a dict.
    Inputs: Graph to scan, list of group A nodes, list of group B nodes
    Output: Returns a pandas dataFrame with group A as columns and group B as rows 
            and the shortest-distance as the values. 
            Else, returns a dict with group A as keys and a list of shortest distances 
            as values (size of list is the length of group B)
    '''
    graphGenes = list(Graph.nodes)
    abc, Uabc = existing(grpA,graphGenes)
    ltp, Ultp = existing(grpB,graphGenes)
    DF={0:ltp}
    for x in abc:
        valList = []
        for y in ltp:
            if (x in graphGenes) and (y in graphGenes): # checks if both genes are part of the overlapping set
                val = nx.astar_path_length(Graph, x,y)
            else:
                val = np.nan
            valList.append(val)
        DF.update({x:valList})
    if dataFrame:
        Data = pd.DataFrame(DF,dtype='float64')
        Data = Data.astype({0:'int'})
        Data = Data.set_index(0)
        return Data
    else:
        return DF
# option: filter out nodes with 1 connection
# store as a dictionary, dict: 
def simplifyGraph(grpA, grpB, Graph, kdeg):
    '''
    Gathers all nodes within kdeg degrees of separation from the given groups.
    Inputs: list of group A nodes, list of group B nodes, graph to simplify, int representing
            the degree of seperation from the key nodes allow.
    Output: subgraph of inputted graph as well as dictionary with each node in the new graph as key
            and # of neighboring nodes in group A, # of neighboring nodes in group B, and 
            list of all neighboring nodes.
    '''
    keyGenes = sorted(set(grpA+grpB))
    i = kdeg
    while i >0:
        # Stores neighbors at each key gene
        neighbors = []
        # Loops through all key genes 
        for gene in keyGenes:
            nodes = [n for n in Graph.neighbors(gene)]
            neighbors.append(nodes)
        neighbors.append(keyGenes)
        neighbors = [item for sublist in neighbors for item in sublist]
        keyGenes = sorted(set(neighbors))
        i-=1
    # Creates subgraph
    newGraph = Graph.subgraph(keyGenes)
    # Iterates through the new nodes and checks if neighbors are in either group A or B
    geneDict = {}
    for gene in newGraph.nodes:
        nodes = [n for n in newGraph.neighbors(gene)]
        grpACount = 0
        grpBCount = 0
        for ngene in nodes:
            if ngene in grpA:
                grpACount +=1
            elif ngene in grpB:
                grpBCount +=1
        geneDict[gene] = (grpACount,grpBCount, nodes)
    return newGraph, geneDict


# Read In Data 
ABC_trans = pd.read_csv('convABC_Genes.txt', sep = '\t')
ABC_trans = ABC_trans.rename(columns = {'From':'TAIR_ID','To':'ENTREZ_ID'})
ABC_orig = pd.read_excel('ABC_Genes.xls', sheet_name = 'Sheet2')
abcAT = sorted(set(list(ABC_trans['TAIR_ID'])))
abcEntrez = sorted(set(list(ABC_trans['ENTREZ_ID'])))
abcAT_orig = sorted(set(list(ABC_orig['TAIR_ID'])))
abc_NC  = nonConvertible(abcAT_orig,abcAT)
LTP_trans = pd.read_csv('convLTP_Genes.txt', sep = '\t')
LTP_trans = LTP_trans.rename(columns = {'From':'TAIR_ID','To':'ENTREZ_ID'})
LTP_orig = pd.read_excel('LTP_Genes.xlsx', sheet_name = 'Sheet1')
ltpAT = sorted(set(list(LTP_trans['TAIR_ID'])))
ltpEntrez = sorted(set(list(LTP_trans['ENTREZ_ID'])))
ltpAT_orig = sorted(set(list(LTP_orig['TAIR_ID'])))
ltp_NC  = nonConvertible(ltpAT_orig,ltpAT)
abc, Uabc = existing(abcEntrez,genesG.nodes)
ltp, Ultp = existing(ltpEntrez,genesG.nodes)
# Run Shortest Path & store as pickle files
simpleShortDistFileName = 'simplifiedShortestDistances.p'
shortDistFileName = 'shortestDistances.p'
simpleG, simpleGDict = simplifyGraph(abc,ltp, genesG,1)
simple_abc_ltp_shortDist = ShortestDistances(genesG,abcEntrez,ltpEntrez, True)
pickle.dump(abc_ltp_shortDist, open(shortDistFileName, 'wb'))
simple_abc_ltp_shortDist = ShortestDistances(simpleG,abcEntrez,ltpEntrez, True)
pickle.dump(simple_abc_ltp_shortDist, open(simpleShortDistFileName, 'wb'))
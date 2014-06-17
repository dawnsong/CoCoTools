import copy
import os, sys
import numpy as np
import networkx as nx
import pickle


#Load networkx graph
#filedir=('C:\\Users\\Robert\\cocomac\\cocomac-tools-rsblume\\applications\\ModhaSingh_PFC\\)
filedir=('/home/despo/jrcohen/gitrepos/graphtheory/newman_partition_test')
filestub=filedir+'G3.pck'
G3=pickle.load(open(filestub, "rb"))


cmat_temp=nx.to_numpy_matrix(G3) # convert nx graph to np matrix adjacency matrix

cmat_upper=np.zeros([cmat_temp.shape[0],cmat_temp.shape[1]]) #initialize NxN zeros matrix

for x in range(cmat_temp.shape[0]):
    cmat_upper[x,x:]=cmat_temp[x,x:] #add only upper quadrant
    cmat_upper[x,x]=np.nan #nan the diagonal

np.save('robB_g3.npy', cmat_upper)    
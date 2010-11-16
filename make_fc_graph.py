"""Load FC flat text files into networkx
"""

import networkx as nx
import numpy as np

# Main script

labels_fname = 'fc_network_names.txt'
edges_fname = 'fc_network_edges.txt'

fc_g = nx.DiGraph() 

with open(labels_fname) as labels: 
    for line in labels:
        num, acronym = line.split()
        num = int(num)
        fc_g.add_node(num, acronym=acronym)

edges = np.loadtxt(edges_fname)

fc_g.add_edges_from(edges)

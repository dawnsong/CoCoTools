import networkx as nx
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import csv
import json
import brainx.modularity as mod



# Functions
#First pass
def make_remove_node_graph_dict(in_graph, labels_list):
	"""builds a dict with the removed node as the keys and the resulting graph as values, that is the key is not in the graph
	Parameters
	-----------
	in_graph : networkx graph
	networkx graph with strings as node labels
	labels_list: list 
	list with the names of for txt files with subgraph nodes
	Returns
	-----------
	graph_dict :   dict with remove nodes as keys and graphs as values"""
	graph_dict=dict()
	for x in range(0,len(labels_list)):
		remove_node=labels_list[x]
		in_graph_remove=in_graph.copy()
		in_graph_remove.remove_node(remove_node)
		graph_dict[remove_node]=in_graph_remove
	return graph_dict

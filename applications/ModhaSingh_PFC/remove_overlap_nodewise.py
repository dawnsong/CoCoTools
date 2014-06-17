#This notebook is for additional analysis in which one node at a time is removed and modularity is computed

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
	for x in range(0,size(labels_list)):
		remove_node=labels_list[x]
		in_graph_remove=in_graph.copy()
		in_graph_remove.remove_node(remove_node)
		graph_dict[remove_node]=in_graph_remove
	return graph_dict

def make_mod_dict_from_separate_dicts(graph_dict,q_dict, part_dict):
	""" used to convert all the work I have already done into the better format graph_dicts"""
	new_mod_dict=dict()
	for keys in graph_dict:
		values_of_new_mod_dict=dict()
		values_of_new_mod_dict['q']=q_dict[keys]
		values_of_new_mod_dict['part']=part_dict[keys]
		new_mod_dict[keys]=values_of_new_mod
	return new_mod_dict

def run_sim_store_as_dict(g):
	"""runs one pass of simulated annealing, stores results in a nice dict"""
	mod_out=dict()
	mod_raw_out=mod.simulated_annealing(g)
	mod_out['q']=mod_raw_out[0].modularity()
	mod_out['part']=mod_raw_out[0].index_as_node_names()
	return mod_out

	
def run_sim_loop_with_graph_dict(graph_dict, filename):
	""""loops through simulated annealing stores as mod_dict"""
	for keys in graph_dict:
		mod_out_dict[keys]=run_sim_store_as_dict(graph_dict[keys])
		out_name=file_name + "save.pck"
		out_file = open(out_name, 'wb')
		pickle.dump(first_pass_mod_out, out_file)
	return mod_out_dict

#Second Pass
def make_mod_subgraph(nodelist, in_graph):
	"""makes subgraph, useful shorthand"""
	subg=in_graph.copy()
	subg.subgraph(nodelist)
	return subg

	
def make_subgraph_dict(mod_out_dict, in_graph):
    """"makes subgraph_dict"""
    subgraph_dict=dict()
    for keys in mod_out_dict:
        if mod_out_dict[keys]['q']>=.30 and size(mod_out_dict[keys]['part'])<5:
            part_subgraph=dict()
            for x in range(0,size(mod_out_dict[keys]['part'])):
                part_nodelist=mod_out_dict[keys]['part'][x]
                part_subgraph[x]=make_mod_subgraph(part_nodelist, in_graph)
            subgraph_dict[keys]=part_subgraph
    return subgraph_dict
	
def run_sim_loop_with_subgraph_dict(subgraph_dict, filename):
    mod_out_dict_subgraph=dict()
    for keys in subgraph_dict:
        mod_out_part_dict=dict()
        for xkeys in subgraph_dict[keys]:
            mod_out_part_dict[xkeys]=run_sim_store_as_dict(subgraph_dict[keys][xkeys])
            #mod_out_dict_subgraph[keys][xkeys]=run_sim_store_as_dict(subgraph_dict[keys][xkeys])
        mod_out_dict_subgraph[keys]=mod_out_part_dict
        out_name=file_name + "save.pck"
        out_file=open(out_name, 'wb')
        pickle.dump(first_pass_mod_out, out_file)
    return mod_out_dict_subgraph
			
				
	
#Inputs
filestub=('C:\\Users\\Robert\\gitrepos\\cocomac-tools-rsblume\\applications\\ModhaSingh_PFC\\modha\\lowest_full_g.pck')
whole_low=pickle.load(open(filestub,'rb'))

infile = open("pfc_remove_first_part.pck",'rb')
pfc_remove_first_part_dict = pickle.load(infile)
infile = open("pfc_remove_first_mod.pck",'rb')
pfc_remove_first_mod_dict = pickle.load(infile)
infile = open("pfc_remove_first_part.pck",'rb')
pfc_remove_first_graph_dict=pickle.load(infile)

#Basic processing on graph for modularity analyses
whole_low=whole_low.to_undirected()
labels_list=['pfc_labels', 'mod_dlpfc_labels', 'mod_vlpfc_labels', 'mod_mpfc_labels', 'lit_dlpfc_labels', 'lit_vlpfc_labels', 'mod_lateral_labels', 'mod_rostral_labels', 'mod_premotor_labels']
labels={x : list(csv.reader(open(x + '.csv')))[0] for x in labels_list}
pfc_graph=whole_low.subgraph(labels['pfc_labels'])

# Begin Script
graph_dict=make_remove_node_graph_dict(pfc_graph, labels['pfc_labels'])
first_pass=make_mod_dict_from_separate_dicts(graph_dict,pfc_remove_first_mod_dict, pfc_remove_first_part_dict)
second_pass_graphs=make_subgraph_dict(first_pass, whole_low)

second_pass_graph_batch=dict()
for x in second_pass_graphs.keys()[:10]:
    second_pass_graph_batch[x]=second_pass_graphs[x]
    

second_pass_mod_out=run_sim_loop_with_subgraph_dict(second_pass_graph_batch, 'subgraph')

#second_pass_mod_out=run_sim_loop_with_subgraph_dict(subgraph_dict, 'subgraph')



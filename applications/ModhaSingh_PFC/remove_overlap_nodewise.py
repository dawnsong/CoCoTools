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
import numpy as np


################ Functions #################
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

def make_mod_dict_from_separate_dicts(graph_dict,q_dict, part_dict):
	""" used to convert all the work I have already done into the better format graph_dicts"""
	new_mod_dict=dict()
	for keys in graph_dict:
		values_of_new_mod_dict=dict()
		values_of_new_mod_dict['q']=q_dict[keys]
		values_of_new_mod_dict['part']=part_dict[keys]
		new_mod_dict[keys]=values_of_new_mod_dict
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
	        mod_out_dict=dict()
		mod_out_dict[keys]=run_sim_store_as_dict(graph_dict[keys])
		out_name=filename + "save.pck"
		out_file = open(out_name, 'wb')
		pickle.dump(mod_out_dict, out_file)
	return mod_out_dict

#Second Pass
def make_mod_subgraph(nodelist, in_graph):
	"""makes subgraph, useful shorthand"""
	subg=in_graph.copy()
	subg=subg.subgraph(nodelist)
	return subg

	
def make_subgraph_dict(mod_out_dict, in_graph):
    """"makes subgraph_dict"""
    subgraph_dict=dict()
    for keys in mod_out_dict:
        if mod_out_dict[keys]['q']>=.30 and len(mod_out_dict[keys]['part'])<5:
            part_subgraph=dict()
            for x in range(0,len(mod_out_dict[keys]['part'])):
                part_nodelist=mod_out_dict[keys]['part'][x]
                part_subgraph[x]=make_mod_subgraph(part_nodelist, in_graph)
            subgraph_dict[keys]=part_subgraph
    return subgraph_dict
	
def run_sim_loop_with_subgraph_dict(subgraph_dict, filename):
    mod_out_dict_subgraph=dict()
    for keys in subgraph_dict:
        mod_out_part_dict=dict()
        for xkeys in subgraph_dict[keys]:
            if len(subgraph_dict[keys][xkeys].nodes())>6:
                mod_out_part_dict[xkeys]=run_sim_store_as_dict(subgraph_dict[keys][xkeys])
            else:
                mod_out_part_dict[xkeys]="n/a"          
            #mod_out_dict_subgraph[keys][xkeys]=run_sim_store_as_dict(subgraph_dict[keys][xkeys])
        mod_out_dict_subgraph[keys]=mod_out_part_dict
        out_name=filename + "save.pck"
        out_file=open(out_name, 'wb')
        pickle.dump(mod_out_dict_subgraph, out_file)
    return mod_out_dict_subgraph
			
				
	
#############Inputs#########################
filestub=("/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/modha/lowest_full_g.pck")
whole_low=pickle.load(open(filestub,'rb'))

infile = open("/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/pfc_remove_first_part.pck",'rb')
pfc_remove_first_part_dict = pickle.load(infile)
infile = open("/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/pfc_remove_first_mod.pck",'rb')
pfc_remove_first_mod_dict = pickle.load(infile)
infile = open("/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/pfc_remove_first_part.pck",'rb')
pfc_remove_first_graph_dict=pickle.load(infile)

#Basic processing on graph for modularity analyses
whole_low=whole_low.to_undirected()
labels_list=['pfc_labels', 'mod_dlpfc_labels', 'mod_vlpfc_labels', 'mod_mpfc_labels', 'lit_dlpfc_labels', 'lit_vlpfc_labels', 'mod_lateral_labels', 'mod_rostral_labels', 'mod_premotor_labels']
labels={x : list(csv.reader(open("/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/" + x + '.csv')))[0] for x in labels_list}
pfc_graph=whole_low.subgraph(labels['pfc_labels'])

nodelist=labels['mod_dlpfc_labels']+labels['mod_vlpfc_labels']+labels['mod_premotor_labels'] + labels['mod_rostral_labels']+ labels['mod_mpfc_labels']



######### Begin Script##############
graph_dict=make_remove_node_graph_dict(pfc_graph, labels['pfc_labels'])
first_pass=make_mod_dict_from_separate_dicts(graph_dict,pfc_remove_first_mod_dict, pfc_remove_first_part_dict)

second_pass_graphs=make_subgraph_dict(first_pass, pfc_graph)
second_pass_graph_batches=dict()

#I have already run all the modularity passes and stored the results in subgraphsave_all_50.pck
#Lets load it:
second_pass_batches=dict()
filename="/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/subgraphsave_all_50.pck"
second_pass_graph_all=pickle.load(open(filename,'rb'))
infile=open(filename, 'rb')
second_pass_graph_all=pickle.load(infile)

filenames = ["/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/subgraphsave.pck", "/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/subgraph2save.pck","/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/subgraph3save.pck"]
for x in filenames:
    infile=open(x, "rb")
    second_pass_batch_cat=pickle.load(infile)
    for keys in second_pass_batch_cat:
       second_pass_batches[keys]= second_pass_batch_cat[keys]
       

#Now determine which second pass Q's >= 0.30
# I initially tried Q of 0.30 but  many right on the cusp...
second_pass_accepted=dict()
for keys in second_pass_graph_all:
    second_pass_accepted_list=[]
    for xkeys in second_pass_graph_all[keys].keys():
        if second_pass_graph_all[keys][xkeys]!='n/a' and second_pass_graph_all[keys][xkeys]['q']>= 0.28:
            second_pass_accepted_list.append(xkeys)
    second_pass_accepted[keys]=second_pass_accepted_list
           

# Now build combined partition dicts with first pass and for those first pass parts that were accepted insert second pass partitions
combined_parts_dict=first_pass.copy()
#first find the graphs that were not modular at first pass
for keys in first_pass:
    if first_pass[keys]['q']<0.3000:
        combined_parts_dict[keys]=['n/a']
        second_pass_accepted[keys]=['n/a']

#now replace modules with submodules for those that made it
for keys in combined_parts_dict:
    combo_first_second_pass_part_list=[]
    first_pass_nparts_remove=range(0,len(first_pass[keys]['part']))
    for values in second_pass_accepted[keys]:
        if values != 'n/a':
            for x,parts in enumerate(second_pass_batches[keys][values]['part']):
                combo_first_second_pass_part_list.append(parts)
            #combo_first_second_pass_part_list.append(second_pass_batches[keys][values]['part']) # needs to be appended
            first_pass_nparts_remove.remove(values)
    for x in first_pass_nparts_remove:
        combo_first_second_pass_part_list.append(first_pass[keys]['part'][x])
    combined_parts_dict[keys]=combo_first_second_pass_part_list
    
    

#Now build adja/correlation matrices with the combined_parts_dict

def _build_relabel_mapping_dict(nodes_list):
    """"builds the mapping dict that is needed to relabel the complete graphs"""
    relabel_mapping_list=zip(range(0,len(nodes_list)),nodes_list)
    relabel_mapping_dict=dict(relabel_mapping_list)
    return relabel_mapping_dict

def _build_complete_graph_from_part(part_list):
    """builds a complete graph from a partition"""
    complete_graph_from_part=nx.complete_graph(len(part_list))
    relabel_map= _build_relabel_mapping_dict(part_list)
    complete_graph_from_part=nx.relabel_nodes(complete_graph_from_part,relabel_map)
    return complete_graph_from_part
    
def _build_one_adj_mat_from_parts(parts_lists,node_list):
    """"build adj matrix from all the partitions for one key in a graph dict"""
    adj_mat_graph=nx.Graph()
    adj_mat_graph.add_nodes_from(nodelist)
    for x in parts_lists:
        part_graph=_build_complete_graph_from_part(x)
        adj_mat_graph.add_edges_from(part_graph.edges())
    adj_mat=nx.to_numpy_matrix(adj_mat_graph,node_list)
    return (adj_mat, adj_mat_graph)


def build_combined_adj_mat(combo_parts_dict,nodelist):
    """ build a adjacency matrix that adds up all the single adj matrices"""
    combined_adj_mat=np.matrix(np.zeros((len(nodelist),len(nodelist))))
    for x in nodelist:
        (single_adj_mat,single_graph)=_build_one_adj_mat_from_parts(combo_parts_dict[x],nodelist)
        combined_adj_mat=combined_adj_mat+single_adj_mat
    return combined_adj_mat
 
combined_adj_mat=build_combined_adj_mat(combined_parts_dict,nodelist)   
    
imshow(combined_adj_mat,interpolation='nearest',cmap=plt.cm.binary)

    


    

    

foo_graph=nx.Graph()
for x in combined_parts_dict[nodelist[0]]:
    test=_build_complete_graph_from_part(x)
    foo_graph.add_edges_from(test.edges())



snork=zip(range(0,len(combined_parts_dict[nodelist[0]][0])),combined_parts_dict[nodelist[0]][0])


foo= np.zeros((len(nodelist),len(nodelist)))

    







second_pass_graph_all=second_pass_batches.copy()

#out_name="subgraphsave_all.pck"
#out_file=open(out_name, 'wb')
#pickle.dump(second_pass_graph_batches,out_file )

#filename="/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/subgraphsave_all.pck"
#second_pass_graph_all=pickle.load(open(filename,'rb'))
#infile=open(filename, 'rb')
#second_pass_graph_all=pickle.load(infile)

#second_pass_graph_batch=dict()
#for x in second_pass_graphs.keys()[50:]:
#    second_pass_graph_batch[x]=second_pass_graphs[x]
    

#second_pass_mod_out=run_sim_loop_with_subgraph_dict(second_pass_graph_batch, 'subgraph4')

#for keys in second_pass_mod_out:
#    second_pass_batches[keys]= second_pass_mod_out[keys]
    
filename="/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/subgraphsave_all_50.pck"
outfile=open(filename, 'wb')
pickle.dump(second_pass_batches, outfile)


#second_pass_mod_out=run_sim_loop_with_subgraph_dict(subgraph_dict, 'subgraph')



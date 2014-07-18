#This script tests whether remove groups of nodes that overlap with other nodes 
#will affect the dorso-ventral or rostro-caudal nature of partitions

###Imports#####
import pickle
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl 
import csv
import brainx.modularity as mod
import numpy as np




########### Functions ################################
#First Pass
def make_remove_overlap_graphs(no_overlap_group, base_graph):
    no_overlap_graph_dict=dict()
    no_overlap_0_g=base_graph.copy()
    no_overlap_1_g=base_graph.copy()
    no_overlap_0_g.remove_nodes_from(no_overlap_group[0])
    no_overlap_1_g.remove_nodes_from(no_overlap_group[1])
    no_overlap_graph_dict=dict(zip([0,1], [no_overlap_0_g,  no_overlap_1_g]))
    return no_overlap_graph_dict

#def run_sim_with_no_overlap_g(no_overlap_g):
#    no_overlap_dict=dict()
#    part=mod.simulated_annealing(no_overlap_g)
#    no_overlap_dict['part']=part[0].index_as_node_names()
#    no_overlap_dict['q']=part[0].modularity()
#    return  no_overlap_dict
    
def run_sim_store_as_dict(g):
	"""runs one pass of simulated annealing, stores results in a nice dict"""
	mod_out=dict()
	mod_raw_out=mod.simulated_annealing(g)
	mod_out['q']=mod_raw_out[0].modularity()
	mod_out['part']=mod_raw_out[0].index_as_node_names()
	return mod_out    

def run_results(no_overlap_graph_dict):
    results_dict=dict()
    for keys in no_overlap_graph_dict:
        results_dict[keys]=run_sim_store_as_dict(no_overlap_graph_dict[keys])
    return results_dict

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
        if mod_out_dict[keys]['q']>=.297 and len(mod_out_dict[keys]['part'])<5:
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


#Private functions to combine data into adj mat
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
    adj_mat_graph.add_nodes_from(node_list)
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


#############Inputs#########################
filestub=("/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/modha/lowest_full_g.pck")
whole_low=pickle.load(open(filestub,'rb'))


#Basic processing on graph for modularity analyses
whole_low=whole_low.to_undirected()
labels_list=['pfc_labels', 'mod_dlpfc_labels', 'mod_vlpfc_labels', 'mod_mpfc_labels', 'lit_dlpfc_labels', 'lit_vlpfc_labels', 'mod_lateral_labels', 'mod_rostral_labels', 'mod_premotor_labels']
labels={x : list(csv.reader(open("/home/rblumenf/gitrepos/CoCoTools/applications/ModhaSingh_PFC/" + x + '.csv')))[0] for x in labels_list}
pfc_graph=whole_low.subgraph(labels['pfc_labels'])

nodelist=labels['mod_dlpfc_labels']+labels['mod_vlpfc_labels']+labels['mod_premotor_labels'] + labels['mod_rostral_labels']+ labels['mod_mpfc_labels']

###############################################
PMC_overlap=dict(zip((0,1), (['6Va', '6Vb', '4c', 'PrCO'],['F4', 'F5', 'PrOM'])))
SMA_overlap=dict(zip((0,1), (['SMAr', 'SMAc'], ['M2-FL', 'M2-HL'])))
FP_overlap=dict(zip((0,1), (['10o', '10m'], ['10v', '10d'])))
###############################################

#Run First Pass
pfc_graph_PMC=make_remove_overlap_graphs(PMC_overlap,pfc_graph)
pfc_graph_SMA=make_remove_overlap_graphs(SMA_overlap,pfc_graph)
pfc_graph_FP=make_remove_overlap_graphs(FP_overlap,pfc_graph)


PMC_results=run_results(pfc_graph_PMC)
SMA_results=run_results(pfc_graph_SMA)
FP_results=run_results(pfc_graph_FP)

#Run Second Pass
PMC_second_pass_g = make_subgraph_dict(PMC_results, pfc_graph)
PMC_second_pass_results= run_sim_loop_with_subgraph_dict(PMC_second_pass_g, 'PMC_remove_overlap')

SMA_second_pass_g = make_subgraph_dict(SMA_results, pfc_graph)
SMA_second_pass_results= run_sim_loop_with_subgraph_dict(SMA_second_pass_g, 'SMA_remove_overlap')

FP_second_pass_g = make_subgraph_dict(FP_results, pfc_graph)
FP_second_pass_results= run_sim_loop_with_subgraph_dict(FP_second_pass_g, 'SMA_remove_overlap')


#Now determine which second pass Q's >= 0.30
# I initially tried Q of 0.30 but  many right on the cusp...
second_pass_accepted=dict()
for keys in PMC_second_pass_results:
    second_pass_accepted_list=[]
    for xkeys in PMC_second_pass_results[keys].keys():
        if PMC_second_pass_results[keys][xkeys]!='n/a' and PMC_second_pass_results[keys][xkeys]['q']>= 0.28:
            second_pass_accepted_list.append(xkeys)
    second_pass_accepted[keys]=second_pass_accepted_list

# Now build combined partition dicts with first pass and for those first pass parts that were accepted insert second pass partitions
combined_parts_dict=PMC_results.copy()
#first find the graphs that were not modular at first pass
for keys in PMC_results:
    if PMC_results[keys]['q']<0.297:
        PMC_results[keys]=['n/a']
        second_pass_accepted[keys]=['n/a']


#now replace modules with submodules for those that made it
for keys in combined_parts_dict:
    combo_first_second_pass_part_list=[]
    first_pass_nparts_remove=range(0,len(PMC_results[keys]['part']))
    for values in second_pass_accepted[keys]:
        if values != 'n/a':
            for x,parts in enumerate(PMC_second_pass_results[keys][values]['part']):
                combo_first_second_pass_part_list.append(parts)
            #combo_first_second_pass_part_list.append(second_pass_batches[keys][values]['part']) # needs to be appended
            first_pass_nparts_remove.remove(values)
    for x in first_pass_nparts_remove:
        combo_first_second_pass_part_list.append(PMC_results[keys]['part'][x])
    combined_parts_dict[keys]=combo_first_second_pass_part_list
    


#PMC 0:
a=[y for x in combined_parts_dict[0] for y in x]
a_mat, a_g=_build_one_adj_mat_from_parts(combined_parts_dict[0],a)

ax=plt.figure().add_subplot(111)
plt.imshow(a_mat,interpolation='nearest'),plt.colorbar()
plt.xticks(range(0,len(a)))
plt.yticks(range(0,len(a)))
ax.set_xticklabels(a)
ax.set_yticklabels(a)

#PMC 1:
b=[y for x in combined_parts_dict[1] for y in x]
b_mat, b_g=_build_one_adj_mat_from_parts(combined_parts_dict[1],b)





a=nx.complete_graph(len(combined_parts_dict[0]))
map_a=_build_relabel_mapping_dict(combined_parts_dict[0][0])



#Needs work here ---> the keys are 0,1 this time
combined_adj_mat=build_combined_adj_mat(combined_parts_dict,nodelist)
_build_complete_graph_from_part(combined_parts_dict[0][0])



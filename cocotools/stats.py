from __future__ import division
import copy

import numpy as np
import scipy.io
import networkx as nx


def random_stats(g, n_rand):
    """Return clustering coeffs and char. path lengths for random graphs.

    Parameters
    ----------
    g : NetworkX DiGraph

    n_rand : number of random graphs to generate

    Returns
    -------
    clust_coeffs
    char_paths

    Notes
    -----
    This routine matches randmio_dir in the Sporns Matlab toolbox, with ITER
    set to 1.
    """
    clust_coeffs = []
    char_paths = []
    A = nx.adjacency_matrix(g)
    n = A.shape[0]
    K = len(x_indices)
    max_attempts = round(n*K/(n*(n-1)))
    for i in range(n_rand):
        R = copy.deepcopy(A)
        x_indices, y_indices = np.nonzero(A)
        att = 0
        while att <= max_attempts:
            while True:
                r1, r2 = np.random.rand(2)
                e1 = np.floor(K*r1)
                e2 = np.floor(K*r2)
                while e1 == e2:
                    e2 = np.floor(K*np.random.rand(1)[0])
                a = x_indices[e1]
                b = y_indices[e1]
                c = x_indices[e2]
                d = y_indices[e2]
                if a not in (c, d) and b not in (c, d):
                    break
            if not (R[a,d] or R[c,b]):
                R[a,d] = R[a,b]
                R[a,b] = 0
                R[c,b] = R[c,d]
                R[c,d] = 0
                y_indices[e1] = d
                y_indices[e2] = b
                break
            att += 1
        clust_coeffs.append(directed_clustering(R))
        char_paths.append(directed_char_path_length(R))
    return clust_coeffs, char_paths
            

def directed_char_path_length(A):
    """Compute the char. path length for a DiGraph.

    This matches charpath in the Sporns Matlab toolbox.

    Parameters
    ----------
    A : NetworkX DiGraph or Numpy matrix
      Graph or adjacency matrix
    """
    if isinstance(A, nx.DiGraph):
        A = nx.adjacency_matrix(A)
    # First get the distance matrix D.
    D = np.eye(A.shape[0])
    n = 1
    n_path_matrix = copy.deepcopy(A) # We don't want to f with A.
    n_path_array = np.array(n_path_matrix)
    L = np.nan_to_num(n_path_array / n_path_array)
    try:
        nonzero = np.nonzero(L)[0][0]
    except IndexError:
        nonzero = False
    else:
        nonzero = True          # If index is 0, we still want True condition.
    while nonzero:
        D += n * L
        n += 1
        n_path_matrix *= A
        n_path_array = np.array(n_path_matrix)
        D_zeros = np.zeros(D.shape)
        D_zeros[D == 0] = 1
        L = np.nan_to_num(n_path_array / n_path_array) * D_zeros
        try:
            nonzero = np.nonzero(L)[0][0]
        except IndexError:
            nonzero = False
        else:
            nonzero = True
    D[D==0] = np.inf
    D -= np.eye(A.shape[0])
    # Now use D to compute the characteristic path length.
    return np.sum(np.sum(D[D!=np.inf]))/len(np.nonzero(D[D!=np.inf])[0])
    

def directed_clustering(g):
    """Compute the clustering coefficient for a DiGraph.

    Edges are considered binary (i.e., unweighted).  All directed triangles
    are used.

    G. Fagiolo, 2007, Physical Review E

    See also clustering_coef_bd in the Sporns Matlab toolbox.
    """
    if isinstance(A, nx.DiGraph):
        A = nx.adjacency_matrix(A)
    S = A + A.transpose()
    K = np.array(S.sum(axis=1)) # Make array for elementwise operations.
    cyc3 = np.array((S**3).diagonal() / 2.0).reshape(A.shape[0], 1)
    CYC3 = K*(K-1) - np.array(2*(A**2).diagonal()).reshape(A.shape[0], 1)
    # If there are zero possible 3-cycles, make the value in CYC3 Inf,
    # so that C = 0.  This is the definition of Rubinov & Sporns,
    # 2010, NeuroImage.
    CYC3[np.where(CYC3==0)] = np.inf
    C = cyc3/CYC3
    return C.mean()


def _compute_directed_closeness(path_lengths, node, num_nodes, all_nodes,
                                direction):
    lengths = 0.0
    for other in all_nodes:
        if node == other:
            continue
        if direction == 'in':
            try:
                lengths += path_lengths[other][node]
            except KeyError:
                return
        elif direction == 'out':
            try:
                lengths += path_lengths[node][other]
            except KeyError:
                return
        else:
            raise ValueError('invalid direction')
    return lengths / (num_nodes - 1)


def directed_closeness(g, direction='in'):
    """Calculate in- or out-closeness for nodes in g.

    Parameters
    ----------
    g : NetworkX DiGraph

    direction : string (optional)
      'in' or 'out'.

    Returns
    -------
    closeness : dict
      Dict mapping nodes to their in- or out-closeness value.
    """
    path_lengths = nx.shortest_path_length(g)
    all_nodes = g.nodes()
    num_nodes = g.number_of_nodes()
    closeness = {}
    for node in all_nodes:
        closeness[node] = _compute_directed_closeness(path_lengths, node,
                                                      num_nodes, all_nodes,
                                                      direction)
    return closeness
        

def compute_graph_of_unknowns(end):
    """Return the inverse of end.

    Because end contains known-absent and known-present edges, the
    returned graph will contain only those edges whose existence is
    unknown.

    Parameters
    ----------
    end : EndGraph
      EndGraph after translated edges have been added.

    Returns
    -------
    g : NetworkX DiGraph
      Graph with edges not in EndGraph.

    Notes
    -----
    end is not modified by this function; a new graph is returned.
    """
    u = nx.DiGraph()
    nodes = end.nodes()
    for source in nodes:
        for target in nodes:
            if source != target and not end.has_edge(source, target):
                u.add_edge(source, target)
    return u
    

def get_top_ten(measures, better='greater'):
    """Returns top ten nodes in order from best to worst.

    Parameters
    ----------
    measures : dict
      Mapping of nodes to values for a particular measure.

    better : string
      Complete this sentence using the word "greater" or "smaller": The
      nodes with the better scores in this dict are the ones with the
      _____ values.

    Returns
    -------
    top_ten : list
      The top ten nodes.

    Notes
    -----
    Nodes corresponding to 10 distinct scores (or as many as there are in
    the graph if there are fewer than 10) are returned; ties are placed in
    brackets.
    """
    top_ten = []
    best_to_worst = sorted(set(measures.values()))
    if better == 'greater':
        best_to_worst.reverse()
    while len(top_ten) < 10:
        try:
            next_best_score = best_to_worst.pop(0)
        except IndexError:
            return top_ten
        next_best_nodes = []
        for key, value in measures.iteritems():
            if value == next_best_score:
                next_best_nodes.append(key)
        if len(next_best_nodes) == 1:
            top_ten.append(next_best_nodes[0])
        else:
            top_ten.append(next_best_nodes)
    return top_ten


def strip_absent_edges(end):
    """Return graph with known-absent edges removed.

    Do not modify the graph passed as input; return a new graph.

    Parameters
    ----------
    end : EndGraph
      EndGraph after translated edges have been added.

    Returns
    -------
    g : NetworkX DiGraph
      Graph with those edges in EndGraph known to be present.  Edge
      attributes are not transferred to this graph.
    """
    g = nx.DiGraph()
    for source, target in end.edges_iter():
        source_ec, target_ec = end[source][target]['ECs']
        ns = ('N', 'Nc', 'Np', 'Nx')
        if source_ec not in ns and target_ec not in ns:
            g.add_edge(source, target)
    return g

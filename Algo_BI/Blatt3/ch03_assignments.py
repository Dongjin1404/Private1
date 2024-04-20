import math
import numpy
import random
import itertools

def generate_string():
    return ''.join([random.choice('ACGT') for i in range(40)])

def sorted_kmers(genome, k):
    """Assignment 1a. 

    Arguments:
        genome (str): A nucleotide sequence.
        k (int): The k-mer length.

    Returns:
        [str, ...]: A sorted list of k-mers.
    """
    kmers = []
    for i in range(len(genome)-k+1):
        kmers.append(genome[i:i+k])
    
    return sorted(kmers)

def all_genomes(kmers):
    """Assignment 1b. 

    arguments:
        kmers (list): a sorted list of k-mers.

    returns:
        [str, ...]: a list/set/generator of all genomes that can be reconstructed from input k-mers.
    """
    def prefix(string):
        return string[:-1]
    
    def suffix(string):
        return string[1:]
    
    def overlap_graph(kmers):
        ograph = {kmer : [] for kmer in kmers}
        for read in kmers:
            for oread in kmers:
                if read != oread and suffix(read) == prefix(oread):
                    ograph[read].append(oread)
        return ograph

    def hamiltonian_paths(graph):
        def dfs(vertex, path):                          # is a recursive function that performs a depth-first search (DFS) on a graph
            if len(path) == len(graph):                 # current path that has been taken to reach the vertex 
                paths.append(path)                      # keeps track of all found Hamiltonian paths
                return                                  #  returns control to the calling function
            for neighbor in graph[vertex]:
                if neighbor not in path:                #avoids visiting the same vertex more than once in the same path
                    dfs(neighbor, path + [neighbor])    #recursive call with current neighbor as new vertex  
                                                        #and the current path plus the neighbor as the new path.
        paths = []
        for vertex in graph:
            dfs(vertex, [vertex])                       # path is represented as a list of vertices
        return paths

    def reconstruct_string(path):
        return path[0] + ''.join([kmer[-1] for kmer in path[1:]])
    
    # Construct the overlap graph
    graph = overlap_graph(kmers)

    # Find all Hamiltonian paths in the graph
    paths = hamiltonian_paths(graph)

    # Reconstruct the full-length strings from the paths
    genomes = [reconstruct_string(path) for path in paths]

    return genomes
#overlap_graph: iterates over all pairs of k-mers and checks if the suffix of one is equal to the prefix of the other; O(n^2 * k)
#each of the n k-mers, the algorithm checks all other n - 1 k-mers for overlap

#hamiltonian_paths: This recursive function iterates over all k-mers and performs a depth-first search (DFS) on the overlap graph; O(n!)
#In the worst case, the algorithm has to explore all possible permutations of the vertices
#These algorithms can become very slow even for relatively small inputs, as n! grows extremely fast with n.

#Finding Hamiltonian paths is a well-known NP-complete problem (nondeterministic polynomial time)
#easy to verify (given a solution, you can check its correctness in polynomial time), 
#but finding a solution is believed to be difficult and can take exponential time in the worst case.

def get_deBruijn_graph(kmers):
    """ Assignment 2a. 

    Arguments:
        kmers (list): A sorted list of k-mers.

    Returns:
        nodes ({nodes}): a set of nodes
        edges ([(node_from, node_to), ...]): a list of edges connecting nodes

    """
    def prefix(string):
        return string[:-1]
    
    def suffix(string):
        return string[1:]
    
    nodes = set()
    edges = []

    for kmer in kmers:
        nodes.add(prefix(kmer))
        nodes.add(suffix(kmer))
        edges.append((prefix(kmer), suffix(kmer)))

    return nodes, edges
    
def eulerian_path(graph):
    edge_graph = [(prefix, suffix) for prefix, suffixes in graph.items() for suffix in suffixes]
    start_vertex = next(vertex for vertex in graph if sum(edge[0] == vertex for edge in edge_graph) > sum(edge[1] == vertex for edge in edge_graph))
    #checks if the start node of the current edge (represented by edge[0]) is equal to the current vertex
    #generates a sequence of boolean values (True or False) indicating whether the start node of each edge is equal to the current vertex.
    # So it calculates the number of edges that start at the current vertex.

    path = []
    stack = [start_vertex]

    while stack:
        vertex = stack[-1]
        if graph[vertex]:
            stack.append(graph[vertex].pop())
        else:
            path.append(stack.pop())
        
    return path[::-1]

def get_genome_deBruijn(kmers):
    """ Assignment 2b: 

    arguments:
        kmers (list): a sorted list of k-mers.

    returns:
        str: a valid genomes that can be reconstructed from input k-mers.
    """
    nodes, edges = get_deBruijn_graph(kmers)
    graph = {node: [] for node in nodes}
    for edge in edges:                          #adds the end node of the current edge to the list of nodes 
        graph[edge[0]].append(edge[1])          #that can be reached from the start node of the current edge
    path = eulerian_path(graph)
    genome = path[0] + ''.join([kmer[-1] for kmer in path[1:]])
    return genome

def get_contigs(kmers):
    nodes, edges = get_deBruijn_graph(kmers)
    graph = {node: [] for node in nodes}        #dictionary stores the adjacency list representation of the De Bruijn graph;
    for edge in edges:                          # each key is a node, and the value associated is a list of all nodes adjacent
        graph[edge[0]].append(edge[1])          #populates dictionary with destination nodes; tuple (node_from, node_to)

    contigs = []
    for node in graph:# if the length of the list of outgoing edges (graph[node]) is not 1 
        if len(graph[node]) != 1 or len([edge for edge in edges if edge[1] == node]) != 1#or the length of the list of incoming edges is not 1
            for next_node in graph[node]:#loop over all nodes that are directly reachable from the current node
                contig = [node] #initializes a new contig with the current node
                while len(graph[next_node]) == 1 and len([edge for edge in edges if edge[1] == next_node]) == 1:
                    contig.append(next_node) #continues as long as the next_node is a 1-in-1-out node
                    next_node = graph[next_node][0]# moves to the next node in the path by following the edge from the current next_node
                contig.append(next_node)#[0] index to get the first node in the list of nodes that the current node (next_node) has an edge to.
                contigs.append(contig)

    return contigs

if __name__ == '__main__':
    # You can use this section to test some functions when you execute this
    # file directly.
    pass

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

def get_deBruijn_graph(kmers): 
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
    #creates list of all edges in deBrujin graph (as a tuple of two nodes)
    #graph is a dictionary where for each node (prefix) there is a list of reachable nodes (suffixes)
    start_vertex = next(vertex for vertex in graph if sum(edge[0] == vertex for edge in edge_graph) > sum(edge[1] == vertex for edge in edge_graph))
    #start vertex is the one that appears more times as a prefix than as a suffix in the edge graph
    #checks if the prefix of the edge is equal to the current vertex = counts how many times the vertex is a prefix
    #checks if the suffix of the edge is equal to the current vertex = counts how many times the vertex is a suffix
    #next() function returns the first vertex that appears more times as a prefix than as a suffix.
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
    for edge in edges:
        graph[edge[0]].append(edge[1])          # adjacency list representation of the de Bruijn graph
    path = eulerian_path(graph)
    genome = path[0] + ''.join([kmer[-1] for kmer in path[1:]])
    return genome

from graphviz import Digraph

def draw_deBruijn_graph(kmers, filename):
    nodes, edges = get_deBruijn_graph(kmers)
    graph = Digraph(format='png')
    for node in nodes:
        graph.node(node)
    for edge in edges:
        graph.edge(*edge)                   #syntax for argument unpacking. It takes a list or tuple and 'unpacks' it into positional arguments in a function call.
    graph.render(filename, view=False)      #saves the graph to a file with the given filename, image file won't be opened after it's created.

s1 = 'AGGCTAGTACGGACTTACGCACAACGCTTGCGGA'
s2 = 'AGGCTAGTACGGACTTGCGCACAAGGCTTGCGGA'

for k in [5, 6]:
    print(f"De Bruijn graph for s1 with k={k}:")
    kmers = sorted_kmers(s1, k)
    draw_deBruijn_graph(kmers, f"s1_k{k}")

    print(f"De Bruijn graph for s2 with k={k}:")
    kmers = sorted_kmers(s2, k)
    draw_deBruijn_graph(kmers, f"s2_k{k}")



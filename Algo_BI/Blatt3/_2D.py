import random

def generate_string(length):
    return ''.join([random.choice('ACGT') for i in range(length)])

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
#contigs : long, contiguous segments of the genome
#a contig is a set of overlapping DNA segments
    
#The reason you visualize the de Bruijn graph as soon as you get a string that is different 
#from the original string is to understand why the reconstructed string is different. 
#The de Bruijn graph and the highlighted contigs can give you insights into the structure of the genome 
#and the overlaps between the k-mers.
    
#When the reconstructed string is different from the original string, it means that 
#there are multiple valid paths through the de Bruijn graph (i.e., the graph is not a simple path).
#This can happen, for example, when there are repeated sequences in the genome. 
#By visualizing the graph and the contigs, you can see these multiple paths and 
#understand how they correspond to the structure of the genome.

import matplotlib.pyplot as plt

import networkx as nx

def visualize_deBruijn_graph(nodes, edges, contigs):
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=1500)

    colors = ['red', 'blue']  # Add more colors if there are more contigs
    for i, contig in enumerate(contigs):
        edge_list = [(contig[j], contig[j+1]) for j in range(len(contig)-1)]
        nx.draw_networkx_edges(G, pos, edgelist=edge_list, edge_color=colors[i], width=2)

    plt.savefig("_2D.png")

# Generate a random string of length 40
genome = generate_string(40)

# Start with k=10 and decrease it
for k in range(10, 0, -1):
    # Generate sorted k-mers
    kmers = sorted_kmers(genome, k)

    # Reconstruct the genome from the k-mers
    reconstructed_genome = get_genome_deBruijn(kmers)

    # If the reconstructed genome is different from the original genome, visualize the de Bruijn graph
    if reconstructed_genome != genome:
        print(f"Reconstructed genome is different from original genome at k={k}")
        print(f"Original genome: {genome}")
        print(f"Reconstructed genome: {reconstructed_genome}")

        nodes, edges = get_deBruijn_graph(kmers)
        contigs = get_contigs(kmers)
        print(contigs)
        visualize_deBruijn_graph(nodes, edges, contigs[:2])
        break


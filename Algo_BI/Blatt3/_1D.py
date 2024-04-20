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
            unused = list(graph.keys() - set(path))
            for neighbor in unused:
                if neighbor in graph[vertex]:           
                    dfs(neighbor, path + [neighbor])    #recursive call with current neighbor as new vertex  
                                                        #and the current path plus the neighbor as the new path.
        paths = []
        for vertex in graph:
            dfs(vertex, [vertex])                       # path is represented as a list of vertices
        return paths

def reconstruct_string(path):
        return path[0] + ''.join([kmer[-1] for kmer in path[1:]])

def all_genomes(kmers):
    """Assignment 1b. 

    arguments:
        kmers (list): a sorted list of k-mers.

    returns:
        [str, ...]: a list/set/generator of all genomes that can be reconstructed from input k-mers.
    """
    
    # Construct the overlap graph
    graph = overlap_graph(kmers)

    # Find all Hamiltonian paths in the graph
    paths = hamiltonian_paths(graph)

    # Reconstruct the full-length strings from the paths
    genomes = [reconstruct_string(path) for path in paths]

    return genomes

import matplotlib.pyplot as plt

# Generate a random string of length 40
genome = generate_string(40)

# Initialize lists to store k values and number of genomes
k_values = []
num_genomes = []

# Start with k=10 and decrease it
for k in range(10, 0, -1):
    # Generate sorted k-mers
    kmers = sorted_kmers(genome, k)

    # Find all genomes that can be reconstructed from the k-mers
    genomes = all_genomes(kmers)

    # Add the k value and the number of genomes to the lists
    k_values.append(k)
    num_genomes.append(len(genomes))
    
    # Print a message indicating that the loop has been processed
    print(f"Processed loop for k={k}")

# Plot the number of possibilities to assemble a string of the correct length against k
plt.plot(k_values, num_genomes)
plt.xlabel('k')
plt.ylabel('Number of genomes')
plt.title('Number of possible genomes for different k values')
plt.savefig("_1D.png")
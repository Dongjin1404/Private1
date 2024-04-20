def hammingdistance(g1, g2):
    if len(g1) != len(g2):
        raise ValueError("Both sequences must have the same length")
    return sum(n1 != n2 for n1, n2 in zip(g1, g2))


def kd_neighbors(kmer, d):
    if d == 0:
        return [kmer]
    if len(kmer) == 1:
        return ["A", "C", "G", "T"]
    neighbors = []
    suffix_neighbors = kd_neighbors(kmer[1:], d)
    for str in suffix_neighbors:
        if hammingdistance(kmer[1:], str) < d:
            for n in "ACGT":
                neighbors.append(n + str)
        else:
            neighbors.append(kmer[0] + str)
    return neighbors


# The termination condition is when the Hamming distance d is 0 or when the length of the kmer is 1.


import time
import matplotlib.pyplot as plt

k_values = list(range(1, 11))  # Test k values from 1 to 10
d = 4  # Keep d constant
runtimes = []

for k in k_values:
    kmer = 'A' * k  # Create a k-mer of length k
    start_time = time.time()
    kd_neighbors(kmer, d)
    end_time = time.time()
    runtimes.append(end_time - start_time)

plt.plot(k_values, runtimes)
plt.xlabel('k')
plt.ylabel('Runtime (seconds)')
plt.title('Runtime of kd_neighbors as a function of k')
plt.savefig('_1A_1.png')

# The runtime of kd_neighbors increases exponentially with k. This is because the number of recursive calls increases exponentially with k.

d_values = list(range(1, 8))  # Test d values from 1 to 10
k = 10  # Keep k constant
runtimes = []

for d in d_values:
    kmer = 'A' * k  # Create a k-mer of length k
    start_time = time.time()
    kd_neighbors(kmer, d)
    end_time = time.time()
    runtimes.append(end_time - start_time)

plt.plot(d_values, runtimes)
plt.xlabel('d')
plt.ylabel('Runtime (seconds)')
plt.title('Runtime of kd_neighbors as a function of d')
plt.savefig('_1A_2.png')

# The runtime of kd_neighbors increases exponentially with d. This is because the number of recursive calls increases exponentially with d.
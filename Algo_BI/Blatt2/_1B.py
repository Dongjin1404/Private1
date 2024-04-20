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
        if hammingdistance(kmer[1:], str) <= d:
            for n in "ACGT":
                neighbors.append(n + str)
        else:
            neighbors.append(kmer[0] + str)
    return neighbors

def motif_enumeration(dna, k, d):
    motifs = set()
    for str in dna:
        for i in range(len(str)-k+1):
            neighbors_set = kd_neighbors(str[i:i+k], d)
            for x in neighbors_set:
                if all(any(hammingdistance(x, s[j:j+k]) <= d for j in range(len(s)-k+1)) for s in dna):
                    motifs.add(x)
    return motifs






#
import time
import matplotlib.pyplot as plt

d_values = list(range(1, 7))  # Test d values from 1 to 7
k = 8  # Keep k constant
dna = ['ATGCATCGGTC', 'GATGAATTGCC', 'CGTAACTTCTG']  # Example DNA sequences
runtimes = []

for d in d_values:
    start_time = time.time()
    motif_enumeration(dna, k, d)
    end_time = time.time()
    runtimes.append(end_time - start_time)

plt.plot(d_values, runtimes)
plt.xlabel('d')
plt.ylabel('Runtime (seconds)')
plt.title('Runtime of motif_enumeration as a function of d')
plt.savefig('_1B_2.png')
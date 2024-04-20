import math
import numpy as np
import random
import itertools


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


def motif_score(motifs, pcount=True):
    score = 0
    k = len(motifs[0])
    for i in range(k):
        column = [motif[i] for motif in motifs]  # For each position, it creates a list of the nucleotides at that position in all motifs
        nucleotide_counts = [column.count(nucleotide) for nucleotide in "ACGT"]
        if pcount:
            nucleotide_counts = [count + 1 for count in nucleotide_counts]  # Add pseudocounts (+1 to each count)
        score += k - max(nucleotide_counts)  # number of non-majority nucleotides
    return score


def motif_profile(motifs, pcount=True):
    k = len(motifs[0])
    profile = np.zeros((4, k))  # 4xk matrix
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            profile["ACGT".index(nucleotide), i] += 1
    if pcount:
        profile += 1  # Add pseudocounts
    for i in range(4):
        for j in range(k):
            profile[i][j] /= len(motifs) + 4 if pcount else len(motifs)
    return profile


def profile_entropy(profile):
    # Calculate the entropy for each position and sum them
    entropy = -np.sum(profile * np.log2(profile + np.finfo(float).eps)) # Add epsilon to avoid log(0) which is undefined
    # If a probability is zero, it becomes np.finfo(float).eps before taking the log, so there's no error or undefined value.
    return entropy


def get_motifs_from_profile(profile, dna):
    motifs = []
    k = len(profile[0])
    for sequence in dna:
        opt_score = -1
        opt_motif = sequence[:k]
        for i in range(len(sequence) - k + 1):
            motif = sequence[i:i+k]
            score = 1
            for j, nucleotide in enumerate(motif): #j is the index of the nucleotide in the motif string + nucleotide
                score *= profile["ACGT".index(nucleotide)][j] # index in the string "ACGT", which gives the row index, j = column index
            if score > opt_score:
                opt_score = score
                opt_motif = motif
        motifs.append(opt_motif)
    return motifs


def randomized_motifsearch(dna, k):
    Motifs = []
    for string in dna:
        start = random.randint(0, len(string)-k) # randomly select a k-mer from each string in dna
        Motifs.append(string[start:start+k])
       
    BestMotifs = Motifs.copy() # creates a copy of the list Motifs

    while True: # infinite loop, keeps iterating until it finds the best motifs 
        Profile = motif_profile(Motifs) # creates a profile matrix from the motifs

        Motifs = get_motifs_from_profile(Profile, dna) # creates the best motifs from the profile matrix

        # Compare scores of motifs and best motifs
        if motif_score(Motifs) < motif_score(BestMotifs):
            BestMotifs = Motifs.copy()
        else:
            return BestMotifs


def read_file(filename):
    with open(filename, 'r') as file:
        dna = [line.strip() for line in file if line.strip()] # Remove empty lines with if line.strip()
    return dna

if __name__ == "__main__":
    for filename in ['implanted_8mers.txt', 'implanted_10mers.txt', 'implanted_15mers.txt']:
        dna = read_file(filename)
        k = int(filename.split('_')[1].split('mers')[0])  # Extract k from filename
        best_motifs = randomized_motifsearch(dna, k)
        print(f"Best motifs for {filename}: {best_motifs}")

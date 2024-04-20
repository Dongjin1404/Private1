import numpy as np

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

motifs = ['AAAACC', 
          'GACAAA', 
          'AGACAA', 
          'CAGGAA',
          'ACAAGG']

implanted_8mers = ['TCAGATGGGCTGCTTGCAGGTTTCTTTTGACCCGGGCCGCGCGGTAGCACCCCGAGGACGCTATCTGAGGGATAC',
                   'CAGTTAAGCAGTTTCCTTGCTCGCCGGAACCCGACTCGCAAGCCAACCGTTTAGTGGAAGGAACCCACCGCGGGG', 
                   'CCCCGGAGGACCTACGACTTGGGAGGGTAATGCACTTTTTCCACACAGTCGGTCCAAAGTCGGGAAGAACTTACC',
                   'CTGCCAGTCCAAGGTATCGTATAGACCGTCAGTAATTGTATACCACGGGGGGTGCAGCTCTGTGCCGGTCGGTGC',
                   'GTTCAGAACAATAGCCCCGGTGCATACGCGTGAAATAAATCCGCTAGCTTCGGTTTTTGCCGAGCAGCTCTTAAT']

# Convert motifs to profile
profile = motif_profile(motifs)

# Get best scoring motifs from each DNA sequence
best_motifs = get_motifs_from_profile(profile, implanted_8mers)

print(best_motifs)

def motif_probability(motif, profile):
    probability = 1
    for j, nucleotide in enumerate(motif):
        probability *= profile["ACGT".index(nucleotide)][j]
    return probability

for motif in best_motifs:
    print(f"The probability of motif {motif} is {motif_probability(motif, profile)}")
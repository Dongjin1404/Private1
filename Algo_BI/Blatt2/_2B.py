def motif_profile(motifs, pcount=True):
    """Task 2b: Make a motif profile.

    Note: it is recommended to use numpy for working with matrices.

    Arguments:
        motifs (list): A list of DNA sequences.
        pcount (bool, optional): Add pseudocounts if True.

    Returns:
        profile (matrix): A matrix where entry ij is the probability of
            nucleotide i being at position j.
    """
    return [[]]

def motif_profile(motifs, pcount=True):
    k = len(motifs[0])
    profile = np.zeros((4, k)) # 4xk matrix
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            profile["ACGT".index(nucleotide), i] += 1
    if pcount:
        profile += 1  # Add pseudocounts
    for i in range(4):
        for j in range(k):
            profile[i][j] /= len(motifs) + 4 if pcount else len(motifs) 
    return profile
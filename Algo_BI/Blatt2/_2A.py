def motif_score(motifs):
    """Task 2a: Score motifs by diff to consensus.

    Arguments:
        motifs (list): A list of DNA sequences.
        pcount (bool, optional): Add pseudocounts if True.

    Returns:
        int: The score.
    """

    score = 0
    k = len(motifs[0])
    for i in range(k):
        column = [motif[i] for motif in motifs] # For each position, it creates a list of the nucleotides at that position in all motifs
        nucleotide_counts = [column.count(nucleotide) for nucleotide in "ACGT"]
        score += len(motifs) - max(nucleotide_counts) # number of non-majority nucleotides
    return score

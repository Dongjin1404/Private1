import random

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
 
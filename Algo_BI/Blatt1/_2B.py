def kmers(text, k):
    no_kmers = []
    for i in range(0, len(text), k):  # incrementing by k at each step.
        no_kmers.append(text[i : i + k])
    return no_kmers

def calculate_hamming_distance(genome1, genome2):
    if len(genome1) != len(genome2):
        raise ValueError(
            "Genomes must be the same length for calculating Hamming distance"
        )

    return sum(nucl1 != nucl2 for nucl1, nucl2 in zip(genome1, genome2))


with open("fna/query.fna", "r") as file:
    query_content = file.read()


def get_genomes(filename):
    with open(filename, "r") as f:
        genomes = f.read().split(">")[1:]  # split the file content by '>', and ignore the first empty string
        return genomes


genomes = get_genomes("fna/genomes.fna")
genomes1_content, genomes2_content, genomes3_content = genomes[:3]  # The syntax for slicing is list[start:stop]


for i, genome in enumerate(
    [genomes1_content, genomes2_content, genomes3_content], start=1
):
    HDs = []
    for percentage in range(10, 101, 10):
        portion_query = query_content[: int(len(query_content) * (percentage / 100))]
        portion_genome = genome[: int(len(genome) * (percentage / 100))]
        HD = calculate_hamming_distance(portion_query, portion_genome)
        HDs.append((percentage, HD))
    print(f"Genome {i}: {HDs}")
# look for the point where the Hamming distance starts to increase for one genome
# and decrease for another.

# starts off being most similar to Genome 1
# For Genome 1, the Hamming distance starts to increase more rapidly after the 40% mark, going from 1267 to 1575.
# For Genome 2, the Hamming distance starts to increase less rapidly after the 40% mark, going from 1649 to 1889.
# suggests that the switch from one genome to the other likely occurs around the 40% mark of the query genome.

# The accuracy of this estimation depends on:
# Analysis of the Hamming distance are done at 10% intervals. So the margin of error is 10%.
# If the genomes are very similar, the Hamming distance might not increase/decrease rapidly enough to be able to distinguish between them.
# If the genomes are very complex and contain many repeated sequences, it might be difficult to identify the switch point.


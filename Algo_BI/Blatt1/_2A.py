with open("fna/query.fna", "r") as file:
    query_content = file.read()


def get_genomes(filename):
    with open(filename, "r") as f:
        genomes = f.read().split(">")[1:]  # split the file content by '>', and ignore the first empty string
        return genomes


genomes = get_genomes("fna/genomes.fna")
genomes1_content, genomes2_content, genomes3_content = genomes[:3] # The syntax for slicing is list[start:stop]


# Calculate the Hamming distance between two genomes
def calculate_hamming_distance(genome1, genome2):
    if len(genome1) != len(genome2):
        raise ValueError("Genomes must be the same length for calculating Hamming distance")

    return sum(nucl1 != nucl2 for nucl1, nucl2 in zip(genome1, genome2))


def main():
    # Calculate the Hamming distance between query and genomes
    hamming_distance_q1 = calculate_hamming_distance(query_content, genomes1_content)
    hamming_distance_q2 = calculate_hamming_distance(query_content, genomes2_content)
    hamming_distance_q3 = calculate_hamming_distance(query_content, genomes3_content)

    print(f"Hamming distance between query and genome1: {hamming_distance_q1}")
    print(f"Hamming distance between query and genome2: {hamming_distance_q2}")
    print(f"Hamming distance between query and genome3: {hamming_distance_q3}")

main()
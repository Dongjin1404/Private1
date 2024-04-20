from _2B import kmers
from _2C import patternCount
from _2I import approximatePatternCount

with open("fna/query.fna", "r") as file:
    query_content = file.read()


def get_genomes(filename):
    with open(filename, "r") as f:
        genomes = f.read().split(">")[
            1:
        ]  # split the file content by '>', and ignore the first empty string
        return genomes


genomes = get_genomes("fna/genomes.fna")
genomes1_content, genomes2_content, genomes3_content = genomes[
    :3
]  # The syntax for slicing is list[start:stop]

tenpercent_query = query_content[-int(len(query_content) * 0.1) :]  # defaults to 0
tenpercent_tenmers_query = kmers(tenpercent_query, 10)

# dictionaries to store the counts
tenmers_counts_genomes1 = {}
tenmers_counts_genomes2 = {}
tenmers_counts_genomes3 = {}

# Add 10-mers and their counts to the dictionaries
for tenmer in tenpercent_tenmers_query:
    tenmers_counts_genomes1[tenmer] = approximatePatternCount(genomes1_content, tenmer, 1)
    tenmers_counts_genomes2[tenmer] = approximatePatternCount(genomes2_content, tenmer, 1)
    tenmers_counts_genomes3[tenmer] = approximatePatternCount(genomes3_content, tenmer, 1)

# Print the counts of each 10-mer
print(" Each 10-mer in Genome 1:")
for tenmer, count in tenmers_counts_genomes1.items():  # key-value pairs
    print(f"{tenmer}: {count}")

print("Each 10-mer in Genome 2:")
for tenmer, count in tenmers_counts_genomes2.items():
    print(f"{tenmer}: {count}")

print("Each 10-mer in Genome 3:")
for tenmer, count in tenmers_counts_genomes3.items():
    print(f"{tenmer}: {count}")

# Print the total number of 10-mers per genome
print("Total number of 10-mers in Genome 1:", sum(tenmers_counts_genomes1.values()))
print("Total number of 10-mers in Genome 2:", sum(tenmers_counts_genomes2.values()))
print("Total number of 10-mers in Genome 3:", sum(tenmers_counts_genomes3.values()))

# Most similar to the query genome, according to this measure, is Genome 1: 5

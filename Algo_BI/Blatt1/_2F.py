import time
import matplotlib.pyplot as plt
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

# lists to store the percentages and runtimes
percentages = []
runtimes = []

# Apply the k-mer procedure for each percentage from 20% to 100%
for percentage in range(20, 101, 10):
    start_time = time.time()

    portion_query = query_content[: int(len(query_content) * (percentage / 100))]

    # Use the kmer_procedure from _2D.py
    portion_tenmers_query = kmers(portion_query, 10)

    # Initialize dictionaries to store the counts
    tenmers_counts_genomes1 = {}
    tenmers_counts_genomes2 = {}
    tenmers_counts_genomes3 = {}

    # Add the 10-mers and their counts to the dictionaries
    for tenmer in portion_tenmers_query:
        tenmers_counts_genomes1[tenmer] = patternCount(genomes1_content, tenmer)
        tenmers_counts_genomes2[tenmer] = patternCount(genomes2_content, tenmer)
        tenmers_counts_genomes3[tenmer] = patternCount(genomes3_content, tenmer)

    # calculate total number of 10-mers per genome
    total_tenmers_genomes1 = sum(tenmers_counts_genomes1.values())
    total_tenmers_genomes2 = sum(tenmers_counts_genomes2.values())
    total_tenmers_genomes3 = sum(tenmers_counts_genomes3.values())

    end_time = time.time()
    runtime = end_time - start_time

    percentages.append(percentage)
    runtimes.append(runtime)
# Plot the runtimes
plt.plot(percentages, runtimes)
plt.xlabel("Percentage of Query Genome")
plt.ylabel("Runtime (seconds)")
plt.title("Runtime of K-mer Procedure")

# Save the plot
plt.savefig("runtime2F.png")

# ONLY KMERS PROCEDURE:
# runtime of your implementation depends on the size of the input: length of the query genome and the number of k-mers searched for.
# As the size of the input increases, the runtime will also increase.
# Whether the implementation is fast enough to process the entire dataset depends on the size of the dataset and the computational resources

# for an operation that needs to be run frequently or where results are needed quickly, a runtime of a few minutes to an hour might be acceptable.
# For larger datasets or less urgent scenarios, a runtime of several hours to a day might be acceptable.

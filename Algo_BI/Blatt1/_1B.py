import matplotlib.pyplot as plt

def calculate_gc_skew(genome):
    skew = [0]  # Starting value
    for nucleotide in genome:
        if nucleotide == "G":
            skew.append(skew[-1] + 1)  # list[-1] the last element of a List.
        elif nucleotide == "C":
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew


with open("fna/halomonas.fna", "r") as file:
    next(file)  # skip the first line
    genome = file.read()

# Calculate GC-Skew
gc_skew = calculate_gc_skew(genome)

# Diagram of the GC-Skew
plt.plot(range(len(gc_skew)), gc_skew)  # x-values, y-values
plt.xlabel("Position")
plt.ylabel("GC Skew")
plt.savefig("GC-Skew1B.png")

# the oriC at the GC SKEW minimum at the position where the reverse half-strand ends and the forward half-strand begins
# The accuracy of this method depends on several factors, including the quality of the sequence data and
# the specific characteristics of the organism from which the genome comes.
# this is a simplified method, additional analyses and experiments are often required

# Find the minimum GC skew
min_gc_skew = min(gc_skew)
origin_of_replication = gc_skew.index(min_gc_skew)

print(f"The origin of replication is approximately located at {origin_of_replication}.")

import matplotlib.pyplot as plt


def calculate_gt_skew(genome):
    skew = [0]  # Starting value
    for nucleotide in genome:
        if nucleotide == "G":
            skew.append(skew[-1] + 1)  # list[-1] the last element of a List.
        elif nucleotide == "T":
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew


with open("fna/halomonas.fna", "r") as file:
    next(file)  # skip the first line
    genome = file.read()

# Calculate GC-Skew
gt_skew = calculate_gt_skew(genome)

# Erstellen Sie das Diagramm
plt.plot(range(len(gt_skew)), gt_skew)
plt.xlabel("Position")
plt.ylabel("GT Skew")
plt.savefig("GT-Skew1C.png")

# the oriC at the minimum at the position where the reverse half-strand ends and the forward half-strand begins
# The accuracy of this method depends on several factors, including the quality of the sequence data and
# the specific characteristics of the organism from which the genome comes.
# this is a simplified method, additional analyses and experiments are often required

# A steadily rising GT skew plot could be due to a number of factors:
# some organisms have a higher proportion of G to T bases
# also be due to the specific region of the genome, may have a higher concentration of G bases due to the presence of certain genes
# (can be influenced by the direction of DNA replication:
# leading strand of DNA (which is synthesized continuously during replication) has a higher proportion of G to T bases
# positive GT skew in the direction of replication and a negative GT skew in the opposite direction.)

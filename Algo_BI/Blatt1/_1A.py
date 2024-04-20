import matplotlib.pyplot as plt

with open("fna/halomonas.fna", "r") as file:
    next(file)  # skip the first line
    halomonas_content = file.read()

# Count the occurrence of each nucleotide
nucleotides = ["A", "C", "G", "T"]
counts = [halomonas_content.count(nucleotide) for nucleotide in nucleotides]

# Print out each nucleotide count
for nucleotide, count in zip(nucleotides, counts):
    print(f"Nucleotide {nucleotide}: {count}")

# Create a bar chart
plt.bar(nucleotides, counts)

# Label the axes
plt.xlabel("Nucleotides")
plt.ylabel("Count")

plt.savefig("bar_chart1A.png")  # Save diagramm as PNG file

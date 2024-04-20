# The annotations attribute of a SeqRecord object in Biopython is a dictionary that contains additional information about the sequence.
# This information is fetched from the database using the Entrez API, and the keys of this dictionary represent the fields of the metadata.

# In the provided code, record.annotations.keys() is used to print the keys of the annotations dictionary,
# which are the fields of the metadata for the fetched sequence.

from Bio import Entrez
from Bio import SeqIO


# Set the email (replace with your email)
Entrez.email = "dongjin@gmx.net"


# Search the Protein database for the first 3 datasets
handle = Entrez.esearch(db="protein", term="cancer", retmax=3)
record = Entrez.read(handle)
ids = record["IdList"]

# For each ID in the list
for id in ids:
    # Fetch the metadata of the dataset
    handle = Entrez.efetch(db="protein", id=id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")

    print(record.annotations.keys())

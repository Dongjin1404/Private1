import mysql.connector
from Bio import Entrez
from Bio import SeqIO
import datetime
import os

# Connect to the MySQL server
db = mysql.connector.connect(
    user="root",
    password="Biologie2511835812",
    host="localhost",
    database="softwareprojekt",
)

# Set the email (replace with your email)
Entrez.email = "dongjin@gmx.net"

# Create a cursor
cursor = db.cursor()

# Check if the metadata_Nucleotides table exists, if not, create it
cursor.execute(
    """
CREATE TABLE IF NOT EXISTS metadata_Nucleotides (
    accession VARCHAR(255),
    molecule_type VARCHAR(255),
    topology VARCHAR(255),
    data_file_division VARCHAR(255),
    date VARCHAR(255),
    sequence_version INT,
    keywords TEXT,
    source VARCHAR(255),
    organism VARCHAR(255),
    taxonomy TEXT,
    ref_text TEXT,
    comment TEXT,
    timestamp DATETIME,
    PRIMARY KEY (accession, timestamp)
)
"""
)

# Create the directory for the downloaded data
if not os.path.exists("Nucleotide_download_data"):
    os.makedirs("Nucleotide_download_data")

# Search the Nucleotide database for the first 10 datasets
handle = Entrez.esearch(db="nucleotide", term="cancer", retmax=10)
record = Entrez.read(handle)
ids = record["IdList"]

# For each ID in the list
for id in ids:
    # Fetch the metadata of the dataset
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    # Insert the metadata into the table
    cursor.execute(
        """
    INSERT INTO metadata_Nucleotides (
        accession,
        molecule_type,
        topology,
        data_file_division,
        date,
        sequence_version,
        keywords,
        source,
        organism,
        taxonomy,
        ref_text,
        comment,
        timestamp
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """,
        (
            record.id,
            record.annotations.get("molecule_type"),
            record.annotations.get("topology"),
            record.annotations.get("data_file_division"),
            record.annotations.get("date"),
            record.annotations.get("sequence_version"),
            ", ".join(record.annotations.get("keywords")),
            record.annotations.get("source"),
            record.annotations.get("organism"),
            ", ".join(record.annotations.get("taxonomy")),
            str(record.annotations.get("references")),
            record.annotations.get("comment"),
            datetime.datetime.now(),  # Current timestamp
        ),
    )

    # Commit the transaction
    db.commit()

    # Get the current timestamp and format it as a string
    timestamp_str = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Write the record to a local file
    with open(f"Nucleotide_download_data/{record.id}_{timestamp_str}.gb", "w") as f:
        SeqIO.write(record, f, "genbank")

# Close the database connection
db.close()

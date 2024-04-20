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

# Set the email address
Entrez.email = "dongjin@gmx.net"

# Create a cursor
cursor = db.cursor()

# Check if the metadata_protein_sequences table exists, if not, create it
cursor.execute("SHOW TABLES LIKE 'metadata_protein_sequences'")
result = cursor.fetchone()
if not result:
    cursor.execute(
        """
    CREATE TABLE metadata_protein_sequences (
        accession VARCHAR(255),
        name VARCHAR(255),
        organism VARCHAR(255),
        sequence_length INT,
        source VARCHAR(255),
        taxonomy TEXT,
        sequence TEXT,
        features TEXT,
        ref_text TEXT,
        timestamp DATETIME,
        PRIMARY KEY (accession, timestamp)
    )
    """
    )

# Create the directory if it doesn't exist
if not os.path.exists("protein_sequence_data"):
    os.makedirs("protein_sequence_data")

# Search the Protein database for the first 10 datasets
handle = Entrez.esearch(db="protein", term="cancer", retmax=10)
record = Entrez.read(handle)
ids = record["IdList"]

# For each ID in the list
for id in ids:
    # Fetch the metadata of the dataset
    handle = Entrez.efetch(db="protein", id=id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")

    # Insert the metadata into the table
    cursor.execute(
        """
    INSERT INTO metadata_protein_sequences (
        accession,
        name,
        organism,
        sequence_length,
        source,
        taxonomy,
        sequence,
        features,
        ref_text,
        timestamp
    ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """,
        (
            record.id,
            record.name,
            record.annotations.get("organism"),
            len(record),
            str(record.annotations.get("source")),
            ", ".join(record.annotations.get("taxonomy")),
            str(record.seq),
            str(record.features),
            str(record.annotations.get("references")),
            datetime.datetime.now(),  # Current timestamp
        ),
    )

    # Commit the transaction
    db.commit()

    # Get the current timestamp and format it as a string
    timestamp_str = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # Write the record to a local file
    with open(f"protein_sequence_data/{record.id}_{timestamp_str}.gb", "w") as f:
        SeqIO.write(record, f, "genbank")

# Close the database connection
db.close()

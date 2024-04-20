# E-utilities is simply another way to search PubMed and the other NCBI databases.
# E-utilities is an API, or Application Programming Interface: a set of rules, protocols, and tools for building software and applications.

import mysql.connector
from Bio import Entrez
from Bio import Medline
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

# Check if the metadata_pubmed table exists, if not, create it
cursor.execute("SHOW TABLES LIKE 'metadata_pubmed'")
result = cursor.fetchone()
if not result:
    cursor.execute(
        """
    CREATE TABLE metadata_pubmed (
        PMID VARCHAR(255),
        Title TEXT,
        Authors TEXT,
        Abstract TEXT,
        PubDate VARCHAR(255),
        Keywords TEXT,
        Timestamp DATETIME,
        PRIMARY KEY (PMID, Timestamp)
    )
    """
    )

# Create a directory for the downloaded files
if not os.path.exists("Pubmed_data"):
    os.makedirs("Pubmed_data")

# Search the PubMed database for "covid" in the title, limit to 10 results
handle = Entrez.esearch(db="pubmed", term="covid[TI]", retmax=10)
record = Entrez.read(handle)

# Get the list of Ids
idlist = record["IdList"]

# For each Id in the list
for id in idlist:
    # Fetch the article
    handle = Entrez.efetch(db="pubmed", id=id, rettype="medline", retmode="text")
    try:
        record = Medline.read(handle)
    except StopIteration:
        print(f"No record found for PubMed ID {id}")
        continue

    # Insert the metadata into the MySQL database
    query = """
    INSERT INTO metadata_pubmed (PMID, Title, Authors, Abstract, PubDate, Keywords, Timestamp)
    VALUES (%s, %s, %s, %s, %s, %s, %s)
    """
    values = (
        id,
        record.get("TI", ""),
        ", ".join(record.get("AU", "")),
        record.get("AB", ""),
        record.get("DP", ""),
        ", ".join(record.get("OT", "")),  # OT field contains keywords
        datetime.datetime.now(),  # Current timestamp
    )
    cursor.execute(query, values)
    db.commit()  # Commit the transaction

    # Get the current timestamp and format it as a string
    timestamp_str = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Save the record locally
    with open(f"Pubmed_data/{id}_{timestamp_str}.txt", "w") as f:
        f.write(str(record))

# Close the database connection
db.close()

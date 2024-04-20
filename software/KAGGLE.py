import os
import mysql.connector
from kaggle.api.kaggle_api_extended import KaggleApi
from datetime import datetime


# Connect to the MySQL server
cnx = mysql.connector.connect(
    user="root",
    password="Biologie2511835812",
    host="localhost",
    database="softwareprojekt",
)

# create a cursor
cursor = cnx.cursor()


# create the metadata table
# We're using TEXT because the tags are stored as a comma-separated string, which can be quite long.
# Check if the metadata table exists, if not, create it

cursor.execute("SHOW TABLES LIKE 'metadata'")
result = cursor.fetchone()
if not result:
    cursor.execute(
        """
    CREATE TABLE IF NOT EXISTS metadata (
        file_name VARCHAR(255),
        last_updated TIMESTAMP,
        website VARCHAR(255),
        subtitle VARCHAR(255),
        tags TEXT,
        total_bytes BIGINT,
        download_count INT,
        description TEXT,
        Timestamp DATETIME,
        PRIMARY KEY (file_name, timestamp)
    )
    """
    )

# authenticate API
api = KaggleApi()
api.authenticate()

# Create a directory for the downloaded files
if not os.path.exists("Kaggle_downloaded_files"):
    os.makedirs("Kaggle_downloaded_files")

timestamp_str = datetime.now().strftime("%Y-%m-%d_%H-%M-%S.%f")
# get the first 20 pages of datasets....2000 datasets....
for page in range(1, 21):
    datasets = api.datasets_list(page=page, max_size=100)
    for ds in datasets:
        if "covid" in ds["title"].lower():
            # get the dataset info
            # Note that the tags field is a list of dictionaries,
            # so we're using a list comprehension to convert it to a string of tag names separated by commas.
            file_name = ds["id"]
            last_updated = ds["lastUpdated"]
            # Parse the ISO 8601 datetime string
            last_updated = datetime.strptime(ds["lastUpdated"], "%Y-%m-%dT%H:%M:%S.%fZ")
            # Format the datetime in the format MySQL expects
            last_updated = last_updated.strftime("%Y-%m-%d %H:%M:%S")
            website = ds["url"]
            subtitle = ds["subtitle"]
            tags = ", ".join(
                tag["name"] for tag in ds["tags"]
            )  # convert list of tags to a string
            total_bytes = ds["totalBytes"]
            download_count = ds["downloadCount"]
            description = ds["description"]
            timestamp = datetime.now()
            # insert the dataset info into the metadata table
            cursor.execute(
                """
            INSERT INTO metadata (file_name, last_updated, website, subtitle, tags, total_bytes, download_count, description, timestamp)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
            """,
                (
                    file_name,
                    last_updated,
                    website,
                    subtitle,
                    tags,
                    total_bytes,
                    download_count,
                    description,
                    datetime.now(),
                ),
            )
            # commit the changes
            cnx.commit()

            # Download the dataset files locally
            print(f"Starting download for {file_name}...")
            try:
                api.dataset_download_files(
                    ds["ref"],
                    path=f"Kaggle_downloaded_files/{file_name}_{timestamp_str}",
                )
            except BaseException as e:  # Catch all types of exceptions
                print(f"An error occurred while downloading the dataset: {e}")
            else:
                print(f"Download completed for {file_name}.")

# close the connection
cursor.close()
cnx.close()

from kaggle.api.kaggle_api_extended import KaggleApi

# authenticate API
api = KaggleApi()
api.authenticate()

# get the first 3 pages of datasets....300 datasets
for page in range(1, 4):
    datasets = api.datasets_list(page=page, max_size=100)
    for ds in datasets:
        # only print the title if it contains 'covid'
        if "covid" in ds["title"].lower():
            # print fields of metadata for each dataset
            for key in ds.keys():
                print(f"{key}, {ds[key]}")
            print("--------------------------------------")

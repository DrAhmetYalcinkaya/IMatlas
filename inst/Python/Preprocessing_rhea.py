import requests
import pandas as pd
import numpy as np
import gzip
import io
import zipfile
import re
from itertools import combinations

class RheaDB:
    """
    This class handles all RheaDB processing. This includes retrieving data and 
    parsing each reaction to get a metabolite-metabolite format.
    """
    def __init__(self, options, log):
        self.options = options
        self.log = log
        self.chebi_df = []
            
    def parse_rhea(self, mapping):
        """
        Main method call to start parsing of the RheaDB. Data is retrieved and 
        filtered for duplicates before writing the dataframe to a file.
        """
        print("Parsing RheaDB")
        with gzip.GzipFile(fileobj=io.BytesIO(self.retrieve_rhea()), mode = "r") as gzipped:
            all = self.parse_gzip(gzipped, mapping)
            df = pd.DataFrame(all, columns = ["From", "To"]).drop_duplicates()
            df.to_csv(f"{self.options['folder']}/Metabolite-metabolite.csv", index = False)
            self.log.write(f"Wrote {len(df.index)} Metabolite-metabolite interactions to a file")

    def retrieve_rhea(self):
        """
        This method downloads retrieves data from RheaDB.
        """
        url = "http://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz"
        return requests.get(url, stream=True).content
    
    def parse_gzip(self, content, mapping):
        """
        This method will parse the gzip provided by RheaDB.
        """
        all = []
        while True:
            line = content.readline().decode("utf-8")
            chebis = list(set(re.findall("(CHEBI:\d+)", line)))
            if len(chebis) > 0:
                hmdbs = list(mapping.iloc[np.where(mapping.index.isin(chebis))]["ID"])
                all += combinations(hmdbs, 2)
            
            if line == "":
                break
        return all
    
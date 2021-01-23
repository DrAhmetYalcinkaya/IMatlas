import requests
import pandas as pd
import numpy as np
import gzip
import io
import zipfile

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
        with gzip.GzipFile(fileobj=io.BytesIO(self.retrieve_rhea()), mode = "r") as gzipped:
            all = self.parse_gzip(gzipped, mapping)
            df = pd.DataFrame(all, columns = ["From", "To"]).drop_duplicates()
            df = df[df['From'] != df['To']]
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
        parse = False
        all = []
        while True:
            line = content.readline().decode("utf-8")
            if "ENTRY" in line:
                parse = True

            if "EQUATION" in line and parse:
                all += self.parse_equation(line, mapping)
                parse = False
            
            if line == "":
                break
        return all
    
    def parse_equation(self, line, mapping):
        """
        This method parses a line in a RheaDB file. Given the line and the 
        mapping towards HMDB, combinations between metabolites are made 
        to obtain metabolite-metabolite interactions.
        """
        all = []
        chebis = " ".join(line.split()[1:]).split("=")
        try:
            b = mapping.loc[self.parse_reaction(chebis[0]), "ID"]
            a = mapping.loc[self.parse_reaction(chebis[1]), "ID"]
            all = [(id1, id2) for id1 in b for id2 in a if str(id1) != "nan" and str(id2) != "nan"]
        except: # some or all keys not found in hmdb
            pass
        return all

    def parse_reaction(self, reactants_string):
        """
        This method parses one part of a reaction to obtain all reactants.
        """
        all = set()
        for reac in reactants_string.split("+"):
            reac = "CHEBI" + reac.split("CHEBI", 1)[1].rstrip("<").strip()
            if "," in reac:
                reac = reac.split(",")[0]
            all.add(reac)
        return list(all)
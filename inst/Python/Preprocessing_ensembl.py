from concurrent.futures import ThreadPoolExecutor
import requests
import pandas as pd
import json
import numpy as np
from time import sleep
import logging
from tqdm import tqdm
import os

class Ensembl:
    """
    This class handles the conversion of Ensembl transcripts to Ensembl proteins as used by StringDB.
    """
    def __init__(self, options):
        logging.info("Connecting with Ensembl DB")

        self.server = "http://rest.ensembl.org"
        self.ext = "/lookup/id/Translation"
        self.headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
        self.responses = []
        self.options = options

    def get_response(self, ids):
        """
        Returns a 2-column list to be converted into a Pandas dataframe 
        containing the transcript and protein identifier
        """
        try:
            req = requests.post(self.server + self.ext, headers=self.headers, params={"expand": True, "format": "condensed"}, 
                                data='{ "ids" : %s }' % str(ids).replace("'", '"'))
            self.json = req.json()
            return [[id, self.json[id]["Translation"]["id"]] for id in self.json.keys() if self.json.get(id) is not None and "Translation" in self.json[id].keys()]
        except:
            logging.warning(f"Encounterd a problem with status code {req.status_code}. Trying again in 60 seconds.")
            sleep(60)
            print(ids)
            self.get_response(ids)

    def set_transcripts(self, df):
        """

        """
        df["Ensembl transcript"] = df["Ensembl transcript"].str.findall(r'(ENST\d+)')
        df = df.explode("Ensembl transcript")
        transcripts = list(df["Ensembl transcript"])
        mapping_loc = str(os.path.realpath(__file__)).replace(r'Python\Preprocessing_ensembl.py', "") + "extdata\Ensembl_Mapping.csv"
        mapping = pd.read_csv(mapping_loc)
        known = list(mapping["Ensembl transcript"])
        self.transcripts = list(set(transcripts) - set(known))
        return df

    def get_transcript_mapping(self):
        """

        """
        logging.info("Start mapping transcripts to Ensembl Proteins")
        mapping_loc = str(os.path.realpath(__file__)).replace(r'Python\Preprocessing_ensembl.py', "") + "extdata\Ensembl_Mapping.csv"
        n = 500
        for i in range(0, len(self.transcripts), n):
            results = self.get_response(self.transcripts[i:i+n])
            conv = pd.DataFrame(results, columns = ["Ensembl transcript", "StringDB"])
            conv.to_csv(mapping_loc, mode="a", header=False, index = False)
            
        logging.info("Completed mapping to Ensembl Proteins")
        conv = pd.read_csv(mapping_loc).drop_duplicates()
        return conv.set_index("Ensembl transcript", drop = False)

    def get_ensembl_proteins(self, df, mapping):
        """

        """
        df.set_index("Ensembl transcript", inplace = True)
        df = df.drop_duplicates()
        df = df.loc[~df.index.duplicated()]
        mapping = mapping.drop_duplicates()
        indexes = mapping.index.intersection(df.index)
        df = pd.concat([df.loc[indexes], mapping.loc[indexes]["StringDB"]], axis=1, sort=True)
        df.reset_index(drop = True, inplace = True)
        logging.info(f"Found {len(df.index)} proteins that were mapped to StringDB identifiers")
        return df

from concurrent.futures import ThreadPoolExecutor
import requests
import pandas as pd
import json
import numpy as np
from time import sleep
import logging

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
        self.submitted_ids = 0
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
            self.submitted_ids += len(ids)
            logging.info(f"Mapped {self.submitted_ids} out of {len(self.transcripts)}", end="\r")
            return [[id, self.json[id]["Translation"]["id"]] for id in self.json.keys() if self.json.get(id) is not None and "Translation" in self.json[id].keys()]
        except:
            logging.warning(f"Encounterd a problem with status code {req.status_code}. Trying again in 60 seconds.")
            sleep(60)
            self.get_response(ids)

    def set_transcripts(self, df):
        """

        """
        df["Ensembl transcript"] = df["Ensembl transcript"].str.findall(r'(ENST\d+)')
        df = df.explode("Ensembl transcript")
        self.transcripts = list(df["Ensembl transcript"])
        return df

    def get_transcript_mapping(self):
        """

        """
        logging.info("Start mapping transcripts to Ensembl Proteins")
        n = 1000
        with ThreadPoolExecutor() as ex:
            futures = [ex.submit(self.get_response, self.transcripts[i:i+n]) 
                        for i in range(0, len(self.transcripts), n)]

            for f in futures:
                self.responses += f.result()

        conv = pd.DataFrame(self.responses, columns = ["Ensembl transcript", "StringDB"])
        conv.set_index("Ensembl transcript", inplace = True, drop = False)
        logging.info("Completed mapping to Ensembl Proteins")
        return conv

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
import pandas as pd
import numpy as np
import logging

class StringDB:
    """ 
    This class handles all processing regarding StringDB, 
    including retrieving data and parsing it.
    """
    def __init__(self, options):
        self.options = options
        self.string_df = []
        logging.info("Start protein-protein interactions")

    def get_stringdb_df(self):
        """
        Retrieving the data from StringDB. It preprocesses the data to obtain correct 
        column names and stripping the organism identifier from each protein.
        """
        url = "https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz"
        self.string_df = pd.read_csv(url, compression="gzip", sep = " ")
        self.string_df.reset_index(drop = True, inplace = True)
        self.string_df["protein1"] = self.string_df["protein1"].str.lstrip("9606.")
        self.string_df["protein2"] = self.string_df["protein2"].str.lstrip("9606.")
        self.string_df.reset_index(drop = True, inplace = True)
        self.string_df.columns = ["From", "To", "Confidence"]
        logging.info(f"Found {len(self.string_df.index)} unique interactions in StringDB")

    def map_string_to_uniprot(self, df):
        """
        This methods maps the dataframe obtained from StringDB to Uniprot 
        entries using the Uniprot dataframe.
        """
        df["StringDB"] = df["StringDB"].str.rstrip(";")
        df.set_index("StringDB", inplace=True)
        self.string_df = self.string_df.loc[(self.in_df(df, "From")) & (self.in_df(df, "To"))]
        self.string_df["From"] = list(df.loc[self.string_df["From"]]["Entry"])
        self.string_df["To"] = list(df.loc[self.string_df["To"]]["Entry"])
        logging.info(f"Mapped {len(self.string_df.index)} interactions containing proteins to Uniprot identifiers")
        return df
    
    def in_df(self, df, column):
        """
        Method that returns which values in the given column are 
        in the index of the given dataframe.
        """
        return self.string_df[column].isin(df.index)

    def extract_protein_interactions(self, df, rows):
        """
        This method extracts protein-protein interactions that 
        """
        df = self.map_string_to_uniprot(df)
        ids = list(df.iloc[rows]["Entry"])
        rows = list(set(np.where(self.string_df["From"].isin(ids))[0]) & set(np.where(self.string_df["To"].isin(ids))[0]))
        logging.info(f"Found {len(rows)} interactions between proteins of interest")
        self.string_df.iloc[rows].to_csv(f"{self.options['folder']}/Protein-protein.csv", index = False)
        self.extract_neighbours(ids)
        logging.info("Done Protein-Protein Interactions")

    def extract_neighbours(self, ids):
        """
        This method find neighbours of the given identifiers and extracts their interactions 
        to a file. It repeats this process twice so two rounds of neighbours are found.
        """
        files = ["Protein-protein_1.csv", "Protein-protein_2.csv"]
        logging.info(f"Start {len(files)} iterations of intermediate protein-protein interactions")

        for f in files:
            rows = list(np.where(self.string_df["From"].isin(ids))[0]) 
            logging.info(f"Found {len(rows)} interactions between proteins of interest")
            self.string_df.iloc[rows].to_csv(f"{self.options['folder']}/{f}", index = False)
            ids = list(set(self.string_df.iloc[rows][["From", "To"]].stack().values))
            logging.info(f"Found {len(ids)} neighbouring proteins")
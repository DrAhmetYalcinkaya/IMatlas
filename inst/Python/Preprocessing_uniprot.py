import pandas as pd
import numpy as np
from Preprocessing_ensembl import Ensembl
from Preprocessing_stringdb import StringDB

class Uniprot:
    """ 
    The Uniprot class handles everything (uniprot) protein related. 
    This includes retrieving data from Uniprot, mapping transcripts 
    to proteins using Ensembl, extracting metadata and extracting 
    protein-protein interactions from StringDB. 
    """
    def __init__(self, options, log):
        self.options = options
        self.log = log
        self.df = []
        self.indexes = []

    def retrieve_uniprot_df(self):
        """
        This mmethod retrieves a dataframe from Uniprot. It includes columns: 
        'Entry', 'Protein names', 'Cofactors', 'EC numbers', 'Transporter proteins', 
        'GO', and 'Ensembl transcripts'. Since only proteins that have a Ensembl mapping
        can be used, the dataframe is filtered for those without mapping.
        """
        self.df = pd.read_csv("https://www.uniprot.org/uniprot/?query=*&format=tab&columns=id,protein%20names,comment(COFACTOR),ec,database(TCDB)," + \
        "go-id,database(Ensembl)&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes", sep="\t")
        self.log.write(f"Found {len(self.df.index)} proteins in the Uniprot database\n")
        self.df = self.df.dropna(subset = ["Ensembl transcript"])
        self.log.write(f"Found {len(self.df.index)} proteins with a known transcript\n")
        self.df.reset_index(drop = True, inplace = True)

    def map_transcripts_to_proteins(self):
        """
        Here, using the Ensembl class, transcripts in the Uniprot dataframe are mapped 
        to Ensembl proteins for usage in StringDB protein-protein interactions. 
        """
        ensembl_db = Ensembl(self.options, self.log)
        self.df = ensembl_db.set_transcripts(self.df)
        mapping = ensembl_db.get_transcript_mapping()
        self.df = ensembl_db.get_ensembl_proteins(self.df, mapping)

    def parse_metadata(self, chebi_mapping):
        """
        This method calls other methods that extract metadata 
        obtained from Uniprot.
        """
        self.extract_transporter()
        self.extract_protein_names() 
        self.extract_ec_numbers() 
        self.extract_cofactor(chebi_mapping)

    def extract_transporter(self):
        """
        This method extracts transporter proteins from Uniprot.
        """
        self.log.write("\nStart extracting transporter proteins\n")
        df = self.df.dropna(subset=["Cross-reference (TCDB)"])
        df = df[["Entry", "Cross-reference (TCDB)"]]
        df.columns = ["ID", "Transporter"]
        df["Transporter"] = df["Transporter"].str.strip(";")
        df.to_csv(f"{self.options['folder']}/Protein_transporter.csv", index = False)
        self.log.write(f"Found {len(df.index)} transporter proteins\n")
        print("Done Transporters")


    def extract_protein_names(self):
        """
        This method extracts protein names from Uniprot.
        """
        df = self.df[["Entry", "Protein names"]]
        names = df["Protein names"].str.split("(", expand = False).str[0]
        synonyms = df["Protein names"].str.split("(", expand = False).str[1].str.rstrip(") ")
        synonyms[synonyms.isnull()] = names[synonyms.isnull()]
        df = pd.concat([df["Entry"], names, synonyms], axis = 1)
        df.columns = ["ID", "Name", "Synonym"]
        df.to_csv(f"{self.options['folder']}/Protein_names.csv", index = False)
        self.log.write(f"written {len(df.index)} protein names to its file\n")
        print("Done Protein Names")

    def extract_ec_numbers(self):
        """
        This method extracts enzyme (EC) proteins from Uniprot.
        """
        df = self.df.dropna(subset=["EC number"])
        column = df["EC number"].str.split(";", expand = False).str[0]
        df = pd.concat([df["Entry"], column], axis = 1)
        df.columns = ["ID", "Number"]
        df.to_csv(f"{self.options['folder']}/Ec_numbers.csv", index = False)
        self.log.write(f"Found {len(df.index)} proteins with a EC number\n")
        print("Done EC Numbers")

    def extract_cofactor(self, conv):
        """
        This method extracts cofactor proteins from Uniprot and metabolites from HMDB.
        It uses ChEBI (provided by uniprot) as a mapping between uniprot and HMDB identifiers.        
        """
        conv.set_index("chebi_id", drop = False, inplace = True)
        df = self.df.dropna(subset=["Cofactor"])
        df.reset_index(drop = True, inplace = True)
        column = df["Cofactor"].str.split("ChEBI:", expand = False).str[1]
        column = column.str.split(";", expand = False).str[0].str.strip()
        df = pd.concat([df["Entry"], column], axis = 1)
        df = df.dropna(subset=["Cofactor"])
        df.reset_index(drop = True, inplace = True)
        cof = conv.loc[df.Cofactor.loc[df.Cofactor.isin(conv.index)], "ID"]
        cof.reset_index(drop = True, inplace = True)
        df = pd.concat([df["Entry"], cof], axis = 1)
        df.columns = ["ID", "Cofactor"]
        df["Cofactor"].replace('', np.nan, inplace=True)
        df.dropna(subset=["Cofactor"], inplace=True)
        df.to_csv(f"{self.options['folder']}/Cofactors.csv", index = False)
        self.log.write(f"Found {len(df.index)} protein-metabolite cofactors combinations\n")

    def parse_protein_interactions(self, gos):
        """
        Here, the StringDB class is used to obtain a StringDB interactions dataframe
        and extracts proteins that are associated with Gene Ontologies that are 
        descendants of the original provided Gene Ontology.
        """
        string_db = StringDB(self.options, self.log)
        string_db.get_stringdb_df()
        self.df["Gene ontology IDs"] = self.df["Gene ontology IDs"].str.split(";", expand = False)
        go_series = self.df["Gene ontology IDs"].explode()
        indexes = go_series.iloc[np.where(go_series.isin(gos))].index      
        string_db.extract_protein_interactions(self.df, indexes)        
        return self.df.iloc[indexes]
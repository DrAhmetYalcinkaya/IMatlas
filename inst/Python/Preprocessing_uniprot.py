import pandas as pd
import numpy as np
import sys
import csv
from Preprocessing_ensembl import Ensembl
from Preprocessing_stringdb import StringDB
import logging

class Uniprot:
    """ 
    The Uniprot class handles everything (uniprot) protein related. 
    This includes retrieving data from Uniprot, mapping transcripts 
    to proteins using Ensembl, extracting metadata and extracting 
    protein-protein interactions from StringDB. 
    """
    def __init__(self, options):
        self.options = options
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
        "go-id,database(Ensembl),database(Reactome)&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes", sep="\t")
        logging.info(f"Found {len(self.df.index)} proteins in the Uniprot database")
        self.df = self.df.dropna(subset = ["Ensembl transcript"])
        logging.info(f"Found {len(self.df.index)} proteins with a known transcript")
        self.df.reset_index(drop = True, inplace = True)
        self.extract_protein_names() 

    def map_transcripts_to_proteins(self):
        """
        Here, using the Ensembl class, transcripts in the Uniprot dataframe are mapped 
        to Ensembl proteins for usage in StringDB protein-protein interactions. 
        """
        ensembl_db = Ensembl(self.options)
        self.df = ensembl_db.set_transcripts(self.df)
        mapping = ensembl_db.get_transcript_mapping()
        self.df = ensembl_db.get_ensembl_proteins(self.df, mapping)

    def parse_metadata(self, chebi_mapping):
        """
        This method calls other methods that extract metadata 
        obtained from Uniprot.
        """
        self.extract_transporter()
        self.extract_ec_numbers() 
        self.extract_reactome()
        self.extract_cofactor(chebi_mapping)

    def extract_reactome(self):
        df = self.df.dropna(subset=["Cross-reference (Reactome)"])
        df = df[["Entry", "Cross-reference (Reactome)"]]
        df["Cross-reference (Reactome)"] = df["Cross-reference (Reactome)"].str.findall(r'(R-HSA-\d+)')
        df = df.explode("Cross-reference (Reactome)")  
        df.columns = ["ID", "Pathway"]
        df.to_csv(f"{self.options['folder']}/Protein_reactome.csv", index = False)
        logging.info(f"Found {len(df.index)} protein pathways")


    def extract_transporter(self):
        """
        This method extracts transporter proteins from Uniprot.
        """
        logging.info("\nStart extracting transporter proteins")
        df = self.df.dropna(subset=["Cross-reference (TCDB)"])
        df = df[["Entry", "Cross-reference (TCDB)"]]
        df.columns = ["ID", "Transporter"]
        df["Transporter"] = df["Transporter"].str.findall(r'(\d+\.[A-Z]+\.\d+\.\d+\.\d+)')
        df = df.explode("Transporter")
        df.to_csv(f"{self.options['folder']}/Protein_transporter.csv", index = False)
        logging.info(f"Found {len(df.index)} transporter proteins")


    def extract_protein_names(self):
        """
        This method extracts protein names from Uniprot.
        """
        df = self.df[["Entry", "Protein names"]]
        
        names = self.df["Protein names"].str.findall(r'.*? \(')
        synonyms = names.str[1].str.rstrip(") (")
        names = names.str[0].str.rstrip(" (")
        names[names.isnull()] = df[names.isnull()]["Protein names"]
        synonyms[synonyms.isnull()] = names[synonyms.isnull()]
        df = pd.concat([df["Entry"], names, synonyms], axis = 1)
        df.columns = ["ID", "Name", "Synonym"]
        df.to_csv(f"{self.options['folder']}/Protein_names.csv", index = False, quoting=csv.QUOTE_ALL)
        logging.info(f"written {len(df.index)} protein names to its file")

    def extract_ec_numbers(self):
        """
        This method extracts enzyme (EC) proteins from Uniprot.
        """
        df = self.df.dropna(subset=["EC number"])
        column = df["EC number"].str.findall(r'(\d+.\d+.\d+.\d+)')
        df = pd.concat([df["Entry"], column], axis = 1)
        df.columns = ["ID", "Number"]
        df = df.explode("Number")
        df.dropna(subset=["Number"], inplace=True)
        df.to_csv(f"{self.options['folder']}/Ec_numbers.csv", index = False)
        logging.info(f"Found {len(df.index)} proteins with a EC number")

    def extract_cofactor(self, conv):
        """
        This method extracts cofactor proteins from Uniprot and metabolites from HMDB.
        It uses ChEBI (provided by uniprot) as a mapping between uniprot and HMDB identifiers.        
        """
        conv.set_index("chebi_id", drop = False, inplace = True)
        df = self.df.dropna(subset=["Cofactor"])
        df["Cofactor"] = df["Cofactor"].str.findall(r'(CHEBI:\d+)')
        df = df.explode("Cofactor")
        df.reset_index(drop = True, inplace = True)
        cof = conv.loc[df.Cofactor.loc[df.Cofactor.isin(conv.index)], "ID"]
        cof.reset_index(drop = True, inplace = True)
        df = pd.concat([df["Entry"], cof], axis = 1)
        df.columns = ["ID", "Cofactor"]
        df["Cofactor"].replace('', np.nan, inplace=True)
        df.dropna(subset=["Cofactor"], inplace=True)
        df.to_csv(f"{self.options['folder']}/Cofactors.csv", index = False)
        logging.info(f"Found {len(df.index)} protein-metabolite cofactors combinations")

    def parse_protein_interactions(self, gos):
        """
        Here, the StringDB class is used to obtain a StringDB interactions dataframe
        and extracts proteins that are associated with Gene Ontologies that are 
        descendants of the original provided Gene Ontology.
        """
        string_db = StringDB(self.options)
        string_db.get_stringdb_df()
        self.df["Gene ontology IDs"] = self.df["Gene ontology IDs"].str.findall(r'(GO:\d+)')
        go_series = self.df["Gene ontology IDs"].explode()
        indexes = go_series.iloc[np.where(go_series.isin(gos))].index      
        string_db.extract_protein_interactions(self.df, indexes)        
        return self.df.iloc[indexes]
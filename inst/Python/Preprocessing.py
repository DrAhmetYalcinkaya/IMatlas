#!/usr/bin/env python

import pip

def import_or_install(packages):
    for package, to_install in packages.items():
        try:
            __import__(package)
        except ImportError:
            pip.main(['install', to_install])  

packages = {"yaml": "pyYaml", "pandas": "pandas", "numpy": "numpy", "zipfile": "zipfile",
    "json": "json", "gzip": "gzip", "requests": "requests", "pathlib": "pathlib",
    "xml": "xml", "re": "re", "concurrent": "concurrent", "io": "io", "itertools": "itertools"
}
import_or_install(packages)

import os
import sys
import yaml
from pathlib import Path
from Preprocessing_hmdb import HMDB
from Preprocessing_go import GODB
from Preprocessing_uniprot import Uniprot
from Preprocessing_rhea import RheaDB
from Preprocessing_stringdb import StringDB

def main(config_path):
    """ 
    This function describes the preprocessing and includes data from
    HMDB, RheaDB, GODB, Uniprot, Ensembl and StringDB. First, a data folder
    is made if it's non-existing and a log file is made. HMDB is parsed, and
    afterwards RheaDB for metabolite-metabolite interactions. GO is queried
    for the given GO identifier to obtain its descendants. These are used
    to extract proteins from Uniprot, which includes mapping Ensembl transcripts
    to Ensembl proteins in order to parse StringDB correctly.
    All data obtained is written to the given data folder.  
    """
    with open(config_path) as file: 
        options = yaml.load(file, Loader=yaml.FullLoader)
        Path(options["folder"]).mkdir(parents=True, exist_ok=True)
    print("Opened Config file.")
    options['folder'] = options['folder'].rstrip("/")

    log = open(f"{options['folder']}/Log_preprocessing.txt", "w", buffering=1)

    hmdb_db = HMDB(options, log)
    rhea_db = RheaDB(options, log)
    go_db = GODB(options, log)
    uniprot = Uniprot(options, log)
    string_db = StringDB(options, log)

    #hmdb_db.parse_hmdb()
    #rhea_db.parse_rhea(hmdb_db.get_chebi_mapping())
    gos = go_db.get_descendants()
    df = uniprot.retrieve_uniprot_df()
    entries = go_db.extract_gos(df) 

    #uniprot.parse_metadata(hmdb_db.get_chebi_mapping())
    df = uniprot.map_transcripts_to_proteins()
    string_db.get_stringdb_df()
    string_db.extract_protein_interactions(df, entries)  
    log.close()

if __name__ == "__main__":
    config_path = "config.yaml"
    if len(sys.argv) > 1:
        config_path = sys.argv[1].lstrip("'").rstrip("'").lstrip('"').rstrip('"')
    main(config_path)
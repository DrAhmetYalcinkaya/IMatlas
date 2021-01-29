def import_or_install(packages):
    for package, to_install in packages.items():
        try:
            __import__(package)
        except ImportError:
            pip.main(['install', to_install]) 

packages = {"yaml": "pyYaml", "pandas": "pandas", "json": "json","requests": "requests", 
            "pathlib": "pathlib", "concurrent": "concurrent", "time": "time", "random": "random", 
            "math": "math", "sys": "sys"}
import_or_install(packages)

import pandas as pd
import requests
import json
import time
import random
import yaml
from math import ceil
from concurrent.futures import ThreadPoolExecutor
import sys
from pathlib import Path

class TextMining:
    """
    This class handles all text mining related searches. 
    """
    def __init__(self, options, log):
        self.options = options
        self.log = log
        self.processed = 0
        self.total = 0

    def request_query(self, search, size, cursor):
        """
        Here are queries for a term determined. The fields in which to search are:
        Abstract, Methods, Results, Tables, Figures and Supplements. Other than that,
        The paper must contain either 'patient' or 'human' in order to limit the results
        to human-related studies. Lastly, the paper must be of the type 'research article'.
        This query has been made with the query builder of EuropePMC.
        """
        query = '(ABSTRACT:"{0}" OR RESULTS:"{0}" OR METHODS:"{0}" OR TABLE:"{0}" OR SUPPL:"{0}" OR FIG:"{0}") AND ("patient" OR "human") AND PUB_TYPE:"Journal Article"'.format(search)
        url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={query}&synonym=true&resultType=idlist&pageSize={size}&cursorMark={cursor}&format=json'      
        print("Current progress: %.2f%%" % ((self.processed + 1) / self.total * 100), end = "\r")  
        return json.loads(requests.get(url).content)

    def search(self, term, n):
        """
        This method performs the search for a single term. Since results are paginated 
        and limited to 1000 results per page, a while-true construction has been made 
        to loop through all pages. For each search, the page identifiers are stored and
        returned once all pages have been processed.
        """
        ids = []
        cursor = "*"
        res = ""
        try:
            while True:
                res = self.request_query(term, 1000, cursor)
                hitcount = int(res["hitCount"])
                ids += [result["id"] for result in res["resultList"]["result"]]
                if hitcount == 0 or cursor == res["nextCursorMark"]:
                    break
                cursor = res["nextCursorMark"] 

        except Exception:
            time.sleep(30)

        self.processed += 1
        return ids

    def search_all(self, list_of_terms):
        """
        Wrapper method for a list of given terms to search for. 
        Returns a dictionary containing a term - pubmed identifiers construction.
        """
        dic = {}
        self.processed = 0
        self.total = len(list_of_terms)
        with ThreadPoolExecutor() as ex:
            futures = [ex.submit(self.search, term, n) 
                        for n, term in enumerate(list_of_terms)]
            for i, f in enumerate(futures):
                ids = f.result()
                    
                if len(ids) > 0:
                    term = list_of_terms[i]
                    dic[term] = ids
            

        return dic


class EBI:
    """
    This class handles all requests that are Gene Ontology (GO) related. 
    """
    def __init__(self, options, log):
        self.options = options
        self.log = log
        self.n = 50
        self.base_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"

    def get_request(self, params, endpoint=""):
        """
        Performs a request to QuickGO using parameters and a endpoint.
        Returns a dictionary of the results.
        """
        url = f"{self.base_url}/{params}/{endpoint}"
        return dict(json.loads(requests.get(url).content))["results"]

    def get_descendants(self, go_id):
        """
        Returns a list of the descendants of a given GO identifier. 
        """
        dic = self.get_request(params = go_id, endpoint = "descendants")[0]
        gos = list(set(dic["descendants"])) + [go_id]
        return gos

    def get_go_names(self, identifiers):
        """
        Returns a dictionary of names : identifiers
        of the given identifiers.
        """
        total = {}
        for i in range(0, len(identifiers), self.n):
            params = ",".join(identifiers[i : i + self.n])
            for r in self.get_request(params = params):
                total[r["name"]] = r["id"]
        return total
    
    def get_ancestors(self, to_search, allowed, names):
        """
        Returns the ancestors for the given GO names. 
        Not all ancestors are part of the same tree as the given GO id 
        in the config, so only certain ancestors are allowed. 
        """
        total = {}
        id_names = {allowed[i]: names[i] for i in range(len(names))}
        to_search = [names[a] for a in to_search]
        for i in range(0, len(to_search), self.n):
            params = ",".join(to_search[i : i + self.n])
            for r in self.get_request(params, endpoint = "ancestors"):
                key = id_names[r["id"]]
                total[key] = [id_names[ances] for ances in r["ancestors"] if ances in allowed] + [key]
        return total


def write_counts(df, columns, prefix, suffix):
    """
    This function writes the value counts of the given column(s).
    This allows for further processing like statistics.
    """
    count_df = df.value_counts(subset = columns).to_frame()
    count_df.reset_index(inplace=True)
    count_df.columns = columns + ["Count"]
    count_df.to_csv(f"{prefix}_textmining_{suffix}.tsv", sep = "\t", index = False)

def get_expanded_df(ebi, df, options):
    """
    Returns a dataframe that includes ancestors
    """
    go_df = pd.read_csv(f"{options['folder']}/Go_names.csv")
    gos = list(go_df["GOID"])
    names = list(go_df["Name"])
    ancestors = ebi.get_ancestors(set(df["Gene Ontology"]), gos, names)
    df["Ancestors"] = ["\t".join(ancestors[go]) for go in df["Gene Ontology"]]

    df = df[["Ancestors", "Metabolite", "Paper ID"]]
    df = df.assign(Ancestors = df.Ancestors.str.split("\t")).explode("Ancestors")
    df.columns = ["Gene Ontology", "Metabolite", "Paper ID"]
    return df

def find_overlap(mets, gos):
    """
    Returns a dataframe that stores the overlap between paper identifiers
    of the GO and metabolite searches. An overlap indicates that both are 
    present in the paper, which indicates an association between the two.
    """
    total = []
    count = 0
    tot = len(gos) * len(mets)
    for go, go_ids in gos.items():
        for met, met_ids in mets.items():
            count += 1
            print(round(count * 100 / tot, 3), end="\r")
            to_loop = set(go_ids) & set(met_ids)
            total += [[go,met,inter] for inter in to_loop]
    return pd.DataFrame(total, columns=["Gene Ontology", "Metabolite", "Paper ID"])

def main(config_path):
    """
    Here the text mining script is called. It uses names obtained from preprocessing
    for both Gene Ontologies and Metabolites. For each, their names are searched using
    EuropePMC, which returns the paper IDs. Overlap between these papers indicate an
    association between two terms. The GO-Metabolite associations are counted and
    written to individual files. Afterwards, an 'extended' dataframe is created. 
    This is done to not only find direct associations, but also implied ones using
    the ancestors of the matched GO names. 
    """
    with open(config_path) as file: 
        options = yaml.load(file, Loader=yaml.FullLoader)
        Path(options["folder"]).mkdir(parents=True, exist_ok=True)
    print("Opened Config file.")
    log = open(f"{options['folder']}/Log_textmining.txt", "w", buffering=1)

    text_mining = TextMining(options, log)
    ebi = EBI(options, log)

    
    gos = list(set(pd.read_csv(f"{options['folder']}/Go_names.csv")["Name"]))
    mets = list(set(pd.read_csv(f"{options['folder']}/Metabolite_name.csv")["name"]))
    df = find_overlap(text_mining.search_all(mets), text_mining.search_all(gos))
    df.to_csv(f"raw_textmining_all.tsv", sep="\t", index = False)

    for prefix in ["raw", "expanded"]:
        write_counts(df, ["Paper ID"], prefix, "papers")
        write_counts(df, ["Metabolite"], prefix, "metabolites")
        write_counts(df, ["Gene Ontology"], prefix, "go")
        write_counts(df, ["Gene Ontology", "Metabolite"], prefix, "counts")
        if prefix == "raw":
            df = get_expanded_df(ebi, df, options)

    df.to_csv(f"extended_textmining_all.tsv", sep="\t", index = False)

if __name__ == "__main__":
    config_path = "config.yaml"
    if len(sys.argv) > 1:
        config_path = sys.argv[1].lstrip("'").rstrip("'").lstrip('"').rstrip('"')
    main(config_path)

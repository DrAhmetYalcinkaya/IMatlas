import json
import requests
import pandas as pd
import numpy as np

class GODB:
    """
    This class handles are Gene Ontology requests using the api from EBI.
    For a given Gene Ontology, descendants and ancestors can be requested
    so all relevant GOs of the given GO in the configuration can be retrieved.
    """
    def __init__(self, options, log):
        self.options = options
        self.go_id = options["GO_ID"]
        self.log = log
        self.gos = []
        self.log.write(f"Start log\n\nGO term: {self.go_id}\n")


    def get_descendants(self):
        """
        This method sets and returns the Gene Ontology identifiers that are 
        descendants of the given term. These GO identifiers are used as a 
        filter for selecting proteins.
        """
        url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}/descendants".format(self.go_id) 
        dic = dict(json.loads(requests.get(url).content))
        desc = list(set(dic["results"][0]["descendants"]))
        self.log.write(f"Number of descendants: {len(desc)}\n")
        print("Gathered Gene Ontologies")
        desc = [a.strip() for a in desc]
        self.gos = desc + [self.go_id]
        return self.gos

    def get_go_names(self):
        """
        This method returns a dictionary of an ID and its name using the api of EBI.
        """
        total = []
        n = 100
        for i in range(0, len(self.gos), n):
            url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}".format(",".join(self.gos[i:i+n])) 
            with requests.get(url) as req:
                res = dict(json.loads(req.content))["results"]
                total += [[r["id"], r["name"]] for r in res]

        return total


    def get_ancestors(self):
        """
        Here, the ancestors of the GO identifiers are retrieved. This method returns a 
        dictionary that contains a mapping between a ID and its ancestors
        """
        total = {}
        n = 100
        for i in range(0, len(self.gos), n):
            url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}/ancestors".format(",".join(self.gos[i:i+n]))
            res = dict(json.loads(requests.get(url).content))["results"]
            for r in res:
                if "ancestors" in r.keys():
                    total[r["id"]] = [ances for ances in r["ancestors"] if ances in self.gos] + [r["id"]]
        return total

    
    def extract_gos(self, df):
        """
        Here all protein-GO and GO-Name associations are written to a file. 
        For each protein, its associated GOs are requested, including the 
        ancestors of the given term. This is due to the hierarchial
        nature of Gene Ontology.
        """
        ancestors = self.get_ancestors()
        tot = []
        for entry, gos in zip(df["Entry"], df["Gene ontology IDs"]):
            all = []
            for go in gos:
                all += ancestors.get(go.strip(), [])
            for i in set(all):
                tot.append([entry, i])

        prot_gos = pd.DataFrame(tot, columns = ["ID", "GOID"])
        prot_gos.to_csv(f"{self.options['folder']}/Protein_gos.csv", index = False)
        pd.DataFrame(self.get_go_names(), columns = ["GOID", "Name"]).to_csv(f"{self.options['folder']}/Go_names.csv", index = False)
        print("Done Gene Ontologies")

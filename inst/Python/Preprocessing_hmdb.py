import re
import xml.etree.cElementTree as ET
import requests
import pandas as pd
import numpy as np
import io
import zipfile

class HMDB:
    """ This class chunks an xml file (zipped) into chunks, so that a DOM parser can be used
    It stores only fields by using regex of the fields given
    """
    def __init__(self, options, log):
        self.files = {}
        self.options = options
        self.log = log        
        self.mapping = []

    def set_variables(self, source, chunk_by, fields, exclude):
        """
        This method sets the variables needed to parse the HMDB XML file in chunks
        It primarly builds a regex expression to be searched in the file and
        sets up files so the results can be written to their respective file.
        """
        self.source = source
        self.chunk_start = "<{}".format(chunk_by)
        self.chunk_end = "</{}".format(chunk_by)
        self.regex = "|".join(["{}>".format(x) for x in fields])
        self.regex = "(" + self.chunk_start + "|" + self.chunk_end + "|" + self.regex + ")"
        self.set_files(fields, exclude)

    def set_files(self, fields, exclude):
        """
        Here, files are created that are being use to write the metadata to, except
        for all fields in 'exclude'Niet d
        """
        for field in fields:
            if field not in exclude:
                self.files[field] = open(f"{self.options['folder']}/Metabolite_{field}.csv", "w", encoding="utf-8")
                self.write(field, ["ID", field])
                self.log.write(f"Created file {self.options['folder']}/Metabolite_{field}.csv\n")

    def chunk(self):
        """
        This method creates a generator by yielding a ElementTree object for each chunk
        """
        print("Start parsing")
        with self.source.open("hmdb_metabolites.xml", "r") as f:
            next(f)
            next(f)
            data = []
            s = f.readline().decode("utf-8", "ignore")
            while s:
                data += [s] if re.search(self.regex, s) is not None else []
                if self.chunk_end in s:
                    yield ET.fromstring("".join(data))
                    data = []   
                s = f.readline().decode("utf-8", "ignore")            

    def write(self, field, vals):
        """
        This method writes the values of a given field to a file.
        """
        self.files[field].write('"' + '","'.join(vals) + '"\n')

    def close(self):
        """
        Here, all files are closed. Is used after processing the XML file.
        """
        for f in self.files.values():
            f.close()

    def get_chebi_mapping(self):
        """
        Returns a dataframe that contains a mapping between HMDB and ChEBI 
        identifiers in order to use with Uniprot
        """
        df = pd.DataFrame(self.mapping, columns = ["ID", "chebi_id"]).drop_duplicates()
        df.set_index("chebi_id", inplace=True, drop = False)
        df.to_csv(f"{self.options['folder']}/Metabolite-chebi.csv", index = False)
        return df

    
    def parse_hmdb(self): 
        """
        
        """
        print("Start downloading HMDB XML")
        url = 'https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip'
        r = requests.get(url, stream=True)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        count = 0
        with z.open("hmdb_metabolites.xml", "r") as source:
            self.set_variables(source = z, chunk_by = "metabolite", 
                        fields = ["accession", "chebi_id", "uniprot_id", 
                                  "class", "kegg_map_id", "super_class", 
                                  "biospecimen", "cellular", "name", "pathway", "term"],
                        exclude = ["term", "kegg_map_id", "chebi_id"])

            for n, chunk in enumerate(self.chunk()):
                print(n, end="\r")
                terms = {x.text for x in chunk.findall("term")}
                if "Biological role" in terms or "Naturally occurring process" in terms:
                    count += self.process_metabolite(chunk)
            self.close()  
            self.log.write(f"Extracted {count} metabolites out of {n} from HMDB\n")
            print("Done Metabolite data")

    
    def process_metabolite(self, chunk):
        """

        """
        accession = chunk.findtext("accession")
        
        self.write("name", [accession, chunk.findtext("name").replace('"', "'")])
        self.write("class", [accession, str(chunk.findtext("class"))])
        self.write("super_class", [accession, str(chunk.findtext("super_class"))])

        for tag in ["accession", "uniprot_id", "biospecimen", "cellular"]:
            for x in chunk.findall(tag):
                self.write(tag, [accession, x.text]) 

        chebi = chunk.findtext("chebi_id")
        if chebi is not None:
            self.mapping.append([accession, "CHEBI:" + chebi])

        chunks = chunk.findall("pathway//")
        for name, kegg in zip(chunks, chunks[1:]):
            if name.tag == "name" and kegg.tag == "kegg_map_id" and kegg.text is not None:
                self.write("pathway", [accession, name.text]) 
        return 1
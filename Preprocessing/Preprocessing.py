#!/usr/bin/env python

"""
This file is the preprocessing file for the immuno-metabolome atlas.
It takes all info from the 'config.yaml' file.  
"""

import re
import pandas as pd
import xml.etree.cElementTree as ET
import yaml
import os
from widediaper import R
import sys

# Open config file
print("Reading config file")
with open("/home/pascal/config.yaml") as file: # Change this line when run locally.
    options = yaml.load(file, Loader=yaml.FullLoader)

def regex(line):
    """
    This function uses regular expressions to find lines that contain data of
    interest. RE (regular expressions) were used instead of XML parsing because of XML parsing
    skipping certain lines due to an unknown bug. Because of the dictionary
    construction, the GO names need to be searched using a second re.search().
    :param line: current line of the file
    :return: When found, a name and its value are returned, otherwise None.
    """
    to_search = {"accession": "<{0}>(.+?)<\/{0}>",
                 "dbReference": '<.*type="STRING" id="(.+?)"\/>',
                 "GO": '<dbReference type="GO" id="(.+?)">',
                 "fullName": "<{0}.*>(.+?)<\/{0}>",
                 "name": '<{0}.*type="primary".*>(.+?)<\/{0}>',
                 "property": '<.*type="protein sequence ID" value="ENSP(.+?)"\/>'}

    if "</entry" in line:
        return "stop"

    for string, value in to_search.items():
        x = re.search(value.format(string), line)
        if x is not None:
            return string, x.group(1)
    x = re.search('<.*type="term" value="(.+?)"\/>', line)
    if x is not None:
        return "GO", x.group(1)

    return None


def add_names(total, dictionary, string):
    """
    This function adds names to the StringDB dictionary for later processing
    :param total: dictionary of all StringDBs so far
    :param dictionary: dictionary of data of current Protein
    :param string: StringDB identifier
    :return: total dictionary for further processing
    """

    if "name" in dictionary.keys():
        name = dictionary["name"][0]
        total[string]["name"] = name

    if "fullName" in dictionary.keys():
        name = dictionary["fullName"][0]
        total[string]["fullName"] = name

    return total


def parse_uniprot_xml(source):
    """
    This function parses the Uniprot XML file. First it retrieves all immune
    proteins by a line of R code. Some files are created to store the
    extracted data in. For each line, if the regex() function finds code,
    it is stored in an temporary dictionary
    :param source: The path where the XML file is located
    :return: dictionary which acts as a translation between StringDB and
    Uniprot info
    """
    print("Start parsing Uniprot XML\n")
    r = R()
    r.load_library("GO.db")
    immune_filter = r('c(GOBPOFFSPRING[["{}"]], "{}")'.format(
        options["GO_ID"], options["GO_ID"]))

    l = immune_filter.split(" ")
    immune_filter = []
    for go in l:
        for go1 in go.split("\r\n"):
            immune_filter.append(go1.lstrip('"').rstrip('"'))

    total = {}
    uniprot_cellular = {}
    immune_set = set()
    synonyms = open(options["folder"] + options["synonyms"], "w")
    go_cellular = open(options["folder"] + options["GO-cellular"], "w")
    go_cellular.write("Uniprot\tGO-id\tCellular\n")
    synonyms.write("from\tsynonym\n")
    go_translation = {}
    with open(options["folder"] + source) as fp:
        dictionary = {}
        for count, line in enumerate(fp):
            print("\rDone {}".format(count), end="\r")
            finds = regex(line)
            if finds is not None and finds != "stop":
                if finds[0] in dictionary:
                    dictionary[finds[0]].append(finds[1])
                else:
                    dictionary[finds[0]] = [finds[1]]
            if finds == "stop":
                uniprot = dictionary["accession"][0]
                stringDBs = set()
                if "dbReference" in dictionary.keys():
                    string = dictionary["dbReference"][0]
                    total[string] = {"Uniprot": uniprot}
                    total = add_names(total, dictionary, string)
                    stringDBs.add(string)

                if "property" in dictionary.keys():
                    for id in dictionary["property"]:
                        string = "9606.ENSP" + id
                        total[string] = {"Uniprot": uniprot}
                        total = add_names(total, dictionary, string)
                        stringDBs.add(string)


                if "GO" in dictionary.keys():
                    l = dictionary["GO"]

                    for go, name in list(zip(l, l[1:]))[::2]:
                        if go in immune_filter:
                            immune_set.add(uniprot)
                            name = name[2:]
                            go_translation[go] = name
                            for string in stringDBs:
                                if "GO-process" in total[string].keys():
                                    total[string]["GO-process"].add(go)
                                else:
                                    total[string]["GO-process"] = set([go])

                        elif name[0] == "C":
                             if uniprot not in uniprot_cellular.keys():
                                 uniprot_cellular[uniprot] = [name[2:]]
                             else:
                                 uniprot_cellular[uniprot].append(name[2:])

                if "name" in dictionary.keys():
                    name = dictionary["name"][0]
                    uniprot = dictionary["accession"][0]
                    synonyms.write("{}\t{}\n".format(uniprot, name))

                dictionary = {}

    for id, locations in uniprot_cellular.items():
        for l in locations:
            go_cellular.write("{}\tNone\t{}\n".format(id, l))

    go_cellular.close()
    synonyms.close()
    print("\nDone parsing Uniprot XML")
    return total, list(immune_set), go_translation

def write_interaction(f, t, output):
    """
    This function writes protein-protein interactions to the output given
    :param f: Dictionary of the 'from' protein
    :param t: Dictionary of the 'to' protein
    :param output: File object of the file to write to
    """
    from_id = f["Uniprot"]
    from_fullname = f["fullName"]
    to_id = t["Uniprot"]
    to_fullname = t["fullName"]
    output.write(
        "{}\t{}\t{}\t{}\tProtein\tProtein\n".format(
            from_id, to_id, from_fullname, to_fullname))

def parse_stringDB(source, trans):
    """
    This function parses the tab-delimited StringDB file. It converts the stringDB identifiers
    to uniprot identifiers and writes these to a temporary file for further processing
    :param source: Location of the StringDB file
    :param f: Output file of the temporary file
    :param trans: Translation dictionary with stringDB as keys and contains uniprot identifiers
    """
    print("Start parsing StringDB")
    f = open(options["folder"] + "string_to_uniprot.tsv", "w")
    f.write("from\tto\tconfidence\n")
    with open(source) as fp:
        for count, line in enumerate(fp):
            print("\rDone {}".format(count), end="\r")
            if count > 0:
                splitted = line.split(" ")
                from_uni = trans.get(splitted[0])
                to_uni = trans.get(splitted[1])

                if from_uni is not None and to_uni is not None:
                    if "Uniprot" in from_uni.keys() and "Uniprot" in \
                            to_uni.keys():
                        f.write("{}\t{}\t{}".format(from_uni["Uniprot"], to_uni[
                            "Uniprot"], splitted[-1]))

    f.close()
    print("\nDone parsing StringDB")


def find_neighbours(df, identifiers):
    """
    This function finds neighbours of the current identifiers (IDs).
    Afterwards, they are combined into a single list and the
    identifiers set is subtracted to get the unique neighbour identifiers.
    :param df: Dataframe with a size of n x 2 containing uniprot identifiers
    :param identifiers: A set of current identifiers
    :return: A set of neighbour identifiers (Uniprot)
    """
    part1 = df[df["from"].isin(identifiers)]["to"].tolist()
    return list(set(part1) - set(identifiers))


def disease(acc, name, elem):
    """
    This function writes to the 'Diseases' file if diseases are available
    :param acc: HMDB accession
    :param name: Metabolite name
    :param elem: Disease XML element
    """
    if len(list(elem)) > 0:
        disease = list(elem)[0].text
        OUTPUTS["Diseases"].write("{}\t{}\t{}\n".format(acc, name, disease))


def biospecimen(acc, name, elem, biospecimen):
    """
    This function writes to the 'Biospecimen' file if biospecimen are available
    :param acc: HMDB accession
    :param name: Metabolite name
    :param elem: Disease XML element
    :param biospecimen: A set of biospecimen as strings.
    """
    if elem.text not in biospecimen:
        biospecimen.add(elem.text)
        OUTPUTS["Biospecimen"].write(
            "{}\t{}\t{}\n".format(acc, name, elem.text))


def pathways(acc, name, elem):
    """
    This function writes to the 'Pathways' file if pathways are available
    :param acc: HMDB accession
    :param name: Metabolite name
    :param elem: Pathway XML element
    """
    for child in list(elem):
        if child.tag == PRESET + "name":
            pathway_name = child.text
        if child.tag == PRESET + "kegg_map_id" and child.text is not None:
            OUTPUTS["Pathways"].write(
                "{}\t{}\t{}\n".format(acc, name, pathway_name))


def protein(acc, name, elem, interacts):
    """
    This function writes to the 'Protein' file if protein are available
    :param acc: HMDB accession
    :param name: Metabolite name
    :param elem: Protein XML element
    """
    for child in list(elem):
        if child.tag == PRESET + "uniprot_id":
            uniprot = child.text
            to_name = id_converter.get(uniprot)  # Find translation for ID
            if to_name is not None:
                string = "{}\t{}\t{}\t{}\tMetabolite\tProtein\n".format(
                    acc, uniprot, name, to_name)

                if uniprot in interacts["Direct"]:
                    HMDB_0.add(acc)
                    INTERACTIONS["m_direct"].write(string)

                if uniprot in interacts["Indirect"]:
                    HMDB_1.add(acc)
                    INTERACTIONS["m_indirect_1"].write(string)

                if uniprot in interacts["Neighbour"]:
                    HMDB_2.add(acc)
                    INTERACTIONS["m_indirect_2"].write(string)

                OUTPUTS["Proteins"].write(string)


def parse_hmdb_xml(interacts):
    """
    This function parses the HMDB xml file. It finds metabolites associated with the protein-protein interactions
    found with the filter given. It also writes metadata about the proteins like diseases, subcellular locations and
    pathways. Metabolites in the xml file not associated with the protein-protein interactions can be still be added
    to the file by setting 'keep_metabolites' to True in the config file. This results in metabolites still be able
    to find in the Rshiny application but will not have any interactions.

    :param interacts: a dictionary containing all protein-protein interactions per
     neighbour set (direct, indirect, second indirect). This is used to find metabolites that are associated with
     these proteins.

    :return: dictionary containing ChEBI to HMDB identifiers.
    """
    print("Start HMDB XML parsing\n")
    source = options["folder"] + options["hmdb"]

    context = ET.iterparse(source, events=("end",))
    count = 0
    chebi_translation = {}
    event, root = next(context)
    output_alt = open(options["folder"] + options["alternative_hmdbs"], "w")
    output_alt.write("Accession\tSecondary\n")
    acc = False
    all_acc = {}
    met_name = False
    for event, elem in context:

        if elem.tag == PRESET + "accession":
            if not acc:
                count += 1
                print("\rParsing accession {}".format(count), end="\r")
                bio_specimen = set()
                accession = elem.text
                acc = True

            else:
                sec_accession = elem.text
                output_alt.write("{}\t{}\n".format(accession, sec_accession))

        elif elem.tag == PRESET + "name" and met_name == False:
            metabolite_name = elem.text.rstrip().lstrip()
            met_name = True
            all_acc[accession] = metabolite_name

        elif elem.tag == PRESET + "disease" and "Diseases" in OUTPUTS.keys():
            disease(accession, metabolite_name, elem)

        elif elem.tag == PRESET + "cellular" and "Cellular" in OUTPUTS.keys():
            OUTPUTS["Cellular"].write(
                "{}\t{}\t{}\n".format(accession, metabolite_name, elem.text))

        elif elem.tag == PRESET + "biospecimen" and "Biospecimen" in \
                OUTPUTS.keys():
            biospecimen(accession, metabolite_name, elem, bio_specimen)

        elif elem.tag == PRESET + "tissue" and "Tissues" in OUTPUTS.keys():
            OUTPUTS["Tissues"].write(
                "{}\t{}\t{}\n".format(accession, metabolite_name, elem.text))

        elif elem.tag == PRESET + "protein":
            protein(accession, metabolite_name, elem, interacts)

        elif elem.tag == PRESET + "pathway" and "Pathways" in OUTPUTS.keys():
            pathways(accession, metabolite_name, elem)

        elif elem.tag == PRESET + "chebi_id":
            chebi_translation[elem.text] = [accession, metabolite_name]

        elif elem.tag == PRESET + "metabolite":
            acc = False
            met_name = False
            elem.clear()

        if count % 10000 == 0:
            elem.clear()
            root.clear()

    for file in list(OUTPUTS.values()) + list(INTERACTIONS.values()):
        file.close()

    output_alt.close()

    if bool(options["keep_metabolites"]):
        keep_metabolites(all_acc)

    print("Done HMDB XML parsing")
    return chebi_translation


def keep_metabolites(all_acc):
    """
    This function adds metabolites to data files when they are not found in the filter. This way, they
    can be found in the Rshiny application, but do not have any connections to other metabolites or proteins.
    :param all_acc: A dictionary with HMDB accessions as key and its name as a value.
    """
    all_accs = set(all_acc.keys())
    for name in ["m_direct", "m_indirect_1", "m_indirect_2"]:
        file = open("{}{}".format(options["folder"], options[name]), "r")
        contents = set()
        for line in file.readlines():
            contents.add(line.split("\t")[0])
        file.close()
        file = open("{}{}".format(options["folder"], options[name]), "a")
        for rest_id in all_accs - contents:
            file.write("{}\t{}\t{}\t{}\tMetabolite\tMetabolite\n".format(
                rest_id, rest_id, all_acc[rest_id], all_acc[rest_id]))
        file.close()


def write_gos(comb, from_go, to_go, file, go_translation):
    """
    this function will write Gene Ontologies associated with proteins.
    :param comb: Data from both the 'from' as well as the 'to' protein (identifiers).
    :param from_go: list with associated GO's of the from-protein
    :param to_go: list with associated GO's of the to-protein
    :param file: file object of the GO file to write to
    :param go_translation: Dictionary of identifier : name pair.
    """
    if from_go is not None:
        for id in from_go:
            file.write("{}\t{}\t{}\n".format(comb[0], id, go_translation.get(id)))

    if to_go is not None:
        for id in to_go:
            file.write("{}\t{}\t{}\n".format(comb[1], id, go_translation.get(
                id)))


def parse_stringDB_with_filter2(direct, neighbours, trans):
    """
    This function filters StringDB again, but using the neighbours found.
    :param direct: Uniprot identifiers within the given GO filter.
    :param neighbours: Iniprot identifiers of the first neighbour set.
    :param trans: A dictionary containing the translation between StringDB and Uniprot
    :return: 2 variables are returned:
     - total: a dictionary containing all identifiers in all
    neighbour sets with their full names as values.
     - interacts: a dictionary containing all protein-protein interactions per
     neighbour set (direct, indirect, second indirect)
    """

    protein_file = open(options["folder"] + options["p_direct"], "w")
    protein_file1 = open(options["folder"] + options["p_indirect_1"], "w")
    protein_file2 = open(options["folder"] + options["p_indirect_2"], "w")

    go_file = open(options["folder"] + options["g_direct"], "w")
    go_file1 = open(options["folder"] + options["g_indirect_1"], "w")
    go_file2 = open(options["folder"] + options["g_indirect_2"], "w")

    for f in [protein_file, protein_file1, protein_file2]:
        f.write("from\tto\talias_a\talias_b\ttype_a\ttype_b\n")

    for f in [go_file, go_file1, go_file2]:
        f.write("Uniprot\tGO-id\tProcess\n")

    total = {}
    count = 0
    interacts = {"Direct" : [], "Indirect": [], "Neighbour": []}

    df = pd.read_csv(options["folder"] + options["stringdb"], sep=" ")

    for f, t, conf in zip(df["protein1"], df["protein2"], df["combined_score"]):
        count += 1
        print("\rDone {}".format(count), end="\r")
        from_uni = trans.get(f)
        to_uni = trans.get(t)
        conf = int(conf)

        if from_uni is not None and to_uni is not None and conf >= options["confidence"]:
            comb = [from_uni["Uniprot"], to_uni["Uniprot"]]
            from_go = from_uni.get("GO-process")
            to_go = to_uni.get("GO-process")
            if set(comb).issubset(direct):
                write_interaction(from_uni, to_uni, protein_file)
                interacts["Direct"] += comb
                write_gos(comb, from_go, to_go, go_file, go_translation)

            if any(x in direct for x in comb):
                write_interaction(from_uni, to_uni, protein_file1)
                write_gos(comb, from_go, to_go, go_file1, go_translation)
                interacts["Indirect"] += comb

            if any(x in neighbours for x in comb):
                write_interaction(from_uni, to_uni, protein_file2)
                write_gos(comb, from_go, to_go, go_file2, go_translation)
                total[from_uni["Uniprot"]] = from_uni["fullName"]
                total[to_uni["Uniprot"]] = to_uni["fullName"]
                interacts["Neighbour"] += comb

    for f in [go_file, go_file1, go_file2, protein_file, protein_file1, protein_file2]:
        f.close()

    if bool(options["keep_proteins"]):
        keep_proteins(total)

    interacts["Direct"] = set(interacts["Direct"])
    interacts["Indirect"] = set(interacts["Indirect"])
    interacts["Neighbour"] = set(interacts["Neighbour"])
    return total, interacts


def keep_proteins(total):
    """
    This function will add proteins to the direct, first indirect and second indirect if they
    were not found using the filter. It will add them having an interaction with themselves, so they will be plotted
    in Rshiny, but without any connection to other proteins.
    :param total: Dictionary with accessions and names.
    """

    all_accs = set(total.keys())
    for name in ["p_direct", "p_indirect_1", "p_indirect_2"]:
        file = open("{}{}".format(options["folder"], options[name]), "r")
        contents = set()
        for line in file.readlines():
            contents.add(line.split("\t")[0])
        file.close()
        file = open("{}{}".format(options["folder"], options[name]), "a")
        for rest_id in all_accs - contents:
            file.write("{}\t{}\t{}\t{}\tProtein\tProtein\n".format(
                rest_id, rest_id, total[rest_id], total[rest_id]))
        file.close()


def read_rxn(file):
    """
    The function reads a single RXN file from RheaDB. It will find all unique combinations of indentifiers
    :param file: The current file (without folder)
    :return: A tuple of all 'before' and 'after' identifiers from the chemical reaction
    """
    from_number = 0
    ids = []
    next_line = False
    with open(options["folder"] + "rxn/" + file, "r") as f:
        for count, line in enumerate(f):
            if next_line:
                from_to = line.lstrip().split(" ")
                from_number = int(from_to[0])
                next_line = False

            if line[0:4] == "RHEA":
                next_line = True

            if line[0:5] == "CHEBI":
                ids.append(line.split(":")[1].rstrip())
    return ids[:from_number], ids[from_number:]


def write_interactions(from_ids, to_ids, file, file1, file2):
    """
    This function writes metabolite-metabolite interactions to files based on
    the neighbour_number.
    :param output: File to which all metabolite-metabolite interactions are written
    :param from_ids: All ids in the 'before' section of a chemical equation
    :param to_ids: All ids in the 'after' section of a chemical equation
    :param filter: Boolean value if a filter value should be used, defaults to true
    """
    for id1 in from_ids:
        hmdb_id1 = chebi_translation.get(id1)
        for id2 in to_ids:
            hmdb_id2 = chebi_translation.get(id2)
            if hmdb_id1 is not None and hmdb_id2 is not None:
                comb = set([hmdb_id1[0], hmdb_id2[0]])
                if comb.issubset(HMDB_0):
                    file.write("{}\t{}\t{}\t{}\tMetabolite\tMetabolite\n".format(
                        hmdb_id1[0], hmdb_id2[0], hmdb_id1[1], hmdb_id2[1]))

                if comb.issubset(HMDB_1):
                    file1.write("{}\t{}\t{}\t{}\tMetabolite\tMetabolite\n".format(
                        hmdb_id1[0], hmdb_id2[0], hmdb_id1[1], hmdb_id2[1]))

                if comb.issubset(HMDB_2):
                    file2.write("{}\t{}\t{}\t{}\tMetabolite\tMetabolite\n".format(
                        hmdb_id1[0], hmdb_id2[0], hmdb_id1[1], hmdb_id2[1]))


def write_metadata():
    """
    This function will write headers of several metadata files. They will be used in other functions
    so the files aren't closed at this point.
    """
    outputs = interactions = {}

    for metabolite_info in options["Metabolite Info"]:
        outputs[metabolite_info] = open("{}{}{}.tsv".format(
            options["folder"], options["metabolite_prefix"], metabolite_info), "w")

    for name in ["m_direct", "m_indirect_1", "m_indirect_2"]:
        interactions[name] = open("{}{}".format(options["folder"], options[name]),
                                  "w")
        interactions[name].write("from\tto\talias_a\talias_b\ttype_a\ttype_b\n")

    for key, value in outputs.items():
        if key == "Proteins":
            value.write("{}\t{}\t{}\t{}\ttype_a\ttype_b\n".format(
                "Accession", "Uniprot", "Metabolite", "Protein"))
        else:
            value.write("{}\t{}\t{}\n".format("Accession", "Metabolite", key))
    return outputs, interactions


def write_metabolite_metabolite():
    """
    This function will write metabolite-metabolite interactions. It uses 3 files to write all the neighbours just like
    the other interactions files. Since RheaDB uses individual files, it reads each file and writes to
    the files respectively.
    """
    file = open(options["folder"] + options["metabolite_metabolite"], "w")
    file1 = open(options["folder"] + options["metabolite_metabolite1"], "w")
    file2 = open(options["folder"] + options["metabolite_metabolite2"], "w")
    for output in [file, file1, file2]:
        output.write("from\tto\talias_a\talias_b\ttype_a\ttype_b\n")

    print("writing intermediates")
    files = os.listdir(options["folder"] + "rxn/")
    for f in files:
        f_ids, t_ids = read_rxn(f)
        write_interactions(f_ids, t_ids, file, file1, file2)

    for output in [file, file1, file2]:
        output.close()



# Parse Uniprot
translation, direct, go_translation = parse_uniprot_xml(options["uniprot"])

# Parse and write StringDB to Uniprot identifiers
parse_stringDB(options["folder"] + options["stringdb"], translation)

# Find neighbours of given filter
df = pd.read_csv(options["folder"] + "string_to_uniprot.tsv", sep="\t")
n1 = find_neighbours(df, direct)
indirect_1 = set(n1 + direct)

# Filter StringDB using neighbours
id_converter, interacts = parse_stringDB_with_filter2(direct, indirect_1, translation)

# Parse HMDB
OUTPUTS, INTERACTIONS = write_metadata()
HMDB_0 = HMDB_1 = HMDB_2 = set()
PRESET = "{http://www.hmdb.ca}"

chebi_translation = parse_hmdb_xml(interacts)

# Parse and write RheaDB
write_metabolite_metabolite()

# Remove temporary file
os.remove(options["folder"] + "string_to_uniprot.tsv")

print("Done Preprocessing...")
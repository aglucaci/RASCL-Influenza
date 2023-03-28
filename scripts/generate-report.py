# =============================================================================
# Generates summary and annotation JSON files
# Results are compiled from statistically significant results
#    and mapped to genomic positions
# =============================================================================

# =============================================================================
# Imports
# =============================================================================
import argparse
import sys
import json
import re
import datetime
import os
import math
import csv
from os import path
from Bio import SeqIO
import BioExt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BioExt.uds import _align_par
from BioExt.scorematrices import (
    DNAScoreMatrix,
    FrequenciesError,
    ProteinScoreMatrix
)
import operator
from collections import defaultdict
from pathlib import Path
import glob

# =============================================================================
# Declares
# =============================================================================

# This is for H3N2
# Load the FASTAs from "../data/reference/H3N2/genome" (To do..)
#ref_genes = [["segment1", ], 
#             ["segment2", ""], [], [], [], [], [], []]


# This is for H3N2
# 1-based gene coordinates from NCBI
# from, to, name, within gene codon offset (to deal with ORF1a/b overlap)
gene_coordinates = [[30, 1730, "HA", 0, "segment4"],
                    [46, 1542, "NP", 0, "segment5"],
                    [119, 391, "PB1_F2", 0, "segment2"],
                    [25, 2298, "PB1", 0, "segment2"],
                    [20, 1429, "NA", 0, "segment6"],
                    [26, 784, "M1", 0, "segment7"], 
                    #[[26, 51] [740, 1007], "M2", 0, "segment7"],
                    [28, 2307, "PB2", 0, "segment1"],
                    [25, 2175, "PA", 0, "segment3"],
                    #[[25, 597], [599,784], "PA_X", 0, "segment3"],
                    #[[27, 56], [529, 864], "NEP", 0, "segment8"],
                    [27, 719, "NS1", 0, "segment3"]
]

score_matrix_ = BioExt.scorematrices.DNA95.load()

# =============================================================================
# Parse commandline arguments
# =============================================================================

arguments = argparse.ArgumentParser(
    description='Summarize selection analysis results.')
arguments.add_argument(
    '-f',
    '--file',
    help='File to process',
    required=True,
    type=str,
    nargs='*')
arguments.add_argument(
    '-p',
    '--pvalue',
    help='p-value to use',
    required=False,
    type=float,
    default=0.05)
arguments.add_argument(
    '-r',
    '--reference',
    help='the key sequence to highlight',
    type=str,
    default='MN908947')
arguments.add_argument(
    '-A',
    '--annotation',
    help='Write a JSON file with site annotations',
    required=True,
    type=str)
arguments.add_argument(
    '-S',
    '--summary',
    help='Write a JSON file here with segment annotations',
    required=True,
    type=str)
arguments.add_argument(
    '-d',
    '--default_tag',
    help='Default name for sequences that have no explicit label',
    required=False,
    type=str,
    default="Reference")

# =============================================================================
# Process commandline arguments
# =============================================================================
print()
print("# --------------------------------------------------------------")

import_settings = arguments.parse_args()
print("# Import settings file:", import_settings.file[0])

results_dir = os.path.dirname(os.path.realpath(import_settings.file[0]))
print("# Results directory:", results_dir)

base_dir = os.path.dirname(os.path.realpath(import_settings.file[0])).split("/")[0: -2]
base_dir = "/".join(base_dir)

data_dir = os.path.join(base_dir, "data")

print("# Data directory:", data_dir)

GENE_KEY = os.path.basename(import_settings.file[0]).split(".")[0]
FLU_TYPE_KEY = import_settings.file[0].split("/")[-2]

print("Gene key:", GENE_KEY)
print("Influenza type key:", FLU_TYPE_KEY)

# =============================================================================
# Load annotation and summary json's
# =============================================================================

annotation_json = None
summary_json = None

if (import_settings.annotation):
    try:
        print("# Opening annotation file:", import_settings.annotation)
        with open(import_settings.annotation) as ann:
            try:
                annotation_json = json.load(ann)
            except BaseException:
                annotation_json = {}
            # end try
        # end with
    except FileNotFoundError as fnf:
        annotation_json = {}
    # end try
# end if

if (import_settings.summary):
    try:
        print("# Opening summary file:", import_settings.summary)
        with open(import_settings.summary) as ann:
            try:
                summary_json = json.load(ann)
            except BaseException:
                summary_json = {}
            # end try
        # end with
    except FileNotFoundError as fnf:
        summary_json = {}
    # end try
# end if

# =============================================================================
# Helper functions
# =============================================================================

def load_reference_genome(ref_genome_dir):
    """
    genome = []
    ref_genome_dir = os.path.join(ref_genome_dir, flu_type, "genome")
    print("# Loading reference data from:", ref_genome_dir)
    print("# Looking for Gene Key:", gene_key)
    gene_segment_map = {"HA": "segment4"}
    ref_genome_dir_file = os.path.join(ref_genome_dir, gene_segment_map[gene_key] + ".fasta")
    record = SeqIO.read(ref_genome_dir_file, "fasta")
    genome.append("genome")
    genome.append(str(record.seq))
    return [genome]
    """

    genome = []

    # Get FASTA Files.
    for _file in glob.glob(os.path.join(ref_genome_dir, "*.fasta")):
        print("# Processing reference segment:", _file) 
        _basename = os.path.basename(_file)
        segment = _basename.split(".")[0]
        #print(segment)
        record = SeqIO.read(_file, "fasta")
        genome.append([segment, str(record.seq)])
 
    #end for
    return genome
#end method

def newick_parser(nwk_str, bootstrap_values, track_tags, json_map):
    global tags
    clade_stack = []
    automaton_state = 0
    current_node_name = ""
    current_node_attribute = ""
    current_node_annotation = ""
    quote_delimiter = None
    name_quotes = {
        "'": 1,
        '"': 1
    }

    def add_new_tree_level():
        new_level = {
            "name": None
        }
        the_parent = clade_stack[len(clade_stack) - 1]
        if ("children" not in the_parent):
            the_parent["children"] = []
        # end if
        clade_stack.append(new_level)
        the_parent["children"].append(clade_stack[len(clade_stack) - 1])
        clade_stack[len(clade_stack) -
                    1]["original_child_order"] = len(the_parent["children"])
    # end nested method

    def finish_node_definition():
        nonlocal current_node_name
        nonlocal current_node_annotation
        nonlocal current_node_attribute
        this_node = clade_stack.pop()
        if (bootstrap_values and "children" in this_node):
            this_node["bootstrap_values"] = current_node_name
        else:
            this_node["name"] = current_node_name
        # end if
        this_node["attribute"] = current_node_attribute
        this_node["annotation"] = current_node_annotation

        try:

            if 'children' not in this_node:
                node_tag = import_settings.default_tag
                if json_map:
                    tn = json_map["branch attributes"]["0"][this_node["name"]]
                else:
                    tn = this_node
                # end if
                nn = tn["original name"] if "original name" in tn else tn["name"]
                for k, v in tags.items():
                    if nn.find(k) >= 0:
                        node_tag = v
                        break
                    # end if
                # end for
            else:
                counts = {}
                node_tag = ""
                for n in this_node['children']:
                    counts[n["tag"]] = 1 + \
                        (counts[n["tag"]] if n["tag"] in counts else 0)
                # end for
                if len(counts) == 1:
                    node_tag = list(counts.keys())[0]
                # end if
            # end if
            this_node["tag"] = node_tag
        except Exception as e:
            print("Exception ", e)
        # end try

        if track_tags is not None:
            track_tags[this_node["name"]] = [
                this_node["tag"], 'children' in this_node]
        # end if

        current_node_name = ""
        current_node_attribute = ""
        current_node_annotation = ""
    # end nested method

    def generate_error(location):
        return {
            'json': None,
            'error':
            "Unexpected '" +
            nwk_str[location] +
            "' in '" +
            nwk_str[location - 20: location + 1] +
            "[ERROR HERE]" +
            nwk_str[location + 1: location + 20] +
            "'"
        }
    # end nested method
    tree_json = {
        "name": "root"
    }

    clade_stack.append(tree_json)

    space = re.compile("\\s")

    for char_index in range(len(nwk_str)):
        try:
            current_char = nwk_str[char_index]
            if automaton_state == 0:
                # look for the first opening parenthesis
                if (current_char == "("):
                    add_new_tree_level()
                    automaton_state = 1
                # end if
            elif automaton_state == 1 or automaton_state == 3:
                # case 1: // name
                # case 3: { // branch length
                # reading name
                if (current_char == ":"):
                    automaton_state = 3
                elif current_char == "," or current_char == ")":
                    try:
                        finish_node_definition()
                        automaton_state = 1
                        if (current_char == ","):
                            add_new_tree_level()
                        # end if
                    except Exception as e:
                        return generate_error(char_index)
                    # end try

                elif (current_char == "("):
                    if len(current_node_name) > 0:
                        return generate_error(char_index)
                    else:
                        add_new_tree_level()
                    # end if

                elif (current_char in name_quotes):
                    if automaton_state == 1 and len(current_node_name) == 0 and len(
                            current_node_attribute) == 0 and len(current_node_annotation) == 0:
                        automaton_state = 2
                        quote_delimiter = current_char
                        continue
                    # end if
                    return generate_error(char_index)
                else:
                    if (current_char == "["):
                        if len(current_node_annotation):
                            return generate_error(char_index)
                        else:
                            automaton_state = 4
                        # end if
                    else:
                        if (automaton_state == 3):
                            current_node_attribute += current_char
                        else:
                            if (space.search(current_char)):
                                continue
                            # end if
                            if (current_char == ";"):
                                char_index = len(nwk_str)
                                break
                            # end if
                            current_node_name += current_char
                        # end if
            elif automaton_state == 2:
                # inside a quoted expression
                if (current_char == quote_delimiter):
                    if (char_index < len(nwk_str - 1)):
                        if (nwk_str[char_index + 1] == quote_delimiter):
                            char_index += 1
                            current_node_name += quote_delimiter
                            continue
                        # end if
                    # end if

                    quote_delimiter = 0
                    automaton_state = 1
                    continue
                else:
                    current_node_name += current_char
                # end if
            elif automaton_state == 4:
                # inside a comment / attribute
                if (current_char == "]"):
                    automaton_state = 3
                else:
                    if (current_char == "["):
                        return generate_error(char_index)
                    # end if
                    current_node_annotation += current_char
                # end if
            # end if
        except Exception as e:
            return generate_error(char_index)
        # end try

    if (len(clade_stack) != 1):
        return generate_error(len(nwk_str) - 1)
    # end if

    if (len(current_node_name)):
        tree_json['name'] = current_node_name
    # end if

    return {
        'json': tree_json,
        'error': None
    }

# end method


def print_distribution(d, title, labels):
    print("#### %s" % title)
    print("| %s | %s |" % labels)
    print("|:---:|:---:|")
    for i in range(len(d)):
        r = d["%d" % i]
        print("| %.4g | %.4g |" % (r[labels[0]], r[labels[1]]))
    # end for
# end method


def extract_distribution(d, labels):
    result = []
    for i in range(len(d)):
        r = d["%d" % i]
        result.append([r[labels[0]], r[labels[1]]])
    # end for
    return result
# end method


def output_record(x):
    global aligned_str
    l = list(x)
    if len(l) == 1:
        aligned_str = l[0]
    # end if
# end nested method


def ignore_record(x):
    pass
# end nested method


def make_report_dict(row, indices):
    result = {}
    for i, t in indices:
        result[t] = row[i]
    # end for
    return result
# end method


def read_labels(json_file):
    tags = {}
    with open(json_file, "r") as cfh:
        try:
            tags = json.load(cfh)
        except BaseException:
            tags = {}
        # end try
    # end with
    return tags
# end method


def get_genomic_annotation(site):
    global gene_coordinates
    genomic_site_coord = -1
    gene = ""
    gene_site = -1
    if len(ref_seq_map):
        genomic_site_coord = ref_seq_map[site]
        if genomic_site_coord < 0:
            gene_site = "Not in SC2 (deletion)"
        else:
            gene = None
            for k in gene_coordinates:
                if k[0] <= genomic_site_coord and k[1] > genomic_site_coord:
                    genomic_site = (
                        (genomic_site_coord + k[3]) - k[0]) // 3
                    gene = k[2]
                    gene_site = genomic_site + 1
                    #gs = "%s %g" % (k[2], genomic_site_coord + 1)
                    break
            # end for
            if gene is None:
                gene = "Not mapped"
            # end if
        # end if
    else:
        gene = "N/A"
    # end if
    return (genomic_site_coord, gene, gene_site)
# end method


def process_cfel(json_file):
    global include_in_annotation, annotation_json, import_settings, site_reports, summary_json, summary_json_key

    print("# Processing CFEL JSON:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    with open(json_file, "r") as cfh:
        cfel = json.load(cfh)
        node_tags = {}
        the_tree = newick_parser(
            cfel["input"]["trees"]['0'],
            False,
            node_tags,
            cfel)['json']
    # end with

    if summary_json is not None:
        omegas = {}
        T = {}
        for k in [[k.split("*")[1], v[0][0]] for k, v in cfel['fits']
                  ['Global MG94xREV']['Rate Distributions'].items()]:

            if k[0] != 'background':  # Check
                test_map[k[0]] = "Test"
            else:
                test_map[k[0]] = "Reference"
            # end if
            omegas[k[0]] = k[1]
            T[k[0]] = 0.
        # end for

        for branch, nt in cfel["tested"]["0"].items():
            info = cfel["branch attributes"]["0"][branch]
            if nt != '':
                T[nt] += info["Global MG94xREV"]
            # end if

            if node_tags.get(branch, 0) != 0:
                node_tags[branch].append(info["Global MG94xREV"])
            else:
                node_tags.update({branch: info["Global MG94xREV"]})
            # end if
        # end for
        summary_json[summary_json_key]['rates'] = {
            'mean-omega': omegas, 'T': T}
    # end if
    beta_indices = []
    p_indices = []
    subs = []

    for i, tag in enumerate(cfel["MLE"]["headers"]):
        if tag[0].find('beta') == 0:
            beta_indices.append([i, re.split('\\(|\\)', tag[0])[1]])
        elif tag[0].find('P-value') == 0:
            p_indices.append([i, re.split('\\(|\\)|for ', tag[0])[1]])
        elif tag[0].find('subs') == 0:
            subs.append([i, re.split('\\(|\\)', tag[0])[1]])
        else:
            pass
        # end if
    # end for

    for i, row in enumerate(cfel["MLE"]["content"]["0"]):
        if annotation_json is not None and len(
                ref_map):  # if this is specified, write everything out
            gs = get_genomic_annotation(i)
            if gs[0] >= 0:
                include_in_annotation[i] = gs[0]
                annotation_json[gs[0]] = {
                    'G': gs[1],
                    'S': gs[2],
                    'index': i,
                    'bCFEL': {
                        'p': row[4],
                        'a': row[0],
                        'b': make_report_dict(row, beta_indices),
                        'p': make_report_dict(row, p_indices),
                        'pp': row[-2],
                        's': make_report_dict(row, subs),
                        'q': row[-3]
                    }
                }
            # end if
        # end if
        if row[-4] <= import_settings.pvalue:
            site_reports[i] = {'cfel': row}
        # end if
    # end for
    return cfel
# end method


def process_relax(json_file):
    global summary_json, summary_json_key

    print("# Processing RELAX JSON:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    with open(json_file, "r") as cfh:
        try:
            relax = json.load(cfh)
            if summary_json is not None:
                relax_d = {}
                #print (summary_json[summary_json_key]['rates']['mean-omega'])
                for r, rr in summary_json[summary_json_key]['rates']['mean-omega'].items():
                    relax_d[r] = []
                    for ignored, rd in relax["fits"]["RELAX alternative"]["Rate Distributions"][test_map[r]].items(
                    ):
                        relax_d[r].append(rd)

                summary_json[summary_json_key]['rates']['relax'] = relax_d
                summary_json[summary_json_key]['relax'] = {
                    'p': relax["test results"]["p-value"],
                    'K': relax["test results"]['relaxation or intensification parameter']}
            # end if
        except BaseException:
            print("Issue loading relax", file=sys.stderr)
        # end try
    # end with
# end method


def traverse_tree_in_order(
        node,
        labels,
        slac_data,
        i,
        parent_tag,
        root):
    node_tag = None

    nn = []

    if node is None:
        return

    nn = root if node["name"] == 'root' else node["name"]

    if nn in slac_data:
        node_tag = slac_data[nn]["codon"][0][i]
        #print (node_tag, parent_tag)
        if (parent_tag != node_tag):
            labels[nn] = node_tag
            labels[node["name"]] = node_tag
    else:
        print("Not in %s" % nn)

    if "children" in node:
        for c in node["children"]:
            traverse_tree_in_order(
                c, labels, slac_data, i, node_tag, root)
        # end for
    # end if
# end inner method


def match_node_names(qry_node, ref_node, mapping):
    #print (qry_node["name"])
    if "children" in qry_node and "children" in ref_node:
        mapping[ref_node["name"]] = qry_node["name"]
        if len(qry_node["children"]) != len(ref_node["children"]):
            # for c in qry_node["children"]:
            #    print (c["name"])
            #print ()
            # for c in ref_node["children"]:
            #    print (c["name"])

            raise Exception("Internal topology mismatch")
        for i, n in enumerate(ref_node["children"]):
            match_node_names(qry_node["children"][i], n, mapping)
    elif "children" in qry_node:
        raise Exception("Topology mismatch")
    elif "children" in ref_node:
        raise Exception("Topology mismatch")
    else:
        if qry_node["name"] != ref_node["name"]:
            raise Exception("Leaf name mismatch")
        # end if
    # end if
# end method


def process_busteds(json_file):
    global summary_json, summary_json_key

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    print("# Processing BUSTED[S] json:", json_file)

    with open(json_file, "r") as cfh:
        busted = json.load(cfh)
        if summary_json is not None:
            if "rates" in summary_json[summary_json_key]:

                summary_json[summary_json_key]['rates']['busted'] = busted["fits"]["Unconstrained model"]["Rate Distributions"]

                summary_json[summary_json_key]['busted'] = {
                    'p': busted["test results"]["p-value"],
                }
        # end if
    # end with
# end method


def def_value():
    return defaultdict(int)
# end method


def process_slac(json_file):
    global summary_json, summary_json_key, include_in_annotation, annotation_json

    print("# Processing SLAC json:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    with open(json_file, "r") as sh:
        slac = json.load(sh)
        compressed_subs = {}
        node_tags = {}
        the_tree = newick_parser(
            slac["input"]["trees"]['0'],
            False,
            node_tags,
            slac)['json']
        root_node = None

        if summary_json is not None:

            for branch, info in (slac["branch attributes"]["0"]).items():
                if branch in node_tags:
                    node_tags[branch].append(info["Global MG94xREV"])
                else:
                    root_node = branch
                # end if
            # end for

            summary_json[summary_json_key]['tree'] = slac["input"]["trees"]["0"]
            summary_json[summary_json_key]['tree_tags'] = node_tags
        # end if

        if len(include_in_annotation):
            for i in include_in_annotation:
                report = annotation_json[include_in_annotation[i]]
                counts_codon_site = {}
                counts_aa_site = {}

                gs = get_genomic_annotation(i)
                if gs[0] >= 0:
                    labels = {}
                    labels[root_node] = slac["branch attributes"]["0"][root_node]["codon"][0][i]
                    traverse_tree_in_order(
                        the_tree, labels, slac["branch attributes"]["0"], i, None, root_node)
                    compressed_subs[gs[0]] = labels
                # end if

                for k in set([k[0] for k in node_tags.values()]):
                    if len(k):
                        counts_codon_site[k] = defaultdict(int)
                        counts_aa_site[k] = defaultdict(int)
                    # end if
                # end for
                #print("node tags:", node_tags)
                if node_tags != {}:
                    for branch, tag in node_tags.items():
                        if len(tag[0]) > 0 and tag[1] == False:
                            codon = slac["branch attributes"]["0"][branch]["codon"][0][i]
                            aa = slac["branch attributes"]["0"][branch]["amino-acid"][0][i]
                            counts_codon_site[tag[0]][codon] += 1
                            counts_aa_site[tag[0]][aa] += 1
                        # end if
                    # end for
                # end if
                report['cdn'] = counts_codon_site
                report['aa'] = counts_aa_site
            # end for
            summary_json[summary_json_key]['subs'] = compressed_subs
            # end for
        # end if
    # end with
# end method


def process_prime(json_file):
    global summary_json, summary_json_key, include_in_annotation

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    print("# Processing PRIME JSON:", json_file)
    with open(json_file, "r") as ph:
        prime = json.load(ph)
        if summary_json is not None:
            h = prime["MLE"]["headers"]
            summary_json[summary_json_key]['prime-properties'] = [h[k]
                                                                  [1].replace('Importance for ', '') for k in range(6, len(h), 3)]
        # end if
        if len(include_in_annotation):
            for i in include_in_annotation:
                report = annotation_json[include_in_annotation[i]]

                prime_info = prime["MLE"]["content"]["0"][i]
                if prime_info:
                    report['prime'] = {
                        'p': [prime_info[k] for k in ([5, ] + list(range(7, len(prime_info), 3)))],
                        'lambda': [prime_info[k] for k in range(6, len(prime_info), 3)]
                    }
                else:
                    report['prime'] = None  # invariable
                # end if
            # end for
        # end if
# end method


def process_fade(json_file):
    global summary_json, summary_json_key, include_in_annotation

    print("# Processing FADE JSON:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    if os.path.exists(json_file):
        with open(json_file, "r") as ph:
            fade = json.load(ph)
            if len(include_in_annotation):
                for i in include_in_annotation:
                    report = annotation_json[include_in_annotation[i]]
                    report['fade'] = {}
                    for residue, info in fade["MLE"]["content"].items():
                        if len(residue) == 1:
                            report['fade'][residue] = {
                                'rate': info["0"][i][1],
                                'BF': info["0"][i][-1]
                            }
                        # end if
                    # end for
                # end for
            # end if
        # end with
    else:
        # Empty FADE file.
        pass
    # end if
# end method


def process_bgm(json_file):
    global summary_json, summary_json_key, include_in_annotation

    print("# Processing BGM JSON:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    with open(json_file, "r") as ph:
        bgm = json.load(ph)
        if summary_json is not None:
            try:
                summary_json[summary_json_key]['bgm'] = bgm["MLE"]["content"]
            except KeyError:
                summary_json[summary_json_key]['bgm'] = []
            # end try
        # end if
    # end with
# end method


def process_fel(json_file, cfel):
    global summary_json, summary_json_key, include_in_annotation, annotation_json, import_settings, site_reports

    print("# Processing FEL JSON:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    with open(json_file, "r") as ffh:
        fel = json.load(ffh)

        for i, row in enumerate(fel["MLE"]["content"]["0"]):
            if i in include_in_annotation:
                annotation_json[include_in_annotation[i]]['bFEL'] = {
                    'a': row[0],
                    'b': row[1],
                    'p': row[4]
                }

            if i in site_reports or row[4] <= import_settings.pvalue and row[1] > row[0]:
                if i in site_reports:
                    site_reports[i]["fel"] = row
                else:
                    site_reports[i] = {
                        'fel': row, 'cfel': cfel["MLE"]["content"]["0"][i]}
                # end if
            # end if
        # end for
    # end with
    return fel
# end method


def process_meme_internal(json_file, fel, cfel):
    global summary_json, summary_json_key, include_in_annotation, annotation_json, import_settings, site_reports

    print("# Processing MEME Internal JSON:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    with open(json_file, "r") as bh:
        meme = json.load(bh)

        for i, row in enumerate(meme["MLE"]["content"]["0"]):
            if i in include_in_annotation:
                annotation_json[include_in_annotation[i]]['bMEME'] = {
                    'p': row[6],
                    'a': row[0],
                    'b+': row[3],
                    'w+': row[4],
                    'b-': row[1],
                    'w-': row[2],
                    'br': row[7]
                }
            # end if

            if i in site_reports or row[6] <= import_settings.pvalue:
                if i in site_reports:
                    #site_reports[i]["meme"] = row
                    site_reports[i].update({"meme": row})
                # else:
                #    site_reports[i] = {'meme' : row,
                #                       'fel'  : fel ["MLE"]["content"]["0"][i],
                #                       'cfel' : cfel ["MLE"]["content"]["0"][i]}
                # end if
            # end if
        # end for
        for n, info in meme["branch attributes"]["0"].items():
            ## checking if the node is in the summary json ##
            if "tree_tags" not in summary_json[summary_json_key]:
                continue
            # end if

            if n in summary_json[summary_json_key]['tree_tags']:
                sig_sites = []
                ## check if the info at the node has the desired key ##
                ## iterating over each site within each node ##
                # print(info)
                #print(n, info["Posterior prob omega class by site"][1])
                if "Posterior prob omega class by site" in info.keys():
                    for pos, post_prob in enumerate(
                            info["Posterior prob omega class by site"][1]):
                        #print(f"CHECkING NODE {n} | SITE {pos+1}, values {j}")
                        ## compute EBF for each site ##
                        # Take the ratio of posterior/prior odds and you have
                        # yourself a Bayes Factor #
                        prior_prob = meme["MLE"]["content"]["0"][pos][4]
                        if post_prob != 0 and post_prob != 1 and prior_prob != 0 and prior_prob != 1:
                            posterior_odds = post_prob / (1 - post_prob)
                            prior_odds = prior_prob / (1 - prior_prob)
                            #print(f"POST -- {posterior_prob} | PRIOR -- {prior_prob}")
                            EBF_update = posterior_odds / prior_odds
                            # print(EBF_update)
                            if EBF_update >= 100:
                                #print(f"MADE IT ----- NODE {n} -- SITE {pos+1} -- EBF {EBF_update}")
                                sig_sites.append(include_in_annotation[pos])
                            # else:
                            #    print(f"DID NOT MAKE IT ----- NODE {n} -- SITE {pos+1} -- EBF {EBF_update}")
                            # end if

                        # end if
                    # end for
                # end if
                summary_json[summary_json_key]['tree_tags'][n].append(
                    sig_sites)
            else:
                # print("")
                pass
                #print ("Node %s is not in SLAC labeled nodeset" % n, file = sys.stderr)
            # end if
        # end for
# end method


def process_meme_full(json_file, fel, cfel):
    global summary_json, summary_json_key, include_in_annotation, annotation_json, import_settings, site_reports

    print("# Processing MEME Full JSON:", json_file)

    if not os.path.exists(json_file) or os.stat(json_file).st_size == 0:
        print("# File is empty or does not exist.")
        return
    # end if

    with open(json_file, "r") as bh:
        full_meme = json.load(bh)
        for i, row in enumerate(full_meme["MLE"]["content"]["0"]):
            if i in include_in_annotation:
                annotation_json[include_in_annotation[i]]['lMEME'] = {
                    'p': row[6],
                    'a': row[0],
                    'b+': row[3],
                    'w+': row[4],
                    'b-': row[1],
                    'w-': row[2],
                    'br': row[7]
                }
            # end if
            if i in site_reports or row[6] <= import_settings.pvalue:
                if i in site_reports:
                    #site_reports[i]["full-meme"] = row
                    site_reports[i].update({"full-meme": row})
                # else:
                #    site_reports[i] = {'full-meme' : row,
                #                       'meme' : meme ["MLE"]["content"]["0"][i],
                #                       'fel'  : fel ["MLE"]["content"]["0"][i],
                #                       'cfel' : cfel ["MLE"]["content"]["0"][i]}
                # ennd if
            # end if
        # end for
        # annotate branches with EBF support
        for n, info in full_meme["branch attributes"]["0"].items():
            if "tree_tags" not in summary_json[summary_json_key]:
                continue
            # end if
            if n in summary_json[summary_json_key]['tree_tags']:
                sig_sites = []
                if "Posterior prob omega class by site" in info.keys():
                    for pos, post_prob in enumerate(
                            info["Posterior prob omega class by site"][1]):
                        #print(f"CHECkING NODE {n} | SITE {pos+1}, values {j}")
                        ## compute EBF for each site ##
                        # Take the ratio of posterior/prior odds and you have
                        # yourself a Bayes Factor #
                        prior_prob = full_meme["MLE"]["content"]["0"][pos][4]
                        if post_prob != 0 and post_prob != 1 and prior_prob != 0 and prior_prob != 1:
                            posterior_odds = post_prob / (1 - post_prob)
                            prior_odds = prior_prob / (1 - prior_prob)
                            #print(f"POST -- {posterior_prob} | PRIOR -- {prior_prob}")
                            EBF_update = posterior_odds / prior_odds
                            # print(EBF_update)
                            if EBF_update >= 100:
                                #print(f"MADE IT ----- NODE {n} -- SITE {pos+1} -- EBF {EBF_update}")
                                sig_sites.append(include_in_annotation[pos])
                            # else:
                            #    print(f"DID NOT MAKE IT ----- NODE {n} -- SITE {pos+1} -- EBF {EBF_update}")
                            # end if
                        # end if
                    # end for
                # end if
                summary_json[summary_json_key]['tree_tags'][n].append(
                    sig_sites)

                # old code
                # for tag, ebf in info.items():
                #    bits = tag.split (" ")
                #    if len (bits) >=4 and ebf>=100:
                #        sig_sites.append (include_in_annotation[int (bits[2]) - 1])
                #summary_json[summary_json_key]['tree_tags'][n].append (sig_sites)
            else:
                #print ("Node %s is not in SLAC labeled nodeset" % n, file = sys.stderr)
                pass
            # end if
        # end for
    # end with
# end method

# =============================================================================
# Main subroutine
# =============================================================================

print()
print("# --------------------------------------------------------------")
print("# Working directory location:", base_dir)
print("# Loading reference genome segments")
flu_reference_genomes_dir = os.path.join(data_dir, "reference", FLU_TYPE_KEY, "genome") 
print("#", FLU_TYPE_KEY, "reference segmented genomes are located in:", flu_reference_genomes_dir)
#print("Gene key:", GENE_KEY)
#print("Influenza type key:", FLU_TYPE_KEY)
ref_genes = load_reference_genome(flu_reference_genomes_dir)
print("")
#print(ref_genes)

print("# --------------------------------------------------------------")
print("# Starting to process all results")
print()

for file_name in import_settings.file:
    print("# Input filename:", file_name)
    summary_json_key = None
    this_file = file_name.split("/")[-1].split(".")[0]
    # print("# This gene (file):", this_file)

    tags = {}
    label_json = os.path.join(results_dir, this_file + ".labels.json")
    print("# Opening:", label_json)
    tags = read_labels(label_json)
    summary_json_key = os.path.basename(this_file)
    if summary_json is not None and (summary_json_key not in summary_json):
        summary_json[summary_json_key] = {}
    # end if
    site_reports = {}
    ref_seq_re = re.compile(import_settings.reference)
    ref_seq_map = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        seq_id = seq_record.description
        if ref_seq_re.search(seq_id):
            ref_seq = str(seq_record.seq).upper()
            aligned_str = None
            for s in ref_genes:
                _align_par(SeqRecord(Seq(s[1]),
                                     id=s[0]),
                           [SeqRecord(Seq(ref_seq),
                                      id="ref")],
                           score_matrix_,
                           False,
                           False,
                           0.6,
                           ignore_record,
                           output_record)
                if (aligned_str is not None):
                    break
                # end if
            # end for
            ref_map = aligned_str.seq.strip('-')
            i = 0
            map_to_genome = []
            while i < len(ref_map):
                #ref_seq_map.append (c)
                if ref_map[i:i + 3] != '---':
                    map_to_genome.append(i)
                # end if
                i += 3
            # end while
            i = 0
            c = 0
            while i < len(ref_seq):
                if ref_seq[i:i + 3] != '---':
                    ref_seq_map.append(
                        map_to_genome[c // 3] + aligned_str.annotations['position'])
                    c += 3
                else:
                    ref_seq_map.append(-1)
                # end if
                i += 3
            # end while
        # end if
    # end for
    if summary_json is not None:
        summary_json[summary_json_key]['map'] = ref_seq_map
    # end if

    include_in_annotation = {}
    test_map = {}

    # Assign JSON Files ---
    print("# Assigning JSON results files")
    CFEL_JSON = os.path.join(results_dir, this_file + ".CFEL.json")
    RELAX_JSON = os.path.join(results_dir, this_file + ".RELAX.json")
    BUSTEDS_JSON = os.path.join(results_dir, this_file + ".BUSTEDS.json")
    SLAC_JSON = os.path.join(results_dir, this_file + ".SLAC.json")
    MEME_FULL_JSON = os.path.join(results_dir, this_file + ".MEME-full.json")
    MEME_JSON = os.path.join(results_dir, this_file + ".MEME.json")
    FEL_JSON = os.path.join(results_dir, this_file + ".FEL.json")
    PRIME_JSON = os.path.join(results_dir, this_file + ".PRIME.json")
    FADE_JSON = os.path.join(results_dir, this_file + ".FADE.json")
    BGM_JSON = os.path.join(results_dir, this_file + ".combined.fas.BGM.json")
    FMM_JSON = os.path.join(results_dir, this_file + ".FMM.json")

    # Process JSON Files ---
    cfel = process_cfel(CFEL_JSON)
    process_relax(RELAX_JSON)
    process_slac(SLAC_JSON)
    process_busteds(BUSTEDS_JSON)
    process_bgm(BGM_JSON)
    fel = process_fel(FEL_JSON, cfel)
    process_fade(FADE_JSON)
    process_prime(PRIME_JSON)
    process_meme_internal(MEME_JSON, fel, cfel)
    process_meme_full(MEME_FULL_JSON, fel, cfel)

# end

# =============================================================================
# Write to file
# =============================================================================

print()
print("# --------------------------------------------------------------")

print("# Writing to file ...", import_settings.annotation)

if annotation_json is not None:
    with open(import_settings.annotation, "w") as ann:
        json.dump(annotation_json, ann, indent=1)
    # end with
# end if

print("# Writing to file ...", import_settings.summary)

if summary_json is not None:
    with open(import_settings.summary, "w") as sm:
        json.dump(summary_json, sm, indent=1)
    # end with
# end if

# =============================================================================
# END OF FILE
# =============================================================================

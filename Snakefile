# Snakefile for Influenza analysis (RASCL-Flu).
# 2023
# @Author: Alexander Lucaci

#----------------------------------------------------------------------
# Imports
#----------------------------------------------------------------------
import os
import sys
import json
import csv
from pathlib import Path

#----------------------------------------------------------------------
# Declares
#----------------------------------------------------------------------
with open("config.json", "r") as in_sc:
  config = json.load(in_sc)
#end with

with open("cluster.json", "r") as in_c:
  cluster = json.load(in_c)
#end with

#----------------------------------------------------------------------
# User settings
#----------------------------------------------------------------------
BASEDIR = os.getcwd()

print("We are operating out of base directory:", BASEDIR)

# Settings
LABEL = config["label"] # Analysis label
QUERY_DATA_DIR = config["queryDataDir"]
BACKGROUND_DATA_DIR = config["backgroundDataDir"]
FILE_ENDING = config["fileEnding"]
REFERENCE_SEQUENCES = config["referenceSequencesDataDir"]
QUERY_FILE = config["queryFile"]
BACKGROUND_FILE = config["backgroundFile"]

print("# Processing query sequences:", QUERY_DATA_DIR, QUERY_FILE)
print("# Processing background sequences:", BACKGROUND_DATA_DIR, BACKGROUND_FILE)

#----------------------------------------------------------------------
# End -- User defined settings 
#----------------------------------------------------------------------

# For debugging or single gene analyses
#genes = ["HA", "NA"]

#genes = ["PB2","PB1_F2","PB1","PA_X","PA","NS1","NP","NEP","NA","M2","M1","HA"]
genes = ["PB2","PB1_F2","PB1","PA_X","PA","NS1","NP","NEP","NA","M2","M1","HA"]

# Set output directory
OUTDIR = os.path.join(BASEDIR, "results", LABEL)

# Create output dir.
Path(os.path.join(BASEDIR,"results")).mkdir(parents=True, exist_ok=True)
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

# Hyphy-analyses
HYPHY_ANALYSES_DIR = config["hyphy-analyses"]
FMM = os.path.join(HYPHY_ANALYSES_DIR, "FitMultiModel", "FitMultiModel.bf")
BUSTEDSMH = os.path.join(HYPHY_ANALYSES_DIR, "BUSTED-MH", "BUSTED-MH.bf")

#----------------------------------------------------------------------
# Rule All 
#----------------------------------------------------------------------
rule all:
    input:
        os.path.join(OUTDIR, QUERY_FILE + ".fa"),
        os.path.join(OUTDIR, BACKGROUND_FILE + ".fa"),
        #expand(os.path.join(OUTDIR, "{GENE}.query.fa"), GENE=genes),
        #expand(os.path.join(OUTDIR, "{GENE}.reference.fa"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.bam"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.msa.OG"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.msa.NS"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.msa.SA"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.compressed.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.query.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.background.bam"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.background.msa.OG"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.background.msa.NS"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.background.msa.SA"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.background.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.background.compressed.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.combined.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.AA.fas"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.combined.fas.raxml.bestTree"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.int.nwk"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.clade.nwk"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.full.nwk"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.combined.fas.BGM.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.SLAC.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME-full.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.PRIME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FADE.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTEDS.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.CFEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.ABSREL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.ABSRELS.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.ABSREL-MH.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.ABSRELS-MH.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTEDS-MH.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FMM.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX-MH.json"), GENE=genes),
        os.path.join(OUTDIR, LABEL + "_summary.json"),
        os.path.join(OUTDIR, LABEL + "_annotation.json")
#end rule -- all

#----------------------------------------------------------------------
# Rules -- Main analysis 
#----------------------------------------------------------------------
rule cleaner_query:
    input:
        input = os.path.join(BASEDIR, "data", QUERY_DATA_DIR, QUERY_FILE)
    output:
        output = os.path.join(OUTDIR, QUERY_FILE + ".fa")
    shell:
       "bash scripts/cleaner.sh {input.input} {output.output}"
#end rule

rule cleaner_background:
    input:
        input = os.path.join(BASEDIR, "data", BACKGROUND_DATA_DIR, BACKGROUND_FILE)
    output:
        output = os.path.join(OUTDIR, BACKGROUND_FILE + ".fa")
    shell:
       "bash scripts/cleaner.sh {input.input} {output.output}"
#end rule

#---------------------------------------------------------------------
# PROCESS QUERY SEQUENCES
#----------------------------------------------------------------------
rule bealign_query:
    input:
        in_genome = rules.cleaner_query.output.output,
        in_gene_RefSeq = os.path.join("data", "reference", REFERENCE_SEQUENCES, "{GENE}" + FILE_ENDING)
    output:
        output = os.path.join(OUTDIR, "{GENE}.query.bam")
    shell:
        "bealign -r {input.in_gene_RefSeq} -m HIV_BETWEEN_F {input.in_genome} {output.output}"
#end rule

rule bam2msa_query:
    input:
        in_bam = rules.bealign_query.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.query.msa.OG")
    shell:
        "bam2msa {input.in_bam} {output.out_msa}"        
#end rule

rule remove_stop_codons_query:
   input:
       input = rules.bam2msa_query.output.out_msa
   output:
       output = os.path.join(OUTDIR, "{GENE}.query.msa.NS")
   shell:
      "hyphy cln Universal {input.input} 'No/No' {output.output}"
#end rule

rule strike_ambigs_query:
   input:
       in_msa = rules.remove_stop_codons_query.output.output
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.query.msa.SA")
   conda: 'environment.yml'
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule

rule tn93_cluster_query:
    params:
        THRESHOLD_QUERY = config["threshold_query"],
        MAX_QUERY = config["max_query"] 
    input:
        in_msa = rules.strike_ambigs_query.output.out_strike_ambigs
    output:
        out_fasta = os.path.join(OUTDIR, "{GENE}.query.compressed.fas"),
        out_json = os.path.join(OUTDIR, "{GENE}.query.json")
    shell:
        "python3 scripts/tn93_cluster.py --input {input.in_msa} --output_fasta {output.out_fasta} --output_json {output.out_json} --threshold {params.THRESHOLD_QUERY} --max_retain {params.MAX_QUERY}"
#end rule

#----------------------------------------------------------------------
# Do the above for background sequences.
#----------------------------------------------------------------------
rule bealign_background:
    input:
        in_genome_background = rules.cleaner_background.output.output,
        in_gene_RefSeq = rules.bealign_query.input.in_gene_RefSeq
    output:
        output = os.path.join(OUTDIR, "{GENE}.background.bam")
    shell:
        "bealign -r {input.in_gene_RefSeq} -m HIV_BETWEEN_F -K {input.in_genome_background} {output.output}"
#end rule 

rule bam2msa_background:
    input:
        in_bam = rules.bealign_background.output.output
    output:
        out_msa = os.path.join(OUTDIR, "{GENE}.background.msa.OG")
    shell:
        "bam2msa {input.in_bam} {output.out_msa}"
#end rule

rule remove_stop_codons_background:
   input:
       input = rules.bam2msa_background.output.out_msa
   output:
       output = os.path.join(OUTDIR, "{GENE}.background.msa.NS")
   shell:
      "hyphy cln Universal {input.input} 'No/No' {output.output}"
#end rule

rule strike_ambigs_background:
   input:
       in_msa = rules.remove_stop_codons_background.output.output
   output:
       out_strike_ambigs = os.path.join(OUTDIR, "{GENE}.background.msa.SA")
   conda: 'environment.yml'
   shell:
      "hyphy scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule

rule tn93_cluster_background:
    params:
        THRESHOLD_background = config["threshold_background"],
        MAX_background = config["max_background"],
    input:
        in_msa = rules.strike_ambigs_background.output.out_strike_ambigs,
        in_gene_RefSeq = rules.bealign_query.input.in_gene_RefSeq
    output:
        out_fasta = os.path.join(OUTDIR, "{GENE}.background.compressed.fas"),
        out_json = os.path.join(OUTDIR, "{GENE}.background.json")
    shell:
        "python3 scripts/tn93_cluster.py --input {input.in_msa} --output_fasta {output.out_fasta} --output_json {output.out_json} --threshold {params.THRESHOLD_background} --max_retain {params.MAX_background} --reference_seq {input.in_gene_RefSeq}"
#end rule

# Combine them, the alignments ----------------------------------------------------
rule combine:
    params:
        THRESHOLD_QUERY = config["threshold_query"]
    input:
        in_compressed_fas = rules.tn93_cluster_query.output.out_fasta,
        in_msa = rules.tn93_cluster_background.output.out_fasta,
	in_gene_RefSeq = rules.bealign_query.input.in_gene_RefSeq
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas")
        #output_csv = os.path.join(OUTDIR, "{GENE}.combined.fas.csv")
    conda: 'environment.yml'
    shell:
        "python3 scripts/combine.py --input {input.in_compressed_fas} -o {output.output} --threshold {params.THRESHOLD_QUERY} --msa {input.in_msa} --reference_seq {input.in_gene_RefSeq}"
#end rule

# Convert to protein
rule convert_to_protein:
    input:
        combined_fas = rules.combine.output.output
    output:
        protein_fas = os.path.join(OUTDIR, "{GENE}.AA.fas")
    conda: 'environment.yml'
    shell:
        "hyphy conv Universal 'Keep Deletions' {input.combined_fas} {output.protein_fas}"
#end rule

# Combined ML Tree
rule raxml:
    params:
        THREADS = PPN
    input:
        combined_fas = rules.combine.output.output
    output:
        combined_tree = os.path.join(OUTDIR, "{GENE}.combined.fas.raxml.bestTree")
    shell:
        "raxml-ng --model GTR --msa {input.combined_fas} --threads {params.THREADS} --tree pars{{3}} --force"
#end rule

rule annotate:
    input:
       in_tree = rules.raxml.output.combined_tree,
       in_compressed_fas = rules.tn93_cluster_query.output.out_fasta
    output:
       out_int_tree = os.path.join(OUTDIR, "{GENE}.int.nwk"),
       out_clade_tree = os.path.join(OUTDIR, "{GENE}.clade.nwk"),
       out_full_tree = os.path.join(OUTDIR, "{GENE}.full.nwk")
    conda: 'environment.yml'
    shell:
       "bash scripts/annotate.sh {input.in_tree} 'REFERENCE' {input.in_compressed_fas} {LABEL} {BASEDIR}"
#end rule 

######################################################################
#---------------------Selection analyses ----------------------------#
######################################################################

rule slac:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.SLAC.json")
    conda: 'environment.yml'
    shell:
        "hyphy   SLAC --alignment {input.in_msa} --samples 0 --tree {input.in_tree} --output {output.output}"
#end rule -- slac

rule bgm:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.combined.fas.BGM.json")
    conda: 'environment.yml'
    shell:
        "hyphy BGM --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- bgm

rule fel:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FEL.json")
    conda: 'environment.yml'
    shell:
        "hyphy  FEL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- fel

rule meme:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME.json")
    conda: 'environment.yml'
    shell:
        "hyphy  MEME --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- MEME

# These are exlcuded from Minimal run (not implemented)
#rule absrel:
#    input:
#        in_msa = rules.combine.output.output,
#        in_tree = rules.annotate.output.out_int_tree
#    output:
#        output = os.path.join(OUTDIR, "{GENE}.ABSREL.json")
#    conda: 'environment.yml'
#    shell:
#        "hyphy ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- absrel

rule busted:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.BUSTEDS.json")
    conda: 'environment.yml'
    shell:
        "hyphy BUSTED --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL} --starting-points 10"
#end rule -- busted

rule relax:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.RELAX.json")
    conda: 'environment.yml'
    shell:
        "hyphy RELAX --alignment {input.in_msa} --models Minimal --tree {input.in_tree_clade} --output {output.output} --test {LABEL} --reference Reference --starting-points 10 --srv Yes"
#end rule -- relax
# End exclusion --

rule prime:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.PRIME.json")
    conda: 'environment.yml'
    shell:
        "hyphy  PRIME --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- prime

rule meme_full:
    input:
        in_msa = rules.combine.output.output,
        in_tree_full = rules.annotate.output.out_full_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.MEME-full.json")
    conda: 'environment.yml'
    shell:
        "hyphy  MEME --alignment {input.in_msa} --tree {input.in_tree_full} --output {output.output} --branches {LABEL}"
#end rule -- meme_full

rule fade:
    input:
        in_msa = rules.convert_to_protein.output.protein_fas,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FADE.json")
    conda: 'environment.yml'
    shell:
        "hyphy FADE --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL}"
#end rule -- fade

# cFEL
rule cfel:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.CFEL.json")
    conda: 'environment.yml'
    shell:
        "hyphy contrast-fel --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branch-set {LABEL} --branch-set Reference"
#end rule -- cfel

# MH Models ---
# aBSREL is left off, I am including here as a comparison for aBSREL-MH
rule absrel:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.ABSREL.json")
    conda: 'environment.yml'
    shell:
        "hyphy ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL}"
#end rule -- absrel

rule absrels:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.ABSRELS.json")
    conda: 'environment.yml'
    shell:
        "hyphy ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL} --srv Yes"
#end rule -- absrel

rule absrelmh:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.ABSREL-MH.json")
    conda: 'environment.yml'
    shell:
        "hyphy ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL} --multiple-hits Double+Triple"
#end rule 

rule absrelsmh:
    input:
        in_msa = rules.combine.output.output,
        in_tree = rules.annotate.output.out_int_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.ABSRELS-MH.json")
    conda: 'environment.yml'
    shell:
        "hyphy ABSREL --alignment {input.in_msa} --tree {input.in_tree} --output {output.output} --branches {LABEL} --multiple-hits Double+Triple --srv Yes"
#end rule 

rule bustedmh:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.BUSTEDS-MH.json")
    conda: 'environment.yml'
    shell:
        #"mpirun -np {PPN} HYPHYMPI {BUSTEDSMH} --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL} --starting-points 10"
        "hyphy {BUSTEDSMH} --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --branches {LABEL} --starting-points 10"
#end rule -- busted

rule fmm:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.FMM.json")
    conda: 'environment.yml'
    shell:
        "hyphy {FMM} --alignment {input.in_msa} --tree {input.in_tree_clade} --output {output.output} --triple-islands Yes"
#end rule -- busted

# RELAX-MH
rule relax_mh:
    input:
        in_msa = rules.combine.output.output,
        in_tree_clade = rules.annotate.output.out_clade_tree
    output:
        output = os.path.join(OUTDIR, "{GENE}.RELAX-MH.json")
    conda: 'environment.yml'
    shell:
        "hyphy RELAX --alignment {input.in_msa} --models Minimal --tree {input.in_tree_clade} --output {output.output} --test {LABEL} --reference Reference --starting-points 10 --srv Yes --multiple-hits Double+Triple"
#end rule -- relax

rule generate_report:
    input:
        expand(os.path.join(OUTDIR, "{GENE}.SLAC.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.CFEL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FADE.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.MEME-full.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.FADE.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.PRIME.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.ABSREL.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.BUSTEDS.json"), GENE=genes),
        expand(os.path.join(OUTDIR, "{GENE}.RELAX.json"), GENE=genes)
    output:
        SUMMARY_JSON = os.path.join(OUTDIR, LABEL + "_summary.json"),
        ANNOTATION_JSON = os.path.join(OUTDIR, LABEL + "_annotation.json")
    conda: 'environment.yml'
    shell:
         "bash scripts/process_json.sh {BASEDIR} {LABEL}"
#end rule generate_report


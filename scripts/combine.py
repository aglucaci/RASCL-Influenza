# Imports -------------------------------------------------------------
import os
import sys
import argparse
import json
import shutil
import csv
import Bio
from Bio import SeqIO

# Declares
# Argparse here
arguments = argparse.ArgumentParser(description='Cluster an MSA with genetic distance (TN93)')
arguments.add_argument('-i', '--input',            help = 'MSA file to process',                                  required = True, type = str )
arguments.add_argument('-o', '--output_fasta',           help = 'Output json file',                                     required = True, type = str)
arguments.add_argument('-m', '--msa',              help = 'Distance threshold for clustering query sequences',    required = True, type = str)
arguments.add_argument('--threshold',              help = 'Distance threshold for clustering query sequences',    required = True, type = float)
arguments.add_argument('-r', '--reference_seq',    help = 'Wuhan reference sequence',               required = False, type = str)
settings = arguments.parse_args()

# Output is {GENE}.combined.fas
# Declares
_ref_seq_name = ""
if settings.reference_seq:
    with open (settings.reference_seq) as fh:
        for l in fh:
            if l[0] == '>':
                _ref_seq_name = l[1:].split (' ')[0].strip()
                break
            #end if
        #end for
    #end with
#end if
           
print ("Reference seq_name %s" % _ref_seq_name)

query_compressed = settings.input #, {GENE}.compressed.fas
combined_msa = settings.output_fasta # which is the combined.fas
tn93_pairwise_calcs = combined_msa + ".csv"

threshold = settings.threshold # query threshold
ref_msa = settings.msa #reference.msa.SA, settings.msa

# Software
# Need to load this from a config.json file.
task_runners = {}
#task_runners['tn93'] = "/usr/local/bin/tn93"
task_runners['tn93'] = "tn93"

# Helper functions
def run_command (exec, arguments, filename, tag):
    #print (colored('Running ... %s\n' % (tag), 'cyan'))
    #print ("\t", colored('Command ... %s\n' % (" ".join ([exec] + arguments)), 'yellow'))
    cmd = " ".join ([exec] + arguments)
    print(cmd)
    result = os.system (cmd)
    if result != 0:
        #raise Exception ('Command exection failed code %s' % result)
        print ('Command exection failed code %s' % result)
        return None
    return os.path.getmtime(filename)
#end method

# Main -------
input_stamp = run_command (task_runners['tn93'], ['-o', tn93_pairwise_calcs, '-s', ref_msa, '-t', "%g" % (threshold*2.0), query_compressed], tn93_pairwise_calcs, "filtering reference sequuences that are closer than %g to any query cluster" % (threshold*2.0))

with open(tn93_pairwise_calcs) as fh:
    reader = csv.reader (fh, delimiter = ',')
    next (reader)
    seqs_to_filter = set ()
    for l in reader:
        seqs_to_filter.add (l[1])
    #end for
    if _ref_seq_name in seqs_to_filter:
        seqs_to_filter.remove (_ref_seq_name)
    #end if
#end with

# Copy query_compressed to output .combined.fas
shutil.copy (query_compressed, combined_msa)

ADD_REF = False
REF_SEQ = ""

# Open .combined.fas to add the close reference sequences to it.
with open (combined_msa, "a+") as fh:
    check_uniq = set ()
    for seq_record in SeqIO.parse(ref_msa, "fasta"):
        if not seq_record.name in seqs_to_filter:
            if seq_record.name == _ref_seq_name:
                print ("\n>%s\n%s" % ("REFERENCE", str(seq_record.seq)), file = fh)
                check_uniq.add ('REFERENCE')
                ADD_REF = True
            else:
                seq_id = seq_record.name
                while seq_id in check_uniq:
                    seq_id = seq_record.name + '_' + ''.join(random.choices ('0123456789abcdef', k = 10))
                #end while
                check_uniq.add (seq_id)
                print ("\n>%s\n%s" % (seq_id, str(seq_record.seq)), file = fh)
            #end if
        #end if
    #end for
#end with

#os.remove (pairwise)
#input_stamp = os.path.getmtime(combined_msa)
  
if ADD_REF == False:
    # Add the reference
    #with open (combined_msa, "a+") as fh:
    #    print ("\n>%s\n%s" % ("REFERENCE", str(REF_SEQ)), file = fh)
    #shutil.copy (settings.reference_seq, combined_msa)
    pass

sys.exit(0)
        
# End of file

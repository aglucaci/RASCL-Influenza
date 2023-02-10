#!/bin/bash
#@Usage: bash scripts/process_json.sh /home/aglucaci/RASCL
#@Usage: bash scripts/process_json.sh /data/shares/veg/aglucaci/RASCL BA.2-01312022


## Declares
BASEDIR=$1
TAG=$2

DATA_DIR="$BASEDIR"/results/"$TAG"

# Static settings
REF_TAG="REFERENCE"
ANNOTATION_JSON="$DATA_DIR"/"$TAG"_annotation.json
SUMMARY_JSON="$DATA_DIR"/"$TAG"_summary.json

for file in "$DATA_DIR"/*.combined.fas; do
   echo ""
   echo python3 scripts/generate-report.py -f $file -A $ANNOTATION_JSON -S $SUMMARY_JSON -r "$REF_TAG"
   python3 scripts/generate-report.py -f $file -A $ANNOTATION_JSON -S $SUMMARY_JSON -r "$REF_TAG"
done

exit 0

# end of file

#!/bin/bash

IN_TREE=$1
REF=$2
IN_ALN=$3
LABEL=$4

# This could be pulled from a config file.
BASEDIR=$5

DIR=`dirname $IN_TREE`
echo "# Base directory: "$DIR
NAME=${IN_TREE##*/}
echo "# Filename: "$NAME

#combined.fas.raxml.bestTree
PREFIX=${NAME%combined.fas.raxml.bestTree}
echo "# Gene prefix: "$PREFIX

PASSING="$DIR"/"$PREFIX"
echo "# Passing in full path prefix: "$PASSING
echo ""

echo "hyphy scripts/annotator.bf $IN_TREE $REF $IN_ALN $LABEL $PASSING"
hyphy scripts/annotator.bf $IN_TREE $REF $IN_ALN $LABEL $PASSING

exit 0

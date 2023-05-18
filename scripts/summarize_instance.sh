#!/bin/bash


if [ "$#" -ne 3 ]; then
    echo "usage: summarize_instance.sh <instance_dir> <root_dir> <data_dir>"
    exit
fi

INSTANCE=$1
ROOT=$2
DATA=$3

python3 $ROOT/python/summarize.py $INSTANCE $DATA 2>&1 | tee $INSTANCE/output_summarize.txt

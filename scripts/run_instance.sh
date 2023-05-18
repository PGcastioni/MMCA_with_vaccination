#!/bin/bash


if [ "$#" -ne 4 ]; then
    echo "usage: ruh_instance.sh <instance_dir> <root_dir> <data_dir> <model_script>"
    exit
fi

INSTANCE=$1
ROOT=$2
DATA=$3
MODEL_SCRIPT=$4

if [ -z "$SLURM_JOB_ID" ]
then
    # running locally
    echo "Running locally"
    export JULIA_NUM_THREADS=2
    julia --project=$ROOT/scripts $ROOT/scripts/$MODEL_SCRIPT $DATA $INSTANCE 2>&1 | tee $INSTANCE/output.txt
else
    # running on SLURM
    echo "Running on SLURM"
    export JULIA_NUM_THREADS=40
    srun -n1 -N1 --cpus-per-task=40 julia --project=$ROOT/scripts $ROOT/scripts/$MODEL_SCRIPT $DATA $INSTANCE 2>&1 | tee $INSTANCE/output.txt
fi

touch $INSTANCE/finished

#!/bin/bash



if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <data_folder> <instance_folder> <config_file> <date>" 
    exit
fi


data=$1
instance=$2
config=$3
date=$4
scripts_dir=$(dirname "$0")

echo $scripts_dir

julia --project=scripts/ $scripts_dir/run_simulations.jl -d $data -c $config -i $instance --start-date $date --end-date $date --export-compartments-time-t 1 --params $scripts_dir/../data/params.csv


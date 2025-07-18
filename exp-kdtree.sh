#!/bin/bash

resultados_file="$HOME/RESULTADOS/resultHauss.txt"

# Define the datasets and their types
datasets=(
    "$HOME/DATASETS/setsHD1/"
    "$HOME/DATASETS/setsHD2/"
    "$HOME/DATASETS/setsHD3/"
  #  "$HOME/DATASETS/setsHD4/"
  #  "$HOME/DATASETS/setsHD5/"
  #  "$HOME/DATASETS/setsHD6/"
)

# Define the distributions
distributions=("R" "G" "M")

# Rango de conjuntos
start=7
end=7

# Primero construimos
for dataset in "${datasets[@]}"; do
    for dist in "${distributions[@]}"; do
        for ((i = start; i <= end; i++)); do
            dataset1="${dataset}set1Conjunto${i}${dist}.csv"
            dataset2="${dataset}set2Conjunto${i}${dist}.csv"
            ./hausdorff "$dataset1" "$dataset2" 0 1
        done
    done
    # ejecutamos consultas
    for dist in "${distributions[@]}"; do
        for ((i = start; i <= end; i++)); do
            dataset1="${dataset}set1Conjunto${i}${dist}.csv"
            dataset2="${dataset}set2Conjunto${i}${dist}.csv"

           # ./hausdorff "$dataset1" "$dataset2" 0 0 >> "$resultados_file" # HDKD1
           # ./hausdorff "$dataset1" "$dataset2" 1 0 >> "$resultados_file" # HDKD2
            ./hausdorff "$dataset1" "$dataset2" 2 0 >> "$resultados_file" # KAMATA
            ./hausdorff "$dataset1" "$dataset2" 3 0 >> "$resultados_file" # HDKD2 -v2
        done
    done
done




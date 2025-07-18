#!/bin/bash

resultados_file="$HOME/RESULTADOS/resultHauss-reales.txt"

# Define the datasets and their types
datasets=(
    "$HOME/DATASETS/setsHDR/"
)

# Define the distributions
distributions=("T")

# Rango de conjuntos
start=1
end=5

# Primero construimos
for dataset in "${datasets[@]}"; do
    for dist in "${distributions[@]}"; do
        for ((i = start; i <= end; i++)); do
            dataset1="${dataset}set1Conjunto${i}${dist}.csv"
            dataset2="${dataset}set2Conjunto${i}${dist}.csv"
            ./hausdorff "$dataset1" "$dataset2" 0 1
        done
    done
done

    # ejecutamos consultas
start=1
end=4
for dataset in "${datasets[@]}"; do
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

./hausdorff $HOME/DATASETS/setsHDR/set1Conjunto5T.csv $HOME/DATASETS/setsHDR/set1Conjunto5T.csv 3 0 >> "$resultados_file" # HDKD2 -v2
./hausdorff $HOME/DATASETS/setsHDR/set1Conjunto5T.csv $HOME/DATASETS/setsHDR/set1Conjunto5T.csv 2 0 >> "$resultados_file" # KAMATA




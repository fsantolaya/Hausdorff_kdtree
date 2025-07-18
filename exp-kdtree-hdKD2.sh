#!/bin/bash

resultados_file="$HOME/RESULTADOS/resultHauss.txt"
ruta_dataset_uniforme="$HOME/DATASETS/setsHD3/UNIFORME/"
ruta_dataset_mixta="$HOME/DATASETS/setsHD3/MIXTA/"
ruta_dataset_gauss="$HOME/DATASETS/setsHD3/GAUSS/"


        
#hdKD2 IN UNIFORM DISTRIBUTION

./hausdorff $ruta_dataset_uniforme"set1Conjunto1R.csv" $ruta_dataset_uniforme"set2Conjunto1R.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_uniforme"set1Conjunto2R.csv" $ruta_dataset_uniforme"set2Conjunto2R.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_uniforme"set1Conjunto3R.csv" $ruta_dataset_uniforme"set2Conjunto3R.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_uniforme"set1Conjunto4R.csv" $ruta_dataset_uniforme"set2Conjunto4R.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_uniforme"set1Conjunto5R.csv" $ruta_dataset_uniforme"set2Conjunto5R.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_uniforme"set1Conjunto6R.csv" $ruta_dataset_uniforme"set2Conjunto6R.csv"  2 >> "$resultados_file"
#./hausdorff $ruta_dataset_uniforme"set1Conjunto7R.csv" $ruta_dataset_uniforme"set2Conjunto7R.csv"  2 >> "$resultados_file"

#hdKD2 IN GAUSS DISTRIBUTION

./hausdorff $ruta_dataset_gauss"set1Conjunto1G.csv" $ruta_dataset_gauss"set2Conjunto1G.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_gauss"set1Conjunto2G.csv" $ruta_dataset_gauss"set2Conjunto2G.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_gauss"set1Conjunto3G.csv" $ruta_dataset_gauss"set2Conjunto3G.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_gauss"set1Conjunto4G.csv" $ruta_dataset_gauss"set2Conjunto4G.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_gauss"set1Conjunto5G.csv" $ruta_dataset_gauss"set2Conjunto5G.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_gauss"set1Conjunto6G.csv" $ruta_dataset_gauss"set2Conjunto6G.csv"  2 >> "$resultados_file"
#./hausdorff $ruta_dataset_gauss"set1Conjunto7G.csv" $ruta_dataset_gauss"set2Conjunto7G.csv"  2 >> "$resultados_file"


#hdKD2 IN MIX DISTRIBUTION

./hausdorff $ruta_dataset_mixta"set1Conjunto4M.csv" $ruta_dataset_mixta"set2Conjunto4M.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_mixta"set1Conjunto5M.csv" $ruta_dataset_mixta"set2Conjunto5M.csv"  2 >> "$resultados_file"
./hausdorff $ruta_dataset_mixta"set1Conjunto6M.csv" $ruta_dataset_mixta"set2Conjunto6M.csv"  2 >> "$resultados_file"
#./hausdorff $ruta_dataset_mixta"set1Conjunto7M.csv" $ruta_dataset_mixta"set2Conjunto7M.csv"  2 >> "$resultados_file"
      
  

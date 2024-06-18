#!/bin/bash

# Asignar argumentos a variables
pileup_con1=$1
pileup_con2=$3
pileup_con3=$4
pileup_pat=$2  
chr=$5
start=$6
end=$7
nom_output=$8

# Convertir start y end a enteros, en caso necesario
start=$(($start))
end=$(($end))

# Imprimir los argumentos recibidos
echo "Argumentos recibidos:"
echo "pileup_con1: $pileup_con1"
echo "pileup_con2: $pileup_con2"
echo "pileup_con3: $pileup_con3"
echo "pileup_pat: $pileup_pat"
echo "chr: $chr"
echo "start: $start"
echo "end: $end"
echo "nom_output: $nom_output"  

# Crear un directorio para los archivos de salida
output_dir="${nom_output}_output"
mkdir -p "$output_dir"

# Función para procesar un archivo pileup
process_pileup() {
    local pileup_file=$1
    local output_prefix=$2

    # Ejecutar AWK para filtrar el archivo pileup 
    echo "Ejecutando AWK para $pileup_file"
    awk -v chr="$chr" -v start="$start" -v end="$end" '$1 == chr && $2 >= start && $2 <= end {print $0}' "$pileup_file" | tr ' ' '\t' > "${output_dir}/${output_prefix}_tab.bed"

    # Verificar si el archivo fue creado y tiene contenido
    if [ ! -s "${output_dir}/${output_prefix}_tab.bed" ]; then
        echo "Error: El archivo ${output_dir}/${output_prefix}_tab.bed no fue creado o está vacío."
        exit 1
    fi
}

# Procesar los archivos pileup
process_pileup "$pileup_con1" "${nom_output}_con1"
process_pileup "$pileup_con2" "${nom_output}_con2"
process_pileup "$pileup_con3" "${nom_output}_con3"
process_pileup "$pileup_pat" "${nom_output}_pat"  

# Imprimir los nombres de los archivos creados
echo "Archivos creados en el directorio $output_dir:"
ls -l "${output_dir}"


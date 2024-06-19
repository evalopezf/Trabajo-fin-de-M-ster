#!/bin/bash
#
#SBATCH -p [nombre_de_la_particion]       # Partición donde se ejecutará el trabajo
#SBATCH --gres=gpu:1                      # Solicita una GPU
#SBATCH --chdir=[directorio_de_trabajo]   # Directorio de trabajo
#SBATCH -J basecalling_twins              # Nombre del trabajo
#SBATCH --cpus-per-task=4                 # Número de CPUs asignadas (máximo 96)
#SBATCH -o [directorio_de_salida_logs]/%x_%j.out  # Archivo de salida para los logs

# Cargar entornos
source $HOME/.bashrc
mamba activate samtools


outdir="directorio_de_salida_para_archivos"
input_pod5="ruta_del_archivo_pod5"
reference="ruta_del_archivo_de_referencia_genomica"
dorado_bin="ruta_del_ejecutable_dorado"
model_r10="ruta_del_modelo_r10_para_basecalling"
mod_base_model="ruta_del_modelo_de_bases_modificadas"
output_bam="ruta_del_archivo_bam_de_salida"
error_log="ruta_del_archivo_de_log_de_errores"
sorted_bam="ruta_del_archivo_bam_ordenado_de_salida"

# Ejecutar Dorado basecaller
CUDA_VISIBLE_DEVICES=0 $dorado_bin basecaller \
    $model_r10 \
    $input_pod5 \
    --reference $reference \
    --modified-bases-models $mod_base_model \
    --secondary no \
    --min-qscore 10 > $output_bam 2> $error_log

# Validación y procesamiento posterior
if [ $? -eq 0 ]; then
    samtools sort $output_bam -o $sorted_bam
    samtools index $sorted_bam
else
    echo "dorado basecaller failed for sample 2023-09_GEMELOS_21-165. Check $error_log for details"
fi

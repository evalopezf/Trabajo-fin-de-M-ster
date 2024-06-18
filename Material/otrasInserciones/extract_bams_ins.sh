#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Uso: $0 chr start"
    exit 1
fi

# Capturar los argumentos
CHR=$1
START=$2
POSITION="${CHR}:${START}"

# Nombres de los archivos de salida
VCF_OUTPUT="only_ins_${CHR}_${START}.vcf"
READS_ID_OUTPUT="reads_ID_SVA_F1_${CHR}_${START}.txt"
READS_ID_LINES_OUTPUT="reads_ID_SVA_F1_lines_${CHR}_${START}.txt"
READSAM_OUTPUT="read_locs_${CHR}_${START}.sam"
SAM_OUTPUT="SVA_F1_ins_${CHR}_${START}.sam"
FASTQ_OUTPUT="SVA_F1_ins_${CHR}_${START}.fastq"
SORTED_BAM_OUTPUT="SVA_F1_sorted_${CHR}_${START}.bam"
INFO_OUTPUT="info_${CHR}_${START}.info"

# Filtrar el VCF para obtener solo las inserciones en la posición especificada
echo "Filtrando inserciones en la posición: $POSITION"
bcftools filter -i 'SVTYPE="INS"' "/home/alumno10/TFM/SVAF1_met/SVA_F1_pat.vcf.gz" -Ov -o $VCF_OUTPUT -r $POSITION

# Extraer los IDs de las inserciones capturadas
echo "Extrayendo IDs de inserciones..."
bcftools query -f '%INFO/RNAMES\n' $VCF_OUTPUT | tr ';' '\n' > $READS_ID_OUTPUT
tr ',' '\n' < $READS_ID_OUTPUT > $READS_ID_LINES_OUTPUT

# Filtrar el BAM usando los IDs extraídos
echo "Filtrando BAM para obtener las lecturas de inserciones..."
samtools view -h "/home/alumno10/TFM/basecalling_SUP/2023-06-02_Genoma_SVAF1_unalign_sorted.bam" | grep -F -f $READS_ID_LINES_OUTPUT > $READSAM_OUTPUT

# Convertir las lecturas SAM a formato FASTQ
echo "Convirtiendo SAM a FASTQ..."
samtools fastq -T MM,ML $READSAM_OUTPUT > $FASTQ_OUTPUT

# Alinear lecturas FASTQ usando Minimap2
echo "Alineando lecturas FASTQ con Minimap2..."
minimap2 -ax map-ont -y "/home/alumno10/TFM/SVAF1_met/consenso/SVAF1_cons_reformatted.fasta" $FASTQ_OUTPUT > $SAM_OUTPUT

# Convertir el archivo SAM a BAM y ordenarlo
echo "Convirtiendo SAM a BAM y ordenando..."
samtools view -bS $SAM_OUTPUT | samtools sort -o $SORTED_BAM_OUTPUT

# Indexar el archivo BAM ordenado
echo "Indexando BAM ordenado..."
samtools index $SORTED_BAM_OUTPUT

# Crear el archivo de información con el contenido de reads_ID_SVA_F1_lines.txt y only_ins.vcf
echo "Creando el archivo de información..."
{
    echo "Contenido de $READS_ID_LINES_OUTPUT:"
    cat $READS_ID_LINES_OUTPUT
    echo
    echo "Contenido de $VCF_OUTPUT:"
    cat $VCF_OUTPUT
} > $INFO_OUTPUT

# Eliminar archivos intermedios
echo "Eliminando archivos intermedios..."
rm $VCF_OUTPUT $READS_ID_OUTPUT $READS_ID_LINES_OUTPUT $READSAM_OUTPUT $SAM_OUTPUT $FASTQ_OUTPUT

echo "Proceso completado. Los resultados están en $SORTED_BAM_OUTPUT, $SORTED_BAM_OUTPUT.bai y $INFO_OUTPUT"


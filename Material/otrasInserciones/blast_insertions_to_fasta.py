import subprocess
import pandas as pd
from Bio import SeqIO

# Paso 1: Ejecutar BLAST y guardar los resultados en formato tabular
def run_blast(query_fasta, db_fasta, results_file):
    blast_cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_fasta,
        "-out", results_file,
        "-outfmt", "6"
    ]
    subprocess.run(blast_cmd, check=True)
    print("BLAST search completed.")

# Paso 2: Extraer los IDs únicos de los resultados de BLAST
def extract_unique_ids(results_file, ids_file):
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    blast_results = pd.read_csv(results_file, sep="\t", names=columns)
    unique_ids = blast_results["sseqid"].unique()
    with open(ids_file, "w") as f:
        for seq_id in unique_ids:
            f.write(f"{seq_id}\n")
    print("Unique IDs extracted.")

# Paso 3: Filtrar el archivo FASTA original
def filter_fasta(input_fasta, ids_file, output_fasta):
    with open(ids_file) as f:
        unique_ids = set(line.strip() for line in f)
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in unique_ids:
                SeqIO.write(record, output_handle, "fasta")
    print(f"Filtered sequences have been saved to {output_fasta}")

def main():
    query_fasta = "/home/alumno10/TFM/SVAF1_met/mast2_seq/mast2_svaf1.fa"  
    db_fasta = "/home/alumno10/TFM/VCF/seq_ins_pat.fasta"  
    results_file = "results_blast.txt"  
    ids_file = "unique_ids.txt"  
    output_fasta = "filtered_sequences.fasta" 
    

    run_blast(query_fasta, db_fasta, results_file)
    

    extract_unique_ids(results_file, ids_file)
    

    filter_fasta(db_fasta, ids_file, output_fasta)

if __name__ == "__main__":
    main()


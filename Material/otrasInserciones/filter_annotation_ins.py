
import requests
import pandas as pd
import os

# Función para obtener información de genes y elementos repetitivos
def get_info(chrom, start, end):
    base_url = f"https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}"
    headers = {"Content-Type": "application/json"}

    gene_response = requests.get(base_url, headers=headers, params={"feature": "gene"})
    if not gene_response.ok:
        print(f"Error fetching gene data for {chrom}:{start}-{end}: {gene_response.text}")
        return [], []

    gene_info = gene_response.json()
    
    repeat_response = requests.get(base_url, headers=headers, params={"feature": "repeat"})
    if not repeat_response.ok:
        print(f"Error fetching repeat data for {chrom}:{start}-{end}: {repeat_response.text}")
        return gene_info, []

    repeat_info = repeat_response.json()

    return gene_info, repeat_info

# Función para analizar una coordenada y calcular la distancia al gen más cercano
def analyze_coordinate(chrom, start, extend=100000):
    gene_info, repeat_info = get_info(chrom, start, start + 1)  
    
    repeat_status = "Si" if repeat_info else "No"

    if not gene_info:
        extended_gene_info, _ = get_info(chrom, start - extend, start + extend)
        if not extended_gene_info:
            gene_status = "Región intergénica"
            nearest_gene = "N/A"
            distance_to_nearest_gene = "N/A"
        else:
            distances = [min(abs(start - gene['start']), abs(start - gene['end'])) for gene in extended_gene_info]
            nearest_gene_index = distances.index(min(distances))
            nearest_gene = extended_gene_info[nearest_gene_index].get('external_name', 'N/A')
            distance_to_nearest_gene = distances[nearest_gene_index]
            gene_status = "Región intergénica"
    else:
        genes = [f"{gene.get('external_name', 'N/A')} (intrón {gene['start']}-{gene['end']})" for gene in gene_info]
        gene_status = ", ".join(genes)
        nearest_gene = gene_info[0].get('external_name', 'N/A')
        distance_to_nearest_gene = 0

    return chrom, start, repeat_status, gene_status, nearest_gene, distance_to_nearest_gene

# Ruta al archivo de resultados BLAST
results_file = r"C:\Users\evalo\OneDrive\Escritorio\Master_Bioinformatica\TFM\results2_blast.txt"

if not os.path.exists(results_file):
    raise FileNotFoundError(f"The file {results_file} does not exist.")
print(f"Reading BLAST results from {results_file}")


posibles = pd.read_csv(results_file, sep="\t", header=None)
print("BLAST results loaded successfully")


posibles.columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]



filtered_df = posibles[posibles["length"] > 2000].copy()


filtered_df.loc[:, "length"] = pd.to_numeric(filtered_df["length"], errors='coerce')
filtered_df.loc[:, "pident"] = pd.to_numeric(filtered_df["pident"], errors='coerce')
filtered_df.loc[:, "evalue"] = pd.to_numeric(filtered_df["evalue"], errors='coerce')
filtered_df.loc[:, "qstart"] = pd.to_numeric(filtered_df["qstart"], errors='coerce')
filtered_df.loc[:, "qend"] = pd.to_numeric(filtered_df["qend"], errors='coerce')

if 'qcovs' in filtered_df.columns:
    filtered_df.loc[:, "qcovs"] = pd.to_numeric(filtered_df["qcovs"], errors='coerce')

if 'qlen' in filtered_df.columns:
    filtered_df.loc[:, "qlen"] = pd.to_numeric(filtered_df["qlen"], errors='coerce')

numeric_columns = filtered_df.select_dtypes(include=['number'])
numeric_columns['sseqid'] = filtered_df['sseqid']

grouped = numeric_columns.groupby('sseqid').mean()

if 'qcovs' in grouped.columns:
    grouped = grouped[grouped["qcovs"] > 95]

grouped_reset = grouped.reset_index()
ids = list(grouped_reset['sseqid'].unique())

# Verificar las coordenadas obtenidas
print("Coordinates obtained from BLAST results:")
print(ids)

split_coords = []
for coord in ids:
    if ':' in coord:
        chrom, start = coord.split(':')
        split_coords.append((chrom, int(start)))

print(f"Total coordinates to analyze: {len(split_coords)}")

results = []
for chrom, start in split_coords:
    print(f"Analyzing {chrom}:{start}")
    result = analyze_coordinate(chrom, start)
    print(f"Result: {result}")
    results.append(result)

df = pd.DataFrame(results, columns=[
    "Chrom", "Start", "Sobre elemento repetitivo", "Región",
    "Gen más cercano", "Distancia al gen más cercano"
])
df['chr:start'] = df['Chrom'] + ':' + df['Start'].astype(str)
df_result = df[[
    'chr:start', "Gen más cercano", "Distancia al gen más cercano",
    "Región", "Sobre elemento repetitivo"
]]

output_file = r"C:\Users\evalo\OneDrive\Escritorio\Master_Bioinformatica\TFM\anotacion_inserciones.csv"

if os.path.exists(output_file):
    os.remove(output_file)

df_result.to_csv(output_file, sep=";", index=None, header=True)




import requests
import pandas as pd

# Función para obtener información de genes y elementos repetitivos
def get_info(chrom, start, end):
    base_url = f"https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}"
    headers = {"Content-Type": "application/json"}

    
    gene_response = requests.get(base_url, headers=headers, params={"feature": "gene"})
    gene_info = gene_response.json()
    
   
    repeat_response = requests.get(base_url, headers=headers, params={"feature": "repeat"})
    repeat_info = repeat_response.json()

    return gene_info, repeat_info

# Función para analizar una coordenada y calcular la distancia al gen más cercano
def analyze_coordinate(chrom, start, extend=100000):
    gene_info, repeat_info = get_info(chrom, start, start + 1)  
    
    # Verificar si está en un elemento repetitivo
    repeat_status = "Si" if repeat_info else "No"

    # Verificar si está en un gen o región intergénica y encontrar el gen más cercano
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
        nearest_gene = gene_info[0].get('external_name', 'N/A')  # Toma el primer gen si está dentro de un gen
        distance_to_nearest_gene = 0  # La coordenada está dentro de un gen

    return chrom, start, repeat_status, gene_status, nearest_gene, distance_to_nearest_gene



# Mostrar el DataFrame
posibles=pd.read_csv("results_blast.txt",sep="\t",header=None)
posibles.columns=[ "qseqid","sseqid" ,"pident", "length" ,"evalue" ,"qcovs" ,"qlen" ,"qstart", "qend"]
filtered_df = posibles[posibles["length"] > 2000]

# Convertir las columnas necesarias a numéricas para poder calcular la media
posibles["length"] = pd.to_numeric(posibles["length"], errors='coerce')
posibles["qlen"] = pd.to_numeric(posibles["qlen"], errors='coerce')
posibles["pident"] = pd.to_numeric(posibles["pident"], errors='coerce')
posibles["evalue"] = pd.to_numeric(posibles["evalue"], errors='coerce')
posibles["qcovs"] = pd.to_numeric(posibles["qcovs"], errors='coerce')
posibles["qstart"] = pd.to_numeric(posibles["qstart"], errors='coerce')
posibles["qend"] = pd.to_numeric(posibles["qend"], errors='coerce')
numeric_columns = filtered_df.select_dtypes(include=['number'])
numeric_columns['sseqid'] = filtered_df['sseqid']
grouped = numeric_columns.groupby('sseqid').mean()
grouped = grouped[grouped["qcovs"] > 95]
grouped_reset = grouped.reset_index()
ids = list(grouped_reset['sseqid'].unique())
split_coords = []
for coord in ids:
        chrom, start = coord.split(':')
        split_coords.append((chrom, int(start)))

# Analizar cada coordenada
results = [analyze_coordinate(chrom, start) for chrom, start in split_coords]

# Crear un DataFrame con los resultados
df = pd.DataFrame(results, columns=["Chrom", "Start", "Sobre elemento repetitivo", "Región", "Gen más cercano", "Distancia al gen más cercano"])
df['chr:start'] = df['Chrom'] + ':' + df['Start'].astype(str)
df_result = df[['chr:start', "Gen más cercano", "Distancia al gen más cercano", "Región", "Sobre elemento repetitivo"]]

df_result.to_csv("anotacion_inserciones.csv",sep=";",index=None,header=True)

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import requests

def download_file(url, filename):
    # Realiza una solicitud GET al URL
    response = requests.get(url)
    response.raise_for_status()  # Esto verificará si hubo un error

    # Guarda el contenido del archivo en un archivo local
    with open(filename, 'wb') as f:
        f.write(response.content)


def fill_missing_combinations(df, all_chromosomes, sva_types, value_column):
    # Crear un DataFrame con todas las combinaciones posibles de cromosomas y tipos de SVA
    all_combinations = pd.DataFrame([(chrom, sva) for chrom in all_chromosomes for sva in sva_types], columns=['chrom', 'sva_type'])
    
    # Realizar un merge con el DataFrame original para asegurar que todas las combinaciones estén presentes
    df_full = all_combinations.merge(df, on=['chrom', 'sva_type'], how='left')
    
    # Rellenar los valores faltantes en la columna de valor con 0
    df_full[value_column] = df_full[value_column].fillna(0)
    
    # Ordenar los cromosomas y tipos de SVA para mantener el orden correcto
    df_full['chrom'] = pd.Categorical(df_full['chrom'], categories=all_chromosomes, ordered=True)
    df_full['sva_type'] = pd.Categorical(df_full['sva_type'], categories=sva_types, ordered=True)
    df_full.sort_values(['chrom', 'sva_type'], inplace=True)
    
    return df_full

# Función para ordenar cromosomas
def ordenar_cromosomas(df, chrom_column):
    chrom_order = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    df[chrom_column] = pd.Categorical(df[chrom_column], categories=chrom_order, ordered=True)
    return df.sort_values(chrom_column)

# Función para crear gráficos
def visualización(df_pivot, title_bar, title_pie, ylabel):
    # Crear una figura con dos subgráficos
    fig = plt.figure(figsize=(15, 9))
    gs = GridSpec(1, 2, width_ratios=[1, 2])  # Relación de anchos para los subgráficos

    # Subgráfico para el gráfico circular (pie plot) a la izquierda
    ax2 = fig.add_subplot(gs[0])
    ax2.annotate('a)', xy=(0.005, 1.08), xycoords='axes fraction', ha='left', fontsize=13, fontweight='bold')

    colors = {
        'SVA_A': '#44a5c2',
        'SVA_B': '#ffae49',
        'SVA_C': '#ff4b49',
        'SVA_D': '#49ff4b',
        'SVA_E': '#494bff',
        'SVA_F': '#a549ff'
    }

    # Calcular el porcentaje de cada tipo de SVA en el genoma completo
    total_counts = df_pivot.sum(axis=0)
    total_counts_sorted = total_counts.reindex(['SVA_A', 'SVA_B', 'SVA_C', 'SVA_D', 'SVA_E', 'SVA_F'])

    # Crear el gráfico circular
    ax2.pie(total_counts_sorted, labels=total_counts_sorted.index, colors=[colors[key] for key in total_counts_sorted.index], 
            autopct='%1.1f%%', startangle=140, counterclock=False, wedgeprops={'edgecolor': 'black', 'linewidth': 1})

    ax2.set_title(title_pie)

    # Subgráfico para el gráfico de barras apiladas a la derecha
    ax1 = fig.add_subplot(gs[1])
    ax1.annotate('b)', xy=(0.005, 1.05), xycoords='axes fraction', ha='left', fontsize=12, fontweight='bold')

    # Graficar cada tipo de SVA
    bottom_values = [0] * len(df_pivot)
    for sva_type in df_pivot.columns:
        ax1.bar(df_pivot.index, df_pivot[sva_type], bottom=bottom_values, color=colors[sva_type], 
                edgecolor='black', linewidth=1, label=sva_type)
        bottom_values = [i+j for i,j in zip(bottom_values, df_pivot[sva_type])]

    ax1.legend(title='SVA Type')

    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel(ylabel)
    ax1.set_title(title_bar)
    ax1.tick_params(axis='x', rotation=90)

    # Mostrar la figura
    plt.tight_layout()
    plt.show()

# Lista de todos los cromosomas esperados y tipos de SVA
all_chromosomes = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
    'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
]
sva_types = ['SVA_A', 'SVA_B', 'SVA_C', 'SVA_D', 'SVA_E', 'SVA_F']

# Leer y procesar el primer conjunto de datos
df1 = pd.read_csv(r"C:\Users\evalo\OneDrive\Escritorio\Master_Bioinformatica\TFM\TFM_organizado\Distrib_SVAs\SVAs_patient.bed", sep="\t", header=None)
df1.columns = ['chrom', 'start', 'end', 'seqid', 'score', 'strand', 'sva_chrom', 'sva_start', 'sva_end', 'sva_type', 'sva_score', 'sva_strand']

# Filtrar datos
df1_unique = df1.drop_duplicates(subset=['sva_chrom', 'sva_start', 'sva_end'])
df1_unique = df1_unique[df1_unique["strand"] == "+"]
df1_unique = df1_unique.drop_duplicates(subset=['chrom', 'start', 'end'])
df1_barplot = df1_unique.groupby(["chrom", "sva_type"]).size().reset_index(name='count')

pattern = 'alt|Un|random'
df1_barplot = df1_barplot[~df1_barplot['chrom'].str.contains(pattern, regex=True)]

df1_barplot = fill_missing_combinations(df1_barplot, all_chromosomes, sva_types, 'count')

df1_pivot = df1_barplot.pivot(index='chrom', columns='sva_type', values='count').fillna(0)

# Crear gráficos para el primer conjunto de datos
visualización(df1_pivot, '', '', "count")

# Leer y filtrar el segundo conjunto de datos
df2 = pd.read_csv(r"C:\Users\evalo\OneDrive\Escritorio\Master_Bioinformatica\TFM\TFM_organizado\Distrib_SVAs\Ditrib_ins\RM_INS_pat.csv", sep=";")
df2['repeat'] = df2['repeat'].fillna('Unknown') 
filtered_data = df2[df2['repeat'].str.contains('SVA')].copy()
df2['class/family'] = df2['class/family'].fillna('Unknown')  
filtered_data = filtered_data[filtered_data['class/family'].str.contains('SVA')].copy()  

split_columns = filtered_data['query sequence'].str.split(':', expand=True)
filtered_data['chr'] = split_columns[0]
filtered_data['start'] = split_columns[1]
filtered_data = filtered_data.drop(columns=['query sequence'])
filtered_data = filtered_data[["chr", "start", "repeat"]]
filtered_data.columns = ["chrom", "start", "sva_type"]

# Agrupar y pivotar datos
df2_barplot = filtered_data.groupby(["chrom", "sva_type"], observed=True).size().reset_index(name='count')
df2_barplot = fill_missing_combinations(df2_barplot, all_chromosomes, sva_types, 'count')

df2_pivot = df2_barplot.pivot(index='chrom', columns='sva_type', values='count').fillna(0)

# Crear gráficos para el segundo conjunto de datos
visualización(df2_pivot, '', '', "count")

# Guardar los DataFrames procesados
df1_barplot.to_csv("distrb_ref.csv", sep="\t", index=None)
df2_barplot.to_csv("distrb_ins.csv", sep="\t", index=None)

 #Descargar y cargar tamaño de cromosomas
url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
filename = 'hg38.chrom.sizes'
download_file(url, filename)
chr_size = pd.read_csv("hg38.chrom.sizes", sep="\t", header=None)
pattern = 'alt|Un|random'
chr_size = chr_size[~chr_size[0].str.contains(pattern, regex=True)]
chr_size.columns = ["chrom", "size"]
chr_size["size(kb)"]=chr_size["size"]/1000


merged_df1 = chr_size.merge(df1_barplot, on="chrom")
merged_df1["Densidad"] = (merged_df1["count"] / merged_df1["size(kb)"] ) * 100


merged_df1=merged_df1[["chrom","sva_type","Densidad"]]
merged_df1 = merged_df1.dropna(subset=['sva_type'])
merged_df1=fill_missing_combinations(merged_df1, all_chromosomes, sva_types, 'Densidad')
df1_pivot_norm = merged_df1.pivot(index='chrom', columns='sva_type', values='Densidad').fillna(0)
visualización(df1_pivot_norm,"","","Densidad")


merged_df2 = chr_size.merge(df2_barplot, on="chrom")
merged_df2["Densidad"] = (merged_df2["count"] / merged_df2["size(kb)"] ) * 100
merged_df2=fill_missing_combinations(merged_df2, all_chromosomes, sva_types, 'Densidad')

merged_df2=merged_df2[["chrom","sva_type","Densidad"]]
merged_df2 = merged_df2.dropna(subset=['sva_type'])

df2_pivot_norm = merged_df2.pivot(index='chrom', columns='sva_type', values='Densidad').fillna(0)
visualización(df2_pivot_norm,"","","Densidad")
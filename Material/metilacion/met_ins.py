import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings

def leer_procesar_archivo(ruta, tipo_base):
    try:
        df = pd.read_csv(ruta, sep='\t', header=None)
        df[10] = df[10] / 100
        df = df[df[3] == tipo_base]
        df_selected = df[[0, 1, 2, 3, 5, 10]]
        df_selected.columns = ['chrom', 'start', 'end', 'tipo', 'strand', 'prob_met']
        return df_selected
    except Exception as e:
        warnings.warn(f"Advertencia: No se pudo leer el archivo {ruta}. Error: {e}", UserWarning)
        return pd.DataFrame()

def procesar_muestra(ruta, ins_name, tipo, muestra_nombre):
    data = []
    for tipo_met in ['m', 'h']:
        met = leer_procesar_archivo(ruta, tipo_met)
        prob_met = met['prob_met'].fillna(0).mean() if not met.empty else 0
        data.append([ins_name, tipo, muestra_nombre, tipo_met, prob_met])
    return data

def metilacion_region_dt(ins_name, tipo):
    data = []
    try:
        rutas = [
            f"/home/alumno10/TFM/met_inser_500/{ins_name}_output/{ins_name}_con_tab.bed",
            f"/home/alumno10/TFM/met_inser_500/{ins_name}_output/{ins_name}_con2_tab.bed",
            f"/home/alumno10/TFM/met_inser_500/{ins_name}_output/{ins_name}_con3_tab.bed",
            f"/home/alumno10/TFM/met_inser_500/{ins_name}_output/{ins_name}_pat_tab.bed"
        ]
        nombres_muestras = ['Control1', 'Control2', 'Control3', 'Paciente']

        for ruta, muestra_nombre in zip(rutas, nombres_muestras):
            data.extend(procesar_muestra(ruta, ins_name, tipo, muestra_nombre))

        df = pd.DataFrame(data, columns=['muestra', 'tipo_ins', 'nombre_muestra', 'tipomet', 'media'])
        return df
    except Exception as e:
        warnings.warn(f"Advertencia: Error procesando datos para {ins_name}: {e}", UserWarning)
        return pd.DataFrame(columns=['muestra', 'tipo_ins', 'nombre_muestra', 'tipomet', 'media'])

def plot_methylation_data(df, filename, figsize=(20, 15)):
    if df.empty:
        print("No hay datos para graficar.")
        return

    data_mean = df.groupby(['tipo_ins', 'nombre_muestra', 'tipomet'])['media'].mean().reset_index()
    data_mean['tipomet'] = pd.Categorical(data_mean['tipomet'], categories=['m', 'h'], ordered=True)

    # Función para añadir valores numéricos encima de las barras
    def add_values(ax):
        for p in ax.patches:
            ax.annotate(f'{p.get_height():.2f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', xytext=(0, 9), textcoords='offset points')

    grid = sns.FacetGrid(data_mean, row='nombre_muestra', col='tipomet', margin_titles=True, aspect=3, height=1.9, sharey=False)
    grid.map_dataframe(sns.barplot, x='tipo_ins', y='media', hue='tipo_ins', palette='muted', dodge=False)

    for ax in grid.axes[:, 1]:
        ax.set_ylim(0, 0.09)
    for ax in grid.axes[:, 0]:
        ax.set_ylim(0, 1.2)

    for ax in grid.axes.flatten():
        add_values(ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)

    grid.set_axis_labels("", "")
    grid.set_titles(col_template="{col_name}", row_template="{row_name}")
    grid.fig.subplots_adjust(hspace=0.3, wspace=0.3)

    grid.axes[0, 0].set_title('a) Metilación (mC)', fontsize=10)
    grid.axes[0, 1].set_title('b) Hidroximetilación (hmC)', fontsize=10)

    grid.fig.text(0.01, 0.5, 'Media', va='center', ha='center', rotation='vertical', fontsize=10)
    grid.fig.text(0.5, 0.01, 'Tipo de Inserción', va='center', ha='center', fontsize=10)

    grid.fig.set_size_inches(figsize)
    grid.savefig(filename)
    plt.close()

def obtener_resultados(consulta, sufix):
    results = []
    for i in consulta:
        ins_name = f"{i}_{sufix}"
        ins_id = f"{i}"
        result = metilacion_region_dt(ins_name, ins_id)
        results.append(result)

    df = pd.concat(results, ignore_index=True)
    for ins in consulta:
        for muestra in ['Control1', 'Control2', 'Control3', 'Paciente']:
            for tipo_met in ['m', 'h']:
                if not ((df['tipo_ins'] == ins) & (df['nombre_muestra'] == muestra) & (df['tipomet'] == tipo_met)).any():
                    df = df.append({'muestra': ins_name, 'tipo_ins': ins, 'nombre_muestra': muestra, 'tipomet': tipo_met, 'media': 0}, ignore_index=True)

    return df

def plot_stacked_bars_for_insertion(df, insertion_type, filename, figsize=(10, 6)):
    df_filtered = df[df['tipo_ins'] == insertion_type]

    categories = ['Control1', 'Control2', 'Control3', 'Paciente']
    types = ['m', 'h']

    data = {category: {t: df_filtered[(df_filtered['nombre_muestra'] == category) & (df_filtered['tipomet'] == t)]['media'].values[0]
                       if not df_filtered[(df_filtered['nombre_muestra'] == category) & (df_filtered['tipomet'] == t)].empty else 0
                       for t in types} for category in categories}

    plt.figure(figsize=figsize)
    bar_width = 0.4
    index = np.arange(len(categories))
    hydroxymethylation = [data[cat]['h'] for cat in categories]
    methylation = [data[cat]['m'] for cat in categories]

    bars_hydro = plt.bar(index, hydroxymethylation, bar_width, label='Hidroximetilación', color='green')
    bars_methyl = plt.bar(index, methylation, bar_width, bottom=hydroxymethylation, label='Metilación', color='lightgreen')

    def add_labels(bars):
        for bar in bars:
            height = bar.get_height()
            plt.annotate(f'{height:.4f}',
                         xy=(bar.get_x() + bar.get_width() / 2, bar.get_y() + height),
                         xytext=(0, 3),  # 3 points vertical offset
                         textcoords="offset points",
                         ha='center', va='bottom')

    add_labels(bars_hydro)
    add_labels(bars_methyl)

    plt.xlabel('Muestras')
    plt.ylabel('Media')
    plt.title(f'Metilación y Hidroximetilación para {insertion_type}')
    plt.xticks(index, categories)
    plt.legend()
    plt.ylim(0, 1.135)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# Obtener y guardar resultados para 'menos500'
consulta = ["chr15:45697481", "chr16:79830903", "chr1:173914366", "chr1:199957671", "chr5:142000088", "chr7:65580692", "chr8:117611591"]
df_results_menos = obtener_resultados(consulta, "menos500")
df_results_menos.to_csv("df_result_menos500.csv", sep='\t', index=None)
plot_methylation_data(df_results_menos, "menos500_plot.png", figsize=(15, 10))

# Obtener y guardar resultados para 'mas500'
df_results_mas = obtener_resultados(consulta, "mas500")
df_results_mas.to_csv("df_result_mas500.csv", sep='\t', index=None)
plot_methylation_data(df_results_mas, "mas500_plot.png", figsize=(15, 10))

plot_stacked_bars_for_insertion(df_results_mas, "chr1:173914366", "mas500_stacked_plot_chr1_173914366.png", figsize=(10, 6))
plot_stacked_bars_for_insertion(df_results_menos, "chr1:173914366", "menos500_stacked_plot_chr1_173914366.png", figsize=(10, 6))


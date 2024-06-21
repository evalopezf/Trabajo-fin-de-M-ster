# Scripts de Análisis de Metilación

Este conjunto de scripts permite obtener los archivos pileup filtrados de diferentes muestras en regiones específicas del genoma, facilitando la comparación de los niveles de metilación. A continuación se proporciona una descripción detallada del script principal y las instrucciones para su uso.

## meth_flank_ins.sh

El script `meth_flank_ins.sh` procesa archivos de secuenciación para generar archivos pileup filtrados en las regiones de interés. Este procesamiento es esencial para comparar los niveles de metilación entre distintas muestras.

## met_ins.sh
Genera las figuras y df para la visualización de resultados

### Uso

Para ejecutar `meth_flank_ins.sh`, es necesario proporcionar una serie de argumentos específicos. En este caso, se utiliza `run_methyl.sh` para lanzar `meth_flank_ins.sh` en paralelo con GNU Parallel.

Descripción de los Argumentos (example_args_list.txt)
- arg1,arg2,arg3,arg4: Ruta a los pileups de las muestras.
- arg5: chrom
- arg6: range-start
- arg7: range-end
- arg8: output name 

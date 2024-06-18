#!/bin/bash

# Nombre del archivo de argumentos
args_file="example_args_list.txt"

# Imprimir los comandos generados
parallel --colsep ' ' -j 4 ./meth_flank_ins.sh {1} {2} {3} {4} {5} {6} {7} {8} :::: $args_file


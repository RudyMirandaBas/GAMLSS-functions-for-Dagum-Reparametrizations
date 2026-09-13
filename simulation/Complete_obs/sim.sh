#!/bin/bash
## ============================================================
## Uso:  ./run_simulation_RDQ1.sh <tau>
## Ejemplo: ./run_simulation_RDQ1.sh 0.5
##
## Tambien se puede correr para varios valores de tau seguidos:
##   for t in 0.1 0.25 0.5 0.75 0.9; do ./run_simulation_RDQ1.sh $t; done
## ============================================================

if [ -z "$1" ]; then
  echo "Error: debes indicar el valor de tau."
  echo "Uso: $0 <tau>"
  echo "Ejemplo: $0 0.5"
  exit 1
fi

TAU="$1"

echo "Lanzando simulacion RDQ1 con tau = ${TAU}"
Rscript RDQ1_simulation.R "${TAU}"

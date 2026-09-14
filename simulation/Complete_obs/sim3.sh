#!/bin/bash
## ============================================================
## Corre los escenarios de Monte Carlo para RD
## para nu = 1.05, 1.1, 1.25, 1.5 (uno despues del otro)
##
## Uso:
##   ./run_all_RD_scenarios.sh
##
## Requiere que MC_RD_simulation.R este en el mismo directorio
## que este script.
## ============================================================

set -e  # si algo falla, se detiene en vez de seguir silenciosamente

NUS=(1.05 1.1 1.25 1.5)

echo "===== Corriendo RD ====="
for nu in "${NUS[@]}"; do
  echo "[RD] nu = ${nu} - iniciando"
  Rscript RD_simulation.R "${nu}"
  echo "[RD] nu = ${nu} - terminado"
done

echo "Todos los escenarios de RD terminaron."

#!/bin/bash
## ============================================================
## Corre los escenarios de Monte Carlo para RD con censura,
## para todas las combinaciones de nu x censura (uno tras otro)
##
## Uso:
##   ./run_all_RD_censored_scenarios.sh
##
## Requiere que MC_RD_censored_simulation.R (y find_delta.R) esten
## en el mismo directorio que este script.
## ============================================================

set -e  # si algo falla, se detiene en vez de seguir silenciosamente

NUS=(1.05 1.1 1.25 1.5)
CENS=(0.1 0.25 0.5)

echo "===== Corriendo RD con censura ====="
for cens in "${CENS[@]}"; do
  for nu in "${NUS[@]}"; do
    echo "[RD] censorship = ${cens}  nu = ${nu} - iniciando"
    Rscript MC_RD_censored_simulation.R "${nu}" "${cens}"
    echo "[RD] censorship = ${cens}  nu = ${nu} - terminado"
  done
done

echo "Todos los escenarios de RD con censura terminaron."

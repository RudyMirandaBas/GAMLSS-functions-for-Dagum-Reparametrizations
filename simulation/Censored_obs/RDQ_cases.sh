#!/bin/bash
## ============================================================
## Corre los escenarios de Monte Carlo con censura:
## primero todo RDQ1, y despues todo RDQ2.
## Cada uno para todas las combinaciones de tau x censura.
##
## Uso:
##   ./run_all_RDQ_censored_scenarios.sh
##
## Requiere que MC_RDQ1_censored_simulation.R,
## MC_RDQ2_censored_simulation.R y find_delta.R esten en el mismo
## directorio que este script.
## ============================================================

set -e  # si algo falla, se detiene en vez de seguir silenciosamente

#TAUS=(0.1 0.25 0.5 0.75 0.9)
TAUS=(0.1)
#CENS=(0.1 0.25 0.5)
CENS=(0.1)

echo "===== Corriendo RDQ1 con censura ====="
for cens in "${CENS[@]}"; do
  for tau in "${TAUS[@]}"; do
    echo "[RDQ1] censorship = ${cens}  tau = ${tau} - iniciando"
    Rscript RDQ1_withH2.R "${tau}" "${cens}"
    echo "[RDQ1] censorship = ${cens}  tau = ${tau} - terminado"
  done
done

#echo "===== Corriendo RDQ2 con censura ====="
#for cens in "${CENS[@]}"; do
#  for tau in "${TAUS[@]}"; do
#    echo "[RDQ2] censorship = ${cens}  tau = ${tau} - iniciando"
#    Rscript RDQ2_withH2.R "${tau}" "${cens}"
#    echo "[RDQ2] censorship = ${cens}  tau = ${tau} - terminado"
#  done
#done

echo "Todos los escenarios de RDQ1 y RDQ2 con censura terminaron."

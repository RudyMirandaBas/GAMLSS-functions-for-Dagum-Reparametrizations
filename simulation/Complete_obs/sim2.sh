#!/bin/bash

set -e  # si algo falla, se detiene en vez de seguir silenciosamente

TAUS=(0.1 0.25 0.5 0.75 0.9)

echo "===== Corriendo RDQ1 ====="
for tau in "${TAUS[@]}"; do
  echo "[RDQ1] tau = ${tau} - iniciando"
  Rscript RDQ1_simulation.R "${tau}"
  echo "[RDQ1] tau = ${tau} - terminado"
done

echo "===== Corriendo RDQ2 ====="
for tau in "${TAUS[@]}"; do
  echo "[RDQ2] tau = ${tau} - iniciando"
  Rscript RDQ2_simulation.R "${tau}"
  echo "[RDQ2] tau = ${tau} - terminado"
done

echo "Todos los escenarios (RDQ1 y RDQ2) terminaron."

#!/bin/bash

python3 attack_on_gs_test.py \
  --count=10 \
  --timeout=60 \
  --size=10 \
  --min_matrix_elem=-10000 \
  --max_matrix_elem=10000 \
  --min_poly_deg=5 \
  --max_poly_deg=15 \
  --min_poly_coef=-10000 \
  --max_poly_coef=10000

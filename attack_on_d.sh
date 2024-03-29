#!/bin/bash

python3 attack_on_d.py \
  --count=10 \
  --timeout=60 \
  --size=10 \
  --min_matrix_elem=0 \
  --max_matrix_elem=100000 \
  --min_poly_deg=5 \
  --max_poly_deg=15 \
  --min_poly_coef=-100000 \
  --max_poly_coef=100000 \
  --poly_deg_bound=20


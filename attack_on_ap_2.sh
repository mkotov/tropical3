#!/bin/bash

python3 attack_on_ap_2.py \
  --count=10 \
  --timeout=60 \
  --size=10 \
  --min_matrix_elem=-100000 \
  --max_matrix_elem=100000 \
  --min_matrix_param=-100000 \
  --max_matrix_param=100000 \
  --min_matrix_step=-100000 \
  --max_matrix_step=100000

#!/bin/bash

python3 attack_on_ap_1_test.py \
  --count=10 \
  --timeout=60 \
  --size=25 \
  --min_matrix_elem=-100000 \
  --max_matrix_elem=100000 \
  --min_matrix_param=-100000 \
  --max_matrix_param=100000

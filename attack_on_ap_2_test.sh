#!/bin/bash

python3 attack_on_ap_2_test.py \
  --count=10 \
  --timeout=60 \
  --size=25 \
  --min_matrix_elem=-10000 \
  --max_matrix_elem=10000

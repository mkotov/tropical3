#!/bin/bash

python3 attack_on_hld_test.py \
  --count=100 \
  --timeout=60 \
  --size=10 \
  --min_matrix_elem=-100000 \
  --max_matrix_elem=100000

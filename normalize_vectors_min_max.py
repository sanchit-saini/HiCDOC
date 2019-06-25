#!/usr/bin/env python3

import argparse
from sklearn.preprocessing import minmax_scale
import lib.parse_matrix as pm

parser = argparse.ArgumentParser(
  description = 'Normalize interaction vectors with min-max'
)
parser.add_argument('-i', required=True, help='Input matrix')
parser.add_argument('-o', required=True, help='Output matrix')
args = parser.parse_args()

vectors = pm.matrix_to_vectors(pm.import_sparse_matrix(args.i))

for chromosome in vectors['interactions']:
  for r, replicate in enumerate(vectors['interactions'][chromosome]):
    for c, column in enumerate(replicate):
      vectors['interactions'][chromosome][r][c] = minmax_scale(column)

pm.export_matrix(pm.vectors_to_sparse_matrix(vectors), args.o)

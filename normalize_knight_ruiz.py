#!/usr/bin/env python3

import argparse
import gcMapExplorer.lib as gmlib
import lib.parse_matrix as pm

# gcMapExplorer fix
# Waiting for pull https://github.com/rjdkmr/gcMapExplorer/pull/5
import numpy as np
def has_zeros(ccmap):
  ccmap.make_readable()
  zero_count = []
  for i in range(ccmap.matrix.shape[0]):
    if np.sum(ccmap.matrix[i]) != 0:
      zero_count.append((ccmap.matrix[i] == 0).sum())
  return np.array(zero_count).any()


parser = argparse.ArgumentParser(
  description = 'Normalize biological biases with Knight-Ruiz'
)
parser.add_argument('-i', required=True, help='Input matrix')
parser.add_argument('-o', required=True, help='Output matrix')
args = parser.parse_args()

matrix = pm.import_sparse_matrix(args.i)
ccmaps = pm.matrix_to_ccmaps(matrix)

ccmaps['interactions'] = {
  chromosome: [
    gmlib.normalizer.normalizeCCMapByKR(
      ccmap,
      percentile_threshold_no_data = 99 if has_zeros(ccmap) else None
    )
    for ccmap in replicates
  ] for chromosome, replicates in ccmaps['interactions'].items()
}

pm.export_matrix(pm.ccmaps_to_sparse_matrix(ccmaps), args.o)

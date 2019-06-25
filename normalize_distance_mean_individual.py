#!/usr/bin/env python3

import argparse
import numpy as np
import gcMapExplorer.lib as gmlib
import lib.parse_matrix as pm

parser = argparse.ArgumentParser(
  description = 'Reduce distance effect with mean contact frequency '
                'per distance for each individual replicate'
)
parser.add_argument('-i', required=True, help='Input matrix')
parser.add_argument('-o', required=True, help='Output matrix')
parser.add_argument('--expected', required=False, help='Output expected values')
args = parser.parse_args()

matrix = pm.import_sparse_matrix(args.i)
ccmaps = pm.matrix_to_ccmaps(matrix)

ccmaps['interactions'] = {
  chromosome: [
    gmlib.normalizer.normalizeCCMapByMCFS(
      ccmap,
      threshold_data_occup = 0,
      stats = 'mean'
    )
    for ccmap in replicates
  ] for chromosome, replicates in ccmaps['interactions'].items()
}

pm.export_matrix(pm.ccmaps_to_sparse_matrix(ccmaps), args.o)

if args.expected:

  vectors = pm.matrix_to_vectors(matrix)
  weak_bins = pm.find_weak_bins(vectors)

  expected = {}

  for chromosome, bins in vectors['bins'].items():

    expected[chromosome] = [
      [0 for i in range(bins)]
      for _ in vectors['replicates']
    ]

    for r, values in enumerate(vectors['interactions'][chromosome]):

      values = np.array(values, float)

      number_of_intersected_weak_columns = [
        len(weak_bins[chromosome][r] & set(range(0, bins - k)))
        for k in range(0, bins)
      ]
      number_of_intersected_weak_rows = [
        len(weak_bins[chromosome][r] & set(range(k, bins)))
        for k in range(0, bins)
      ]
      number_of_crosses = [
        [
          j - i
          for j in weak_bins[chromosome][r]
          for i in weak_bins[chromosome][r]
          if j >= i
        ].count(k)
        for k in range(0, bins)
      ]

      number_of_intersected_weak_bins = [
        number_of_intersected_weak_columns[k]
        + number_of_intersected_weak_rows[k]
        - number_of_crosses[k]
        for k in range(0, bins)
      ]

      for k in range(0, bins):
        diagonal = np.diagonal(values, k)
        length = len(diagonal) - number_of_intersected_weak_bins[k]
        expected[chromosome][r][k] = sum(diagonal)/length

  pm.export_diagonal(dict(
      entries = expected,
      bins = vectors['bins'],
      resolution = vectors['resolution'],
      replicates = vectors['replicates'],
      comments = vectors['comments']
  ), args.expected, name='distance')

#!/usr/bin/env python3

import argparse
import numpy as np
from sklearn.neighbors import RadiusNeighborsRegressor
import lib.parse_matrix as pm

parser = argparse.ArgumentParser(
  description = 'Reduce distance effect with a radius-neighbors regression '
                'on all combined replicates'
)
parser.add_argument('-i', required=True, help='Input matrix')
parser.add_argument('-o', required=True, help='Output matrix')
parser.add_argument('--expected', required=False, help='Output expected values')
args = parser.parse_args()

vectors = pm.matrix_to_vectors(pm.import_sparse_matrix(args.i))
ignored = pm.find_weak_bins(vectors)

expected = {
  chromosome: [] for chromosome in vectors['interactions']
}

for chromosome, bins in vectors['bins'].items():

  id_values = []

  for r, replicate in enumerate(vectors['interactions'][chromosome]):
    for v, vector in enumerate(replicate):
      if v in ignored[chromosome][r]:
        continue
      for c, cell in enumerate(vector):
        if c in ignored[chromosome][r]:
          continue
        if c > v:
          break
        id_values += [(v-c, cell)]

  xs, ys = np.transpose(id_values)

  rnr = RadiusNeighborsRegressor(
    radius = 10,
    weights = 'distance'
  ).fit(xs.reshape(-1, 1), ys)

  expected[chromosome] = rnr.predict(
   np.arange(0, bins).reshape(-1, 1)
  ).tolist()

  for r, replicate in enumerate(vectors['interactions'][chromosome]):
    values = np.array(replicate, float)

    for i, e in enumerate(expected[chromosome]):
      np.fill_diagonal(values[i:], np.diagonal(values[i:])/e)
      if i != 0:
        np.fill_diagonal(values[:,i:], np.diagonal(values[:,i:])/e)

    vectors['interactions'][chromosome][r] = values.tolist()

pm.export_matrix(pm.vectors_to_sparse_matrix(vectors), args.o)

if args.expected:
  pm.export_diagonal(dict(
    entries = {chromosome: [values] for chromosome, values in expected.items()},
    bins = vectors['bins'],
    resolution = vectors['resolution'],
    replicates = [','.join(vectors['replicates'])],
    comments = vectors['comments']
  ), args.expected, name='distance')

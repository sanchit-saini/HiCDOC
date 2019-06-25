#!/usr/bin/env python3

import argparse
import math
import numpy as np
from sklearn.metrics import silhouette_samples
import lib.parse_matrix as pm
from lib.constrained_k_means import cop_kmeans, l2_distance

parser = argparse.ArgumentParser(description = 'Detect compartments using '
                                               'constrained k-means')
parser.add_argument('-i', required = True, help = 'Input matrix')
parser.add_argument('-o', required = True, help = 'Output compartments')
parser.add_argument('-k', type = int, default = 2,
                    help = 'Number of compartments')
parser.add_argument('--distances', nargs = '+',
                    help = 'Output distances to centroids. '
                           'One file required per compartment')
parser.add_argument('--concordance', help = 'Output concordance')
parser.add_argument('--silhouette', help = 'Output Silhouette')
args = parser.parse_args()

vectors = pm.matrix_to_vectors(pm.import_sparse_matrix(args.i))
vectors = pm.filter_vectors(vectors)

# replicates = ['1.1', '2.2', '2.3', '3.1', '1.2', '2.1']
# indices = [[0, 4], [1, 2, 5], [3]]
indices = [
  [
    index for index, r in enumerate(vectors['replicates'])
    if r.split('.')[0] == condition
  ] for condition in sorted(set(r.split('.')[0] for r in vectors['replicates']))
]

# joint_columns
# {
#   chromosome: [        # condition 1
#                 [...],   # column 0 of replicate 1
#                 [...],   # column 0 of replicate 2
#                 [...],   # column 1 of replicate 1
#                 [...],   # column 1 of replicate 2
#                 ...
#               ],
#               [        # condition 2
#                 [...],   # column 0 of replicate 1
#                 [...],   # column 0 of replicate 2
#                 [...],   # column 0 of replicate 3
#                 [...],   # column 1 of replicate 1
#                 [...],   # column 1 of replicate 2
#                 [...],   # column 1 of replicate 3
#                 ...
#               ],
#               ...
# }

joint_columns = {
  chromosome: [
    [
      vectors['interactions'][chromosome][r][v]
      for v in range(bins)
      for r in indices[condition]
    ] for condition in range(len(indices))
  ] for chromosome, bins in vectors['bins'].items()
}

# must_link
# {
#   chromosome: [        # condition 1
#                 (0,1),   # column 0 of replicates 1 and 2
#                 (2,3),   # column 1 of replicates 1 and 2
#                 ...
#               ],
#               [        # condition 2
#                 (0,1),   # column 0 of replicates 1 and 2
#                 (0,2),   # column 0 of replicates 1 and 3
#                 (3,4),   # column 1 of replicates 1 and 2
#                 (3,5),   # column 1 of replicates 1 and 3
#                 ...
#               ],
#               ...
# }

must_link = {
  chromosome: [
    [
      pair for i in range(
        0, bins*len(indices[condition]), len(indices[condition])
      )
      for pair in [(i, i+j) for j in range(1, len(indices[condition]))]
    ] for condition in range(len(indices))
  ] for chromosome, bins in vectors['bins'].items()
}

compartments = {
  chromosome: [
    [] for _ in vectors['replicates']
  ] for chromosome in vectors['bins']
}

if args.distances:
  distances = [{
    chromosome: [
      [] for _ in vectors['replicates']
    ] for chromosome in vectors['bins']
  } for i in range(args.k)]

if args.concordance:
  concordance = {
    chromosome: [
      [] for _ in vectors['replicates']
    ] for chromosome in vectors['bins']
  }

if args.silhouette:
  silhouette = {
    chromosome: [
      [] for _ in vectors['replicates']
    ] for chromosome in vectors['bins']
  }

kmeans = {
  chromosome: [
    {
      'clusters': [],
      'centroids': []
    } for condition in indices
  ] for chromosome in vectors['interactions']
}

for chromosome in joint_columns:

  # Detect compartments
  for condition in range(len(indices)):

    clusters, centroids = cop_kmeans(
      dataset = joint_columns[chromosome][condition],
      k = args.k,
      ml = must_link[chromosome][condition]
    )

    kmeans[chromosome][condition]['clusters'] = clusters
    kmeans[chromosome][condition]['centroids'] = centroids

  # Make compartments in each condition correspond
  # The centroids that are closest to each other are assumed
  # to be of the same compartment
  for condition in range(1, len(indices)):

    choices = sorted([
      (l2_distance(p, q), a, b)
      for a, p in enumerate(kmeans[chromosome][0]['centroids'])
      for b, q in enumerate(kmeans[chromosome][condition]['centroids'])
    ])

    correspondence = [-1] * args.k

    while choices:
      a = choices[0][1]
      b = choices[0][2]
      correspondence[b] = a
      i = 0
      while i < len(choices):
        if choices[i][1] == a or choices[i][2] == b:
          choices.pop(i)
        else:
          i += 1

    kmeans[chromosome][condition]['clusters'] = [
      correspondence[i] for i in kmeans[chromosome][condition]['clusters']
    ]

    kmeans[chromosome][condition]['centroids'] = [
      kmeans[chromosome][condition]['centroids'][correspondence[i]]
      for i in range(args.k)
    ]

  for condition in range(len(indices)):

    clusters = np.array(kmeans[chromosome][condition]['clusters']).reshape(
      -1, len(indices[condition])
    ).transpose().tolist()

    centroids = kmeans[chromosome][condition]['centroids']

    for i, index in enumerate(indices[condition]):
      compartments[chromosome][index] = clusters[i]

      if args.distances:
        for c, centroid in enumerate(centroids):
          distances[c][chromosome][index] = [
            l2_distance(vector, centroid)
            for vector in vectors['interactions'][chromosome][index]
          ]

      # Concordance is computed between the first 2 centroids
      if args.concordance:
        min_value, max_value = [
          math.log(
            (l2_distance(centroid, centroids[0]) + 1e-10)
            / (l2_distance(centroid, centroids[1]) + 1e-10)
          ) for centroid in centroids
        ]
        a, b = [-1, 1]

        concordance[chromosome][index] = [
          (b - a) * (
            (
              math.log(
                (l2_distance(vector, centroids[0]) + 1e-10)
                / (l2_distance(vector, centroids[1]) + 1e-10)
              ) - min_value
            ) / (max_value - min_value)
          ) + a
          for vector in vectors['interactions'][chromosome][index]
        ]

      # Add filtered regions to the results
      for removed in sorted(vectors['removed'][chromosome]):
        compartments[chromosome][index] = (
          compartments[chromosome][index][:removed]
          + [None]
          + compartments[chromosome][index][removed:]
        )
        if args.distances:
          for c, centroid in enumerate(centroids):
            distances[c][chromosome][index] = (
              distances[c][chromosome][index][:removed]
              + [None]
              + distances[c][chromosome][index][removed:]
            )
        if args.concordance:
          concordance[chromosome][index] = (
            concordance[chromosome][index][:removed]
            + [None]
            + concordance[chromosome][index][removed:]
          )

    if args.silhouette:
      coefficients = silhouette_samples(
        np.vstack([
          vectors['interactions'][chromosome][index]
          for index in indices[condition]
        ]),
        np.array(clusters).flatten()
      ).reshape(len(indices[condition]), -1).tolist()

      for i, index in enumerate(indices[condition]):
        silhouette[chromosome][index] = coefficients[i]

        # Add filtered regions to the results
        for removed in sorted(vectors['removed'][chromosome]):
          silhouette[chromosome][index] = (
            silhouette[chromosome][index][:removed]
            + [None]
            + silhouette[chromosome][index][removed:]
          )

pm.export_diagonal(dict(
    entries = compartments,
    bins = vectors['bins'],
    resolution = vectors['resolution'],
    replicates = vectors['replicates'],
    comments = vectors['comments']
), args.o)

if args.distances:
  for i, f in enumerate(args.distances):
    pm.export_diagonal(dict(
        entries = distances[i],
        bins = vectors['bins'],
        resolution = vectors['resolution'],
        replicates = vectors['replicates'],
        comments = vectors['comments']
    ), f, name = 'distance')

if args.concordance:
  pm.export_diagonal(dict(
      entries = concordance,
      bins = vectors['bins'],
      resolution = vectors['resolution'],
      replicates = vectors['replicates'],
      comments = vectors['comments']
  ), args.concordance)

if args.silhouette:
  pm.export_diagonal(dict(
      entries = silhouette,
      bins = vectors['bins'],
      resolution = vectors['resolution'],
      replicates = vectors['replicates'],
      comments = vectors['comments']
  ), args.silhouette)

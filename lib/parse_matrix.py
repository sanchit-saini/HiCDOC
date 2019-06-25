# This library is designed for intrachromosomal matrix manipulation

import itertools
import re
import numpy as np
from functools import reduce
import gcMapExplorer.lib as gmlib

# Create a sparse matrix dictionary from a file
#
# chromosome    position 1    position 2    replicate 1.1    ...
#
# {
#   interactions: {
#               (chromosome, position 1, position 2): [interaction 1, ...],
#               ...
#             },
#   sizes: {
#            'chromosome 1': 12780000,
#             ...
#          },
#   resolution: 10000,
#   replicates: ['1.1', '1.2', '2.1', '2.2'],
#   comments: ['# tissue: heart', '# normalization: cyclic loess']
# }
def import_sparse_matrix(file, header=True, diagonal=False):

  interactions = {}
  sizes = {}
  replicates = []
  comments = []
  positions = set()
  position_indices = (1,1) if diagonal else (1,2)

  with open(file) as f:

    for line in f:

      line = line.strip()

      if line.startswith('#'):
        comments += [re.sub(r'^#\s*', '', line)]
        continue

      if header:
        replicates = [
          re.sub(r'replicate\s*', '', i)
          for i in line.split('\t')[position_indices[1]+1:]
        ]
        header = False
        continue

      line = line.split('\t')
      chromosome = str(line[0])

      position_1, position_2 = [int(line[i]) for i in position_indices]

      def format(i):
        if i == 'None':
          return None
        try:
          return float(i)
        except:
          return i

      values = [format(i) for i in line[position_indices[1]+1:]]

      if chromosome not in sizes:
        sizes[chromosome] = 0

      if not replicates:
        replicates = ['1.'+str(i) for i in range(len(interactions))]

      positions |= set([position_1, position_2])

      sizes[chromosome] = max(
        position_2, sizes[chromosome]
      )

      if diagonal or sum(values) > 0:
        interactions[(chromosome, position_1, position_2)] = values

  resolution = min(abs(i - j) for i in positions for j in positions if i != j)

  for chromosome in sizes:
    sizes[chromosome] += resolution

  return dict(
    interactions = interactions,
    sizes = sizes,
    resolution = resolution,
    replicates = replicates,
    comments = comments
  )

# Create a diagonal dictionary from a matrix
#
# chromosome    position    replicate 1.1    ...
#
# {
#   entries: {
#             'chromosome 1': [
#                [position value, position value, ...],  # replicate 1
#                [position value, position value, ...],  # replicate 2
#                ...
#              ],
#              ...
#           },
#   sizes: {
#            'chromosome 1': 12780000,
#             ...
#          },
#   resolution: 10000,
#   replicates: ['1.1', '1.2', '2.1', '2.2'],
#   comments: ['# tissue: heart', '# normalization: cyclic loess']
# }
def matrix_to_diagonal(matrix):

  entries = {
    chromosome: [
      [
        matrix['interactions'][(chromosome, position, position)][i]
        if (chromosome, position, position) in matrix['interactions']
        else None
        for position in range(0, size, matrix['resolution'])
      ] for i in range(len(matrix['replicates']))
    ] for chromosome, size in matrix['sizes'].items()
  }

  bins = {
    chromosome: size // matrix['resolution']
    for chromosome, size in matrix['sizes'].items()
  }

  diagonal = dict(
    entries = entries,
    bins = bins,
    resolution = matrix['resolution'],
    replicates = matrix['replicates'],
    comments = matrix['comments']
  )

  if 'removed' in matrix:
    diagonal['removed'] = {
      chromosome: [
        position // matrix['resolution']
        for position in matrix['removed'][chromosome]
      ] for chromosome in matrix['removed']
    }

  return diagonal

# Fill a sparse matrix with zeros to creata a full matrix
def sparse_to_full_matrix(matrix):

  regions = [
    (chromosome, position_1, position_2)
    for chromosome, size in matrix['sizes'].items()
    for position_1 in range(0, size + 1, matrix['resolution'])
    for position_2 in range(0, size + 1, matrix['resolution'])
    if position_1 <= position_2
  ]

  interactions = {
    region: [
      matrix['interactions'][region]
      if region in matrix['interactions']
      else [0 for _ in matrix['replicates']]
    ]
    for region in regions
  }

  full_matrix = dict(
    interactions = interactions,
    sizes = matrix['sizes'],
    resolution = matrix['resolution'],
    replicates = matrix['replictates'],
    comments = matrix['comments']
  )

  if 'removed' in matrix:
    full_matrix['removed'] = matrix['removed']

  return full_matrix

# Remove zero lines from a full matrix to create a sparse matrix
def full_to_sparse_matrix(matrix):

  sparse_matrix = matrix.copy()

  for region, values in list(sparse_matrix['interactions'].items()):
    if type(values) is not list:
      values = [values]
    if sum(values) == 0:
      del sparse_matrix['interactions'][region]

  return sparse_matrix

# Split a matrix of multiple replicates
# into a list of single replicate sparse matrices
# [
#   {
#     interactions: {
#                    (chromosome, position 1, position 2): interaction,
#                    ...
#                   },
#     sizes: {
#              'chromosome 1': 12780000,
#               ...
#            },
#     resolution: 10000,
#     replicates: ['1.1']
#   },
#   ...
# ]
def split_matrix(matrix):

  matrices = [dict(
    interactions = {},
    sizes = matrix['sizes'],
    resolution = matrix['resolution'],
    replicates = [replicate],
    comments = matrix['comments']
  ) for replicate in matrix['replicates']]

  if 'removed' in matrix:
    for i in range(len(matrix['replicates'])):
      matrices[i]['removed'] = matrix['removed']

  for region, values in matrix['interactions'].items():
    for i, value in enumerate(values):
      if value != 0:
        matrices[i]['interactions'][region] = [value]

  return matrices

# Join single replicate sparse matrices
# into a sparse matrix of multiple replicates
def join_matrices(*matrices):

  interactions = {}
  replicates = []

  for i, matrix in enumerate(matrices):
    for region, values in matrix['interactions'].items():
      if region not in interactions:
        interactions[region] = [0 for _ in matrices]
      interactions[region][i] = values[0]
    replicates += matrices[i]['replicates']

  matrix = dict(
    interactions = interactions,
    sizes = matrices[0]['sizes'],
    resolution = matrices[0]['resolution'],
    replicates = replicates,
    comments = matrices[0]['comments']
  )

  if 'removed' in matrices[0]:
    matrix['removed'] = matrices[0]['removed']

  return matrix

# Write a matrix to a file
# # comments
# chromosome    position 1    position 2    replicate 1.1    ...
def export_matrix(matrix, file, header=True):

  with open(file, 'w') as output:

    if matrix['comments']:
      output.write(
        '\n'.join(
          ['# ' + comment for comment in matrix['comments']]
        ) + '\n'
      )

    if header:
      output.write(
        '\t'.join(
          ['chromosome', 'position 1', 'position 2']
          + ['replicate ' + i for i in matrix['replicates']]
        ) + '\n'
      )

    for region, values in sorted(matrix['interactions'].items()):
      if type(values) is not list:
        values = [values]
      if sum(values) > 0:
        output.write(
          '\t'.join(
          [
            str(i) for i in
            list(region) +
            [str(i)[-2:] == '.0' and str(i)[:-2] or str(i) for i in values]
          ]
        ) + '\n')

# Write a diagonal to a file
# # comments
# chromosome    position    replicate 1.1    ...
def export_diagonal(diagonal, file, header=True, name='position'):

  with open(file, 'w') as output:

    if diagonal['comments']:
      output.write(
        '\n'.join(
          ['# ' + comment for comment in diagonal['comments']]
        ) + '\n'
      )

    if header:
      output.write(
        '\t'.join(
          ['chromosome', name]
          + ['replicate ' + i for i in diagonal['replicates']]
        ) + '\n'
      )

    for chromosome in diagonal['entries']:
      for i in range(diagonal['bins'][chromosome]):
        output.write('\t'.join([
          str(a) for a in [
            chromosome,
            i*diagonal['resolution'],
          ] + [
            replicate[i]
            for replicate in diagonal['entries'][chromosome]
          ]
        ])+'\n')

# Create a vectors dictionary from a matrix dictionary
# {
#   interactions: {
#                   chromosome: [                            # 3D numpy array
#                                 [                            # replicate 1
#                                  [cell 1, cell 2, ...],       # bin 1
#                                  [cell 1, cell 2, ...],       # bin 2
#                                  [cell 1, cell 2, ...],       # bin 3
#                                  ...
#                                ],
#                                ...
#                              ],
#                    ...
#                 },
#   bins: {
#           chromosome 1: 2180,
#            ...
#         },
#   resolution: 10000,
#   replicates: ['1.1', '1.2', '2.1', '2.2']
# }
def matrix_to_vectors(matrix):

  bins = {
    chromosome: size // matrix['resolution']
    for chromosome, size in matrix['sizes'].items()
  }

  interactions = {
    chromosome: [
      [
        [
          0 for _ in range(bins)
        ] for _ in range(bins)
      ] for _ in matrix['replicates']
    ] for chromosome, bins in bins.items()
  }

  for chromosome in bins:
    for bin_1 in range(bins[chromosome]):
      for bin_2 in range(bin_1, bins[chromosome]):
        position_1 = bin_1 * matrix['resolution']
        position_2 = bin_2 * matrix['resolution']
        values = (
          matrix['interactions'][(chromosome, position_1, position_2)]
          if (chromosome, position_1, position_2) in matrix['interactions']
          else [0 for _ in matrix['replicates']]
        )
        for replicate, value in enumerate(values):
          interactions[chromosome][replicate][bin_1][bin_2] = value
          interactions[chromosome][replicate][bin_2][bin_1] = value

  vectors = dict(
    interactions = interactions,
    bins = bins,
    resolution = matrix['resolution'],
    replicates = matrix['replicates'],
    comments = matrix['comments']
  )

  if 'removed' in matrix:
    vectors['removed'] = {
      chromosome: [
        position // matrix['resolution']
        for position in matrix['removed'][chromosome]
      ] for chromosome in matrix['removed']
    }

  return vectors

# Convert vectors dictionary to matrix
def vectors_to_sparse_matrix(vectors):

  sizes = {
    chromosome: bins * vectors['resolution']
    for chromosome, bins in vectors['bins'].items()
  }

  interactions = {}

  for chromosome, replicates in vectors['interactions'].items():
    for bin_1 in range(vectors['bins'][chromosome]):
      for bin_2 in range(bin_1, vectors['bins'][chromosome]):
        values = [replicate[bin_1][bin_2] for replicate in replicates]
        if sum(values) > 0:
          position_1 = bin_1 * vectors['resolution']
          position_2 = bin_2 * vectors['resolution']
          region = (chromosome, position_1, position_2)
          interactions[region] = values

  matrix = dict(
    interactions = interactions,
    sizes = sizes,
    resolution = vectors['resolution'],
    replicates = vectors['replicates'],
    comments = vectors['comments']
  )

  if 'removed' in vectors:
    matrix['removed'] = {
      chromosome: [
        bin * vectors['resolution']
        for bin in vectors['removed'][chromosome]
      ] for chromosome in vectors['removed']
    }

  return matrix

def matrix_to_ccmaps(matrix):

  bins = {
    chromosome: size // matrix['resolution']
    for chromosome, size in matrix['sizes'].items()
  }

  interactions = {
    chromosome: [] for chromosome in matrix['sizes']
  }

  matrices = split_matrix(matrix)

  for single_matrix in matrices:
    positions_1, positions_2, values = ({
      chromosome: [] for chromosome in single_matrix['sizes']
    } for i in range(3))

    for region, value in sorted(single_matrix['interactions'].items()):
      value = value[0]
      if value:
        chromosome = region[0]
        positions_1[chromosome] += [region[1]]
        positions_2[chromosome] += [region[2]]
        values[chromosome] += [value]

    for chromosome in single_matrix['sizes']:
      interactions[chromosome] += [
        gmlib.importer.gen_map_from_locations_value(
          positions_1[chromosome],
          positions_2[chromosome],
          values[chromosome]
        )
      ]

  ccmaps = dict(
    interactions = interactions,
    bins = bins,
    resolution = matrix['resolution'],
    replicates = matrix['replicates'],
    comments = matrix['comments']
  )

  if 'removed' in matrix:
    ccmaps['removed'] = {
      chromosome: [
        position // matrix['resolution']
        for position in matrix['removed'][chromosome]
      ] for chromosome in matrix['removed']
    }

  return ccmaps

def ccmaps_to_sparse_matrix(ccmaps):

  sizes = {
    chromosome: bins * ccmaps['resolution']
    for chromosome, bins in ccmaps['bins'].items()
  }

  interactions = {}

  for chromosome in ccmaps['interactions']:
    for r, ccmap in enumerate(ccmaps['interactions'][chromosome]):
      ccmap.make_readable()
      for bin_1, column in enumerate(ccmap.matrix):
        for bin_2, value in enumerate(column):
          if bin_1 > bin_2:
            continue
          position_1 = bin_1 * ccmaps['resolution']
          position_2 = bin_2 * ccmaps['resolution']
          if (chromosome, position_1, position_2) not in interactions:
            interactions[(chromosome, position_1, position_2)] = [
              0 for _ in ccmaps['replicates']
            ]
          interactions[(chromosome, position_1, position_2)][r] = value

  matrix = dict(
    interactions = interactions,
    sizes = sizes,
    resolution = ccmaps['resolution'],
    replicates = ccmaps['replicates'],
    comments = ccmaps['comments']
  )

  if 'removed' in ccmaps:
    matrix['removed'] = {
      chromosome: [
        bin * ccmaps['resolution']
        for bin in ccmaps['removed'][chromosome]
      ] for chromosome in ccmaps['removed']
    }

  return full_to_sparse_matrix(matrix)

# Find weak rows and columns (of sum <= threshold) in a vectors dictionary
def find_weak_bins(vectors, threshold=0):

  weak = {
    chromosome: [set() for _ in vectors['replicates']]
    for chromosome in vectors['interactions']
  }

  for chromosome, values in vectors['interactions'].items():
    for r, replicate in enumerate(values):
      for v, vector in enumerate(replicate):
        if sum(vector) <= threshold:
          weak[chromosome][r] |= {v}

  return weak

# Remove weak rows and columns (of sum <= threshold) from a vectors dictionary
# Weak rows and columns found in one replicate are removed from all replicates
def filter_vectors(vectors, threshold=0):

  interactions = vectors['interactions'].copy()

  weak_bins = find_weak_bins(vectors, threshold)
  removed = {
    chromosome: reduce(set.union, replicates)
    for chromosome, replicates in weak_bins.items()
  }

  bins = {
    chromosome: vectors['bins'][chromosome] - len(removed[chromosome])
    for chromosome in vectors['bins']
  }

  for chromosome in removed:
    for bin in reversed(sorted(removed[chromosome])):
      interactions[chromosome] = np.delete(
        np.delete(
          interactions[chromosome],
          bin, 2
        ),
        bin, 1
      ).tolist()

  return dict(
    interactions = interactions,
    removed = removed,
    bins = bins,
    resolution = vectors['resolution'],
    replicates = vectors['replicates'],
    comments = vectors['comments']
  )

# Remove weak rows and columns (of sum <= threshold) from a matrix dictionary
# Weak rows and columns found in one replicate are removed from all replicates
def filter_matrix(matrix, threshold=0):
  return vectors_to_sparse_matrix(
    filter_vectors(matrix_to_vectors(matrix), threshold)
  )

# Fill removed rows and columns with 0, in a vectors dictionary
def refill_vectors(vectors):

  interactions = vectors['interactions'].copy()

  for chromosome, removed in vectors['removed'].items():
    for bin in sorted(removed):
      interactions[chromosome] = np.insert(
        np.insert(
          interactions[chromosome],
          bin, 0, 2
        ),
        bin, 0, 1
      ).tolist()

  bins = {
    chromosome: (
      vectors['bins'][chromosome] + len(vectors['removed'][chromosome])
      for chromosome in vectors['bins']
    )
  }

  return dict(
    interactions = interactions,
    bins = vectors['bins'],
    resolution = vectors['resolution'],
    replicates = vectors['replicates'],
    comments = vectors['comments']
  )

# Fill removed rows and columns with 0, in a matrix dictionary
def refill_matrix(matrix):
  return vectors_to_sparse_matrix(refill_vectors(matrix_to_vectors(matrix)))

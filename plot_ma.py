#!/usr/bin/env python3

import argparse
import numpy as np
import lib.parse_matrix as pm

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from statsmodels.nonparametric.smoothers_lowess import lowess

parser = argparse.ArgumentParser(description = 'Plot MA plots')
parser.add_argument('-i', nargs=2, required=True,
                    help='Input raw matrix and normalized matrix')
parser.add_argument('-p', required=True, help='Output figure prefix')
args = parser.parse_args()

raw_vectors = pm.matrix_to_vectors(pm.import_sparse_matrix(args.i[0]))
normalized_vectors = pm.matrix_to_vectors(pm.import_sparse_matrix(args.i[1]))

resolution = str(normalized_vectors['resolution'] // 1000) + 'k'
comments = '\n'.join(normalized_vectors['comments'])

for chromosome in raw_vectors['interactions']:

  ma_plots = [
    [
      [] for row in raw_vectors['replicates']
    ]
    for line in raw_vectors['replicates']
  ]

  for line in range(len(raw_vectors['replicates'])):
    for row in range(0, line):
      ma_plots[line][row] = np.array([
        (
          (raw_vectors['interactions'][chromosome][line][i][j]
           + raw_vectors['interactions'][chromosome][row][i][j])/2,
          (raw_vectors['interactions'][chromosome][line][i][j]
           - raw_vectors['interactions'][chromosome][row][i][j])
        )
        for i in range(raw_vectors['bins'][chromosome])
        for j in range(i, raw_vectors['bins'][chromosome])
      ])
    for row in range(line, len(raw_vectors['replicates'])):
      ma_plots[line][row] = np.array([
        (
          (normalized_vectors['interactions'][chromosome][line][i][j]
           + normalized_vectors['interactions'][chromosome][row][i][j])/2,
          (normalized_vectors['interactions'][chromosome][line][i][j]
           - normalized_vectors['interactions'][chromosome][row][i][j])
        )
        for i in range(normalized_vectors['bins'][chromosome])
        for j in range(i, normalized_vectors['bins'][chromosome])
      ])

  # MA is now of form
  # [                                                 ]
  #    [cond_rep_i                             ], ...
  #       [cond_rep_j                   ], ...
  #           [ (x,y), (x,y), ...  ]

  fig, axs = plt.subplots(
    nrows = len(raw_vectors['replicates']),
    ncols = len(raw_vectors['replicates']),
    figsize = (24, 34)
  )
  plt.subplots_adjust(
    left  = 0.06,
    right = 0.98,
    top = 0.97,
    bottom = 0.36
  )

  plt.figtext(
    x = 0.06,
    y = 0.33,
    s = 'resolution: ' + resolution + '\n'
    + 'chromosome: ' + chromosome + '\n' + comments,
    fontsize = 25,
    color = 'black',
    fontfamily = 'Open Sans Condensed',
    verticalalignment = 'top',
    horizontalalignment = 'left',
    multialignment = 'left'
  )

  fig.legend(
    [
      Line2D([0], [0], color = '#BE143C', lw = 2),
      Line2D([0], [0], color = '#143CAA', lw = 2)
    ],
    ['Before', 'After'],
    loc = (0, 0),
    bbox_to_anchor = (0.86, 0.306),
    prop = {
      'family': 'Open Sans Condensed',
      'size': 25
    },
    borderpad = 0,
    fancybox = False,
    edgecolor = 'white'
  )

  for i in range(0, len(raw_vectors['replicates'])):
    for j in range(0, len(raw_vectors['replicates'])):

      ax = axs[i][j]
      ax.spines['top'].set_visible(False)
      ax.spines['right'].set_visible(False)

      if i == 0:
        ax.set_xlabel(
          raw_vectors['replicates'][j],
          fontsize = 20,
          family = 'Open Sans Condensed'
        )
        ax.xaxis.set_label_coords(0.5, 1.24)
      if j == 0:
        ax.set_ylabel(
          raw_vectors['replicates'][i],
          fontsize = 20,
          family = 'Open Sans Condensed'
        )
        ax.yaxis.set_label_coords(-0.24, 0.5)

      if i != j:
        values = ma_plots[i][j]
        xs, ys = np.transpose(values)
        ax.scatter(xs, ys, color = 'gray')
        span = np.ptp(xs)
        fitted_xs, fitted_ys = np.transpose(np.unique(
          lowess(ys, xs, delta = 0.01*span),
          axis = 0
        ))
        color = '#BE143C' if i > j else '#143CAA'
        ax.plot(
          fitted_xs,
          fitted_ys,
          color = color,
          linewidth = 3
        )
      else:
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.tick_params(
          axis = 'both',
          which = 'both',
          length = 0
        )

  plt.savefig(
    args.p + '_' + resolution + '_chr' + chromosome + '.png',
    format = 'png',
    dpi = 300
  )

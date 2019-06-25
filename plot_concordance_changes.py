#!/usr/bin/env python3

import argparse
import math
import numpy as np
import lib.parse_matrix as pm

from plotly.io import write_image
import plotly.graph_objs as go

parser = argparse.ArgumentParser(
  description = 'Plot distance between centroids in both conditions'
)
parser.add_argument('-i', nargs = 2, required = True,
                    help = 'Input detected compartments and conrdance')
parser.add_argument('-p', required = True, help = 'Output figure prefix')
args = parser.parse_args()

compartments = pm.matrix_to_diagonal(
  pm.import_sparse_matrix(args.i[0], diagonal = True)
)

concordance = pm.matrix_to_diagonal(
  pm.import_sparse_matrix(args.i[1], diagonal = True)
)

resolution = str(compartments['resolution'] // 1000) + 'k'
comments = '<br>'.join(compartments['comments'])

# replicates = ['1.1', '2.2', '2.3', '3.1', '1.2', '2.1']
# conditions = ['1', '2', '3']
# indices = [[0, 4], [1, 2, 5], [3]]
conditions = sorted(set(r.split('.')[0] for r in compartments['replicates']))
indices = [
  [
    index for index, r in enumerate(compartments['replicates'])
    if r.split('.')[0] == condition
  ] for condition in conditions
]

def draw_cross(x, y, color):

  x_mean, y_mean = np.mean(x), np.mean(y)
  x_std, y_std = np.var(x)**0.5, np.var(y)**0.5

  covar = np.cov(x, y, ddof = 0)
  eigenvalues, eigenvectors = np.linalg.eigh(covar)
  # x_var, y_var = eigenvalues
  eigenvector = eigenvectors[0] if x_std < y_std else eigenvectors[1]
  theta = np.arctan(eigenvector[1] / eigenvector[0])

  return [
    dict(
      type = 'line',
      xref = 'x',
      yref = 'y',
      x0 = x_mean
      + x_std * math.cos(math.radians(0)) * math.cos(theta)
      - y_std * math.sin(math.radians(0)) * math.sin(theta),
      y0 = y_mean
      + y_std * math.sin(math.radians(0)) * math.cos(theta)
      + x_std * math.cos(math.radians(0)) * math.sin(theta),
      x1 = x_mean
      + x_std * math.cos(math.radians(180)) * math.cos(theta)
      - y_std * math.sin(math.radians(180)) * math.sin(theta),
      y1 = y_mean
      + y_std * math.sin(math.radians(180)) * math.cos(theta)
      + x_std * math.cos(math.radians(180)) * math.sin(theta),
      line = dict(
        color = color,
        width = 1
      )
    ),
    dict(
      type = 'line',
      xref = 'x',
      yref = 'y',
      x0 = x_mean
      + x_std * math.cos(math.radians(90)) * math.cos(theta)
      - y_std * math.sin(math.radians(90)) * math.sin(theta),
      y0 = y_mean
      + y_std * math.sin(math.radians(90)) * math.cos(theta)
      + x_std * math.cos(math.radians(90)) * math.sin(theta),
      x1 = x_mean
      + x_std * math.cos(math.radians(270)) * math.cos(theta)
      - y_std * math.sin(math.radians(270)) * math.sin(theta),
      y1 = y_mean
      + y_std * math.sin(math.radians(270)) * math.cos(theta)
      + x_std * math.cos(math.radians(270)) * math.sin(theta),
      line = dict(
        color = color,
        width = 1
      )
    )
  ]

def annotate(x, y, i, color):

  x_mean, y_mean = np.mean(x), np.mean(y)
  x_std, y_std = np.var(x)**0.5, np.var(y)**0.5

  covar = np.cov(x, y, ddof = 0)
  eigenvalues, eigenvectors = np.linalg.eigh(covar)
  # x_var, y_var = eigenvalues
  eigenvector = eigenvectors[0] if x_std < y_std else eigenvectors[1]
  theta = np.arctan(eigenvector[1] / eigenvector[0])

  return dict(
    text = i,
    font = dict(
      family = 'Open Sans Condensed',
      size = 22,
      color = color
    ),
    align = 'left',
    showarrow = False,
    x = x_mean
    + x_std/10 * math.cos(theta),
    y = y_mean
    + y_std/10 * math.cos(theta),
    xref = 'x',
    yref = 'y',
    xanchor = 'left',
    yanchor = 'bottom',
  )

for chromosome in compartments['bins']:

  xs = list(range(compartments['bins'][chromosome]))

  for a, b in [
    (a, b) for a in range(len(indices)) for b in range(a+1, len(indices))
  ]:

    compartment_changes = [
      i for i in range(compartments['bins'][chromosome])
      if (compartments['entries'][chromosome][indices[a][0]][i]
      != compartments['entries'][chromosome][indices[b][0]][i])
    ]

    c = np.array(concordance['entries'][chromosome], float)

    layout = go.Layout(
      margin = go.layout.Margin(
        l = 100,
        r = 320,
        b = 1300,
        t = 200
      ),
      showlegend = False,
      annotations = [
        dict(
          text = 'resolution: ' + resolution + '<br>' +
          'chromosome: ' + chromosome + '<br>'
          + 'conditions: ' + conditions[a] + ' vs ' + conditions[b] + '<br>'
          + comments,
          font = dict(
            family = 'Open Sans Condensed',
            size = 35
          ),
          align = 'left',
          showarrow = False,
          x = 0,
          y = 0,
          xref = 'paper',
          yref = 'paper',
          yshift = -120,
          xanchor = 'left',
          yanchor = 'top'
        ),
        dict(
          text = 'Concordance' + '<br>'
          + 'in condition '+conditions[b],
          font = dict(
            family = 'Open Sans Condensed',
            size = 35
          ),
          align = 'center',
          showarrow = False,
          x = 0,
          y = 1,
          xref = 'x',
          yref = 'paper',
          yshift = 100,
          xanchor = 'center',
          yanchor = 'top',
        ),
        dict(
          text = 'Concordance' + '<br>'
          + 'in condition '+conditions[a],
          font = dict(
            family = 'Open Sans Condensed',
            size = 35
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'middle',
        ),
        *[
          annotate(
            np.tile(c[indices[a], i], len(indices[b])),
            np.repeat(c[indices[b], i], len(indices[a])),
            i,
            'rgb(190, 20, 60)' if i in compartment_changes
            else 'rgb(20, 60, 170)'
          )
          for i in xs
        ]
      ],
      shapes = [
        bar for i in xs for bar in draw_cross(
          np.tile(c[indices[a], i], len(indices[b])),
          np.repeat(c[indices[b], i], len(indices[a])),
          'rgb(190, 20, 60)' if i in compartment_changes
          else 'rgb(20, 60, 170)'
        )
      ],
      xaxis = dict(
        range = [
          (np.min(c[indices[a]])*32)/31,
          (np.max(c[indices[a]])*32)/31
        ],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 20
        )
      ),
      yaxis = dict(
        range = [
          (np.min(c[indices[b]])*32)/31,
          (np.max(c[indices[b]])*32)/31
        ],
        scaleanchor = 'x',
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 20
        )
      )
    )

    figure = go.Figure(data = [], layout = layout)

    write_image(
      figure,
      args.p + '_' + resolution + '_chr'
      + chromosome + '_' + conditions[a] + 'vs' + conditions[b] + '.png',
      width = 2480,
      height = 3508
    )

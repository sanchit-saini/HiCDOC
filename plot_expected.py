#!/usr/bin/env python3

import argparse
import lib.parse_matrix as pm

from plotly.io import write_image
import plotly.graph_objs as go

parser = argparse.ArgumentParser(
 description = 'Plot distance regressions'
)
parser.add_argument('-i', required = True, help = 'Input expected values')
parser.add_argument('-p', required = True, help = 'Output figure prefix')
args = parser.parse_args()

expected = pm.matrix_to_diagonal(
  pm.import_sparse_matrix(args.i, diagonal = True)
)

resolution = str(expected['resolution'] // 1000) + 'k'
comments = '<br>'.join(expected['comments'])

for chromosome, bins in expected['bins'].items():

  entries = expected['entries']

  if not type(entries[chromosome][0]) == list:
    entries[chromosome] = [entries[chromosome]]

  max_x = (len(entries[chromosome][0]) + 2) * expected['resolution']
  max_y = (max([v for e in entries[chromosome] for v in e]) * 16)/15

  layout = go.Layout(
    margin = go.layout.Margin(
      l = 150,
      r = 300,
      b = 2200,
      t = 150
    ),
    showlegend = False,
    annotations = [
      dict(
        text = 'resolution: ' + resolution + '<br>'
          + 'chromosome: ' + chromosome + '<br>' + comments,
        font = dict(
          family = 'Open Sans Condensed',
          size = 27
        ),
        align = 'left',
        showarrow = False,
        x = 0,
        y = 0,
        yshift = -100,
        xanchor = 'left',
        yanchor = 'top',
      ),
      dict(
        text = 'Genomic distance',
        font = dict(
          family = 'Open Sans Condensed',
          size = 27
        ),
        align = 'left',
        showarrow = False,
        x = max_x,
        y = 0,
        xref = 'x',
        yref = 'y',
        xshift = 20,
        xanchor = 'left',
        yanchor = 'bottom',
      ),
      dict(
        text = 'Expected<br>interaction proportion',
        font = dict(
          family = 'Open Sans Condensed',
          size = 27
        ),
        align = 'center',
        showarrow = False,
        x = 0,
        y = max_y,
        xref = 'x',
        yref = 'y',
        yshift = 20,
        xanchor = 'center',
        yanchor = 'bottom',
      )
    ],
    xaxis = dict(
      domain = [0, 0],
      range = [0, max_x],
      ticklen = 8,
      tickfont = dict(
        family = 'Open Sans Condensed',
        size = 19
      ),
    ),
    yaxis = dict(
      domain = [0, 0],
      range = [0, max_y],
      ticklen = 8,
      tickfont = dict(
        family = 'Open Sans Condensed',
        size = 19
      ),
    )
  )

  data = [
    go.Scatter(
      x = [
        i * expected['resolution']
        for i in range(0, len(entries[chromosome][0]))
      ],
      y = e,
      mode = 'lines',
      line = dict(
        width = 2,
        color = 'rgb(190, 20, 60)'
      ),
      xaxis = 'x',
      yaxis = 'y',
    ) for e in entries[chromosome]
  ]

  figure = go.Figure(data = data, layout = layout)

  write_image(
    figure,
    args.p + '_' + resolution + '_chr' + chromosome + '.png',
    width = 2480,
    height = 3508
  )

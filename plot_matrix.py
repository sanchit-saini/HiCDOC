#!/usr/bin/env python3

import argparse
import lib.parse_matrix as pm

from plotly.io import write_image
import plotly.graph_objs as go

parser = argparse.ArgumentParser(
  description = 'Plot measure under matrix'
)
parser.add_argument('-i', required = True, help = 'Input matrix')
parser.add_argument('-p', required = True, help = 'Output figure prefix')
parser.add_argument('--measure', help = 'Input measure')
parser.add_argument('--name', default = '', help = 'Measure name')
args = parser.parse_args()

vectors = pm.matrix_to_vectors(pm.import_sparse_matrix(args.i))
if args.measure:
  measure = pm.matrix_to_diagonal(
    pm.import_sparse_matrix(args.measure, diagonal = True)
  )

resolution = str(vectors['resolution'] // 1000) + 'k'
comments = '<br>'.join(vectors['comments'])

for chromosome in vectors['interactions']:

  for r, values in enumerate(vectors['interactions'][chromosome]):

    replicate = vectors['replicates'][r]

    annotations = [
      dict(
        text = 'resolution: ' + resolution + '<br>' +
        'chromosome: ' + chromosome + '<br>' +
        'replicate: ' + replicate + '<br>' + comments,
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
        yshift = -50,
        xanchor = 'left',
        yanchor = 'top',
      )
    ]

    yaxes = dict(
      yaxis = dict(
        domain = [0.1 if args.measure else 0, 1],
        scaleanchor = 'x',
        autorange = 'reversed',
        visible = False
      )
    )

    data = [
      go.Heatmap(
        z = values,
        showscale = False,
        colorscale = [
           [0, 'rgba(255, 255, 255, 0)'],
           [1/256, 'rgba(255, 250, 250, 1)'],
           [1/128, 'rgba(255, 240, 240, 1)'],
           [1/64, 'rgba(255, 230, 230, 1)'],
           [1/32, 'rgba(255, 220, 220, 1)'],
           [1/16, 'rgba(255, 190, 190, 1)'],
           [1/8, 'rgba(255, 160, 160, 1)'],
           [1/4, 'rgba(255, 80, 80, 1)'],
           [1/2, 'rgba(255, 50, 50, 1)'],
           [1, 'rgba(255, 0, 0, 1)']
         ]
      )
    ]

    if args.measure:
      annotations += [dict(
        text = args.name,
        font = dict(
          family = 'Open Sans Condensed',
          size = 35
        ),
        align = 'left',
        showarrow = False,
        x = len(values)-1,
        y = 0,
        xref = 'x',
        yref = 'y2',
        xshift = 30,
        xanchor = 'left',
        yanchor = 'middle',
      )]
      yaxes = dict(
        **yaxes,
        yaxis2 = dict(
          domain = [0, 0.1],
          ticklen = 10,
          tickfont = dict(
            family = 'Open Sans Condensed',
            size = 20
          )
        )
      )
      data += [
        go.Scatter(
          x = list(range(len(values))),
          y = measure['entries'][chromosome][r],
          mode = 'lines',
          line = dict(
            width = 2,
            color = 'rgb(20, 60, 170)'
          ),
          xaxis = 'x',
          yaxis = 'y2',
        )
      ]

    layout = go.Layout(
      margin = go.layout.Margin(
        l = 100 if args.measure else 70,
        r = 250 if args.measure else 70,
        b = 950,
        t = 70
      ),
      showlegend = False,
      annotations = annotations,
      xaxis = dict(
        domain = [0, 1],
        range = [0, len(values)-1],
        visible = False
      ),
      **yaxes
    )

    figure = go.Figure(data = data, layout = layout)

    write_image(
      figure,
      args.p + '_' + resolution + '_chr' + chromosome + '_' + replicate + '.png',
      width = 2480,
      height = 3508
    )

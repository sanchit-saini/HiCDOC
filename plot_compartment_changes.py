#!/usr/bin/env python3

import argparse
import numpy as np
import lib.parse_matrix as pm

from plotly.io import write_image
import plotly.graph_objs as go

parser = argparse.ArgumentParser(
  description = 'Plot compartment changes and measures'
)
parser.add_argument('-i', required = True, help = 'Input detected compartments')
parser.add_argument('-p', required = True, help = 'Output figure prefix')
parser.add_argument('--concordance', help = 'Input concordance')
parser.add_argument('--distances', nargs = '+',
                    help = 'Input distances to centroids. '
                           'One file per centroid')
parser.add_argument('--silhouette', help = 'Input Silhouette')
args = parser.parse_args()

compartments = pm.matrix_to_diagonal(
  pm.import_sparse_matrix(args.i, diagonal = True)
)

if args.concordance:
  concordance = pm.matrix_to_diagonal(
    pm.import_sparse_matrix(args.concordance, diagonal = True)
  )

if args.distances:
  distances = [
    pm.matrix_to_diagonal(
      pm.import_sparse_matrix(f, diagonal = True)
    )
    for f in args.distances
  ]

if args.silhouette:
  silhouette = pm.matrix_to_diagonal(
    pm.import_sparse_matrix(args.silhouette, diagonal = True)
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

    compartment_changes_labels = []

    i = 0
    while i < len(compartment_changes):
      start = i
      end = i
      increase = True
      while increase and end < len(compartment_changes) - 1:
        increase = False
        max_step = 3 if compartment_changes[end] > 99 else 2
        for step in range(1, max_step+1):
          if compartment_changes[end + 1] == compartment_changes[end] + step:
            end += 1
            increase = True
            break

      if start == end:
        label = str(compartment_changes[start])
        position = compartment_changes[start]
      else:
        label = str(compartment_changes[start]) + '-' + str(compartment_changes[end])
        position = compartment_changes[start] + (
          compartment_changes[end] - compartment_changes[start]
        )/2

      compartment_changes_labels += [(position, label)]
      i = end + 1


    layout = go.Layout(
      margin = go.layout.Margin(
        l = 120,
        r = 500,
        b = 400,
        t = 120
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
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 0,
          y = 0,
          xref = 'paper',
          yref = 'paper',
          yshift = -70,
          xanchor = 'left',
          yanchor = 'top'
        ),
        dict(
          text = 'Compartment<br>in condition '+conditions[a],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        dict(
          text = 'Compartment<br>in condition '+conditions[b],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y2',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        dict(
          text = 'Concordance<br>in condition '+conditions[a],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y3',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'middle',
        ),
        dict(
          text = 'Concordance<br>in condition '+conditions[b],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y4',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'middle',
        ),
        dict(
          text = 'Silhouette coefficient<br>in condition '+conditions[a],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y5',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        dict(
          text = 'Silhouette coefficient<br>in condition '+conditions[b],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y6',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        dict(
          text = 'Euclidian distance' + '<br>'
          + 'to first centroid' + '<br>'
          + 'in condition '+conditions[a],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y7',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        dict(
          text = 'Euclidian distance' + '<br>'
          + 'to second centroid' + '<br>'
          + 'in condition '+conditions[a],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y8',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        dict(
          text = 'Euclidian distance' + '<br>'
          + 'to first centroid' + '<br>'
          + 'in condition '+conditions[b],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y9',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        dict(
          text = 'Euclidian distance' + '<br>'
          + 'to second centroid' + '<br>'
          + 'in condition '+conditions[b],
          font = dict(
            family = 'Open Sans Condensed',
            size = 45
          ),
          align = 'left',
          showarrow = False,
          x = 1,
          y = 0,
          xref = 'paper',
          yref = 'y10',
          xshift = 30,
          xanchor = 'left',
          yanchor = 'bottom',
        ),
        *[
          dict(
            text = label,
            font = dict(
              family = 'Open Sans Condensed',
              size = 35
            ),
            align = 'left',
            showarrow = False,
            x = i,
            y = 1,
            xref = 'x',
            yref = 'paper',
            yshift = 15,
            xanchor = 'center',
            yanchor = 'bottom',
          ) for i, label in compartment_changes_labels
        ]
      ],
      shapes = [
        *[
          dict(
            type = 'line',
            xref = 'x',
            yref = 'paper',
            x0 = i,
            y0 = 0,
            x1 = i,
            y1 = 1,
            line = dict(
              color = 'rgb(190, 20, 60)',
              width = 2,
              dash = 'dash'
            )
          ) for i in compartment_changes
        ]
      ],
      xaxis = dict(
        domain = [0, 1],
        visible = False,
        range = [0, xs[-1]]
      ),
      yaxis = dict(
        domain = [0.96, 1],
        dtick = 1,
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        ),
        showgrid = False,
        zeroline = False
      ),
      yaxis2 = dict(
        domain = [0.9, 0.94],
        dtick = 1,
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        ),
        showgrid = False,
        zeroline = False
      ),
      yaxis3 = dict(
        domain = [0.77, 0.87],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      ),
      yaxis4 = dict(
        domain = [0.65, 0.75],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      ),
      yaxis5 = dict(
        domain = [0.53, 0.63],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      ),
      yaxis6 = dict(
        domain = [0.41, 0.51],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      ),
      yaxis7 = dict(
        domain = [0.3, 0.39],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      ),
      yaxis8 = dict(
        domain = [0.2, 0.29],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      ),
      yaxis9 = dict(
        domain = [0.1, 0.19],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      ),
      yaxis10 = dict(
        domain = [0, 0.09],
        ticklen = 10,
        tickfont = dict(
          family = 'Open Sans Condensed',
          size = 25
        )
      )
    )

    data = [
      go.Scatter(
        x = xs,
        y = compartments['entries'][chromosome][indices[a][0]],
        mode = 'lines',
        line = dict(
          width = 2,
          color = 'rgb(20, 60, 170)'
        ),
        xaxis = 'x',
        yaxis = 'y',
      ),
      go.Scatter(
        x = xs,
        y = compartments['entries'][chromosome][indices[b][0]],
        mode = 'lines',
        line = dict(
          width = 2,
          color = 'rgb(20, 60, 170)'
        ),
        xaxis = 'x',
        yaxis = 'y2',
      )
    ]

    if args.concordance:
      c = np.array(concordance['entries'][chromosome], float)
      data += [
        *[
          go.Scatter(
            x = xs,
            y = ys,
            mode = 'lines',
            line = dict(
              width = 2,
              color = 'rgb(20, 60, 170)'
            ),
            xaxis = 'x',
            yaxis = 'y3',
          ) for ys in c[indices[a]].tolist()
        ],
        *[
          go.Scatter(
            x = xs,
            y = ys,
            mode = 'lines',
            line = dict(
              width = 2,
              color = 'rgb(20, 60, 170)'
            ),
            xaxis = 'x',
            yaxis = 'y4',
          ) for ys in c[indices[b]].tolist()
        ],
      ]

    if args.silhouette:
      s = np.array(silhouette['entries'][chromosome], float)
      data += [
        *[
          go.Scatter(
            x = xs,
            y = ys,
            mode = 'lines',
            line = dict(
              width = 2,
              color = 'rgb(20, 60, 170)'
            ),
            xaxis = 'x',
            yaxis = 'y5',
          ) for ys in s[indices[a]].tolist()
        ],
        *[
          go.Scatter(
            x = xs,
            y = ys,
            mode = 'lines',
            line = dict(
              width = 2,
              color = 'rgb(20, 60, 170)'
            ),
            xaxis = 'x',
            yaxis = 'y6',
          ) for ys in s[indices[b]].tolist()
        ],
      ]

    if args.distances:
      d = [np.array(d['entries'][chromosome], float) for d in distances]
      data += [
        *[
          go.Scatter(
            x = xs,
            y = ys,
            mode = 'lines',
            line = dict(
              width = 2,
              color = 'rgb(20, 60, 170)'
            ),
            xaxis = 'x',
            yaxis = 'y7',
          ) for ys in d[0][indices[a]].tolist()
        ],
        *[
          go.Scatter(
            x = xs,
            y = ys,
            mode = 'lines',
            line = dict(
              width = 2,
              color = 'rgb(20, 60, 170)'
            ),
            xaxis = 'x',
            yaxis = 'y9',
          ) for ys in d[0][indices[b]].tolist()
        ]
      ]
      if len(args.distances) > 1:
        data += [
          *[
            go.Scatter(
              x = xs,
              y = ys,
              mode = 'lines',
              line = dict(
                width = 2,
                color = 'rgb(20, 60, 170)'
              ),
              xaxis = 'x',
              yaxis = 'y8',
            ) for ys in d[1][indices[a]].tolist()
          ],
          *[
            go.Scatter(
              x = xs,
              y = ys,
              mode = 'lines',
              line = dict(
                width = 2,
                color = 'rgb(20, 60, 170)'
              ),
              xaxis = 'x',
              yaxis = 'y10',
            ) for ys in d[1][indices[b]].tolist()
          ]
        ]

    figure = go.Figure(data = data, layout = layout)

    write_image(
      figure,
      args.p + '_' + resolution + '_chr'
      + chromosome + '_' + conditions[a] + 'vs' + conditions[b] + '.png',
      width = 7016,
      height = 4960
    )

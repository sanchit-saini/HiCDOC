#!/usr/bin/env python3

import argparse
import lib.parse_matrix as pm

parser = argparse.ArgumentParser(
  description = 'Join multiple replicates into one matrix file'
)
parser.add_argument(
  '-i',
  nargs = '+',
  required = True,
  help = 'Input matrices'
)
parser.add_argument(
  '-o',
  required = True,
  help = 'Output matrix'
)
parser.add_argument(
  '--inputs-have-headers',
  dest = 'header',
  action = 'store_true',
  help = 'Whether the input matrices have a header line. '
         'Omit if they do not have a header line'
)
parser.add_argument(
  '--replicates',
  nargs = '+',
  required = True,
  help = 'condition.replicate value for each input matrix. '
         'Example: 1.1 1.2 2.1 2.2'
)
parser.add_argument(
  '--comments',
  nargs = '+',
  default = '',
  help = 'Comments, in quotes, separated by spaces. '
         'Example: "tissue: muscle" "days: 12"'
)
args = parser.parse_args()

matrices = [
  pm.import_sparse_matrix(
    input,
    header = args.header
  ) for input in args.i
]
matrix = pm.join_matrices(*matrices)

matrix['replicates'] = args.replicates
matrix['comments'] = args.comments

pm.export_matrix(matrix, args.o)

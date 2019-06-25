#!/bin/bash

PARAMS=""
while (( "$#" )); do
  # split --option=something into 2 variables
  [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
  case "$1" in
    -i)
      input="$2"
      shift 2
      ;;
    -d)
      outdir="${2%/}"
      shift 2
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS \"$1\""
      shift
      ;;
  esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

scriptdir="$(dirname "$(readlink -f "$0")")"
mkdir -p "$outdir"

printf "\n\e[1;32mNormalizing technical biases with cyclic loess\e[0m\n"

"$scriptdir"/normalize_cyclic_loess.r \
  -i "$input" \
  -o "$outdir"/normalized.tsv

printf "\n\e[1;32mNormalizing biological biases with Knight-Ruiz\e[0m\n"

"$scriptdir"/normalize_knight_ruiz.py \
  -i "$outdir"/normalized.tsv \
  -o "$outdir"/normalized.tsv

printf "\n\e[1;32mNormalizing distance effect with combined RNR\e[0m\n"

"$scriptdir"/normalize_distance_rnr_combined.py \
  -i "$outdir"/normalized.tsv \
  -o "$outdir"/normalized.tsv

printf "\n\e[1;32mDetecting compartments\e[0m\n"

"$scriptdir"/detect_constrained_k_means.py \
  -i "$outdir"/normalized.tsv \
  -o "$outdir"/compartments.tsv \
  --concordance "$outdir"/concordance.tsv \
  --silhouette "$outdir"/silhouette.tsv \
  --distances "$outdir"/distance1.tsv "$outdir"/distance2.tsv

printf "\n\e[1;32mPlotting compartment changes\e[0m\n"

"$scriptdir"/plot_compartment_changes.py \
  -i "$outdir"/compartments.tsv \
  -p "$outdir"/compartments \
  --concordance "$outdir"/concordance.tsv \
  --silhouette "$outdir"/silhouette.tsv \
  --distances "$outdir"/distance1.tsv "$outdir"/distance2.tsv

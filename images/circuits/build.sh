#!/usr/bin/env bash
# Render every circuit drawing in this folder (each *.tex but style.tex) to a PNG
# beside its source. Needs pdflatex (with circuitikz and quantikz) and pdftoppm.
# Usage: images/circuits/build.sh [name ...]   (no names: render them all)
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="$(mktemp -d)"
trap 'rm -r "$build_dir"' EXIT

if [[ $# -gt 0 ]]; then
  names=("$@")
else
  names=()
  for source in "$here"/*.tex; do
    name="$(basename "$source" .tex)"
    [[ "$name" == style ]] || names+=("$name")
  done
fi

for name in "${names[@]}"; do
  (cd "$here" && pdflatex -interaction=nonstopmode -halt-on-error \
    -output-directory "$build_dir" "$name.tex" >"$build_dir/$name.out") || {
    echo "LaTeX failed on $name.tex; the end of its log:" >&2
    tail -n 25 "$build_dir/$name.log" >&2
    exit 1
  }
  pdftoppm -png -r 300 -singlefile "$build_dir/$name.pdf" "$here/$name"
  echo "rendered $name.png"
done

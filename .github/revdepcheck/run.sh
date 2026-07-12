#!/usr/bin/env bash
set -euo pipefail

set +e
timeout -s TERM -k 5m 350m Rscript .github/revdepcheck/run.R
code=$?
set -e

case "$code" in
  124|130|137|143)
    echo "again=true" >> "$GITHUB_OUTPUT"
    exit 0
    ;;
  0)
    echo "again=false" >> "$GITHUB_OUTPUT"
    exit 0
    ;;
  *)
    echo "again=false" >> "$GITHUB_OUTPUT"
    exit "$code"
    ;;
esac

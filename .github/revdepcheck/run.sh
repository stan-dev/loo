#!/usr/bin/env bash
set -euo pipefail

set +e
timeout -s INT -k 5m 350m Rscript .github/revdepcheck/run.R
code=$?
set -e

if [ "$code" = "124" ] || [ "$code" = "130" ] || [ "$code" = "137" ]; then
  echo "again=true" >> "$GITHUB_OUTPUT"
  exit 0
fi

echo "again=false" >> "$GITHUB_OUTPUT"
exit "$code"

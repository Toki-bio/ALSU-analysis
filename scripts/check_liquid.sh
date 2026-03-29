#!/usr/bin/env bash
set -euo pipefail

# Quick check for Liquid tokens and common placeholder-style braces in Markdown files.
# Exit non-zero if matches are found so hooks/CI can fail fast.

md_files=$(git ls-files '*.md')
if [ -z "$md_files" ]; then
  exit 0
fi

# Patterns: Liquid openers and common placeholder tokens (e.g. {N}, {LOCALE})
liquid_re='\\{\\%|\\{\\{'
placeholder_re='\\{[A-Z0-9_:\\-]+\\}'

matches=$(echo "$md_files" | xargs -r grep -n -E "${liquid_re}|${placeholder_re}" || true)
if [ -n "$matches" ]; then
  echo "Potential Liquid/template tokens found in Markdown files:" >&2
  echo "$matches" >&2
  echo >&2
  echo "Wrap literal template code in {% raw %}...{% endraw %} or replace placeholders with safe tokens (e.g. %N% or <<N>>)." >&2
  exit 1
fi

exit 0

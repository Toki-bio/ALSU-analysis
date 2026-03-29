#!/usr/bin/env bash
set -euo pipefail

# Configure repository to use .githooks (run once locally)
git config core.hooksPath .githooks
echo "Configured git hooks path to .githooks."

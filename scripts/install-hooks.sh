#!/bin/sh
# Point git at the versioned hooks in .githooks/.
#
# core.hooksPath is LOCAL config -- it is not carried by a clone -- so every fresh
# clone must run this once:
#
#     sh scripts/install-hooks.sh
#
# See .githooks/pre-commit for what is enforced and why.
set -eu

root=$(git rev-parse --show-toplevel)
cd "$root"

git config core.hooksPath .githooks
chmod +x .githooks/* 2>/dev/null || true   # no-op on Windows checkouts

printf 'core.hooksPath = %s\n' "$(git config --get core.hooksPath)"
printf 'installed hooks: %s\n' "$(ls .githooks | tr '\n' ' ')"

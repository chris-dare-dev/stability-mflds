#!/bin/sh
# Configure this clone's LOCAL git settings.  Run once per clone:
#
#     sh scripts/setup-clone.sh
#
# Both settings live in .git/config, which `git clone` does NOT carry, so a fresh clone
# starts without them.  Neither can be committed.
#
#   1. core.hooksPath -> .githooks
#      Activates .githooks/pre-commit, which keeps working artifacts (planning docs,
#      .claude/, figures/) out of the repository.
#
#   2. filter.obsidian-strip
#      This repo sits inside an Obsidian vault that stamps YAML frontmatter into
#      CLAUDE.md, README.md and docs/*.md.  The clean filter strips it on `git add`, so
#      the working tree keeps the metadata for Obsidian while git stores the clean,
#      LF-normalized version.  See .gitattributes and scripts/strip-obsidian-frontmatter.py.
#
# Without (2) the pre-commit hook refuses any document that still carries a frontmatter
# block, so a forgotten setup fails loudly instead of polluting the history.
set -eu

root=$(git rev-parse --show-toplevel)
cd "$root"

git config core.hooksPath .githooks
chmod +x .githooks/* 2>/dev/null || true   # no-op on Windows checkouts

# Prefer python3 where both exist; Git for Windows ships `python`.
py=python
if command -v python3 >/dev/null 2>&1; then
	py=python3
fi

git config filter.obsidian-strip.clean "$py scripts/strip-obsidian-frontmatter.py"
git config filter.obsidian-strip.smudge cat

# Re-run the filter over already-tracked content so the index matches the new rules.
git add --renormalize . >/dev/null 2>&1 || true

printf 'core.hooksPath               = %s\n' "$(git config --get core.hooksPath)"
printf 'filter.obsidian-strip.clean  = %s\n' "$(git config --get filter.obsidian-strip.clean)"
printf 'filter.obsidian-strip.smudge = %s\n' "$(git config --get filter.obsidian-strip.smudge)"
printf 'hooks: %s\n' "$(ls .githooks | tr '\n' ' ')"
printf 'index renormalized.\n'

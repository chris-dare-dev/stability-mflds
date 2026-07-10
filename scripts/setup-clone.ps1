# Configure this clone's LOCAL git settings  (PowerShell equivalent of setup-clone.sh).
#
#     powershell -ExecutionPolicy Bypass -File scripts\setup-clone.ps1
#
# Sets core.hooksPath (activating .githooks/pre-commit) and the obsidian-strip clean
# filter (keeping vault frontmatter out of git while the working tree keeps it).  Both
# live in .git/config, which `git clone` does not carry, so run this once per clone.
#
# The hooks are POSIX sh; Git for Windows ships the shell that runs them.
# See .githooks/pre-commit and .gitattributes for what is enforced and why.

$ErrorActionPreference = 'Stop'

$root = (git rev-parse --show-toplevel)
Set-Location $root

git config core.hooksPath .githooks

# `python` is what Git for Windows puts on PATH; prefer python3 if it resolves.
$py = 'python'
if (Get-Command python3 -ErrorAction SilentlyContinue) { $py = 'python3' }

git config filter.obsidian-strip.clean "$py scripts/strip-obsidian-frontmatter.py"
git config filter.obsidian-strip.smudge cat

# Re-run the filter over already-tracked content so the index matches the new rules.
git add --renormalize . 2>$null

Write-Output "core.hooksPath               = $(git config --get core.hooksPath)"
Write-Output "filter.obsidian-strip.clean  = $(git config --get filter.obsidian-strip.clean)"
Write-Output "filter.obsidian-strip.smudge = $(git config --get filter.obsidian-strip.smudge)"
Write-Output "hooks: $((Get-ChildItem .githooks | Select-Object -ExpandProperty Name) -join ' ')"
Write-Output "index renormalized."

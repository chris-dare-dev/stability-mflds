# Point git at the versioned hooks in .githooks/  (PowerShell equivalent of install-hooks.sh).
#
# core.hooksPath is LOCAL config -- it is not carried by a clone -- so every fresh
# clone must run this once:
#
#     powershell -ExecutionPolicy Bypass -File scripts\install-hooks.ps1
#
# The hooks are POSIX sh; Git for Windows ships the shell that runs them.
# See .githooks/pre-commit for what is enforced and why.

$ErrorActionPreference = 'Stop'

$root = (git rev-parse --show-toplevel)
Set-Location $root

git config core.hooksPath .githooks

Write-Output "core.hooksPath = $(git config --get core.hooksPath)"
Write-Output "installed hooks: $((Get-ChildItem .githooks | Select-Object -ExpandProperty Name) -join ' ')"

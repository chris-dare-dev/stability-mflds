@echo off
rem Bridge Windows -> WSL Macaulay2 (E10-M4 / G16).
rem
rem The oracle (bridgeland_stability/oracle/m2.py) discovers M2 via the
rem BRIDGELAND_M2 env var and invokes it as  M2 --script <windows-temp-path>.
rem Macaulay2 has no native Windows build; it lives in the WSL Debian distro
rem (apt package `macaulay2`).  This shim forwards the invocation, translating
rem any argument that names an existing Windows file into its /mnt/... WSL
rem path via wslpath.  Everything else (flags like --script) passes through.
rem
rem Opt-in usage (the suite's @requires_m2 tests skip without it):
rem   set BRIDGELAND_M2=<repo>\scripts\m2-wsl.cmd
rem   pytest -q
setlocal enabledelayedexpansion
set "ARGS="
for %%A in (%*) do (
  set "ARG=%%~A"
  if exist "%%~A" (
    for /f "usebackq delims=" %%P in (`wsl -d Debian -- wslpath -a "%%~fA"`) do set "ARG=%%P"
  )
  set ARGS=!ARGS! "!ARG!"
)
wsl -d Debian -- M2 %ARGS%

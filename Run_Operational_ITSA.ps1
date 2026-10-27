$ErrorActionPreference = 'Stop'
$taskProject = 'C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Simulation'
if ((Resolve-Path -LiteralPath $PSScriptRoot).Path -ne $taskProject) { throw 'Script must remain in Simulation.' }
Set-Location -LiteralPath $taskProject
if ((Get-Location).Path -ne $taskProject) { throw 'Wrong directory.' }
if (!(Test-Path .\.venv\Scripts\python.exe) -or !(Test-Path .\itsa_recovery_v5\run_analysis.py)) { throw 'Environment or V5 source missing.' }
& .\.venv\Scripts\python.exe -B .\itsa_recovery_v5\run_analysis.py --stage pilot
if ($LASTEXITCODE -ne 0) { throw "V5 stopped with exit code $LASTEXITCODE. Review processing.log and run_state.json." }

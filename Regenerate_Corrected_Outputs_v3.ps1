$ErrorActionPreference = 'Stop'
$taskProject = 'C:\Users\Arlan\OneDrive\Desktop\Jhonjhon\Remote Sensing Research\Simulation'
if ((Resolve-Path -LiteralPath $PSScriptRoot).Path -ne $taskProject) { throw 'Script must remain in Simulation.' }
Set-Location -LiteralPath $taskProject
if ((Get-Location).Path -ne $taskProject) { throw 'Wrong workspace.' }
$taskPython = Join-Path $taskProject '.venv\Scripts\python.exe'
$taskEntry = Join-Path $taskProject 'itsa_recovery_v3\run_analysis.py'
if (!(Test-Path -LiteralPath $taskPython) -or !(Test-Path -LiteralPath $taskEntry)) { throw 'Python or V3 entry point missing.' }
$taskActive = Join-Path $taskProject 'Corrected_outputs'
if (Test-Path -LiteralPath $taskActive) {
    if ((Resolve-Path -LiteralPath $taskActive).Path -ne "$taskProject\Corrected_outputs") { throw 'Unexpected output path.' }
    $taskArchiveRoot = Join-Path $taskProject 'Output_Archive'
    $taskArchive = Join-Path $taskArchiveRoot ('Corrected_outputs_' + (Get-Date -Format 'yyyyMMdd_HHmmss_ffff'))
    if (![IO.Path]::GetFullPath($taskArchive).StartsWith("$taskProject\Output_Archive\", [StringComparison]::OrdinalIgnoreCase)) { throw 'Unsafe archive path.' }
    if (Test-Path -LiteralPath $taskArchive) { throw 'Archive already exists.' }
    New-Item -ItemType Directory -Path $taskArchiveRoot -Force | Out-Null
    Copy-Item -LiteralPath $taskActive -Destination $taskArchive -Recurse
    foreach ($taskFile in Get-ChildItem -LiteralPath $taskActive -File -Recurse) {
        $taskRelative = $taskFile.FullName.Substring($taskActive.Length).TrimStart('\')
        $taskCopy = Join-Path $taskArchive $taskRelative
        if (!(Test-Path -LiteralPath $taskCopy)) { throw "Archive missing $taskRelative" }
        if ((Get-FileHash -LiteralPath $taskFile.FullName -Algorithm SHA256).Hash -ne (Get-FileHash -LiteralPath $taskCopy -Algorithm SHA256).Hash) { throw "Archive mismatch: $taskRelative" }
    }
    if ((Resolve-Path -LiteralPath $taskActive).Path -ne "$taskProject\Corrected_outputs") { throw 'Output path changed.' }
    Remove-Item -LiteralPath $taskActive -Recurse -Force
    Write-Host "Previous outputs preserved and verified: $taskArchive"
}
& $taskPython -B $taskEntry --stage pilot
if ($LASTEXITCODE -ne 0) { throw "V3 failed with exit code $LASTEXITCODE. Archived results remain available." }

param(
    [string]$CsvPath
)

$ErrorActionPreference = "Stop"

$rscript = Join-Path $env:USERPROFILE "anaconda3\envs\rstudio\Scripts\Rscript.exe"
$fishDir = $PSScriptRoot

if (-not $CsvPath) {
    $csvMatch = Get-ChildItem -Path "J:\OneDrive -*\*\salmon\salmon2025-2\salmon2025-2_result\speed_heading_data2025_2_trout.csv" -ErrorAction Stop |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if (-not $csvMatch) {
        throw "Could not locate speed_heading_data2025_2_trout.csv"
    }
    $CsvPath = $csvMatch.FullName
}

$dataDir = Split-Path -Parent $CsvPath
$csvName = Split-Path -Leaf $CsvPath

$workRoot = Join-Path $env:TEMP "salmonBehav2025_run"
$fishLink = Join-Path $workRoot "fish_link"
$dataLink = Join-Path $workRoot "data_link"
$outDir = Join-Path $workRoot "out"

New-Item -ItemType Directory -Force -Path $workRoot | Out-Null
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

if (Test-Path $fishLink) {
    cmd /c rmdir "$fishLink" | Out-Null
}
if (Test-Path $dataLink) {
    cmd /c rmdir "$dataLink" | Out-Null
}

New-Item -ItemType Junction -Path $fishLink -Target $fishDir | Out-Null
New-Item -ItemType Junction -Path $dataLink -Target $dataDir | Out-Null

$scriptPath = Join-Path $fishLink "salmonBehav2025.R"
$csvPath = Join-Path $dataLink $csvName

& $rscript $scriptPath $csvPath $outDir
if ($LASTEXITCODE -ne 0) {
    throw "R script failed with exit code $LASTEXITCODE"
}

$finalSummary = Join-Path $dataLink "salmonBehav2025_summary.csv"
$finalPairwise = Join-Path $dataLink "salmonBehav2025_pairwise_tests.csv"

Copy-Item -LiteralPath (Join-Path $outDir "salmonBehav2025_summary.csv") -Destination $finalSummary -Force
Copy-Item -LiteralPath (Join-Path $outDir "salmonBehav2025_pairwise_tests.csv") -Destination $finalPairwise -Force

Write-Host "Completed:"
Write-Host (Join-Path $dataDir "salmonBehav2025_summary.csv")
Write-Host (Join-Path $dataDir "salmonBehav2025_pairwise_tests.csv")

$ErrorActionPreference = "Stop"

$deploymentDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoDir = Split-Path -Parent $deploymentDir
$distributionDir = Join-Path $repoDir "dist\TranscriptScan"
$versionLine = (Get-Content -LiteralPath (Join-Path $deploymentDir "VERSION.txt") -Raw).Trim()

if ($versionLine -notmatch '^Transcript Scan (?<Version>\d+\.\d+\.\d+)$') {
    throw "VERSION.txt must contain 'Transcript Scan X.Y.Z'."
}

$version = $Matches.Version
$releaseName = "TranscriptScan-$version-Windows-x64"
$distRoot = [System.IO.Path]::GetFullPath((Join-Path $repoDir "dist"))
$releaseDir = [System.IO.Path]::GetFullPath((Join-Path $distRoot $releaseName))
$zipPath = [System.IO.Path]::GetFullPath((Join-Path $distRoot "$releaseName.zip"))
$expectedPrefix = $distRoot.TrimEnd('\') + '\'

if (-not $releaseDir.StartsWith($expectedPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to package outside the repository dist folder: $releaseDir"
}

$requiredFiles = @(
    "TranscriptScan.exe",
    "README_TRANSCRIPT_SCAN.txt",
    "THIRD_PARTY_NOTICES.txt",
    "VERSION.txt",
    "multiple_sequence_blast_template.xlsx"
)

if (-not (Test-Path -LiteralPath (Join-Path $distributionDir "_internal") -PathType Container)) {
    throw "Build output is missing dist\TranscriptScan\_internal. Run build_transcript_scan.ps1 first."
}
foreach ($name in $requiredFiles) {
    if (-not (Test-Path -LiteralPath (Join-Path $distributionDir $name) -PathType Leaf)) {
        throw "Build output is missing required file: $name"
    }
}

if (Test-Path -LiteralPath $releaseDir) {
    Remove-Item -LiteralPath $releaseDir -Recurse -Force
}
New-Item -ItemType Directory -Path $releaseDir | Out-Null

Copy-Item -LiteralPath (Join-Path $distributionDir "_internal") -Destination $releaseDir -Recurse
foreach ($name in $requiredFiles) {
    Copy-Item -LiteralPath (Join-Path $distributionDir $name) -Destination $releaseDir
}

$unexpected = Get-ChildItem -LiteralPath $releaseDir -Recurse -Force | Where-Object {
    $_.Name -eq "TranscriptScanData" -or
    $_.Name -eq "settings.json" -or
    $_.Extension -eq ".log" -or
    $_.Name -like "phase*_manual*"
}
if ($unexpected) {
    throw "Release staging contains generated user/test data: $($unexpected.FullName -join ', ')"
}

if (Test-Path -LiteralPath $zipPath) {
    Remove-Item -LiteralPath $zipPath -Force
}
Compress-Archive -LiteralPath $releaseDir -DestinationPath $zipPath -CompressionLevel Optimal

$zipFile = Get-Item -LiteralPath $zipPath
$zipHash = Get-FileHash -LiteralPath $zipPath -Algorithm SHA256
Write-Output "RELEASE_VERSION=$version"
Write-Output "RELEASE_DIR=$releaseDir"
Write-Output "ZIP_PATH=$zipPath"
Write-Output "ZIP_SIZE_BYTES=$($zipFile.Length)"
Write-Output "ZIP_SHA256=$($zipHash.Hash)"

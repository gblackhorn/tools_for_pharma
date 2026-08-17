$ErrorActionPreference = "Stop"

$deploymentDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoDir = Split-Path -Parent $deploymentDir
Set-Location -LiteralPath $repoDir

python -m PyInstaller --noconfirm --clean deployment\transcript_scan.spec
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller build failed with exit code $LASTEXITCODE."
}

$distributionDir = Join-Path $repoDir "dist\TranscriptScan"
Copy-Item -LiteralPath "deployment\README_TRANSCRIPT_SCAN.txt" -Destination $distributionDir -Force
Copy-Item -LiteralPath "deployment\VERSION.txt" -Destination $distributionDir -Force
Copy-Item -LiteralPath "deployment\THIRD_PARTY_NOTICES.txt" -Destination $distributionDir -Force

Write-Host "Built portable app: $distributionDir"
Write-Host "Run packaged verification:"
Write-Host "  Start-Process -FilePath '$distributionDir\TranscriptScan.exe' -ArgumentList '--self-test' -Wait -PassThru"

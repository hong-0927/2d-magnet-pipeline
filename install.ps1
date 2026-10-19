# One-shot setup for Windows: virtual environment, dependencies, environment check.
# Usage:  powershell -ExecutionPolicy Bypass -File install.ps1
#         powershell -ExecutionPolicy Bypass -File install.ps1 -Venv C:\path\to\venv
param(
    [string]$Venv = (Join-Path $PSScriptRoot '.venv')
)
$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot

# The Microsoft Store stub answers as "python" and then opens the store; use py when present.
$py = 'py'
if (-not (Get-Command py -ErrorAction SilentlyContinue)) { $py = 'python' }
& $py -c "import sys; assert sys.version_info >= (3, 9), sys.version" 2>$null
if ($LASTEXITCODE -ne 0) { throw "[error] Python >= 3.9 required (py -3 or python3)" }

Write-Host "[1/3] venv -> $Venv"
& $py -m venv $Venv

$activate = Join-Path $Venv 'Scripts\Activate.ps1'
if (-not (Test-Path $activate)) { throw "[error] venv activation script missing: $activate" }
& $activate

Write-Host "[2/3] dependencies"
& python -m pip install --upgrade pip | Out-Null
& python -m pip install -r (Join-Path $PSScriptRoot 'requirements.txt')
if ($LASTEXITCODE -ne 0) { throw "[error] pip install failed" }

Write-Host "[3/3] doctor"
& python (Join-Path $PSScriptRoot 'run_2d_magnet.py') doctor
if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Warning "doctor reported WARN lines. They are environment hints, not code errors."
    Write-Host "         Symptom -> fix table: knowledge/KNOWN_ISSUES.md"
} else {
    Write-Host ""
    Write-Host "[OK] environment clean. Next: README.md 'Quick start'."
}

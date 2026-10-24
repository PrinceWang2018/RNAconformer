# Build the RCV Windows installer.
#   Requirements: 64-bit Python 3.12 (python.org or `winget install Python.Python.3.12`),
#                 Inno Setup 6 (`winget install JRSoftware.InnoSetup`) for the installer step.
#   Run from the repository root in PowerShell:   powershell -ExecutionPolicy Bypass -File packaging\build_windows.ps1
param([string]$Python = "py -3.12", [switch]$SkipInstaller)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent $PSScriptRoot)

Write-Host "== creating build environment (.venv-build)"
Invoke-Expression "$Python -m venv .venv-build"
$py = ".venv-build\Scripts\python.exe"
& $py -m pip install --upgrade pip
& $py -m pip install -r packaging\requirements-desktop.txt
& $py -m pip install --no-deps -e .

Write-Host "== self-check (engine + ViennaRNA)"
& $py -c "import RNA, numpy, scipy, webview; from rcv import __version__; print('RCV', __version__, 'ViennaRNA', RNA.__version__)"

Write-Host "== PyInstaller"
& $py -m PyInstaller packaging\rcv.spec --noconfirm --distpath dist --workpath build\pyinstaller

Write-Host "== smoke test of the frozen engine (windowed exe: output goes to %APPDATA%\RCV\logs)"
$p = Start-Process -FilePath "dist\RCV\RCV.exe" -ArgumentList "engine","ss","--sequence","GGCACUUCGGUGCC" -Wait -PassThru
if ($p.ExitCode -ne 0) { throw "frozen engine self-test failed (exit code $($p.ExitCode))" }
$p = Start-Process -FilePath "dist\RCV\RCV.exe" -ArgumentList "convert","examples\viewer_demo\HIV1_TAR_NMR_1ANR.pdb","-o","build\selftest.rcv" -Wait -PassThru
if ($p.ExitCode -ne 0 -or -not (Test-Path "build\selftest.rcv")) { throw "frozen .rcv self-test failed (exit code $($p.ExitCode))" }

Write-Host "== help pages (user guide, plugin guide, what is new)"
& $py packaging\make_docs.py build\docs_html
if ($LASTEXITCODE -ne 0) { throw "make_docs.py failed" }

if (-not $SkipInstaller) {
  $iscc = "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe"
  if (-not (Test-Path $iscc)) { $iscc = "${env:ProgramFiles}\Inno Setup 6\ISCC.exe" }
  if (-not (Test-Path $iscc)) { $iscc = "${env:LOCALAPPDATA}\Programs\Inno Setup 6\ISCC.exe" }
  if (Test-Path $iscc) {
    $version = & $py -c "from rcv import __version__; print(__version__)"
    & $iscc "/DAppVersion=$version" packaging\windows\RCV.iss
    Write-Host "== installer written to dist\installer"
  } else {
    Write-Warning "Inno Setup 6 not found: the portable app is in dist\RCV (zip it to share)."
  }
}

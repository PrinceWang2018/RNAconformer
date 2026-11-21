# Build the RCV MSIX package for the Microsoft Store (and for testing it on your own PC).
#   Requirements: the PyInstaller build in dist\RCV and the help pages in build\docs_html (run
#                 packaging\build_windows.ps1 first, -SkipInstaller is enough), and the Windows 10/11 SDK
#                 (MakeAppx.exe, MakePri.exe; `winget install Microsoft.WindowsSDK.10.0.26100`).
#   Run from the repository root in PowerShell:
#     powershell -File packaging\windows\build_msix.ps1 `
#       -IdentityName "12345ZixiangWang.RCV" -Publisher "CN=XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX" `
#       -PublisherDisplayName "Zixiang Wang"
#   The three values come from Microsoft Partner Center > the app > Product management > Product identity and must
#   match it exactly. Output: dist\msix\RCV-<version>-Windows-x64.msix (unsigned: upload it to Partner Center,
#   the Store signs it).
#
#   To install the package on your own PC before submitting it, sign it with a test certificate whose subject is
#   the -Publisher value:  -CertificatePfx test.pfx -CertificatePassword ...  (trust the certificate once in
#   "Local Machine > Trusted People"), then double-click the .msix. Without -Publisher a local test identity
#   (CN=RCV Test) is used.
param(
  [string]$IdentityName = "RCV.RNAConformationViewer",
  [string]$Publisher = "CN=RCV Test",
  [string]$PublisherDisplayName = "Zixiang Wang",
  [string]$Version = "",
  [string]$CertificatePfx = "",
  [string]$CertificatePassword = ""
)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path -Parent (Split-Path -Parent $PSScriptRoot))

if (-not (Test-Path "dist\RCV\RCV.exe")) { throw "dist\RCV\RCV.exe not found: run packaging\build_windows.ps1 first" }
if (-not (Test-Path "build\docs_html\index.html")) { throw "build\docs_html not found: run packaging\build_windows.ps1 first" }

# ---- version: a.b.c.0 (the Store needs the fourth number to be 0)
if (-not $Version) {
  $Version = (Get-Content "src\rcv\__init__.py" | Select-String "__version__ = '([^']+)'").Matches[0].Groups[1].Value
}
$nums = @(($Version -replace '[^0-9.].*$', '').Split('.') | Where-Object { $_ -ne '' } | ForEach-Object { [int]$_ })
while ($nums.Count -lt 3) { $nums += 0 }
$msixVersion = "{0}.{1}.{2}.0" -f $nums[0], $nums[1], $nums[2]
Write-Host "== RCV $Version -> MSIX version $msixVersion ($IdentityName, $Publisher)"

# ---- Windows SDK tools (newest SDK first)
function Find-SdkTool([string]$name) {
  $roots = @("${env:ProgramFiles(x86)}\Windows Kits\10\bin", "${env:ProgramFiles}\Windows Kits\10\bin")
  foreach ($r in $roots) {
    if (-not (Test-Path $r)) { continue }
    $hit = Get-ChildItem $r -Directory | Where-Object { $_.Name -match '^10\.' } | Sort-Object { [version]$_.Name } -Descending |
      ForEach-Object { Join-Path $_.FullName "x64\$name" } | Where-Object { Test-Path $_ } | Select-Object -First 1
    if ($hit) { return $hit }
  }
  $cmd = Get-Command $name -ErrorAction SilentlyContinue
  if ($cmd) { return $cmd.Source }
  throw "$name not found: install the Windows 10/11 SDK (winget install Microsoft.WindowsSDK.10.0.26100)"
}
$makeappx = Find-SdkTool "makeappx.exe"
$makepri = Find-SdkTool "makepri.exe"

# ---- package folder: the program, the help pages, licences, the manifest and the icons
$stage = "build\msix\RCV"
if (Test-Path "build\msix") { Remove-Item "build\msix" -Recurse -Force }
New-Item -ItemType Directory -Force $stage | Out-Null
Copy-Item "dist\RCV\*" $stage -Recurse
Copy-Item "build\docs_html" "$stage\docs" -Recurse
Copy-Item "LICENSE" "$stage\LICENSE.txt"
Copy-Item "packaging\THIRD-PARTY-NOTICES.txt" $stage
Copy-Item "packaging\windows\msix\Assets" "$stage\Assets" -Recurse

$esc = { param($s) [System.Security.SecurityElement]::Escape($s) }
$manifest = (Get-Content "packaging\windows\msix\AppxManifest.xml.in" -Raw -Encoding UTF8).
  Replace('@IDENTITY_NAME@', (& $esc $IdentityName)).
  Replace('@PUBLISHER@', (& $esc $Publisher)).
  Replace('@PUBLISHER_DISPLAY@', (& $esc $PublisherDisplayName)).
  Replace('@VERSION@', $msixVersion)
[System.IO.File]::WriteAllText((Resolve-Path $stage).Path + "\AppxManifest.xml", $manifest, (New-Object System.Text.UTF8Encoding $false))

# ---- resources.pri: lets Windows pick the icon of the right size (scale-100/200/400, target sizes 16-256)
Write-Host "== MakePri"
& $makepri createconfig /cf build\msix\priconfig.xml /dq en-US /pv 10.0.0 /o | Out-Null
if ($LASTEXITCODE -ne 0) { throw "MakePri createconfig failed" }
& $makepri new /pr $stage /cf build\msix\priconfig.xml /mn "$stage\AppxManifest.xml" /of "$stage\resources.pri" /o
if ($LASTEXITCODE -ne 0) { throw "MakePri new failed" }

# ---- pack
New-Item -ItemType Directory -Force "dist\msix" | Out-Null
$out = "dist\msix\RCV-$Version-Windows-x64.msix"
Write-Host "== MakeAppx -> $out"
& $makeappx pack /d $stage /p $out /o
if ($LASTEXITCODE -ne 0) { throw "MakeAppx pack failed" }

# ---- optional: sign for testing on this PC (the Store signs the submitted package itself)
if ($CertificatePfx) {
  $signtool = Find-SdkTool "signtool.exe"
  $sargs = @("sign", "/fd", "SHA256", "/f", $CertificatePfx)
  if ($CertificatePassword) { $sargs += @("/p", $CertificatePassword) }
  & $signtool @sargs $out
  if ($LASTEXITCODE -ne 0) { throw "signtool failed (the certificate subject must equal -Publisher)" }
  Write-Host "== signed with $CertificatePfx"
}
$mb = [math]::Round((Get-Item $out).Length / 1MB, 1)
Write-Host "== MSIX package written: $out ($mb MB)"

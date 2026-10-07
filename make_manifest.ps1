Set-Location D:\理论体系\repo
$lines = Get-ChildItem -Recurse -File |
  Where-Object { $_.Name -ne 'MANIFEST.sha256' -and $_.Name -ne 'MANIFEST.tmp' -and $_.Name -ne 'make_manifest.ps1' } |
  Sort-Object FullName |
  Get-FileHash -Algorithm SHA256 |
  ForEach-Object {
    $rel = $_.Path.Substring((Get-Location).Path.Length + 1) -replace '\\','/'
    "$($_.Hash.ToLower())  $rel"
  }
[IO.File]::WriteAllLines("$pwd\MANIFEST.sha256", $lines, (New-Object System.Text.UTF8Encoding($false)))
Write-Host "完成:共 $((Get-Content MANIFEST.sha256).Count) 行"

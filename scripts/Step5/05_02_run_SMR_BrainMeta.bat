@echo off
setlocal enabledelayedexpansion

REM =========================================================
REM Step 5.2
REM Run SMR analysis using BrainMeta cis-eQTL data
REM =========================================================

REM Project root: this BAT file is stored in scripts\Step5\
set PROJECT_ROOT=%~dp0..\..
for %%I in ("%PROJECT_ROOT%") do set PROJECT_ROOT=%%~fI

set SMR=%PROJECT_ROOT%\tools\smr\smr.exe
set GWAS=%PROJECT_ROOT%\results\Step5\01_prepare_gwas\Step5_ADHD2022_for_SMR.txt
set LD=%PROJECT_ROOT%\data\Step5\LD_reference\data_maf0.01_rs_ref
set SNP=%PROJECT_ROOT%\results\Step2\Step2_panel_snps_rsid.txt
set EQTL_DIR=%PROJECT_ROOT%\data\Step5\BrainMeta
set OUTDIR=%PROJECT_ROOT%\results\Step5\02_smr_raw

if not exist "%OUTDIR%" mkdir "%OUTDIR%"

echo ========================================================
echo Step 5.2: Run SMR with BrainMeta cis-eQTL
echo ========================================================
echo PROJECT: %PROJECT_ROOT%
echo SMR    : %SMR%
echo GWAS   : %GWAS%
echo LD     : %LD%
echo SNP    : %SNP%
echo EQTL   : %EQTL_DIR%
echo OUTDIR : %OUTDIR%
echo.

if not exist "%SMR%" (
  echo ERROR: SMR executable not found.
  pause
  exit /b 1
)

if not exist "%GWAS%" (
  echo ERROR: GWAS summary file not found.
  pause
  exit /b 1
)

if not exist "%LD%.bed" (
  echo ERROR: LD reference BED file not found.
  pause
  exit /b 1
)

if not exist "%LD%.bim" (
  echo ERROR: LD reference BIM file not found.
  pause
  exit /b 1
)

if not exist "%LD%.fam" (
  echo ERROR: LD reference FAM file not found.
  pause
  exit /b 1
)

if not exist "%SNP%" (
  echo ERROR: SNP list file not found.
  pause
  exit /b 1
)

for %%c in (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22) do (
  echo Running chromosome %%c ...

  if exist "%EQTL_DIR%\BrainMeta_cis_eQTL_chr%%c.besd" (
    "%SMR%" ^
      --bfile "%LD%" ^
      --gwas-summary "%GWAS%" ^
      --beqtl-summary "%EQTL_DIR%\BrainMeta_cis_eQTL_chr%%c" ^
      --extract-snp "%SNP%" ^
      --diff-freq 0.2 ^
      --out "%OUTDIR%\Step5_SMR_chr%%c"
  ) else (
    echo WARNING: BrainMeta_cis_eQTL_chr%%c.besd not found, skipped.
  )

  echo.
)

echo All chromosomes finished.
pause
endlocal
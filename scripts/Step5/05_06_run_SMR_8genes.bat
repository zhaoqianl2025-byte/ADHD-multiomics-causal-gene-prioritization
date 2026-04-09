@echo off
chcp 65001 >nul
setlocal

REM =========================================================
REM Step 5.6
REM Run SMR plot jobs for 8 selected genes
REM =========================================================

REM This script is located in scripts\Step5\
set PROJECT_ROOT=%~dp0..\..
for %%I in ("%PROJECT_ROOT%") do set PROJECT_ROOT=%%~fI

set SMR=%PROJECT_ROOT%\tools\smr\smr.exe
set GWAS=%PROJECT_ROOT%\results\Step5\01_prepare_gwas\Step5_ADHD2022_for_SMR.txt
set BFILE=%PROJECT_ROOT%\data\Step5\EU1000QC\EU1000QC
set GENE_LIST=%PROJECT_ROOT%\data\Step5\glist-hg19.txt
set EQTL_DIR=%PROJECT_ROOT%\data\Step5\BrainMeta
set OUTROOT=%PROJECT_ROOT%\results\Step5\05_plots\genes

if not exist "%OUTROOT%" mkdir "%OUTROOT%"

echo ============================================
echo Step 5.6: Running SMR plots for 8 selected genes
echo ============================================
echo PROJECT  : %PROJECT_ROOT%
echo SMR      : %SMR%
echo GWAS     : %GWAS%
echo BFILE    : %BFILE%
echo GENE LIST: %GENE_LIST%
echo EQTL DIR : %EQTL_DIR%
echo OUTROOT  : %OUTROOT%
echo.

REM ================= 检查 =================
if not exist "%SMR%" (
  echo ERROR: SMR executable not found.
  pause & exit /b 1
)

if not exist "%GWAS%" (
  echo ERROR: GWAS summary file not found.
  pause & exit /b 1
)

if not exist "%BFILE%.bed" (
  echo ERROR: BFILE BED not found.
  pause & exit /b 1
)

if not exist "%GENE_LIST%" (
  echo ERROR: gene list file not found.
  pause & exit /b 1
)

REM =========================================================
REM FUNCTION: run one gene
REM =========================================================
set i=0

call :run_gene 01_TMEM161B-AS1 TMEM161B-AS1 ENSG00000247828.8 chr5
call :run_gene 02_AC091826.2   AC091826.2   ENSG00000271904.1 chr5
call :run_gene 03_LSM6         LSM6         ENSG00000164167.12 chr4
call :run_gene 04_AC097372.2   AC097372.2   ENSG00000279845.1 chr4
call :run_gene 05_KIFC2        KIFC2        ENSG00000167702.13 chr8
call :run_gene 06_REELD1       REELD1       ENSG00000250673.2 chr4
call :run_gene 07_PTPRF        PTPRF        ENSG00000142949.17 chr1
call :run_gene 08_ZNF34        ZNF34        ENSG00000196378.13 chr8

echo.
echo ============================================
echo All jobs finished.
echo ============================================
pause
exit /b

REM =========================================================
REM Subroutine
REM =========================================================
:run_gene
set /a i+=1

set FOLDER=%1
set SYMBOL=%2
set PROBE=%3
set CHR=%4

echo [%i%/8] Running %SYMBOL% ...

set OUTDIR=%OUTROOT%\%FOLDER%
if not exist "%OUTDIR%" mkdir "%OUTDIR%"

pushd "%OUTDIR%"

"%SMR%" ^
  --bfile "%BFILE%" ^
  --gwas-summary "%GWAS%" ^
  --beqtl-summary "%EQTL_DIR%\BrainMeta_cis_eQTL_%CHR%" ^
  --probe %PROBE% ^
  --probe-wind 1000 ^
  --gene-list "%GENE_LIST%" ^
  --plot ^
  --out plot_%SYMBOL%

popd
echo.
exit /b
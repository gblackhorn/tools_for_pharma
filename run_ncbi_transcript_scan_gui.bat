@echo off
cd /d "%~dp0"
python -m tools_for_pharma.oligo.ncbi_blast --gui
if errorlevel 1 pause

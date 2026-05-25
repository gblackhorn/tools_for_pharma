@echo off
cd /d "%~dp0"
python qpcr_ref_qc_extract.py --gui
if errorlevel 1 pause

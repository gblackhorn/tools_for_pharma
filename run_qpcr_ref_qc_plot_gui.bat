@echo off
cd /d "%~dp0"
python -m tools_for_pharma.qpcr.ref_qc_plot --gui
if errorlevel 1 pause

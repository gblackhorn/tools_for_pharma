@echo off
cd /d "%~dp0"
python qpcr_ref_qc_plot.py --gui
if errorlevel 1 pause

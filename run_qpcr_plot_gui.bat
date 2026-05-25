@echo off
cd /d "%~dp0"
python qPCR_plot_excel_table.py --gui
if errorlevel 1 pause

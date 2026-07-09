@echo off
cd /d "%~dp0"
python -m tools_for_pharma.plotting.IC50 --gui
if errorlevel 1 pause

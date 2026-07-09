@echo off
cd /d "%~dp0"
python -m tools_for_pharma.plotting.curve --gui
if errorlevel 1 pause

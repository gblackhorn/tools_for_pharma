@echo off
cd /d "%~dp0"
python -m tools_for_pharma.oligo.off_target --gui
if errorlevel 1 pause

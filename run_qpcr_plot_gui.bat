@echo off
cd /d "%~dp0"
python -m tools_for_pharma.qpcr.plot --gui
if errorlevel 1 pause

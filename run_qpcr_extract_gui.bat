@echo off
cd /d "%~dp0"
python -m tools_for_pharma.qpcr.extract --gui
if errorlevel 1 pause

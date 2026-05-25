@echo off
cd /d "%~dp0"
python oligo_utils.py --gui
if errorlevel 1 pause

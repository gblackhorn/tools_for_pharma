@echo off
cd /d "%~dp0"
python util_oligo_offtarget.py --gui
if errorlevel 1 pause

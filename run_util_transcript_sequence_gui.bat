@echo off
cd /d "%~dp0"
python util_transcript_sequence.py --gui
if errorlevel 1 pause

Transcript Scan 1.0.0
=====================

Requirements
------------
- Windows 10 or Windows 11, 64-bit
- No Python installation is required
- Network access to NCBI is needed only when a requested transcript is not cached

Starting the app
----------------
1. Keep the complete TranscriptScan folder together.
2. Double-click TranscriptScan.exe.
3. On first use, enter your own contact email for NCBI transcript requests.

Local data
----------
The app creates a TranscriptScanData folder beside TranscriptScan.exe:

- settings.json stores the NCBI contact email entered by the user.
- transcript_cache stores public transcript FASTA records for reuse.
- logs stores diagnostic logs.

Oligo sequences are not written to the transcript cache. Local transcript scans
do not send oligo sequences to NCBI. Only public transcript accessions are sent
when a transcript needs to be downloaded or refreshed.

Updating the app
----------------
Preserve the TranscriptScanData folder when replacing the application files so
the saved email and transcript cache remain available.

Troubleshooting
---------------
If the app reports a failure, review:

TranscriptScanData\logs\transcript_scan.log

If Windows blocks the executable, contact your IT department rather than
disabling company security controls.

Transcript Scan 1.1.0
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

Single-sequence targets
-----------------------
The single-sequence workflow can use an exact NM/XM/NR/XR transcript accession,
a pasted transcript sequence, or one local FASTA/text file. Pasted and local-file
targets do not use NCBI. Genomic NC accessions are not whole-transcript records
and are not accepted in transcript-accession mode.

After reviewing a result, choose Edit and run again to return to the populated
form. These form values stay only until the app is closed. Refresh from NCBI is
a one-time action and resets after a successful run.

Local data
----------
The app creates a TranscriptScanData folder beside TranscriptScan.exe:

- settings.json stores the NCBI contact email entered by the user.
- transcript_cache stores public transcript FASTA records for reuse.
- logs stores diagnostic logs.

Pasted target sequences and selected local files are not copied into settings or
the transcript cache.

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

# tools_for_pharma

Utilities for pharma-oriented oligo sequence work and qPCR report processing.

## Repository Layout

```text
tools_for_pharma/
  shared/        Generic helpers used by multiple workflows
  oligo/         Oligo/off-target/transcript sequence tools
  qpcr/          qPCR extraction, plotting, and reference-gene QC tools
data/examples/  Development and example input workbooks
outputs/plots/  Generated plot examples
tests/fixtures/ Test data, when tests are added
```

## Logseq Summary Export

Keep the local Logseq Markdown mirror under `Logseq_markdown`, with source notes
in its `journals` and `pages` subfolders. These source files are the editable
source of truth; do not manually maintain the generated summary documents.

After refreshing the mirror, first perform a read-only comparison:

```powershell
python -m tools_for_pharma.logseq_export --source Logseq_markdown --check
```

The check lists source paths that were added, modified, deleted, or changed
between ID-only and substantive content. It does not print note contents. Exit
status `0` means both exports are current, while status `1` means an update is
needed.

Update only the stale journal or page export with:

```powershell
python -m tools_for_pharma.logseq_export --source Logseq_markdown --update
```

Omitting `--update` retains the original update behavior. Each affected summary
is regenerated completely and replaced atomically in `Logseq_markdown/export`:

- `Logseq_journals_for_summary.md`
- `Logseq_pages_for_summary.md`

The exporter removes Logseq `id::` UUID properties, removes SharePoint and
OneDrive destinations, masks standalone DNA/RNA strings of at least 10 bases,
preserves source boundaries, and lists ID-only files in an appendix. Re-run
`--check` after updating to verify that both summaries are synchronized.

## MSH3 Transcript Panel Comparison

Download the requested versioned Ensembl MSH3 mature cDNAs and compare each one
against the exact RefSeq reference `NM_002439.5`:

```powershell
python -m tools_for_pharma.oligo.transcript_panel --output-dir outputs\msh3_transcript_comparison
```

The output directory contains separate cDNA and CDS FASTA files, one combined
multi-record cDNA FASTA, a retrieval manifest, a reference-anchored comparison
summary, difference blocks, conserved coordinate blocks, and JSON data for a
review workbook. The command stops instead of silently substituting a newer
Ensembl transcript version when an exact requested version is not available
from the live Ensembl release.

The batch files stay in the repo root so they remain easy to double-click. Each
launcher changes into the repo root and runs the matching Python module with
`python -m ...`.

## Oligo Tools

Open the combined oligo GUI:

```text
run_oligo_gui.bat
```

Or run it from PowerShell:

```powershell
python -m tools_for_pharma.oligo.app --gui
```

### Off-Target Analysis Preparation

This workflow starts from an oligo antisense sequence. It extracts the selected
antisense region, usually positions **2-18**, and computes the complementary
sense sequence in **5'->3'** orientation.

Process one antisense sequence:

```powershell
python -m tools_for_pharma.oligo.off_target "AUGCUACGGAUCUAGCUAGCU"
```

Use a non-default antisense region:

```powershell
python -m tools_for_pharma.oligo.off_target "AUGCUACGGAUCUAGCUAGCU" --start 3 --end 19
```

Process many antisense sequences from an Excel/CSV table:

```powershell
python -m tools_for_pharma.oligo.off_target --input oligos.xlsx --column antisense
```

Open the off-target table GUI:

```text
run_util_oligo_offtarget_gui.bat
```

### Transcript Sequence Extraction

This workflow starts from a FASTA/plain-text transcript. It extracts a 1-based
inclusive transcript range and outputs the matched **SS** and **AS** strands for
oligo design.

Extract one transcript range:

```powershell
python -m tools_for_pharma.oligo.transcript_sequence --transcript-file transcript.fasta --start 120 --end 140
```

Match multiple ranges from a table with `start` and `end` columns:

```powershell
python -m tools_for_pharma.oligo.transcript_sequence --transcript-file transcript.fasta --range-table ranges.xlsx
```

Open the transcript sequence GUI:

```text
run_util_transcript_sequence_gui.bat
```

### NCBI Transcript / BLAST Checks

For AS oligo checks against a specific transcript, provide the AS sequence and
an NM/XM/NR/XR accession. The tool fetches the transcript through NCBI EFetch and
scans for the AS reverse-complement target.

Open the private local transcript-scan GUI:

```text
run_ncbi_transcript_scan_gui.bat
```

`run_ncbi_blast_gui.bat` remains as a legacy alias, but this GUI performs a
local transcript scan and does not submit oligo sequences to BLAST.

The GUI asks for the user's NCBI contact email on first use and saves it only on
that computer. The opening screen allows the saved email to be changed later.

The first dialog offers two workflows. **Single sequence and one transcript**
accepts one AS or SS sequence, one exact-version RefSeq accession such as
`NM_002439.5`, and these scan-region presets:

- Full sequence (selected by default)
- Seed, positions 2-8
- Core, positions 2-18

The single-sequence workflow shows five closest transcript windows per selected
region by default. These windows are not filtered by the maximum-mismatch
setting, and different region lengths are ranked separately. Results appear
directly in a scrollable text window with **Copy all** and **New scan** buttons;
this workflow does not create an Excel workbook.

**Excel sequence table** retains the batch workflow. It lets you choose AS or SS
sequence/name columns and a transcript source: a `target_accession` column, one
RefSeq accession for all rows, a local transcript FASTA/text file, or a private
versioned-accession panel. Single-target results are saved beside the input as
`<input filename>_ncbi_transcript_scan_results.xlsx`; panel results use
`<input filename>_private_transcript_panel_results.xlsx`.

Excel-table GUI workbooks begin with `input_queries`, then the compact
`comparison_results` sheet, followed by the unchanged `local_transcript_scan`
technical-detail sheet. `comparison_results` contains one best summary for each
query, target, and selected region, including the region start/end, result
status, qualifying-site count, best transcript coordinates, the two sequences
in query orientation, and human-readable differences. Mismatch positions are
1-based coordinates in the complete entered sequence: for example, the first
base of `seed:2-8` is reported as position `2`, not position `1`. A query with no
qualifying match remains in the summary with its closest transcript window.

GUI workflows share `TranscriptScanData\transcript_cache` beside the source or
packaged app. Each exact accession version is downloaded once and reused on
later runs. Single-sequence mode automatically downloads the transcript when it
is not already cached and can explicitly refresh it when requested. Cache-only
offline mode remains available for Excel private-panel scans. The cache contains
public transcript FASTA records only; oligo sequences remain local and are not
written into the cache. `TranscriptScanData\settings.json` stores the contact
email, and `TranscriptScanData\logs` contains rotating diagnostic logs.

#### Build a portable Windows app

The one-folder distribution includes Python and the required libraries, so the
recipient does not need to install Python. Build it on 64-bit Windows from the
repository root:

```powershell
python -m pip install -r deployment\requirements-transcript-scan.txt
powershell -NoProfile -ExecutionPolicy Bypass -File deployment\build_transcript_scan.ps1
```

Share the complete `dist\TranscriptScan` folder, or ZIP that folder without
moving `TranscriptScan.exe` out of it. The app creates its writable
`TranscriptScanData` subfolder on first use. Build and end-user instructions are
in `deployment\README_TRANSCRIPT_SCAN.txt`.

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-sequence "AUGCUACGGAUCUAGCUAGCU" --target-accession NM_000000.0
```

For an SS/sense sequence that is already in transcript orientation, use
`--ss-sequence`. The local transcript scan compares SS directly against the
transcript instead of reverse-complementing it:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --ss-sequence "UAGCUAGCUAGAUCCGUAGCA" --target-accession NM_000000.0
```

For a single AS sequence and target accession, omitting `--output` prints a
readable quick-look summary directly in the terminal. Add `--output` when you
want to save the same local scan results as CSV. If no local matches are found
within the default `--max-mismatches 3`, the quick-look output automatically
shows the 10 closest transcript windows so you can still inspect the nearest
sites:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-sequence "AUGCUACGGAUCUAGCUAGCU" --target-accession NM_000000.0 --output transcript_scan.csv
```

You can also compare against a local FASTA/plain transcript:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-sequence "AUGCUACGGAUCUAGCUAGCU" --target-file transcript.fasta --max-mismatches 3
```

The local scanner accepts exactly one transcript record per target file. It
rejects multi-record FASTA files instead of concatenating their records.

The local scan output includes both the transcript window in transcript
orientation and `transcript_match_as_5to3`, which is reverse-complemented back to
AS orientation so it can be compared directly with your AS sequence. If you want
both a CSV file and the terminal summary, add `--terminal`; if you prefer raw CSV
printed to the terminal, add `--stdout-csv`.

To choose a different number of closest windows, add `--closest N` anywhere in
the command. The normal scan still uses `--max-mismatches`, and the terminal
also shows the closest transcript windows without applying that cutoff:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-sequence "AUGCUACGGAUCUAGCUAGCU" --target-accession NM_000000.0 --closest 10
```

For oligo risk review, you can scan full AS plus custom subregions:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-table as_sequences.xlsx --as-column antisense --as-name-column oligo_id --target-file transcript.fasta --scan-region full --scan-region seed:2-8 --scan-region core:2-18 --result-workbook as_review.xlsx
```

For a privacy-preserving local panel scan, repeat an exact-version RefSeq
accession. NCBI EFetch receives only these public accessions; guide sequences
remain on the local computer:

```powershell
$guide = Read-Host "Enter AS sequence"
python -m tools_for_pharma.oligo.ncbi_blast --as-sequence $guide --private-panel --target-accession NM_000041.4 --target-accession NM_001302688.2 --cache-dir .ncbi_transcript_cache --max-mismatches 0
```

Targets can also come from a text, CSV, or Excel table. The default accession
column names are `target_accession`, `accession`, or `refseq`:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-table guides.xlsx --as-column antisense --as-name-column oligo_id --private-panel --target-table transcript_targets.xlsx --target-column target_accession --result-workbook private_panel_results.xlsx
```

Use `--download-targets-only` to retrieve and validate the public references
before entering any guide. After the one-record FASTA files are cached, add
`--offline` to prohibit all NCBI requests during scanning. Private panel
CLI private-panel workbooks contain `input_queries`, `transcript_targets`,
`local_transcript_scan`, `query_target_summary`, and `run_metadata`. The GUI
uses its compact `comparison_results` sheet instead of the redundant
`query_target_summary`. Every guide-target pair receives a match, no-match, or
target-error result.

The transcript-scan GUI also supports AS or SS tables, private accession panels,
target tables, and offline cache-only scanning.

For broader NCBI BLAST URL API searches, use `--blast` or `--blast-only`:

**Privacy warning:** these options transmit the input oligo sequence(s) to the
remote NCBI BLAST service. Local transcript scans do not transmit the oligo.

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-sequence "AUGCUACGGAUCUAGCUAGCU" --blast-only --database refseq_rna --blast-output blast_hits.csv
```

Batch BLAST can read multiple AS sequences from FASTA/plain text or from an
Excel/CSV table:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-file as_sequences.fasta --blast-only --database refseq_rna --blast-output blast_hits.csv
python -m tools_for_pharma.oligo.ncbi_blast --as-table as_sequences.xlsx --as-column antisense --as-name-column oligo_id --blast-only --database refseq_rna --blast-output blast_hits.csv
```

The same BLAST-only workflow accepts SS/sense inputs with `--ss-file`,
`--ss-table`, `--ss-column`, and `--ss-name-column`.

Short oligo queries are submitted as multi-FASTA batches instead of one BLAST job
per sequence. The default batch cap is 1,000 total bases per BLAST request.
The output CSV includes the BLAST RID and the query ID for each hit.

For batch work, the preferred output is an Excel result workbook:

```powershell
python -m tools_for_pharma.oligo.ncbi_blast --as-table as_sequences.xlsx --as-column antisense --as-name-column oligo_id --blast-only --database refseq_rna --result-workbook as_blast_results.xlsx
```

Remote BLAST workbooks contain `input_queries`, `local_transcript_scan`,
`blast_hits_raw`, `blast_hits_filtered`, `blast_batches`, and `run_metadata`.
Local-only workbooks omit the BLAST sheets. Without explicit CSV or workbook
paths, local scans write `<input>_ncbi_transcript_scan_results.xlsx` and remote
BLAST runs write `<input>_ncbi_blast_results.xlsx`. Use `--cache-dir` to reuse
fetched NM/XM transcript FASTA files across runs.

NCBI asks API users to include `tool` and `email`, avoid contacting BLAST more
than once every 10 seconds, and avoid polling a single RID more than once per
minute. The tool uses safer defaults: at least 15 seconds between NCBI requests
and at least 75 seconds between status checks for the same RID. The GUI asks for
and locally saves the current user's email. CLI commands that need NCBI network
access require `--email user@example.com`; cached/offline local scans do not.

## qPCR Table Extraction And Plotting

These tools turn qPCR Excel report tables into plot-ready data, then make bar
plots from the reviewed extracted sheet.

Use the extraction GUI:

```text
run_qpcr_extract_gui.bat
```

After choosing a workbook and worksheet, the GUI shows a column-mapping review.
It auto-selects Group, Compound ID, Sample ID, Animal ID, the individual plotted
value, bar mean, SEM, and optional sample size. Two-row headers are displayed as
combined labels such as `APOE | Mean CT`, so the final group `Mean` is not
confused with a gene-level `Mean CT`. Adjust any dropdown before continuing.
Known schemas remain automatic, including both `MEAN RQ` with `Relative to
control group` and the shorter `Mean` with `Normalized RQ`. If no sample-size
column exists, it is inferred from the number of extracted individual values.
Before writing the `plotdata-...` sheet, a preview reports the number of bars
and individual values and warns when a provided mean differs materially from
the mean calculated from its individual values.

Or run extraction from PowerShell:

```powershell
python -m tools_for_pharma.qpcr.extract -i "data/examples/qpcr/BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet qPCR
```

Use the plotting GUI:

```text
run_qpcr_plot_gui.bat
```

Or run plotting from PowerShell:

```powershell
python -m tools_for_pharma.qpcr.plot -i "data/examples/qpcr/BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet "plotdata-qPCR" --plot both
```

Plot modes:

- `split`: one plot per reference source
- `combined`: one grouped plot with all reference sources together
- `both`: create both styles

New extractions keep one row per animal in the `plotdata-...` sheet. Along with
the reviewed group `MEAN RQ` and `SEM`, the sheet includes `Sample ID`,
`Animal ID`, and `Individual RQ`. The individual value is taken from `Relative
to control group` for each reference-gene table and from `Geomean` for the
aggregate table. The qPCR plotter overlays these values as jittered dots on the
mean +/- SEM bars. Older extracted sheets without `Individual RQ` remain
compatible and are plotted as bars without dots.

By default, plots are saved beside the Excel file in a subfolder based on the
workbook name. Existing generated examples have been moved to `outputs/plots/`.

Generic grouped-bar plotting now lives in `tools_for_pharma.plotting.bar`, so it
can be used outside qPCR workflows. The old qPCR command still works as a
compatibility entrypoint.

For a simple two-column table where the first column is `Group` and the second
column contains values like `0.72 +/- 0.13` or `0.72 ± 0.13`, make a grouped bar
plot directly:

```text
run_simple_group_plot_gui.bat
```

Or run it from PowerShell:

```powershell
python -m tools_for_pharma.qpcr.simple_group_plot -i "group_plot.xlsx" --title "MSH3 remaining on D33 relative to baseline in Liver" --y-label "Remaining relative to baseline"
python -m tools_for_pharma.plotting.bar -i "group_plot.xlsx" --title "MSH3 remaining on D33 relative to baseline in Liver" --y-label "Remaining relative to baseline"
```

Labels such as `G1-baseline`, `G1-2mpk D33`, and `G1-5mpk D33` are grouped under
`G1`; the text after the hyphen becomes the bar label in the legend.

The same tool also supports wider tables with `Dose`, `Group`, and multiple
`Time-...` columns. For example, with columns such as `Dose (mpk)`, `Group`,
`Time-baseline`, `Time-D8`, and `Time-D29`, the default mode creates:

- one plot with timepoints on the x-axis and compound+dose bars
- one plot per compound comparing doses across time
- one plot per dose comparing compounds across time

To create only the all-variable plot:

```powershell
python -m tools_for_pharma.qpcr.simple_group_plot -i "group_plot.xlsx" --plot-mode all-variables
```

## Generic Curve Interpolation Plotting

Use the curve plotter for dose-response or inhibition-rate tables with a
positive concentration/dose column and an inhibition/response column. By
default, it fits a smooth 4-parameter logistic curve, the common sigmoidal
dose-response model used for IC50-style analysis in biopharma, and marks
`IC50`, `IC75`, and `IC90`. The plot also includes a curve summary with response
range, log-dose AUC, Hill slope, R-squared, and the local slope at each marked IC
value.

Input data should usually be a wide table: the first column is concentration,
and every following column is one compound's inhibition rate.

```csv
Concentration (nM),AD-001,AD-002
0.1,5,3
1,28,18
10,61,49
100,92,88
```

Values can also include an inline error value. The mean is used for the curve,
and the value after `+/-` or the plus-minus symbol is drawn as an error bar:

```csv
Concentration (nM),AD-001,AD-002
0.1,5.00+/-0.50,3.00+/-0.20
1,28.0+/-1.2,18.0+/-0.9
10,61.0+/-3.4,49.0+/-2.1
100,92.0+/-2.8,88.0+/-3.0
```

If the first row only contains numbers, the tool treats it as data instead of a
header row. The first column is still concentration, and later columns are named
`Compound 1`, `Compound 2`, and so on.

```csv
0.1,5,3
1,28,18
10,61,49
100,92,88
```

```text
run_IC50_plot_gui.bat
run_curve_plot_gui.bat
```

Or run it from PowerShell:

```powershell
python -m tools_for_pharma.plotting.IC50 -i "inhibition_curve.xlsx" --sheet "Sheet1"
```

By default, the title comes from the input file name, output files are saved
beside the input file, and the marked values are `IC50`, `IC75`, and `IC90`.
Choose different IC markers with `--ic`:

```powershell
python -m tools_for_pharma.plotting.IC50 -i "inhibition_curve.csv" --ic "IC50,IC80,IC90"
```

For many compounds, use one image per compound so the IC summary remains
readable:

```powershell
python -m tools_for_pharma.plotting.IC50 -i "inhibition_curve.csv" --plot-mode single
```

To create both the combined overview and one detailed image per compound:

```powershell
python -m tools_for_pharma.plotting.IC50 -i "inhibition_curve.csv" --plot-mode both
```

The IC50-specific summary panel is added by `tools_for_pharma.plotting.IC50`.
The lower-level `tools_for_pharma.plotting.curve` module stays summary-free for
general curve plotting. When a combined IC50 plot has more than three curves,
the detailed summary panel is omitted from that combined image; use
`--plot-mode single` or `--plot-mode both` for the per-compound summaries.

To use the older point-to-point log-dose interpolation instead of 4PL fitting:

```powershell
python -m tools_for_pharma.plotting.IC50 -i "inhibition_curve.csv" --fit-method interpolation
```

## qPCR Reference-Gene QC

This exploratory QC workflow checks whether reference genes look stable across
groups before relying on them for normalization. It is separate from the main
`plotdata-...` MEAN RQ workflow.

Use the reference-QC extraction GUI:

```text
run_qpcr_ref_qc_extract_gui.bat
```

Or run extraction from PowerShell:

```powershell
python -m tools_for_pharma.qpcr.ref_qc_extract -i "data/examples/qpcr/BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet qPCR
```

Use the reference-QC plotting GUI:

```text
run_qpcr_ref_qc_plot_gui.bat
```

Or run plotting from PowerShell:

```powershell
python -m tools_for_pharma.qpcr.ref_qc_plot -i "data/examples/qpcr/BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet refqc-qPCR
```

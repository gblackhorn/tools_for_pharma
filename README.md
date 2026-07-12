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

Use the GUI for an Excel AS table:

```text
run_ncbi_blast_gui.bat
```

The GUI lets you choose the input workbook, AS sequence/name columns, and the
transcript source: use a `target_accession` column, type one RefSeq accession for
all rows, or choose a local transcript FASTA/text file. Results are saved beside
the input workbook as `<input filename>_ncbi_blast_results.xlsx`.

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

For broader NCBI BLAST URL API searches, use `--blast` or `--blast-only`:

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

The workbook contains `input_queries`, `local_transcript_scan`,
`blast_hits_raw`, `blast_hits_filtered`, `blast_batches`, and `run_metadata`.
If you use `--as-file` or `--as-table` without CSV output paths, the tool writes
`<input>_ncbi_blast_results.xlsx` by default. Use `--cache-dir` to reuse fetched
NM/XM transcript FASTA files across runs.

NCBI asks API users to include `tool` and `email`, avoid contacting BLAST more
than once every 10 seconds, and avoid polling a single RID more than once per
minute. The tool uses safer defaults: at least 15 seconds between NCBI requests
and at least 75 seconds between status checks for the same RID. The default
contact email is `da.guo@argobiopharma.com`; pass `--email` to override it.

## qPCR Table Extraction And Plotting

These tools turn qPCR Excel report tables into plot-ready data, then make bar
plots from the reviewed extracted sheet.

Use the extraction GUI:

```text
run_qpcr_extract_gui.bat
```

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

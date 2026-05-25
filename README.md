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

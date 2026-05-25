# tools_for_pharma

Simple utilities for pharma-oriented oligo sequence work.

## Utility Files

Use these two user-facing scripts first:

- `util_oligo_offtarget.py`: prepare oligo antisense region information for
  in silico off-target analysis.
- `util_transcript_sequence.py`: extract SS/AS oligo design sequences from long
  NCBI-style transcript sequences.

`oligo_utils.py` is a combined entrypoint for all workflows. The `oligo_*`
and `util_...` names are the consistent names to use.

## Off-Target Analysis Preparation

This workflow starts from an oligo antisense sequence. It extracts the selected
antisense region, usually positions **2-18**, and computes the complementary
sense sequence in **5'->3'** orientation. For the default region, the output
columns are:

- `normalized_antisense`: cleaned full input sequence, uppercase RNA with `T`
  converted to `U`
- `antisense_2-18_5to3`: antisense positions 2 through 18, still written 5'->3'
- `sense_18-2_5to3`: reverse-complement sense sequence aligned to antisense
  positions 18 through 2, written 5'->3'
- `error`: row-level validation message, blank when processing succeeded

The batch table workflow appends these columns to the end of the original table
by default. It does not create a separate output file unless you provide
`--output`.

Process one antisense sequence:

```powershell
python util_oligo_offtarget.py "AUGCUACGGAUCUAGCUAGCU"
```

Use a non-default antisense region:

```powershell
python util_oligo_offtarget.py "AUGCUACGGAUCUAGCUAGCU" --start 3 --end 19
```

Process many antisense sequences from an Excel/CSV table:

```powershell
python util_oligo_offtarget.py --input oligos.xlsx --column antisense
```

Save to a separate file only when needed:

```powershell
python util_oligo_offtarget.py --input oligos.xlsx --column antisense --output oligos_checked.xlsx
```

Open the off-target table GUI:

```powershell
python util_oligo_offtarget.py --gui
```

Or double-click:

```text
run_util_oligo_offtarget_gui.bat
```

## Transcript Sequence Extraction

This workflow starts from a long FASTA/plain-text transcript, such as a sequence
copied from NCBI. It extracts a 1-based inclusive transcript range and outputs
the matched **SS** and **AS** strands for oligo design.

Extract one transcript range:

```powershell
python util_transcript_sequence.py --transcript-file transcript.fasta --start 120 --end 140
```

For multiple ranges, provide an Excel/CSV table with `start` and `end` columns.
Column matching is case-insensitive. The output table fills `SS matched`,
`AS matched`, `matched length (nt)`, and `match error`. Rows with blank `start`
or `end` are skipped without an error:

```powershell
python util_transcript_sequence.py --transcript-file transcript.fasta --range-table ranges.xlsx
```

Open the transcript sequence GUI:

```powershell
python util_transcript_sequence.py --gui
```

Or double-click:

```text
run_util_transcript_sequence_gui.bat
```

For the combined GUI, double-click:

```text
run_oligo_gui.bat
```

## qPCR Table Extraction And Plotting

These scripts turn qPCR Excel report tables into plot-ready data, then make
bar plots from the reviewed extracted sheet.

`qPCR_extract_excel_table.py` reads a qPCR worksheet, including tables embedded
inside a larger report sheet. It fills merged-cell values, reads saved Excel
formula results as values, skips NRT/NTC/QC rows, and extracts:

- `Group`
- `Compound ID`
- `Reference source`
- `Sample size (n)`
- `MEAN RQ`
- `SEM`

By default, extraction appends the output to the same workbook in a new sheet
named `plotdata-[original sheet name]`. For example, extracting from sheet
`qPCR` creates `plotdata-qPCR`. Close the workbook in Excel before extraction,
because Excel may block writing the appended sheet.

Use the extraction GUI:

```powershell
python qPCR_extract_excel_table.py --gui
```

The GUI lets you choose the Excel file, choose the source worksheet, and confirm
the new `plotdata-...` sheet name.

Or double-click:

```text
run_qpcr_extract_gui.bat
```

Use extraction from the command line:

```powershell
python qPCR_extract_excel_table.py -i "BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet qPCR
```

`qPCR_plot_excel_table.py` reads a reviewed `plotdata-...` sheet and creates
plots. The selected plot sheet must start with `plotdata-` and contain the
required qPCR plotting columns. If `--sheet` is omitted, the plot script uses
the last worksheet whose name starts with `plotdata-`.

Plot modes:

- `split`: one plot per reference source
- `combined`: one grouped plot with all reference sources together
- `both`: create both styles

By default, plots are saved beside the Excel file in a subfolder with the same
name as the workbook stem. For
`BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx`, the
default plot folder is:

```text
BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20
```

The folder is created automatically if it does not exist.
Each plot is saved in both `.png` and `.svg` format.

Use the plotting GUI:

```powershell
python qPCR_plot_excel_table.py --gui
```

The GUI lets you choose the Excel file, choose the `plotdata-...` sheet, and
choose the plot mode. It saves plots to the default folder described above.

Or double-click:

```text
run_qpcr_plot_gui.bat
```

Use plotting from the command line:

```powershell
python qPCR_plot_excel_table.py -i "BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet "plotdata-qPCR" --plot both
```

You can omit `--sheet` when the workbook already has the correct `plotdata-...`
sheet:

```powershell
python qPCR_plot_excel_table.py -i "BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --plot both
```

## qPCR Reference-Gene QC

This exploratory QC workflow checks whether reference genes look stable across
groups before relying on them for normalization. It is separate from the main
`plotdata-...` MEAN RQ workflow.

`qpcr_ref_qc_extract.py` extracts each reference gene's sample-level `Mean CT`,
then calculates group-level mean CT and SEM. The output is appended to the same
workbook in a sheet named `refqc-[original sheet name]`, such as `refqc-qPCR`.
It can follow simple Excel formulas used by these qPCR reports, including direct
cell links and `AVERAGE(...)`, so it can recover CT values even when cached
formula results are missing.

Use the reference-QC extraction GUI:

```text
run_qpcr_ref_qc_extract_gui.bat
```

Or from the command line:

```powershell
python qpcr_ref_qc_extract.py -i "BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet qPCR
```

`qpcr_ref_qc_plot.py` reads a `refqc-...` sheet and creates a Mean CT QC plot
with one point and SEM error bar per group/reference gene. The plot is saved as
both `.png` and `.svg` in the default plot folder beside the workbook.

Use the reference-QC plotting GUI:

```text
run_qpcr_ref_qc_plot_gui.bat
```

Or from the command line:

```powershell
python qpcr_ref_qc_plot.py -i "BWS-2a ICV #10-Brain-HTT1a-3内参geomean-qPCR report-2026-05-20.xlsx" --sheet refqc-qPCR
```

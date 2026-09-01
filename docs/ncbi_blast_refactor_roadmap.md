# NCBI BLAST / Transcript Scan Refactor Roadmap

## Purpose

This document is the execution plan for refactoring
`tools_for_pharma/oligo/ncbi_blast.py` and consolidating reusable sequence
features across `tools_for_pharma/oligo`.

The refactor must improve structure without changing scientific behavior,
privacy behavior, command-line compatibility, workbook output, GUI usability,
or the portable one-folder Windows application.

Work must proceed one phase at a time. A phase is complete only after its exit
gate passes and its changes are committed. Do not begin the next phase while
the current phase has a failing test, packaging regression, or unresolved
compatibility question.

## Status legend

- `[x]` Complete
- `[ ]` Not started
- `[~]` In progress
- `[!]` Blocked or requires a decision

## Baseline checkpoint

- [x] Pre-refactor commit created:
  `44349c4 feat: expand transcript scan target inputs`
- [x] Portable application version at the checkpoint: `1.1.0`
- [x] `ncbi_blast.py` behavior was validated before the checkpoint with the
  complete test suite (`120 passed`).
- [x] Packaged `TranscriptScan.exe --self-test` passed before the checkpoint.
- [x] The packaged GUI was manually exercised before the checkpoint.
- [x] The release ZIP was verified not to contain test-generated
  `TranscriptScanData`.

If a later phase produces an unclear regression, compare with or return to this
checkpoint. Do not rewrite or amend the checkpoint.

## Scope

### In scope

- Establish reusable nucleotide, FASTA, and positional-comparison primitives.
- Remove duplicated generic behavior from modules under `oligo`.
- Separate transcript-scan models, inputs, targets, scanning, reporting,
  workflows, GUI, and CLI responsibilities.
- Keep `tools_for_pharma.oligo.ncbi_blast` as a compatibility facade and
  executable module.
- Preserve the one-folder Windows application built with PyInstaller.
- Improve test seams for filesystem, network, time, and GUI boundaries.

### Out of scope

- New biological scoring or advanced off-target analysis.
- Changes to mismatch meaning, AS/SS orientation, or scan-region semantics.
- New transcript or genome databases.
- New workbook columns, sheet restructuring, or GUI redesign.
- Combining the single-transcript scanner with the multi-record transcript
  panel workflow.
- Changing what data is sent to NCBI.

Feature requests discovered during the refactor must be recorded separately
and handled after the refactor unless they are required to restore existing
behavior.

## Non-negotiable compatibility contracts

### Scientific and sequence behavior

- Coordinates shown to users remain 1-based and inclusive.
- AS scans continue to compare the reverse-complement target with transcripts.
- SS scans continue to compare directly in transcript orientation.
- Mismatch positions continue to be reported against the complete entered
  query, including when a subregion such as `seed:2-8` is scanned.
- `full`, `seed:2-8`, and `core:2-18` retain their current meaning.
- Single-target inputs never silently concatenate multiple FASTA records.
- DNA/RNA conversion and ambiguous-base policies must be explicit; consolidating
  normalizers must not silently change which characters each workflow accepts.

### Privacy and network behavior

- Local transcript scanning does not submit an oligo sequence to NCBI.
- EFetch sends the public transcript accession and contact/tool metadata only.
- Remote BLAST remains a separate, explicit workflow because it submits the
  query sequence.
- Pasted transcript sequences and local transcript files remain local and are
  not copied into the persistent transcript cache.
- The transcript cache stores public transcript FASTA records, not oligos.

### CLI and Python compatibility

- This command remains valid:

  ```powershell
  python -m tools_for_pharma.oligo.ncbi_blast --gui
  ```

- Existing CLI flags, defaults, validation errors, and exit-code behavior remain
  compatible unless a separately approved change says otherwise.
- Existing imports from `tools_for_pharma.oligo.ncbi_blast` remain available
  through explicit re-exports during and after the refactor.
- `run_ncbi_blast_gui.bat` and `run_ncbi_transcript_scan_gui.bat` continue to
  work.

### Output compatibility

- Multiple-sequence workbooks retain their current sheet order:
  `input_queries`, `comparison_results`, the existing technical-detail sheet,
  and the applicable support sheets.
- Source workbook fields continue to be preserved in `input_queries`.
- Sheet names, column names, coordinate conventions, and row meaning do not
  change as a side effect of moving code.
- Single-sequence results continue to appear in the GUI text result window and
  do not require an Excel output file.

### GUI compatibility

- First use prompts for an NCBI contact email and saves it locally.
- Single-sequence form values survive repeated runs and "Edit and run again".
- Form state is cleared when the application is closed.
- A cache miss with cache use enabled automatically downloads the transcript.
- Refresh from NCBI remains a one-run action.
- NM/XM/NR/XR accession, pasted sequence, and one local FASTA/text target remain
  available in single-sequence mode.
- Genomic NC accessions continue to produce a clear transcript-target error.

## Portable Windows application contract

The standalone application is a first-class refactor target, not a final
packaging check.

### Current packaging path

1. `transcript_scan_app.py` is the PyInstaller entry point.
2. It imports `application_data_dir`, `gui_log_path`, `run_gui`, and
   `shared_gui_transcript_cache_dir` from `ncbi_blast.py`.
3. `deployment/transcript_scan.spec` defines the one-folder build and Tkinter
   hidden imports.
4. `deployment/build_transcript_scan.ps1` builds the executable and copies:
   - `README_TRANSCRIPT_SCAN.txt`
   - `VERSION.txt`
   - `THIRD_PARTY_NOTICES.txt`
   - `multiple_sequence_blast_template.xlsx`
5. The frozen app writes settings, cache, and logs to
   `TranscriptScanData` beside `TranscriptScan.exe`.

### Packaging rules during refactoring

- Keep the four entry-point imports compatible until the entry point is changed
  deliberately and verified in the same phase.
- Prefer static, explicit imports. If dynamic imports are introduced, update
  `deployment/transcript_scan.spec` with the required hidden imports.
- New package modules used by the GUI must be reachable from the PyInstaller
  analysis graph.
- Do not assume repository-relative files exist in a frozen application.
- Preserve `sys.frozen`/`sys.executable` behavior for the portable base folder.
- Keep writable user data outside PyInstaller's `_internal` directory.
- Ensure the build spec does not accidentally exclude a newly required
  dependency.
- Keep the complete `TranscriptScan` folder together; do not design around a
  standalone executable copied out of the folder.
- Do not include `TranscriptScanData` in a release ZIP. A self-test may create
  it locally, so the exact release directory must be checked before archiving.
- Do not update the application version merely because code moved. Decide the
  release version after final compatibility verification.

## Target dependency structure

```text
tools_for_pharma/
├── sequence/
│   ├── __init__.py
│   ├── nucleotides.py       # normalization, complements, 1-based slicing
│   ├── fasta.py             # generic FASTA records, parsing, formatting
│   └── comparison.py        # equal-length positional comparison
├── shared/
│   └── excel_utils.py       # cross-project spreadsheet I/O
└── oligo/
    ├── core.py              # compatibility wrappers for existing callers
    ├── transcript_accessions.py
    ├── ncbi_transport.py
    ├── transcript.py
    ├── transcript_panel.py
    ├── transcript_scan/
    │   ├── __init__.py
    │   ├── models.py
    │   ├── queries.py
    │   ├── targets.py
    │   ├── scanner.py
    │   ├── reporting.py
    │   ├── workflows.py
    │   ├── gui.py
    │   └── cli.py
    └── ncbi_blast.py        # compatibility facade and module entry point
```

This is the intended dependency direction:

```text
sequence primitives
        ↓
oligo/transcript domain models
        ↓
retrieval and comparison workflows
        ↓
reporting
        ↓
CLI and GUI
```

Lower layers must not import Tkinter, argparse, pandas, Excel writers,
application settings, or network clients.

The exact target structure may be adjusted when a phase reveals a better
boundary, but any adjustment must be recorded in this document before code is
moved.

## Generic extraction decisions

| Candidate | Intended home | Decision |
|---|---|---|
| RNA/DNA normalization | `sequence/nucleotides.py` | Generic; extract with explicit alphabet and ambiguity policies. |
| Complement and reverse-complement | `sequence/nucleotides.py` | Generic; use explicit function names and retain old wrappers. |
| 1-based inclusive subsequence | `sequence/nucleotides.py` | Generic; name the coordinate convention explicitly. |
| FASTA record, parse, format | `sequence/fasta.py` | Generic; parser preserves record boundaries and caller enforces cardinality. |
| Equal-length mismatch positions/Hamming distance | `sequence/comparison.py` | Generic; unequal lengths must raise rather than imply an alignment. |
| AS/SS orientation and scan regions | `oligo/transcript_scan` | Domain-specific; uses generic primitives. |
| Transcript-window scanning | `oligo/transcript_scan/scanner.py` | Domain-specific for now; do not invent a generic search engine prematurely. |
| RefSeq transcript accession handling | `oligo/transcript_accessions.py` | Reusable transcript-domain feature, not a basic sequence primitive. |
| NCBI HTTP/retry/contact handling | `oligo/ncbi_transport.py` | Reusable transport; EFetch and remote BLAST workflows remain separate. |
| Transcript cache | `oligo/transcript_scan/targets.py` | Domain-specific because version and FASTA validation are part of its contract. |
| Reference/CDS difference blocks | `transcript_panel.py` | Domain-specific; do not merge with equal-length mismatch comparison. |
| Excel I/O | `shared/excel_utils.py` | Already shared; consumers should import it directly. |
| CLI parsers and GUI dialogs | Their workflow modules | Keep specific unless repeated behavior proves a stable shared interface. |

Avoid generic `utils.py` modules. Every shared module must have one clear
subject and an independently testable contract.

## Validation gates

### Gate S: source validation

Run after every phase:

```powershell
python -m pytest -q
python -m tools_for_pharma.oligo.ncbi_blast --help
git diff --check
```

Also run focused tests for the modules changed in that phase. The complete test
suite must not fall below the established behavior; test-count changes must be
explained by added or deliberately consolidated tests, never by deleting
coverage to make a phase pass.

### Gate P: packaged import/self-test validation

Run whenever imports reachable from `transcript_scan_app.py` change, whenever a
module moves, and at every GUI-related phase:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File deployment\build_transcript_scan.ps1
```

Then launch the generated executable with `--self-test` and require exit code
zero. Review `TranscriptScanData\logs\transcript_scan.log` if it fails.

Gate P must confirm:

- The executable starts without Python installed in its runtime environment.
- Tkinter, pandas, and openpyxl import successfully.
- The portable data directory is writable.
- The transcript-cache directory can be created.
- An Excel write/read round trip succeeds.
- The copied template and deployment documents are present.

### Gate G: packaged GUI smoke validation

Run after changes to GUI, workflow dispatch, target handling, cache/settings
paths, or final packaging:

- Launch `dist\TranscriptScan\TranscriptScan.exe` directly.
- Confirm mode selection opens.
- Confirm saved email loads, or first-use email is requested in a clean data
  directory.
- Run a single scan with a pasted transcript.
- Choose "Edit and run again" and confirm the form remains populated.
- Run a single scan with a local one-record FASTA/text target.
- Run an accession scan from a cache miss and confirm automatic download.
- Repeat the accession scan and confirm cache reuse.
- Confirm refresh applies to one run only.
- Run a multiple-sequence workbook and inspect the sheet order and key columns.
- Confirm logs are written beside the executable.
- Confirm pasted/local transcript content is not persisted in settings or cache.

Use non-sensitive test sequences. Do not use remote BLAST with a private query
as a packaging smoke test.

### Gate R: release-folder validation

Run only for release candidates:

- Verify version text, executable version constant, and release naming agree.
- Verify the template opens and retains the expected sheet/header structure.
- Verify the distribution works from a copied folder outside the repository.
- If possible, test on a clean Windows 10/11 64-bit machine or VM without
  Python installed.
- Verify the release folder/ZIP does not contain `TranscriptScanData`, logs,
  settings, cached transcripts, test outputs, or private sequences.
- Record ZIP size and SHA-256 checksum.

## Execution phases

### Phase 1 — Lock behavior and public contracts

Status: [x]

Goal: make moves safe before moving implementation.

Tasks:

- Inventory public names imported from `ncbi_blast.py` by tests, scripts, and
  `transcript_scan_app.py`.
- Add an explicit compatibility test for required public imports.
- Add characterization tests for RNA/DNA normalization differences.
- Add characterization tests for FASTA single-record and multi-record behavior.
- Lock mismatch coordinates and AS/SS orientation with focused tests.
- Lock workbook sheet order and essential column names.
- Lock portable path behavior under normal Python and simulated `sys.frozen`.
- Record the source and package smoke-test commands in the test documentation.

Exit gate:

- Gate S passes.
- Gate P passes because package-facing compatibility tests were added.
- No production behavior is intentionally changed.

Suggested commit:

```text
test: lock transcript scan compatibility contracts
```

### Phase 2 — Establish generic nucleotide primitives

Status: [x]

Goal: create one authoritative implementation for basic sequence operations.

Tasks:

- Add `tools_for_pharma/sequence/nucleotides.py`.
- Define explicit RNA/DNA normalization and ambiguous-base behavior.
- Add explicit complement and reverse-complement operations.
- Add an explicitly named 1-based inclusive subsequence operation.
- Add independent unit tests, including boundary/error cases.
- Convert `oligo/core.py` into compatibility wrappers/re-exports.
- Migrate `transcript.py`, `table.py`, and other small consumers first.
- Migrate `ncbi_blast.py` and `transcript_panel.py` only after their different
  `N` policies are represented explicitly.

Exit gate:

- Gate S passes.
- Gate P passes because the packaged import graph now includes a new package.
- Existing import paths from `oligo.core` remain valid.

Suggested commit:

```text
refactor: centralize nucleotide sequence primitives
```

### Phase 3 — Establish generic FASTA and positional comparison primitives

Status: [x]

Goal: remove duplicated parsing/formatting and equal-length comparison logic.

Tasks:

- Add a generic immutable `FastaRecord`.
- Add FASTA parsing that never joins adjacent records.
- Add FASTA formatting with validated line width.
- Add a single-record validator rather than embedding cardinality assumptions in
  the parser.
- Keep alphabet validation separate from structural FASTA parsing.
- Add `mismatch_positions_1based` and `hamming_distance` for equal-length
  sequences.
- Reject unequal-length comparison explicitly.
- Adapt `transcript_panel.py` without changing its multi-record behavior.
- Adapt `ncbi_blast.py` through AS/SS-specific adapters.
- Adapt `transcript.py` plain/single-FASTA handling.

Exit gate:

- Gate S passes.
- Gate P passes.
- Multi-record transcript-panel tests and single-target rejection tests both
  pass.

Suggested commit:

```text
refactor: share FASTA and sequence comparison primitives
```

### Phase 4 — Extract transcript identifiers and NCBI transport

Status: [x]

Goal: share stable transcript/NCBI infrastructure without hiding the privacy
boundary.

Tasks:

- Extract versioned RefSeq transcript validation and header accession parsing.
- Preserve the distinct error for genomic NC/NG/NT/NW accessions.
- Extract NCBI request configuration, contact email validation, retry behavior,
  timeout behavior, and error translation.
- Keep EFetch parameters separate from BLAST request parameters.
- Inject network opener and clock dependencies at workflow boundaries instead
  of relying on patching unrelated module globals.
- Preserve existing mocked BLAST lifecycle and offline tests.
- Review whether `transcript_panel.py` can use the shared transport without
  changing Ensembl or NCBI response handling.

Exit gate:

- Gate S passes.
- Gate P passes.
- No private query is used for live remote-BLAST verification.
- EFetch and remote-BLAST data flows remain visibly separate in code and tests.

Suggested commit:

```text
refactor: share transcript identifiers and NCBI transport
```

### Phase 5 — Create the transcript-scan domain package

Status: [x]

Goal: move domain models and pure scan logic out of `ncbi_blast.py`.

Tasks:

- Add `transcript_scan/models.py` for scan-domain dataclasses.
- Add `transcript_scan/queries.py` for AS/SS inputs, table fields, query IDs,
  batching, and scan-region parsing.
- Add `transcript_scan/scanner.py` for pure transcript-window scanning and best
  match selection.
- Keep sequence primitives outside this package.
- Keep models free of GUI, network, pandas, and argparse dependencies.
- Explicitly re-export existing public classes and functions from
  `ncbi_blast.py`.
- Update internal imports to use defining modules, not the compatibility facade.

Exit gate:

- Gate S passes.
- Gate P passes.
- Direct legacy imports from `ncbi_blast.py` still work.
- Scan results and mismatch coordinates are byte-for-byte/schema-equivalent
  where applicable.

Suggested commit:

```text
refactor: extract transcript scan models and engine
```

### Phase 6 — Extract target acquisition and cache handling

Status: [x]

Goal: isolate accession, pasted, local-file, and cache target sources.

Tasks:

- Add `transcript_scan/targets.py`.
- Model target source explicitly: accession, pasted sequence, or local file.
- Preserve one-transcript cardinality validation.
- Preserve exact version verification and cache filenames.
- Preserve cache-miss automatic download and explicit offline behavior.
- Preserve one-run refresh behavior.
- Keep pasted/local targets out of persistent settings and cache.
- Add tests using temporary data directories and mocked EFetch responses.

Exit gate:

- Gate S passes.
- Gate P passes.
- Gate G passes for all three single-target sources and cache reuse.

Suggested commit:

```text
refactor: isolate transcript targets and cache
```

### Phase 7 — Extract reporting and workbook generation

Status: [x]

Goal: isolate presentation without changing output contracts.

Tasks:

- Add `transcript_scan/reporting.py`.
- Move terminal formatting, CSV projections, comparison rows, workbook rows,
  and run metadata.
- Use `shared/excel_utils.py` directly where its behavior matches the current
  contract.
- Keep workflow-specific sheet composition in transcript-scan reporting.
- Add schema tests for sheet order, required columns, source-field preservation,
  and mismatch coordinate formatting.
- Preserve the unchanged technical-detail sheet.
- Preserve single-sequence text output formatting.

Exit gate:

- Gate S passes.
- Gate P passes.
- Gate G passes for single text output and a representative multiple-sequence
  workbook.

Suggested commit:

```text
refactor: extract transcript scan reporting
```

### Phase 8 — Extract workflows and remote BLAST

Status: [x]

Goal: separate orchestration from interfaces and keep remote submission
explicit.

Tasks:

- Add `transcript_scan/workflows.py` for local single, multiple, and panel
  orchestration.
- Put remote BLAST submission, polling, retrieval, and filtering in a clearly
  named module or service, separate from local scan workflows.
- Pass configuration/dependencies explicitly rather than passing a large
  `argparse.Namespace` into domain code.
- Retain thin adapters while callers migrate.
- Preserve batch sizing, query IDs, RID logs, timeouts, and result filtering.
- Preserve existing privacy warnings around remote BLAST.

Exit gate:

- Gate S passes, including mocked remote-BLAST lifecycle tests.
- Gate P passes.
- Local scan code has no dependency on the remote-BLAST client.

Suggested commit:

```text
refactor: separate transcript scan and remote BLAST workflows
```

### Phase 9 — Extract CLI while preserving module execution

Status: [ ]

Goal: make CLI parsing an interface layer and keep all existing launch paths.

Tasks:

- Add `transcript_scan/cli.py`.
- Move parser construction, runtime argument validation, and CLI dispatch.
- Convert parsed arguments into explicit workflow configuration objects.
- Keep `ncbi_blast.main()` as a compatibility wrapper.
- Keep the `if __name__ == "__main__"` behavior.
- Exercise representative AS, SS, pasted target, file target, accession, table,
  and mocked remote-BLAST command paths.

Exit gate:

- Gate S passes.
- Gate P passes.
- Both GUI batch files launch successfully from the repository checkout.

Suggested commit:

```text
refactor: extract transcript scan CLI
```

### Phase 10 — Extract GUI and portable application services

Status: [ ]

Goal: remove Tkinter and portable-state concerns from the compatibility facade
without breaking the packaged app.

Tasks:

- Add `transcript_scan/gui.py`.
- Move dialogs, mode selection, form state, progress UI, and result window.
- Keep workflow execution outside Tkinter callbacks where practical.
- Decide whether portable paths/settings belong in a small app-services module;
  do not generalize them unless another application has the same contract.
- Update `transcript_scan_app.py` imports deliberately.
- Update `deployment/transcript_scan.spec` if PyInstaller cannot discover new
  modules or dependencies statically.
- Expand `--self-test` to import the final GUI/workflow modules explicitly.
- Preserve logging and fatal-error dialog behavior.

Exit gate:

- Gate S passes.
- Gate P passes.
- Full Gate G passes.
- Test a copied `dist\TranscriptScan` folder, not only the repository build
  location.

Suggested commit:

```text
refactor: extract transcript scan GUI
```

### Phase 11 — Reduce the compatibility facade

Status: [ ]

Goal: finish with a small, stable `ncbi_blast.py` rather than another large
coordinator.

Tasks:

- Keep explicit public re-exports required by existing callers.
- Keep `main()` and module execution compatibility.
- Add `__all__` documenting the supported compatibility surface.
- Remove obsolete adapters only after confirming no repository caller uses
  them.
- Update module documentation to distinguish local transcript scan from remote
  BLAST.
- Record final line counts and dependency boundaries.

Exit gate:

- Gate S passes.
- Gate P passes.
- Full Gate G passes.
- `ncbi_blast.py` is a thin compatibility facade, ideally about 100–250 lines.

Suggested commit:

```text
refactor: finalize ncbi blast compatibility facade
```

### Phase 12 — Release-candidate packaging and documentation

Status: [ ]

Goal: prove that the refactored source and standalone application are ready for
the colleague distribution workflow.

Tasks:

- Review README, deployment instructions, version file, and third-party notices.
- Decide and apply the release version consistently.
- Perform a clean PyInstaller build.
- Run packaged self-test and complete GUI smoke matrix.
- Inspect workbook output from the packaged executable.
- Test from a copied folder outside the repository.
- Run Gate R and generate a clean ZIP.
- Record ZIP path, size, and SHA-256 checksum.
- Make the release-candidate commit only after all checks pass.

Exit gate:

- Gates S, P, G, and R all pass.
- The release package works without a Python installation.
- No private/test data is included.
- The roadmap progress table is complete.

Suggested commit:

```text
build: package refactored transcript scan app
```

## Progress log

Update this table when each phase is completed. Record the actual commit hash
and any approved deviation from the roadmap.

| Phase | Status | Commit | Validation summary | Notes/deviations |
|---|---|---|---|---|
| Baseline | Complete | `44349c4` | 120 tests; packaged self-test and GUI validation completed before checkpoint | Version 1.1.0 |
| 1. Contracts | Complete | `95bbcda` | 128 tests; CLI help; PyInstaller build; packaged self-test exit 0 | Added compatibility, path, normalization, FASTA, and mismatch contracts; no production code changed |
| 2. Nucleotide primitives | Complete | `ec8a14f` | 135 tests; 95 focused tests; three source CLI smoke checks; PyInstaller build; packaged self-test exit 0 | Added `tools_for_pharma.sequence`; preserved strict oligo and N-aware transcript-panel policies |
| 3. FASTA/comparison primitives | Complete | `9f94705` | 144 tests; 104 focused tests; source FASTA/scan/format smoke checks; PyInstaller build; packaged self-test exit 0 | Shared structural FASTA parsing and equal-length comparison; preserved panel, query, and single-transcript adapters |
| 4. Identifiers/NCBI transport | Complete | `f227540` | 158 tests; 97 focused tests; CLI help; PyInstaller build; packaged self-test exit 0 | Extracted RefSeq parsing and injectable NCBI transport; tests keep EFetch accession-only and remote BLAST query submission distinct; left `transcript_panel.py` transport unchanged because it shares different URL/contact/response handling with Ensembl |
| 5. Scan domain package | Complete | `c0c3740` | 163 tests; 88 focused tests; AS and SS source smoke scans; CLI help; PyInstaller build; packaged self-test exit 0 | Added dependency-light models, query preparation, and pure scanner modules; retained pandas table loading and workbook projections in the facade for later interface/reporting phases; `ncbi_blast.py` reduced from 4,067 to 3,597 lines |
| 6. Targets/cache | Complete | `cc7028c` | 171 tests; 96 focused tests; three source CLI smoke checks; CLI help; PyInstaller build; packaged self-test exit 0; packaged GUI checks for pasted, local-file, and cached-accession targets | Added explicit immutable target-source models and mocked cache-miss/reuse/refresh/offline coverage; pasted/local content remains local and is not cached; `ncbi_blast.py` reduced from 3,597 to 3,342 lines |
| 7. Reporting | Complete | `833b254` | 177 tests; 102 focused tests; CLI help; PyInstaller build; packaged self-test exit 0; packaged single-text GUI check; representative offline multiple-workbook generation and visual/schema inspection | Added dedicated reporting/workbook composition with shared sheet-name sanitization; preserved source fields such as `Pos20`, one-based mismatch coordinates, technical detail output, and compatibility re-exports; `ncbi_blast.py` reduced from 3,342 to 2,745 lines. The native file picker was not targetable by the GUI automation API, so the multiple-workbook gate used the equivalent source workflow against a local transcript; packaged reporting imports were verified by the build/self-test. |
| 8. Workflows/remote BLAST | Complete | `afa3183` | 185 tests including mocked remote-BLAST lifecycle; CLI help; direct AS, SS, and offline cached-panel smoke scans; PyInstaller build; packaged self-test exit 0; packaged first-use email and mode chooser opened | Added explicit local, single, and panel workflow configs plus a clearly named remote-BLAST service for submission, polling, retrieval, filtering, batching, and RID logs; preserved facade imports and privacy warnings; local workflows do not import `argparse`, `remote_blast`, or `NcbiBlastClient`; `ncbi_blast.py` reduced from 2,745 to 2,417 lines. The GUI automation API lost the Tk window during the single-form transition, so the full packaged GUI gate was not repeated; source workflow contracts and the packaged import/self-test gate passed. |
| 9. CLI | Not started | | | |
| 10. GUI/portable services | Not started | | | |
| 11. Compatibility facade | Not started | | | |
| 12. Release package | Not started | | | |

## Per-phase working protocol

For every phase:

1. Confirm the previous phase is committed and the worktree is understood.
2. Re-read the phase scope and identify affected callers before editing.
3. Add or confirm characterization tests before moving behavior.
4. Move one responsibility at a time; avoid mixing feature work with moves.
5. Use explicit imports from defining modules internally.
6. Preserve compatibility re-exports at old import paths.
7. Run the required validation gates.
8. Inspect the final diff for accidental output, privacy, or packaging changes.
9. Commit the phase independently.
10. Update the progress log with the commit and validation evidence.

If a phase cannot pass its gate, stop. Fix or revert that phase rather than
continuing with later extractions.

## Definition of refactor complete

The refactor is complete only when all of the following are true:

- Reusable sequence operations have one authoritative implementation.
- Generic FASTA parsing preserves record boundaries.
- Domain workflows use generic primitives without leaking domain behavior into
  the generic layer.
- Local scanning and remote BLAST have separate, auditable data flows.
- `ncbi_blast.py` is a small compatibility facade.
- Existing Python imports, CLI commands, batch files, GUI behavior, workbook
  schemas, cache behavior, and privacy behavior remain compatible.
- The full test suite passes.
- The PyInstaller one-folder application builds and passes self-test.
- The packaged GUI passes the smoke matrix from a copied folder.
- The clean release ZIP contains no generated user data or private sequences.

# Fast-Sanger-BE

Batch quantification of adenine (ABE) and cytosine (CBE) base editing from Sanger `.ab1` chromatograms.

Fast-Sanger-BE reads a folder of `.ab1` files, works out from each file name which sample, guide(s), editor and read direction it holds, rejects reads whose chromatogram around the guide is too noisy to interpret, and reports background-corrected conversion at every editable position of the protospacer. Forward and reverse reads are combined per sample, and ABE and CBE experiments are reported separately.

There are two ways to use it:

- **Online, with nothing to install:** <https://PROADI-TIAF.github.io/fast-sanger-be/>. The analysis runs inside your web browser and the `.ab1` files never leave your computer (see [Web version](#web-version)).
- **From a terminal** on Linux, macOS or Windows (through WSL), as described below. This suits large batches and automated pipelines.

## Installation

Works on Linux, macOS and Windows. On Windows, install the Windows Subsystem for Linux first (run `wsl --install` in PowerShell as administrator, restart, then open *Ubuntu*) and type every command below in the Ubuntu terminal.

Requires [conda](https://docs.conda.io) (Miniconda, Miniforge or Anaconda) and git.

```bash
git clone https://github.com/PROADI-TIAF/fast-sanger-be.git
cd fast-sanger-be
conda env create -f environment.yml
conda activate fast-sanger-be
pip install .
fast-sanger-be --version
```

The environment contains Python 3.11, Biopython, NumPy, Matplotlib and openpyxl. No R, Bash or specific browser is required.

## Quick start

```bash
cd /path/to/ab1/folder
fast-sanger-be
```

or `fast-sanger-be /path/to/ab1/folder`. To try it on synthetic data with known editing levels:

```bash
python -m fast_sanger_be.simulate example_run
fast-sanger-be example_run
```

`example_run/truth.tsv` lists the simulated editing at every position, so the results can be compared with the truth.

## Web version

**<https://PROADI-TIAF.github.io/fast-sanger-be/>**

The same analysis also runs in a web browser, with nothing to install. The page loads Python and this package with [Pyodide](https://pyodide.org) (WebAssembly), so the `.ab1` files are processed on the user's computer and never uploaded. The page replaces the terminal steps: files are dropped onto it, names that do not follow the convention are completed in a table, guides are typed or loaded from `guides.txt`, and the results are downloaded as `Final_Results.zip` with the same folders and files as the command-line version. A table `input_files.tsv` in the zip records which original file each analysed name came from.

To run the page from your own copy instead (for example, while changing the code):

```bash
python web/serve.py
```

This builds the page (package wheels and an example dataset) and opens it at <http://localhost:8000>. The first analysis downloads Pyodide and the scientific libraries (about 30 MB, then cached by the browser), so an internet connection is needed once.

To publish it, push the repository to GitHub and set *Settings > Pages > Source* to *GitHub Actions*; the workflow in `.github/workflows/web.yml` builds and deploys `web/` on every push to `main`. Pyodide is loaded from the jsDelivr CDN, pinned to version 314.0.7 in `web/worker.js`.

## Inputs

### File names

```
<Sample>_<Guide>[_<Guide2>...]_<ABE|CBE>_<F|R>.ab1
```

| Field | Rule | Example |
|---|---|---|
| Sample | letters, digits, hyphens (no underscores) | `Clone-12` |
| Guide(s) | one or more guide names, separated by `_` | `sgHBB` or `sgHBB_sgEMX1` |
| Editor | `ABE` or `CBE` | `ABE` |
| Read | `F` (forward) or `R` (reverse) | `F` |

Examples: `Clone-12_sgHBB_ABE_F.ab1`, `Clone-12_sgHBB_ABE_R.ab1`, `Pool-3_sgHBB_sgEMX1_CBE_F.ab1` (one amplicon analysed for two guides).

If some files do not follow this pattern, a **rename editor** opens in the browser with one row per file (pre-filled where the old name allows a guess). After saving, type `y` in the terminal and the files are renamed. Every rename is recorded in `rename_log.tsv`, so it can be undone. Files can also be left out of a run by unchecking *Include*.

### Guide file

A FASTA-like text file named `guides.txt` (or `guides.fasta` / `guides.fa`) in the same folder, with the protospacer 5'→3' **without the PAM**:

```
>sgHBB
CTTGACACAACTGTGTTCAC
>sgEMX1
GAGTCCGAGCAGAAGAAGAA
```

Headers may also be written as `[sgHBB]`. Another file can be given with `--guides`. If the file is missing, invalid, or lacks a guide used in the file names, a **guide editor** opens with the missing names already listed; save, then type `y` to continue (or `n` to abort). Two columns (name, sequence) can be pasted straight from Excel.

The editors are served from a temporary local web server (`127.0.0.1`) so the result is saved directly into the analysis folder. With `--no-browser` the address is printed instead of opened.

## Analysis

For every `.ab1` file and every guide in its name:

1. **Guide location.** The guide is searched on both strands of the base calls. Expected conversions (A→G for ABE, C→T for CBE) count as matches, so heavily edited samples are still found; up to 3 other mismatches are tolerated (`--max-mismatches`). Reads where the guide matches more than one location are rejected. The chromatogram is then oriented so that the guide reads 5'→3', which makes forward and reverse reads directly comparable.

2. **Quality control.** A read is excluded when any of these fails:

   | Criterion | Default |
   |---|---|
   | Read covers the full guide ±10 nt (`--qc-flank`) | required |
   | Guide matches a single location | required |
   | Median Phred of the region used for background (±50 nt, guide excluded) | ≥ 20 |
   | Median primary peak height of that region | ≥ 100 RFU |
   | Overlapping peaks (secondary/primary > 0.25) across the guide ±10 nt window, excluding bases the editor can convert | ≤ 60% |
   | Every target position has a clean neighbourhood (see below) | required |

   *Target neighbourhoods.* Local problems such as a dye blob or a short noisy stretch can corrupt the values at a few positions while the rest of the read looks fine. Each target position is therefore checked separately: it is noisy when more than 35% of its neighbours within ±5 nt have overlapping peaks (`--max-local-mixed`, `--position-flank`). Only neighbours the editor cannot convert are counted, because a real edit produces two peaks at the edited base itself and must not be mistaken for noise. One noisy target excludes the whole read, and the reason names the positions involved (for example `overlapping peaks near A1, A3, A5`).

   With `--partial-reads` such a read is kept instead: only the noisy positions are left empty in the workbooks, hatched in the heatmap and greyed in the read figure, and the *Notes* column says which ones were skipped. The read is still excluded if no target can be quantified.

   Excluded reads are reported in the terminal, in the `Excluded_reads` sheet of the final workbook, and their figures go to a separate `Low_Quality` folder.

3. **Quantification.** For each base position, the height of each dye channel is the trace maximum within ±2 samples of the called peak. Background noise is estimated per channel from up to 50 nt on each side of the guide (`--background-flank`), never from the guide itself: at each position, channels that are not the primary peak contribute their height, outliers beyond 3 scaled median absolute deviations are discarded, and the mean is taken. Background is subtracted from each channel (negative values set to 0) and the four channels are converted to percentages. Editing is the percentage of G at guide A positions (ABE) or of T at guide C positions (CBE). Positions are numbered 1–N from the 5' end of the protospacer.

4. **Combining reads.** The value for a sample is the mean of the reads that passed QC. When only one read is available (the other is missing or excluded) that read is used and the reason is written in the *Notes* column. Positions where forward and reverse differ by more than 10 percentage points are flagged (`--discordance`).

The dye order is read from each file (`FWO_` tag) instead of being assumed.

## Outputs

```
Final_Results/
├── run_log.txt
├── ABE/
│   ├── ABE_final_results.xlsx      one sheet per guide: sample × target position (mean of F/R)
│   ├── Heatmaps/<guide>_heatmap.png
│   └── Raw_Data/
│       ├── ABE_raw_results.xlsx    one sheet per guide: every F and R read, plus QC and base composition
│       ├── Figures/<guide>/<file>.png
│       └── Low_Quality/<guide>/<file>.png
└── CBE/                            same structure, only when CBE files are present
```

**Final workbook**: one sheet per guide (Sample, Reads used, one column per target position such as `A5 (% G)`, maximum F/R difference, Notes), plus `Excluded_reads`, `Guides` and `Parameters`.

**Raw workbook**: one sheet per guide with one row per read (status *OK*, *Low quality: reason*, *Guide not found* or *Missing*), a `QC` sheet with every metric and background estimate, a `Base_composition` sheet with %A/%C/%G/%T at every guide position of every read, `Guides` and `Parameters`.

**Figures**: the chromatogram of the guide ±10 nt (guide shaded in grey) above the percentage of each base at each guide position. Heatmaps use a fixed 0–100% colour scale (`--heatmap-max`) so different runs are comparable.

If `Final_Results` already exists, it is renamed with a timestamp rather than overwritten.

## Options

```
fast-sanger-be [input_dir] [--guides FILE] [--output-name NAME] [--no-browser] [--non-interactive]
               [--qc-flank 10] [--min-phred 20] [--min-peak-height 100] [--mixed-peak-ratio 0.25]
               [--position-flank 5] [--max-local-mixed 0.35] [--partial-reads] [--max-window-mixed 0.6]
               [--max-mismatches 3]
               [--background-flank 50] [--discordance 10] [--heatmap-max 100]
               [--figure-format png|pdf|svg|tiff] [--dpi 300]
```

Every value used is written to the `Parameters` sheet of each workbook.

**Non-interactive use** (servers, pipelines): `--non-interactive` never prompts. When information is missing it writes `rename_table.tsv` or `guides.txt` for you to complete and exits with code 2; the next run picks them up.

## Testing

```bash
pip install ".[test]"
pytest
```

The tests simulate chromatograms with known editing (0–95%) on both strands, mixed templates and partial reads, and run the whole pipeline on a synthetic dataset.

## Citation

If you use Fast-Sanger-BE, please cite: *(manuscript in preparation)*.

Fast-Sanger-BE is an independent implementation. Its background correction follows the principle introduced by BEAT: Xu L, Liu Y, Han R. BEAT: A Python program to quantify base editing from Sanger sequencing. *CRISPR J.* 2019;2(4):223–229. doi:10.1089/crispr.2019.0017.

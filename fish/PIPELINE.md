# Magnetic Field Experiment Analysis Pipeline

## Overview

This document summarizes the current analysis workflow under
`SpikeTools-for-IntanLabview/fish`.

The recommended workflow is MATLAB-first:

1. MATLAB: generate a unified CSV from DeepLabCut `.mat` session files.
2. MATLAB: run the main report and statistical analysis with `report_magnet_stats.m`.
3. Optional R: run `magnet_circular_mixed_model.R` for additional circular mixed-model summaries.

`report_magnet_stats.m` is the primary analysis entry point after CSV generation.
The R script is useful as a secondary or confirmatory analysis, but it is not the
default report workflow.

---

## Step 1: Generate the CSV in MATLAB

**Script:** `magnetStatsCSV2025b.m`  
**Main dependency:** `magnetBehav2.m`

### Usage

```matlab
% Salmon
[ALL, results_stats] = magnetStatsCSV2025b("salmon");

% Trout
[ALL, results_stats] = magnetStatsCSV2025b("trout");
```

If the argument is omitted, the default is currently `salmon`:

```matlab
[ALL, results_stats] = magnetStatsCSV2025b();
```

### What the script does

1. Selects `mat_files_salmon` or `mat_files_trout` based on the `species` argument.
2. Calls `magnetBehav2()` for each `.mat` session.
3. Builds one combined table with three row types:
   - `sample`: frame-level time series
   - `window`: per-fish window summaries
   - `intensity`: per-window aggregate rows
4. Adds pairing keys such as `PairKey_Fish` and `PairKey_FishWin`.
5. Adds `OrderInSession` for each cohort.
6. Adds `AgeGroup` from the fourth column of `mat_files`.

### Main output

- `speed_heading_data2025_salmon-all.csv`
- `speed_heading_data2025_trout-all.csv`

### Important CSV columns

| Column | Meaning |
|---|---|
| `Level` | `sample`, `window`, or `intensity` |
| `Session` | source `.mat` file name |
| `Condition` | e.g. `beringSea`, `okhotskSea`, `controlSea1`, `microT10`, `microT100`, `controlT1` |
| `Cohort` | cohort identifier |
| `AgeGroup` | e.g. `2mo`, `5mo` |
| `Fish` | fish identifier used for paired analysis |
| `Track` | track index within session |
| `OrderInSession` | session order within a cohort |
| `Time` | seconds from session start |
| `HeadingMag` | head direction relative to magnetic north |
| `cosMag`, `sinMag` | circular-model components |
| `Speed`, `TurnRate`, `Curvature` | activity metrics |
| `Intensity_uT` | assigned magnetic-field intensity |
| `PairKey_Fish` | fish pairing key across conditions |
| `PairKey_FishWin` | fish-plus-window pairing key |

---

## Step 2: Main MATLAB Report Analysis

**Primary script:** `magnet_report_code/report_magnet_stats.m`  
**Core implementation:** `magnet_report_code/report_magnet_stats_core.m`

Before running the report, add the report-code folder to the MATLAB path:

```matlab
addpath(fullfile(pwd, "magnet_report_code"))
```

If MATLAB is not currently in the `fish` folder, use the full path instead:

```matlab
addpath("C:\Users\stakahas.DESKTOP-JMVVEDC\Documents\SpikeTools-for-IntanLabview\fish\magnet_report_code")
```

### Basic usage

```matlab
OUT = report_magnet_stats(csvPath, conditions, winStartSec, winLenSec, ...
    "outdir", outdir, ...
    "doLME", true);
```

`winStartSec` is the start of the analysis window.  
`winLenSec` is the window length. For example, `300, 600` analyzes 300-900 s.

### Salmon sea-condition example

```matlab
csvPath = "speed_heading_data2025_salmon-all.csv";

OUT = report_magnet_stats( ...
    csvPath, ...
    ["controlsea1", "beringsea", "okhotsksea"], ...
    300, 600, ...
    "outdir", "magnet_report_salmon_sea_300_900", ...
    "doLME", true);
```

### Salmon intensity-condition example

```matlab
csvPath = "speed_heading_data2025_salmon-all.csv";

OUT = report_magnet_stats( ...
    csvPath, ...
    ["controlt1", "microt10", "microt100"], ...
    300, 600, ...
    "outdir", "magnet_report_salmon_intensity_300_900", ...
    "doLME", true);
```

### Trout example

Use the same report commands with:

```matlab
csvPath = "speed_heading_data2025_trout-all.csv";
```

### What `report_magnet_stats.m` does

The report workflow uses `Level == "sample"` rows from the CSV and then:

1. Filters rows by the requested time window.
2. Filters rows to the requested conditions.
3. Keeps only fish that are paired across all requested conditions.
4. Computes additional metrics such as:
   - `Distance10s`
   - 10 s tortuosity
5. Aggregates per fish and condition.
6. Runs paired tests for continuous metrics:
   - `Speed`
   - `Distance10s`
   - `TurnRate`
   - `Curvature`
   - `Tortuosity`
7. Runs paired tests for heading concentration `R`.
8. Runs paired tests on fish-mean `cosMag` and `sinMag`.
9. Runs circular pairwise tests on fish mean headings.
10. Runs Rayleigh uniformity tests across fish means.
11. Optionally fits mixed models for activity metrics, `cosMag`, and `sinMag`.
12. Produces report plots and summary files.

### Main helper functions

`report_magnet_stats.m` delegates most work to functions in `magnet_report_code`,
including:

- `report_magnet_stats_core.m`
- `per_fish_cond_means.m`
- `per_fish_cond_circular.m`
- `fishmeans_circ_pairwise_tests.m`
- `rayleigh_uniformity_across_fish.m`
- `target_direction_vtest.m`
- `fitglme_safe.m`
- `fitlm_safe.m`
- `draw_rose_pdf_all.m`
- `draw_R_visualization.m`
- `draw_fish_means_rose_multi.m`
- `draw_paired_raincloud_allmetrics.m`
- `draw_circ_pairwise_with_raw.m`
- `plot_circ_pairwise_overview.m`
- `draw_lme_anova_tables.m`
- `summarize_results.m`

### Important outputs

The exact files depend on options, but the report directory typically includes:

- `summary.txt`
- `rose_headingMag.pdf`
- fish-mean rose plots
- paired metric plots
- circular pairwise plots
- LME ANOVA table PDFs

The returned MATLAB struct `OUT` also contains the key tables:

- `OUT.pairwise.ttest_cont`
- `OUT.pairwise.ttest_R`
- `OUT.pairwise.ttest_proj`
- `OUT.circPair`
- `OUT.uniformityRayleigh`
- `OUT.lme.continuous`
- `OUT.lme.cos`
- `OUT.lme.sin`
- `OUT.groupR`
- `OUT.files`
- `OUT.summary`

---

## Step 3: Optional R Mixed-Model Analysis

**Script:** `magnet_circular_mixed_model.R`

Use the R script when you want an additional command-line circular mixed-model
summary, bootstrap circular mean-difference intervals, or R-based output tables.
This step is optional and should be interpreted alongside the MATLAB report.

### Command-line usage

```powershell
Rscript magnet_circular_mixed_model.R <csv_path> <out_dir> [window_start] [window_end] [condition_filter]
```

### Example

```powershell
Rscript magnet_circular_mixed_model.R `
  speed_heading_data2025_salmon-all.csv `
  r_salmon_sea_300_600 `
  300 `
  600 `
  "controlsea1,beringsea,okhotsksea"
```

Important difference from the MATLAB report:

- In `report_magnet_stats.m`, `300, 600` means 300-900 s.
- In `magnet_circular_mixed_model.R`, `300, 600` means 300-600 s.

### Typical R outputs

- `circular_model_input.csv`
- `circular_mean_difference_ci.csv`
- `cos_mixed_model_summary.txt`
- `sin_mixed_model_summary.txt`
- `cos_mixed_model_anova.csv`
- `sin_mixed_model_anova.csv`
- `individual_alignment_rows.csv`
- `individual_alignment_summary.csv`
- `bayesian_circular_summary.txt`
- `README_circular_models.txt`

---

## Recommended Execution Order

### Salmon

1. Generate the CSV:

```matlab
[ALL, results_stats] = magnetStatsCSV2025b("salmon");
```

2. Run the main MATLAB sea-condition report:

```matlab
addpath(fullfile(pwd, "magnet_report_code"))

OUT_sea = report_magnet_stats( ...
    "speed_heading_data2025_salmon-all.csv", ...
    ["controlsea1", "beringsea", "okhotsksea"], ...
    300, 600, ...
    "outdir", "magnet_report_salmon_sea_300_900", ...
    "doLME", true);
```

3. Run the main MATLAB intensity report:

```matlab
OUT_intensity = report_magnet_stats( ...
    "speed_heading_data2025_salmon-all.csv", ...
    ["controlt1", "microt10", "microt100"], ...
    300, 600, ...
    "outdir", "magnet_report_salmon_intensity_300_900", ...
    "doLME", true);
```

4. Optionally run the R script for secondary analyses.

### Trout

1. Generate the CSV:

```matlab
[ALL, results_stats] = magnetStatsCSV2025b("trout");
```

2. Run the same report commands with:

```matlab
"speed_heading_data2025_trout-all.csv"
```

---

## Notes and Caveats

1. `report_magnet_stats.m` lowercases condition names internally, so lowercase
   condition names are recommended in calls.
2. `report_magnet_stats.m` uses a start time plus window length. The R script
   uses start time plus end time.
3. The report workflow keeps only fish that are present in all requested
   conditions, so `OUT.N_pairs_common` is an important quality-control value.
4. `doLME=true` requests mixed models. If mixed-model fitting fails, the report
   code attempts safer fallbacks where possible.
5. Keep the `magnet_report_code` folder on the MATLAB path when using
   `report_magnet_stats.m`.
6. Non-ASCII OneDrive paths can be troublesome for some R and MATLAB commands.
   A short local path or junction can be safer when needed.

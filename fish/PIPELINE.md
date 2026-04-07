# Magnetic Field Experiment Analysis Pipeline

## Overview

This document summarizes the current analysis workflow under `SpikeTools-for-IntanLabview/fish`.

The pipeline has two main stages and one optional descriptive step:

1. MATLAB: generate a unified CSV from DeepLabCut `.mat` session files
2. Optional R summary: descriptive behavioral statistics
3. R mixed-model analysis: circular heading, activity, age-group, and order effects

---

## Step 1: MATLAB CSV Generation

**Script:** `magnetStatsCSV2025b.m`  
**Dependency:** `magnetBehav2.m`

### Usage

```matlab
% Salmon
[ALL, results_stats] = magnetStatsCSV2025b("salmon");

% Trout
[ALL, results_stats] = magnetStatsCSV2025b("trout");
```

If the argument is omitted, the default is:

```matlab
[ALL, results_stats] = magnetStatsCSV2025b();
```

which currently means `salmon`.

### What the script does

1. Selects `mat_files_salmon` or `mat_files_trout` based on the `species` argument.
2. Calls `magnetBehav2()` for each session.
3. Builds one combined table with three row types:
   - `sample`: frame-level time series
   - `window`: per-fish window summaries
   - `intensity`: per-window aggregate rows
4. Adds paired-analysis keys such as `PairKey_Fish` and `PairKey_FishWin`.
5. Adds `OrderInSession` for each cohort.
6. Adds `AgeGroup` from the fourth column of `mat_files`.

### Main output

- `speed_heading_data2025_salmon-all.csv`
- `speed_heading_data2025_trout-all.csv`

The output file name is selected automatically from the `species` argument.

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
| `Speed`, `Distance10s`, `TurnRate`, `Curvature` | activity metrics |
| `PairKey_Fish` | fish pairing key across conditions |
| `PairKey_FishWin` | fish-plus-window pairing key |

---

## Step 2: Optional Descriptive Statistics

**Script:** `salmonBehav2025.R`  
**Wrapper:** `run_salmonBehav2025.ps1`

### Typical usage

```powershell
.\run_salmonBehav2025.ps1 -CsvPath "path\to\speed_heading_data2025_salmon-all.csv"
```

or directly:

```powershell
Rscript salmonBehav2025.R <csv_path> <out_dir>
```

### Main outputs

- `salmonBehav2025_summary.csv`
- `salmonBehav2025_pairwise_tests.csv`

This step is descriptive and optional. The main inferential workflow is in Step 3.

---

## Step 3: R Mixed-Model Analysis

**Script:** `magnet_circular_mixed_model.R`

### Required R packages

- `lme4`
- `lmerTest` (optional but preferred)
- `bpnreg`

### Command-line usage

```powershell
Rscript magnet_circular_mixed_model.R <csv_path> <out_dir> [window_start] [window_end] [condition_filter]
```

### Arguments

| Argument | Meaning |
|---|---|
| `csv_path` | CSV produced by `magnetStatsCSV2025b.m` |
| `out_dir` | output directory |
| `window_start` | start time in seconds |
| `window_end` | end time in seconds |
| `condition_filter` | comma-separated lowercase condition names |

### Example runs

#### Salmon Sea conditions

```powershell
$rscript = "C:\Program Files\R\R-4.5.3\bin\Rscript.exe"
$script  = "J:\OneDrive - 同志社大学\Documents\MATLAB\SpikeTools-for-IntanLabview\fish\magnet_circular_mixed_model.R"
$csv     = "J:\...\speed_heading_data2025_salmon-all.csv"

& $rscript $script $csv "J:\...\r_300_1000" 300 1000 "controlsea1,beringsea,okhotsksea"
& $rscript $script $csv "J:\...\r_600_900" 600 900 "controlsea1,beringsea,okhotsksea"
```

#### Salmon Intensity conditions

```powershell
& $rscript $script $csv "J:\...\t_300_1000" 300 1000 "controlt1,microt10,microt100"
& $rscript $script $csv "J:\...\t_600_900" 600 900 "controlt1,microt10,microt100"
```

#### Trout

Use the same commands with:

```powershell
$csv = "J:\...\speed_heading_data2025_trout-all.csv"
```

### What the script does

The script now includes all of the following:

1. Circular pairwise heading comparisons with bootstrap CI
2. Mixed models for `cosMag` and `sinMag`
3. Bayesian circular model output
4. Activity pairwise tests for available metrics
5. Cohort-group summaries
6. Age-group summaries
7. Within-fish alignment and shift analyses
8. Shift clustering and cluster plot
9. General `OrderInSession` activity analyses
10. Subset-specific `OrderInSession` models for:
   - `Sea x 2mo`
   - `Sea x 5mo`
   - `Intensity x 2mo`
   - `Intensity x 5mo`

### Important output files

Core files:

- `circular_model_input.csv`
- `circular_mean_difference_ci.csv`
- `results_interpretation.txt`
- `README_circular_models.txt`

Within-fish outputs:

- `within_fish_reference_correlations.csv`
- `within_fish_shift_summary.csv`
- `within_fish_shift_distribution.csv`
- `within_fish_shift_cohort_summary.csv`
- `within_fish_shift_age_summary.csv`
- `within_fish_shift_pair_correlations.csv`
- `within_fish_shift_pair_cohort_correlations.csv`
- `within_fish_shift_pair_age_correlations.csv`
- `within_fish_shift_cluster_assignments.csv`
- `within_fish_shift_cluster_summary.csv`
- `within_fish_shift_cluster_stats.csv`
- `within_fish_shift_cluster_centers.csv`
- `within_fish_shift_clusters.png`

Activity and model outputs:

- `activity_pairwise_tests.csv`
- `activity_mixed_model_terms_summary.csv`
- `order_activity_pairwise_tests.csv`
- `order_activity_mixed_model_terms_summary.csv`
- `age_group_pairwise_effects.csv`
- `age_group_activity_pairwise_tests.csv`
- `age_mixed_model_terms_summary.csv`
- `subset_order_activity_pairwise_tests.csv`
- `subset_order_activity_mixed_model_terms_summary.csv`
- `subset_order_activity_model_meta.csv`

Subset-specific `OrderInSession` model files:

- `subset_order_sea_2mo_<metric>_anova.csv`
- `subset_order_sea_5mo_<metric>_anova.csv`
- `subset_order_intensity_2mo_<metric>_anova.csv`
- `subset_order_intensity_5mo_<metric>_anova.csv`

plus matching `_coefficients.csv` and `_summary.txt` files.

---

## Recommended Execution Order

### Salmon

1. Generate CSV in MATLAB:

```matlab
[ALL, results_stats] = magnetStatsCSV2025b("salmon");
```

2. Run Sea analyses:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" `
  "J:\OneDrive - 同志社大学\Documents\MATLAB\SpikeTools-for-IntanLabview\fish\magnet_circular_mixed_model.R" `
  "J:\...\speed_heading_data2025_salmon-all.csv" `
  "J:\...\r_300_1000" `
  300 `
  1000 `
  "controlsea1,beringsea,okhotsksea"
```

3. Run Intensity analyses:

```powershell
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" `
  "J:\OneDrive - 同志社大学\Documents\MATLAB\SpikeTools-for-IntanLabview\fish\magnet_circular_mixed_model.R" `
  "J:\...\speed_heading_data2025_salmon-all.csv" `
  "J:\...\t_300_1000" `
  300 `
  1000 `
  "controlt1,microt10,microt100"
```

### Trout

1. Generate CSV in MATLAB:

```matlab
[ALL, results_stats] = magnetStatsCSV2025b("trout");
```

2. Run the same R commands, but with:

```powershell
"J:\...\speed_heading_data2025_trout-all.csv"
```

---

## Notes and Caveats

1. `condition_filter` should be lowercase because the R script lowercases condition names internally.
2. Some historical output folder names may not match the exact numeric window in the arguments. The arguments are the true source of the analysis window.
3. Non-ASCII OneDrive paths can still be troublesome for some R and MATLAB commands. Junctions or short proxy paths are safer when needed.
4. `OrderInSession` effects in Sea subsets should be interpreted together with the subset-specific outputs, not only the global activity models.
5. The most direct tests for procedural artifacts are now the subset-specific order models in `Sea x AgeGroup` and `Intensity x AgeGroup`.


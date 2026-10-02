```text
8888888b.                         d8b 8888888b.                   
888  "Y88b                        Y8P 888  "Y88b                  
888    888                            888    888                  
888    888  .d88b.  88888b.d88b.  888 888    888  .d88b.  88888b. 
888    888 d88""88b 888 "888 "88b 888 888    888 d8P  Y8b 888 "88b
888    888 888  888 888  888  888 888 888    888 88888888 888  888
888  .d88P Y88..88P 888  888  888 888 888  .d88P Y8b.     888  888
8888888P"   "Y88P"  888  888  888 888 8888888P"   "Y8888  888  888
```

**DomiDen** (Dominance analysis of Denudation rates) is a set of R scripts that
estimates the relative importance (general dominance, LMG) of groups of
environmental predictors, such as topography, climate and seismicity, on
catchment denudation rates. It uses an ensemble of random OLS models ranked by AIC.

## Key features

- **Any number of variables.** Use as many predictors as your tables contain,
  grouped into any number of classes. You only declare them once in the `config`
  block of `scripts/03_setup_config.R`.
- **Only one decision per predictor.** Decide whether the variable suits a
  **log**, **arcsin(sqrt)** or **logit** transformation. DomiDen does the rest.
- **Any number of datasets.** Each table in `data/` is analysed with the same
  model ensemble and shown side by side.
- **No absolute paths.** Everything is relative to the project root.

## Folder structure

All paths are relative to the project root `DomiDen/`, detected with the `here` package.

DomiDen/
├── DomiDen.Rproj
├── README.md
├── scripts/
│ ├── 01_DA_functions.R <- dominance analysis functions
│ ├── 02_model_generator.R <- random model ensemble generator
│ ├── 03_setup_config.R <- packages, config, data loading
│ └── 04_run_and_plot.R     <- RUN THIS: analysis, tables, figures
├── data/    <- your input tables (.xlsx / .csv)
├── results/ <- CSV outputs (created automatically)
└── figures/ <- PDF figures (created automatically)

## Requirements

R >= 4.x and these packages:

```r
install.packages(c("here", "relaimpo", "dplyr", "tidyr", "tibble",
                   "ggplot2", "readxl", "knitr"))
```

## How to run

1. Download or clone the `DomiDen` folder and keep its structure unchanged.
2. Open `DomiDen.Rproj` in RStudio. From a terminal, `cd DomiDen` instead.
3. Put your data tables in `data/`. Each row is a catchment, with one response
   column and one column per predictor.
4. Edit **only** the `config` block (section 2) in `scripts/03_setup_config.R`.
5. Open `scripts/04_run_and_plot.R` and run it (Source button in RStudio), or use
   `source("scripts/04_run_and_plot.R")` in the console, or
   `Rscript scripts/04_run_and_plot.R` from the terminal in the `DomiDen/` folder.
   Scripts 01-03 are loaded automatically; you never run them separately.
6. Find the outputs in `results/` and `figures/`.

## Configuring your analysis

Everything is set in the `config` list of `scripts/03_setup_config.R`.

| Setting | What to put |
|---|---|
| `response_var` | Column name of the response (e.g. `"DR_tkm2y"`) |
| `response_transform` | `"log"`, `"sqrt"` or `"none"` |
| `n_random_models`, `n_top_models`, `seed` | Ensemble size, models kept after AIC ranking, random seed |
| `pool` | Named list: class name -> vector of predictor column names |
| `Positive_to_Log` | Strictly positive predictors -> `log(x)` |
| `Arcsin_Sqrt_Proportions` | Proportions (0-1 or 0-100) -> `asin(sqrt(x))` |
| `Proportions_to_Logit` | Proportions strictly between 0 and 1 -> `log(x/(1-x))` |
| `datasets` | `"label" = "file name in data/"`, one line per dataset |
| `dataset_colors`, `class_colors`, `predictor_labels` | Optional plot styling (`NULL` for automatic) |

### Using any number of variables

DomiDen does not limit the number of variables. To add a predictor:

1. Add its column name to a class in `pool` (or create a new class).
2. Add the same name to **one** transformation list.

That is all. Models, tables, figure axes, class separators and axis limits adapt
automatically. The only judgement you make is whether each variable is suited to
a log, arcsin or logit transformation. Example:

```r
pool = list(
  Climate = c("MAP", "MAT_K", "NewVariable"),   # <- added
  Biota   = c("NDVI")                           # <- new class
),
Positive_to_Log         = c("MAP", "MAT_K", "NewVariable"),
Arcsin_Sqrt_Proportions = c("NDVI"),
```

DomiDen checks the config before running. It stops with a clear message if a
predictor is missing from the transformation lists, appears in two lists, or is
absent from a data table.

## Input requirements

- The response is strictly positive if `response_transform = "log"`.
- Column names are identical across all datasets.
- Each predictor is in exactly one transformation list.
- Log: values > 0. Arcsin(sqrt): proportions 0-1 or 0-100 (auto-detected).
  Logit: proportions strictly between 0 and 1.

## Outputs

- `results/DominanceResults_byClass.csv` and `results/DominanceResults_byPredictor.csv`:
  mean, SD, SE, 95% CI, median, min, max and number of models.
- `figures/Fig_DomiDen_byClass.pdf` and `figures/Fig_DomiDen_byPredictor.pdf`.

## Troubleshooting

- "cannot open file" / "File not found": run `here::here()`. It must print the path
  ending in `DomiDen`. If it doesn't, create an empty file named `.here` inside
  `DomiDen/`, or open the `.Rproj`. Also check that the file names in
  `config$datasets` match the files in `data/` exactly.
- Script names in `source(here("scripts", ...))` must match the file names exactly,
  including case.
- "Dataset ... lacks column(s)": a name in `pool` is spelled differently in your table.

## Citation and license

MIT License

Copyright (c) 2026 YOUR NAME

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
# ==============================================================================
# DomiDen | 03_setup_config.R
# Packages, USER CONFIGURATION and data loading.
#
# >>> This is the ONLY file you normally need to edit. <<<
# Edit section 2 (config). Sections 1, 3 and 4 run as they are.
#
# VERSION 1.0 | Last update 01.10.2026
# ==============================================================================


# ==============================================================================
# 1. PACKAGES AND FUNCTIONS
# ==============================================================================
required_pkgs <- c("here", "relaimpo", "dplyr", "tidyr", "tibble",
                   "ggplot2", "readxl", "knitr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them with: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))", call. = FALSE)
}

library(here)
library(relaimpo)
library(ggplot2)
library(readxl)
library(tibble)
library(dplyr)
library(tidyr)

message("Project root: ", here::here())

source(here("scripts", "01_DA_functions.R"))     # dominance analysis
source(here("scripts", "02_model_generator.R"))  # random model ensemble


# ==============================================================================
# 2. USER CONFIGURATION  (edit this block)
# ==============================================================================
config <- list(
  
  # --- 2a. Response variable --------------------------------------------------
  response_var       = "DR_tkm2y",  # column name of the response in every table
  response_transform = "log",       # "log" (needs values > 0), "sqrt" or "none"
  
  # --- 2b. Ensemble settings --------------------------------------------------
  n_top_models    = 1000,    # number of best models (lowest AIC) that are kept
  n_random_models = 10000,   # number of random models generated before ranking
  seed            = 123,     # random seed, for reproducible results
  
  # --- 2c. Predictor pool -----------------------------------------------------
  # Any number of classes, and any number of predictors per class.
  # Names must match the column names in your tables exactly.
  # The class names (left side) are used as labels in tables and figures;
  # underscores are shown as spaces. The class order here is the order in the plots.
  pool = list(
    Topography           = c("Ksn", "LocalRelief", "Slope",
                             "Grad_prc90", "Grad_prc10"),
    Seismic              = c("PGA50pc50yr", "Seis_Energy_Mean"),
    Land_cover           = c("Barren_Vegetated", "glacier_fr"),
    Substrate_properties = c("Erodibility_lithological", "permeability",
                             "Prop_UnconsolidatedDeposits", "SoilThickness_m"),
    Climate              = c("MAP", "AnnualMax_daily_pr", "MAT_K", "AridityIndex")
  ),
  
  # --- 2d. Transformations ----------------------------------------------------
  # You only decide whether each predictor suits a log, arcsin or logit
  # transformation. Every predictor in the pool must be in exactly ONE list
  # below. Leave a list as c() if you do not need it.
  
  # Strictly positive values (> 0)   -> log(x)
  Positive_to_Log = c(
    "Slope", "Ksn", "LocalRelief",
    "PGA50pc50yr", "Seis_Energy_Mean",
    "Grad_prc90", "Grad_prc10",
    "Barren_Vegetated",
    "Erodibility_lithological", "SoilThickness_m", "permeability",
    "MAP", "AnnualMax_daily_pr", "MAT_K", "AridityIndex"
  ),
  
  # Proportions, 0 to 1 (0-100 is detected) -> asin(sqrt(x)); zeros are allowed
  Arcsin_Sqrt_Proportions = c(
    "glacier_fr", "Prop_UnconsolidatedDeposits"
  ),
  
  # Proportions strictly between 0 and 1 -> log(x / (1 - x))
  Proportions_to_Logit = c(),
  
  # --- 2e. Datasets (REQUIRED) -------------------------------------------------
  # One line per dataset:  "label" = "file name in data/"
  # The label is mandatory. It is used in every table, facet title and legend,
  # so choose something short and meaningful. Labels must be unique.
  # Example:
  datasets = c("Meteoric 10Be" = "df_DR_Meteoric.csv",
               "In-situ 10Be"  = "df_DR_InSitu.csv",
               "Gauging"       = "df_ER_Gauging.csv"),
  # datasets = c(
  #   "Dataset 1" = "dataset_1.csv",
  #   "Dataset 2" = "dataset_2.csv",
  #   "Dataset 3" = "dataset_3.csv"
  # ),
  
  # --- 2f. Optional plot styling (NULL = automatic) -----------------------------
  # Names must match the dataset labels / pool class names / column names above.
  # Example:
   dataset_colors = c("Meteoric 10Be" = "#009590", "In-situ 10Be" = "#DEC0A3", "Gauging" = "grey90"),
  #dataset_colors   = NULL,
  class_colors     = NULL,   # e.g. 
  #class_colors  =  c(Topography = "#8c510a", Climate = "#92c5de"),
  predictor_labels = NULL    # e.g. 
  #predictor_labels =c(LocalRelief = "Local Relief", Ksn = "ksn")
)


# ==============================================================================
# 3. CONFIG CHECKS  (stops early with a clear message)
# ==============================================================================
all_predictors <- unlist(config$pool, use.names = FALSE)
transform_lists <- list(
  Positive_to_Log         = config$Positive_to_Log,
  Arcsin_Sqrt_Proportions = config$Arcsin_Sqrt_Proportions,
  Proportions_to_Logit    = config$Proportions_to_Logit
)
all_transformed <- unlist(transform_lists, use.names = FALSE)

no_transform <- setdiff(all_predictors, all_transformed)
if (length(no_transform) > 0)
  stop("In the pool but in no transformation list: ",
       paste(no_transform, collapse = ", "), call. = FALSE)

not_in_pool <- setdiff(all_transformed, all_predictors)
if (length(not_in_pool) > 0)
  stop("In a transformation list but not in the pool: ",
       paste(not_in_pool, collapse = ", "), call. = FALSE)

twice <- unique(all_transformed[duplicated(all_transformed)])
if (length(twice) > 0)
  stop("Predictor in more than one transformation list: ",
       paste(twice, collapse = ", "), call. = FALSE)


# --- Dataset labels must exist, be non-empty and be unique --------------------
ds_labels <- names(config$datasets)
if (is.null(ds_labels) || any(is.na(ds_labels)) || any(trimws(ds_labels) == ""))
  stop("Every entry in config$datasets needs a label, e.g.\n",
       '  datasets = c("Dataset 1" = "file1.xlsx", "Dataset 2" = "file2.csv")',
       call. = FALSE)
if (anyDuplicated(ds_labels))
  stop("Duplicated dataset labels: ",
       paste(unique(ds_labels[duplicated(ds_labels)]), collapse = ", "),
       call. = FALSE)


# ==============================================================================
# 4. LOAD DATA, then GENERATE MODELS
# ==============================================================================
read_table_any <- function(file) {
  path <- here("data", file)
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) readxl::read_excel(path)
  else if (ext == "csv")         read.csv(path, check.names = FALSE)
  else stop("Unsupported file type (use .xlsx, .xls or .csv): ", file, call. = FALSE)
}

datasets <- lapply(config$datasets, read_table_any)
names(datasets) <- names(config$datasets)

for (nm in names(datasets)) {
  miss <- setdiff(c(config$response_var, all_predictors), names(datasets[[nm]]))
  if (length(miss) > 0)
    stop("Dataset '", nm, "' lacks column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
}

Models_ensmbl <- generate_models_ensemble(config$pool,
                                          n_random_models = config$n_random_models,
                                          seed = config$seed)
message("Setup finished: ", length(datasets), " datasets, ",
        length(all_predictors), " predictors, ",
        length(config$pool), " classes.")
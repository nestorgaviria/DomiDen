# ==============================================================================
# DomiDen | 01_DA_functions.R
# Dominance analysis functions (class-level and predictor-level LMG).
#
# These functions are called by 04_run_and_plot.R. You do not need to edit this
# file unless you want to change the analysis itself.
#
# Exported functions:
#   - Mdls_dominance_Analysis_OLS_AIC()
#       Class-level dominance (groups = pool classes with 2+ predictors).
#   - Mdls_dominance_Analysis_OLS_AIC_predictor()
#       Predictor-level dominance (groups = NULL).
#   - robust_logit_transform()
#       Robust logit transform for proportions in [0,1] or [0,100].
#
# Inputs:
#   - df: data.frame with response and predictor columns.
#   - model_list: list of predictor sets from 02_model_generator.R.
#   - dataset_name: character label (from config$datasets).
#   - cfg: config list from 03_setup_config.R (response, pool, transforms, etc.).
#
# Outputs:
#   A list with:
#     $models            : screening results for all models (AIC, BIC, R2, etc.)
#     $winners           : top n_top_models by AIC
#     $dominance         : LMG results per model (class or predictor level)
#     $dominance_summary : aggregated LMG by (Dataset, Class) or
#                          (Dataset, Class, Predictor), with mean, SD, SE, 95% CI
#     $logit_transform_info, $arcsin_sqrt_transform_info : diagnostics
#
# VERSION 1.0 | Last update 01.10.2026
# ==============================================================================

Mdls_dominance_Analysis_OLS_AIC <- function(df, model_list, dataset_name, cfg) {
  
  # START TIMER
  start_time <- Sys.time()
  print(paste(">>> Processing:", dataset_name, "| N =", nrow(df)))
  print(paste("    Started at:", start_time))
  
  # =========================================================
  # 1. DATA TRANSFORMATION LAYER
  # =========================================================
  
  # A. Target Transformation
  y_col <- cfg$response_var
  if(cfg$response_transform == "log") {
    df$target_y <- log(df[[y_col]])
  } else {
    df$target_y <- df[[y_col]]
  }
  
  # B. Predictor Transformation (Log & Logit)
  log_vars <- cfg$Positive_to_Log
  logit_vars <- cfg$Proportions_to_Logit
  arcsin_sqrt_vars <- if (is.null(cfg$Arcsin_Sqrt_Proportions)) {
    character(0)
  } else {
    cfg$Arcsin_Sqrt_Proportions
  }
  
  overlap_log_logit <- intersect(log_vars, logit_vars)
  overlap_log_sqrt   <- intersect(log_vars, arcsin_sqrt_vars)
  overlap_logit_sqrt <- intersect(logit_vars, arcsin_sqrt_vars)
  
  if (length(overlap_log_logit) > 0 ||
      length(overlap_log_sqrt) > 0 ||
      length(overlap_logit_sqrt) > 0) {
    stop(
      "Predictors must belong to only one transformation group. Overlaps: ",
      paste(
        c(overlap_log_logit, overlap_log_sqrt, overlap_logit_sqrt),
        collapse = ", "
      )
    )
  }
  
  # Apply Log (with safety)
  for(col in log_vars) {
    if(col %in% names(df)) {
      vals <- df[[col]]
      if(any(vals <= 0, na.rm=TRUE)) {
        vals <- ifelse(vals <= 0, min(vals[vals>0], na.rm=TRUE)/2, vals)
      }
      df[[paste0("log_", col)]] <- log(vals)
    }
  }
  
  # =========================================================
  # Apply robust, data-driven logit transformations
  # =========================================================
  
  logit_transform_info <- list()
  
  for (col in logit_vars) {
    
    if (col %in% names(df)) {
      
      logit_out <- robust_logit_transform(
        x = df[[col]],
        variable_name = col
      )
      
      df[[paste0("logit_", col)]] <- logit_out$transformed
      
      # Store transformation metadata for reproducibility and diagnostics
      logit_transform_info[[col]] <- data.frame(
        Variable = col,
        Epsilon_Lower = logit_out$epsilon_lower,
        Epsilon_Upper = logit_out$epsilon_upper,
        N_Zero = logit_out$n_zero,
        N_One = logit_out$n_one,
        Converted_From_Percent = logit_out$converted_from_percent,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      warning(
        "Logit variable '", col,
        "' is listed in cfg$Proportions_to_Logit but is absent from df."
      )
    }
  }
  
  # Make an inspectable table; retained in the output object below.
  logit_transform_info <- if (length(logit_transform_info) > 0) {
    do.call(rbind, logit_transform_info)
  } else {
    data.frame(
      Variable = character(0),
      Epsilon_Lower = numeric(0),
      Epsilon_Upper = numeric(0),
      N_Zero = integer(0),
      N_One = integer(0),
      Converted_From_Percent = logical(0),
      stringsAsFactors = FALSE
    )
  }
  
  # =========================================================
  # Apply square-root transformations to zero-containing proportions
  # =========================================================
  
  arcsin_sqrt_transform_info <- list()
  
  for (col in arcsin_sqrt_vars) {
    
    if (col %in% names(df)) {
      
      vals <- as.numeric(df[[col]])
      vals_obs <- vals[!is.na(vals)]
      
      if (length(vals_obs) == 0) {
        stop("Arcsin Square-root transform failed for '", col, "': all values are NA.")
      }
      
      converted <- FALSE
      
      if (max(vals_obs) > 1) {
        if (max(vals_obs) <= 100 && min(vals_obs) >= 0) {
          vals <- vals / 100
          converted <- TRUE
        } else {
          stop(
            "Arcsin Square-root transform failed for '", col,
            "': expected values in [0,1] or [0,100]."
          )
        }
      }
      
      if (any(vals < 0 | vals > 1, na.rm = TRUE)) {
        stop(
          "Arcsin Square-root transform failed for '", col,
          "': values must lie in [0,1]."
        )
      }
      
      df[[paste0("arcsin_sqrt_", col)]] <- asin(sqrt(vals))

      arcsin_sqrt_transform_info[[col]] <- data.frame(
        Variable = col,
        N_Zero = sum(vals == 0, na.rm = TRUE),
        N_One = sum(vals == 1, na.rm = TRUE),
        Converted_From_Percent = converted,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      warning(
        "Square-root variable '", col,
        "' is listed in cfg$Arcsin_Sqrt_Proportions but is absent from df."
      )
      
    }
  }
  
  arcsin_sqrt_transform_info <- if (length(arcsin_sqrt_transform_info) > 0) {
    do.call(rbind, arcsin_sqrt_transform_info)
  } else {
    data.frame(
      Variable = character(0),
      N_Zero = integer(0),
      N_One = integer(0),
      Converted_From_Percent = logical(0),
      stringsAsFactors = FALSE
    )
  }
  
  df <- as.data.frame(df)
  
  # =========================================================
  # 2. PHASE 1: SCREENING LOOP (OLS + AIC)  
  # =========================================================
  
  print(paste("    >>> PHASE 1: Screening", length(model_list), "models using OLS + AIC..."))
  
  df_screening_results <- data.frame()
  n_failed <- 0  # Track failures
  
  # Loop through ALL models to get metrics
  for (m_name in names(model_list)) {
    tryCatch({
      raw_preds <- model_list[[m_name]]
      final_preds <- c()
      
      # Name Mapping Logic
      for (p in raw_preds) {
        if (p %in% log_vars && paste0("log_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("log_", p))
        } else if (p %in% logit_vars && paste0("logit_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("logit_", p))
        } else if (p %in% arcsin_sqrt_vars && paste0("arcsin_sqrt_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("arcsin_sqrt_", p))
        } else {
          final_preds <- c(final_preds, p)
        }
      }
      
      # Fit OLS regression
      model_data <- df[, c("target_y", final_preds), drop = FALSE]
      model_data <- na.omit(model_data)
      
      # Check if we have enough data after NA removal
      if(nrow(model_data) < length(final_preds) + 2) {
        warning(paste("Model", m_name, "skipped: insufficient data after NA removal"))
        n_failed <- n_failed + 1
        next
      }
      
      # Scale predictors
      model_data[, final_preds] <- scale(model_data[, final_preds])
      
      # Fit OLS
      formula_str <- paste("target_y ~", paste(final_preds, collapse = " + "))
      ols_model <- lm(as.formula(formula_str), data = model_data)
      
      res_row <- data.frame(
        Dataset     = dataset_name,
        Model_ID    = m_name,
        AIC         = AIC(ols_model),
        BIC         = BIC(ols_model),
        R2_adj      = summary(ols_model)$adj.r.squared,
        RMSE        = sqrt(mean(residuals(ols_model)^2)),
        LogLik      = as.numeric(logLik(ols_model)),
        N_obs       = nrow(model_data),  # Track sample size
        N_preds     = length(final_preds)
      )
      
      df_screening_results <- rbind(df_screening_results, res_row)
      
    }, error = function(e) {
      warning(paste("Model", m_name, "failed:", e$message))
      n_failed <- n_failed + 1
    })
  }
  
  print(paste("    >>> Screening complete.", nrow(df_screening_results), 
              "models succeeded,", n_failed, "failed."))
  
  # =========================================================
  # 3. SELECTION: PICK THE WINNERS  
  # =========================================================
  
  n_select <- ifelse(is.null(cfg$n_top_models), 50, cfg$n_top_models)
  
  df_winners <- df_screening_results %>%
    arrange(AIC) %>%  # Lower AIC = better
    slice_head(n = n_select)
  
  print(paste("    >>> SELECTION: Selected top", nrow(df_winners), "models based on AIC"))
  
  # =========================================================
  #  4. PHASE 2: DOMINANCE ANALYSIS ON WINNERS
  # =========================================================
  print("    >>> PHASE 2: Running Dominance Analysis on Winners...")
  
  df_dominance_results <- data.frame()
  
  for(i in 1:nrow(df_winners)) {
    
    m_name <- df_winners$Model_ID[i]  
    raw_preds <- model_list[[m_name]]
    
    # Apply transformations
    final_preds <- c()
    
    for (p in raw_preds) {
      if (p %in% log_vars) {
        final_preds <- c(final_preds, paste0("log_", p))
      } else if (p %in% logit_vars) {
        final_preds <- c(final_preds, paste0("logit_", p))
      } else if (p %in% arcsin_sqrt_vars) {
        final_preds <- c(final_preds, paste0("arcsin_sqrt_", p))
      } else {
        final_preds <- c(final_preds, p)
      }
    }
    
    # Prepare data for this specific model
    model_data <- df[, c("target_y", final_preds), drop = FALSE]
    model_data <- na.omit(model_data)
    
    # Scale predictors (critical for fair dominance comparison)
    model_data[, final_preds] <- scale(model_data[, final_preds])
    
    # Fit OLS model for dominance analysis
    formula_str <- paste("target_y ~", paste(final_preds, collapse = " + "))
    ols_model <- lm(as.formula(formula_str), data = model_data)
    
    # Create groups for THIS specific model AND mapping for singles
    model_groups <- list()
    var_to_class_map <- list()  # Track which class each variable belongs to
    
    for(class_name in names(cfg$pool)) {
      class_vars_raw <- cfg$pool[[class_name]]
      
      # Find which transformed predictors belong to this class
      class_vars_transformed <- c()
      for(var in class_vars_raw) {
        transformed_var <- NULL
        
        if (paste0("log_", var) %in% final_preds) {
          transformed_var <- paste0("log_", var)
        } else if (paste0("logit_", var) %in% final_preds) {
          transformed_var <- paste0("logit_", var)
        } else if (paste0("arcsin_sqrt_", var) %in% final_preds) {
          transformed_var <- paste0("arcsin_sqrt_", var)
        } else if (var %in% final_preds) {
          transformed_var <- var
        }
        
        if(!is.null(transformed_var)) {
          class_vars_transformed <- c(class_vars_transformed, transformed_var)
          var_to_class_map[[transformed_var]] <- class_name  #Store mapping
        }
      }
      
      # Only include as group if it has 2+ predictors
      if(length(class_vars_transformed) >= 2) {
        model_groups[[class_name]] <- class_vars_transformed
      }
    }
    
    # Run dominance analysis
    tryCatch({
      lmg_result <- relaimpo::calc.relimp(
        ols_model, 
        type = "lmg",
        rela = FALSE,
        groups = if(length(model_groups) > 0) model_groups else NULL,
        groupnames = if(length(model_groups) > 0) names(model_groups) else NULL
      )
      
      # Extract results
      result_names <- names(lmg_result$lmg)
      result_values <- as.numeric(lmg_result$lmg)
      
      # Process each result
      for(j in seq_along(result_names)) {
        result_name <- result_names[j]
        result_value <- result_values[j]
        
        # Determine if this is a group or individual predictor
        if(result_name %in% names(model_groups)) {
          # It's a group
          class_name <- result_name
          n_in_class <- length(model_groups[[result_name]])
        } else {
          # It's an individual predictor - map back to its class
          class_name <- var_to_class_map[[result_name]]
          if(is.null(class_name)) {
            # Fallback: couldn't map (shouldn't happen)
            class_name <- result_name
            n_in_class <- 0
          } else {
            n_in_class <- 1
          }
        }
        
        # Store result
        group_importance <- data.frame(
          Model_ID = m_name,
          Dataset = dataset_name,
          Class = class_name,
          LMG_Importance = result_value * 100,
          Model_R2 = summary(ols_model)$r.squared,
          N_Predictors_Total = length(final_preds),
          N_Predictors_InClass = n_in_class,
          Is_Group = n_in_class >= 2,  # NEW: Flag for grouped vs individual
          stringsAsFactors = FALSE
        )
        
        df_dominance_results <- rbind(df_dominance_results, group_importance)
      }
      
    }, error = function(e) {
      warning(paste("Dominance analysis failed for model", m_name, ":", e$message))
    })
  }
  
  # =========================================================
  # 6. FINAL REPORTING
  # =========================================================
  
  end_time <- Sys.time()
  duration <- round(difftime(end_time, start_time, units = "mins"), 2)
  print(paste("    Finished at:", end_time))
  print(paste("    Total Execution Time:", duration, "minutes"))
  
  # Return results
  return(list(
    models = df_screening_results,
    winners = df_winners,
    dominance = df_dominance_results,
    dominance_summary = df_dominance_results %>%
      group_by(Dataset, Class) %>%
      summarise(
        Mean_Importance   = mean(LMG_Importance),
        SD_Importance     = sd(LMG_Importance),
        SE_Importance     = sd(LMG_Importance) / sqrt(n()),
        N_Models          = n(),
        Median_Importance = median(LMG_Importance),
        Min_Importance    = min(LMG_Importance),
        Max_Importance    = max(LMG_Importance),
        N_AsGroup         = sum(Is_Group),
        N_AsSingle        = sum(!Is_Group),
        .groups = "drop"
      ) %>%
      mutate(
        CI_lower = ifelse(N_Models > 1,
                          Mean_Importance - qt(0.975, N_Models - 1) * SE_Importance,
                          NA_real_),
        CI_upper = ifelse(N_Models > 1,
                          Mean_Importance + qt(0.975, N_Models - 1) * SE_Importance,
                          NA_real_)
      ) %>%
      # Reorder columns so CI comes right after SE
      select(Dataset, Class, Mean_Importance, SD_Importance, SE_Importance,
             CI_lower, CI_upper, Median_Importance, Min_Importance, Max_Importance,
             N_Models, N_AsGroup, N_AsSingle),
    logit_transform_info = logit_transform_info,
    arcsin_sqrt_transform_info = arcsin_sqrt_transform_info
    
    
  ))
}

# ≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
# *******************************************************************************
# *********   Function to run predictor-level dominance analysis   **************
# *********   on top AIC-selected models from an ensemble          **************
# *********   (predictor-level LMG via relaimpo)                   **************
# ≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
# version 08.09.2026
#
# This function:
# - Transforms response and predictors as per cfg
# - Screens all models in model_list using OLS + AIC
# - Selects top n_top_models (or default 50) by AIC
# - Runs predictor-level LMG dominance (groups = NULL) on winners
# - Returns one row per predictor per model in $dominance
# - Summarizes by (Dataset, Class, Predictor) in $dominance_summary



Mdls_dominance_Analysis_OLS_AIC_predictor <- function(df, model_list, dataset_name, cfg) {
  
  # START TIMER
  start_time <- Sys.time()
  print(paste(">>> Processing (predictor-level):", dataset_name, "| N =", nrow(df)))
  print(paste("    Started at:", start_time))
  
  # =========================================================
  # 1. DATA TRANSFORMATION LAYER
  # =========================================================
  
  # A. Target Transformation
  y_col <- cfg$response_var
  if (cfg$response_transform == "log") {
    df$target_y <- log(df[[y_col]])
  } else {
    df$target_y <- df[[y_col]]
  }
  
  # B. Predictor Transformation (Log & Logit)
  log_vars <- cfg$Positive_to_Log
  logit_vars <- cfg$Proportions_to_Logit
  arcsin_sqrt_vars <- if (is.null(cfg$Arcsin_Sqrt_Proportions)) {
    character(0)
  } else {
    cfg$Arcsin_Sqrt_Proportions
  }

  overlap_log_logit <- intersect(log_vars, logit_vars)
  overlap_log_sqrt   <- intersect(log_vars, arcsin_sqrt_vars)
  overlap_logit_sqrt <- intersect(logit_vars, arcsin_sqrt_vars)
  
  if (length(overlap_log_logit) > 0 ||
      length(overlap_log_sqrt) > 0 ||
      length(overlap_logit_sqrt) > 0) {
    stop(
      "Predictors must belong to only one transformation group. Overlaps: ",
      paste(
        c(overlap_log_logit, overlap_log_sqrt, overlap_logit_sqrt),
        collapse = ", "
      )
    )
  }
  
  # Apply Log (with safety)
  for (col in log_vars) {
    if (col %in% names(df)) {
      vals <- df[[col]]
      if (any(vals <= 0, na.rm = TRUE)) {
        vals <- ifelse(vals <= 0, min(vals[vals > 0], na.rm = TRUE) / 2, vals)
      }
      df[[paste0("log_", col)]] <- log(vals)
    }
  }
  
  # =========================================================
  # Apply robust, data-driven logit transformations
  # =========================================================
  
  logit_transform_info <- list()
  
  for (col in logit_vars) {
    
    if (col %in% names(df)) {
      
      logit_out <- robust_logit_transform(
        x = df[[col]],
        variable_name = col
      )
      
      df[[paste0("logit_", col)]] <- logit_out$transformed
      
      # Store transformation metadata for reproducibility and diagnostics
      logit_transform_info[[col]] <- data.frame(
        Variable = col,
        Epsilon_Lower = logit_out$epsilon_lower,
        Epsilon_Upper = logit_out$epsilon_upper,
        N_Zero = logit_out$n_zero,
        N_One = logit_out$n_one,
        Converted_From_Percent = logit_out$converted_from_percent,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      warning(
        "Logit variable '", col,
        "' is listed in cfg$Proportions_to_Logit but is absent from df."
      )
      
    }
  }
  
  # Make an inspectable table; retained in the output object below.
  logit_transform_info <- if (length(logit_transform_info) > 0) {
    do.call(rbind, logit_transform_info)
  } else {
    data.frame(
      Variable = character(0),
      Epsilon_Lower = numeric(0),
      Epsilon_Upper = numeric(0),
      N_Zero = integer(0),
      N_One = integer(0),
      Converted_From_Percent = logical(0),
      stringsAsFactors = FALSE
    )
  }
  
  # =========================================================
  # Apply square-root transformations to zero-containing proportions
  # =========================================================
  
  arcsin_sqrt_transform_info <- list()
  
  for (col in arcsin_sqrt_vars) {
    
    if (col %in% names(df)) {
      
      vals <- as.numeric(df[[col]])
      vals_obs <- vals[!is.na(vals)]
      
      if (length(vals_obs) == 0) {
        stop("Arcsin Square-root transform failed for '", col, "': all values are NA.")
      }
      
      converted <- FALSE
      
      if (max(vals_obs) > 1) {
        if (max(vals_obs) <= 100 && min(vals_obs) >= 0) {
          vals <- vals / 100
          converted <- TRUE
        } else {
          stop(
            "Arcsin Square-root transform failed for '", col,
            "': expected values in [0,1] or [0,100]."
          )
        }
      }
      
      if (any(vals < 0 | vals > 1, na.rm = TRUE)) {
        stop(
          "Arcsin Square-root transform failed for '", col,
          "': values must lie in [0,1]."
        )
      }
      
      df[[paste0("arcsin_sqrt_", col)]] <- asin(sqrt(vals))
      
      arcsin_sqrt_transform_info[[col]] <- data.frame(
        Variable = col,
        N_Zero = sum(vals == 0, na.rm = TRUE),
        N_One = sum(vals == 1, na.rm = TRUE),
        Converted_From_Percent = converted,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      warning(
        "Square-root variable '", col,
        "' is listed in cfg$Arcsin_Sqrt_Proportions but is absent from df."
      )
      
    }
  }
  
  arcsin_sqrt_transform_info <- if (length(arcsin_sqrt_transform_info) > 0) {
    do.call(rbind, arcsin_sqrt_transform_info)
  } else {
    data.frame(
      Variable = character(0),
      N_Zero = integer(0),
      N_One = integer(0),
      Converted_From_Percent = logical(0),
      stringsAsFactors = FALSE
    )
  }
  
  df <- as.data.frame(df)
  
  # =========================================================
  # 2. PHASE 1: SCREENING LOOP (OLS + AIC)
  # =========================================================
  
  print(paste("    >>> PHASE 1: Screening", length(model_list), "models using OLS + AIC..."))
  
  df_screening_results <- data.frame()
  n_failed <- 0
  
  for (m_name in names(model_list)) {
    tryCatch({
      raw_preds <- model_list[[m_name]]
      final_preds <- c()
      
      # Name Mapping Logic
      for (p in raw_preds) {
        if (p %in% log_vars && paste0("log_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("log_", p))
        } else if (p %in% logit_vars && paste0("logit_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("logit_", p))
        } else if (p %in% arcsin_sqrt_vars && paste0("arcsin_sqrt_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("arcsin_sqrt_", p))
        } else {
          final_preds <- c(final_preds, p)
        }
      }
      
      # Fit OLS regression
      model_data <- df[, c("target_y", final_preds), drop = FALSE]
      model_data <- na.omit(model_data)
      
      # Check if we have enough data after NA removal
      if (nrow(model_data) < length(final_preds) + 2) {
        warning(paste("Model", m_name, "skipped: insufficient data after NA removal"))
        n_failed <- n_failed + 1
        next
      }
      
      # Scale predictors
      model_data[, final_preds] <- scale(model_data[, final_preds])
      
      # Fit OLS
      formula_str <- paste("target_y ~", paste(final_preds, collapse = " + "))
      ols_model <- lm(as.formula(formula_str), data = model_data)
      
      res_row <- data.frame(
        Dataset     = dataset_name,
        Model_ID    = m_name,
        AIC         = AIC(ols_model),
        BIC         = BIC(ols_model),
        R2_adj      = summary(ols_model)$adj.r.squared,
        RMSE        = sqrt(mean(residuals(ols_model)^2)),
        LogLik      = as.numeric(logLik(ols_model)),
        N_obs       = nrow(model_data),
        N_preds     = length(final_preds),
        stringsAsFactors = FALSE
      )
      
      df_screening_results <- rbind(df_screening_results, res_row)
      
    }, error = function(e) {
      warning(paste("Model", m_name, "failed:", e$message))
      n_failed <- n_failed + 1
    })
  }
  
  print(paste("    >>> Screening complete.", nrow(df_screening_results),
              "models succeeded,", n_failed, "failed."))
  
  
  # =========================================================
  # 3. SELECTION: PICK THE WINNERS
  # =========================================================
  
  n_select <- ifelse(is.null(cfg$n_top_models), 50, cfg$n_top_models)
  
  df_winners <- df_screening_results %>%
    arrange(AIC) %>%
    slice_head(n = n_select)
  
  print(paste("    >>> SELECTION: Selected top", nrow(df_winners), "models based on AIC"))
  
  
  # =========================================================
  # 4. PHASE 2: PREDICTOR-LEVEL DOMINANCE ANALYSIS ON WINNERS
  # =========================================================
  print("    >>> PHASE 2: Running predictor-level Dominance Analysis on Winners...")
  
  df_dominance_results <- data.frame()
  
  for (i in 1:nrow(df_winners)) {
    
    m_name <- df_winners$Model_ID[i]
    raw_preds <- model_list[[m_name]]
    
    # Apply transformations
    final_preds <- c()
    
    for (p in raw_preds) {
      if (p %in% log_vars) {
        final_preds <- c(final_preds, paste0("log_", p))
      } else if (p %in% logit_vars) {
        final_preds <- c(final_preds, paste0("logit_", p))
      } else if (p %in% arcsin_sqrt_vars) {
        final_preds <- c(final_preds, paste0("arcsin_sqrt_", p))
      } else {
        final_preds <- c(final_preds, p)
      }
    }
    
    # Prepare data for this specific model
    model_data <- df[, c("target_y", final_preds), drop = FALSE]
    model_data <- na.omit(model_data)
    
    # Scale predictors (critical for fair dominance comparison)
    model_data[, final_preds] <- scale(model_data[, final_preds])
    
    # Fit OLS model for dominance analysis
    formula_str <- paste("target_y ~", paste(final_preds, collapse = " + "))
    ols_model <- lm(as.formula(formula_str), data = model_data)
    
    # Build mapping from transformed predictor -> class
    var_to_class_map <- list()
    
    for (class_name in names(cfg$pool)) {
      class_vars_raw <- cfg$pool[[class_name]]
      
      for (var in class_vars_raw) {
        transformed_var <- NULL
        
        if (paste0("log_", var) %in% final_preds) {
          transformed_var <- paste0("log_", var)
        } else if (paste0("logit_", var) %in% final_preds) {
          transformed_var <- paste0("logit_", var)
        } else if (paste0("arcsin_sqrt_", var) %in% final_preds) {
          transformed_var <- paste0("arcsin_sqrt_", var)
        } else if (var %in% final_preds) {
          transformed_var <- var
        }
        
        if (!is.null(transformed_var)) {
          var_to_class_map[[transformed_var]] <- class_name
        }
      }
    }
    
    # Run predictor-level dominance analysis (NO grouping)
    tryCatch({
      lmg_result <- relaimpo::calc.relimp(
        ols_model,
        type = "lmg",
        rela = FALSE,
        groups = NULL,        # always predictor-level
        groupnames = NULL
      )
      
      # Extract predictor-level results
      pred_names <- names(lmg_result$lmg)
      pred_values <- as.numeric(lmg_result$lmg)
      
      for (j in seq_along(pred_names)) {
        pred_name <- pred_names[j]
        pred_value <- pred_values[j]
        
        # Map predictor to class
        class_name <- var_to_class_map[[pred_name]]
        if (is.null(class_name)) {
          class_name <- NA_character_
        }
        
        pred_row <- data.frame(
          Model_ID           = m_name,
          Dataset            = dataset_name,
          Predictor          = pred_name,
          Class              = class_name,
          LMG_Importance     = pred_value * 100,
          Model_R2           = summary(ols_model)$r.squared,
          N_Predictors_Total = length(final_preds),
          stringsAsFactors   = FALSE
        )
        
        df_dominance_results <- rbind(df_dominance_results, pred_row)
      }
      
    }, error = function(e) {
      warning(paste("Dominance analysis failed for model", m_name, ":", e$message))
    })
  }
  
  
  # =========================================================
  # 5. FINAL REPORTING
  # =========================================================
  
  end_time <- Sys.time()
  duration <- round(difftime(end_time, start_time, units = "mins"), 2)
  print(paste("    Finished at:", end_time))
  print(paste("    Total Execution Time:", duration, "minutes"))
  
  # Summarize by (Dataset, Class, Predictor)
  dominance_summary <- df_dominance_results %>%
    group_by(Dataset, Class, Predictor) %>%
    summarise(
      Mean_Importance   = mean(LMG_Importance),
      SD_Importance     = sd(LMG_Importance),
      SE_Importance     = sd(LMG_Importance) / sqrt(n()),
      N_Models          = n(),
      Median_Importance = median(LMG_Importance),
      Min_Importance    = min(LMG_Importance),
      Max_Importance    = max(LMG_Importance),
      .groups = "drop"
    ) %>%
    mutate(
      CI_lower = ifelse(N_Models > 1,
                        Mean_Importance - qt(0.975, N_Models - 1) * SE_Importance,
                        NA_real_),
      CI_upper = ifelse(N_Models > 1,
                        Mean_Importance + qt(0.975, N_Models - 1) * SE_Importance,
                        NA_real_)
    ) %>%
    dplyr::select(Dataset, Class, Predictor, Mean_Importance, SD_Importance, SE_Importance,
           CI_lower, CI_upper, Median_Importance, Min_Importance, Max_Importance, N_Models)
  
  # Return results
  return(list(
    models            = df_screening_results,
    winners           = df_winners,
    dominance         = df_dominance_results,
    dominance_summary = dominance_summary,
    logit_transform_info = logit_transform_info,
    arcsin_sqrt_transform_info = arcsin_sqrt_transform_info
    
    
  ))
}

# =========================================================
# Robust logit transformation for bounded proportions
# =========================================================
# - Accepts proportions in [0, 1] or percentages in [0, 100]
# - Converts percentages to proportions only if values exceed 1
# - Replaces exact 0 with half the smallest observed positive value
# - Replaces exact 1 with 1 - half the smallest distance below 1
# - Keeps observed interior values unchanged
# - Stops clearly if values lie outside valid bounds
#
# Returns:
#   list(
#     transformed = numeric vector,
#     epsilon_lower = numeric,
#     epsilon_upper = numeric,
#     converted_from_percent = logical,
#     n_zero = integer,
#     n_one = integer
#   )

robust_logit_transform <- function(x, variable_name = deparse(substitute(x))) {
  x <- as.numeric(x)
  x_obs <- x[!is.na(x)] # Isolate non-missing for bounds/checks
  
  if (length(x_obs) == 0) stop("Logit transform failed for '", variable_name, "': all values are NA.")
  
  rng <- range(x_obs)
  converted <- FALSE
  
  # Auto-detect percentages and scale to proportions
  if (rng[2] > 1) {
    if (rng[2] <= 100 && rng[1] >= 0) {
      x <- x / 100; x_obs <- x_obs / 100; rng <- rng / 100 # Keep all synced
      converted <- TRUE
    } else {
      stop("Failed for '", variable_name, "': must be in [0,1] or [0,100]. Range: [", rng[1], ", ", rng[2], "].")
    }
  }
  
  if (rng[1] < 0 || rng[2] > 1) stop("Failed for '", variable_name, "': values outside [0, 1].")
  
  n_zero <- sum(x_obs == 0)
  n_one  <- sum(x_obs == 1)
  eps_lower <- 0
  eps_upper <- 0
  
  # Adjust exact zeros (using which() safely ignores NAs)
  if (n_zero > 0) {
    if (n_zero == length(x_obs)) stop("Failed for '", variable_name, "': all values are 0 (no variation).")
    eps_lower <- min(x_obs[x_obs > 0]) / 2
    x[which(x == 0)] <- eps_lower
  }
  
  # Adjust exact ones
  if (n_one > 0) {
    if (n_one == length(x_obs)) stop("Failed for '", variable_name, "': all values are 1 (no variation).")
    eps_upper <- min(1 - x_obs[x_obs < 1]) / 2
    x[which(x == 1)] <- 1 - eps_upper
  }
  
  # Return (qlogis is R's built-in logit transform: log(x / (1 - x)))
  list(
    transformed = qlogis(x), 
    epsilon_lower = eps_lower,
    epsilon_upper = eps_upper,
    converted_from_percent = converted,
    n_zero = n_zero,
    n_one = n_one
  )
}

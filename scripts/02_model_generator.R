# ==============================================================================
# DomiDen | 02_model_generator.R
# Random model ensemble generator.
#
# This function is called by 03_setup_config.R to create a single ensemble of
# predictor combinations that is then used for all datasets in 04_run_and_plot.R.
# You do not need to edit this file unless you want to change how models are built.
#
# Exported function:
#   - generate_models_ensemble(pool, n_random_models, seed)
#
# Inputs:
#   - pool: named list from config$pool (class -> predictor names).
#   - n_random_models: number of random models to generate.
#   - seed: random seed for reproducibility.
#
# Behaviour:
#   - Generates models with 2–6 predictors.
#   - Predictors can come from 1 or multiple classes; not all classes must appear.
#   - No limit on the number of predictors per class within a model.
#   - Candidate predictors are weighted by inverse squared usage frequency, so
#     frequently used predictors become less likely over time.
#   - Duplicate predictor combinations are discarded; only unique models are kept.
#
# Outputs:
#   A named list: each element is a character vector of predictor names for one model.
#
# VERSION 1.0 | Last update 01.10.2026
# ==============================================================================

generate_models_ensemble <- function(pool, n_random_models, seed) {
  set.seed(seed)
  model_list <- list()
  
  # Track unique model signatures
  used_signatures <- character(0)
  
  # Helper function to create signature
  make_signature <- function(vars) {
    paste(sort(vars), collapse = "|")
  }
  
  # Helper function to add model (only if unique)
  add_model <- function(name, vars) {
    sig <- make_signature(vars)
    if(!(sig %in% used_signatures)) {
      model_list[[name]] <<- vars
      used_signatures <<- c(used_signatures, sig)
      return(TRUE)
    }
    # ROLLBACK: undo the counts accumulated during failed construction
    global_counts[vars] <<- global_counts[vars] - 1
    return(FALSE)
  }
  
  # 1. SETUP
  all_preds <- unlist(pool, use.names = FALSE)
  pred_class_map <- data.frame(
    predictor = all_preds,
    class = rep(names(pool), lengths(pool)),
    stringsAsFactors = FALSE
  )
  
  # Global Usage Tracker
  global_counts <- setNames(rep(0, length(all_preds)), all_preds)
  class_names <- names(pool)

  print(paste(" Generating", n_random_models, "Interaction Models..."))
  
  # =========================================================
  # STEP 1: GENERATE BALANCED RANDOM MODELS (Interactions)
  # =========================================================
  
  successful_models <- 0
  max_attempts <- n_random_models * 5  # Allow more attempts
  attempt <- 0
  
  while(successful_models < n_random_models && attempt < max_attempts) {
    attempt <- attempt + 1
    
    # A. Determine Size (Range 2 to 6)
    target_size <- sample(2:6, 1, prob = c(0.15, 0.2, 0.3, 0.2, 0.15))
    #target_size <- sample(1:3, 1, prob = c(0.35, 0.5, 0.15))
    
    
    current_model <- c()
    current_classes <- c()
    
    class_names <- names(pool)
    
    # B. Fill remaining slots (if target_size = 6, we need 1 more)
    while(length(current_model) < target_size) {
      
      candidates <- setdiff(all_preds, current_model)
      
      if(length(candidates) == 0) break
      
      # Weights
      cand_usage <- global_counts[candidates]
      weights <- 1 / (1 + cand_usage)^2 
      
      selected <- sample(candidates, 1, prob = weights)
      current_model <- c(current_model, selected)
      
      selected_cls <- pred_class_map$class[match(selected, pred_class_map$predictor)]
      current_classes <- c(current_classes, selected_cls)
      global_counts[selected] <- global_counts[selected] + 1
    }
    
    # Try to add model (only if unique)
    if(add_model(paste0("Mdl_", successful_models + 1), current_model)) {
      successful_models <- successful_models + 1
    }
  }
  # 
  if(successful_models < n_random_models) {
    warning(paste("Only generated", successful_models, "unique models out of",
                  n_random_models, "requested"))
  }
  
  # =========================================================
  # DIAGNOSTICS
  # =========================================================
  usage_stats <- as.data.frame(global_counts)
  colnames(usage_stats) <- "Times_Selected"
  
  print("--- Balance Check (Top 5 & Bottom 5) ---")
  print(head(usage_stats[order(usage_stats$Times_Selected), , drop=FALSE], 5))
  print(tail(usage_stats[order(usage_stats$Times_Selected), , drop=FALSE], 5))
  
  print(paste("\n>>> Total unique models generated:", length(model_list)))
  
  return(model_list)
}

require(ggplot2)
require(dplyr)
require(tidyr)
library(digest)
library(cmdstanr)
library(posterior)
require(rethinking)

# Change these paths to your Rtools gcc and g++ executables:
Sys.setenv(
  CC = "C:/rtools45/x86_64-w64-mingw32.static.posix/bin/gcc.exe",
  CXX = "C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe"
)


load_and_preprocess_rMT_data <- function (fn = './data/motor_thresholds.csv'){
  x <- read.csv(fn)
  xx <- x %>% dplyr::select(-all_of(c('education_years', 'age', 'Sex', 'MT_fig8_Lt', 'MT_H7_Lt', 'MT_rf_Lt', 'MT_fig8_Rt', 'MT_H7_Rt', 'MT_rf_Rt')))

  xx.long <- xx %>% pivot_longer(cols = !c(pcode),
                                names_to = c('coil','side'),
                                names_pattern = '^MT_(.*?)_(.*?)_V$',
                                values_to = 'value')

  xx.long$coil[xx.long$coil == 'rf'] <- 'RF'
  xx.long$coil[xx.long$coil == 'h7'] <- 'H7'
  xx.long$side[xx.long$side == 'lt'] <- 'Left'
  xx.long$side[xx.long$side == 'rt'] <- 'Right'

  return(xx.long)
}

Vmax <- function(coil){
  return(ifelse(coil=='fig8',1670,1700))
}
adjust_intensity <- function (coil, machine_setting){
  return(machine_setting * Vmax(coil) / Vmax('RF'))
}



load_and_preprocess_MEP_data <- function(fn = 'data/cleaned_MEP_data.csv', 
                                         possible_percent = seq(from = 60, to = 150, by = 10),
                                         filter_percent = NULL,
                                         prepare_for_modelling = TRUE,
                                         intensity_varname = 'adjusted_MT_cat'){
  emg <- read.csv(fn)

  # Filter instances were data was invalid: subject 303, subject 301, fig8 coil, right side
  p2p <- emg %>%
    filter(toUse == 1 & patient_code != 303) %>% filter (!(patient_code == 301 & coil == 'fig8' & side == 'rt')) %>%
    dplyr::select(patient_code, coil, side, block, repetition, machine_setting, percent_threshold, PeakToPeak)

                                          
  p2p <- adjust_rMT_for_all(p2p)


  p2p$adjusted_MT_cat <- possible_percent[argmin_edges(p2p$adjusted_MT, possible_percent)]
  p2p$percent_threshold_cat <- possible_percent[argmin_edges(p2p$percent_threshold, possible_percent)]
                                        
  # adjust machine setting
  p2p <- p2p %>% mutate(raw_Intensity = adjust_intensity(coil, machine_setting)) %>%
                 mutate(Intensity = Vmax(coil)*raw_Intensity/(Vmax('RF'))) %>%
                 mutate(Vmax = 100*Vmax(coil)/Vmax('RF')) 
  

  if (prepare_for_modelling){  
    p2p <- prepare_realworld_data_for_modelling(p2p, intensity_varname =  intensity_varname)
    if (!is.null(filter_percent)){
      p2p <- p2p %>% filter(Irel >= filter_percent[1] & Irel <= filter_percent[2] )
    }
    }
  return(p2p)
}

prepare_realworld_data_for_modelling <- function(p2p, keep_cols = c('subj','coil','side','Intensity','Y','Vmax'), intensity_varname = 'adjusted_MT_cat'){
  data <- p2p %>% rename(
    subj = patient_code,
    Intensity = !!rlang::sym(intensity_varname),
    Y = PeakToPeak
  ) %>% select(all_of(keep_cols))

  
  return(data)
}

adjust_rMT <- function(p2p, fixed_percent = 110)
{
  #Find the closest MSO threshold in RF corresponding to the MSO if H7 in the fixed MT:
  sides <- unique(p2p$side)
  subjects <- unique(p2p$patient_code)
  coils <- unique(p2p$coil)
  p2p$adjusted_MT <- p2p$percent_threshold
  p2p$adjusted_MT_min_dist <- 0

  save_h7_mso <- list()

  for (subj in subjects){
    for (sid in sides){
      cat(sprintf("%s, %s:",subj, sid))
      #find the H7 MSO corresponding to fixed percent:
      
      h7_mso_df <- p2p %>%
                  filter(patient_code == subj & coil == 'h7' & side == sid &
                            percent_threshold == fixed_percent)
      h7_mso <- h7_mso_df$machine_setting
      if (is.null(h7_mso)){
        next
      }
      if (length(h7_mso) > 1){
        h7_mso_val <- min(h7_mso)
      }
      else {
        h7_mso_val <- h7_mso
      }
      
      #extract MSOs and percent_threshold for RF:
      rf_mso <- p2p %>% filter(patient_code == subj & coil == 'rf' & side == sid) 
      
      #Compute the distance between the observed RF machine_setting and the H7 machine_setting which correspond to fixed percent_thrshold
      #This will find the closest corresponding MSO
      dist <- sqrt((rf_mso$machine_setting - h7_mso_val)^2)
      #Find the closest percent threshold:
      
      argmin <- (dist == min(dist))
      #And retrieve the corresponding MSO, if there are more than one take the minimal on
      rf_mso_val <- unique(rf_mso$machine_setting[argmin])
      
      print("Closest RF MSOs:")
      print(rf_mso[argmin,]$machine_setting)
     
      
      #Update an adjusted column:
      subset <- (p2p$patient_code == subj & p2p$coil == 'rf' & p2p$side == sid & p2p$machine_setting %in% rf_mso_val)
      subset_orig_fixed <- (p2p$patient_code == subj & p2p$coil == 'rf' & p2p$side == sid & p2p$adjusted_MT == fixed_percent)
      if (sum(subset_orig_fixed)>0) {p2p[subset_orig_fixed,]$adjusted_MT <- NA}
      p2p[subset,]$adjusted_MT <- fixed_percent
      original_MT <- p2p[subset,]$percent_threshold[1]
      
      #log the minimal distace. If the distance is too large these records should be excluded:
      p2p[subset,]$adjusted_MT_min_dist <- min(dist)
      
      
      
      cat(sprintf("rMT=%d, H7 MSO =  %d, original RF = %d ,picked RF MSO = %d, Min distance= %d\n", fixed_percent, h7_mso_val, original_MT, rf_mso_val, min(dist)))
      
    }
  }

  p2p$coil[p2p$coil == 'rf'] <- 'RF'
  p2p$coil[p2p$coil == 'h7'] <- 'H7'
  p2p$side[p2p$side == 'lt'] <- 'Left'
  p2p$side[p2p$side == 'rt'] <- 'Right'

  return(p2p)
}

adjust_rMT_for_all <- function(p2p)
{
  #For each stimulus intensity in RF find the closest corresponding adjusted rMT based on H7 coil
  
  compute_dist_matrix <- function(rf_tab, h7_tab){
    calc_dist <- function (x,y){
      return (sqrt((x - y)^2))
    }
      
    dist <- dist_matrix <- outer(rf_tab$stim_output, h7_tab$stim_output, calc_dist)
    colnames(dist) <- h7_tab$percent_threshold
    rownames(dist) <- rf_tab$percent_threshold

    return (dist)
  }

  get_adjusted_rMT_by_dist <- function(original_rf_percent, dist_matrix){
    #Find the row corresponding to original_rf_percent
    row_index <- which(rownames(dist_matrix) == original_rf_percent)
    if (length(row_index) == 0){
      return (NA)
    }
    #Find the column index of the minimal distance in that row
    col_index <- which(dist_matrix[row_index, ] == min(dist_matrix[row_index, ]))
    # if (length(col_index) > 1){
      # compute the weighted mean of the corresponding percent thresholds
    corresponding_percents <- as.numeric(colnames(dist_matrix)[col_index]) 
    adjusted_rMT <- round(mean(corresponding_percents))
    return (list(adjusted_rMT = adjusted_rMT, dist = as.vector(dist_matrix[row_index, col_index])))
    # }
    # else
    #   #Retrieve the corresponding percent threshold from column names
    #   adjusted_rMT <- as.numeric(colnames(dist_matrix)[col_index])
    #   return (adjusted_rMT)
  }
fit_linear_model <- function(h7_tab){
    #Fit a linear model predicting H7 percent_threshold from RF stim_output
    fit <- lm(percent_threshold ~ stim_output, data = h7_tab)
    return (fit)
  }

  get_adjusted_rMT_by_fit <- function (original_rf_percent, fit_h7){
    
    pred_percent_rf <- predict(fit_h7, newdata = data.frame(stim_output = original_rf_percent))
    r_squared <- summary(fit_h7)$r.squared
    return (list(adjusted_rMT = round(pred_percent_rf), r_squared = r_squared))
  }

  get_adjusted_rMT_exact_or_fit <- function(original_rf_percent, original_rf_mso, dist_matrix, fit_h7){
    #Find the row corresponding to original_rf_percent
    row_index <- which(rownames(dist_matrix) == original_rf_percent)
    if (length(row_index) == 0){
      return (NA)
    }
    #Find the column index of the minimal distance in that row
    col_index <- which(dist_matrix[row_index, ] == 0)
    if (!is.null(col_index) && length(col_index) == 1){
      corresponding_percents <- as.numeric(colnames(dist_matrix)[col_index]) 
      adjusted_rMT <- round(mean(corresponding_percents))
      r_squared <- 1
    } else {
 
      adjusted_rMT_res <- get_adjusted_rMT_by_fit(original_rf_mso, fit_h7)
      adjusted_rMT <- adjusted_rMT_res$adjusted_rMT    
      r_squared <- adjusted_rMT_res$r_squared
    }
    return (list(adjusted_rMT = adjusted_rMT, r_squared = r_squared))
  }

  sides <- unique(p2p$side)
  subjects <- unique(p2p$patient_code)
  coils <- unique(p2p$coil)
  p2p$adjusted_MT <- p2p$percent_threshold
  p2p$adjusted_MT_min_dist <- NA
  p2p$adjusted_MT_r_squared <- 1



  for (subj in subjects){
    for (sid in sides){

      cat(sprintf("%s, %s:",subj, sid))
      #find the H7 MSO corresponding to fixed percent:
      
      rf_tab <- p2p %>%
                  filter(patient_code == subj & coil == 'rf' & side == sid) %>%
                  group_by(percent_threshold) %>%
                  summarise(stim_output = mean(machine_setting))
      
      h7_tab <- p2p %>%
                  filter(patient_code == subj & coil == 'h7' & side == sid) %>%
                  group_by(percent_threshold) %>%
                  summarise(stim_output = mean(machine_setting))
      
      fit_h7 <- fit_linear_model(h7_tab)

      dist <- compute_dist_matrix(rf_tab, h7_tab)
       
     #Compute the distance between the observed RF machine_setting and all H7 machine_settings
      subset <- which(p2p$patient_code == subj & p2p$coil == 'rf' & p2p$side == sid)
      for (i in subset){
        #adjusted_rMT_res <- get_adjusted_rMT_exact_or_fit(p2p$percent_threshold[i], p2p$machine_setting[i], dist, fit_h7)
        adjusted_rMT_res <- get_adjusted_rMT_exact_or_fit(p2p$percent_threshold[i], p2p$machine_setting[i], dist, fit_h7)
        p2p$adjusted_MT[i] <- adjusted_rMT_res$adjusted_rMT
        p2p$adjusted_MT_r_squared[i] <- adjusted_rMT_res$r_squared
        
      }
    }
  }

  p2p$coil[p2p$coil == 'rf'] <- 'RF'
  p2p$coil[p2p$coil == 'h7'] <- 'H7'
  p2p$side[p2p$side == 'lt'] <- 'Left'
  p2p$side[p2p$side == 'rt'] <- 'Right'

  return(p2p)
}

factorize <- function (df, vnames){
  for (v in vnames){
    df[[v]] <- as.factor(df[[v]])
  }
  return (df)
}

argmin_edges <- function(a, b) {
  sapply(a, function(x) {
    # candidates in b greater than x
    candidates <- which(b >= x)
    if (length(candidates) == 0) {
      return(NA_integer_)  # no edge greater than x
    } else {
      # among those, pick the one with smallest b value
      return(candidates[which.min(b[candidates])])
    }
  })
}


# Function to verify file integrity against stored hash
verify_file_hash <- function(filepath, checksum_file = "CHECKSUMS.md") {
  if (!file.exists(filepath)) {
    stop("File does not exist: ", filepath)
    return(FALSE)
  }
  
  if (!file.exists(checksum_file)) {
    warning("Checksum file not found: ", checksum_file)
    cat("Please verify the file manually or regenerate checksums.\n")
    return(TRUE)  # Allow execution but warn user
  }
  
  # Compute current hash
  current_hash <- compute_file_hash(filepath)
  
  # Read checksums file
  checksums <- readLines(checksum_file)
  
  # Extract filename from path
  filename <- basename(filepath)
  
  # Find the line with this file's hash
  hash_line <- grep(filename, checksums, value = TRUE)
  
  if (length(hash_line) == 0) {
    warning("No checksum found for: ", filename)
    cat("Please verify the file manually or regenerate checksums.\n")
    return(TRUE)  # Allow execution but warn user
  }
  
  # Extract expected hash (format: "hash  filename" or "- `hash` - filename")
  expected_hash <- sub(".*`([a-f0-9]{64})`.*", "\\1", hash_line)
  
  # If markdown format didn't match, try simple format
  if (nchar(expected_hash) != 64) {
    expected_hash <- sub("^([a-f0-9]{64}).*", "\\1", hash_line)
  }
  
  # Compare hashes
  if (current_hash == expected_hash) {
    cat("✓ File integrity verified:", filename, "\n")
    return(TRUE)
  } else {
    cat("✗ Hash mismatch for:", filename, "\n")
    cat("  Expected:", expected_hash, "\n")
    cat("  Got:     ", current_hash, "\n")
    return(FALSE)
  }
}

# ==== Modelling Functions ====
make_sampler_params_list <- function(iter_sampling = 1000,
                               iter_warmup = 1000,
                               chains = 4, 
                               parallel_chains = 4,
                               adapt_delta = 0.999, 
                               max_treedepth = 20,
                               refresh = 500){
  params <-
    list(
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      chains = chains,
      parallel_chains = parallel_chains,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = refresh
    ) 
}


make_priors_list <- list(
 
IOcurve9 = function(

  # Population means
  prior_mu_theta_mean = 50,
  prior_mu_theta_sd   = 10,

  prior_mu_slope_mean = 0.15,
  prior_mu_slope_sd   = 0.05,

  prior_mu_A_mean = 1,
  prior_mu_A_sd   = 0.2,

  # Population scale priors
  prior_sigma_theta_loc   = 0,
  prior_sigma_theta_scale = 1,

  prior_sigma_slope_loc   = 0,
  prior_sigma_slope_scale = 0.03,

  prior_sigma_A_loc   = 0,
  prior_sigma_A_scale = 0.3,

  prior_sigma_loc   = 0,
  prior_sigma_scale = 0.3,

  # Coil/side effects on θ
  prior_delta_theta_coil_loc   = 0,
  prior_delta_theta_coil_scale = 0.3,

  prior_delta_theta_side_loc   = 0,
  prior_delta_theta_side_scale = 0.3,

  # Coil/side effects on logA
  prior_delta_logA_coil_loc   = 0,
  prior_delta_logA_coil_scale = 0.1,

  prior_delta_logA_side_loc   = 0,
  prior_delta_logA_side_scale = 0.1,

  # NEW: coil/side effects on slope
  prior_delta_slope_coil_loc   = 0,
  prior_delta_slope_coil_scale = 0.3,

  prior_delta_slope_side_loc   = 0,
  prior_delta_slope_side_scale = 0.3,

  # Upper bound for all sigmas
  sigma_ubound = 2,

  # Cap for logA_prior and logA
  logA_cap = 2
)
{
  list(

    # Population means
    prior_mu_theta_mean = prior_mu_theta_mean,
    prior_mu_theta_sd   = prior_mu_theta_sd,

    prior_mu_slope_mean = prior_mu_slope_mean,
    prior_mu_slope_sd   = prior_mu_slope_sd,

    prior_mu_A_mean = prior_mu_A_mean,
    prior_mu_A_sd   = prior_mu_A_sd,

    # Population scale priors
    prior_sigma_theta_loc   = prior_sigma_theta_loc,
    prior_sigma_theta_scale = prior_sigma_theta_scale,

    prior_sigma_slope_loc   = prior_sigma_slope_loc,
    prior_sigma_slope_scale = prior_sigma_slope_scale,

    prior_sigma_A_loc   = prior_sigma_A_loc,
    prior_sigma_A_scale = prior_sigma_A_scale,

    prior_sigma_loc   = prior_sigma_loc,
    prior_sigma_scale = prior_sigma_scale,

    # Coil/side effects on θ
    prior_delta_theta_coil_loc   = prior_delta_theta_coil_loc,
    prior_delta_theta_coil_scale = prior_delta_theta_coil_scale,

    prior_delta_theta_side_loc   = prior_delta_theta_side_loc,
    prior_delta_theta_side_scale = prior_delta_theta_side_scale,

    # Coil/side effects on logA
    prior_delta_logA_coil_loc   = prior_delta_logA_coil_loc,
    prior_delta_logA_coil_scale = prior_delta_logA_coil_scale,

    prior_delta_logA_side_loc   = prior_delta_logA_side_loc,
    prior_delta_logA_side_scale = prior_delta_logA_side_scale,

    # NEW: coil/side effects on slope
    prior_delta_slope_coil_loc   = prior_delta_slope_coil_loc,
    prior_delta_slope_coil_scale = prior_delta_slope_coil_scale,

    prior_delta_slope_side_loc   = prior_delta_slope_side_loc,
    prior_delta_slope_side_scale = prior_delta_slope_side_scale,

    # Upper bound for all sigmas
    sigma_ubound = sigma_ubound,

    # Cap for logA_prior and logA
    logA_cap = logA_cap
  )
}

)


dims_from_index_strings <- function(names) {
  if (length(names) == 0L) return(integer(0))

  
}


strip_brackets <- function(x){
  y <- sub(pattern = "\\[(.*)\\]", replacement = "", x = x)
  return(y)
}

get_parameters_names <- function(draws, clean = TRUE){
  pnames <- unique(sapply(colnames(draws), strip_brackets))
  pnames_clean <- pnames[!pnames %in% c(".chain", ".iteration", ".draw", "lp__")]
  if (clean){
    return(pnames_clean)
  } else {
    return(pnames)
  }
}

extract_bracket_matrix <- function(x) {
    
    #browser()
    # Find bracketed parts
    m <- gregexpr("\\[(.*?)\\]", x)
    inside <- regmatches(x, m)
    
    # Flatten list and remove empty entries
    inside <- inside[sapply(inside, length) > 0]
    
    # If no bracketed numbers found → return 1×1 matrix with 1
    if (length(inside) == 0) {
      return(matrix(1, nrow = 1, ncol = 1))
    }
    
    # Remove brackets and split into numeric vectors
    numeric_rows <- lapply(inside, function(s) {
      s <- gsub("^\\[|\\]$", "", s)   # remove [ ]
      as.numeric(strsplit(s, ",")[[1]])
    })
    
    # All rows should have equal length
    max_len <- max(lengths(numeric_rows))
    
    # Build output matrix
    out <- matrix(NA_real_, nrow = length(numeric_rows), ncol = max_len)
    for (i in seq_along(numeric_rows)) {
      out[i, seq_along(numeric_rows[[i]])] <- numeric_rows[[i]]
    }
  
    return(out)
}

draws_to_array <- function(draws, parameters = NULL, exact = FALSE){
  n_iter <- length(draws$.iteration)
  out <- list()
  
  if (is.null(parameters)){
    parameters <- get_parameters_names(draws)
  }

  

  for (parameter in parameters){
    cat(sprintf('Converting %s...\n',parameter))
    
    param_colnames1 <- grep(paste0("^",parameter,"$"),colnames(draws), value = TRUE)
    param_colnames2 <- grep(paste0("^",parameter,"\\["),colnames(draws), value = TRUE)
    param_colnames <- c(param_colnames1, param_colnames2)
    # check if the parameter name contains brackets:
    
    index_matrix <- extract_bracket_matrix(param_colnames)
    
    #browser()
    if (is.null(index_matrix) || all(is.na(index_matrix))){
      index_matrix <- matrix(1, nrow = 1, ncol = 1)
      n_dims <- 0
      arr <- array(NA, dim = c(n_iter,1))
    }
    else{
      n_dims <- dim(index_matrix)[2]
      arr <- array(NA, dim = c(n_iter, apply(index_matrix, 2, max)))
    }

    for (i in seq_along(param_colnames)){
       indices <- index_matrix[i,]
      
       if (length(indices) == 1){
         # for 1D parameters, this is much faster:
        arr[,indices] <- draws[[param_colnames[[i]]]]
       } else {
        # for multi-dimensional parameters, need to use do.call with `[<-`:
        idx <- as.list(index_matrix[i, ])
        # Build full index: first dimension is always ':', then indices
        full_idx <- c(list(seq_len(n_iter)), idx)
        # Assign using do.call
        arr <- do.call(`[<-`, c(list(arr), full_idx, list(value = draws[[param_colnames[[i]]]])))
       }
    }
    out[[parameter]] <- arr
  }
    
  return(out)
  
}

array_to_df <- function(arr, levels_list) {
  
  # Number of posterior draws (first dimension)
  n_draws <- dim(arr)[1]
  
  # Factor names
  factor_names <- names(levels_list)
  
  # Number of factors
  n_factors <- length(levels_list)
  
  # Check dimensions
  if (length(dim(arr)) != n_factors + 1) {
    stop("dim(arr) must be length(levels_list) + 1")
  }

  # Build all combinations of factor levels
  factor_grid <- expand.grid(levels_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  
  # Total number of combinations
  n_combinations <- nrow(factor_grid)
  
  # Repeat factor grid for each draw
  df <- factor_grid[rep(seq_len(n_combinations), each = n_draws), ]
  
  if (is.null(dim(df))) # Handle case with single factor
    {
      df <- data.frame(df)
      names(df) <- factor_names
  }

  # Add Value column
  df$Value <- as.vector(arr)
  
  # Convert factor columns to factors
  df[factor_names] <- lapply(df[factor_names], factor)
  
  # Reorder columns: Value first
  df <- df[, c("Value", factor_names)]
  
  return(tibble(df))
  }

compute_contrasts <- function(arr, contrasts_list){
  out <- matrix(NA, nrow = dim(arr)[1], ncol = length(contrasts_list))
  for (cont_index in seq_along(contrasts_list)){
    contrast <- contrasts_list[[cont_index]]
    if (dims(contrast) == 1){
      A <- contrasts_list[[cont_index]][1]
      B <- contrasts_list[[cont_index]][2]
      out[,cont_index] <- arr[,A] - arr[,B]
    } else if (dims(contrast) == 2){
      A1 <- contrasts_list[[cont_index]][1,1]
      A2 <- contrasts_list[[cont_index]][2,1]
      B1 <- contrasts_list[[cont_index]][1,2]
      B2 <- contrasts_list[[cont_index]][2,2]

      out[, cont_index] <- arr[, A1, A2] - arr[, B1, B2]
    }
    
  }
  return(list(arr = out, contrast_names = names(contrasts_list)))
}

plot_contrasts <- function(contrasts_df, contrast_name, title, ci = 0.89){
  
  g <- ggplot(contrasts_df, aes(x = Value, fill = .data[[contrast_name]])) +
    geom_density(alpha = 0.5) +
    labs(
      title = title,
      x = "Difference",
      y = "Density"
    ) +
    theme_minimal(base_size = 14)

  print(g)
  
  
  info <- contrasts_df %>%
    group_by(.data[[contrast_name]]) %>%
    summarise(
      median = median(Value),
      lower = HPDI(Value, ci)[1],
      upper = HPDI(Value, ci)[2],
      area_above_zero = mean(Value > 0)
    )
  
  cat(sprintf("\n%s:\n", title))
  print(info)

  invisible(info)

}

prior_predictive <- function(mod,
                             data_list_template,
                             priors_list,
                             sampler_params = make_sampler_params_list(),
                             ci = 0.99) {
  # data_list_template: same structure as for posterior fit, but N can be 0
  # prior_list: named list of prior hyperparameters, see make_prior_list below
  


  
  data_list <- c(
    data_list_template,
    priors_list,
    list(use_likelihood = 0L)  # ensure we skip likelihood
  )
  

  fit_prior <- mod$sample(
    data = data_list,
    iter_warmup = sampler_params$iter_warmup,
    iter_sampling = sampler_params$iter_sampling,
    chains = sampler_params$chains,
    parallel_chains = sampler_params$parallel_chains,
    max_treedepth = sampler_params$max_treedepth,
    adapt_delta = sampler_params$adapt_delta,
    refresh = sampler_params$refresh,
    #init = 1,
    fixed_param = FALSE
  )
    
  draws <- fit_prior$draws(format = "df")
  # Extract the prior:

  Y_prior <- draws_to_array(draws, parameters = c("Y_prior"))$Y_prior


  df_hpdi <- NULL
  # Compute the HPDI of with credible interval of alpha for Y_prior:
  Y_prior_hpdi <- t(apply(Y_prior, 2, HPDI, prob = ci))

  attr(data_list_template,"coils")[3] <- 'T360° Coil Array'

  df_hpdi <- data.frame(
  Intensity = data_list_template$Intensity,
  coil = factor(data_list_template$coil, levels = 1:data_list_template$J, labels = attr(data_list_template, "coils")),
  side = factor(data_list_template$side, levels = 1:data_list_template$K, labels = attr(data_list_template, "sides")),
  Y = data_list_template$Y,
  lower = Y_prior_hpdi[, 1],
  upper = Y_prior_hpdi[, 2]
  )

  # Plot with shaded HPDI
  g <- ggplot(df_hpdi, aes(x = Intensity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = coil), 
              alpha = 0.3) +
  geom_point(aes(x = Intensity, y= Y), color = 'black') +
  facet_grid(side ~ coil) +
  labs(
    title = sprintf("Prior Predictive Distribution with %d%% HPDI",ci*100),
    x = "Stimulator Intensity (% of MSO)",
    y = "MEP Amplitude (peak-to-peak)"
  ) +
  theme_minimal()

  print(g)

  #Plot general distributions of Y and Yprior
  df_long <- data.frame(
    Y = c(as.vector(Y_prior), data_list_template$Y),
    type = as.factor(
              c(
                rep("Y_prior", length(as.vector(Y_prior))),
                rep("Y_observed", length(data_list_template$Y))
              )
            )
  )

  max_data <- max(data_list_template$Y)

  g2 <- ggplot(df_long, aes(x = Y, fill = type)) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Prior Predictive Distribution vs Observed Data",
      x = "MEP Amplitude (peak-to-peak)",
      y = "Density"
    ) + xlim(0, max_data) 
  
  print(g2)

    # You can inspect distribution of alpha_cs, mu_100, mu_110, mu_120, etc.
  return(list(fit = fit_prior,
               draws = draws,
               Y_prior = Y_prior,
               df_hpdi = df_hpdi))
}


# Posterior predictive check with numeric summaries
posterior_predictive_check_numeric <- function(fit, data, n_samples = 1000) {
  library(dplyr)
  
  # Extract observed data
  Y_obs <- data$Y
  
  # Extract posterior predictive samples
  Y_rep <- fit$draws(variables = "Y_rep", format = "matrix")
  
  # Sample if we have more draws than needed
  if (nrow(Y_rep) > n_samples) {
    Y_rep <- Y_rep[sample(nrow(Y_rep), n_samples), ]
  }
  
  # 1. Bayesian p-values for summary statistics
  cat("=== Bayesian p-values (posterior predictive p-values) ===\n")
  
  # Mean
  T_obs_mean <- mean(Y_obs)
  T_rep_mean <- apply(Y_rep, 1, mean)
  ppp_mean <- mean(T_rep_mean >= T_obs_mean)
  cat(sprintf("Mean: p = %.3f (observed = %.3f)\n", ppp_mean, T_obs_mean))
  
  # Median
  T_obs_median <- median(Y_obs)
  T_rep_median <- apply(Y_rep, 1, median)
  ppp_median <- mean(T_rep_median >= T_obs_median)
  cat(sprintf("Median: p = %.3f (observed = %.3f)\n", ppp_median, T_obs_median))
  
  # SD
  T_obs_sd <- sd(Y_obs)
  T_rep_sd <- apply(Y_rep, 1, sd)
  ppp_sd <- mean(T_rep_sd >= T_obs_sd)
  cat(sprintf("SD: p = %.3f (observed = %.3f)\n", ppp_sd, T_obs_sd))
  
  # IQR
  T_obs_iqr <- IQR(Y_obs)
  T_rep_iqr <- apply(Y_rep, 1, IQR)
  ppp_iqr <- mean(T_rep_iqr >= T_obs_iqr)
  cat(sprintf("IQR: p = %.3f (observed = %.3f)\n", ppp_iqr, T_obs_iqr))
  
  # Min/Max
  T_obs_min <- min(Y_obs)
  T_rep_min <- apply(Y_rep, 1, min)
  ppp_min <- mean(T_rep_min <= T_obs_min)
  cat(sprintf("Min: p = %.3f (observed = %.3f)\n", ppp_min, T_obs_min))
  
  T_obs_max <- max(Y_obs)
  T_rep_max <- apply(Y_rep, 1, max)
  ppp_max <- mean(T_rep_max >= T_obs_max)
  cat(sprintf("Max: p = %.3f (observed = %.3f)\n\n", ppp_max, T_obs_max))
  
  # 2. Percentage of observations within posterior predictive intervals
  cat("=== Coverage: Observed data within posterior predictive intervals ===\n")
  
  # For each observation, compute 50%, 89%, 95% intervals from Y_rep
  pct_in_50 <- mean(sapply(1:ncol(Y_rep), function(i) {
    interval <- quantile(Y_rep[, i], probs = c(0.25, 0.75))
    Y_obs[i] >= interval[1] & Y_obs[i] <= interval[2]
  }))
  
  pct_in_89 <- mean(sapply(1:ncol(Y_rep), function(i) {
    interval <- quantile(Y_rep[, i], probs = c(0.055, 0.945))
    Y_obs[i] >= interval[1] & Y_obs[i] <= interval[2]
  }))
  
  pct_in_95 <- mean(sapply(1:ncol(Y_rep), function(i) {
    interval <- quantile(Y_rep[, i], probs = c(0.025, 0.975))
    Y_obs[i] >= interval[1] & Y_obs[i] <= interval[2]
  }))
  
  cat(sprintf("50%% intervals: %.1f%% of observations (expected: 50%%)\n", pct_in_50 * 100))
  cat(sprintf("89%% intervals: %.1f%% of observations (expected: 89%%)\n", pct_in_89 * 100))
  cat(sprintf("95%% intervals: %.1f%% of observations (expected: 95%%)\n\n", pct_in_95 * 100))
  
  # 3. Quantile comparison
  # cat("=== Quantile Comparison ===\n")
  # quantiles <- c(0.05, 0.25, 0.50, 0.75, 0.95)
  # obs_quantiles <- quantile(Y_obs, quantiles)
  
  # # Compute quantiles for each posterior sample
  # rep_quantiles_mean <- apply(Y_rep, 1, function(x) quantile(x, quantiles))
  
  # # Summarize across posterior samples (using numeric indexing)
  # rep_quantiles_summary <- apply(rep_quantiles_mean, 1, function(x) {
  #   c(mean = mean(x), sd = sd(x), 
  #     lower = quantile(x, 0.025), upper = quantile(x, 0.975))
  # })
  
  # quantile_df <- data.frame(
  #   Quantile = paste0(quantiles * 100, "%"),
  #   Observed = sprintf("%.3f", obs_quantiles),
  #   Predicted_Mean = sprintf("%.3f", rep_quantiles_summary["mean", ]),
  #   Predicted_95CI = sprintf("[%.3f, %.3f]", 
  #                            rep_quantiles_summary["lower", ],
  #                            rep_quantiles_summary["upper", ])
  # )
  # print(quantile_df, row.names = FALSE)
  
  # 4. Return results as list
  invisible(list(
    ppp_values = data.frame(
      Statistic = c("Mean", "Median", "SD", "IQR", "Min", "Max"),
      Observed = c(T_obs_mean, T_obs_median, T_obs_sd, T_obs_iqr, T_obs_min, T_obs_max),
      p_value = c(ppp_mean, ppp_median, ppp_sd, ppp_iqr, ppp_min, ppp_max)
    ),
    coverage = data.frame(
      Interval = c("50%", "89%", "95%"),
      Observed_Coverage = c(pct_in_50, pct_in_89, pct_in_95) * 100,
      Expected_Coverage = c(50, 89, 95)
    ),
    #quantiles = quantile_df,
    Y_rep = Y_rep
  ))
}

plot_IOcurve <- function(fit, data, n_points = 200, ci = 0.89) {

  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)

  # Extract posterior draws
  draws <- as_draws_df(fit)

  # Dimensions
  I <- max(data$subj)
  J <- max(data$coil)
  K <- max(data$side)

  coils <- attr(data, "coils")
  coils[1] <- 'figure-of-8'
  coils[3] <- 'T360° coil-array'
  sides <- attr(data, "sides")

  # Intensity grid
  Intensity_grid <- seq(min(data$Intensity),
                        max(data$Intensity),
                        length.out = n_points)

  extract_array <- function(draws, name, dim) {
  # Stan outputs names like: name.1.1.1
  idx <- grep(paste0("^", name, "\\."), names(draws))

  if (length(idx) == 0) {
    stop(paste("No variables found matching", name))
  }

  arr <- as.matrix(draws[, idx])

  # reshape into draws × dim1 × dim2 × dim3
  array(arr, dim = c(nrow(draws), dim))
}
  
  # Extract arrays
  p_arrs <- draws_to_array(draws, c("theta_raw", "slope_raw", "A_raw"))

  theta_arr <- p_arrs$theta_raw
  slope_arr <- p_arrs$slope_raw
  A_arr     <- p_arrs$A_raw

  

  # Average over subjects → population-level parameters
  theta_mean <- apply(theta_arr, c(1,3,4), mean)   # draws × J × K
  slope_mean <- apply(slope_arr, c(1,3,4), mean)
  A_mean     <- apply(A_arr,     c(1,3,4), mean)

  

  # Logistic function
  logistic_fun <- function(I, A, slope, theta) {
    A / (1 + exp(-slope * (I - theta)))
  }


  # Compute posterior curves
  curve_df <- map_dfr(1:J, function(j) {
    map_dfr(1:K, function(k) {
      map_dfr(1:length(Intensity_grid), function(ii) {
        Ival <- Intensity_grid[ii]

        mu_draws <- logistic_fun(
          I = Ival,
          A = A_mean[, j, k],
          slope = slope_mean[, j, k],
          theta = theta_mean[, j, k]
        )

        tibble(
          coil = coils[j],
          side = sides[k],
          Intensity = Ival,
          mu = mu_draws
        )
      })
    })
  })

  # Summarize posterior
  # alpha <- (1 - ci) / 2

  curve_summary <- curve_df %>%
    group_by(coil, side, Intensity) %>%
    summarise(
      mu_mean = mean(mu),
      mu_low  = HPDI(mu, ci)[1],
      mu_high = HPDI(mu, ci)[2],
      .groups = "drop"
    )
  
  
  font_name = "Times New Roman"
  font_size = 11
  # Plot
  
  p <- ggplot(curve_summary,
              aes(x = Intensity, y = mu_mean, color = factor(coil))) +
    geom_ribbon(aes(ymin = mu_low, ymax = mu_high,
                    fill = factor(coil), color = factor(coil)),
                alpha = 0.25 ) +
    geom_line(linewidth = 1.3) +
    facet_wrap(~side, nrow = 1, labeller = labeller(side = c("Left" = "Left Hemisphere", "Right" = "Right Hemisphere"))) +
    scale_color_manual(values = c("#3399FF", "#66CC66", "#FF3399")) +
    scale_fill_manual(values = c("#3399FF", "#66CC66", "#FF3399")) +
   labs(
    x = "Stimulation Intensity (% MSO)",
    y = "Predicted MEP Amplitude (mV)",
    title = sprintf("Population IO Curves (%d%% HPDI)", round(ci * 100)),
    color = "Coil Type",
    fill = "Coil Type"
  ) +
    theme_classic(base_size = 10, base_family = "Times New Roman") +
theme(
  strip.background = element_blank(),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black", family = "Times New Roman"),
  axis.title = element_text(size = 12, face = "bold", color = "black", family = "Times New Roman"),
  axis.text = element_text(size = 12, color = "black", family = "Times New Roman"),
  strip.text = element_text(size = 12, face = "bold", color = "black", family = "Times New Roman"),
  legend.title = element_text(size = 12, face = "bold", color = "black", family = "Times New Roman"),
  legend.text = element_text(size = 12, color = "black", family = "Times New Roman")
  # plot.margin = margin(t = 10, r = 10, b = 10, l = 15, unit = "pt")
)




  print(p)
  invisible(p)
}

simulate_coil_side_effects <- function(
  fit,
  data,
  sigma_sim = NULL,
  refit_fun = NULL,   # function that refits the model: function(data_list) -> fit
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  library(posterior)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)

  # --- Basic structure from data list ---

  N <- data$N
  I <- data$I
  J <- data$J
  K <- data$K

  subj      <- data$subj
  coil      <- data$coil
  side      <- data$side
  Intensity <- data$Intensity

  # --- Extract posterior draws ---

  draws <- as_draws_df(fit)

  get_scalar_median <- function(name) {
    if (!name %in% names(draws))
      stop(paste("Parameter", name, "not found in fit"))
    median(draws[[name]])
  }

  get_indexed_median <- function(name, dim_len) {
    idx <- grep(paste0("^", name, "\\["), names(draws))
    if (length(idx) != dim_len)
      stop(paste("Expected", dim_len, "columns for", name, "but found", length(idx)))
    mat <- as.matrix(draws[, idx])
    apply(mat, 2, median)
  }

  # Population means
  mu_theta <- get_scalar_median("mu_theta")
  mu_slope <- get_scalar_median("mu_slope")
  mu_logA  <- get_scalar_median("mu_logA")

  # Scales
  sigma_theta <- get_scalar_median("sigma_theta")
  sigma_slope <- get_scalar_median("sigma_slope")
  sigma_logA  <- get_scalar_median("sigma_logA")

  # Observation noise
  sigma_post <- get_scalar_median("sigma")
  if (is.null(sigma_sim)) sigma_sim <- sigma_post

  # Coil/side deltas
  delta_theta_coil <- get_indexed_median("delta_theta_coil", J)
  delta_theta_side <- get_indexed_median("delta_theta_side", K)

  delta_slope_coil <- get_indexed_median("delta_slope_coil", J)
  delta_slope_side <- get_indexed_median("delta_slope_side", K)

  delta_logA_coil <- get_indexed_median("delta_logA_coil", J)
  delta_logA_side <- get_indexed_median("delta_logA_side", K)

  # --- Construct "true" coil/side parameters for Scenario 1 (effects as estimated) ---

  # For now we ignore subject-level variation in the simulation (can be added later)
  theta_true_1 <- matrix(NA_real_, nrow = J, ncol = K)
  slope_true_1 <- matrix(NA_real_, nrow = J, ncol = K)
  logA_true_1  <- matrix(NA_real_, nrow = J, ncol = K)

  for (j in 1:J) {
    for (k in 1:K) {
      theta_true_1[j, k] <-
        mu_theta +
        delta_theta_coil[j] +
        delta_theta_side[k]

      slope_true_1[j, k] <-
        mu_slope +
        delta_slope_coil[j] +
        delta_slope_side[k]

      logA_true_1[j, k] <-
        mu_logA +
        delta_logA_coil[j] +
        delta_logA_side[k]
    }
  }

  # --- Construct "true" parameters for Scenario 0 (no coil/side effects) ---

  theta_true_0 <- matrix(mu_theta, nrow = J, ncol = K)
  slope_true_0 <- matrix(mu_slope, nrow = J, ncol = K)
  logA_true_0  <- matrix(mu_logA,  nrow = J, ncol = K)

  # --- Logistic IO function ---

  logistic_fun <- function(I, A, slope, theta) {
    A / (1 + exp(-slope * (I - theta)))
  }

  # --- Generate Y_sim under Scenario 1 (effects as estimated) ---

  Y_sim1 <- numeric(N)
  for (n in 1:N) {
    j <- coil[n]
    k <- side[n]

    A_jk    <- exp(logA_true_1[j, k])
    slope_jk <- slope_true_1[j, k]
    theta_jk <- theta_true_1[j, k]

    mu_curve <- logistic_fun(Intensity[n], A_jk, slope_jk, theta_jk)
    mu_log   <- log(mu_curve + 1e-9)

    Y_sim1[n] <- rlnorm(1, meanlog = mu_log, sdlog = sigma_sim)
  }

  # --- Generate Y_sim under Scenario 0 (no coil/side effects) ---

  Y_sim0 <- numeric(N)
  for (n in 1:N) {
    j <- coil[n]
    k <- side[n]

    A_jk    <- exp(logA_true_0[j, k])
    slope_jk <- slope_true_0[j, k]
    theta_jk <- theta_true_0[j, k]

    mu_curve <- logistic_fun(Intensity[n], A_jk, slope_jk, theta_jk)
    mu_log   <- log(mu_curve + 1e-9)

    Y_sim0[n] <- rlnorm(1, meanlog = mu_log, sdlog = sigma_sim)
  }

  # --- Build new data lists ---

  data_sim1 <- data
  data_sim1$Y <- Y_sim1

  data_sim0 <- data
  data_sim0$Y <- Y_sim0

  # --- Optionally refit the model to simulated data ---

  fit_sim1 <- fit_sim0 <- NULL

  if (!is.null(refit_fun)) {
    fit_sim1 <- refit_fun(data_sim1)
    fit_sim0 <- refit_fun(data_sim0)
  }

  list(
    data_sim1 = data_sim1,
    data_sim0 = data_sim0,
    fit_sim1  = fit_sim1,
    fit_sim0  = fit_sim0,
    true_params = list(
      scenario1 = list(
        theta = theta_true_1,
        slope = slope_true_1,
        logA  = logA_true_1
      ),
      scenario0 = list(
        theta = theta_true_0,
        slope = slope_true_0,
        logA  = logA_true_0
      ),
      sigma_sim = sigma_sim
    )
  )
}


plot_simulated_logistic <- function(sim = NULL, df = NULL, true_params = NULL, u_alpha_subj = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)

  if (!is.null(sim)) {
    df <- sim$data
    true_params <- sim$true_params
    u_alpha_subj <- sim$u_alpha_subj
  }

  # Extract factor levels
  coils <- sort(unique(df$coil))
  sides <- sort(unique(df$side))
  subjs <- sort(unique(df$subj))


  # ------------------------------------------------------------
  # (1) Compute panel-level α, β, γ from the data frame
  #     (they are constant within each coil × side)
  # ------------------------------------------------------------
  param_grid <- df %>%
    group_by(coil, side) %>%
    summarise(
      alpha = unique(alpha_pop),
      beta  = unique(beta),
      gamma = unique(gamma),
      subtitle = sprintf("α = %.2f, β = %.2f, γ = %.2f",
                         unique(alpha), unique(beta), unique(gamma)),
      .groups = "drop"
    ) %>% select(-c(alpha, beta, gamma))

  # ------------------------------------------------------------
  # (2) Build subject-specific curves (grey)
  # ------------------------------------------------------------
  Irel_grid <- seq(min(df$Irel), max(df$Irel), length.out = 200)
  
  compute_mu <- function (param = 'mu', Irel, coil, side, subj, true_params, u_alpha_subj) {
    res <- compute_mu_and_params(Irel = Irel,
                                      coil = coil,
                                      side = side,
                                      subj = subj,
                                      true_params = true_params,
                                      u_alpha_subj = u_alpha_subj)
    return(res[[param]])
  }

  subj_curves <- expand.grid(
    Irel = Irel_grid,
    coil = coils,
    side = sides,
    subj = subjs) %>%
    mutate(
    mu = mapply(compute_mu, param = 'mu', Irel = Irel,
                                      coil = coil,
                                      side = side,
                                      subj = subj,
                                      MoreArgs = list(true_params = true_params,
                                      u_alpha_subj = u_alpha_subj)))
    
                                     
  # ------------------------------------------------------------
  # (3) Population-level curves (blue)
  # ------------------------------------------------------------
  pop_curves <- expand.grid(
    Irel = Irel_grid,
    coil = coils,
    side = sides) %>%
    mutate(
    mu = mapply(compute_mu, param = 'mu_pop', Irel = Irel,
                                      coil = coil,
                                      side = side,
                                      MoreArgs = list(subj = NULL,
                                        true_params = true_params,
                                        u_alpha_subj = u_alpha_subj)))
                      

  # ------------------------------------------------------------
  # (4) Attach subtitles to the observed data
  # ------------------------------------------------------------
  df <- df %>%
    left_join(param_grid, by = c("coil", "side"))

  # ------------------------------------------------------------
  # (5) Plot
  # ------------------------------------------------------------
  ggplot() +
    # subject curves
    geom_line(
      data = subj_curves,
      aes(x = Irel, y = mu, group = subj),
      color = "grey70",
      alpha = 0.6
    ) +
    # population curve
    geom_line(
      data = pop_curves,
      aes(x = Irel, y = mu),
      color = "blue",
      linewidth = 1.2
    ) +
    # simulated data
    geom_jitter(
      data = df,
      aes(x = Irel, y = Y),
      width = 0.3,
      height = 0.0,
      color = "grey30",
      alpha = 0.1,
      size = 0.5
    ) +
    facet_grid(
      side ~ coil,
      labeller = labeller(
        coil = label_value,
        side = label_value
      )
    ) +
    # Add α, β, γ subtitles inside each panel
    geom_text(
      data = param_grid,
      aes(
        x = -Inf, y = Inf,
        label = subtitle
      ),
      hjust = -0.1, vjust = 1.2,
      size = 3.5,
      color = "black"
    ) +
    labs(
      title = "Simulated TMS Data with True Logistic IO Curves",
      x = "Relative Intensity (% rMT)",
      y = "MEP Amplitude (peak-to-peak)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(size = 12, face = "bold")
    )

}

plot_meps <- function(df, Y_prior = NULL, Y_rep = NULL, subplot = FALSE) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)

  # df must contain: subj, coil, side, Y
  stopifnot(all(c("subj", "coil", "side", "Y") %in% names(df)))

  N <- nrow(df)

  # Helper: flatten posterior predictive arrays
  flatten_ppc <- function(mat, label) {
    if (is.null(mat)) return(NULL)
    if (!is.matrix(mat)) stop("Y_prior / Y_rep must be matrices: n_iter × N")

    tibble(
      value = as.vector(mat),
      iter  = rep(seq_len(nrow(mat)), each = ncol(mat)),
      obs   = rep(seq_len(ncol(mat)), times = nrow(mat)),
      type  = label
    ) %>%
      left_join(df %>% mutate(obs = row_number()), by = "obs")
  }

  # Build long-format data
  df_obs <- df %>%
    mutate(type = "Observed") %>%
    rename(value = Y)

  df_prior <- flatten_ppc(Y_prior, "Prior predictive")
  df_rep   <- flatten_ppc(Y_rep,   "Posterior predictive")

  df_all <- bind_rows(df_obs, df_prior, df_rep)

  # If subplot = FALSE → single density overlay
  if (!subplot) {
    p <- ggplot(df_all, aes(x = value, color = type, fill = type)) +
      geom_density(alpha = 0.25) +
      theme_minimal(base_size = 14) +
      labs(
        x = "MEP amplitude",
        y = "Density",
        title = "Observed vs Prior/Posterior Predictive Distributions"
      )
    return(p)
  }

  # If subplot = TRUE → facet by side ~ coil
  p <- ggplot(df_all, aes(x = value, color = type, fill = type)) +
    geom_density(alpha = 0.25) +
    facet_grid(side ~ coil) +
    theme_minimal(base_size = 14) +
    labs(
      x = "MEP amplitude",
      y = "Density",
      title = "Observed vs Prior/Posterior Predictive by Coil × Side"
    )

  return(p)
}

plot_data_distribution <- function(df) {
  df <- df %>% filter(Irel >= 90 & Irel <= 120) 
  df <- factorize(df,c('coil','side','Irel')) %>% filter()
  ggplot(df, aes(x = Y)) +
    geom_histogram(alpha = 0.5) +
    facet_grid(coil + side ~ Irel) +
    labs(
      title = "Distribution of Simulated MEP Amplitudes",
      x = "MEP Amplitude (peak-to-peak)",
      y = "Density"
    ) +
    theme_minimal(base_size = 14)
}

fit_TMS_model <- function(mod,
                          data_list,
                          priors_list,
                          sampler_params = make_sampler_params_list(),
                          out_file = NULL) {
  stan_data <- c(
    data_list,
    priors_list,
    list(use_likelihood = 1L)
  )

  fit <- mod$sample(
    data = stan_data,
    iter_sampling = sampler_params$iter_sampling,
    iter_warmup = sampler_params$iter_warmup,
    chains = sampler_params$chains,
    parallel_chains = sampler_params$parallel_chains,
    adapt_delta = sampler_params$adapt_delta,
    max_treedepth = sampler_params$max_treedepth,
    refresh = sampler_params$refresh#,
  #   init = function() list(
  #   alpha_0 = 10,
  #   beta_0 = 1,
  #   gamma_0 = 0,
  #   sigma_y = 0.5
  # )
)

  if (!is.null(out_file)) {
    fit$save_object(file = out_file)
  }

  return(fit)
}

show_diagnostics <- function(fit){

  
  # 4.2.1. Sigmas
 vars <- c("sigma", "sigma_theta", "sigma_slope", "sigma_logA",
   "mu_theta", "mu_slope", "mu_logA")
   #"delta_slope_coil","delta_slope_side","delta_theta_coil","delta_logA_coil","delta_logA_coil","delta_logA_side")
  summ <- check_convergence(fit, params = vars, show_plots = TRUE)

  return(summ)
}

show_results <- function(fit, data, ci = 0.95){
  results <- draws_to_array(fit$draws(format = 'df'), parameters = c('delta_slope_coil','delta_slope_side','delta_theta_coil','delta_logA_coil','delta_theta_side','delta_logA_side'))

  contrasts_out <- list()

  for (p in names(results)) {
    if (grepl('coil', p)) {
      contrasts_list <- list("fig8 - H7" = c(1,2), "fig8 - T360°-Coil-Array" = c(1,3), "H7 - T360°-Coil-Array" = c(2,3))
      contrast <- 'coil_contrast'
    }
    if (grepl('side', p)) {
      contrasts_list <- list("Left-Right" = c(1,2))
      contrast <- 'side_contrast'
    }
    contrasts_arr <- compute_contrasts(results[[p]], contrasts_list = contrasts_list)
    levels_list <- list()
    levels_list[[1]] <- names(contrasts_list)
    names(levels_list)[1] <- contrast
    contrasts_df <- array_to_df(contrasts_arr$arr, levels_list = levels_list)
    contrasts_out[[contrast]] <- contrasts_df
    names(contrasts_out[[contrast]])[2] <- 'contrast'
    plot_contrasts(contrasts_df, contrast, title = p, ci = ci)
  }
  


  plot_IOcurve(fit, data_actual, ci = ci)

  return(list(raw_res = results, contrasts = contrasts_out))
}


check_convergence <- function(fit, params = NULL, show_plots = FALSE) {
  if (!inherits(fit, "CmdStanMCMC")) {
    stop("fit must be a CmdStanMCMC object")
  }
  
  cat("\n====================\n Convergence Summary\n====================\n")
  
  # Extract draws
  draws <- fit$draws(variables = params)
  
  # Posterior package summary
  summ <- posterior::summarise_draws(draws)
  print(summ[, c("variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail")])
  
  # Rhat check
  cat("\n--- R-hat diagnostics ---\n")
  bad_rhat <- summ$variable[summ$rhat > 1.01]
  if (length(bad_rhat) == 0) {
    cat("All R-hat < 1.01 (good)\n")
  } else {
    cat("Parameters with R-hat > 1.01:\n")
    print(bad_rhat)
  }
  
  # ESS check
  cat("\n--- Effective Sample Size (ESS) ---\n")
  low_ess <- summ$variable[summ$ess_bulk < 400]
  if (length(low_ess) == 0) {
    cat("All ESS_bulk > 400 (good)\n")
  } else {
    cat("Parameters with ESS_bulk < 400:\n")
    print(low_ess)
  }
  
  # Divergences
  cat("\n--- Divergences ---\n")
  div <- fit$diagnostic_summary()$num_divergent
  cat("Number of divergent transitions:", div, "\n")
  
  # Treedepth
  cat("\n--- Max treedepth ---\n")
  td <- fit$diagnostic_summary()$num_max_treedepth
  cat("Transitions hitting max treedepth:", td, "\n")
  
  # Energy Bayesian fraction of missing information (E-BFMI)
  cat("\n--- Energy diagnostics (E-BFMI) ---\n")
  ebfmi <- fit$diagnostic_summary()$ebfmi
  print(ebfmi)
  if (any(ebfmi < 0.3)) {
    cat("Warning: E-BFMI < 0.3 indicates poor exploration\n")
  }
  
  # Optional traceplots
  if (show_plots) {
    cat("\n--- Traceplots ---\n")
    g <- bayesplot::mcmc_trace(draws, pars = params)
    print(g)
  }
  
  invisible(summ)
}


diagnose_fit <- function(params = params, fit = NULL, fit_file = NULL) {
  require(posterior)
  require(bayesplot)
  if (!is.null(fit_file)) {
    fit <- readRDS(fit_file)
  }
  if (!is.null(fit)){
    check_convergence(fit = fit, show_plots = TRUE, params = params)
  }
  else{
    warning("No valid fit object!")
  }
}



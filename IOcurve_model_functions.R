require(ggplot2)
require(dplyr)
require(tidyr)


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

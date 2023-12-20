#
# helper functions for all of the ESM calculation/simulation functions.
#

library(dplyr)
library(tidyr)
library(data.table)



#'Move to analysis functions?
make_err_conv_func_log <- function(.coefs, study, use_SOC){
  study_coef <- .coefs[sapply(names(.coefs),function(x) grepl(study, x))]
  if ('(Intercept)' %in% names(.coefs)){
    intercept <- .coefs['(Intercept)']}
  if ('log(pred_err)' %in% names(.coefs)){
    intercept <- study_coef
    .coef <- .coefs['log(pred_err)']
  }else{
    .coef <- study_coef
  } 
  
  primary_func <-  function(soil_sample, x) {
    exp(log((1/sqrt(sum_weights_mass(soil_sample, x))))*.coef + intercept) }
  
  
  if (use_SOC) {
    return(function(soil_sample, x){
      primary_func(soil_sample, x) * calc.linear_SOC_pct(soil_sample, x)
      
      
    })
  }else{
    return (primary_func)}
}




calc.linear_SOC_pct<- function (soil_sample, mass_to_quantify){
  soil_sample%>% 
    filter( lag(Cum_Min_Soil_g_cm2, default = 0) <= mass_to_quantify) %>%
    slice_max(Cum_Min_Soil_g_cm2) %>%
    pull(SOC_pct) %>% 
    max()
}





#'Calculate the spline precision statistic for the soil sample, based on the 
#'interpolation mass
sum_weights_mass <- function(soil_sample, mass_to_quantify){
  soil_sample <- ungroup(soil_sample)
  
  below_weights <- filter(soil_sample, Cum_Min_Soil_g_cm2 <= mass_to_quantify) %>% 
    summarize(x = sum(1/ ((Cum_Min_Soil_g_cm2 - mass_to_quantify)**2))) %>%
    pull(x)
  
  above_weights <- filter(soil_sample, Cum_Min_Soil_g_cm2 >= mass_to_quantify) %>% 
    summarize(x = sum(1/ ((Cum_Min_Soil_g_cm2 - mass_to_quantify)**2))) %>% 
    pull(x)
  
  if (length(above_weights) == 0){
    above_weights <- 0
  }
  if (length(below_weights) == 0){
    below_weights <- 0
  }
  
  return( (sqrt(below_weights)+ sqrt(above_weights))**2 )}





load_study <- function(name){
  path <- file.path('Source_Data', paste(name, 'csv', sep = '.'))
  
  data <-read.csv(path) %>% mutate(Rep = 1, study = name, ID = as.character(ID)) 
  
  if ('site' %in% colnames(data)){
    data <- mutate(data, site = as.character(site))
    data <- mutate(data, ID = paste(site, ID, sep = '_'))
  }
  if ('Ref_ID' %in% colnames(data)){
    data <- mutate(data, Ref_ID = as.character(Ref_ID))
  }
  if ('prof' %in% colnames(data)){
    data <- select(data, -c(prof))
  }
  return(data %>% calc.cumulative_masses() )
}


make_zerodepth_row <- function(soil_sample){
  return(data.frame(ID = first(soil_sample$ID),Lower_cm = 0, Upper_cm = 0,
                    Cum_Min_Soil_g_cm2 = 0, Cum_Soil_g_cm2 = 0,
                    Cum_SOC_g_cm2 = 0, Cum_SOM_g_cm2 = 0,
                    SOC_g_cm2 = 0, SOM_g_cm2 = 0, Soil_g_cm2 = 0
  ))
  
}



#'Calculate layer and cumulative masses for soil data.
#'
calc.cumulative_masses <- function(soil_data, override = F){
  # Calculate soil mass in each interval
  if ('Cum_Min_Soil_g_cm2' %in% colnames(soil_data) & !override){
    return(soil_data)
  }else{ 
    return (
      soil_data %>% mutate(Soil_g_cm2 = (Lower_cm - Upper_cm) * BD_g_cm3) %>%
        #SOC and SOM mass of intervals
        mutate(SOC_g_cm2 = (SOC_pct/100)*Soil_g_cm2,
               SOM_g_cm2 = (SOM_pct/100)*Soil_g_cm2) %>%
        mutate(Min_Soil_g_cm2 = Soil_g_cm2 - SOM_g_cm2) %>%  #mineral mass of intervals
        group_by(ID, Rep) %>%
        mutate(Cum_Soil_g_cm2 = cumsum(Soil_g_cm2),
               Cum_SOC_g_cm2 = cumsum(SOC_g_cm2),
               Cum_SOM_g_cm2 = cumsum(SOM_g_cm2),
               Cum_Min_Soil_g_cm2 = cumsum(Min_Soil_g_cm2))
    )}
}




make_zerodepth_row <- function(soil_sample){
  return(data.frame(ID = first(soil_sample$ID),Lower_cm = 0, Upper_cm = 0,
                    Cum_Min_Soil_g_cm2 = 0, Cum_Soil_g_cm2 = 0,
                    Cum_SOC_g_cm2 = 0, Cum_SOM_g_cm2 = 0,
                    SOC_g_cm2 = 0, SOM_g_cm2 = 0, Soil_g_cm2 = 0
  ))
  
}

#' aggregate a dataframe of soil_data to a given depth.
agg_to_depth <- function(soil_data, bottom_depth, top_depth){
  
  to_depth <- soil_data %>% 
    filter(Lower_cm <= bottom_depth & Upper_cm >= top_depth )
  
  depths <- to_depth$Lower_cm -to_depth$Upper_cm
  soc_weights <- depths * to_depth$BD_g_cm3
  
  return(
    to_depth %>% 
      summarize(
        SOC_pct = weighted.mean(SOC_pct, soc_weights), 
        SOM_pct = weighted.mean(SOM_pct, soc_weights), 
        BD_g_cm3 = weighted.mean(BD_g_cm3, depths),
        Upper_cm = min(Upper_cm),
        Lower_cm = max(Lower_cm),
        #Ref_ID = first(Ref_ID), 
        #ID = first(ID)
      )
    
  )
  
}

#'aggregate the soil data to several depths, to simulate soil-sampling
#'for instance: if new_depths = c(30, 50) simulates soil sampling to the depth of 
#'30cm and 50cm 
agg_to_depths <- function(soil_data, new_depths){
  new_depths <- new_depths[order(new_depths)]
  if (!all(new_depths %in% soil_data$Lower_cm)){
    stop(gettextf("Error, can only aggregate to soil depths found in existing sample"))
  }
  
  new_tops <- dplyr::lag(new_depths) %>% replace_na(0)
  
  new_soil_data <-data.frame()
  for (i in seq(1:length(new_depths))){
    
    new_soil_data <- rbind(
      new_soil_data, 
      agg_to_depth(soil_data, new_depths[i], new_tops[i])) %>% 
      arrange(Lower_cm)
  }
  
  return (new_soil_data)
}

#' Turn a dataframe of soil cores split into depths into a dataframe where these depths have been 
#' aggregated to give the results for a hypothetical sampling event with fewer soil depths
df_agg_to_depths <- function(soil_data, new_depths){
  
  return (soil_data %>% group_by(ID, Rep, Ref_ID) %>% 
            group_modify(~ agg_to_depths(.,  
                                         new_depths = new_depths
            ), .keep = T
            ))
}


dif_matrix <- function(v){
  m <- matrix(outer(v, v, "-"), nrow = length(v))
  diag(m) <- NA
  return(m)}

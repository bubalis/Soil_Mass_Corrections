#
# Linear interpolation methods for ESM corrections

#source('util.R')






#' Run a linear mass correction as described by Fowler et al.
#' This version has a few useful modifications:
#' 1: an adjustment factor has been added. The adjustment factor is an estimate of the ratio
#' of the SOC content of the lowest layer before or including the depth of quantification
#' and the SOC content of the soil mass in the depth-adjustment layer.
#'  Running with adjust_factor = 1 gives the
#' same results as the Fowler et. al method.
#' Adjust_factor == 'standard_correction' gives a version that weights the fixed depth and
#' Mass-corrected approaches differently, giving more weight to the mass-corrected approach
#' when the final soil mass measured starts closer to the depth of quantification.
#'
#' 2: This version allows including soil samples that are all or partially below the depth of quantification.
#' This allows for simulating more mass-correction sampling/testing approaches.
#'
#' @param intial_soil : a dataframe describing the soil at t = 0
#' @param sample: a dataframe describing the soil at t =1
#' both dataframes should have the following columns:
#' SOC_pct : SOC percent content of the layer
#' SOM_pct: SOM percent of they layer
#' BD_g_cm3: Bulk density of the layer
#' Upper_cm : depth of the start of the layer
#' Lower_cm : depth of the bottom of the layer
#'
#' @param quantification_depth : depth, in cm that we are trying to get a mass-correct SOC estimate for
#' @param adjust_factor : adjustment factor. See above.
#'
linear_mass_correction <- function(intial_soil, sample, quantification_depth =30,
                                   adjust_factor = 1
){


  if (is.null(adjust_factor) | adjust_factor == 'standard_correction' ){
    adjust_factor <- .5
  }

  sample <- calc.cumulative_masses(sample)
  intial_soil <- calc.cumulative_masses(intial_soil)


  mass1 <- intial_soil %>%
    filter(Lower_cm == quantification_depth) %>%
    pull(Cum_Min_Soil_g_cm2)

  mass2 <- sample %>%
    filter(Lower_cm == quantification_depth) %>%
    pull(Cum_Min_Soil_g_cm2)

  new_depth <- quantification_depth * mass1/mass2

  delta_d <- new_depth - quantification_depth


  soc_to_depth <- sample %>%
    filter(Lower_cm == quantification_depth) %>%
    pull(Cum_SOC_g_cm2) %>%
    max()

  if (soc_to_depth == -Inf){soc_to_depth <- 0}
  # Calculate the corrected total SOC stock Eq.5


  soc_adjust <- sample %>% filter(Upper_cm <= new_depth) %>%
    last() %>%
    summarize(delta_d  * SOC_pct/100 * BD_g_cm3 * adjust_factor) %>%
    ungroup() %>% pull()

  #print(soc_upper_layers)
  #print(last_layer_soc)

  SOC_Stock_g_cm2 <- soc_to_depth + soc_adjust

  # Calculate the corrected total soil mass
  Min_mass_soil_g_cm2 <- intial_soil %>%
    filter(Lower_cm == quantification_depth) %>%
    pull(Cum_Min_Soil_g_cm2)


  # Collect the corrected soil values
  temp = as.data.frame(
    list("Min_mass_soil_g_cm2" = Min_mass_soil_g_cm2,
         "Cum_SOC_g_cm2"= SOC_Stock_g_cm2,
         'depth' = new_depth))

  return(temp)
}


linear_interpolation_CUM_soc <- function(soil_sample, interpolation_mass){
  above_layers <- soil_sample %>%
    filter(Cum_Min_Soil_g_cm2 <= interpolation_mass)

  soc_to_depth <-  above_layers %>%
    pull(Cum_SOC_g_cm2) %>%
    max()

  if (soc_to_depth == -Inf){soc_to_depth <- 0}


  mass_to_depth <- above_layers %>%
    pull(Cum_Min_Soil_g_cm2) %>%
    max()
  if (mass_to_depth == -Inf){mass_to_depth <- 0}

  next_layer_SOC <- soil_sample %>%
    filter(Cum_Min_Soil_g_cm2 > interpolation_mass) %>%
    pull(SOC_pct) %>%
    first()

  if (length(next_layer_SOC) ==0){
    next_layer_SOC <- above_layers %>% pull(SOC_pct) %>% last()
  }

  #partial mass is denominated in Total Soil, not mineral soil,
  #so that SOC_pct can be multiplied against it
  partial_mass <- (interpolation_mass - mass_to_depth) / (1 - (next_layer_SOC/100))

  return(soc_to_depth + partial_mass * next_layer_SOC /100)


}



group_run_MC <- function(filtered_FD, adjustment_factor = 1,
                         quantification_depth = 30){
  t2.samps <- filter(filtered_FD, ID != Ref_ID)
  t1.samps <- filter(filtered_FD, ID == Ref_ID)

  fn <- function(sample){
    intial_soil <- subset(t1.samps,
                          Ref_ID==first(sample$Ref_ID) &
                            Rep==first(sample$Rep))
    return (
      linear_mass_correction( intial_soil, sample,
                              adjust_factor = adjustment_factor,
                              quantification_depth = quantification_depth))
  }


  return (t2.samps %>%
            group_by(ID, Ref_ID, Rep) %>%
            group_modify(~fn(.), .keep = T))}

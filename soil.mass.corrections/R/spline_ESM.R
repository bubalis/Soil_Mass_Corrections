
#
#spline methods for Equivalent Soil Mass
# Mostly drawn from Von Haden et al 2020
#

#source('util.R')
library(openxlsx)

#'Data validation from the VanHaden et al paper
validate.fd.data <- function(raw_fd){
  # Check input for basic errors ----
  is_error <- FALSE
  required_colnames <- c("ID", "Ref_ID", "Rep", "Upper_cm", "Lower_cm",
                         "SOC_pct", "SOM_pct", "BD_g_cm3")
  if (!all(required_colnames %in% colnames(raw_FD))){
    is_error <- TRUE
    stop(gettextf("Missing or misspelled column name(s)"))
  }
  if (any(is.na(raw_FD$ID) | is.na(raw_FD$Ref_ID) | is.na(raw_FD$Rep))){
    is_error <- TRUE
    stop(gettextf("Missing ID, Rep, or Ref_ID values"))
  }
  if (any(is.na(raw_FD$Upper_cm) | is.na(raw_FD$Lower_cm))) {
    is_error <- TRUE
    stop(gettextf("Missing Upper_cm or Lower_cm values"))
  }
  if (any(!is.numeric(raw_FD$Upper_cm) | !is.numeric(raw_FD$Lower_cm))) {
    is_error <- TRUE
    stop(gettextf("Non-numeric Upper_cm or Lower_cm values"))
  }
  if (any(!is.numeric(raw_FD$SOC_pct) | !is.numeric(raw_FD$BD_g_cm3) |
          !is.numeric(raw_FD$SOM_pct))){
    is_error <- TRUE
    stop(gettextf("Non-numeric SOC_pct, BD_g_cm3, or SOM_pct values"))
  }
  # Stop script if there is an error with the input file ----
  if (is_error){
    warning(gettextf("ESM script has failed due to error(s) listed above"))
  }

}

#'Preprocessing for the data provided with the VanHaden Paper for ESM analysis
preprocess.fd.data <- function(raw_FD){
  #Process from Van Haden et al (2022) for prepping raw soil core data for use
  #in the ESM calculations
  min_core_length_cm <- 0
  # Remove extra columns and add FD type
  reduced_FD <- subset(raw_FD, select=c(ID, Rep, Ref_ID, Upper_cm, Lower_cm,
                                        SOC_pct, SOM_pct, BD_g_cm3))
  reduced_FD$Type <- "FD"

  # Order columns
  reduced_FD <- reduced_FD[, c("Type", "ID", "Rep", "Ref_ID", "Upper_cm",
                               "Lower_cm", "SOC_pct", "SOM_pct", "BD_g_cm3")]

  # Remove cores that do not have a surface (0 cm) sample
  filtered_FD <- reduced_FD %>% group_by(ID, Rep) %>%
    filter(min(Upper_cm) == 0) %>%
    drop_na %>%
    arrange(ID, Rep, Upper_cm, Lower_cm) #



  #Remove any samples that are below a zone of non-contiguity
  all_contiguous_FD <- filtered_FD %>%
    group_by(ID, Rep) %>%
    filter(all(Upper_cm==dplyr::lag(Lower_cm, default=FALSE)))

  noncontiguous_FD <- filtered_FD %>%
    group_by(ID, Rep) %>%
    filter(any(Upper_cm!=dplyr::lag(Lower_cm, default=FALSE)))

  removed_noncontiguous_FD <- noncontiguous_FD %>%
    group_by(ID, Rep) %>%
    filter(Upper_cm < Upper_cm[which(Upper_cm!=dplyr::lag(Lower_cm,
                                                          default=FALSE))])

  filtered_FD <- rbind(all_contiguous_FD, removed_noncontiguous_FD)

  # Remove cores that do not have a sample deeper than min_core_length_cm
  filtered_FD <- filtered_FD %>%
    group_by(ID, Rep) %>%
    filter(!max(Lower_cm) < min_core_length_cm)

  # Throw an error if there are no observations left in dataset
  if (nrow(filtered_FD)==0){
    stop(gettextf("No observations remaining in dataset"))
  }
  return(filtered_FD)
}


hyman_extrapolating_spline <- function(xin, yin, xout, deriv = 0){
  fun <- yman_extrapolating_splinefun(xin, yin)

  return (fun(xout, deriv))

}

hyman_extrapolating_splinefun <- function(xin, yin){
  fun <- splinefun(xin, yin, method = 'hyman')
  fun_to_return <- function(xout, deriv){
    yout <- rep(NA, length(xout))
    yout[xout <= max(xin)] <- fun(xout[xout <= max(xin)], deriv = deriv)
    deriv_at_max <- fun(max(xin), deriv = 1)
    if (deriv == 0){
      yout[xout> max(xin)] <- fun(max(xin)) + (xout[xout > max(xin)] - max(xin)) * deriv_at_max

    }else if (deriv == 1){
      yout[xout> max(xin)] <- deriv_at_max

    }
    return (yout)
  }
  return(fun_to_return)
}


#' calculate estimates for cumulative SOC and other properties using the ESM approach
#' From Van Haden et al (2022), "Soil's Dirty Little Secret"  doi:10.1111/GCB.15124
#' @param filtered_FD A dataframe containing the fixed-depth samples
#'  each row of the dataframe represents a single layer of a soil core
#'  the dataframe must contain the following columns:
#'
#'  ID : a sample ID for the soil core.
#'  Ref_ID : the sample ID of the core used as reference for the ESM estimate.
#'  if the soil core is original, the Ref_ID is the same as the ID, otherwise it is the
#'  core taken from same location at the baseline
#'  Upper_cm : the starting depth of the soil layer
#'  Lower_cm : the ending depth of the soil layer
#'  SOC_pct : the organic carbon content (%) of the soil layer
#'  SOM_pct : the organic matter content (%) of the soil layer
#'  BD_g_cm3 : the bulk density (g/cm^3) of the soil layer.
#'
#' @param extrapolation A bool: should the method try to extrapolate outside the range for which there is data?
#' @param ESM_depths A vector. the depth layers to calculate cumulative masses to
#' @param return_fd A bool: should the fixed-depth data also be returned? Default is false
#' @returns A list: Cumulative_ESM: the dataframe with the cumulative soil properties added to it
#' Cumulative_FD : a dataframe of cumulative values to different depths using the Fixed-depth approach
calc_ESM_VanHaden <- function(filtered_FD, extrapolation, ESM_depths_cm, return_fd = F){

  # Add in zero masses at zero cm (for interpolation of first interval)
  modified_FD <- filtered_FD %>%
    group_by(ID, Ref_ID, Rep) %>%
    summarise() %>%
    mutate(Type = "FD", Upper_cm=0, Lower_cm=0,
           SOC_pct=0, SOM_pct=0, BD_g_cm3=0) %>%
    bind_rows(filtered_FD, .) %>%
    arrange(ID, Ref_ID, Rep, Upper_cm, Lower_cm)

  cumulative_FD <- calc.cumulative_masses(modified_FD)

  # Begin ESM-based calculations
  cumulative_ESM <- data.frame()

  for (i in 1:nrow(distinct(cumulative_FD, ID, Ref_ID, Rep))){
    current_vals <- distinct(cumulative_FD, ID, Ref_ID, Rep)[i,]
    current_Rep <- subset(cumulative_FD, ID==current_vals$ID &
                            Ref_ID==current_vals$Ref_ID &
                            Rep==current_vals$Rep)

    if (current_vals$ID != current_vals$Ref_ID){
      # Subset the reference set of values
      current_refs <- subset(cumulative_FD, ID==current_vals$Ref_ID)

      # Average the reference values (in case of multiple values per depth)
      current_refs_mean <- current_refs %>% group_by(Upper_cm, Lower_cm) %>%
        filter(Lower_cm %in% ESM_depths_cm) %>%
        mutate_at(vars(-Upper_cm, -Lower_cm, -Cum_Min_Soil_g_cm2),
                  function(x) x = NA) %>%
        mutate(Cum_Min_Soil_g_cm2 = mean(Cum_Min_Soil_g_cm2, na.rm=TRUE)) %>%
        summarise_all(mean) %>%
        mutate(ID=current_vals$ID, Ref_ID=current_vals$Ref_ID,
               Rep=current_vals$Rep, Type="ESM")

      #Determine whether extrapolation outside of maximum mass occurs
      if (extrapolation == FALSE) {
        # Remove references where mineral mass is greater than the sample max
        # Completely avoids extrapolation outside of spline model
        current_refs_filtered <-
          current_refs_mean[which(current_refs_mean$Cum_Min_Soil_g_cm2
                                  <= max(current_Rep$Cum_Min_Soil_g_cm2)),]
      } else {
        # Remove references that have a depth greater than sample max
        # Extraoplates only to the maximum depth of the samples
        current_refs_filtered <-
          current_refs_mean#[which(current_refs_mean$Lower_cm <=
        #          max(current_Rep$Lower_cm)),]
      }

      # Interpolate SOC and SOM using cubic spline models
      current_refs_filtered$Cum_SOC_g_cm2 <-
        spline(x=current_Rep$Cum_Min_Soil_g_cm2,
               y=current_Rep$Cum_SOC_g_cm2,
               xout=current_refs_filtered$Cum_Min_Soil_g_cm2,
               method="hyman")$y

      current_refs_filtered$Cum_SOM_g_cm2 <-
        spline(x=current_Rep$Cum_Min_Soil_g_cm2,
               y=current_Rep$Cum_SOM_g_cm2,
               xout=current_refs_filtered$Cum_Min_Soil_g_cm2,
               method="hyman")$y

      # Calculate non-cumulative masses
      current_refs_final <- current_refs_filtered %>%
        group_by(ID, Ref_ID, Rep) %>%
        mutate(Min_Soil_g_cm2 =
                 Cum_Min_Soil_g_cm2-dplyr::lag(Cum_Min_Soil_g_cm2, default=0),
               SOC_g_cm2 = Cum_SOC_g_cm2-dplyr::lag(Cum_SOC_g_cm2, default=0),
               SOM_g_cm2 = Cum_SOM_g_cm2-dplyr::lag(Cum_SOM_g_cm2, default=0))

      current_refs_final$Soil_g_cm2 <-
        current_refs_final$Min_Soil_g_cm2 + current_refs_final$SOM_g_cm2

      current_refs_final$BD_g_cm3 <-
        current_refs_final$Soil_g_cm2/
        (current_refs_final$Lower_cm-current_refs_final$Upper_cm)

      current_refs_final$SOC_pct <-
        current_refs_final$SOC_g_cm2/current_refs_final$Soil_g_cm2*100

      current_refs_final$SOM_pct <-
        current_refs_final$SOM_g_cm2/current_refs_final$Soil_g_cm2*100

      current_refs_final$Cum_Soil_g_cm2 <-
        cumsum(current_refs_final$Soil_g_cm2)

      current_ESM <- data.frame(current_refs_final)
      cumulative_ESM <- rbind(cumulative_ESM, current_ESM)
    }
  }

  # Post-processing cleanup and output ----
  # Remove zero masses at zero depth ESM and FD
  cumulative_FD <- subset(cumulative_FD, !(Upper_cm == 0 & Lower_cm == 0))
  cumulative_ESM <- subset(cumulative_ESM, !(Upper_cm == 0 & Lower_cm == 0))

  # Add NA for FD reference ID
  cumulative_FD$Ref_ID <- NA

  # Re-order column names in ESM dataset
  cumulative_ESM <- cumulative_ESM[colnames(cumulative_FD)]
  if (return_fd) {
    return (list('ESM' = cumulative_ESM, "FD" = cumulative_FD))
  }else{return (cumulative_ESM)}
}




run_ESM_VanHaden <- function(raw_FD, extrapolation, ESM_depths_cm){
  validate.fd.data(raw_fd)
  filtered_FD <- preprocess.fd.data(raw_fd)

  return (calc_ESM_VanHaden(filtered_FD, extrapolation, ESM_depths_cm))

}

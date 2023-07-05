library(dplyr)
library(tidyr)
library(data.table)
library(openxlsx)

library(soilDB)

library(sf)




#Convert the Vanhaden Format sample_dt into the Fowler format
VanHadenFD.from.Fowler.dt <- function(sample_dt){
  as_fd <- sample_dt %>% 
    mutate(Upper_cm = depth - depth_inc, SOM_pct = soc * 1.9, BD = BD * 10) %>%
    rename(SOC_pct = soc, BD_g_cm3 = BD, Lower_cm = depth) %>%
    select(ID, Rep, Ref_ID, Upper_cm, Lower_cm, SOC_pct, SOM_pct, BD_g_cm3, interval,
           scenario )
  
  
  
  
  return(as_fd)
}


#' aggregate a dataframe of soil_data to a given depth.
agg_to_depth <- function(soil_data, depth_to_agg){
  
  to_depth <- soil_data %>% filter(Lower_cm <= depth_to_agg)
  depths <- to_depth$Lower_cm -to_depth$Upper_cm
  soc_weights <- depths * to_depth$BD_g_cm3
  agged <- to_depth %>% summarize(SOC_pct = weighted.mean(SOC_pct, soc_weights), 
                                  SOM_pct = weighted.mean(SOM_pct, soc_weights), 
                                  BD_g_cm3 = weighted.mean(BD_g_cm3, depths),
                                  Upper_cm = min(Upper_cm),
                                  Lower_cm = max(Lower_cm),
                                  #Ref_ID = first(Ref_ID), 
                                  #ID = first(ID)
                                  )
  
  return (agged)
  
}

#'aggregate the soil data to several depths, to simulate soil-sampling
#'for instance: if new_depths = c(30, 50) simulates soil sampling to the depth of 
#'30cm and 50cm 
agg_to_depths <- function(soil_data, new_depths){
  if (!all(new_depths %in% soil_data$Lower_cm)){
    stop(gettextf("Error, can only aggregate to soil depths found in existing sample"))
  }
  new_tops <- lag(new_depths) %>% replace_na(0)
  
  new_soil_data <-data.frame()
  for (i in seq(1:length(new_depths))){
    filtered <- filter(soil_data, Upper_cm >= new_tops[i])
    new_soil_data <- rbind(new_soil_data, 
                           agg_to_depth(filtered, new_depths[i]))
    
  }
  return (new_soil_data)
}

#' Turn a dataframe of soil cores split into depths into a dataframe where these depths have been 
#' aggregated to give the results for a hypothetical sampling event with fewer soil depths
df_agg_to_depths <- function(soil_data, new_depths){
  
  return (soil_data %>% group_by(ID, Rep, Ref_ID) %>% 
            group_modify(~ agg_to_depths(.,  new_depths = new_depths), .keep = T 
          ))
  
}


get_SSurgo_spline_ESM<- function(lat, lon){
  DT_sf <- data.table(longitude= lon,
                      latitude = lat) %>%
    st_as_sf(coords = c("longitude", "latitude"), 
             crs = 4326, agr = "constant")
  
  res <- SDA_spatialQuery(
    DT_sf,
    what = "mupolygon",
    geomIntersection = T,
    #byFeature = T,
    db = "SSURGO",
    idcol = "gid",
    query_string = FALSE,
    as_Spatial = T
  ) %>% data.frame()
  
  vals <- paste(unique(res$mukey), sep = ' , ' , collapse = ",")
  
  q2 <- paste0("SELECT compname, om_h, om_l, om_r, dbovendry_l,  hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h,  dbovendry_r, dbovendry_h, mukey, comppct_r,  hzname  FROM component 
    INNER JOIN chorizon ch ON ch.cokey = component.cokey
    WHERE component.mukey IN (", vals, ")"
  )
  
  r2 <- SDA_query(q2)
  
  #here we set up data to make a spline on the SSurgo data
  # filter to the most common component (in the future, incorporate all components for uncertainty)
  # explore uncertainty by using the high and low values in addition to the means? 
  dt <- r2 %>% filter(comppct_r == max(r2$comppct_r)) %>% arrange(hzdept_r) %>% 
    mutate(cumulative_soil_mass = cumsum(dbovendry_r * (hzdepb_r - hzdept_r) )) %>% 
    mutate(cumulative_soc = cumsum(dbovendry_r * (hzdepb_r - hzdept_r) * om_r/1.9/100 )) %>% 
    mutate(cumulative_min_soil = cumulative_soil_mass - (cumulative_soc * 1.9))
  
  fitted.spline <-  splinefun(x=c(0, dt$cumulative_min_soil), 
                              y=c(0, dt$cumulative_soc),
                              
                              method="hyman")
  return (fitted.spline)  
}


#' Return a mass-correction factor for a given soil.
#' The mass-correction factor reflects the best-guess of the ratio of SOC in the
#' depth change layer to equivalent mass and the SOC in the whole 30-cm profile. 
#' 
Ssurgo_spline_MC <- function(cum_min_mass_t0,
                             cum_min_soil_t1, 
                             cum_SOC_to_depth_t1, 
                             fitted.spline){
  
  xout <- c(cum_min_soil_t0, cum_min_soil_t1)
  
  y <-  fitted.spline(xout, deriv = 0)
  # SOC % of mass to measurement depth, reference soil:
  soc_content_depth_layer_2 <- y[2] / cum_min_soil_t0
  
  # SOC % of mass in depth change layer, reference soil:
  soc_content_change_layer <- (y[2]- y[1]) / (cum_min_soil_t1 - cum_min_soil_t0)
  
  adjustment_factor <- soc_content_change_layer / soc_content_depth_layer_2
  print(adjustment_factor)
  
  return(adjustment_factor)
}

#'Run a mass-correction on single-depth measurements.
#'Optional adjustment factor which is an estimate for:
#'ratio of 
#'carbon percent in the adjustment_layer / carbon percent in total depth of core.
#'
mass_correction_single_depth <- function(cum_min_soil_t0,
                            cum_min_soil_t1, 
                            cum_SOC_to_depth_t1, 
                            adjustment_factor){
  
  mc_factor <- 1 + ( (cum_min_soil_t0 /cum_min_soil_t1) -1)* adjustment_factor
  print(mc_factor)
  return(cum_SOC_to_depth_t1 * mc_factor)
  
}





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
    filter(min(Upper_cm) == 0) %>% drop_na %>% # Remove samples that have NAs
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


#'Calculate layer and cumulative masses for soil data.
#'
calc.cumulative_masses <- function(soil_data){
  # Calculate soil mass in each interval
  return (
  soil_data %>% mutate(Soil_g_cm2 = (Lower_cm - Upper_cm) * BD_g_cm3) %>%
    #SOC and SOM mass of intervals
    mutate(SOC_g_cm2 = (SOC_pct/100)*Soil_g_cm2,
           SOM_g_cm2 = (SOM_pct/100)*Soil_g_cm2) %>%
    #mineral mass of intervals
    mutate(Min_Soil_g_cm2 = Soil_g_cm2 - SOM_g_cm2) %>%
    group_by(ID, Rep) %>%
    mutate(Cum_Soil_g_cm2 = cumsum(Soil_g_cm2),
           Cum_SOC_g_cm2 = cumsum(SOC_g_cm2),
           Cum_SOM_g_cm2 = cumsum(SOM_g_cm2),
           Cum_Min_Soil_g_cm2 = cumsum(Min_Soil_g_cm2))
  )
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
        current_refs_mean[which(current_refs_mean$Lower_cm <=
                                  max(current_Rep$Lower_cm)),]
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
    begin_last_depth <- sample %>% filter(Upper_cm < quantification_depth) %>% pull(Upper_cm) %>%
      max() 
    adjust_factor <- quantification_depth / (quantification_depth + begin_last_depth) * .5
  }
  
  sample <- calc.cumulative_masses(sample)
  intial_soil <- calc.cumulative_masses(intial_soil)
  
  
  mass1 <- intial_soil %>% filter(Lower_cm == quantification_depth) %>% pull(Cum_Min_Soil_g_cm2)
  mass2 <- sample %>% filter(Lower_cm == quantification_depth) %>% pull(Cum_Min_Soil_g_cm2)
  
  new_depth <- quantification_depth * mass1/mass2
  
  delta_d <- new_depth - quantification_depth
  
  
  soc_to_depth <- sample %>% filter(Lower_cm == quantification_depth) %>% 
    pull(Cum_SOC_g_cm2) %>% max()
  if (soc_to_depth == -Inf){soc_upper_layers <-0}
  # Calculate the corrected total SOC stock Eq.5
  
  
  soc_adjust <- sample %>% filter(Upper_cm <= new_depth) %>% last() %>%
    summarize(delta_d  * SOC_pct/100 * BD_g_cm3 * adjust_factor) %>%
    ungroup() %>% pull()
  
  #print(soc_upper_layers)
  #print(last_layer_soc)
  
  SOC_Stock_g_cm2 <- soc_to_depth + soc_adjust
  
  # Calculate the corrected total soil mass
  Min_mass_soil_g_cm2 <- intial_soil %>% filter(Lower_cm == quantification_depth) %>% 
    pull(Cum_Min_Soil_g_cm2)
  
  
  # Collect the corrected soil values
  temp = as.data.frame(list("Min_mass_soil_g_cm2" =        Min_mass_soil_g_cm2,
                            "Cum_SOC_g_cm2"= SOC_Stock_g_cm2, 
                            'depth' = new_depth))
  
  return(temp)
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
  
  
  return (t2.samps %>% group_by(ID, Ref_ID, Rep) %>% group_modify(~fn(.), .keep = T))}
  
dif_matrix <- function(v){
  m <- matrix(outer(v, v, "-"), nrow = length(v))
  diag(m) <- NA
  return(m)}
  

run_SSurgo_Mass_Corr <- function(lat, lon, data){
  
  
  fitted.spline <- get_SSurgo_spline_ESM(lat = lat,
                                         lon = lon)
  
  one_depth <- data %>% df_agg_to_depths(c(30)) %>% calc.cumulative_masses()
  Ids <- one_depth %>% filter(ID != Ref_ID) %>% pull(ID) 
  
  use_ssurgo_for_MC <- function(ID_){
    .Ref_ID <- filter(one_depth, ID == ID_) %>% pull(Ref_ID) %>% first()
    cum_min_soil_t0 <- filter(one_depth, ID == .Ref_ID) %>% pull(Cum_Min_Soil_g_cm2)
    cum_min_soil_t1 <- filter(one_depth, ID == ID_) %>% pull(Cum_Min_Soil_g_cm2)
    cum_soc_to_depth_t1 <- filter(one_depth, ID == ID_) %>% pull(Cum_SOC_g_cm2)
    
    adj_factor <- Ssurgo_spline_MC( cum_min_soil_t0,
                                    cum_min_soil_t1,
                                    cum_soc_to_depth_t1,
                                    fitted.spline)
    print(adj_factor)
    return(mass_correction_single_depth(cum_min_soil_t0,
                           cum_min_soil_t1,
                           cum_soc_to_depth_t1,
                           adj_factor))
    
  }
  
  
  
  res <- sapply(Ids, use_ssurgo_for_MC)
  MC_SSurgo <- data.frame(Cum_SOC_g_cm2 =res, ID = names(res),
                          method = 'Mass Correction, SSurgo',
                          sample_depths = 30, Rep = 1, Cum_SOC_g_cm2_baseline = NA,
      
                                        soc_change = NA)
  return (MC_SSurgo)
}


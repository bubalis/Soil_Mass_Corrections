


#
#SSurgo-Based Methods for Equivalent Soil-Mass Estimates
#
#
#source('util.R')
library(soilDB)
library(sf)
library(R.cache)


nan_or_null <- function(x){
  if (length(x)>0){return (is.na(x))
  }else{
    return(is.null(x))}
}



#'Retrieve Soil Profile data from ssurgo, given either a muname or lat and lon
get_ssurgo_depth_dat <- function(lat = NULL, lon = NULL, muname =NULL){
  if (!nan_or_null(lat) & !nan_or_null(lon)){
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
      #idcol = "gid",
      query_string = FALSE,
      as_Spatial = T
    ) %>% data.frame()


  }else if (!nan_or_null(muname)){
    q <-  paste0("SELECT muname, mukey FROM mapunit
    WHERE muname = '",muname, "'" )
    res <- SDA_query(q)
  }else{
    stop('Must pass either a lat and lon or a muname to retreive SSURGO data')}

  vals <- paste(unique(res$mukey), sep = ' , ' , collapse = ",")

  q2 <- paste0("SELECT compname, om_h, om_l, om_r, dbovendry_l,  hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h,  dbovendry_r, dbovendry_h, mukey, comppct_r,  hzname  FROM component
    INNER JOIN chorizon ch ON ch.cokey = component.cokey
    WHERE component.mukey IN (", vals, ")"
  )

  r2 <- SDA_query(q2)

  #here we set up data to make a spline on the SSurgo data
  # filter to the most common component

  mukey_to_keep <- r2 %>%
    dplyr::filter(!is.na(comppct_r) & !is.na(dbovendry_r) &
                    comppct_r == max(r2$comppct_r, na.rm = T)) %>%
    dplyr::pull(mukey) %>%
    table()  %>%
    names() %>%
    dplyr::first()

  return( r2 %>%
            dplyr::filter(comppct_r == max(r2$comppct_r, na.rm = T) &
                            mukey == mukey_to_keep & !is.na(dbovendry_r)) %>%
            arrange(hzdept_r) %>%
            mutate(
              cumulative_soil_mass = cumsum(dbovendry_r * (hzdepb_r - hzdept_r) ),
              cumulative_soc = cumsum(dbovendry_r * (hzdepb_r - hzdept_r) * om_r/1.9/100 ),
              hzdepb_r = hzdepb_r - min(hzdept_r),
              hzdept_r = hzdept_r - min(hzdept_r)) %>%
            mutate(
              cumulative_min_soil = cumulative_soil_mass - (cumulative_soc * 1.9)
            )

  )}


get_ssurgo_depth_dat <- addMemoization(get_ssurgo_depth_dat)


#'Parameterize the exponential decay SOC accumulation model based on
#'the soil profile data from SSURGO
get_SSurgo_exp_decay_params <- function(lat, lon, muname){
  d <- get_ssurgo_depth_dat(lat,lon, muname)
  xdata <-d$hzdepb_r
  ydata <- d$cumulative_soc
  m <- tryCatch({
    nls(ydata ~ -(a/b) * (exp(-b * xdata) -1),
        start =  list(a = .1, b = .1), control = list(maxiter = 5000))},
    error = function(e){nls(ydata ~ -(a/b) * (exp(-b * xdata) -1),
                            start =  list(a = .5, b = .1), control = list(maxiter = 500))})
  a <- summary(m)$coef[1,1]
  b <- summary(m)$coef[2,1]
  return (list(a = a, b =b))

}



get_SSurgo_spline_ESM<- function(lat, lon, muname){
  dt <- get_ssurgo_depth_dat(lat,lon, muname)

  fitted.spline <-  splinefun(x=c(0, dt$cumulative_min_soil),
                              y=c(0, dt$cumulative_soc),

                              method="hyman")
  return (fitted.spline)
}


#'Extrapolate the cumulative SOC value at depth x
#'Based on the decay parameters a and b.
decay_SOC_model <-function(x, a, b){
  -(a/b) * (exp(-b * x)-1)}


decay_ESM <- function(cum_min_soil_t0,
                      cum_min_soil_t1,
                      cum_SOC_to_depth_t1,
                      a, b,
                      depth_of_quant = 30){

  adj_depth <- cum_min_soil_t0/ cum_min_soil_t1 * depth_of_quant
  quant_depths <- c(depth_of_quant, adj_depth)
  yout <- decay_SOC_model(quant_depths, a,b)
  soc_ratio <- cum_SOC_to_depth_t1 / yout[1]

  soc_estimate <- soc_ratio * yout[2]

  return(soc_estimate)
}


#' Return a mass-correction factor for a given soil.
#' The mass-correction factor reflects the best-guess of the ratio of SOC in the
#' depth change layer to equivalent mass and the SOC in the whole 30-cm profile.
#'
Ssurgo_spline_MC <- function(cum_min_soil_t0,
                             cum_min_soil_t1,
                             cum_SOC_to_depth_t1,
                             fitted.spline){

  xout <- c(cum_min_soil_t0, cum_min_soil_t1)

  y <-  fitted.spline(xout, deriv = 0)
  # SOC % of mass to measurement depth, reference soil:

  soc_content_depth_layer_2 <- y[2] / cum_min_soil_t1

  # SOC % of mass in depth change layer, reference soil:
  soc_content_change_layer <- (y[2]- y[1]) / (cum_min_soil_t1 - cum_min_soil_t0)

  adjustment_factor <- soc_content_change_layer / soc_content_depth_layer_2
  #print(adjustment_factor)

  return(adjustment_factor)
}


mass_correction_single_depth <- function(cum_min_soil_t0,
                                         cum_min_soil_t1,
                                         cum_SOC_to_depth_t1,
                                         adjustment_factor){

  mc_factor <- 1 + ( (cum_min_soil_t0 /cum_min_soil_t1) -1)* adjustment_factor
  return(cum_SOC_to_depth_t1 * mc_factor)

}



run_SSurgo_Mass_Corr <- function(lat = NULL,
                                 lon = NULL, muname = NULL, data,
                                 depth_of_estimate = 30 ){


  fitted.spline <- get_SSurgo_spline_ESM(lat = lat,
                                         lon = lon,
                                         muname = muname)

  print(depth_of_estimate)

  one_depth <- data %>% df_agg_to_depths(c(depth_of_estimate)) %>%
    calc.cumulative_masses(override = T)

  Ids <- one_depth %>% filter(ID != Ref_ID) %>% pull(ID)


  use_ssurgo_for_MC_spline <- function(ID_){
    .Ref_ID <- filter(one_depth, ID == ID_) %>%
      pull(Ref_ID) %>%
      first()

    cum_min_soil_t0 <- filter(one_depth, ID == .Ref_ID) %>%
      pull(Cum_Min_Soil_g_cm2)

    cum_min_soil_t1 <- filter(one_depth, ID == ID_) %>%
      pull(Cum_Min_Soil_g_cm2)

    cum_soc_to_depth_t1 <- filter(one_depth, ID == ID_) %>%
      pull(Cum_SOC_g_cm2)

    adj_factor <- Ssurgo_spline_MC( cum_min_soil_t0,
                                    cum_min_soil_t1,
                                    cum_soc_to_depth_t1,
                                    fitted.spline)
    #print(adj_factor)
    return(mass_correction_single_depth(cum_min_soil_t0,
                                        cum_min_soil_t1,
                                        cum_soc_to_depth_t1,
                                        adj_factor))

  }





  res <- sapply(Ids, use_ssurgo_for_MC_spline)
  MC_SSurgo <- data.frame(Cum_SOC_g_cm2 =res, ID = names(res),
                          method = 'Mass Correction, SSurgo_spline',
                          sample_depths = depth_of_estimate, Rep = 1, Cum_SOC_g_cm2_baseline = NA,
                          soc_change = NA)

  exp_decay_factors <- get_SSurgo_exp_decay_params(lat, lon, muname)

  use_ssurgo_for_MC_EXP <- function(ID_){
    .Ref_ID <- filter(one_depth, ID == ID_) %>%
      pull(Ref_ID) %>% first()

    cum_min_soil_t0 <- filter(one_depth, ID == .Ref_ID) %>%
      pull(Cum_Min_Soil_g_cm2)

    cum_min_soil_t1 <- filter(one_depth, ID == ID_) %>%
      pull(Cum_Min_Soil_g_cm2)

    cum_soc_to_depth_t1 <- filter(one_depth, ID == ID_) %>%
      pull(Cum_SOC_g_cm2)

    return(decay_ESM(cum_min_soil_t0,
                     cum_min_soil_t1,
                     cum_soc_to_depth_t1,
                     a = exp_decay_factors$a, b = exp_decay_factors$b,
                     depth_of_quant = depth_of_estimate))

  }
  res2 <- sapply(Ids, use_ssurgo_for_MC_EXP)

  MC_ssurgo_2 <- data.frame(Cum_SOC_g_cm2=res2, ID = names(res2),
                            method = 'Mass Correction, SSurgo EXP',
                            sample_depths = depth_of_estimate, Rep = 1,
                            Cum_SOC_g_cm2_baseline = NA,
                            soc_change = NA)

  MC_ssurgo_3 <- data.frame(Cum_SOC_g_cm2=(res + res2)/2, ID = names(res2),
                            method = 'Mass Correction, SSurgo avg',
                            sample_depths = depth_of_estimate, Rep = 1,
                            Cum_SOC_g_cm2_baseline = NA,

                            soc_change = NA)

  return (rbind(MC_SSurgo, MC_ssurgo_2, MC_ssurgo_3))
}


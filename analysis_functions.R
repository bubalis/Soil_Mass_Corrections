library(devtools)

library(stringr)
load_all('soil.mass.corrections')

#' Make a single reference curve for all locations
#' Using the universal reference soil used in Fowler, 2023
make_standardized_spline <- function(){

  source('data_sim_fowler.R')

  dt_s2_t0 <- calc_soil (soc= soc_profile, soc_30= 1.4, soc_linear=FALSE,
                         BD_30=1.5, BD_linear=FALSE )

  sp_std <- splinefun(dt_s2_t0$Min_mass_soil_MG_ha/100, dt_s2_t0$SOC_Stock_MG_ha/100, method = 'hyman')
  return(sp_std)
}

sp_std <- make_standardized_spline()


#' Method for ESM estimate using a universal reference curve.
#' Using the universal reference soil used in Fowler, 2023
std_spline_mc <- function(data,
                                 depth_of_estimate = 30 ){



  one_depth <- data %>% df_agg_to_depths(c(depth_of_estimate)) %>%
    calc.cumulative_masses(override = T)

  Ids <- one_depth %>% filter(ID != Ref_ID) %>% pull(ID)


  spline_function <- function(ID_){
    .Ref_ID <- filter(one_depth, ID == ID_) %>%
      pull(Ref_ID) %>%
      first()

    cum_min_soil_t0 <- filter(one_depth, ID == .Ref_ID) %>%
      pull(Cum_Min_Soil_g_cm2)

    cum_min_soil_t1 <- filter(one_depth, ID == ID_) %>%
      pull(Cum_Min_Soil_g_cm2)

    cum_soc_to_depth_t1 <- filter(one_depth, ID == ID_) %>%
      pull(Cum_SOC_g_cm2)


    pred_known <- sp_std(cum_min_soil_t1)
    adj_factor <-  cum_soc_to_depth_t1 /pred_known


    #print(adj_factor)
    return(sp_std(cum_min_soil_t0) * adj_factor )

  }





  res <- sapply(Ids, spline_function)
  out <- data.frame(Cum_SOC_g_cm2 =res, ID = names(res),
                          method = 'Standardized Spline',
                          sample_depths = depth_of_estimate, Rep = 1, Cum_SOC_g_cm2_baseline = NA,
                          soc_change = NA)
  }



site_metadata <- read.csv('Source_data/site_metadata.csv')

mut_Ref_ID <- function(id){
  return (paste('Y5', substring(id, 4), sep = '_'))
}




#' Generate parameters for the beta distribution based on a mean and a variance
estBetaParams <- function(mu, var){

  alpha <- abs(( ((1-mu) /var) - 1/mu) * mu^2)
  beta <- ((1-mu)/var) * mu * (1-mu)
  return (list(alpha = alpha, beta = beta))
}




#' Automatically generate all depths to analyze with
#' the standard spline-based ESM methods.
#' Returns a list of every possible combination of 2 or more depths that includes the
#' the depth of quantification.
#' @param  depths - a vector of depths for core segements.
#' @param depth_to_quantify - the depth of interest for quantification.
select_ESM_depths <- function(depths, depth_to_quantify = 30){
  out <- sapply(seq(1, length(depths)),
                function(i){combn(depths, i, simplify = F)},
                simplify = 'list') %>% unlist(recursive = F)

  return(Filter(function(x){depth_to_quantify %in% x}, out ))
}


#' Automatically generate all depths to analyze with
#' linear mass-correction methods
#' Returns a list of every possible combination of 2 depths where one depth is
#' the depth of quantification
#'
#' below the depth of quantification,
#'
#' and the other is the depth of quantification.
#' @param  depths - a vector of depths for core segements.
#' @param depth_to_quantify - the depth of interest for quantification.
select_MC_depths <- function(depths, depth_to_quantify){
  #depths_to_choose <- depths[depths <= depth_to_quantify]

  return(
    append(
      Filter(function(x){depth_to_quantify %in% x},
             combn(depths, 2, simplify = F)),
      c(depth_to_quantify))
  )

}

#' Automatically generate all depths to analyze with
#' simulated a simulated split-core methods.
#' Returns a list of every possible combination of 2 depths where
#' one depth is <= the depth of quantification,
#' and the other depth is >= the depth of quantification
#'
#' @param  depths - a vector of depths for core segements.
#' @param depth_to_quantify - the depth of interest for quantification.
select_single_SOC_depths <- function(depths, depth_to_quantify){
  out <- combn(depths, 2, simplify = F)

  filter_fun <- function(x){
    cond1 <- max(abs(x - depth_to_quantify)) <= 30
    cond2 <- sum(x >= depth_to_quantify) >  0
    cond3 <- sum(x<= depth_to_quantify) > 0
    cond4 <- abs(x[1]-x[2]) <= 30

    return(all(cond1, cond2, cond3))
  }

  return (Filter(filter_fun, out ))
}


#'Return the string representing the maximum number of sample depths for
#'a given subset of the data.
get_max_depths <- function(data){
  return(data$sample_depths[which.max(nchar(data$sample_depths))])
}

#'Get the metadata/geographic data for a study
get_site_meta <- function(.study, raise.error = F){
  meta <- site_metadata %>%
    filter(study == .study)

  if (nrow(meta) == 0){
    message <- paste('Could not find metadata for', .study)
    if (raise.error){stop(message)
      }else{print(message)}
   }
  return(meta)
}


#' Read in data from the van-haden ESM paper
read.VanHaden.data <- function(){
  input_file_name <- "Source_Data/Example_datasets.xlsx"
  # Name of the sheet on the spreadsheet that contains the data
  input_file_sheet <- "a_temporal_paired"
  return(read.xlsx(input_file_name, sheet=input_file_sheet))
}


load_coefs <-function(){
  if (!file.exists('Results/error_model.rds')){return(list(-6, 1))

  }else {
      m = readRDS(file.path('Results', 'error_model.rds'))}
  return(coef(m))
}


make_error_propogation_func <- function(){
  coefs <- load_coefs()
  intercept <- coefs[[1]]
  .coef <- coefs[[2]]

  error_prop_func <- function(soil_data,x){
    err_stat <- exp(log((1/sqrt(sum_weights_mass(soil_data, x))))*.coef + intercept)
    return(err_stat * calc.linear_SOC_pct(soil_data, x))  }

  return(error_prop_func)
}

error_prop_func <- tryCatch(make_error_propogation_func(),
                            error = {function(x){return(function(x, y){NULL})}})

#use a re-scaled beta distribution to simulate errors
sim_benchmark_empirical_errs_beta <- function(soil_data,
                                              yest,
                                              mass_of_quant,

                                              nsims = 100){


  soil_data <- soil_data %>% bind_rows(make_zerodepth_row(.)) %>%
    arrange(Lower_cm) %>% ungroup()



  var_of_est <- error_prop_func(soil_data, mass_of_quant)**2

  #get parameters to rescale the distribution
  above_SOC <- soil_data %>% filter(Cum_Min_Soil_g_cm2 < mass_of_quant) %>%
    pull(Cum_SOC_g_cm2) %>%
    max()
  below_SOC <- soil_data %>% filter(Cum_Min_Soil_g_cm2 > mass_of_quant) %>%
    pull(Cum_SOC_g_cm2) %>%
    min()


  distance_from_last <- yest - above_SOC
  interval <- below_SOC - above_SOC

  param <- estBetaParams(distance_from_last/interval, #rescaled from 0-1
                         var_of_est/interval**2 # rescaled variance to the interval
  )

  scaled_back <- rbeta(nsims, param$alpha, param$beta)*interval + above_SOC

  return (scaled_back )
}


#'Simulate the benchmark "correct" SOC to the equivalent mass using the generated empirical error factor.
#' @param soil_data - a dataframe representing the soil measurements, to all depths
#' @param yest - the estimated cumulative SOC to the equivalent mass, generated by using a
#' spline-based ESM on all available depth-data
#' @param err_func a function to calcuate the standard deviation of the the errors, based on the total observation weight.
#' @param nsims number of synthetic "True" values to generate.
#'
sim_benchmark_empirical_errs <- function(soil_data, yest, mass_of_quant, err_prop_func, nsims = 200){
  soil_data <- soil_data %>% bind_rows(make_zerodepth_row(.)) %>%
    arrange(Lower_cm)

  err_sd <- error_prop_func(soil_data, mass_of_quant)

  return(yest + rnorm(nsims, 0, err_sd))
}


#'Simulate the benchmark "correct" SOC to the equivalent mass using the generated empirical error factor.
#'Return the res_data dataframe with the sythetic "TRUE" observations added in as the benchmark value.
#' @param FD_Data - a dataframe representing all raw soil data
#' @param res_data - a dataframe with the results from simulating different methods.
#' @param depth_of_est the depth of the reference mass
#' @param err_func a function to calcuate the standard deviation of the the errors, based on the total observation weight.
#' @param nsims number of synthetic "True" values to generate.
sim_all_benchmarks_empirical_errs <- function(FD_data, res_data, depth_of_est,
                                              nsims = 200,
                                              err_prop_func = error_prop_func
){
  out <- data.frame()
  FD_data <- calc.cumulative_masses(FD_data, override = T)



  for (i in seq(1, nrow(distinct(res_data, ID, Ref_ID)))){

    .ID <- distinct(res_data, ID, Ref_ID)[i,] %>%
      pull(ID) %>%
      first()


    .Ref_ID <- distinct(res_data, ID, Ref_ID)[i,] %>%
      pull(Ref_ID) %>%
      first()

    sim_data <- filter(FD_data, ID == .ID) %>%
      calc.cumulative_masses()

    sub_data <- filter(res_data, ID == .ID & Ref_ID == .Ref_ID)
    max_depths <- sub_data %>%
      filter(nchar(sample_depths) == max(nchar(sample_depths))) %>%
      pull(sample_depths) %>%
      first()

    yest <- sub_data%>%
      filter(sample_depths == max_depths) %>%
      pull(Cum_SOC_g_cm2) %>% first()

    mass_of_quant <- FD_data %>%
      filter( ID == .Ref_ID & Lower_cm == depth_of_est) %>%
      pull(Cum_Min_Soil_g_cm2) %>%
      first()

    simmed <- sim_benchmark_empirical_errs_beta(sim_data, yest,
                                                mass_of_quant,
                                                #err_prop_func = err_func_prop_func,
                                                nsims = nsims)

    #print(.Ref_ID)
    #print(yest)
    new <- data.frame(SOC_benchmark = simmed, ID = .ID,
                      Ref_ID = .Ref_ID,
                      simulation.N = seq(1:nsims),
                      benchmark_no_errors = yest)

    out <- rbind(out, new)}
  return (out) }





#'compare the ESM and Mass-correction methods on a dataset.
#'@param fd_data a dataframe: data in the format used by VanHaden
#'@param ESM_depths a list of vectors: depth combinations to use the ESM method on (if null, these are auto-selected)
#'@param MC_depths a list of vectors: depth combinations to use the Mass-correction methods on. (if null, these are auto-selected)
#'@param split_core_single_SOC_depths a list of vectors: depth combinations to simulate a 2-depth core that is weighed
#'and then combined into a single SOC sample
#'@param depth_of_estimate: numeric: the depth we are trying get an accurate estimate of SOC for.
compare.methods <- function(fd_data,  ESM_depths = NULL, MC_depths = NULL,
                            split_core_single_SOC_depths = NULL,
                            depth_of_estimate = 30, lat = NULL, lon = NULL,
                            muname = NULL){

  #if depth combinations are not specified, generate them automatically:
  if (is.null(ESM_depths)){
    ESM_depths <- select_ESM_depths(unique(fd_data$Lower_cm), depth_of_estimate)}

  if (is.null(MC_depths)){
    MC_depths <- select_MC_depths(unique(fd_data$Lower_cm), depth_of_estimate)}

  if (is.null(split_core_single_SOC_depths)){
    split_core_single_SOC_depths <- select_single_SOC_depths(unique(fd_data$Lower_cm),
                                                             depth_of_estimate)}

  comparisons <- fd_data %>% select(ID, Rep, Ref_ID) %>% distinct()

  #fixed depth results
  single_depth_agged <- df_agg_to_depths(fd_data, c(depth_of_estimate))
  FD_res <-  single_depth_agged  %>%
    mutate(Cum_SOC_g_cm2 = Lower_cm * BD_g_cm3 * SOC_pct/100,
           method = 'Fixed Depth',
           sample_depths = as.character(depth_of_estimate)
    )

  mass_corr_avg_res <- single_depth_agged %>%
    group_run_MC(adjustment_factor = 'standard_correction',
                 quantification_depth = depth_of_estimate) %>%
    mutate(
      method = "Linear Average",
      sample_depths = as.character(depth_of_estimate))

  mass_corr_fowler_weights <- single_depth_agged %>%
    group_run_MC(adjustment_factor = 17/32,
                 quantification_depth = depth_of_estimate) %>%
    mutate(method = "Linear Average, Fowler Weights",
           sample_depths = as.character(depth_of_estimate))

  all_res <- rbind(FD_res,
                   mass_corr_avg_res,
                   mass_corr_fowler_weights)


  for (depth_vals in ESM_depths){
    print(depth_vals)
    sample <- df_agg_to_depths(fd_data %>%
                                 filter(ID != Ref_ID), depth_vals)

    initial <- df_agg_to_depths(fd_data %>%
                                  filter(ID == Ref_ID),
                                c(depth_of_estimate))

    data <- rbind(sample, initial)

    ESM_res <- data %>%
      calc_ESM_VanHaden(T, depth_vals) %>%
      filter(Lower_cm == depth_of_estimate) %>%
      mutate(sample_depths = paste(depth_vals, collapse = ', '),
             method = 'ESM Spline') %>% filter(ID != Ref_ID)
    all_res <- rbind(all_res, ESM_res)

    if (length(depth_vals) ==2){ # for 2-depth cores, try fitting a decay model exactly.

      ESM_decay_res <- run_decay_ESM(data, depth_of_estimate) %>%
        mutate(sample_depths = paste(depth_vals, collapse = ', '),
               method = 'Exponential Decay', Type = 'decay_ESM')

      all_res  <- rbind(all_res, ESM_decay_res)
    }}

  for (depth_vals in MC_depths){
    to_depths <-  rbind(df_agg_to_depths(fd_data %>% filter(ID != Ref_ID),
                                         depth_vals),
                        df_agg_to_depths(fd_data %>% filter(ID == Ref_ID),
                                         c(depth_of_estimate)))


    mass_correction_res <- to_depths %>%
      group_run_MC(adjustment_factor = 1,
                   quantification_depth = depth_of_estimate) %>%
      mutate(method = "Linear",
             sample_depths = paste(depth_vals,
                                   collapse = ', '))

    all_res <- rbind(all_res, mass_correction_res)
  }

  for (depth_vals in split_core_single_SOC_depths){
    res <- group_run_joined_sample(fd_data, depth_vals, depth_of_estimate) %>%
      mutate(method = 'Two-depth combined sample',
             sample_depths = paste(depth_vals, collapse = ', '))
    all_res <- rbind(all_res, res)

  }

  if (!all(is.na(c(lat, lon, muname)))){
    ssurgo_res <- run_SSurgo_Mass_Corr(lat = lat, lon = lon,
                                       muname = muname, data = fd_data,
                                       depth_of_estimate = depth_of_estimate) %>%
      mutate(sample_depths = as.character(sample_depths))

    std_spline_est <- std_spline_mc(fd_data, depth_of_estimate) %>%
      mutate(sample_depths = as.character(sample_depths))

    all_res <- rbind(all_res, ssurgo_res, std_spline_est)
  }

  all_res <-all_res %>% ungroup() %>%
    select(ID, Rep, Upper_cm, Lower_cm, Cum_SOC_g_cm2,method, sample_depths) %>%
    merge(comparisons)

  baseline_data <- all_res %>%
    ungroup() %>%
    filter( method == 'Fixed Depth' & Ref_ID == ID) %>%
    rename( Cum_SOC_g_cm2_baseline = Cum_SOC_g_cm2) %>%
    select(Ref_ID, Rep, Cum_SOC_g_cm2_baseline)

  res_changes <- all_res %>%
    filter(ID != Ref_ID) %>%
    merge(baseline_data) %>%
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline) %>%
    select(ID, Rep, Ref_ID, Cum_SOC_g_cm2, method, sample_depths,
           Cum_SOC_g_cm2_baseline,
           soc_change)

  return (res_changes)
}




#'Calculate the mass change to depth of quantification between the soils being quantified
#'And its reference soil
mass_change <- function(FD_data, quant_depth){
  cum_FD <- calc.cumulative_masses(FD_data)

  return (cum_FD %>%
            filter(ID != Ref_ID & Lower_cm == quant_depth) %>%
            merge(cum_FD %>% filter(Lower_cm == quant_depth) %>%
                    select(ID, Cum_Min_Soil_g_cm2),
                  by.x = 'Ref_ID', by.y = 'ID') %>%
            mutate(
            mass.change = Cum_Min_Soil_g_cm2.x -  Cum_Min_Soil_g_cm2.y) %>%
            select(Ref_ID, ID, mass.change))

}



#'Summarize all comparisons made, using an error simulation function to generate
#'Synthetic "TRUE" observations
#' Return the results dataframe with the synthetic true observations and their implied errors
#' added in.
comparison_summarize_sim <- function(res, FD_data, depth_of_estimate, nsims,
                                     simulation_function = sim_all_benchmarks_empirical_errs ){
  res <- res %>%
    arrange(ID) %>%
    tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>%
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline )

  simulated_benchmarks <- simulation_function( FD_data, res, depth_of_estimate, nsims)

  mass_changes <- bind_rows(
    lapply(unique(FD_data$ID),
           function(.ID){mass_change(mutate(FD_data, Ref_ID = .ID), depth_of_estimate)}))

  res <- res %>%
    merge(simulated_benchmarks, by = c('ID', 'Ref_ID')) %>%
    mutate(
      dif_from_benchmark = Cum_SOC_g_cm2 - SOC_benchmark,
      soc_change_benchmark = SOC_benchmark- Cum_SOC_g_cm2_baseline) %>%
    left_join(mass_changes)

  return(res)
}



#' Run comparisons on the Fixed-depth core data
run_comparison <- function(FD_data, ESM_depths = NULL, MC_depths = NULL,
                           split_core_single_SOC_depths = NULL,
                           depth_of_estimate = 30, study_name, .site = NULL
                           ){

  site_meta <- get_site_meta(study_name)
  all_out <- data.frame()


  #if a specific site is not specified, and each comparison is already specified
  if (is.null(.site) & ('Ref_ID' %in% colnames(FD_data))){
    lat <- site_meta %>% pull(lat) %>% first()
    lon <- site_meta %>% pull(lon) %>% first()
    muname <- site_meta %>% pull(muname) %>% first()

    res <- compare.methods(FD_data, ESM_depths = ESM_depths, MC_depths = MC_depths,
                           split_core_single_SOC_depths = split_core_single_SOC_depths,
                           lat = lat, lon = lon, muname = muname) %>%
      mutate(study_name = study_name)

    all_out <- res %>% mutate(1)
    #all_out <- comparison_summarize_sim(res, FD_data, 0, depth_of_estimate) %>%
    #  mutate(weight = 1)

  }else{
    if ('site' %in% colnames(FD_data)){
      #if a site is listed in the dataframe, only run comparisons between plots at the same site
      for (site_name in unique(FD_data$site)){

        new_res <- run_comparison(FD_data %>% filter(site == site_name) %>% select(-c(site)),
                                  ESM_depths = ESM_depths,
                                  MC_depths = MC_depths,
                                  split_core_single_SOC_depths = split_core_single_SOC_depths,
                                  depth_of_estimate = depth_of_estimate,
                                  study_name = study_name, .site = site_name)

        all_out <- rbind(all_out, new_res)
      }


    }else{
      #if a site argument was passed, get location data for specific site
      if (!is.null(.site)){
        site_dat <- filter(site_meta, site == .site)}

      else{site_dat <- site_meta}

      lat <- site_dat %>%
        pull(lat) %>%
        first()

      lon <- site_dat %>%
        pull(lon) %>%
        first()

      muname <- site_dat %>%
        pull(muname) %>%
        first()


      for (.ID in unique(FD_data$ID)){
        sub_fd <- FD_data %>% mutate(Ref_ID = .ID)
        res <- compare.methods(sub_fd,
                               ESM_depths = ESM_depths,
                               MC_depths = MC_depths,
                               split_core_single_SOC_depths = split_core_single_SOC_depths,
                               depth_of_estimate = depth_of_estimate,
                               lat = lat, lon = lon, muname = muname) %>%
          mutate(study_name = study_name) %>%
          #comparison_summarize_sim(sub_fd, .ID, depth_of_estimate) %>%
          mutate(weight = 1/length(unique(FD_data$ID)))

        all_out <- rbind(all_out, res)

      }   }}


    return (all_out)
  }


#'Run a before-after comparison on data, where both are present.
time_series_comparison <- function(fd_data, ESM_depths = NULL, MC_depths = NULL,
                                   split_core_single_SOC_depths = NULL,
                                   study_name, depth_of_estimate = 30){
  if (is.null(ESM_depths)){
    ESM_depths <- select_ESM_depths(unique(fd_data$Lower_cm), depth_of_estimate)}

  if (is.null(MC_depths)){
    MC_depths <- select_MC_depths(unique(fd_data$Lower_cm), depth_of_estimate)}

  if (is.null(split_core_single_SOC_depths)){
    split_core_single_SOC_depths <- select_single_SOC_depths(unique(fd_data$Lower_cm), depth_of_estimate)}

  out <- run_comparison(fd_data, ESM_depths = ESM_depths,
                        MC_depths = MC_depths,
                        split_core_single_SOC_depths = split_core_single_SOC_depths,
                        study_name = study_name)

  res <- out$res
  summary <- out$summary

  #summarize the between-treatment comparisons
  #Compare the calculated differences to the benchmark (ESM on all available data)

  calced_difs <- res %>% filter(ID != 'all_t1') %>% group_by(method, sample_depths) %>%
    reframe(changes = dif_matrix(soc_change))

  as.df <- as.data.frame(calced_difs$changes) %>%
    mutate(method = calced_difs$method, sample_depths = calced_difs$sample_depths)
  calced_difs <- as.df

  max_depths <- as.df %>%
    filter(nchar(sample_depths) == max(nchar(res$sample_depths))) %>%
    pull(sample_depths) %>% first()

  benchmark_difs <- as.df %>% filter(sample_depths == max_depths) %>%
    select(-c(method, sample_depths))


  differences_summary <- data.frame()

  dist <- distinct(calced_difs, method, sample_depths)
  for ( i in seq(1: nrow(dist))){
    .sample_depths <- dist[i, 'sample_depths']
    .method <- dist[i, 'method']
    subset <- calced_difs %>% filter(.sample_depths == sample_depths, .method == method
    )
    difs <- select(subset, -c(method, sample_depths)) - benchmark_difs
    perc_difs <- difs / benchmark_difs
    RMSE <- sqrt(mean(colMeans(difs **2, na.rm = T)))
    RMSE_perc <- sqrt(mean(colMeans(perc_difs **2, na.rm = T)))
    mean_absolute_perc_err <- mean(colMeans(abs(perc_difs), na.rm = T))
    mean_absolute_err <- mean(colMeans(abs(difs), na.rm = T))
    .bias <- mean(difs[upper.tri(difs)])

    differences_summary <- rbind(differences_summary,
                                 data.frame(method = .method, sample_depths = .sample_depths,
                                            RMSE_differences = RMSE, RMSE_differences_perc = RMSE_perc,
                                            mean_absolute_perc_err_dif = mean_absolute_perc_err,
                                            difs_err = mean_absolute_err,
                                            bias = .bias))
  }
  summary <- merge(summary, differences_summary, by = c('method', 'sample_depths') )

  return (list('summary' = summary, 'res' = res %>% mutate('study' = study_name)))}


#'Add in randomly generated benchmark or "True" values.
add_reference_changes <- function(.study_name, out_data, nsims=50,
                                  load_fd_data_func = load_study, depth_of_quant = 30){
  print(.study_name)
  sub_data <- filter(out_data, study_name == .study_name)
  FD_data <- load_fd_data_func(.study_name)
  IDs <- unique(sub_data$ID)

  FD_benchmark <- sub_data %>% filter(method == 'Fixed Depth') %>%
    rename(FD_val = Cum_SOC_g_cm2) %>%
    select(ID, Ref_ID, FD_val)
  sub_data <- sub_data %>% merge(FD_benchmark)

  out <- comparison_summarize_sim(sub_data, FD_data,  depth_of_quant, nsims = nsims)

  return( out %>% mutate(FD_error = FD_val - SOC_benchmark) %>%
            mutate(relative_error = abs(dif_from_benchmark / FD_error)) %>% select(-c(FD_val)))
}


#'Summarize a grouped results dataframe to give summary statistics
summarize_res <- function(grouped_res){
  if (!('weight' %in% colnames(grouped_res) )){
    grouped_res$weight <- 1
  }
  return(grouped_res %>% summarise(
    RMSE = sqrt(weighted.mean(dif_from_benchmark**2, weight)),
    RMPSE = sqrt(weighted.mean((dif_from_benchmark /SOC_benchmark )**2, weight )),
    bias = weighted.mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1), weight),
    #Tons carbon overestimated per ton mass.change in fixed-depth sample
    bias_rate = weighted.mean(dif_from_benchmark / mass.change, weight),
    MAPE_change = weighted.mean(abs(dif_from_benchmark/soc_change_benchmark), weight),
    Max_err =  max(abs(dif_from_benchmark)),
    mean_abs_err = weighted.mean(abs(dif_from_benchmark), weight)
  ))

}



#'Combine a list of 2 or more dataframes by binding their rows.
#'Respond to a memory error by combining them in a csv then loading them.
rbind_memoryerror <- function(li_of_dfs){
  out<- tryCatch({bind_rows(li_of_dfs)},
                 error = function(cond){
                   lapply(li_of_dfs, function(x){write.csv(x,'temp.csv', append = T)})
                   out <- read.csv('temp.csv')
                   file.remove('temp.csv')
                   return(out)
                 }

  )
  return (out)
}

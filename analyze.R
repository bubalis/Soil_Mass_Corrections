
source('ESM_and_mass_correction_functions.R')

site_metadata <- read.csv('Source_data/site_metadata.csv')


get_site_meta <- function(.study){
  return (site_metadata %>% filter(study == .study))
  
}


#' Read in data from the van-haden ESM paper
read.VanHaden.data <- function(){
  input_file_name <- "Source_Data/Example_datasets.xlsx"
  # Name of the sheet on the spreadsheet that contains the data
  input_file_sheet <- "a_temporal_paired"
  return(read.xlsx(input_file_name, sheet=input_file_sheet))
}



#'compare the ESM and Mass-correction methods on a dataset.
#'@param fd_data a dataframe: data in the format used by VanHaden
#'@param ESM_depths a list of vectors: depth combinations to use the ESM method on
#'@param MC_depths a list of vectors: depth combinations to use the Mass-correction methods on.
#'@param depth_of_estimate: numeric: the depth we are trying get an accurate estimate of SOC for.
compare.methods <- function(fd_data,  ESM_depths, MC_depths,
                            depth_of_estimate = 30, lat = NULL, lon = NULL,
                            muname = NULL){
  
  
  
  comparisons <- fd_data %>% select(ID, Rep, Ref_ID) %>% distinct()
  
  #fixed depth results
  single_depth_agged <- df_agg_to_depths(fd_data, c(depth_of_estimate))
  FD_res <-  single_depth_agged  %>% 
    mutate(Cum_SOC_g_cm2 = Lower_cm * BD_g_cm3 * SOC_pct/100,
           method = 'Single Depth',
           sample_depths = as.character(depth_of_estimate)
    )
  
  mass_corr_avg_res <- single_depth_agged %>% 
    group_run_MC(adjustment_factor = 'standard_correction',  quantification_depth = depth_of_estimate) %>% 
    mutate(method = "Mass Correction, SOC_corrected?", 
           sample_depths = as.character(depth_of_estimate))
  
  mass_corr_fowler_weights <- single_depth_agged %>% 
    group_run_MC(adjustment_factor = 17/32,  quantification_depth = depth_of_estimate) %>% 
    mutate(method = "Mass Correction, Fowler Weights", 
           sample_depths = as.character(depth_of_estimate))
  
  all_res <- rbind(FD_res, mass_corr_avg_res, mass_corr_fowler_weights)
  
  
  for (depth_vals in ESM_depths){
    ESM_res <-  rbind(df_agg_to_depths(fd_data %>% filter(ID != Ref_ID), depth_vals),
                          df_agg_to_depths(fd_data %>% filter(ID == Ref_ID), c(depth_of_estimate))
                                 )%>% 
      calc_ESM_VanHaden(T, depth_vals) %>% 
      filter(Lower_cm == depth_of_estimate) %>% #only keep the estimate for the depth we are calculating
      mutate(sample_depths = paste(depth_vals, collapse = ', '),
             method = 'ESM')
    all_res <- rbind(all_res, ESM_res)
  }
  
  for (depth_vals in MC_depths){
    print(depth_vals)
    to_depths <-  rbind(df_agg_to_depths(fd_data %>% filter(ID != Ref_ID), depth_vals),
                        df_agg_to_depths(fd_data %>% filter(ID == Ref_ID), c(depth_of_estimate))
    )
    
    mass_correction_res <- to_depths %>% 
      group_run_MC(adjustment_factor = 1, quantification_depth = depth_of_estimate) %>% 
      mutate(method = "Mass Correction", 
             sample_depths = paste(depth_vals, collapse = ', '))
    
   
    
    all_res <- rbind(all_res, mass_correction_res)
    
  }
  if (!all(is.na(c(lat, lon, muname)))){
    ssurgo_res <- run_SSurgo_Mass_Corr(lat = lat, lon = lon,
                                       muname = munam, data = fd_data,
                                       depth_of_estimate = depth_of_estimate)
    all_res <- rbind(all_res, ssurgo_res)
  }
  
  
  
  all_res <-all_res %>% ungroup() %>%select(ID, Rep, Upper_cm, Lower_cm, 
                                            Cum_SOC_g_cm2,
                                            method, sample_depths) %>% merge(comparisons)
  
  baseline_data <- all_res %>%ungroup() %>% filter( method == 'Single Depth' & Ref_ID == ID) %>% 
    rename( Cum_SOC_g_cm2_baseline = Cum_SOC_g_cm2) %>%
    select(Ref_ID, Rep, Cum_SOC_g_cm2_baseline)
  
  res_changes <- all_res %>% filter(ID != Ref_ID) %>% merge(baseline_data) %>%
    mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline) %>% 
    select(ID, Rep, Cum_SOC_g_cm2, method, sample_depths, Cum_SOC_g_cm2_baseline,
           soc_change)
  
  return (res_changes)
}


mut_Ref_ID <- function(id){
  return (paste('Y5', substring(id, 4), sep = '_'))
}


#'Calculate the mass change to depth of quantification
mass_change <- function(FD_data, quant_depth){
  cum_FD <- calc.cumulative_masses(FD_data) 
  
  return (cum_FD %>% filter(ID != Ref_ID & Lower_cm == quant_depth) %>% 
    merge(cum_FD %>% filter(Lower_cm == quant_depth) %>% select(ID, Cum_Min_Soil_g_cm2), 
          by.x = 'Ref_ID', by.y = 'ID') %>% 
    mutate(mass.change = Cum_Min_Soil_g_cm2.x -  Cum_Min_Soil_g_cm2.y) %>% 
    select(ID, mass.change))
  
}


#'Compare all individual treatment/core data to every possible ref
run_comparison <- function(FD_data, 
                          ESM_depths, MC_depths, 
                          depth_of_estimate = 30, 
                          study_name){
  
  site_meta <- get_site_meta(study_name)
  all_out <- data.frame()
  
  if ('Ref_ID' %in% colnames(FD_data)){
      lat <- site_meta %>% pull(lat) %>% first()
      lon <- site_meta %>% pull(lon) %>% first()
      muname <- site_meta %>% pull(lon) %>% first()
      
      res <- compare.methods(FD_data, ESM_depths = ESM_depths, MC_depths = MC_depths, 
                             lat = lat, lon = lon, muname = muname)
    
      
      res <- res %>% arrange(ID) %>% 
      tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
      mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline, Ref = 0 )
    res <- res %>% merge( mass_change(FD_data, depth_of_estimate), by = 'ID')
    all_out <- res
  }else{
    if ('site' %in% colnames(FD_data)){
      #if a site is listed in the data, only run comparisons between plots at the same site
      
      for (.site in unique(FD_data$site)){
        site_dat <- site_meta %>% filter(site == .site)
        lat <- site_dat %>% pull(lat) %>% first()
        lon <- site_dat %>% pull(lon) %>% first()
        muname <- site_dat %>% pull(lon) %>% first()
        new_res <- run_comparison(FD_data %>% filter(site == .site,) %>% select(-c(site)),
                                  ESM_depths = ESM_depths,
                                  MC_depths = MC_depths,
                                  depth_of_estimate = depth_of_estimate,
                                  study_name = study_name, lat = lat, lon = lon,
                                  muname = muname)$res
        all_out <- rbind(all_out, new_res)
      }
    }else{
  for (.ID in unique(FD_data$ID)){
    FD_data <- FD_data %>% mutate(Ref_ID = .ID) 
    
    res<- compare.methods(FD_data, ESM_depths = ESM_depths,
                          MC_depths = MC_depths, 
                          depth_of_estimate = depth_of_estimate)
    
    
    
    res <- res %>% arrange(ID) %>% 
      tidyr::fill(Cum_SOC_g_cm2_baseline, .direction = 'down') %>% 
      mutate(soc_change = Cum_SOC_g_cm2 - Cum_SOC_g_cm2_baseline, Ref = .ID )
    res <- res %>% merge( mass_change(FD_data, depth_of_estimate), by = 'ID')
    
    all_out <- rbind(all_out, res)
  }
  }}
  
  max_depths <- all_out %>% filter(nchar(sample_depths) == max(nchar(all_out$sample_depths))) %>% 
    pull(sample_depths) %>% first()
  
  all_out <- all_out %>% merge(all_out %>% filter(sample_depths == max_depths) 
                       %>% select(ID, soc_change, Ref) %>% 
                         rename(soc_change_benchmark = soc_change)) %>% 
    mutate(dif_from_benchmark = soc_change - soc_change_benchmark) %>% 
    left_join(mass_change(FD_data, depth_of_estimate)) %>% 
    mutate(study = study_name, weight = 1/ (nrow(distinct(FD_data, ID)) -1))
  
  
  summary<- all_out %>% group_by(method, sample_depths) %>% 
    summarize(RMSE_from_benchmark = sqrt(mean(dif_from_benchmark**2)), 
              Bias.from.benchmark = mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1)))
  return (list(summary = summary, res = all_out))
}


#Run comparisons

filtered_FD <- read.VanHaden.data() %>% preprocess.fd.data()
#filtered_FD <- filtered_FD %>% mutate(Ref_ID = sapply(ID, set_Ref_ID ))



vanHoden_res <- compare.methods(filtered_FD, 
                                ESM_depths = list(c(10, 30, 50, 100),  c(10, 30), c(30, 50),
                                                  c(10, 30, 50)),
                                MC_depths = list(c(10,30), c(30)))


#run the vanHoden data backwards, so that the "after" soil is the 
#reference for the "before" soil
reversed <- filtered_FD %>% mutate (Ref_ID = sapply(Ref_ID, mut_Ref_ID))                            
vanHoden_rev <- compare.methods(reversed, 
                                ESM_depths = list(c(10, 30, 50, 100),  c(10, 30), c(30, 50),
                                                  c(10, 30, 50)),
                                MC_depths = list(c(10,30), c(30)))

vanHoden_all <- rbind(vanHoden_res, vanHoden_rev)

vanHoden_all <- vanHoden_all %>%  merge(vanHoden_all %>% filter(method == 'Single Depth') 
                                        %>% select(ID, Rep, soc_change) 
                                        %>% rename(FD_error = soc_change)) %>% 
  mutate(relative_error = soc_change / abs(FD_error) ) %>% select( -c(FD_error)) %>% 
  mutate(error = soc_change)

vanHoden_all %>% group_by(method, sample_depths) %>% summarize(RMSE = sqrt(mean(soc_change**2)))


#compare the results of different methods on the fowler simulated data, and 
#variations on it
source("data_sim.R")
fowler_sim_data <- simulate.soil.Fowler(soc_profile)


#folwer_dt < - fowler_sim_data$fowler_dt
#sample_dt <- fowler_sim_data$sample_dt

t_dt <- fowler_sim_data$t_dt
soc_delt_actual <-  fowler_sim_data$true_delta_soc

soc_delt_actual <- soc_delt_actual %>% 
  mutate(true_change_soc =true_change_soc/100,
         Ref_ID = ifelse(grepl('t0', ID) & grepl('t2', Ref_ID), paste(Ref_ID, 'x', sep = '_'),
                         Ref_ID)
  ) %>% filter(!(grepl('t1', ID) & grepl('t2', Ref_ID))) %>% filter(Ref_ID != ID) %>%
  select(ID, true_change_soc)





t_dt <- t_dt %>% filter(scenario != 's1')

y0_rev <- t_dt %>% filter(period %in% c('t0', 't1')) %>%
  mutate(Ref_ID = paste(scenario, 't1', sep = '_'))

#reversed sample_dt, with t0 as the ending point, and t2 as the starting point
y0_rev2 <- t_dt %>% filter(period %in% c('t0', 't2')) %>%
  mutate(period = paste(period, 'x', sep = '_')) %>%
  mutate(Ref_ID = paste(scenario, 't2_x', sep = '_'),
         ID = paste(scenario, period, sep = '_'))


fowler.all <- data.frame()
for (sample_data in list(t_dt, y0_rev, y0_rev2))
   {
     ESM_depths <- list(c(10,20,30),  c(15,30), c(24,30), c(24,30, 36),
                        c(30, 36), c(15, 30,50), c(10, 30, 50), c(30, 50), 
                        c(10, 30, 50, 70), c(15, 30, 45, 60, 75))
     MC_depths <- list(c(10,20,30),  c(15,30), c(24,30), c(10, 30), c(30))
     fowler.sim     <- compare.methods(sample_data, ESM_depths, MC_depths)
     fowler.sim <- fowler.sim %>% 
            distinct(method, ID, sample_depths, .keep_all = T)  %>% merge(soc_delt_actual) %>%
            mutate(error = soc_change - true_change_soc)
     
     fowler.all <- rbind(fowler.all, fowler.sim)
     
}


fowler.all <- fowler.all %>% merge(fowler.all %>% filter(method == 'Single Depth') 
                                   %>% select(ID, Rep, error) 
                                   %>% rename(FD_error = error),
                                   by = c('ID', 'Rep')) %>% 
  mutate(relative_error = error / abs(FD_error) ) %>% select( -c(FD_error))

fowler.all <- fowler.all %>% distinct(ID, sample_depths, soc_change, .keep_all = TRUE)


fowler.all %>% group_by(method, sample_depths) %>% summarize(RMSE = sqrt(mean(error**2)))



time_series_comparison <- function(FD_data, ESM_depths, MC_depths, lat, lon,
                                   muname, study_name, depth_of_estimate = 30){
  out <- run_comparison(FD_data, ESM_depths = ESM_depths,
                                MC_depths = MC_depths)
  
  
  res <- out$res
  summary <- out$summary
  
  
  
  #summarize the between-treatment comparisons
  #Compare the calculated differences to the benchmark (ESM on all available data)
  calced_difs <- res %>% filter(ID != 'all_t1') %>% group_by(method, sample_depths) %>% 
    reframe(changes = dif_matrix(soc_change))
  
  as.df <- as.data.frame(calced_difs$changes) %>% 
    mutate(method = calced_difs$method, sample_depths = calced_difs$sample_depths)
  calced_difs <- as.df
  
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
  
  return (list('summary' = summary, 'res' = res %>% mutate('study' = study_name)))
}







#franz_FD <- read.csv('Source_data/Franzluebbers_2013.csv')  %>% mutate(Rep = 1)
#muname <- "Cecil fine sandy loam, 2 to 6 percent slopes"  
#ESM_depths <- list(c(20, 40), c(20,40, 60), c(20, 40, 60, 90), c(40, 60), 
#                  c(20, 40 ,60, 90, 120), c(20, 40 ,60, 90, 120, 150))
#MC_depths <- list(c(20, 40), c(40))


#franz_summary <- run_comparison(franz_FD, lat= NULL, lon = NULL,
#                                           muname = muname, ESM_depths = ESM_depths,
#                                           MC_depths = MC_depths, depth_of_estimate = 40)


sainju_FD <- read.csv('Source_data/Sainju_2013.csv') %>% mutate(Rep = 1)
muname <- 'Dooley sandy loam, 0 to 4 percent slopes'
ESM_depths = list(c(15, 30), c(15,30, 60), c(15, 30, 60, 90), c(30, 60), c(7.5, 15, 30, 60),
                  c(7.5, 15, 30, 60, 90, 120), c(7.5, 15,30))
MC_depths = list(c(15, 30), c(30), c(7.5, 15, 30))

sainju_summary <- run_comparison(sainju_FD, ESM_depths = ESM_depths,
                                 MC_depths = MC_depths, study_name = 'Sainju_2013'
                                            )


mishra_fd <- read.csv('Source_data/mishra_2010.csv') %>% mutate(Rep = 1)
study_name = 'mishra_2010'
mishra_ESM <- list(c(10, 20, 30, 40), c(20, 30), c(10,20,30),
                    c(30, 40), c(20,30, 40))
mishra_MC <- list(c(20, 30), c(10, 20, 30), c(10, 30), c(30))

mishra_fd <- mishra_fd %>% mutate(ID = paste(ID, site, sep = '_'))


for (.site in names(mishra_munames)){
  mishra_res <- run_comparison(mishra_fd,
        ESM_depths = mishra_ESM, MC_depths = mishra_MC, 
                                    tudy_name = study_name)}

vanDoren_fd <- read.csv('Source_data/Van_doren_1986.csv') %>% mutate(Rep = 1)

vanDoren_sp_4_t <- run_comparison(vanDoren_fd %>% filter(!grepl('t1', ID) %>% select(-c(Ref_ID))), 
                                           ESM_depths = ESM_depths, 
                                           MC_depths = MC_depths, 
                                  study_name = 'Van_doren_1986.csv')


devine_fd <- read.csv('Source_Data/Devine_2014.csv') %>% mutate(Rep = 1)

devine_res <-run_comparison(devine_fd,ESM_depths = list(c(5, 15,30, 50, 100),
                                                        c(5, 15, 30), 
                                                          c(15, 30), c(15,30, 50),
                                        c(15, 30, 50), c(15, 30, 50, 100), c(30,50)),
                                        MC_depths = list(c(5, 15, 30), c(30), c(15, 30)),
                                        study_name = 'Devine_2014'
                                        )

blanco_fd <- read.csv('Source_Data/Blanco_Canqui(2008).csv') %>% mutate(Rep = 1)

blanco_res <- run_comparison(blanco_fd, 
                             ESM_depths = list(c(5, 10, 30), c(10, 30),  c(30, 50),
                                                   c(10, 30, 50), c(10, 30, 50, 60),
                                                   c(5, 10, 30, 50, 60)
                                                   ),
                      MC_depths = list(c(30), c(10, 30)),
                     study_name = 'Blanco_Canqui(2008)')
  
  




venterea_fd <- read.csv('Source_data/Venterea_2006.csv') %>% mutate(Rep = 1)


venterea_s4_t <- run_comparison(venterea_fd %>% filter(years == 2000) %>% select(-c(Ref_ID)),
                                           ESM_depths = list(c(5, 10, 20, 30, 45, 60),
                                                             c(30, 45), c(20, 30),
                                                             c(20, 30, 45), c(30, 60),
                                                             c(10, 30, 45), c(5, 10, 30),
                                                             c(10, 20, 30), c(5, 10, 20, 30),
                                                             c(5, 10, 20, 30, 45),
                                                             c(5, 10, 20, 30, 45, 60)
                                            
                                                             ),
                                           MC_depths = list(c(30), c(20, 30), c(10, 20, 30), 
                                                            c(5, 10, 20, 30)),
                                           study_name = 'Venterea_2006'
                                           
                                           )

chatterjee_fd <- read.csv('Source_data/chatterjee_2009.csv') %>% mutate(Rep = 1)

 

chj_res <- run_comparison(subset, 
                                ESM_depths = list(c(5, 10, 30), c(10, 30),  c(30, 50),
                                                   c(10, 30, 50), c(10, 30, 50, 60),
                                                   c(5, 10, 30, 50, 60)),
                                   MC_depths = list(c(30), c(10, 30), c(5,10,30)),
                                   study_name = 'chatterjee_2009'
                                   )
    



blanco_fd2 <- read.csv('Source_data/Blanco-canqui_2011.csv') %>% mutate(Rep = 1)



blanc_res_2 <- run_comparison(subset, ESM_depths = list(c(30,40), c(20,30), c(20,30, 40), c(30,60),
                                                   c(20,30,60), c(20, 30,40, 60),
                                                   c(15, 20, 30,40, 60),
                                                   c(15,30), c(10, 30),
                                                   c(20, 30, 40, 60, 80),
                                                   c(15, 20, 30, 40, 60, 80),
                                                   c(2.5, 5, 10, 15, 20, 30, 40, 60, 80, 100)
  ),
  MC_depths = list(c(30), c(10, 30), c(20, 30), c(15,30)),
   study = 'Blanco_Canqui_2011'
  )
  
  


sanford_st_res <- run_comparison(sanford_fd %>% filter(ID != Ref_ID & ID != 'all_t0') %>% 
                                            select(-c(Ref_ID)),
                                            muname = NULL,
                                            ESM_depths = list(c(15, 30), c(15,30, 60), 
                                                        c(15, 30, 60, 90), c(30, 60)),
                                            MC_depths = list(c(15, 30), c(30)),
                                            study_name = 'Sanford_2012'
  
                                                                                      )


poffenbarger_fd <- read.csv('Source_data/poffenbarger_2020.csv') %>% mutate(Rep = 1)




 
poffenbarger_res <- run_comparison(subset, 
                        ESM_depths = list(c(15, 30), 
                        c(30, 60,), c(15,30, 60), c(15,30,60, 90) ),
                         MC_depths = list(c(15, 30), c(30)), 
                        study_name = 'poffenbarger_2020')
                         
                         

out_s4t <- rbind(res_NAEW$res, res_WARS$res, res_NWARS$res, 
      res_DEL$res, res_Cos$res, res_Hoyt$res, blanco_res %>% mutate(Ref = 0, weight = 1), 
      devine_res$res, vanDoren_sp_4_t$res,
      sainju_summary$res, sanford_st_res$res, venterea_s4_t$res, chj_res, blanc_res_2)


summarize_res <- function(grouped_res){
  grouped_res %>% summarise(RMSE = sqrt(weighted.mean(dif_from_benchmark**2, weight)), 
                            RMPSE = sqrt(weighted.mean((dif_from_benchmark /soc_change_benchmark )**2, weight )), 
                            bias = weighted.mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1), weight), 
                            MAPE = weighted.mean(abs(dif_from_benchmark/soc_change_benchmark), weight), ,
                            Max_err =  max(abs(dif_from_benchmark)),
                            mean_abs_err = mean(abs(dif_from_benchmark)))
}


cross_fold_weights_res <- data.frame()
for (study_name in unique(out_s4t$study)){
  training_data <- out_s4t %>% filter(study != study_name)
  sumed_dat <- training_data %>% filter(!grepl("WL", Ref) &  !grepl('WL', ID) & !grepl('Forest', ID) & !grepl('Forest', Ref)) %>%
    filter(sample_depths == 30) %>% group_by(method) %>% 
    summarise(bias = weighted.mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1), weight))
  single_d_bias <- sumed_dat %>% filter(method == 'Single Depth') %>% 
    pull('bias') %>% first()
  linear_MC_bias <- sumed_dat %>% filter(method == 'Mass Correction') %>% 
    pull('bias') %>% first()
   MC_weight <- abs(single_d_bias)/(abs(single_d_bias) + abs(linear_MC_bias))
   sd_weight <- 1-MC_weight
   testing_data <- out_s4t %>% filter(study == study_name) %>% 
     filter(!grepl("WL", Ref) &  !grepl('WL', ID) & !grepl('Forest', ID) & !grepl('Forest', Ref))
   testing_data %>% filter(method %in% c('Single Depth', 'Mass Correction') & sample_depths == 30) %>% 
     pivot_wider(id_cols = c(ID, Ref), names_from = method, values_from = soc_change) %>% 
     mutate(soc_change = `Mass Correction` * MC_weight + `Single Depth` * sd_weight ) %>% 
     merge(testing_data %>% select(ID, Ref, soc_change_benchmark, mass.change, weight) %>% distinct(ID, Ref, .keep_all = T)) %>%
     mutate(dif_from_benchmark = soc_change - soc_change_benchmark, study = study_name,
            method = 'Mass Correction, empirical weights') %>%
     rbind(cross_fold_weights_res ) -> cross_fold_weights_res
   
}




### calculate time-series differences on real data
ESM_depths <- list(c(5, 10, 15,20, 25, 30,35, 40, 45), c(30, 45), c(25, 30, 35), c(15,30, 45), 
                   c(20, 30, 40), c(30, 40),
                   c(10, 20, 30), c(15, 30), c(20, 30))
MC_depths <- list(c(30), c(15,30), c(20, 30), c(25, 30), c(10, 20, 30))
muname <- 'Wooster silt loam, 2 to 6 percent slopes'

vanDoren_res <- time_series_comparison(vanDoren_fd, ESM_depths, MC_depths, lat = NULL,
                                       lon = NULL, muname = muname, study_name = 'van_doren_1986')




sanford_res <- time_series_comparison(sanford_fd, 
                                      ESM_depths = list(c(15, 30), c(15,30, 60), 
                                                   c(15, 30, 60, 90), c(30, 60)),
                                      MC_depths = list(c(15, 30), c(30)), 
                                      lat =  43.29586, lon = -89.38068, 
                                      muname = NULL,
                                      study_name = 'sanford_2012'
)

venterea_res <- time_series_comparison(venterea_fd,  
                                    ESM_depths = list(
                                      c(30, 45), c(20, 30),
                                               c(20, 30, 45), c(30, 60),
                                               c(10, 30, 45), c(5, 10, 30),
                                               c(10, 20, 30), c(5, 10, 20, 30),
                                               c(5, 10, 20, 30, 45),
                                               c(5, 10, 20, 30, 45, 60)),
                          MC_depths = list(c(30), c(20, 30), c(10, 20, 30), 
                 c(5, 10, 20, 30)),
                 lat = NULL, lon = NULL, 
                 muname = 'Waukegan silt loam, 0 to 2 percent slopes',
                 study = 'Venterea_2006.csv')


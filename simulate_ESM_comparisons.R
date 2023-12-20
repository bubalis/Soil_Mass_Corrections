options(warn = -1)
library(stringr)
source('analysis_functions.R')






#Run comparisons for the simulated change data from vanHoden et al
filtered_FD <- read.VanHaden.data() %>%
  preprocess.fd.data()

vanHoden_res <- compare.methods(filtered_FD)

#run the vanHoden data backwards, so that the "after" soil is the
#reference for the "before" soil
reversed <- filtered_FD %>% mutate (Ref_ID = sapply(Ref_ID, mut_Ref_ID))
vanHoden_rev <- compare.methods(reversed )

vanHoden_all <- rbind(vanHoden_res, vanHoden_rev)

vanHoden_all <- vanHoden_all %>%
  merge(vanHoden_all %>%
          filter(method == 'Fixed Depth') %>%
          select(ID, Rep, soc_change) %>%
          rename(FD_error = soc_change)) %>%
  mutate(relative_error = soc_change / abs(FD_error) ) %>%
  select( -c(FD_error)) %>%
  mutate(error = soc_change)

vanHoden_all %>%
  group_by(method, sample_depths) %>%
  summarize(RMSE = sqrt(mean(soc_change**2))) %>%
  print(n = 21)

#compare the results of different methods on the fowler simulated data, and
#variations on it
source("data_sim_fowler.R")
fowler_sim_data <- simulate.soil.Fowler(soc_profile)


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
     ESM_depths <- list(c(10,20,30),  c(15,30), c(24,30), c(24,30, 36), c(20,30),
                        c(30, 36), c(15, 30,50), c(10, 30, 50), c(30, 50),
                        c(10, 30, 50, 70), c(15, 30, 45, 60, 75))
     MC_depths <- list(c(10,20,30),  c(15,30), c(24,30), c(10, 30), c(30))
     split_core_single_SOC_depths <- list(c(20, 30), c(20, 40), c(24,30),
                                          c(30, 36), c(24, 36), c(30, 40))

     fowler.sim     <- compare.methods(sample_data, ESM_depths, MC_depths,
                                       split_core_single_SOC_depths)
     fowler.sim <- fowler.sim %>%
            distinct(method, ID, sample_depths, .keep_all = T)  %>% merge(soc_delt_actual) %>%
            mutate(error = soc_change - true_change_soc)

     fowler.all <- rbind(fowler.all, fowler.sim)

}

fowler.all <- fowler.all %>% merge(fowler.all %>% filter(method == 'Fixed Depth')
                                   %>% select(ID, Rep, error)
                                   %>% rename(FD_error = error),
                                   by = c('ID', 'Rep')) %>%
  mutate(relative_error = error / abs(FD_error) ) %>% select( -c(FD_error))

fowler.all <- fowler.all %>%
  distinct(ID, sample_depths, soc_change, .keep_all = TRUE) %>%
  mutate(
    weight = ifelse(grepl('s2', ID), 10,
                    ifelse(grepl('s3|s4', ID), 1,
                           ifelse(grepl('s5|s6', ID), 2, NA))))

fowler.all %>% group_by(method, sample_depths) %>% summarize(RMSE = sqrt(weighted.mean(error**2, weight)))




studies <- c('Sainju_2013',
             'mishra_2010',
            'Devine_2014',
            'Blanco_Canqui(2008)',
            'chatterjee_2009',
            'Blanco-Canqui_2011',
            'poffenbarger_2020',
            'Yang_1999')

run_test <- function(){
  venterea_fd <- load_study('Venterea_2006')

  #
  run_comparison(venterea_fd %>%
                                    filter(years == 2000) %>%
                                    select(-c(Ref_ID)),
                                  study_name = 'Venterea_2006'
                                 )
}

results <- lapply(studies,
              function(x){run_comparison(load_study(x) %>% select(-one_of('Ref_ID')),
                                                      study_name = x
                                                      )})

#these two studies have both after and before measurements
vanDoren_fd <- load_study('Van_doren_1986')
#the t1 measurements are not good to worrk with. Also, we don't run all combinations for
#vanDoren because there are way too many.
results[[9]] <- run_comparison(vanDoren_fd %>%
                              filter(!grepl('t1', ID)) %>% select(-c(Ref_ID)),
                               ESM_depths = list(c(10, 20, 30), c(10, 30),
                                                 c(20,30), c(15,30),
                                                 c(30,40), c(30, 45),
                                                 c(25,30), c(25,30,35),
                                                 c(20,30,40), c(5,15,30),
                                                 c(5, 15,30,45),
                                                 c(10,20,30,40),
                                                 c(5,10,15,20,25,30,35,40,45),
                                                 c(30,35), c(30)
                                                 ),
                               study_name = 'Van_doren_1986')

venterea_fd <- load_study('Venterea_2006')

#
results[[10]] <- run_comparison(venterea_fd %>%
                                  filter(years == 2000) %>%
                                  select(-c(Ref_ID)),
                                  study_name = 'Venterea_2006'
                               )

results[[11]]<- run_comparison(venterea_fd %>%
                                 filter(years == 2005) %>%
                                 select(-c(Ref_ID)),
                                 study_name = 'Venterea_2006'
                               )



out_s4t <- bind_rows(results)

out_s4t <- out_s4t %>%
  distinct(ID, Ref_ID, study_name, method, sample_depths, .keep_all = T)


out_s4t %>% #select( ID, Ref_ID, Rep, Cum_SOC_g_cm2, method, sample_depths,
  #        Cum_SOC_g_cm2_baseline, soc_change,  study_name) %>%
  #distinct() %>%
  write.csv( file.path('Results', 'all_data_space_for_time_comparison.csv'))



library(sf)
library(parallel)
source('analysis_functions.R')

potash_data <- read.csv('Source_Data/Potash_et_al_data/measurements.csv') %>%
  rename(Upper_cm = sample_depth_min,
         Lower_cm = sample_depth_max,
         SOC_pct = SOCc,
         BD_g_cm3 = BD,
         ID = location_id) %>%
  mutate(SOM_pct = SOC_pct * 1.9,
         Rep = 1,
         Ref_ID = '1') %>%
  mutate(ID = paste(site, ID, sep = '_')) %>%
  calc.cumulative_masses() %>%
  mutate(study_name = site)


locs <- read.csv('Source_Data/Potash_et_al_data/locations.csv')

locs <- st_as_sf(locs, coords = c('X', 'Y'), remove = F) %>%
  mutate(ID = paste(site, location_id, sep = '_'))

d <- st_distance(locs)

closest <- apply(d, 2, function(arr){sort(arr, index.return = T)$ix[2:11]})



fun <- function(i){
  out <- data.frame()
  .ID <- locs[i,] %>% pull(ID)
  lat <- locs[i,] %>% pull(Y)
  lon <- locs[i,] %>% pull(X)
  sub_data <- filter(potash_data, ID == .ID)
  for (i2 in 1:3){
    ref_I <- closest[i2, i]
    .Ref_ID <- locs[ref_I,] %>% pull(ID)
    compare_data <-  filter(potash_data, ID == .Ref_ID)
    if (length(setdiff(compare_data$Lower_cm, sub_data$Lower_cm))==0  ){

      run_data <- bind_rows(sub_data, compare_data) %>%
        mutate(Ref_ID =.Ref_ID)

      res <- tryCatch({compare.methods(run_data, lat = lat, lon = lon)},
                      error = function(cond){data.frame()})
      out <- bind_rows(out, res)}
  }

  return(out)
}

cl <- makeCluster(detectCores() -1)
clusterExport(cl, c('fun', 'locs', 'potash_data', 'closest'), envir = environment())
clusterEvalQ(cl, {library(dplyr)})
clusterEvalQ(cl, {library(tidyr)})
clusterEvalQ(cl, {source('analysis_functions.R')})
#clusterExport(cl, c('check_errors', 'spline_error_within', 'sum_weights_mass',  'sum_weights_mass_2',
#                    'hyman_extrapolating_splinefun', "make_zerodepth_row", "calc.linear_SOC_pct",
#                    'df_agg_to_depths', 'agg_to_depths', 'agg_to_depth',
#                    'calc.cumulative_masses'))

out <- parLapply(cl, 1:nrow(locs), fun) %>% bind_rows()


save_path <- file.path('Results', 'potash_data_comparisons.csv')
out <- mutate(out, study_name = sub('_.*', "", ID))
write.csv(out, save_path)
















#vanDoren_res <- time_series_comparison(vanDoren_fd,
#                                       study_name = 'Van_doren_1986')

#venterea_res <- time_series_comparison(venterea_fd,
#                 study_name = 'Venterea_2006.csv')





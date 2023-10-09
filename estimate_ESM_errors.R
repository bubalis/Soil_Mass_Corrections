library(tidyr)
library(dplyr)
library(ggplot2)

source('ESM_and_mass_correction_functions.R')




spline_error_within <- function(soil_sample, depth_to_quantify = 30){
    masked_val <- soil_sample %>% filter(Lower_cm == depth_to_quantify) 
    
    mass_to_quantify <- masked_val$Cum_Min_Soil_g_cm2
    
    true_SOC <- masked_val$Cum_SOC_g_cm2
    
    soil_sample <- soil_sample %>% filter(Lower_cm != depth_to_quantify)
    
    sp <- splinefun(x = soil_sample$Cum_Min_Soil_g_cm2, 
                    y = soil_sample$Cum_SOC_g_cm2,
                    method = 'hyman')
    est_y <- sp(mass_to_quantify, deriv = 0)
    est_deriv <- sp(mass_to_quantify, deriv = 1)
    
    
    total_weight_mass <- sum_weights_mass(soil_sample, mass_to_quantify)
    #total_weight_depth <- sum_weights_depth(soil_sample, depth_to_quantify )
    error <- est_y -  true_SOC 
    
    return (c(true_SOC, error, est_deriv, total_weight_mass ))
}

samp_data <- read.csv('Source_Data/Potash_et_al_data/measurements.csv') %>% 
  rename(Upper_cm = sample_depth_min, Lower_cm = sample_depth_max,
         SOC_pct = SOCc, BD_g_cm3 = BD, ID = location_id) %>% 
  mutate(SOM_pct = SOC_pct * 1.9, Rep = 1, Ref_ID = '1') %>% calc.cumulative_masses()

other_samps <- bind_rows(sapply(c('Sainju_2013',
                                  'mishra_2010',
                                  'Devine_2014',
                                  'Blanco_Canqui(2008)',
                                  'Venterea_2006',
                                  'chatterjee_2009',
                                  'Blanco-canqui_2011',
                                  'poffenbarger_2020',
                                  'Yang_1999'),
                                load_study ))


van_doren <- load_study('Van_doren_1986') %>%
  filter(!grepl('t1', ID))


other_samps <- bind_rows(other_samps, van_doren)
other_samps$uniq_id <- other_samps$ID

samp_data <- bind_rows(samp_data, other_samps)
samp_data[is.na(samp_data$study), 'study'] <- samp_data[is.na(samp_data$study), 'site']
samp_data <- samp_data %>% mutate( uniq_id = paste(site, ID, sep = '_'))
samp_data$ID <- samp_data$uniq_id




check_errors <- function(.ID, samp_data, depth = 30){
  #print(.ID)
  soil_sample <- filter(samp_data, uniq_id == .ID)
  if ((depth %in% soil_sample$Lower_cm) & (nrow(soil_sample) >3)){
    depths_to_choose_from <- soil_sample %>% filter(Lower_cm !=  depth) %>% pull(Lower_cm)
    
    depth_combos <- lapply(seq(2, length(depths_to_choose_from)), function(x){combn(depths_to_choose_from, x, simplify =  F)}) %>% 
      unlist(recursive = F)
    
    
    .names <- sapply(depth_combos, function(x){paste(x, collapse = ', ')})
    all_res <- lapply(depth_combos, 
                     FUN = function(depths){
                       spline_error_within(soil_sample %>% 
                       filter(Lower_cm %in% c(depths, depth)), depth )} ) %>% 
       data.frame() %>% t()
    colnames(all_res) <- c('True_SOC', 'error', 'estimated_deriv', 'precision_mass')
    rownames(all_res) <- 1:nrow(all_res)
      
      all_res <- all_res %>% data.frame() %>% mutate(depths_tested = .names, depth = depth, ID = .ID, 
             pred_err = 1/sqrt(precision_mass))
      
    return(all_res)
  
  }else{return(data.frame())}
}


post_process <- function(results, depth = 30) {
  
    
    
    return(
      results %>%
      merge(distinct(samp_data, ID,  study), by = 'ID' )
  ) %>% filter(!is.na(error))
  
}

analyze_for_depth <- function(depth, samp_data){
  ids <-   samp_data %>% filter(Lower_cm == depth) %>% pull(uniq_id) %>% unique()
  fun <-  function(x){check_errors(x, samp_data, depth)}
  cl <- makeCluster(detectCores())
  clusterExport(cl, c('fun', 'samp_data', 'depth'), envir = environment())
  clusterEvalQ(cl, {library(dplyr)})
  clusterExport(cl, c('check_errors', 'spline_error_within', 'sum_weights_mass'))
  
  res <-  parLapply(cl, ids, fun) %>% bind_rows()
  
  stopCluster(cl)
  return( res )
}



df <- bind_rows(analyze_for_depth(30, samp_data),
                analyze_for_depth(25, samp_data),
                analyze_for_depth(35, samp_data),
                analyze_for_depth(40, samp_data),
                analyze_for_depth(20, samp_data)
)




df <- df %>% filter(!is.na(error)) %>% merge(distinct(samp_data, ID, study, by = 'ID')) %>% 
  arrange(pred_err) %>% 
  mutate(pred_err_bin1 = round(pred_err), 
         pred_err_bin2 = round(pred_err /2)*2,
         pred_err_bin3 = round(pred_err /5)*5, 
                  pred_err_bin4 = round(pred_err*2)/2,
         even_err_bin = round(row_number() /20 + .5)
         ) 


write.csv(df, 'esm_errors_finished.csv')


####



grouped <- df %>% group_by(even_err_bin) %>% summarize(mean_sq_err = mean(error**2), mean_pred_err = mean(pred_err))

m1 <- lm(error ~ pred_err, data = df)
m2 <- lm(abs(error) ~ 0 + pred_err, data = df)



d <- filter(df, pred_err < 10)
cor1 <- cor(df$pred_err, df$error)
cor2 <- cor(df$pred_err, abs(df$error))
cor3 <- cor(grouped$mean_pred_err**2, grouped$mean_sq_err)
m3 <- lm(mean_sq_err ~ 0 + pred_err_bin**2, data = grouped )

#ggplot(data = df, aes(x = pred_err, y = error, color = 'site')) + geom_point()

#ggplot(data = df, aes(x = pred_err, y = error, color = 'depth')) + geom_point()


ggplot(data = df, aes(x = pred_err, y = abs(error)), color = as.factor(study)) + geom_point()

ggplot(data = df, aes(x = pred_err, y = abs(error), color = depth)) + geom_point()

ggplot(data = grouped %>% filter(pred_err_bin1 <25), aes(x = pred_err_bin1, y = sqrt(mean_sq_err))) + geom_point()




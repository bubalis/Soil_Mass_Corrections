library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)
#library(truncnorm)

source('ESM_and_mass_correction_functions.R')


sigmoid <- function(x){
  1/ (1 + exp(x*-1))
}

logit <- function(p){
  log(p/ (1-p))
}

site_metadata <- read.csv('Source_data/site_metadata.csv')


get_reference_samps <- function(samp_data){
  samp_data %>% filter(ID == 'NA_CC_CT')
}
#soil_sample <- filter(samp_data, ID == 'IL-MC_16') %>% df_agg_to_depths(c(30, 45, 60, 75)) %>% calc.cumulative_masses()
# calc.linear_SOCpct <- function(soil_sample, depth_to_quantify){
#   depth_difs <- soil_sample$Lower_cm - depth_to_quantify
#   above_ind <- which.max(1/(depth_difs *-1))
#   below_ind <- which.max(1/(depth_difs))
#   
#   if (depth_difs[above_ind]>0 ){
#     above_C <- NA
#   }else{above_C <- soil_sample$SOC_pct[above_ind]}
#   if (depth_difs[below_ind]<0 ){
#     below_C <- NA
#   }else{below_C <- soil_sample$SOC_pct[below_ind]}
#   
#   weights <- 1/abs(depth_difs[c(above_ind, below_ind)])
#   
#   
#   
#   return( weighted.mean(c(above_C, below_C), 
#                                 1/abs(depth_difs[c(above_ind, below_ind)]),
#                                 na.rm = T))}



calc.linear_SOCpct_mass <- function(soil_sample, mass_to_quantify){
  depth_difs <- soil_sample$Cum_Min_Soil_g_cm2- mass_to_quantify
  above_ind <- which.max(1/(depth_difs *-1))
  below_ind <- which.max(1/(depth_difs))
  
  if (depth_difs[above_ind]>0 ){
    above_C <- NA
  }else{above_C <- soil_sample$SOC_pct[above_ind]}
  if (depth_difs[below_ind]<0 ){
    below_C <- NA
  }else{below_C <- soil_sample$SOC_pct[below_ind]}
  
  weights <- 1/abs(depth_difs[c(above_ind, below_ind)])
  
  
  
  return( weighted.mean(c(above_C, below_C), 
                        1/abs(depth_difs[c(above_ind, below_ind)]),
                        na.rm = T))}



spline_error_within <- function(soil_sample, depth_to_quantify = 30){
  depths <- soil_sample$Lower_cm[soil_sample$Lower_cm != depth_to_quantify]  
  masked_val <- soil_sample %>% 
      filter(Lower_cm == depth_to_quantify) 
    
    mass_to_quantify <- masked_val$Cum_Min_Soil_g_cm2
    true_SOC <- masked_val$Cum_SOC_g_cm2
    
    soil_sample <- soil_sample %>% 
      df_agg_to_depths(depths) %>% 
      calc.cumulative_masses() %>% 
      ungroup() %>% 
      bind_rows(make_zerodepth_row(.)) %>%
      arrange(Lower_cm)
    
    sp <- hyman_extrapolating_splinefun(
                    xin = soil_sample$Cum_Min_Soil_g_cm2, 
                    yin =  soil_sample$Cum_SOC_g_cm2
                    )
    
    est_y <- sp(mass_to_quantify, deriv = 0)
    
    #if (mass_to_quantify > max(soil_sample$Cum_Min_Soil_g_cm2)){
    #  est_deriv <- sp(max(soil_sample$Cum_Min_Soil_g_cm2), 1)
    #}
    #else{
    est_deriv <- sp(mass_to_quantify, deriv = 1)#}
    
    total_weight_mass <- sum_weights_mass(soil_sample, mass_to_quantify)
    weight_mass_alt <- sum_weights_mass_2(soil_sample, mass_to_quantify)
    total_weight_SOC <- sum_weights_mass(mutate(soil_sample, Cum_Min_Soil_g_cm2 = Cum_SOC_g_cm2,
                                                ), est_y)
    
    linear_deriv <- calc.linear_SOC_pct(soil_sample, mass_to_quantify)
    
    error <- est_y -  true_SOC 
    
    return (c(true_SOC, error, est_deriv, total_weight_mass, total_weight_SOC, linear_deriv,
              weight_mass_alt))
}



check_errors <- function(.ID, samp_data, depth = 30, debug = F){
  #print(.ID)
  soil_sample <- filter(samp_data, uniq_id == .ID)
  if ((depth %in% soil_sample$Lower_cm) & (nrow(soil_sample) >3)){
    depths_to_choose_from <- soil_sample %>% 
      filter(Lower_cm !=  depth) %>% 
      pull(Lower_cm)
    
    depth_combos <- lapply(seq(2, length(depths_to_choose_from)), 
                    function(x){combn(depths_to_choose_from, x, simplify =  F)}) %>% 
      unlist(recursive = F)
    
   
    if (debug){
      fun <- function(depths){ print(depths)
        dat_to_run <- df_agg_to_depths(soil_sample, c(depth, depths)) %>% 
          calc.cumulative_masses()
        spline_error_within(dat_to_run,
                              depth)}
      
    }else{
      fun <-  function(depths){
        dat_to_run <- df_agg_to_depths(soil_sample, c(depth, depths)) %>% 
          calc.cumulative_masses()
        spline_error_within(dat_to_run, 
                            depth )}
    }
    
    .names <- sapply(depth_combos, function(x){paste(x, collapse = ', ')})
    all_res <- lapply(depth_combos, 
                     FUN = fun ) %>% 
       data.frame() %>% t()
    colnames(all_res) <- c('True_SOC', 'error', 'estimated_deriv', 
                           'precision_mass', 'precision_SOC', 'linear_deriv', 'precision_mass_alt')
    
    rownames(all_res) <- 1:nrow(all_res)
      
      all_res <- all_res %>% data.frame() %>% 
        mutate(depths_tested = .names, depth = depth, ID = .ID, 
             pred_err = 1/sqrt(precision_mass))
      
    return(all_res)
  
  }else{return(data.frame())}
}




analyze_for_depth <- function(depth, samp_data, debug = F){
  print(depth)
  ids <-   samp_data %>% filter(Lower_cm == depth) %>% pull(uniq_id) %>% unique()
  fun <-  function(x){check_errors(x, samp_data, depth, debug = debug)}
  if (!debug){
  cl <- makeCluster(detectCores())
  clusterExport(cl, c('fun', 'samp_data', 'depth'), envir = environment())
  clusterEvalQ(cl, {library(dplyr)})
  clusterEvalQ(cl, {library(tidyr)})
  clusterExport(cl, c('check_errors', 'spline_error_within', 'sum_weights_mass',  'sum_weights_mass_2',
                      'hyman_extrapolating_splinefun', "make_zerodepth_row", "calc.linear_SOC_pct",
                      'df_agg_to_depths', 'agg_to_depths', 'agg_to_depth',
                      'calc.cumulative_masses'))
  
  res <-  parLapply(cl, ids, fun) %>% bind_rows()
  
  stopCluster(cl)
  }else{
    res <- list()
    for (i in ids){
      print(i)
      out <- fun(i)
      res[[i]] <- out
    }
   res <- bind_rows(res)
  }
  
  
  return( res )
}

samp_data <- read.csv('Source_Data/Potash_et_al_data/measurements.csv') %>% 
  rename(Upper_cm = sample_depth_min, 
         Lower_cm = sample_depth_max,
         SOC_pct = SOCc, 
         BD_g_cm3 = BD, 
         ID = location_id) %>% 
  mutate(SOM_pct = SOC_pct * 1.9, 
         Rep = 1, 
         Ref_ID = '1') %>% 
  calc.cumulative_masses()

other_samps <- bind_rows(sapply(c('Sainju_2013',
                                  'mishra_2010',
                                  'Devine_2014',
                                  'Blanco_Canqui(2008)',
                                  'Venterea_2006',
                                  'chatterjee_2009',
                                  'Blanco-canqui_2011',
                                  'poffenbarger_2020',
                                  'Yang_1999',
                                  'Van_doren_1986'),
                                load_study ))



samp_data <- bind_rows(samp_data, other_samps)

samp_data[is.na(samp_data$study), 'study'] <- samp_data[is.na(samp_data$study), 'site']
samp_data <- samp_data %>% mutate( uniq_id = paste(site, ID, sep = '_'), 
                                   site = replace_na(site, ''))
samp_data$ID <- samp_data$uniq_id


df <- bind_rows(analyze_for_depth(30, samp_data),
                analyze_for_depth(25, samp_data),
                analyze_for_depth(35, samp_data),
                analyze_for_depth(40, samp_data),
                analyze_for_depth(20, samp_data),
                analyze_for_depth(22.5, samp_data))


### To do: 
### Estimate squared error relative to error statistic
### To estimate errors: 
### 1: Calculate the Error statistic
### Normalize the interpolation interval to 0-1
### Normalize the variance to the same interval
### Convert variance to a re-paramterized beta distribution:
### alpha = mean*variance
### beta = variance * (1-mean)
### Use that to create errors/confidence intervals.

gen_beta_params<- function(.mean, .var){
  return(c(.mean * 1/.var, (1-.mean)/.var))
}

df <- df %>% filter(!is.na(error)) %>% 
  merge(distinct(samp_data, site, ID, study, by = 'ID')) %>% 
  arrange(pred_err) %>% 
  mutate(even_err_bin = round(row_number() /20 + .5),
         rel_err = error/ True_SOC,
         pred_err_SOC = pred_err * linear_deriv) 

df <- df %>%
  mutate(pred_err_alt = sqrt(1/precision_mass_alt) ) %>%
  
  mutate(pred_err_SOC_alt =  pred_err_alt *linear_deriv) %>%
  mutate(
         pred_err_SOC = pred_err * linear_deriv
         ) %>% 
  mutate(error_relative = error / linear_deriv) %>% 
  mutate(log_err_rel = log(abs(error)) /linear_deriv )


results_path <- file.path('Results', 'esm_errors_empirical.csv')
write.csv(df, results_path)



grouped <- df %>% 
  group_by(even_err_bin) %>% 
  summarize(mean_sq_err = mean(error**2), mean_pred_err = mean(pred_err))



df$err_pred_adj <- df$pred_err * df$estimated_deriv




m.1 <- lm(abs(error) ~ 0 + pred_err:study, data = df  #%>% filter(pred_err < 15) 
          )
m.2 <- lm(abs(error) ~ 0 + pred_err_SOC:study, data = df %>% filter(pred_err < 15) )
m.3 <- lm(log(abs(error)) ~ log(pred_err):study, data = df)
m.4 <- lm(log(abs(error)) ~   log(pred_err_alt):study, data = df )
m.5 <- lm(log(abs(error)) ~   0 +log(pred_err_SOC):study + study, data = df )
m.9 <- lm(log(abs(error_relative)) ~  pred_err:study + study  , data = df %>% filter(pred_err <15))
m.6 <- lm(log(abs(error)) ~   log(pred_err_SOC_alt):study, data = df )
m.7 <- lm(abs(error)  ~ 0 + pred_err_alt:study, data =df)
m.8 <- lm(abs(error) ~ 0 +pred_err_SOC_alt:study, data =df)
m.10 <- lm(abs(error_relative) ~ 0+pred_err:study +study, data = df )
m.12 <- lm(log(abs(error_relative))~  log(pred_err), data = df)
m.13 <- lm(log(abs(error_relative))~ 0 + study + log(pred_err), data = df %>% filter(pred_err <10))

coefs <- coef(m.2)
names(coefs) <- sapply(names(coefs), function(x) str_replace(x, 'pred_err_SOC:study', ''))

m1 <- lm(error ~ pred_err, data = df)
m2 <- lm(abs(error) ~ 0 +pred_err, data = df %>% filter(pred_err < 10 & estimated_deriv >0))

    
m4 <- lm(abs(error)~ 0 + pred_err_SOC, data = df %>% filter(pred_err < 10 & df$estimated_deriv >0))

m5 <- lm(log(abs(error)) ~ log(pred_err), data = df)
m6 <- lm(log(abs(error) ) ~ log(pred_err_SOC), data = df %>% filter(pred_err < 10 & df$estimated_deriv >0 & !grepl('WL|Forest', ID)))
m7 <- lm(log(abs(error))~pred_err_SOC,  data = df  )
m8 <- lm(sqrt(abs(error) )~0 + pred_err, data = df  )
m9 <- lm(sqrt(abs(error))~ 0 + pred_err_SOC, data = df %>% filter(pred_err <10))

m10 <- lm(log(abs(error)) ~   log(pred_err_alt), data = df )
m11 <- lm(log(abs(error)) ~   log(pred_err_SOC), data = df )
m12 <- lm(sqrt(abs(error))~ 0 + pred_err_SOC_alt, data = df )

d <- filter(df, pred_err < 10)
cor1 <- cor(df$pred_err, df$error)
cor2 <- cor(df$pred_err, abs(df$error))
cor3 <- cor(grouped$mean_pred_err**2, grouped$mean_sq_err)
cor4 <- cor(df$err_pred_adj, abs(df$error))

m3 <- lm(mean_sq_err ~ 0 + pred_err_bin**2, data = grouped )

ggplot(data = df %>% filter(pred_err <10), aes(x = err_pred_adj, y = abs(error), color = 'site')) + geom_point()

#ggplot(data = df, aes(x = pred_err, y = error, color = 'depth')) + geom_point()


ggplot(data = df, aes(x = pred_err, y = abs(error_relative), 
       color = as.factor(depth)
       )
       ) +geom_point(alpha = .3) +
  scale_x_continuous(name = 'Predicted Error Statistic (unitless)', 
                     trans = 'log'
                     ) +
  scale_y_continuous(name = 'Absolute Error, Actual (g C/ cm2', trans = 'log'
                     ) + facet_wrap(~as.factor(study))

ggplot(data = df, aes(x = pred_err, y = abs(error), color = depth)) +  geom_point(alpha = .3) +
  scale_x_continuous(name = 'Predicted Error Statistic (unitless)', trans = 'log',
                     breaks =c(3, 8, 20)) +
  scale_y_continuous(name = 'Absolute Error, Actual (g C/ cm2)', trans = 'log',
                     breaks = c(.00001, .0001, .001, .01, .1, 1))

ggplot(data = df, aes(x = pred_err, y = abs(rel_err)* 100, color = depth)) + geom_point(alpha = .4) +
  scale_x_continuous(name = 'Predicted Error Statistic (unitless)', trans = 'log',
                     breaks =c(3, 8, 20)) +
  scale_y_continuous(name = '% Error', trans = 'log',
                     breaks = c(.0001, .001, .01, 1, 10, 100))




ggplot(data = grouped %>% filter(pred_err_bin1 <25), aes(x = pred_err_bin1, y = sqrt(mean_sq_err))) + geom_point()



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
  calc_portion <- function(x){exp(log(x)*.coef + intercept  )}
  if (use_SOC) {
    return(function(soil_sample, x){
      primary_func(soil_sample, x) * calc.linear_SOC_pct(soil_sample, x)
     
      
    })
  }else{
  return (primary_func)}
}
err_conv_func <- function(soil_sample, x){
  .coef <- .9795 
  intercept <- -6.2387
  err_stat<- exp(log((1/sqrt(sum_weights_mass(soil_sample, x))))*.coef + intercept)
  return(err_stat * calc.linear_SOC_pct(soil_sample, x))
}


err_conv_func_log <- function(x){
  exp(log(x)* .2417 - 4.4394079)
}

err_conv_func_sq <- function(x){
  (x * .0133)**2
}

make_err_conv_func <- function(model, study){
  .coefs <- coef(model)
  .coef <- .coefs[sapply(names(.coefs),function(x) grepl(study, x))]
  return(function(x)x*.coef)
}





#samp_data$site <- replace_na(samp_data$site, "") 
#samp_data <- samp_data %>% left_join(site_metadata %>% select(site, study, N.samples)) %>%
#  mutate(N.samples = replace_na(N.samples, 1))



make_counterfactual_df <- function(soil_sample, depth_to_quantify = 30,
                                   err_conv_func){
  
  masked_val <- soil_sample %>% 
    filter(Lower_cm == depth_to_quantify) 
  
  depths <- soil_sample$Lower_cm[soil_sample$Lower_cm != depth_to_quantify]
  
  
  soil_sample <- soil_sample %>% 
    df_agg_to_depths(depths) %>% 
    calc.cumulative_masses() %>% 
    ungroup() %>% 
    bind_rows(make_zerodepth_row(.)) %>%
    arrange(Lower_cm)
  
  
  mass_to_quantify <- masked_val$Cum_Min_Soil_g_cm2
  true_SOC <- masked_val$Cum_SOC_g_cm2
  
  
  xin <- soil_sample$Cum_Min_Soil_g_cm2
  
  sp <- hyman_extrapolating_splinefun(
    xin = xin, 
    yin =  soil_sample$Cum_SOC_g_cm2
    
  )
  
  
  xout <- c(seq(1, ceiling(max(soil_sample$Cum_Min_Soil_g_cm2))+ 5, 1 ), soil_sample$Cum_Min_Soil_g_cm2, 
            masked_val$Cum_Min_Soil_g_cm2)
  
  xout <- xout[order(xout)]
  
  y_est <- sp(xout, 0)
  
  
  
  
  se <- sapply(xout, function(x){err_conv_func(soil_sample, x)})
  se <- replace_na(se, 0)                                    
  vars_of_est <- se**2 
  
  
  
 
  # 
   masses <- soil_sample$Cum_Min_Soil_g_cm2
   above <- sapply(xout, function(x){max(masses[masses < x]) } )
   above[1] <- 0
  # 
   above_SOC_func <- function(i){
     out <- (soil_sample %>% filter(Cum_Min_Soil_g_cm2 == above[i]) %>% pull(Cum_SOC_g_cm2) %>% first() )
     if (length(out) == 0) {out <- 0}
     return(out)
   }
  # 
   above_SOC <- sapply(seq(length(xout)), above_SOC_func)
  # 
   interval_func<- function(i){
     out <- soil_sample %>% mutate(interval = lead(Cum_SOC_g_cm2) - (Cum_SOC_g_cm2)  ) %>%
       filter(Cum_Min_Soil_g_cm2 == above[i]) %>% pull(interval) %>% max()
     if (length(out) == 0) {out <- NA}
     return(out)
     } 
  # 
  interval <- sapply(seq(length(xout)), interval_func) %>% unlist()
  # 
  # 
  # 
  distance <- y_est - above_SOC
  distance[distance == 0 & !(xout %in% xin)] <- y_est[distance == 0 & !(xout %in% xin)] - .000000000001
  # 
  # 
  # 
  params <- estBetaParams(distance/interval, vars_of_est/(interval**2))
  # 
  lo <- qbeta(.025, params$alpha, params$beta) * interval + above_SOC
  hi <- qbeta(.975, params$alpha, params$beta) * interval + above_SOC
  lo[xout %in% xin] = y_est[xout %in% xin]
  hi[xout %in% xin ] = y_est[xout %in% xin]
  lo[is.na(lo)] = above_SOC[is.na(lo)]
  lo[replace_na(lo > y_est, F)] <- y_est[replace_na(lo > y_est, F)]
  hi[replace_na (hi < y_est, F) ] <- y_est[replace_na (hi < y_est, F) ]
  # 
  # 
  
  
  
  return(data.frame(Cum_Min_Soil_g_cm2 = xout, 
                    Cum_SOC_g_cm2 = y_est,
                    low = lo,
                    hi = hi,
                    ID = first(soil_sample$ID),
                    se = se
         ))
  
}




van_doren_samps <- samp_data %>% 
  filter(grepl('doren', study) & !grepl('t1', ID)) %>% 
  calc.cumulative_masses()

devine_samps <- samp_data %>% 
  filter(grepl('Devine', study)) %>%
  calc.cumulative_masses()

mish_samps <- samp_data %>% 
  filter(grepl('mishra', study)) %>% 
  calc.cumulative_masses()

vent_samps <- samp_data %>% filter(grepl('terea', study)) %>%
  calc.cumulative_masses(
    
  )

#study_data <- vent_samps 
for (study_data  in list(van_doren_samps, devine_samps, mish_samps, vent_samps)){
  
  
out_of_conf_int <- c()
#difs_se <- c()
depths <- unique(study_data$Lower_cm)
ho_depths <- depths[10<depths & depths <55 & depths!= max(depths) ]
name <- first(study_data$study)
error_function <- err_conv_func #make_err_conv_func_log(coef(m.12), first(study_data$study), T)

out_dir <- file.path('Figures', 'cross_validations', name)

if (!dir.exists(out_dir)){dir.create(out_dir)}

for (ho_depth in ho_depths){

  study_data <-  study_data %>% mutate(held_out = Lower_cm == ho_depth)
  

sim_d <- study_data %>% group_by(ID) %>% 
  group_map(function(x, ...) make_counterfactual_df(x, ho_depth, error_function), .keep = T)

sim_d <- do.call(rbind, sim_d)
merg <- study_data %>% filter(held_out) %>% 
  merge(sim_d, by = c('ID', 'Cum_Min_Soil_g_cm2'), suffixes = c('true', 'est'))
out_of_conf_int <- append(out_of_conf_int , merg %>% summarize(out_of_int = sum(Cum_SOC_g_cm2true <low | Cum_SOC_g_cm2true >hi )))
#difs_se <- append(difs_se, (merg$Cum_SOC_g_cm2true - merg$Cum_SOC_g_cm2est)/merg$se)

g <- ggplot(data = sim_d, aes(x = -Cum_Min_Soil_g_cm2 * 100, y = Cum_SOC_g_cm2* 100)) + 
  geom_smooth(stat = 'identity', color = 'black', size = .25) + 
  geom_ribbon(aes(ymin = low*100, ymax = hi*100), fill ='blue', alpha = .3)  + 
  
  geom_point(data = filter(study_data, !held_out), 
             aes(x = -Cum_Min_Soil_g_cm2*100, y = Cum_SOC_g_cm2*100), color = 'blue') + 
  
  geom_point(data = filter(study_data, held_out), 
             aes(x = -Cum_Min_Soil_g_cm2*100, y = Cum_SOC_g_cm2*100), color = 'red')  +
  xlab('Cumulative Mineral Soil, T/ha') +
  scale_y_continuous(name ='Cumulative SOC, T/ha', breaks = ) +
  coord_flip() + facet_wrap(~ID)
  


for (.id in unique(sim_d$ID)){
  max_mass <- (ceiling(max(sim_d$Cum_Min_Soil_g_cm2) * 1.01 /20) * 20) * 100 
  sub_dat <- filter(sim_d, ID == .id)
plot <- ggplot(data = sub_dat, aes(x = -Cum_Min_Soil_g_cm2 *100, y = Cum_SOC_g_cm2 * 100)) + 
  geom_smooth(stat = 'identity', color = 'black') + 
  
  geom_ribbon(aes(ymin = low*100, ymax = hi*100), fill ='blue', alpha = .3) + 
  
  geom_point(data = filter(study_data, !held_out & ID == .id), 
        aes(x = -Cum_Min_Soil_g_cm2 * 100, y = Cum_SOC_g_cm2 *100), color = 'black') + 
  
  geom_point(data = filter(study_data, held_out & ID == .id), 
             aes(x = -Cum_Min_Soil_g_cm2 * 100, y = Cum_SOC_g_cm2 * 100),
             color = 'red')  +
  
  ylab('Cumulative SOC, T/ha') + 
  scale_x_continuous(name = 'Cumulative Mineral Soil, T/ha', 
                     breaks = seq(max_mass * -1, 0, 2000 ),
                     labels = seq(max_mass *-1, 0, 2000 ) * -1) +
  coord_flip()+ 
  ggtitle(sprintf('%s with %s cm depth heldout', .id, ho_depth))
print(plot)

ggsave(file.path(out_dir, sprintf('%s_%s.png', .id, ho_depth)))
}
}
print(out_of_conf_int)
#print(sqrt(mean(difs_se**2)))
}



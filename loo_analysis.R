library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(metafor)
library(cowplot)
#library(truncnorm)

source('analysis_functions.R')

site_metadata <- read.csv('Source_data/site_metadata.csv')

locs <- read.csv('Source_Data/Potash_et_al_data/locations.csv')%>%
  mutate(ID = paste(site, location_id, sep ='_')) %>%
  rename(lon = X, lat = Y) %>% mutate(muname = NA)



get_soil_metadata <- function(.ID,.site, .study){
  if (.site == .study){
    data <- filter(locs, ID == .ID & site == .site)
    return(list(muname = NULL, lat = data$lat, lon = data$lon))
  }else{
    data <- filter(site_metadata, study == .study)

  if (.site != ''){
    data <- filter(data, site == .site)
  }
  return(list(muname = data$muname, lat = data$lat, lon = data$lon))
  }
  }



get_reference_samps <- function(samp_data){
  samp_data %>% filter(ID == 'NA_CC_CT')
}





#'Leave-one-out interpolation of Cum_Mineral_Soil, Cum_SOC spline
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
    #weight_mass_alt <- sum_weights_mass_2(soil_sample, mass_to_quantify)
    total_weight_SOC <- sum_weights_mass(mutate(soil_sample, Cum_Min_Soil_g_cm2 = Cum_SOC_g_cm2,
                                                ), est_y)

    linear_deriv <- calc.linear_SOC_pct(soil_sample, mass_to_quantify)

    error <- est_y -  true_SOC

    return (c(true_SOC, error, est_deriv, total_weight_mass, total_weight_SOC, linear_deriv
              #weight_mass_alt
              ))
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
                           'precision_mass', 'precision_SOC', 'linear_deriv')

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
  clusterEvalQ(cl, {library(devtools)})
  clusterEvalQ(cl, {load_all('soil.mass.corrections')})
  clusterExport(cl, c('check_errors', 'spline_error_within', 'sum_weights_mass',
                      #'sum_weights_mass_2',
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

check_errors_linear <- function(.ID, samp_data, depth){

  soil_sample <- filter(samp_data, uniq_id == .ID)
  if ((depth %in% soil_sample$Lower_cm) & (nrow(soil_sample) >3)){
    depths <- soil_sample$Lower_cm

    samp_held_out <- df_agg_to_depths(soil_sample, depths[depths != depth]) %>%
      calc.cumulative_masses(override = T)
    mass_to_quantify <- filter(soil_sample, Lower_cm == depth) %>%
      pull(Cum_Min_Soil_g_cm2)
    true_SOC <- filter(soil_sample, Lower_cm == depth) %>%
      pull(Cum_SOC_g_cm2)

    soc_range <- filter(samp_held_out, Lower_cm > depth) %>%
      pull(Cum_SOC_g_cm2)  %>% min() -  filter(samp_held_out, Lower_cm < depth) %>%
      pull(Cum_SOC_g_cm2) %>%  c(0) %>% max()


  estimate <- linear_interpolation_CUM_soc(samp_held_out, mass_to_quantify)

  return(data.frame(error = estimate - true_SOC, ID = .ID,
                    depth = depth, soc_range = soc_range))
  }else{
    return(data.frame())
  }

}

check_errors_ssurgo <- function(.ID, samp_data, depth){
  soil_sample <- filter(samp_data, uniq_id == .ID)
  if ((depth %in% soil_sample$Lower_cm) & (nrow(soil_sample) >3)){

      soil_meta <- get_soil_metadata(.ID,  first(soil_sample$site), first(soil_sample$study))
      fitted.spline <- get_SSurgo_spline_ESM(lat = soil_meta$lat,
                                             lon = soil_meta$lon,
                                             muname = soil_meta$muname)


      true_depth <- filter(soil_sample, Lower_cm == depth)

      soil_sample$depth_dif <- abs(soil_sample$Lower_cm - depth)
      #test_samp <- filter(soil_sample$Lower_cm > depth)
      test_samp <- filter(soil_sample, Lower_cm != depth & Lower_cm > 10) %>%
        arrange(depth_dif) %>%
        first()

      adj_factor <- Ssurgo_spline_MC(true_depth$Cum_Min_Soil_g_cm2,
                                      test_samp$Cum_Min_Soil_g_cm2,
                                      test_samp$Cum_SOC_g_cm2,
                                      fitted.spline)

      res_spline <- mass_correction_single_depth(true_depth$Cum_Min_Soil_g_cm2,
                                   test_samp$Cum_Min_Soil_g_cm2,
                                   test_samp$Cum_SOC_g_cm2,
                                   adj_factor)

      exp_decay_factors <- get_SSurgo_exp_decay_params(soil_meta$lat, soil_meta$lon, soil_meta$muname)

      decay_res <- decay_ESM(true_depth$Cum_Min_Soil_g_cm2,
                test_samp$Cum_Min_Soil_g_cm2,
                test_samp$Cum_SOC_g_cm2,
                a = exp_decay_factors$a, b = exp_decay_factors$b,
                depth_of_quant = depth)

      estimate <- mean(res_spline, decay_res)

      return(data.frame(error = estimate - true_depth$Cum_SOC_g_cm2, ID = .ID,
                        depth = depth, ssurgo_depth = test_samp$Lower_cm))

  }else{
    return(data.frame())
  }


}


analyze_linear_depth <- function(depth, samp_data, debug = F){
  ids <-samp_data %>% filter(Lower_cm == depth) %>% pull(uniq_id) %>% unique()
  fun <-  function(x){check_errors_linear(x, samp_data, depth)}
  if (!debug){
    cl <- makeCluster(detectCores())
    clusterExport(cl, c('fun', 'samp_data', 'depth',
                        'check_errors_linear'), envir = environment())
    clusterEvalQ(cl, {library(dplyr)})
    clusterEvalQ(cl, {library(tidyr)})
    clusterEvalQ(cl, {library(devtools)})
    clusterEvalQ(cl, {load_all('soil.mass.corrections')})
    clusterExport(cl, 'check_errors')


    res <-  parLapply(cl, ids, fun) %>% bind_rows()
  }


  else{
    res <- list()
    for (i in ids){
      print(i)
      out <- fun(i)
      res[[i]] <- out
    }
    res <- bind_rows(res)
  }

  return(res)
}


analyze_depth_ssurgo <- function(depth, samp_data, debug = F){
  ids <-samp_data %>% filter(Lower_cm == depth) %>% pull(uniq_id) %>% unique()
  fun <-  function(x){check_errors_ssurgo(x, samp_data, depth)}
  if (!debug){
    cl <- makeCluster(detectCores())
    clusterExport(cl, c('fun', 'samp_data', 'depth', 'locs',
                        'check_errors_ssurgo', 'site_metadata',  'get_soil_metadata'),
                  envir = environment())
    clusterEvalQ(cl, {library(dplyr)})
    clusterEvalQ(cl, {library(tidyr)})
    clusterEvalQ(cl, {library(devtools)})
    clusterEvalQ(cl, {load_all('soil.mass.corrections')})


    res <-  parLapply(cl, ids, fun) %>% bind_rows()
  }


  else{
    res <- list()
    for (i in ids){
      print(i)
      out <- fun(i)
      res[[i]] <- out
    }
    res <- bind_rows(res)
  }

  return(res)

}


samp_data <- read.csv('Source_Data/Potash_et_al_data/measurements.csv') %>%
  dplyr::rename(Upper_cm = sample_depth_min,
         Lower_cm = sample_depth_max,
         SOC_pct = SOCc,
         BD_g_cm3 = BD,
         ID = location_id) %>%
  mutate(SOM_pct = SOC_pct * 1.9,
         Rep = 1,
         Ref_ID = '1')

other_samps <- bind_rows(sapply(c('Sainju_2013',
                                  'mishra_2010',
                                  'Devine_2014',
                                  'Blanco_Canqui(2008)',
                                  'Venterea_2006',
                                  'chatterjee_2009',
                                  'Blanco-Canqui_2011',
                                  'poffenbarger_2020',
                                  'Yang_1999',
                                  'Van_doren_1986',
                                  'Franzluebbers_2013'),
                                load_study ))




samp_data <- bind_rows(samp_data, other_samps)

vh_data <- read.VanHaden.data() %>%
  filter(!grepl('Y5', ID)) %>% # in the Van-Haden Data the Y5 data is simulated.
  mutate(study = 'Von_Haden_2020')

samp_data <- bind_rows(samp_data, vh_data)


samp_data[is.na(samp_data$study), 'study'] <- samp_data[is.na(samp_data$study), 'site']
samp_data <- samp_data %>% mutate( uniq_id = paste(site, ID, sep = '_'),
                                   site = replace_na(site, ''))
samp_data$ID <- samp_data$uniq_id

samp_data <- calc.cumulative_masses(samp_data, override = T)


df <- lapply(c(15, 20, 22.5, 25,30, 35, 40,45),
             function(x) analyze_for_depth(x, samp_data)) %>%
  bind_rows() %>%
  rename(linear_error = error)





linear_res <- lapply(c(15, 20, 22.5, 25,30, 35, 40,45),
                     function(x) analyze_linear_depth(x, samp_data)) %>%
  bind_rows() %>%
  rename(linear_error = error)


ssurgo_res <- lapply(c(15, 20, 22.5, 25,30, 35, 40,45),
                     function(x) analyze_depth_ssurgo(x, samp_data %>% filter(study != "Franzluebbers_2013" & study != 'Von_Haden_2020'), T)) %>%
  bind_rows() %>%
  rename(ssurgo_error = error)






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
  mutate(
         pred_err_SOC = pred_err * linear_deriv
         ) %>%
  mutate(error_relative = error / linear_deriv) %>%
  mutate(log_err_rel = log(abs(error)) /linear_deriv )


results_path <- file.path('Results', 'esm_errors_empirical.csv')
linear_res_path <- file.path('Results', 'linear_errors_empirical.csv')
ssurgo_res_path <- file.path('Results', 'ssurgo_errors_loo.csv')

write.csv(linear_res, linear_res_path)

write.csv(ssurgo_res, ssurgo_res_path)

write.csv(df, results_path)
df$max_depth_tested <- str_split(df$depths_tested, ',|, ') %>%
  lapply(function(x) max(as.numeric(x))) %>% unlist()




loo_data <-  df %>%
  group_by(ID, depth) %>%
  group_modify(~filter(.x, nchar(depths_tested) == max(nchar(.x$depths_tested)) )) %>%
  filter(max_depth_tested > depth) %>%
  ungroup %>%
  distinct(error, True_SOC, pred_err, .keep_all = T) %>%
  merge(linear_res) %>%
  mutate(error_difference = abs(linear_error) - abs(error),
         error_frac = error / soc_range) %>%
  merge(ssurgo_res)

loo_data$site_type <- 'Other'
loo_data[grepl('NT', loo_data$ID) | grepl('No Till', loo_data$ID),
         'site_type'] <- 'NT'
loo_data[grepl('WL', loo_data$ID) | grepl('Forest', loo_data$ID),
         'site_type'] <- 'Forest'
loo_data[loo_data$site %in% c('IL-BR', 'IL-RD', 'NE'), 'site_type'] <- 'NT'



loo_data %>% filter((grepl('NT', ID) | grepl('No Till', ID)) &
                !(grepl('t1', ID) | grepl('2000', ID)) ) %>%
  group_by(depth) %>%
  summarize(bias = mean(error),
            RMSE = sqrt(mean(error**2)),
            Mean.Abs.Err = mean(abs(error)),
            count = n())


loo_data %>% filter((grepl('NT', ID) | grepl('No Till', ID)) &
                !(grepl('t1', ID) | grepl('2000', ID)) ) %>%
  filter(depth %in% c(20,30,40)) %>%
  ggplot(aes(x = error, fill = factor(depth))) +
  geom_histogram() +
  facet_wrap(~depth) +
  geom_vline(xintercept = 0, linetype = 'dashed')

loo_data %>% filter(grepl('WL', ID) | grepl('Forest', ID)) %>%
  group_by(depth) %>%
  summarize(bias = mean(error),
            RMSE = sqrt(mean(error**2)),
            Mean.Abs.Err = mean(abs(error)),
            count = n())


mode.stat <- function(x) {
  t <- table(x)
  names(t)[ which.max(t) ]
}


summary1 <- loo_data %>%
  filter(max_depth_tested > depth) %>%
  filter() %>%
  group_by(depth, study) %>%
  summarize(bias = mean(error),
            RMSE = sqrt(mean(error**2)),
            Mean.Abs.Err = mean(abs(error)),
            count = n(),
            perc_greater = sum(error > 0) / n(),
            variance = var(error),
            err_stat = mean(pred_err),
            depths_used = mode.stat(depths_tested),
            bias_linear = mean(linear_error),
            variance_linear = var(linear_error),
            bias_ssurgo = mean(ssurgo_error),
            variance_ssurgo = var(ssurgo_error)
            ) %>%
  mutate(var_est = variance / count,
         var_est_linear = variance_linear/count,
         var_est_ssurgo = variance_ssurgo/count)

.sum <- loo_data %>%
  group_by(ID, depth) %>%
  group_modify(~filter(.x, nchar(depths_tested) == max(nchar(.x$depths_tested)) ))%>%
  distinct(error, True_SOC, pred_err, .keep_all = T) %>%
  filter(max_depth_tested > depth) %>%
  group_by(depth) %>%
  summarize(bias = mean(error),
            RMSE = sqrt(mean(error**2)),
            Mean.Abs.Err = mean(abs(error)))



rma_res_esm <- data.frame()
rma_res_line <- data.frame()
rma_res_ssurg <- data.frame()

parse.rma <- function(m){
  data.frame(depth = .depth,
             bias = coef(m),
             se = m$se,
             pval = m$pval,
             tau = sqrt(m$tau2))
}

rma_res_esm <- data.frame()
rma_res_line <- data.frame()
rma_res_ssurg <- data.frame()

for (.depth in c(15, 20,30,40,45)){
  m.esm <- rma(bias, var_est, data = filter(summary1, depth == .depth))
  m.line <- rma(bias_linear, var_est_linear, data = filter(summary1, depth == .depth))
  m.ssurg <- rma(bias_ssurgo, var_est_ssurgo, data = filter(summary1, depth == .depth))

  rma_res_esm <- bind_rows(rma_res_esm, parse.rma(m.esm))


  rma_res_line <- bind_rows(rma_res_line, parse.rma(m.line))

  rma_res_ssurg <- bind_rows(rma_res_ssurg, parse.rma(m.ssurg))
}

formt.strs <- function(rma.res){
  rma.res %>% mutate(
    str_format = sprintf('mean= %s \np=%s\ntau =%s ',  round(bias *100, 2),
                         round(pval, 3), round(tau*100, 3)
    ))}

rma_res_esm <- formt.strs(rma_res_esm)
rma_res_line <- formt.strs(rma_res_line)
rma_res_ssurg <- formt.strs(rma_res_ssurg)



s <- summary1%>%
  ungroup %>%
  group_by(depth) %>%
  summarize(
    tau = sqrt(max(0, var(bias) - mean(var_est))),
    se  = sqrt(1/sum(1/(max(0, var(bias) - mean(var_est))+variance)  ))
  )


interp.theme <- theme(strip.text.x = element_text(size = 15),
                      axis.text.x = element_text(size = 12),
                      axis.title.x = element_text(size = 12),
                      plot.title = element_text(size = 18))

 plot_sp<- loo_data %>%
  filter(depth %in% c(15, 20,30, 45)) %>%
  ggplot(aes(x = error*100, fill = factor(depth))) +
  geom_histogram() +
  facet_wrap(~depth) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('Error, Tons C/ha') + ylab('') +
  ggtitle('Spline Interpolation') +
  scale_fill_discrete(name = "Depth Left Out") +
  geom_text(x = -8, y = 100,
            data = rma_res_esm %>% filter(depth %in% c(15,20,30,45)),
            aes(label = str_format)) +
          guides(fill = F) +
     interp.theme

plot_lin <- loo_data %>%
  filter(depth %in% c(15, 20,30, 45)) %>%
  ggplot(aes(x = linear_error*100, fill = factor(depth))) +
  geom_histogram() +
  facet_wrap(~depth) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('Error, Tons C/ha') + ylab('') +
  ggtitle('Linear Interpolation') +
  scale_fill_discrete(name = "Depth Left Out") +
  geom_text(x = -20, y = 80,
            data = rma_res_line %>% filter(depth %in% c(15,20,30,45)),
            aes(label = str_format)) +
  guides(fill = F) +
  interp.theme

plot_ssurg <- loo_data %>%
  filter(depth %in% c(15, 20,30, 45)) %>%
  ggplot(aes(x = ssurgo_error*100, fill = factor(depth))) +
  geom_histogram() +
  facet_wrap(~depth) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('Error, Tons C/ha') + ylab('') +
  ggtitle('SSurgo Method') +
  scale_fill_discrete(name = "Depth Predicted") +
  geom_text(x = 30, y = 120,
            data = rma_res_ssurg %>% filter(depth %in% c(15,20,30,45)),
            aes(label = str_format)) +
  guides(fill = F) +
  interp.theme

p_err <- plot_grid(plot_sp, plot_lin, plot_ssurg, nrow = 1)
ggsave(file.path('Figures', 'loo_interp.png'), p_err,
       height = 5, width = 13, units ='in' )

loo_data %>%
  filter(depth %in% c(15, 20,30, 45)) %>%
  ggplot(aes(x = error_difference*-100, fill = factor(depth))) +
  geom_histogram() +
  facet_wrap(~depth) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('Reduction of Error, Spline vs Linear') + ylab('') +
  ggtitle('Comparison of Errors between Linear and ESM in Leave-One-Out Cross-Validation') +
  scale_fill_discrete(name = "Depth Left Out") +
  geom_text(x = -7, y = 100,
            data = rma_res_esm %>% filter(depth %in% c(15,20,30,45)),
            aes(label = str_format))

rma_data_all <- rma(bias, var_est, data = summary1)



grouped <- df %>%
  group_by(even_err_bin) %>%
  summarize(mean_sq_err = mean(error**2), mean_pred_err = mean(pred_err))



df$err_pred_adj <- df$pred_err * df$estimated_deriv




##
## Lets be honest, I tested a lot of different error models.
## This is not that important to any specific result.

#m.1 <- lm(abs(error) ~ 0 + pred_err:study, data = df  #%>% filter(pred_err < 15)
#          )
#m.2 <- lm(abs(error) ~ 0 + pred_err_SOC:study, data = df %>% filter(pred_err < 15) )
#m.3 <- lm(log(abs(error)) ~ log(pred_err):study, data = df)
#m.4 <- lm(log(abs(error)) ~   log(pred_err_alt):study, data = df )
#m.5 <- lm(log(abs(error)) ~   0 +log(pred_err_SOC):study + study, data = df )
#m.9 <- lm(log(abs(error_relative)) ~  pred_err:study + study  , data = df %>% filter(pred_err <15))
#m.6 <- lm(log(abs(error)) ~   log(pred_err_SOC_alt):study, data = df )
#m.7 <- lm(abs(error)  ~ 0 + pred_err_alt:study, data =df)
#m.8 <- lm(abs(error) ~ 0 +pred_err_SOC_alt:study, data =df)
#m.10 <- lm(abs(error_relative) ~ 0+pred_err:study +study, data = df )
m.12 <- lm(log(abs(error_relative))~  log(pred_err), data = df)
#m.13 <- lm(log(abs(error_relative))~ 0 + study + log(pred_err), data = df %>% filter(pred_err <10))


#m1 <- lm(error ~ pred_err, data = df)
#m2 <- lm(abs(error) ~ 0 +pred_err, data = df %>% filter(pred_err < 10 & estimated_deriv >0))


#m4 <- lm(abs(error)~ 0 + pred_err_SOC, data = df %>% filter(pred_err < 10 & df$estimated_deriv >0))

#m5 <- lm(log(abs(error)) ~ log(pred_err), data = df)
#m6 <- lm(log(abs(error) ) ~ log(pred_err_SOC), data = df %>% filter(pred_err < 10 & df$estimated_deriv >0 & !grepl('WL|Forest', ID)))
#m7 <- lm(log(abs(error))~pred_err_SOC,  data = df  )


#m8 <- lm(sqrt(abs(error) )~0 + pred_err, data = df  )
#m9 <- lm(sqrt(abs(error))~ 0 + pred_err_SOC, data = df %>% filter(pred_err <10))

#m10 <- lm(log(abs(error)) ~   log(pred_err_alt), data = df )
m11 <- lm(log(abs(error)) ~   log(pred_err_SOC), data = df )
#m12 <- lm(sqrt(abs(error))~ 0 + pred_err_SOC_alt, data = df )

#d <- filter(df, pred_err < 10)
#cor1 <- cor(df$pred_err, df$error)
#cor2 <- cor(df$pred_err, abs(df$error))
#cor3 <- cor(grouped$mean_pred_err**2, grouped$mean_sq_err)
#cor4 <- cor(df$err_pred_adj, abs(df$error))



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




saveRDS(m.12,file = file.path('Results', 'error_model.RDS'))


err_conv_func <- function(soil_sample, x){

  .coef <- coef(m.12)[[2]]
  intercept <- coef(m.12)[[1]]
  err_stat <- exp(log((1/sqrt(sum_weights_mass(soil_sample, x))))*.coef + intercept)
  return(err_stat * calc.linear_SOC_pct(soil_sample, x))
}




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



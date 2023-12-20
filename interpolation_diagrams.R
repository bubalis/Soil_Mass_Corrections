
source('analysis_functions.R')
library(ggplot2)

m <- readRDS('Results/error_model.rds')


err_conv_func <- function(soil_sample, x){

  .coef <- coef(m)[[2]]
  intercept <- coef(m)[[1]]
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




van_doren_samps <- load_study('Van_doren_1986') %>%
  filter(!grepl('t1', ID)) %>%
  calc.cumulative_masses()

devine_samps <- load_study('Devine_2014') %>%
  calc.cumulative_masses()

mish_samps <- load_study('mishra_2010') %>%
  calc.cumulative_masses()

vent_samps <- load_study('Venterea_2006') %>%
  calc.cumulative_masses(
  )

#study_data <- vent_samps
for (study_data  in list(van_doren_samps, devine_samps, mish_samps, vent_samps)){
  study_data$Ref_ID <- 1

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

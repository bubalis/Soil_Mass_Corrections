#####################################
#  Functions for simulating soils for ESM/mass corrections from Fowler et al
#  Author: Ames Fowler              #
#  Date: 8/10/2022                  # 
#####################################
# Guess initial values

library(tidyverse)
library(data.table)

depth = seq(0,80,1)


sample_soil_dube <- function(intervals, soil){
  
  lagged_intervals <- lag(intervals)
  lagged_intervals[1] <- 0
  
  soil$ints <- NA
  soil$ints[1] <- 0
  for (i in seq(1, length(intervals))){
    soil <- soil %>% 
      mutate(
  ints = replace(ints, is.na(ints) & depth > lagged_intervals[i] & depth<= intervals[i],
                i))
    
  }
  soil_sum <- soil %>% filter(!is.na(ints)) %>%
    group_by(ints) %>% 
    summarise(depth_inc = sum(depth_inc),
              Min_mass_soil_MG_ha_inc  = sum(T_soil_mass_Mg_ha_inc*(1-1.9*soc/100)),
              T_soil_mass_Mg_ha_inc = sum(T_soil_mass_Mg_ha_inc),
              SOC_Stock_MG_ha_inc = sum(SOC_Stock_MG_ha_inc),
              BD = T_soil_mass_Mg_ha_inc/depth_inc/10^3,
              soc = SOC_Stock_MG_ha_inc/T_soil_mass_Mg_ha_inc*100,
              scenario = first(scenario),
              period= first(period),
              interval = paste(intervals, collapse = ', ')) %>% ungroup %>% 
    mutate(depth = cumsum(depth_inc),
           Min_mass_soil_MG_ha = cumsum(Min_mass_soil_MG_ha_inc),
           SOC_Stock_MG_ha = cumsum(SOC_Stock_MG_ha_inc),
           T_soil_mass_Mg_ha = cumsum(T_soil_mass_Mg_ha_inc),
           soc_30 = sum(BD*soc)/sum(BD))
  
  return (soil_sum)
}


sample_soil <- function(interval, soil){
  
  # Split soils to relavent parts.
  
  soil0 <- soil[1,]%>% 
    mutate( ints = 0, scenario = first(scenario),
            period= first(period),
            interval = interval)
  soil1 <- soil[2:31,]
  soil3 <- soil[32:81,]
  
  #build intervals if sample depth is less than 30 
  if((30/interval)>1){
    soil1 <- soil1 %>% 
      mutate( ints = cut(depth, breaks =(30/interval)) %>% as.integer())
  }else{
    soil1 <- soil1 %>% 
      mutate( ints = 1)
  }
  
  ## group by intervals (ints) and sum mass values 
  ## 
  soil_sum <- soil1 %>% 
    group_by(ints) %>% 
    summarise(depth_inc = sum(depth_inc),
              Min_mass_soil_MG_ha_inc  = sum(T_soil_mass_Mg_ha_inc*(1-1.9*soc/100)),
              T_soil_mass_Mg_ha_inc = sum(T_soil_mass_Mg_ha_inc),
              SOC_Stock_MG_ha_inc = sum(SOC_Stock_MG_ha_inc),
              BD = T_soil_mass_Mg_ha_inc/depth_inc/10^3,
              soc = SOC_Stock_MG_ha_inc/T_soil_mass_Mg_ha_inc*100,
              scenario = first(scenario),
              period= first(period),
              interval = interval) %>% ungroup %>% 
    mutate(depth = cumsum(depth_inc),
           Min_mass_soil_MG_ha = cumsum(Min_mass_soil_MG_ha_inc),
           SOC_Stock_MG_ha = cumsum(SOC_Stock_MG_ha_inc),
           T_soil_mass_Mg_ha = cumsum(T_soil_mass_Mg_ha_inc),
           soc_30 = sum(BD*soc)/sum(BD))
  
  ## add back in the zero sample for plotting. 
  sum_soil_out <- rbind(soil0 %>% dplyr::select(names(soil_sum)), soil_sum)
  
}


calc.standard.soil.profile <- function(
    # SOC stock percentages from Jobbagy & Jackson 2000
  target_values = c(41, 23, 15, 12, 9),
  # depth vector for Jobbagy & Jackson 2000
    depths = seq(20, 100, 20) ){
 
  
   #note: this is perhaps a starting point for 
  
  v = c(1,.8,-.04)
  
 
  
  
  
 
    
    #function to calculate SOC curve 
    fn <- function(v) {
      soc_o = v[1]
      soc_inf = v[2]
      k = v[3]
      vectorc <- soc_inf+(soc_o-soc_inf)*exp(k*depths)
      
      vectorc <- vectorc/sum(vectorc)*100
      return <- sum(abs(target_values - vectorc))
    }
  
  
  
  #R solver for SOC curve parameters 
  solve_out <- optim(v, fn)
  
  # get SOC curve perameters
  soc_o = solve_out$par[1]
  soc_inf = solve_out$par[2]
  k = solve_out$par[3]
  
  #Build the hypotheical soil profile
  depth = seq(0,80,1)
  soc <- soc_inf+(soc_o-soc_inf)*exp(k*depth)
  soc_0_ratio <- soc[1]/soc[2] 
  
  return (list('SOC' = soc, 'soc_0_ratio' = soc_0_ratio))
}

standard_vals <- calc.standard.soil.profile()
soc_profile <- standard_vals$SOC
soc_0_ratio <-standard_vals$soc_0_ratio




calc_soil <- function(soc, soc_30, soc_linear=FALSE, BD_30, BD_linear ){
  depth_vals = seq(0,80,1)
  
  
  if(soc_linear == FALSE){
    # scale the soc profile curve to soc_30 average
    soc_scale <- soc_30/mean(soc[2:31])
    soc <- soc * soc_scale
  }else{
    soc = soc_30
  }
  
  if(BD_linear == FALSE){
    # build the BD sieres with =+/- 10% 
    BD = BD_30+seq(-1,1,2/29)*BD_30*.1
    W = BD/30
    soc[2:31] = mean(W)/W*soc[2:31]  
    # with changing BD there is transient mass. 
    # To ensure the same 30m mass percentage the SOC 
    # percentage needs to be weighted up with less BD
    soc[1] = soc[2]*ifelse(soc_linear==FALSE,soc_0_ratio,1)
    BD = c(BD[1], BD, rep(BD[30], 50))  
  }else{
    BD = BD_30#1.5    
  }  
  
  dt<- as.data.frame(cbind(depth_vals, BD, soc, soc_linear, BD_linear)) %>% rename(depth = depth_vals)
  
  # Calculate all the soil values from depth, BD, and SOC 
  return (
  calculate.soil.values(dt))
  
  ##temp_dt <- dt %>% mutate(depth_inc = depth - lag(depth) %>% replace_na(0),
  #                         soc_30 = sum(BD[2:31]*soc[2:31])/sum(BD[2:31]),
  #                         T_soil_mass_Mg_ha_inc = depth_inc*BD*10^3/10,
  #                         Min_mass_soil_MG_ha_inc = T_soil_mass_Mg_ha_inc*(1-1.9*soc/100),
  #                         T_soil_mass_Mg_ha   = T_soil_mass_Mg_ha_inc %>% cumsum,
  #                         Min_mass_soil_MG_ha = Min_mass_soil_MG_ha_inc %>% cumsum,
  #                         SOC_Stock_MG_ha_inc = T_soil_mass_Mg_ha_inc*soc/100,
  #                         SOC_Stock_MG_ha     = SOC_Stock_MG_ha_inc%>% cumsum)
  #return (temp_dt)
}

#'
calculate.soil.values <- function(simmed_soil_profile){
  return(
    simmed_soil_profile %>% 
       mutate(depth_inc = depth - lag(depth) %>% replace_na(0),
               soc_30 = sum(BD[2:31]*soc[2:31])/sum(BD[2:31]),
               T_soil_mass_Mg_ha_inc = depth_inc*BD*10^3/10,
               Min_mass_soil_MG_ha_inc = T_soil_mass_Mg_ha_inc*(1-1.9*soc/100),
               T_soil_mass_Mg_ha   = T_soil_mass_Mg_ha_inc %>% cumsum,
               Min_mass_soil_MG_ha = Min_mass_soil_MG_ha_inc %>% cumsum,
               SOC_Stock_MG_ha_inc = T_soil_mass_Mg_ha_inc*soc/100,
               SOC_Stock_MG_ha     = SOC_Stock_MG_ha_inc%>% cumsum)
    
  )
  
}


# Funct. expand soil ---- 


#'From fowler et al. expand a soil to simulate a reduction in bulk density.
#'
#'This function could be improved in two ways: 
#'1: it needs to be changed so that it behaves correctly when simulating and INCREASE in BD
#'2 : it needs to be changed so that soc_new does not need to be passed as an argument,
#'rather that soc level to 30 cm is automatically kept constant.
expand_soil <- function(dt, BD_new, soc_new){
  
  #explained soil depth by type 
  # BD_new = 1.1
  BD_mean = mean(dt$BD[2:31])
  #calculate the depth of the old soil that fills 30cm with new BD average.
  DA  = 30 *BD_new/BD_mean  
  
  if( dt$BD_linear %>% first == 0){
    # calc percent change. 
    prct_change = BD_mean/BD_new - 1 
    # build the new BD sequence 
    BD_seq = BD_new*(1-seq(prct_change, -prct_change, -prct_change*2/(DA)))
    # add in the zero point and values greater than 30   
    BD_seq_out <- c(BD_seq[1]-(BD_seq[2] -BD_seq[1]), BD_seq,rep(BD_seq[length(BD_seq)], ceiling(80-DA-1)))
    
  }else{
    # If not varying with depth... all is constant. 
    BD_seq= rep(BD_new, 81)
    BD_seq_out= rep(BD_new, 81)
  }
  
  # calculate the depth series of teh expanded 1cm increments, now some value 1+cm for each. 
  Depth_out = c(0,dt$BD[2:length(dt$BD)]/BD_seq_out[2:length(dt$BD)]) %>% cumsum
  
  # interpolate BD values for new profile (back to 1:30 depth by 1cm)
  temp <- approx(x = Depth_out, y = BD_seq_out, xout = dt$depth, method = "linear",
                 yright =BD_seq[length(BD_seq)], ties = mean, na.rm = TRUE) %>% 
    as.data.frame()
  
  names(temp) <- c("depth", "BD")
  # mean(temp$BD[2:31]) ## for checking
  
  
  # interpolate SOC values for new profile (back to 1:30 depth by 1cm)    
  temp2 <- approx(x = Depth_out[1:81], y = dt$soc, xout = dt$depth, method = "linear",
                  yright =BD_seq[length(BD_seq)], ties = mean, na.rm = TRUE) %>% 
    as.data.frame()
  names(temp2) <- c("depth", "soc")
  
  ### merge all
  temp = merge(temp, temp2) %>% cbind( dt[,c("soc_linear", "BD_linear", "scenario")])
  
  #sum(temp$soc[2:31]*temp$BD[2:31])/sum(temp$BD[2:31]) ###check average SOC 
  
  #correct SOC to fixed level
  if(first(temp$soc_linear) ==0){
    
    ## find the net change in SOC stock from the expanded soil from the 30cm Base line 
    net_soc_stock = soc_new*mean(temp$BD[2:31])*30 - sum(temp$soc[2:31]*temp$BD[2:31]) 
    
    ## weight by SOC stock by SOC distribution and 
    ## multiple by BD to convert back to SOC concentration 
    net_soc = net_soc_stock*temp$soc[2:31]/mean(temp$soc[2:31])/30/temp$BD[2:31]
    
    ## add the increase in SOC and fix the zero value 
    temp$soc[2:31] = c((net_soc))+temp$soc[2:31] #,rep(0,20)
    temp$soc[1] = temp$soc[2]*ifelse(first(temp$soc_linear)==0,soc_0_ratio,1)
    
    # sum(temp$soc[2:31]*temp$BD[2:31])/sum(temp$BD[2:31])  ## check average 
  }else{
    temp$soc = soc_new    
  }
  
  #recalucate the soil variables. 
  return( calculate.soil.values(temp))
  
  #temp_dt <- temp %>% mutate(depth_inc = depth - lag(depth) %>% replace_na(0),
  #                           soc_30 = sum(BD[2:31]*soc[2:31])/sum(BD[2:31]),
  #                           T_soil_mass_Mg_ha_inc = depth_inc*BD*10^3/10,
  #                           Min_mass_soil_MG_ha_inc = T_soil_mass_Mg_ha_inc*(1-1.9*soc/100),
  #                           T_soil_mass_Mg_ha   = T_soil_mass_Mg_ha_inc %>% cumsum,
  #                           Min_mass_soil_MG_ha = Min_mass_soil_MG_ha_inc %>% cumsum,
  #                           SOC_Stock_MG_ha_inc = T_soil_mass_Mg_ha_inc*soc/100,
  #                           SOC_Stock_MG_ha     = SOC_Stock_MG_ha_inc%>% cumsum)
  
  
}

#' Take a simulated soil with data for each 1-inch increment and alter it so that 
#' it is composed of multiple genetic horizons that are completely homogenous
#' @param simmed_soil simulated soil produced by function calc_soil
#' @param horizon_depths a vector of integer values of the bottoms of the horizon depths.
simulate.homogenous.horizons <- function(simmed_soil, horizon_depths){
  lagged_depths <- lag(horizon_depths)
  lagged_depths[1] <- 0
  
  for (i in seq(1, length(horizon_depths))){
    simmed_soil <- simmed_soil %>% 
      mutate(
        soc = replace(soc,  depth > lagged_depths[i] & depth<= horizon_depths[i],
        simmed_soil %>% filter(depth > lagged_depths[i] & depth<= horizon_depths[i])%>%
                        summarize(weighted.mean(soc, BD))) %>% unlist() ,
        BD = replace(BD, depth > lagged_depths[i] & depth<= horizon_depths[i],
            simmed_soil %>% filter(depth > lagged_depths[i] & depth<= horizon_depths[i])%>%
         pull(BD) %>% mean()))
}
return (calculate.soil.values(simmed_soil))
}

smoothstep <- function(x, edge0, edge1){
  x <- clamp((x-edge0) /  (edge1-edge0))
  return((3*x**2 - 2* x**3))
  
  
}


clamp <-function(x, lowerlimit = 0, upperlimit =1){
  if (x> upperlimit){return(upperlimit)
  }else if(x< lowerlimit){return(lowerlimit)
    }else{return(x)}
}

sim.horiz.smoothstep <- function(horizon){
  min_soc <- min(horizon$soc)
  max_soc <- max(horizon$soc)
  
  ys <- sapply(horizon$depth, 
  function(x){smoothstep(x, min(horizon$depth), max(horizon$depth)) }) %>% rev()
  
  soil_c <- ys * (max_soc-min_soc) + min_soc
  return(soil_c*  weighted.mean(horizon$soc, horizon$BD) /weighted.mean(soil_c, horizon$BD))  
}


simulate.smoothstep.horizons <- function(simmed_soil, horizon_depths){
  lagged_depths <- lag(horizon_depths)
  lagged_depths[1] <- 0
   
  
   for (i in seq(1, length(horizon_depths))){
    simmed_soil <- simmed_soil %>% 
      mutate(
        soc = replace(soc,  depth > lagged_depths[i] & depth<= horizon_depths[i],
                      simmed_soil %>% filter(depth > lagged_depths[i] & depth<= horizon_depths[i])%>%
                        sim.horiz.smoothstep()))
                        
                        
  
   }
  return (calculate.soil.values(simmed_soil))
}

get_correct_SOC_change <- function(intial_soil, new_soil, depth_to = 30){
  xout <- intial_soil %>% filter(depth == depth_to) %>% pull(Min_mass_soil_MG_ha)
  interp <- approx(new_soil$Min_mass_soil_MG_ha, new_soil$SOC_Stock_MG_ha, xout = xout )
  return(interp$y)
}

get_true_soc_change <- function(t_dt){
  all_res <- data.frame()
  for (ref_point in c('t0', 't1', 't2')){
  ref_points <- t_dt %>% filter(period ==  ref_point)
   new_measures <- t_dt %>% filter(period != ref_point)
   distinct.groups <- new_measures %>% distinct(scenario, period)
   out <- c()
   names_ <- c()
   ref_IDs = c()
   for (i in seq(1, nrow(distinct.groups))){
     grp <- distinct.groups[i,]
     new_soil_dat <- new_measures %>% filter(scenario == grp$scenario & period == grp$period)
     baseline <- ref_points %>% filter(scenario == grp$scenario)
     baseline_soc <- baseline %>% filter(depth == 30) %>% pull(SOC_Stock_MG_ha)
     ref_IDs <- ref_IDs %>% append( baseline %>% mutate(ID = paste(scenario, period, sep = '_')) 
                                    %>% pull(ID) %>% first())
     
     true_delta_SOC <- get_correct_SOC_change(baseline, new_soil_dat) - baseline_soc
     out <- append(out, true_delta_SOC)
     names_ <- append(names_, paste(grp$scenario, grp$period, sep = '_'))
   }
  
  all_res <- rbind(all_res, data.frame(true_change_soc = out, ID =  names_,
                                       Ref_ID = ref_IDs))
  
}
  return(all_res)
}

simulate.soil.Fowler <-function(soc_profile){
  dt_s1_t0 <- calc_soil (soc= soc_profile, soc_30= 1.4, soc_linear=TRUE, BD_30=1.5, BD_linear=TRUE ) %>% mutate(scenario = "s1", period= "t0")
  dt_s2_t0 <- calc_soil (soc= soc_profile, soc_30= 1.4, soc_linear=FALSE, BD_30=1.5, BD_linear=FALSE ) %>% mutate(scenario = "s2", period= "t0")
  #homogenous horizons:
  dt_s3_t0 <- simulate.homogenous.horizons(dt_s2_t0, c(20,40, 60,80)) %>% mutate(scenario = "s3", period= "t0")
  dt_s4_t0 <- simulate.homogenous.horizons(dt_s2_t0, c(15, 30, 50, 80)) %>% mutate(scenario = "s4", period= "t0")
  #smoothstep horizons:
  dt_s5_t0 <- simulate.smoothstep.horizons(dt_s2_t0, c(20,40, 60,80)) %>% mutate(scenario = "s5", period= "t0")
  dt_s6_t0 <- simulate.smoothstep.horizons(dt_s2_t0, c(15, 30, 50, 80)) %>% mutate(scenario = "s6", period= "t0")
  
  
  # Expand base soil t1
  dt_s1_t1 <- expand_soil(dt_s1_t0, BD_new = 1.3, soc_new = 1.5) %>% mutate(scenario = "s1", period= "t1")
  dt_s2_t1 <- expand_soil(dt_s2_t0, BD_new = 1.3, soc_new = 1.5) %>% mutate(scenario = "s2", period= "t1")
  dt_s3_t1 <- expand_soil(dt_s3_t0, BD_new = 1.3, soc_new = 1.5) %>% mutate(scenario = "s3", period= "t1")
  dt_s4_t1 <- expand_soil(dt_s4_t0, BD_new = 1.3, soc_new = 1.5) %>% mutate(scenario = "s4", period= "t1")
  dt_s5_t1 <- expand_soil(dt_s5_t0, BD_new = 1.3, soc_new = 1.5) %>% mutate(scenario = "s5", period= "t1")
  dt_s6_t1 <- expand_soil(dt_s6_t0, BD_new = 1.3, soc_new = 1.5) %>% mutate(scenario = "s6", period= "t1")
  
  # Expand base soil t2
  dt_s1_t2 <- expand_soil(dt_s1_t1, BD_new = 1.1, soc_new = 1.6) %>% mutate(scenario = "s1", period= "t2")
  dt_s2_t2 <- expand_soil(dt_s2_t1, BD_new = 1.1, soc_new = 1.6) %>% mutate(scenario = "s2", period= "t2")
  dt_s3_t2 <- expand_soil(dt_s3_t0, BD_new = 1.1, soc_new = 1.6) %>% mutate(scenario = "s3", period= "t2")
  dt_s4_t2 <- expand_soil(dt_s4_t0, BD_new = 1.1, soc_new = 1.6) %>% mutate(scenario = "s4", period= "t2")
  dt_s5_t2 <- expand_soil(dt_s5_t0, BD_new = 1.1, soc_new = 1.6) %>% mutate(scenario = "s5", period= "t2")
  dt_s6_t2 <- expand_soil(dt_s6_t0, BD_new = 1.1, soc_new = 1.6) %>% mutate(scenario = "s6", period= "t2")
  
  # Combine all expanded base soil and name.  
  t_list <- list(dt_s1_t0, dt_s2_t0, dt_s1_t1, dt_s2_t1, dt_s1_t2, dt_s2_t2,
                 dt_s3_t0, dt_s4_t0, dt_s3_t1, dt_s4_t1, dt_s3_t2, dt_s4_t2,
                 dt_s5_t0, dt_s6_t0, dt_s5_t1, dt_s6_t1, dt_s5_t2, dt_s6_t2
                 
                 )
  names(t_list) <- c("dt_s1_t0", "dt_s2_t0",  "dt_s1_t1", "dt_s2_t1", "dt_s1_t2", "dt_s2_t2",
                     'dt_s3_t0', 'dt_s4_t0', 'dt_s3_t1', 'dt_s4_t1', 'dt_s3_t2', 'dt_s4_t2',
                     'dt_s5_t0', 'dt_s6_t0', 'dt_s5_t1', 'dt_s6_t1', 'dt_s5_t2', 'dt_s6_t2'
                     )
  
  
  t_dt <- rbindlist(t_list,use.names=TRUE)
  
  
  soc_delt_actual <- get_true_soc_change(t_dt)
  soc_delt_actual <- soc_delt_actual %>% filter(!(grepl('t2', Ref_ID) & grepl('t1', ID))) %>% 
    mutate(ID = ifelse(grepl('t2', Ref_ID) & grepl('t0', ID), paste(ID, 'x', sep = '_'), ID)) %>%
    select(ID, Ref_ID, true_change_soc)
  
  t_dt <- t_dt %>% mutate(Lower_cm = depth, Upper_cm = lag(depth), 
                  SOM_pct = soc * 1.9, 
                  ID = paste(scenario, period, sep = '_' )) %>% 
                  rename(SOC_pct = soc, BD_g_cm3 = BD) %>% 
    filter(!is.na(Upper_cm)) %>% select(ID, SOM_pct, SOC_pct,
                                        BD_g_cm3, Upper_cm, Lower_cm, period, 
                                        scenario) %>% mutate(Rep = 1) %>% 
                                        filter(Lower_cm >0)
  
  ref_keys <- t_dt %>% filter(period == 't0') %>% 
    select(scenario,  ID) %>% rename( Ref_ID = ID  ) %>% distinct()
  
  t_dt <- merge(t_dt, ref_keys)
  
  
  # # Sample all soils and combine 
  # sample_10cm <- map(t_list,function(H) sample_soil(interval = 10, H))
  # sample_15cm <- map(t_list,function(H) sample_soil(interval = 15, H))
  # sample_30cm <- map(t_list,function(H) sample_soil(interval = 30, H))
  # sample_4_depth <- map(t_list, function(H) sample_soil_dube(intervals = c(10, 30,50, 70), H))
  # #sample_30_50 <- sample_30_50 <- map(t_list, function(H) sample_soil_dube(intervals = c(30,50), H))
  # sample_close <- map(t_list, function(H) sample_soil_dube(intervals = c(24, 30, 36), H))
  # 
  # sample_list <- c(sample_10cm, sample_15cm, sample_30cm, sample_4_depth,
  #                  sample_close)
  # sample_dt <- rbindlist(sample_list)
  # 
  # sample_dt <- sample_dt %>% mutate(ID = paste(scenario, period, 
  #                                               sep = '_'),
  #                                   Rep = 1)
  # 
  # #added by me: add a unique ID and Ref ID to  
  # ref_keys <- sample_dt %>% filter(period == 't0') %>% 
  #   select(scenario, interval, ID) %>% rename( Ref_ID = ID  ) %>% distinct()
  # 
  # sample_dt <- merge(sample_dt, ref_keys)
  #              
  # 
  # 
  # 
  # return (list('sample_dt' = sample_dt %>% filter(depth_inc >0),
  # 'true_delta_soc' =   soc_delt_actual)
  #              )
  
  return(list('t_dt' = t_dt,
              'true_delta_soc' =   soc_delt_actual)
  )
  
}




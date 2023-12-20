##
## Using the 2-depth method described by Wendt & Hauser (2013)
##



#source('util.R')



#' Simulate a 2-depth sample that is split into two, weighted separately,
#' then mixed to acheive the correct equivalent mass, then has a single SOC test conducted on it.
#' As described in Wendt & Hauser (2013).
sim_joined_sample <- function(soil_sample, initial_soil, depth_of_quant = 30){
  mass_of_quantification <- initial_soil %>% filter(Lower_cm == depth_of_quant) %>%
    pull(Cum_Soil_g_cm2) %>%
    first()
  #print(mass_of_quantification)
  # Simulate a 2-depth sample that is split into two, weighted separately,
  # then mixed to acheive the correct equivalent mass, then has a single SOC test conducted on it
  masses <- soil_sample %>% pull(Soil_g_cm2)

  weights <- c(masses[1],   mass_of_quantification - masses[1] ) / masses[1]
  if (masses[1] > mass_of_quantification){
    weights <- c(1,0)
  }

  pct <- weighted.mean(soil_sample$SOC_pct, weights)
  #print(pct)
  yout <-  pct * mass_of_quantification /100
  print(yout)
  return (                soil_sample %>%
                            first() %>%
                            select(-c(ID, Rep, Ref_ID, Cum_SOC_g_cm2)) %>%
                            mutate(
                              Upper_cm = 0,
                              Lower_cm = depth_of_quant) %>%
                            mutate(Cum_SOC_g_cm2 =yout)
  )

}



group_run_joined_sample <- function(soil_data, depths, depth_of_quant = 30){
  if (length(depths) != 2) {stop('Can only simulate 2-depth cores if depths argument is length 2')}
  t2.samps <- filter(soil_data, ID != Ref_ID) %>%
    df_agg_to_depths(depths) %>%
    calc.cumulative_masses()

  t1.samps <- filter(soil_data, ID == Ref_ID) %>%
    calc.cumulative_masses()

  fn <- function(sample){
    intial_soil <- subset(t1.samps,
                          Ref_ID==first(sample$Ref_ID) &
                            Rep==first(sample$Rep))
    return (
      sim_joined_sample(sample, intial_soil,
                        depth_of_quant = depth_of_quant))

  }
  return(t2.samps %>%
           group_by(ID, Ref_ID, Rep) %>%
           group_modify(~fn(.), .keep = T))
}


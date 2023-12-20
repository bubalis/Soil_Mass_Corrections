
#
# ESM on two-depth cores using a exponential decay function
#


#source('util.R')






#'get decay parameters for the exponential decay model of SOC stock with depth
#'Give the exact fit for 2-depth data.
solve_decay_params_2depth <- function (x1, x2,  y1,y2){
  func <- function(b){(exp(-b * x2) -1) / (exp(-b* x1) -1) * y1 - y2 }

  b <- uniroot(func, lower = -1, upper = 1)$root

  a <- y1/ (exp(-b* x1 ) -1) * -b
  return(list( 'a' =a, 'b' = b))
}



#For a two-depth sample, use the exponential decay function to calculate the equivalent soil mass.
#Because there are only two observations and two parameters, one version of the curve will fit perfectly.
#This can also be used for a single depth sample

decay_ESM_two_depth <- function(sample, intial_soil, depth_of_quant){
  sample <- calc.cumulative_masses(sample)
  intial_soil <- calc.cumulative_masses(intial_soil)

  #if there is only one depth in the sample
  #add a zero row
  if (nrow(sample) == 1){
    row1 <- make_zerodepth_row(sample)
    row2 <- sample %>% filter(Lower_cm == max(Lower_cm))
  }else if(nrow(sample == 2)){
    row1 <- sample %>% filter(Upper_cm == 0)
    row2 <- sample %>% filter(Lower_cm == max(Lower_cm))
  }else{
    stop('2-depth decay fitting only works with 2 depths')
  }

  mass1 <- intial_soil %>% filter(Lower_cm == depth_of_quant) %>% pull(Cum_Min_Soil_g_cm2)
  mass2 <- sample %>% filter(Lower_cm == depth_of_quant) %>% pull(Cum_Min_Soil_g_cm2)

  bd_at_quant_mass <- sample %>% filter(Cum_Min_Soil_g_cm2 < depth_of_quant) %>%
    pull(BD_g_cm3) %>% last()


  if (length(bd_at_quant_mass) == 0){
    bd_at_quant_mass <- sample %>% pull(BD_g_cm3) %>% first()
  }


  decay_params <- solve_decay_params_2depth(row1 %>% pull(Lower_cm),
                                            row2 %>% pull(Lower_cm),
                                            row1 %>% pull(Cum_SOC_g_cm2),
                                            row2 %>% pull(Cum_SOC_g_cm2))
  a <- decay_params$a
  b <- decay_params$b

  depth_adj <- (mass1  - mass2) / bd_at_quant_mass

  yout <- -(a/b) * (exp(-b * (depth_of_quant + depth_adj))-1)

  out <- cbind(data.frame(Cum_SOC_g_cm2 =yout),
               data.frame(row2) %>%
                 select(-c(ID, Rep, Ref_ID, Cum_SOC_g_cm2)))

  return(out)
}


run_decay_ESM <- function(filtered_FD, depth_of_quant){
  t2.samps <- filter(filtered_FD, ID != Ref_ID)
  t1.samps <- filter(filtered_FD, ID == Ref_ID)

  #print(nrow(t2.samps))
  #print(nrow(t1.samps))

  fn <- function(sample){
    #print(sample)
    intial_soil <- subset(t1.samps,
                          Ref_ID==first(sample$Ref_ID) &
                            Rep==first(sample$Rep))
    return (
      decay_ESM_two_depth(sample, intial_soil,
                          depth_of_quant = depth_of_quant))
  }


  return (t2.samps %>%
            group_by(ID, Ref_ID, Rep) %>%
            group_modify(~fn(.), .keep = T)
  )}


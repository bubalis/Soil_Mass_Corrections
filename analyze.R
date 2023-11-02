options(warn = -1) 
#library(constrKriging)
library(stringr)
source('analysis_functions.R')


load_coefs <- function(filepath){
  temp <- read.csv(filepath)
  
  setNames(temp$x, temp$X)
}

m.12.coefs <- load_coefs('model_coefficients.csv')


#error_prop_func <- function(soil_sample, mass_of_est){
#  1/sqrt(sum_weights_mass(soil_sample, mass_of_est)) * calc.linear_SOCpct_mass(soil_sample, mass_of_est) * .0015
#  
#}




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
source("data_sim.R")
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
            'Blanco-canqui_2011',
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
results <- NULL

out_s4t <- out_s4t %>% 
  distinct(ID, Ref_ID, study_name, method, sample_depths, .keep_all = T)


out_s4t %>% #select( ID, Ref_ID, Rep, Cum_SOC_g_cm2, method, sample_depths, 
  #        Cum_SOC_g_cm2_baseline, soc_change,  study_name) %>% 
  #distinct() %>% 
  write.csv( file.path('Results', 'all_data_space_for_time_comparison.csv'))


out_s4t <- rbind_memoryerror(lapply(unique(out_s4t$study_name), function(x){
  add_reference_changes(x, out_s4t %>% distinct(ID, Ref_ID, sample_depths, method, study_name, .keep_all = T))}))









ag_res <- out_s4t %>% 
  filter(!grepl('WL', ID) & ID!= "Forest Succession" & !grepl('WL', Ref_ID) & Ref_ID != "Forest Succession") 

summary_2 <- ag_res %>% group_by(study_name, method, sample_depths) %>% summarize_res()

ag_res %>%  
  group_by(method, sample_depths) %>% 
  summarize_res() %>% 
  filter(sample_depths %in% c('30, 40', '30, 45', '30, 50', '30, 60', '20, 30', '30'))


filtered_res <- filter(ag_res, str_count(sample_depths, ',') < 2)

filtered_res <- filtered_res %>% 
  mutate(strategy = paste(method, sample_depths, ', '))

strategies <- distinct(filtered_res, method, sample_depths, strategy) %>% 
  mutate(other_depth = str_remove(sample_depths, '30,|, 30'))  %>% 
  mutate(other_depth = ifelse(!grepl(',', sample_depths), 0, as.numeric(other_depth)))


fit_data2 <-  filtered_res %>%
  mutate(strategy = as.factor(strategy)) %>% 
  mutate(strategy = relevel(strategy, "Fixed Depth 30 , "))
  
m.15 <- lmer(log(relative_error) ~ strategy + (1|study_name), data = fit_data2)
m.8 <- lm(log(relative_error) ~ strategy, data = fit_data2)

fit_data <- filtered_res %>% 
  group_by(ID, Ref_ID, study_name, strategy) %>% 
  summarize_res() %>% 
  ungroup() %>% 
  mutate(strategy = as.factor(strategy)) %>% 
  mutate(strategy = relevel(strategy, 'Fixed Depth 30 , '))



m3 <- lm(RMSE ~  as.factor(strategy) + as.factor(study_name), data = fit_data)
m4 <- lm(abs(bias) ~ as.factor(strategy) + as.factor(study_name), data = fit_data)
m3.1 <- lmer(RMSE ~ as.factor(strategy) + (1 | as.factor(study_name), data = fit_data))


#'Upscale factor accounts for 
reg_results_dataframe <- function(model, upscale_factor = 200){
  coefs  <- coef(summary(model))
  d <- data.frame(coef. = coefs[,'Estimate'], se =  coefs[, "Std. Error"] * sqrt(upscale_factor) ) 
  d <- d %>% 
    mutate(strategy = rownames(d)) %>%
    mutate(strategy = str_remove( strategy, 'strategy')) %>%
    filter(!(grepl('interecept',strategy) |(grepl('study_name', strategy) ) )) %>%
    mutate(strategy = gsub('as.factor\\(strategy\\)','', strategy)) %>% merge(strategies) %>%
    mutate(other_depth = str_remove(sample_depths, '30,|, 30'))
      
  return(d)
}

df_logratio <- reg_results_dataframe(m.15)

data <-  fit_data  %>% merge(strategies) %>%
              filter(other_depth >9 & other_depth <71)
                     
ggplot(data =data,  aes(x = factor(other_depth), y = bias, fill = method, color = method)) + 
  geom_boxplot( position = position_dodge()) + 
  geom_hline(yintercept = 0)


df_rmse <- reg_results_dataframe(m3)
df_bias <- reg_results_dataframe(m4)

data_1 <- filter(df_rmse, !grepl(',', other_depth)) %>% 
  mutate(other_depth = ifelse(!grepl(',', sample_depths), 0, as.numeric(other_depth)))

data_2 <- filter(df_bias, !grepl(',', other_depth)) %>% 
    mutate(other_depth = ifelse(!grepl(',', sample_depths), 0, as.numeric(other_depth)))

df_logratio <- filter(df_logratio, !grepl(',', other_depth)) %>% 
  mutate(other_depth = ifelse(!grepl(',', sample_depths), 0, as.numeric(other_depth)))

df_logratio <- df_logratio %>% mutate(rel_err_fd =  exp(coef.), lo = exp(coef. - se * 1.96),
                                      hi = exp(coef. + se *1.96))




df_logratio %>% filter((other_depth >9 & other_depth < 71  & !grepl('Decay', method)) | 
                    (method %in% c('Linear', 'Mass Correction, SSurgo avg', 'Fixed Depth', 'Linear Average, Fowler Weights') & 
                       other_depth == 0) ) %>% 
  mutate(method = factor(method, 
                         levels = c('ESM Spline', 'Two-depth combined sample', 
                                    'Linear Average, Fowler Weights', "Linear",
                                    'Mass Correction, SSurgo avg'))) %>%
  ggplot(aes(x = other_depth, y= rel_err_fd, group = method, color = method)) + 
  geom_point(position = position_dodge(width = 2.5)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), 
                position = position_dodge(width = 2.5)) +
  ylab("Relative Error (fraction), compared with fixed-depth approach") + 
  ggtitle('Relative Error of sampling-strategies\nand quantification methods for quantifying carbon to 30 cm') +
  scale_x_continuous(name = 'Depth of Other Sample', 
                     breaks = c(0, 10, 15, 20,25,35, 40,45, 50, 60, 70), 
                     labels = list('Single Depth', 10, 15,20,25,35,40,45, 50, 60, 70)
  )

ggplot(data = data_1 %>% filter(other_depth > 3 & other_depth <30), 
       aes(x = other_depth, y= coef., group = method, color = method)) + 
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = coef. - se * 1.96, ymax = coef. + se* 1.96 ), 
                position = position_dodge(width = 1)) +
  xlab('Depth of Other Sample')+ ylab("Relative Error, g/cm2") + 
  ggtitle(
'Relative Error of 2-depth sampling-strategies
and quantification methods for quantifying carbon to 30 cm')

ggplot(data = data_1 %>% filter(other_depth >30), 
       aes(x = other_depth, y= coef., group = method, color = method)) + 
  geom_point(position = position_dodge(width = 2.5)) +
  geom_errorbar(aes(ymin = coef. - se, ymax = coef. + se ), 
                position = position_dodge(width = 2.5)) +
  xlab('Depth of Other Sample')+ ylab("Relative Error, g/cm2") + 
  ggtitle(
'Relative Error of 2-depth sampling-strategies
and quantification methods for quantifying carbon to 30 cm')



#I actually could plot this data raw? 
ggplot(data = data_1 %>% filter(other_depth ==0 & !grepl('SOC', method) & !grepl('ESM', method))  ,  
       aes(x = 0, y= coef., group = method, color = method)) + 
  geom_point(position = position_dodge(width = .8)) +
  geom_errorbar(aes(ymin = coef. - se, ymax = coef. + se ), position = position_dodge(width = .8)) + 
  ylab("Relative Error, g/cm2") + 
  ggtitle(
'Relative Error of different methods for quantifying soil carbon on a single 30 cm sample')


data_1 %>% filter((other_depth >9 & other_depth < 71  & !grepl('Decay', method)) | 
                    (method %in% c('Linear', 'Mass Correction, SSurgo avg', 'Fixed Depth', 'Linear Average, Fowler Weights') & 
                  other_depth == 0) ) %>% 
  mutate(method = factor(method, 
        levels = c('ESM Spline', 'Two-depth combined sample', 
                   'Linear Average, Fowler Weights', "Linear",
                    'Mass Correction, SSurgo avg'))) %>%
ggplot(aes(x = other_depth, y= coef., group = method, color = method)) + 
  geom_point(position = position_dodge(width = 2.5)) +
  geom_errorbar(aes(ymin = coef. - se * 1.96, ymax = coef. + se * 1.96 ), 
                position = position_dodge(width = 2.5)) +
   ylab("Relative Error, g/cm2, compared with fixed-depth approach") + 
  ggtitle('Relative Error of sampling-strategies\nand quantification methods for quantifying carbon to 30 cm') +
  scale_x_continuous(name = 'Depth of Other Sample', 
                     breaks = c(0, 10, 15, 20,25,35, 40,45, 50, 60, 70), 
                     labels = list('Single Depth', 10, 15,20,25,35,40,45, 50, 60, 70)
                   )

ggplot(data = fit_data %>% merge(strategies) %>% 
         filter((other_depth <30 & other_depth >9)),
       aes(x = factor(method), y = bias, fill = method)) + geom_boxplot() + 
  xlab('Quantification Method') + ylab('bias, gC / cm2') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ggtitle('Average Bias for each point quantified, 2-depths cores, 30 cm lower depth')



ggplot(data = fit_data %>% merge(strategies) %>% 
         filter((other_depth <30 & other_depth >9)),
       aes(x = factor(method), y = bias_rate, fill = method)) + geom_boxplot() + 
  xlab('Quantification Method') + ylab('Bias rate, Tons C overestimated per ton increase in fixed-depth sample') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ggtitle('Average bias rate of each point quantified,  2-depths cores, 30 cm lower depth')


ggplot(data = fit_data %>% merge(strategies) %>%  filter((other_depth <30 & other_depth >9) ), 
      aes(x = factor(method), y = RMSE, fill = method)) + 
  geom_boxplot() + xlab('Quantification Method') + ylab('RMSE, gC / cm2')+ 
  ggtitle('RMSE for each point quantified, 2-depths cores, 30 cm lower depth')


ggplot(data = fit_data %>% merge(strategies) %>%  
         filter(other_depth >30 & other_depth <75 | method == 'Fixed Depth'), 
       aes(x = factor(method), y = RMSE, fill = method)) + geom_boxplot() + 
  xlab('Quantification Method') + ylab('RMSE, gC / cm2') +
  ggtitle('RMSE of each point quantified, 2-depth cores lower depth below 30 cm')


ggplot(data = fit_data %>% merge(strategies) %>%  filter(other_depth >30 & other_depth <75 ), 
       aes(x = factor(method), y = bias, fill = method)) + geom_boxplot() + 
  xlab('Quantification Method') + ylab('Bias, gC / cm2') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ggtitle('Average bias of each point quantified, 2-depth cores, lower depth below 30 cm')


ggplot(data = fit_data %>% merge(strategies) %>%  filter(other_depth >30 & other_depth <75 ), 
       aes(x = factor(method), y = bias_rate, fill = method)) + geom_violin() + 
  xlab('Quantification Method') + ylab('Bias rate, Tons C overestimated per ton increase in fixed-depth sample') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ggtitle('Average bias rate of each point quantified, 2-depth cores, lower depth beneath 30 cm')


ggplot(data = fit_data %>% merge(strategies) %>% filter(other_depth ==0 & !grepl('ESM', method)), 
       aes(x = factor(method), y = bias, fill = method)) + geom_violin() + 
  xlab('Quantification Method') + ylab('Bias, gC / cm2') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ggtitle('Average bias of each point quantified, single-depth cores')


ggplot(data = fit_data %>% merge(strategies) %>% filter(other_depth ==0 & !grepl('ESM', method)), 
       aes(x = factor(method), y = RMSE, fill = method)) + geom_boxplot() + 
  xlab('Quantification Method') + ylab('RMSE, gC / cm2') + 
  ggtitle('Average RMSE of points quantified, g C /cm2, single depth cores')


ggplot(data = fit_data %>% merge(strategies) %>% filter(other_depth ==0 & !grepl('ESM', method)), 
       aes(x = factor(method), y = bias_rate, fill = method)) + geom_boxplot() + 
  xlab('Quantification Method') + ylab('Bias rate, Tons C overestimated per ton increase in fixed-depth sample') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  ggtitle('Average bias rate of each point quantified, single-depth cores')

# 
# 
# ggplot(data = data_2 %>% filter(other_depth > 8 & other_depth <30), 
#        aes(x = other_depth, y= coef., group = method, color = method)) + 
#   geom_point(position = position_dodge(width = 1)) +
#   geom_errorbar(aes(ymin = coef. - se, ymax = coef. + se ), position = position_dodge(width = 1)) +
#   geom_hline(yintercept = 0 , color = 'black') +
#   xlab('Depth of Other Sample')+ ylab("Relative Bias, g/cm2") + 
#   ggtitle('Relative Bias of 2-depth sampling-strategies\nand quantification methods for quantifying carbon to 30 cm')
# 
# ggplot(data = data_2 %>% filter(other_depth >30), 
#        aes(x = other_depth, y= coef., group = method, color = method)) + 
#   geom_point(position = position_dodge(width = 2)) +
#   geom_errorbar(aes(ymin = coef. - se, ymax = coef. + se ), position = position_dodge(width = 2)) +
#   geom_hline(yintercept = 0, color = 'black') +
#   xlab('Depth of Other Sample')+ ylab("Relative Bias, g/cm2") + 
#   ggtitle('Relative Bias of 2-depth sampling-strategies\nand quantificaiton methods for quantifying carbon to 30 cm')
# 
# ggplot(data = data_2 %>% filter(other_depth >9 & other_depth < 71),
#        aes(x = other_depth, y= coef., group = method, color = method)) + 
#   geom_point(position = position_dodge(width = 2)) +
#   geom_errorbar(aes(ymin = coef. - se, ymax = coef. + se ), position = position_dodge(width = 2)) +
#   geom_hline(yintercept = 0 , color = 'black') +
#   xlab('Depth of Other Sample')+ ylab("Relative Bias, g/cm2") + 
#   ggtitle('Relative Bias of 2-depth sampling-strategies\nand quantification methods for quantifying carbon to 30 cm')
# 
# 
# #I actually could plot this data raw? 
# ggplot(data = data_2 %>% filter(other_depth ==0 & !grepl('SOC', method)),  
#        aes(x = other_depth, y= coef., group = method, color = method)) + 
#   geom_point(position = position_dodge(width = .8)) +
#   geom_errorbar(aes(ymin = coef. - se, ymax = coef. + se ), position = position_dodge(width = .8)) + 
#   geom_hline(yintercept = 0, color = 'black') +
#   ylab("Relative Bias, g/cm2") + 
#   ggtitle('Relative Bias of different methods for quantifying soil carbon on a single 30 cm sample')




bias_difs <- data.frame()
error_difs <- data.frame()
row_names <- c()

for (i in 1:nrow(strategies)){
  .method1 <- strategies[i, 'method']
  .sample_depths1 <- strategies[i, 'sample_depths']
  
  err_row <- c()
  bias_row <- c()
  col_names <- c()
  for (i2 in 1:nrow(strategies)){
    .method2 <- strategies[i2, 'method']
    .sample_depths2 <- strategies[i2, 'sample_depths']
  
  studies1 <- filtered_res%>% filter(method == .method1 & sample_depths == .sample_depths1) %>%
    pull(study_name) %>% unique()
  studies2 <-  filtered_res%>% filter(method == .method2 & sample_depths == .sample_depths2) %>%
    pull(study_name) %>% unique()
  
  studies <- intersect(studies1, studies2)
  
  if (length(studies) >0 ){
  sum_data <- filtered_res %>% filter(study_name %in% studies & method %in% c(.method1, .method2) 
                                      & sample_depths %in% c(.sample_depths1, .sample_depths2)) %>% 
    group_by(method, sample_depths) %>%
    summarize_res()
  ref_data <- sum_data %>% filter(.sample_depths1 == sample_depths & .method1 == method)
  compare_data <- sum_data %>% filter(sample_depths == .sample_depths2 & method == .method2)
  
  .err_difs <- compare_data$RMSE - ref_data$RMSE
  .bias_difs <- abs(compare_data$bias) - abs(ref_data$bias)
  new_name <- paste(.method2, .sample_depths2, sep = ', ')
  
  col_names <- append(col_names, new_name)
  err_row <- append(err_row, .err_difs)
  bias_row <- append(bias_row, .bias_difs)
  }
  }
  
  names(err_row) <- col_names
  names(bias_row) <- col_names
  error_difs <- bind_rows(error_difs, err_row)
  bias_difs <- bind_rows(bias_difs, bias_row)
  row_names <- append(row_names, paste(.method1, .sample_depths1, sep = ', '))
  }
rownames(error_difs) <- row_names
rownames(bias_difs) <-row_names


ggplot(ag_res %>% filter(sample_depths %in% c('30, 40', '30, 45', '30, 50', '30, 60')), 
       aes(mass.change, dif_from_benchmark, color = as.factor(method))) + 
  geom_violin(alpha = .25) + geom_hline(yintercept  = 0) + 
  facet_wrap(~sample_depths)

ggplot(ag_res %>% filter(sample_depths %in% c('10, 30', '20, 30', '15, 30')), 
       aes(mass.change, dif_from_benchmark, color = as.factor(method))) + 
  geom_point(alpha = .25) + geom_hline(yintercept  = 0) + 
  facet_wrap(~sample_depths)


s<- ag_res %>% filter(sample_depths == '30, 45') %>% pull(study_name) %>% unique()
dat_to_plot <- ag_res %>% filter(
  study_name %in% s & ((sample_depths == '30, 45' & grepl('ESM', method)) | 
  (sample_depths == '30' & 
     method %in% c('Mass Correction', 'Mass Correction, SSurgo avg', 'Fixed Depth' ))))

ggplot(dat_to_plot, aes(x = mass.change, y = dif_from_benchmark*100, color = as.factor(method))) + geom_point(alpha = .5) + 
  geom_hline(yintercept = 0)


cross_fold_weights_res <- data.frame()
for (study in unique(out_s4t$study_name)){
  training_data <- out_s4t %>% filter(study_name != study)
  sumed_dat <- training_data %>% filter(!grepl("WL", Ref) &  !grepl('WL', ID) & !grepl('Forest', ID) & !grepl('Forest', Ref)) %>%
    filter(sample_depths == 30) %>% group_by(method) %>% 
    
    summarise(bias = weighted.mean(dif_from_benchmark * ifelse(mass.change >0, 1, -1), weight))
  
  single_d_bias <- sumed_dat %>% filter(method == 'Fixed Depth') %>% 
    pull('bias') %>% first()
  
  linear_MC_bias <- sumed_dat %>% filter(method == 'Linear') %>% 
    pull('bias') %>% first()
   MC_weight <- abs(single_d_bias)/(abs(single_d_bias) + abs(linear_MC_bias))
   sd_weight <- 1-MC_weight
   testing_data <- out_s4t %>% filter(study == study_name) %>% 
     filter(!grepl("WL", Ref) &  !grepl('WL', ID) & !grepl('Forest', ID) & !grepl('Forest', Ref))
   testing_data %>% filter(method %in% c('Single Depth', 'Mass Correction') & sample_depths == '30') %>% 
     pivot_wider(id_cols = c(ID, Ref), names_from = method, values_from = soc_change) %>% 
     mutate(soc_change = `Mass Correction` * MC_weight + `Single Depth` * sd_weight ) %>% 
     merge(testing_data %>% select(ID, Ref, soc_change_benchmark, mass.change, weight) %>% distinct(ID, Ref, .keep_all = T)) %>%
     mutate(dif_from_benchmark = soc_change - soc_change_benchmark, study = study_name,
            method = 'Mass Correction, empirical weights') %>%
     rbind(cross_fold_weights_res ) -> cross_fold_weights_res
   
}


vanDoren_res <- time_series_comparison(vanDoren_fd,
                                       study_name = 'van_doren_1986')

venterea_res <- time_series_comparison(venterea_fd,  
                 study_name = 'Venterea_2006.csv')





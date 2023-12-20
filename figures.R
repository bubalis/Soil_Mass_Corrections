library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(cowplot)

source('analysis_functions.R')


#' Get the legend of a ggplot plot
get.legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

add_sites <- function(grp, .name){
  .name <- first(.name) %>% pull()
  sub_sites <- filter(site_metadata, study == .name) %>%
    pull(site)

  grp$site <- NA
  for (.site in sub_sites){
    if (length(.site) >0){
      print(.site)

      grp <- grp %>%
        mutate(site = ifelse(grepl(.site, grp$ID), .site, site))
    }
  }
  print(grp$site %>% unique)
  return (grp)

}

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

first_of_str <- function(str, split_char = '_'){
  str_split(str, split_char)[[1]][1]
}

out_s4t <- read.csv(file.path('Results', 'all_data_space_for_time_comparison.csv'))

out_s4t <- out_s4t %>%
  group_by(study_name) %>%
  group_modify(add_sites)

out_potash <- read.csv(file.path('Results', 'potash_data_comparisons.csv'))
out_potash$study_name <- str_split(out_potash$ID, '_') %>% sapply(function(x){x[1]})
sub_potash <- filter(potash_data, ID %in% unique(out_potash$ID) | ID %in% unique(out_potash$Ref_ID)  )

get_fd_data <- function(id){filter(sub_potash, study_name == id)}


out_potash <- bind_rows(lapply(unique(out_potash$study_name), function(x){
  add_reference_changes(x, out_potash, load_fd_data_func = get_fd_data, nsims = 50)})) %>%
  mutate(site = study_name)

out_s4t <- rbind_memoryerror(lapply(unique(out_s4t$study_name), function(x){
    add_reference_changes(x, out_s4t, nsims = 50 )}))

ag_res <- out_s4t %>%
  filter(!grepl('WL', ID) & ID!= "Forest Succession" & !grepl('WL', Ref_ID) &
           Ref_ID != "Forest Succession") %>% mutate(site = ifelse(nchar(site) == 0,
                                                     study_name, site))







get_IDs_for_analysis <- function(res_data, depths1, depths2){
  IDs1 <- res_data %>% filter(sample_depths == depths1) %>% pull(ID) %>% unique()
  IDs2 <- res_data %>% filter(sample_depths == depths2) %>% pull(ID) %>% unique()

  intersect(IDs1, IDs2)}




summary_potash <- out_potash %>%
  group_by(method, sample_depths, site, simulation.N) %>%
  mutate(weight = 1) %>%
  summarize_res() %>%
  mutate(data_type = 'Single Core Data')


summary_fieldtrials <- ag_res %>%
  group_by(site, method, sample_depths, simulation.N) %>%
  summarize_res() %>%
  mutate(data_type = 'Field-Level Aggregated Data')


onedepth_sum <-  rbind(summary_potash, summary_fieldtrials) %>%
  filter(sample_depths == '30') %>%
  mutate(method = factor(method)) %>%
  filter(method != 'ESM Spline') %>%
  mutate(method_group = case_when(
    method %in% c('Fixed Depth', 'Linear') ~ "Raw Data Methods",
    method %in% c('Linear Average', 'Linear Average, Fowler Weights') ~ 'Averaging Methods',
    grepl('SSurgo', method) ~ 'SSurgo-based methods'
  )) %>%
  mutate(
  method_group = factor(method_group,
                        levels = c('Raw Data Methods', 'Averaging Methods',
                                   'SSurgo-based methods'))) %>%
  mutate(method = factor(method, levels = c("Fixed Depth",   "Linear",
                                             "Linear Average",
                                             "Linear Average, Fowler Weights",
                                             "Mass Correction, SSurgo avg",
                                            "Mass Correction, SSurgo EXP",
                                          "Mass Correction, SSurgo_spline")))


ggplot(onedepth_sum, aes(x = method_group, y= bias * 100, color = method#, fill = method#,
                           #shape = data_type
                           )) +
    geom_boxplot(position = position_dodge(.75) ) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_y_continuous(name = 'Bulk-Density Bias, T C/ha ', breaks = c(-2,0,2,4),
                     label = c(-2,0,2,4), limits = c(-2, 4.5)) + xlab('') +
  ggtitle('Bias of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples, by Site/Study and Simulation')

  ggplot(onedepth_sum, aes(x = method_group, y= bias_rate, color = method, fill = method)) +
    geom_point(position = position_dodge(.5) ) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_y_continuous(name = 'Bias Rate: Tons C / Tons Mineral Soil',
                       breaks = c(-.01, 0, .01, .02),
                       label = c(-.1, 0, .1, .2),
                       limits = c(-.01, .02)) +
    xlab('') +
    ggtitle('Bias, as a proportion of change in mass to 30 cm, \n of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples,by Site/Study and Simulation')


  ggplot(onedepth_sum, aes(x = method_group, y= RMSE * 100, color = method)) +
    geom_boxplot(position = position_dodge(.75) ) +
    scale_y_continuous(name = 'RMSE, T C/ha', breaks = 0:5, labels = 0:5) +
    xlab('') +
    ggtitle('Error of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples, by Site/Study and Simulation')


  ggplot(onedepth_sum, aes(x = method_group, y = MAPE_change * 100, color = method, fill = method)) +
    geom_boxplot(position = position_dodge(.5) ) +
    #geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_y_continuous(trans = 'log', breaks = c(10, 100, 1000, 2000)) +
    xlab('') +
    ggtitle('Bias, as a proportion of change in mass to 30 cm, \n of Different Methods of Estimating Soil Organic Carbon Change \nfrom Single-Depth Samples, by Site/Study')

#' Plot the data, restricting to the IDs that have the sample depths represented in
#' Sample_depths1 and sample_depths2
plot_for_present_in_all <- function(sample_depths1, sample_depths2){

    ids_potash <- get_IDs_for_analysis(out_potash, sample_depths1, sample_depths2)

    ids_other <- get_IDs_for_analysis(out_s4t, sample_depths1, sample_depths2)

  sum2.potash <- out_potash %>%
    filter(ID %in% ids_potash) %>%
    mutate(weight =1) %>%
    group_by(site, method, sample_depths, simulation.N) %>%
    summarize_res()

  sum2.other <- out_s4t %>%
    filter(ID %in% ids_other) %>%
    group_by(site, method, sample_depths, simulation.N) %>%
    summarize_res()


  sum2_to_plot <- bind_rows(sum2.potash, sum2.other)




  p <- sum2_to_plot %>% filter(sample_depths %in% c(sample_depths1, '30', sample_depths2)) %>%
    filter(sample_depths!= '30' | grepl('SSurgo', method) | method %in% c('Fixed Depth',  'Linear')) %>%
    ggplot(aes(x = sample_depths, y = RMSE *100, color = method)) +
    geom_boxplot(position = position_dodge(.75), alpha = .5) +
    #geom_hline(yintercept = 0, linetype = 'dashed')
    ylab('RMSE T C/ha') + xlab('Depths sampled') +
    ggtitle('RMSE by sampling depth and quantification method \n aggregated by site/study and simulation')
  print(p)

  sum2_to_plot %>%
    filter(sample_depths %in% c(sample_depths1, '30', sample_depths2)) %>%
    filter(sample_depths!= '30' | grepl('SSurgo', method) | method %in% c('Fixed Depth',  'Linear')) %>%
    ggplot(aes(x = sample_depths, y = bias *100, color = method)) +
    geom_boxplot(position = position_dodge(.75)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    ylab('Bulk Density bias, T C/ha') + xlab('Depths sampled') +
    ggtitle('Bias by sampling depth and quantification method \n aggregated by site/study and simulation')

  }



plot_for_present_in_all('15, 30', '30, 45')

plot_for_present_in_all('15, 30', '30, 60')


all_out <- bind_rows(out_potash, ag_res)

standard_summary_stats <- function(grped){
  summarize(grped, bias = mean(dif_from_benchmark), RMSE = sqrt(mean(dif_from_benchmark**2)))
}

NT_comparisons <- ag_res %>% filter((grepl('NT', ID) &
                                      !grepl('NT', Ref_ID) &
                                      !grepl('2000', ID) ) | grepl('No Till', ID)
                                      )

NT.bias_table <- data.frame()
NT.rmse_table <- data.frame()
for (below_depth in c('30, 40', '30, 45', '30, 50', '30, 60') ){
  IDs <- get_IDs_for_analysis(NT_comparisons, '30', below_depth)
  sum.dat <- NT_comparisons %>% filter(ID %in% IDs) %>% group_by(method, sample_depths, study_name) %>%
    standard_summary_stats() %>%

    filter((method %in% c('ESM Spline', 'Two-depth combined sample') & sample_depths == below_depth) |
             method == 'Mass Correction, SSurgo_spline' |method == 'Fixed Depth') %>%
    ungroup %>%
    select(-c(sample_depths))

    NT.bias_table <- sum.dat %>%
      select(-c(RMSE)) %>%
      pivot_wider(names_from = method, values_from = bias) %>%
      mutate(two_depths = below_depth) %>%
      bind_rows(NT.bias_table)

    NT.rmse_table <- sum.dat %>%
            select(-c(bias)) %>%
            pivot_wider(names_from = method, values_from = RMSE) %>%
            mutate(two_depths = below_depth) %>%
            bind_rows(NT.rmse_table)

}

NT_plot_data <- data.frame()
for (below_depth in c('30, 40', '30, 45', '30, 50', '30, 60') ){
  IDs <- get_IDs_for_analysis(NT_comparisons, '30', below_depth)
   NT_plot_data <- NT_comparisons %>% filter(ID %in% IDs) %>%
     filter((method %in% c('ESM Spline', 'Two-depth combined sample') & sample_depths == below_depth) |
              (method %in% c('Fixed Depth', 'Linear', 'Mass Correction, SSurgo avg') &
                 sample_depths == '30')) %>%
    mutate(comp = paste(ID, Ref_ID, sep = '_'),
           compare_depths = below_depth) %>%
    bind_rows(NT_plot_data)



}

NT_plot_data <- mutate(NT_plot_data,
                       sample_type = ifelse(sample_depths == '30',
                                            'Single Depth', 'Two Depth'),
                       method = factor(method,
                                       levels = c('ESM Spline',
                                                  'Two-depth combined sample',
                                                  "Mass Correction, SSurgo avg",
                                                  'Fixed Depth',
                                                  'Linear'
                                                  )))
NT_plot_data <- NT_plot_data %>%
  mutate(x = match(method, levels(NT_plot_data$method))) %>%
  mutate(x = ifelse(x>2, x+1, x))

NT_plt_summed <- NT_plot_data %>%
  group_by(compare_depths, method) %>%
  summarize(.dif_from_benchmark = mean(dif_from_benchmark), x = first(x),
            RMSE = sqrt(mean(dif_from_benchmark**2))
            ) %>%
  rename(dif_from_benchmark = .dif_from_benchmark)

NT_plt_summed <- NT_plt_summed %>%
  rbind(mutate(NT_plt_summed, x = x -.3),
        mutate(NT_plt_summed, x = x + .3))

out_dir <- file.path('Figures')
plot_NT <- ggplot(data = NT_plot_data,
           aes(x = x, y = dif_from_benchmark * 100,
               fill = method, color = method)
           ) +
      geom_boxplot(position =  position_dodge(width = .75), outlier.shape = NA,
                   aes(dodge =comp)) +
      #facet_grid(compare_depths ~ comp + method, scales = "free", space = "free_x") +
      geom_hline(yintercept = 0, linetype = 'dashed') +
      geom_smooth(stat = 'identity', data = NT_plt_summed,
                  aes(y = dif_from_benchmark*100,  x = x, group = method), color = 'black') +
      scale_y_continuous('Error, Tons Carbon per ha', limits = c(-12.5, 10),
                   breaks = c(-10, -5,  0,  5,10)) +
  facet_wrap(~compare_depths) +
  scale_x_continuous('Sample Type',
                     breaks = c(1.5, 5), labels = c('Two Depth', 'One Depth')) +


  ggtitle(
'Summarization of Comparisons No-Till vs Tillage')

leg <- cowplot::get_legend(plot_NT)
plot_NT <- plot_NT +
  theme(strip.text.x = element_text(size = 15),
        axis.text.y =element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size =12))

print(plot_NT)
ggsave(file.path(out_dir, 'no_till_comparisons.png'))


grouped_plot_data <- data.frame()
for (below_depth in c('30, 40', '30, 45', '30, 50', '30, 60') ){
  IDs <- get_IDs_for_analysis(all_out, '30', below_depth)
  grouped_plot_data <- all_out %>%
  filter(
            (method %in% c('ESM Spline', 'Two-depth combined sample') & sample_depths == below_depth) |
           (method %in% c('Fixed Depth', 'Linear', 'Mass Correction, SSurgo avg') &
            sample_depths == '30')) %>%
    filter(ID %in% IDs) %>%
    mutate(comp = paste(ID, Ref_ID, sep = '_'),
           compare_depths = below_depth) %>%
    bind_rows(grouped_plot_data)

  }

grouped_plot_data <- mutate(grouped_plot_data,
                            sample_type = ifelse(sample_depths == '30',
                                                 'Single Depth', 'Two Depth'),
                            method = factor(method,
                                            levels = c('ESM Spline',
                                                       'Two-depth combined sample',
                                                       "Mass Correction, SSurgo avg",
                                                       'Fixed Depth',
                                                       'Linear'
                                            )))
grouped_plot_data <- grouped_plot_data %>%
  mutate(x = match(method, levels(grouped_plot_data$method))) %>%
  mutate(x = ifelse(x>2, x+1, x))

Grp_plt_BD_summed <- grouped_plot_data %>%
  group_by(compare_depths, method) %>%
  summarize(.dif_from_benchmark = mean(dif_from_benchmark * ifelse( mass.change > 0, 1,-1)),
            x = first(x), RMSE = sqrt(mean(dif_from_benchmark**2))) %>%
  rename(dif_from_benchmark = .dif_from_benchmark)

Grp_plt_BD_summed <- Grp_plt_BD_summed %>%
  rbind(mutate(Grp_plt_BD_summed, x = x - .3),
        mutate(Grp_plt_BD_summed, x = x + .3))

plot_All_BD_bias <- ggplot(data = grouped_plot_data,
       aes(x = x, y = dif_from_benchmark * 100 * ifelse( mass.change > 0, 1,-1),
           fill = method, color = method
       )) +
  geom_boxplot(position =  position_dodge(width = .8),  alpha = .5,
               outlier.shape = NA,
               aes(dodge = factor(site))) +
  #facet_grid(compare_depths ~ comp + method, scales = "free", space = "free_x") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_smooth(stat = 'identity', data = Grp_plt_BD_summed,
              aes(y = dif_from_benchmark*100,  x = x, group = method), color = 'black') +
  scale_y_continuous('Error in favor of increased BD, Tons C/ha', limits = c(-5, 10),
                     breaks = c(-5, 0,  5, 10)) +
  facet_wrap(~compare_depths) +
  scale_x_continuous('Sample Type',
                     breaks = c(1.5, 5), labels = c('Two Depth', 'One Depth')) +
  ggtitle(
    'All comparisons, grouped by site, error in favor of > bulk density') +
  theme(strip.text.x = element_text(size = 15),
        axis.text.y =element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size =12))

print(plot_All_BD_bias)
ggsave(file.path(out_dir, 'all_errors_by_site_BD.png'))





Grp_plt_summed <- grouped_plot_data %>%
  group_by(compare_depths, method) %>%
  summarize(.dif_from_benchmark = mean(dif_from_benchmark),
            x = first(x), RMSE = sqrt(mean(dif_from_benchmark**2))) %>%
  rename(dif_from_benchmark = .dif_from_benchmark)

Grp_plt_summed <- Grp_plt_summed %>%
  rbind(mutate(Grp_plt_summed, x = x - .3),
        mutate(Grp_plt_summed, x = x + .3))


plot_All <- ggplot(data = grouped_plot_data,
                           aes(x = x, y = dif_from_benchmark * 100,
                               fill = method, color = method
                           )) +
  geom_boxplot(position =  position_dodge(width = .8),  alpha = .5,
               outlier.shape = NA,
               aes(dodge = factor(site))) +
  #facet_grid(compare_depths ~ comp + method, scales = "free", space = "free_x") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_smooth(stat = 'identity', data = Grp_plt_summed,
              aes(y = dif_from_benchmark*100,  x = x, group = method), color = 'black') +
  scale_y_continuous('ESM Error, Tons C/ha', limits = c(-5, 10),
                     breaks = c(-5, 0,  5, 10)) +
  facet_wrap(~compare_depths) +
  scale_x_continuous('Sample Type',
                     breaks = c(1.5, 5), labels = c('Two Depth', 'One Depth')) +
  ggtitle(
    'Bias') +
  theme(strip.text.x = element_text(size = 15),
        axis.text.y =element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size =12))

print(plot_All)
ggsave(file.path(out_dir, 'all_errors_by_site.png'))


p <- plot_grid(plot_NT +guides(color = F, fill = F) + ggtitle('(A) No-Till vs Tillage'),
               plot_All_BD_bias +guides(color = F, fill = F) + ggtitle('(B) All comparisons, error towards greater bulk density') )
print(p)
ggsave(file.path('Figures', 'combined_for_poster.png'), p, width = 10, height = 5, units = 'in')

p2 <- plot_grid(leg)
ggsave(file.path('Figures', 'comb_legend.png'), p2, width =2, height =2, units = 'in')



#
#Beyond here is mostly nonsense.
#
#

all_out[is.na(all_out$weight), 'weight'] <- 1

bias_table <- data.frame()
rmse_table <- data.frame()
for (below_depth in c('30, 40', '30, 45', '30, 50', '30, 60') ){
  IDs <- get_IDs_for_analysis(all_out, '30', below_depth)
  sum.dat <- all_out %>% filter(ID %in% IDs) %>% group_by(method, sample_depths, study_name) %>%
    summarize_res %>%
    filter((method %in% c('ESM Spline', 'Two-depth combined sample') & sample_depths == below_depth) |
             method == "Mass Correction, SSurgo avg"    |method == 'Fixed Depth') %>% ungroup

  bias_table <- sum.dat %>%
    select(method, bias, study_name) %>%
    pivot_wider(names_from = method, values_from = bias) %>%
    mutate(two_depths = below_depth) %>%
    bind_rows(bias_table)

  rmse_table <- sum.dat %>%
    select(method, RMSE, study_name) %>%
    pivot_wider(names_from = method, values_from = RMSE) %>%
    mutate(two_depths = below_depth) %>%
    bind_rows(rmse_table)

}


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


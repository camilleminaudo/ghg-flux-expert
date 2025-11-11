
# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2023"
# https://github.com/camilleminaudo/restore4cs-scripts
# ---

# --- Description of this script
# This script


# my questions:
# 1- how different are CO2 fluxes estimates between expert timeseries selection and "blind" automation
# 2- how different are CH4 fluxes estimates between expert diffusion/ebullition separation and "blind" automation
# 3- how do experts usually select timeseries?




rm(list = ls()) # clear workspace
cat("/014") # clear console

load_preprocessed <- T

# ---- Directories ----

# You have to make sure this is pointing to the write folder on your local machine
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 

datapath <- paste0(dropbox_root,"/GHG/RAW data")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/RESTORE4Cs/data/Harmonized_GHG")
# plots_path <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/meetings_presentations/2024_PPNW_Girona"
plots_path <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/RESTORE4Cs/GHG_expert_vs_automated/plots/"
results_path <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/RESTORE4Cs/GHG_expert_vs_automated/results"


# ---- packages ----
library(tidyverse)
library(lubridate)
library(zoo)
library(ggplot2)
library(egg)
library(goFlux)
require(dplyr)
require(purrr)
require(msm)
require(data.table)
require(tools)
require(pbapply)

library(ggExtra)
library(ggpubr)
library(scales)


# ---- Functions ----


files.sources = list.files(path = "C:/Projects/myGit/aquaGHG/R/", full.names = T)
for (f in files.sources){source(f)}

get_df_exceedance <- function(relative_diff){
  relative_diff <- relative_diff[!is.na(relative_diff)]
  df_exceedance <- NULL
  for(thresh in sort(unique(relative_diff))){
    n_less <- length(which(abs(relative_diff) < thresh))
    df_exceedance <- rbind(df_exceedance,
                           data.frame(t = thresh,
                                      n = n_less,
                                      p = n_less/length(relative_diff)*100))
  }
  return(df_exceedance)
}

load_this <- function(mylist){
  # load corresponding files
  mydata_all <- NULL
  for(id in mylist){
    message("Loading data for ",id)
    
    myauxfile <- auxfile[which(auxfile$UniqueID==id),]
    gas <- unique(myauxfile$gas_analiser)
    
    setwd(RData_path)
    if(gas== "LI-COR"){
      gs_suffix <- "LI-7810"
    } else {
      gs_suffix <- gas
    }
    load(file = paste0(myauxfile$subsite,"_",gs_suffix,".RData"))
    mydata <- mydata[,c("POSIX.time", 
                        "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                        "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
    
    mydata <- mydata[which(mydata$POSIX.time>=myauxfile$start.time & 
                             mydata$POSIX.time<=myauxfile$start.time+myauxfile$duration),]
    if(dim(mydata)[1]>0){
      mydata$UniqueID <- myauxfile$UniqueID
      mydata$Etime <- as.numeric(mydata$POSIX.time) - min(as.numeric(mydata$POSIX.time))
      
      # how many times assessed?
      ind_which <- which(table_draws$UniqueID==id)
      mydata$n_assessed <- length(ind_which)
      # how many times flagged as anomalous?
      mydata$n_flagged_co2 <- length(which(table_draws$duration_expert_co2[ind_which]<10))
      mydata$n_flagged_ch4 <- length(which(table_draws$duration_expert_ch4[ind_which]<10))
      # defining flag score
      mydata$flag_score_char_co2 <- paste0(mydata$n_flagged_co2,"/",mydata$n_assessed)
      mydata$flag_score_num_co2 <- round(mydata$n_flagged_co2/mydata$n_assessed*100)/100
      mydata$flag_score_char_ch4 <- paste0(mydata$n_flagged_ch4,"/",mydata$n_assessed)
      mydata$flag_score_num_ch4 <- round(mydata$n_flagged_ch4/mydata$n_assessed*100)/100
      
      
      mydata_all <- rbind(mydata_all, mydata)
    }
    
    rm(mydata)
  }
  
  return(mydata_all)
}

# a function to identify outliers
getOutliers_n <- function(y){
  
  # outlier detection based on inter-quartile range
  Q_0.25 <- as.numeric(quantile(y, probs = c(0.25), na.rm = T))
  Q_0.75 <- as.numeric(quantile(y, probs = c(0.75), na.rm = T))
  IQR <- Q_0.75-Q_0.25
  
  delta <- 3
  flag_outlier_low <- y < Q_0.25 - delta * IQR
  flag_outlier_upp <- y > Q_0.75 + delta * IQR
  
  n <- sum(flag_outlier_low) + sum(flag_outlier_upp)
  return(n)
}

standardize.it <- function(x, xmin, xmax){
  x_stand <- (x-xmin)/xmax
}

# a function to get statistics for a group of incubations
get_stats_incubs <- function(mylist){
  # load corresponding files
  mystats_all <- NULL
  for(id in mylist){
    message("Loading data for ",id)
    
    myauxfile <- auxfile[which(auxfile$UniqueID==id),]
    gas <- unique(myauxfile$gas_analiser)
    
    setwd(RData_path)
    if(gas== "LI-COR"){
      gs_suffix <- "LI-7810"
    } else {
      gs_suffix <- gas
    }
    load(file = paste0(myauxfile$subsite,"_",gs_suffix,".RData"))
    mydata <- mydata[,c("POSIX.time", 
                        "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                        "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
    
    mydata <- mydata[which(mydata$POSIX.time>=myauxfile$start.time & 
                             mydata$POSIX.time<=myauxfile$start.time+myauxfile$duration),]
    mydata$Etime <- as.numeric(mydata$POSIX.time) - min(as.numeric(mydata$POSIX.time))
    if(dim(mydata)[1]>0){
      d_dCdt_co2 <- get_dCdt_density(dataframe = mydata, gastype = "CO2dry_ppm")
      d_dCdt_ch4 <- get_dCdt_density(dataframe = mydata, gastype = "CH4dry_ppb")
      
      # detrended_co2 <- data.frame(time = d_dCdt_co2[[2]]$time,
      #                             conc = d_dCdt_co2[[2]]$conc,
      #                             trend = d_dCdt_co2[[2]]$concsmooth,
      #                             flag = abs(d_dCdt_co2[[2]]$conc - d_dCdt_co2[[2]]$concsmooth)/d_dCdt_co2[[2]]$concsmooth > 0.5)
      # 
      # # ggplot(detrended_co2, aes(time, conc, colour = flag))+geom_point()+
      # #   theme_article()
      # 
      # n_flag_CO2 <- sum(detrended_co2$flag)
      
      
      t.win <- 30
      C0_co2 = min(d_dCdt_co2[[2]]$concsmooth[d_dCdt_co2[[2]]$time<t.win])
      Cf_co2 = max(d_dCdt_co2[[2]]$concsmooth[d_dCdt_co2[[2]]$time>max(d_dCdt_co2[[2]]$time)-t.win])
      delta_co2 = Cf_co2-C0_co2
      
      C0_ch4 = min(d_dCdt_ch4[[2]]$concsmooth[d_dCdt_ch4[[2]]$time<t.win])
      Cf_ch4 = max(d_dCdt_ch4[[2]]$concsmooth[d_dCdt_ch4[[2]]$time>max(d_dCdt_ch4[[2]]$time)-t.win])
      delta_ch4 = Cf_ch4-C0_ch4
      
      bubbles <- find_bubbles(time = mydata$Etime, conc = mydata$CH4dry_ppb, window.size = 30)
      if(is.null(bubbles)){n_bubbles = 0} else {n_bubbles <- dim(bubbles)[1]}
      
      mystats <- data.frame(UniqueID = myauxfile$UniqueID,
                            mean_co2 = mean(mydata$CO2dry_ppm),
                            var_co2 = var(mydata$CO2dry_ppm),
                            sd_co2 = sd(mydata$CO2dry_ppm),
                            CV_co2 = sd(mydata$CO2dry_ppm)/mean(mydata$CO2dry_ppm),
                            n_over_3IQR_co2 = getOutliers_n(y = mydata$CO2dry_ppm),
                            delta_co2 = delta_co2,
                            var_dCdt_co2 = var(d_dCdt_co2[[2]]$dydt),
                            p90_dCdt_co2 = as.numeric(quantile(d_dCdt_co2[[2]]$dydt, 0.9)),
                            # n_flag_co2 = n_flag_CO2,
                            
                            mean_ch4 = mean(mydata$CH4dry_ppb),
                            var_ch4 = var(mydata$CH4dry_ppb),
                            sd_ch4 = sd(mydata$CH4dry_ppb),
                            CV_ch4 = sd(mydata$CH4dry_ppb)/mean(mydata$CH4dry_ppb),
                            n_over_3IQR_ch4 = getOutliers_n(y = mydata$CH4dry_ppb),
                            delta_ch4 = delta_ch4,
                            var_dCdt_ch4 = var(d_dCdt_ch4[[2]]$dydt),
                            p90_dCdt_ch4 = as.numeric(quantile(d_dCdt_ch4[[2]]$dydt, 0.9)),
                            n_bubbles = n_bubbles)
      
      mystats_all <- rbind(mystats_all, mystats)
    }
    rm(mystats)
  }
  return(mystats_all)
}



plot_incubations_slcts <- function(list_ids, var){
  mydata_sel <- load_this(mylist = list_ids)
  table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
  if(var == "CO2dry_ppm"){
    table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_co2-(table_draws_sel$start.time_auto)
    table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_co2-(table_draws_sel$start.time_auto)
  } else {
    table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
    table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto) 
  }
  
  n_users <- length(unique(table_draws_sel$userID))
  n_assess <- length(unique(table_draws_sel$uniqAssessID))
  
  table_draws_sel$ymins <- NA
  # for (i in seq(1, length(table_draws_sel$uniqAssessID))){
  #   table_draws_sel$ymins[i] <- mean(mydata_sel[[var]][which(mydata_sel$Etime > table_draws_sel$Etime_start[i] & mydata_sel$Etime < table_draws_sel$Etime_stop[i])], na.rm = T)
  # }
  # delta <- (max(mydata_sel[[var]]) - min(mydata_sel[[var]]))/n_assess/4
  
  table_draws_sel$ymins <- table_draws_sel$delta <- NA
  
  for(id in unique(table_draws_sel$UniqueID)){
    
    n_assess_id <- length(unique(table_draws_sel$uniqAssessID[table_draws_sel$UniqueID==id]))
    
    range_y <- range(mydata_sel[[var]][mydata_sel$UniqueID==id])
    ymins_vect <- seq(range_y[1],range_y[2], length.out = length(unique(table_draws_sel$uniqAssessID[table_draws_sel$UniqueID==id])))
    
    table_draws_sel$ymins[table_draws_sel$UniqueID==id] <- ymins_vect
    
    table_draws_sel$delta <- 0.5*(range_y[2]-range_y[1])/n_assess_id
    
  }
  
  
  
  table_draws_sel$label <- paste0("ID = ",table_draws_sel$draw_corrected)
  mydata_sel$label <- table_draws_sel$label[match(mydata_sel$UniqueID, table_draws_sel$UniqueID)]
  
  plt_bottom <- ggplot()+
    geom_vline(data = table_draws_sel, aes(xintercept=Etime_start, group = UniqueID, colour = userID), alpha=0.3, size=1)+
    geom_vline(data = table_draws_sel, aes(xintercept=Etime_stop, group = UniqueID, colour = userID), alpha=0.3, size=1)+
    geom_segment(data = table_draws_sel, aes(x = Etime_start, xend= Etime_stop, y = userID, colour = userID), alpha=0.9, size=4)+
    theme_article()+
    xlab("Elapsed time [secs]")+
    ylab("")+
    scale_color_brewer(palette = "Dark2")
  if (length(unique(table_draws_sel$UniqueID))>1){
    plt_bottom <- plt_bottom+
      facet_wrap(label~., scales = "free_x")
  }
  
  
  
  if(var == "CO2dry_ppm"){
    
    plt_top <- ggplot()+
      geom_vline(data = table_draws_sel, aes(xintercept=Etime_start, group = UniqueID, colour = userID), alpha=0.3, size=1)+
      geom_vline(data = table_draws_sel, aes(xintercept=Etime_stop, group = UniqueID, colour = userID), alpha=0.3, size=1)+
      # geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
      #                                       ymin = ymins-delta, ymax = ymins+delta, fill=userID), alpha=0.2)+
      geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
      theme_article()+
      xlab("Elapsed time [secs]")+
      ylab(expression(""*CO[2~~dry]~" [ppm]"))+
      scale_color_brewer(palette = "Dark2")
    
    if (length(unique(table_draws_sel$UniqueID))>1){
      plt_top <- plt_top +
        facet_wrap(label~., scales = "free_x")
    }
    
    
    
  } else {
    
    plt_top <- ggplot()+
      geom_vline(data = table_draws_sel, aes(xintercept=Etime_start, group = UniqueID, colour = userID), alpha=0.3, size=1)+
      geom_vline(data = table_draws_sel, aes(xintercept=Etime_stop, group = UniqueID, colour = userID), alpha=0.3, size=1)+
      # geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
      #                                       ymin = ymins-delta, ymax = ymins+delta, fill=userID), alpha=0.2)+
      geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
      theme_article()+
      xlab("Elapsed time [secs]")+
      ylab(expression(""*CH[4~~dry]~" [ppb]"))+
      scale_color_brewer(palette = "Dark2")
    
    if (length(unique(table_draws_sel$UniqueID))>1){
      plt_top <- plt_top +
        facet_wrap(label~., scales = "free_x")
    }
    
  }
  
  myplt <- ggarrange(plt_top, plt_bottom, align = "hv", ncol = 1, heights = c(0.6,0.4), legend = "none", labels = c("a","b"))
 
  return(myplt)
}

# ---- Loading data ----

setwd(dirname(results_path))

if (load_preprocessed){
  load("results/results_ghg_experts_preprocessed.Rdata")
} else {
  
  
  # loading auxfile
  auxfile <- read.csv(file = "auxfile.csv")
  auxfile$start.time <- as.POSIXct(auxfile$start.time, tz = 'UTC') 
  
  load("results/results_ghg_experts_recalculated.Rdata")
  
  
  
  # load_fs <- function(path, pattern){
  #   
  #   list_f <- list.files(path = path, pattern = pattern, full.names = T)
  #   list_f <- list_f[grep(pattern = "csv", x = list_f)]
  #   
  #   isF <- T
  #   for (f in list_f){
  #     df_tmp <- read.csv(f)
  #     if(isF){
  #       isF <- F
  #       df_out <- df_tmp
  #     } else {
  #       df_out <- rbind(df_out,df_tmp)
  #     }
  #   }
  #   return(df_out)
  # }
  # 
  # table_all <- load_fs(path = results_path, pattern = "BLIND_vs_EXPERT_co2_ch4_fluxes_")
  # table_ebull <- load_fs(path = results_path, pattern = "ch4_ebullition")
  # table_ebull$ebullition[which(table_ebull$ebullition>table_ebull$total_estimated)] <- table_ebull$total_estimated[which(table_ebull$ebullition>table_ebull$total_estimated)]
  # table_ebull$ebullition[which(table_ebull$ebullition < 0)] <- 0
  # table_draws <- load_fs(path = results_path, pattern = "table_draw")
  # 
  # 
  # 
  # # clean names
  # table_all <- table_all %>%
  #   mutate(username = str_replace_all(username, fixed("Camille Minaudo"), "Camille"),
  #          username = str_replace_all(username, fixed("Dmorant"), "DMorant"),
  #          username = str_replace_all(username, fixed("Jorge Montes"), "Jorge"),
  #          username = str_replace_all(username, fixed("Katrin Attermeyer"), "Katrin"),
  #          username = str_replace_all(username, fixed("Miguel Cabrera"), "Miguel"))
  # 
  # table_ebull <- table_ebull %>%
  #   mutate(username = str_replace_all(username, fixed("Camille Minaudo"), "Camille"),
  #          username = str_replace_all(username, fixed("Dmorant"), "DMorant"),
  #          username = str_replace_all(username, fixed("Jorge Montes"), "Jorge"),
  #          username = str_replace_all(username, fixed("Katrin Attermeyer"), "Katrin"),
  #          username = str_replace_all(username, fixed("Miguel Cabrera"), "Miguel"))
  # 
  # table_draws <- table_draws %>%
  #   mutate(username = str_replace_all(username, fixed("Camille Minaudo"), "Camille"),
  #          username = str_replace_all(username, fixed("Dmorant"), "DMorant"),
  #          username = str_replace_all(username, fixed("Jorge Montes"), "Jorge"),
  #          username = str_replace_all(username, fixed("Katrin Attermeyer"), "Katrin"),
  #          username = str_replace_all(username, fixed("Miguel Cabrera"), "Miguel"))
  
  
  # statistics on all incubations
  stats_all <- get_stats_incubs(mylist = unique(table_draws$UniqueID))
  names(stats_all)
  
  stats_all$flux_co2 <- CO2_flux.auto$best.flux[match(stats_all$UniqueID, CO2_flux.auto$UniqueID)]
  stats_all$flux_ch4 <- CH4_flux.auto$best.flux[match(stats_all$UniqueID, CH4_flux.auto$UniqueID)]
  stats_all$ebull_ch4 <- CH4_flux.auto$ebullition.flux[match(stats_all$UniqueID, CH4_flux.auto$UniqueID)]
  stats_all$diffus_ch4 <- CH4_diff_flux.man$best.flux[match(stats_all$UniqueID, CH4_diff_flux.man$UniqueID)]
  
  
  # ---- Retrieve timestamp processing for table_draw and making unique IDs ----
  
  # this was done in recalculate_fluxes_with_aquaGHG.R
  
  
  save.image(file = "results/results_ghg_experts_preprocessed.Rdata")
  
}



# ---- overview database ----

message("overview of in-situ database")

p_duration <- ggplot(auxfile, aes(duration))+geom_density(alpha=0.5)+
  xlab("Duration [secs]")+
  # ylab("Water depth [cm]")+
  theme_article()

p_depth <- ggplot(auxfile, aes(water_depth))+geom_density(alpha=0.5)+
  xlab("Water depth [cm]")+
  theme_article()


p_mean_sd_co2 <- ggplot(stats_all, aes(mean_co2, sd_co2))+geom_point(alpha=0.5)+
  xlab("Average CO2 [ppm]")+
  ylab("Standard deviation CO2 [ppm]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()

p_mean_sd_ch4 <- ggplot(stats_all, aes(mean_ch4, sd_ch4))+geom_point(alpha=0.5)+
  xlab("Average CH4 [ppb]")+
  ylab("Standard deviation CH4 [ppb]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()

p_dist_co2 <- ggMarginal(p_mean_sd_co2, type="density")
p_dist_ch4 <- ggMarginal(p_mean_sd_ch4, type="density")

p_overview_db <- ggpubr::ggarrange(p_duration, p_depth, p_dist_co2, p_dist_ch4, ncol = 2, nrow = 2,labels = c("a","b","c","d"))
p_overview_db

# ggsave(plot = p_overview_db, filename = "Fig1_overview_database_ghg_expert.jpeg", path = plots_path, 
#        width = 6, height = 6, dpi = 300, units = 'in', scale = 1.0)




message("contributions by experts")
table(table_draws$username)

tab_users <- as.data.frame(table(table_draws$username))
tab_users$userID <- paste0("user",seq(1,dim(tab_users)[1]))

tab_users$userID <- factor(x = tab_users$userID, levels = tab_users$userID[order(tab_users$Freq, decreasing = T)])

table_draws$userID <- tab_users$userID[match(table_draws$username, tab_users$Var1)]

p_users <- ggplot(tab_users, aes(userID, Freq))+geom_point()+theme_article()+
  # ggtitle(paste0("n (users) = ",dim(tab_users)[1]))+
  xlab("User ID")+
  ylab("Number of incubations")+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

p_users



n <- length(unique(table_draws$UniqueID))
n_inspections <- dim(table_draws)[1]
table_n <- data.frame(id = names(table(table_draws$UniqueID)),
                      n = as.numeric(table(table_draws$UniqueID)))


p_overview <- ggplot(table_n, aes(n))+geom_histogram()+theme_article()+xlab("Repetitions")+ylab("Number of incubations")+
  # ggtitle(paste0("n (unique incubations) = ",n, ", n (inspections) = ",n_inspections))+
  scale_x_continuous(breaks = breaks_pretty())
p_overview



table_n <- data.frame(id = names(table(table_draws$UniqueID)),
                      n = as.numeric(table(table_draws$UniqueID)),
                      key = "by_incubation")


table_n_userID <- data.frame(id = names(table(paste(table_draws$draw,table_draws$userID,sep = "-"))),
                             n = as.numeric(table(paste(table_draws$draw,table_draws$userID,sep = "-"))),
                             key = "by_user")


table_n_both <- rbind(table_n, table_n_userID)


p_overview_userID <- ggplot(table_n_both, aes(n, fill = key))+
  geom_histogram(position = "jitter")+theme_article()+xlab("Repetitions")+ylab("Number of incubations")+
  scale_x_continuous(breaks = breaks_pretty())
p_overview_userID



length(which(table_n_userID$n>1))

table_multiple_assess_by_user <- table_n_userID[which(table_n_userID$n==2),]

find_dash <- table_multiple_assess_by_user$id %>% str_locate_all("-") %>% sapply(., function(x) x[,2])

length(unique(str_sub(table_multiple_assess_by_user$id, start = 1, end = find_dash-1)))


p_overview_contribs <- ggpubr::ggarrange(p_users, p_overview, ncol = 2, nrow = 1, align = "h",labels = c("a","b"))
p_overview_contribs
ggsave(plot = p_overview_contribs, filename = "Fig1_overview_contributions_ghg_expert.jpeg", path = plots_path,
       width = 6, height = 3, dpi = 300, units = 'in', scale = 1.2)


# example of experts' selection
table_draws$duration_ratio_co2 <- table_draws$duration_expert_co2/(table_draws$end.time_auto - table_draws$start.time_auto) 
table_draws$duration_ratio_ch4 <- table_draws$duration_expert_ch4/(table_draws$end.time_auto - table_draws$start.time_auto) 

# table(table_draws$UniqueID[table_draws$duration_ratio_co2<.5])
# list_ids <- names(table(table_draws$UniqueID[table_draws$duration_ratio_co2<.5])[which(table(table_draws$UniqueID[table_draws$duration_ratio_co2<.5])>5)])


list_ids <- names(table(table_draws$UniqueID[table_draws$duration_ratio_co2>.8])[which(table(table_draws$UniqueID[table_draws$duration_ratio_co2>.8])>5)])



# selected incubations for draft
# list_ids <- c("s1-cu-a2-16-o-d-11:58","s3-va-a1-4-o-d-10:04")
list_ids <- "s3-va-a1-4-o-d-10:04"
p_exemple_co2 <- plot_incubations_slcts(list_ids = list_ids[6:10], var = "CO2dry_ppm")
p_exemple_co2

ggsave(plot = p_exemple_co2, filename = "Fig2_exemple_CO2_selection_experts.jpeg", path = plots_path, 
       width = 4, height = 5, dpi = 300, units = 'in', scale = 1.0)


p_exemple_ch4 <- plot_incubations_slcts(list_ids = list_ids, var = "CH4dry_ppb")
p_exemple_ch4

ggsave(plot = p_exemple_ch4, filename = "FigS1_exemple_CH4_selection_experts.jpeg", path = plots_path, 
       width = 4, height = 5, dpi = 300, units = 'in', scale = 1.0)






# ---- Flagged anomalous ----

setwd(results_path)
load("table_draws_fixed.RData")

threshold_anomalous <- 10 # seconds or n observations

ggplot(table_draws_all_fixed, aes(duration_expert_co2 ))+
  geom_histogram(alpha=0.5)+theme_article()+geom_vline(xintercept = threshold_anomalous)
ggplot(table_draws_all_fixed, aes(duration_expert_ch4))+
  geom_histogram(alpha=0.5)+theme_article()+geom_vline(xintercept = threshold_anomalous)

# table_all$uniqAssessID <- paste(table_all$timestamp_processing, table_all$UniqueID, sep="/")
# table_ebull$uniqAssessID <- paste(table_ebull$timestamp_processing, table_ebull$UniqueID, sep="/")

table_flagged_co2 <- data.frame(variable = "CO2",
                                uniqAssessID = unique(table_draws_all_fixed$uniqAssessID[which(
                                  table_draws_all_fixed$duration_expert_co2 < threshold_anomalous)]))
table_flagged_co2$UniqueID <- table_draws_all_fixed$UniqueID[match(table_flagged_co2$uniqAssessID, table_draws_all_fixed$uniqAssessID)]

table_flagged_co2$n_assessed <- table_flagged_co2$n_flagged <- NA
for(i in seq(1,length(table_flagged_co2$uniqAssessID))){
  table_flagged_co2$n_assessed[i] <- length(which(table_draws_all_fixed$UniqueID==table_flagged_co2$UniqueID[i]))
  table_flagged_co2$n_flagged[i] <- length(which(table_flagged_co2$UniqueID==table_flagged_co2$UniqueID[i]))
}
table_flagged_co2$flagged_rate <- table_flagged_co2$n_flagged/table_flagged_co2$n_assessed


table_flagged_ch4 <- data.frame(variable = "CH4",
                                uniqAssessID = unique(table_draws_all_fixed$uniqAssessID[which(
                                  table_draws_all_fixed$duration_expert_ch4 < threshold_anomalous)]))
table_flagged_ch4$UniqueID <- table_draws_all_fixed$UniqueID[match(table_flagged_ch4$uniqAssessID, table_draws_all_fixed$uniqAssessID)]

table_flagged_ch4$n_assessed <- table_flagged_ch4$n_flagged <- NA
for(i in seq(1,length(table_flagged_ch4$uniqAssessID))){
  table_flagged_ch4$n_assessed[i] <- length(which(table_draws_all_fixed$UniqueID==table_flagged_ch4$UniqueID[i]))
  table_flagged_ch4$n_flagged[i] <- length(which(table_flagged_ch4$UniqueID==table_flagged_ch4$UniqueID[i]))
}
table_flagged_ch4$flagged_rate <- table_flagged_ch4$n_flagged/table_flagged_ch4$n_assessed


table_flagged <- rbind(table_flagged_co2, table_flagged_ch4)

table_flagged_nassessed_more_than_3 <- table_flagged[which(table_flagged$n_assessed>=3),]

which(table_flagged_nassessed_more_than_3$flagged_rate==1)


# identify incubations flagged for both CO2 and CH4
list_both <- table_flagged$UniqueID[duplicated(paste0(table_flagged$variable, table_flagged$UniqueID))]

tab_flagged_plt <- table_flagged#[!duplicated(table_flagged$UniqueID),]
tab_flagged_plt$variable[tab_flagged_plt$UniqueID %in% list_both] <- "CO2 and CH4"


# ggplot(tab_flagged_plt, 
#        aes(n_assessed, flagged_rate, colour = variable))+geom_jitter()+
#   xlab("n")+ylab("Flagging rate")+
#   theme_article()+
#   ylim(c(0,1))+
#   scale_colour_viridis_d(begin = 0.2, end = .9, option = "A", direction = -1)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
#   scale_x_continuous(breaks = breaks_pretty())


plt_flaggingRate <- ggplot(tab_flagged_plt,
                           aes(as.factor(n_assessed), flagged_rate, colour = variable))+geom_jitter(alpha = 0.9, size=3)+
  xlab("Counts of incubation assessment")+
  ylab("Flagging rate")+
  theme_article()+
  ylim(c(0,1))+
  scale_colour_viridis_d(begin = 0.2, end = .9, option = "A", direction = -1)


ggplot(tab_flagged_plt)+
  # geom_point(aes(reorder(UniqueID, flagged_rate), n_assessed))+
  geom_point(aes(reorder(UniqueID, n_assessed), n_flagged))+
  # coord_flip()+
  # xlab("Incubation ID")+ylab("Flagging rate")+
  theme_article()+
  ylim(c(0,10))+
  scale_colour_viridis_d(begin = 0.2, end = .9, option = "A", direction = -1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))#+facet_grid(.~variable, scales = "free_x")


ggsave(plot = plt_flaggingRate, filename = "flaggingRate.jpeg", path = plots_path,
       width = 6, height = 4, dpi = 300, units = 'in', scale = 1)



tab_flagged_plt[which(tab_flagged_plt$UniqueID=="s4-da-a2-10-o-d-08:29"),]





ggplot(table_flagged[duplicated(table_flagged$UniqueID),], #http://127.0.0.1:18381/graphics/74e43cbf-6549-4826-8127-6bcf615f7086.png
       aes(n_assessed, flagged_rate, colour = variable))+geom_jitter(size=3)+theme_article()+
  xlab("n(assessed)")+ylab("Flagged rate")+
  scale_x_continuous(breaks = breaks_pretty())


tmpstmp_flag_co2 <- unique(table_draws_all_fixed$uniqAssessID[which(table_draws_all_fixed$duration_expert_co2 < threshold_anomalous)])
tmpstmp_flag_ch4 <- unique(table_draws_all_fixed$uniqAssessID[which(table_draws_all_fixed$duration_expert_ch4 < threshold_anomalous)])

incubs_flagged_co2 <- unique(table_draws_all_fixed$UniqueID[table_draws_all_fixed$duration_expert_co2<threshold_anomalous])
incubs_flagged_ch4 <- unique(table_draws_all_fixed$UniqueID[table_draws_all_fixed$duration_expert_ch4<threshold_anomalous])

tmp <- rbind(data.frame(variable = "co2",
                        ID = incubs_flagged_co2),
             data.frame(variable = "ch4",
                        ID = incubs_flagged_ch4))

tmp_flag_both <- tmp[duplicated(tmp$ID), ]
tmp_flag_NOTboth <- tmp[ ! tmp$ID %in% tmp_flag_both$ID, ]

# what is special about these flagged incubations?
stats_all$co2flagged <- "valid"
stats_all$co2flagged[which(stats_all$UniqueID%in%incubs_flagged_co2)] <- "flagged"
stats_all$ch4flagged <- "valid"
stats_all$ch4flagged[which(stats_all$UniqueID%in%incubs_flagged_ch4)] <- "flagged"


p_mean_sd_co2_f <- ggplot(stats_all, aes(mean_co2, sd_co2, colour = co2flagged))+geom_point(alpha=0.5)+
  xlab("Average CO2 [ppm]")+
  ylab("Standard deviation CO2 [ppm]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  scale_colour_viridis_d(begin = 0.2, end = .9, option = "A", direction = 1)+
  theme(legend.title = element_blank())+ggtitle(paste0("n flagged = ", length(incubs_flagged_co2)))

p_mean_sd_ch4_f <- ggplot(stats_all, aes(mean_ch4, sd_ch4, colour = ch4flagged))+geom_point(alpha=0.5)+
  xlab("Average CH4 [ppb]")+
  ylab("Standard deviation CH4 [ppb]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  scale_colour_viridis_d(begin = 0.2, end = .9, option = "A", direction = 1)+
  theme(legend.title = element_blank())+ggtitle(paste0("n flagged = ", length(incubs_flagged_ch4)))

# p_dist_co2 <- ggMarginal(p_mean_sd_co2, type="density", groupFill = T)
p_flagged <- ggpubr::ggarrange(p_users, p_overview, p_mean_sd_co2_f, p_mean_sd_ch4_f, common.legend = T, legend = "right")

ggsave(plot = p_flagged, filename = "flagged.jpeg", path = plots_path,
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.0)




# plot incubations with the highest disagreements between experts
list_ids <- tail(incubs_flagged_co2, 9)

mydata_sel <- load_this(mylist = list_ids)

p_flagged_CO2 <- ggplot()+
  # geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
  #                                       ymin = -Inf, ymax = Inf, fill=userID), alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  ylab("CO2_dry [ppm]")+
  facet_wrap(UniqueID~., scales = "free")
p_flagged_CO2
ggsave(plot = p_flagged_CO2, filename = "examples_flagged_CO2.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.0)


plot_incubations_slcts(list_ids, var = "CO2dry_ppm")



# ---- remove from the tables the incubations flagged as anomalous ----

# first we remove from the tables the incubations flagged as anomalous
table_draws <- table_draws[ ! table_draws$uniqAssessID %in% unique(table_flagged$uniqAssessID), ]




# ---- How do expert trim time series? ----

# ggplot(table_draws, aes(end.time_auto-end.time_expert_ch4))+geom_density()+theme_article()


p_less_than_50perc_co2 <- length(which(table_draws$duration_ratio_co2 < .5))/dim(table_draws)[1]*100
p_less_than_50perc_ch4 <- length(which(table_draws$duration_ratio_ch4 < .5))/dim(table_draws)[1]*100


p_more_than_90perc_co2 <- length(which(table_draws$duration_ratio_co2 > .9))/dim(table_draws)[1]*100
p_more_than_90perc_ch4 <- length(which(table_draws$duration_ratio_ch4 > .9))/dim(table_draws)[1]*100


# how many trimmed time series start after half the original measurements?
length(which((table_draws$start.time_expert_co2 - table_draws$start.time_auto) > ((table_draws$end.time_auto - table_draws$start.time_auto)/2)))/dim(table_draws)[1]*100
length(which((table_draws$start.time_expert_ch4 - table_draws$start.time_auto) > ((table_draws$end.time_auto - table_draws$start.time_auto)/2)))/dim(table_draws)[1]*100


p_density_duration <- ggplot(table_draws)+
  geom_density(aes(duration_expert_co2/(end.time_auto - start.time_auto), fill="CO2"), color=alpha("black", 0.1), alpha=.5)+
  geom_density(aes(duration_expert_ch4/(end.time_auto - start.time_auto), fill="CH4"), color=alpha("black", 0.1), alpha=.5)+
  scale_fill_manual(values = c("#fc8d62", "#8da0cb"), name = "Variable")+
  xlab(expression(duration[~expert~selection]~~"/ " * duration[~original~data]))+
  theme_article()+theme(legend.position = c(0.2,0.8))
p_density_duration




row_blank <- seq(0,1,0.01)
mymatco2 <- mymatch4 <- NULL
table_draws$code_row_co2 <- table_draws$code_row_ch4 <- NA
for(i in seq(1,dim(table_draws)[1])){
  t_start_co2 <- (table_draws$start.time_expert_co2[i] - table_draws$start.time_auto[i])/(table_draws$end.time_auto[i] - table_draws$start.time_auto[i])
  t_stop_co2 <- (table_draws$end.time_expert_co2[i] - table_draws$start.time_auto[i])/(table_draws$end.time_auto[i] - table_draws$start.time_auto[i])
  
  myrow <- 0*row_blank
  myrow[which(row_blank>=t_start_co2 & row_blank<=t_stop_co2)] <- 3
  code_row_co2 <- 3
  if((t_stop_co2 - t_start_co2) < .5){myrow[which(row_blank>=t_start_co2 & row_blank<=t_stop_co2)] <- 2; code_row_co2 <- 2}
  if(t_start_co2 > .5){myrow[which(row_blank>=t_start_co2 & row_blank<=t_stop_co2)] <- 4; code_row_co2 <- 4}
  if(t_stop_co2 < .5){myrow[which(row_blank>=t_start_co2 & row_blank<=t_stop_co2)] <- 1; code_row_co2 <- 1}
  
  mymatco2 <- rbind(mymatco2,myrow)
  table_draws$code_row_co2[i] <- code_row_co2
  
  
  t_start_ch4 <- (table_draws$start.time_expert_ch4[i] - table_draws$start.time_auto[i])/(table_draws$end.time_auto[i] - table_draws$start.time_auto[i])
  t_stop_ch4 <- (table_draws$end.time_expert_ch4[i] - table_draws$start.time_auto[i])/(table_draws$end.time_auto[i] - table_draws$start.time_auto[i])
  
  myrow <- 0*row_blank
  myrow[which(row_blank>=t_start_ch4 & row_blank<=t_stop_ch4)] <- 3
  code_row_ch4 <- 3
  if((t_stop_ch4 - t_start_ch4) < .5){myrow[which(row_blank>=t_start_ch4 & row_blank<=t_stop_ch4)] <- 2; code_row_ch4 <- 2}
  if(t_start_ch4 > .5){myrow[which(row_blank>=t_start_ch4 & row_blank<=t_stop_ch4)] <- 4; code_row_ch4 <- 4}
  if(t_stop_ch4 < .5){myrow[which(row_blank>=t_start_ch4 & row_blank<=t_stop_ch4)] <- 1; code_row_ch4 <- 1}
  
  mymatch4 <- rbind(mymatch4,myrow)
  table_draws$code_row_ch4[i] <- code_row_ch4
}

# image(t(mymatco2[order((table_draws$start.time_expert_co2 - table_draws$start.time_auto)/(table_draws$end.time_auto - table_draws$start.time_auto)),]),
#       col=c("#DCDCDC", "#1b9e77"))

# order by duration:
# mymatco2_ord <- mymatco2[order((table_draws$end.time_expert_co2 - table_draws$start.time_expert_co2)/(table_draws$end.time_auto - table_draws$start.time_auto)),]

# order by starting time
# mymatco2_ord <- mymatco2[order(table_draws$code_row_co2, (table_draws$start.time_expert_co2-table_draws$start.time_auto)/(table_draws$end.time_auto-table_draws$start.time_auto)),]
mymatco2_ord <- mymatco2[order((table_draws$start.time_expert_co2-table_draws$start.time_auto)/(table_draws$end.time_auto-table_draws$start.time_auto)),]

rownames(mymatco2_ord) <- as.character(seq(1,dim(mymatco2)[1]))
longData_co2<-reshape2::melt(mymatco2_ord)
longData_co2$value_str <- "void"
longData_co2$value_str[longData_co2$value==1] <- "1st half"
longData_co2$value_str[longData_co2$value==2] <- "less than half"
longData_co2$value_str[longData_co2$value==3] <- "more than half"
longData_co2$value_str[longData_co2$value==4] <- "2nd half"

longData_co2$value_str <- factor(x = longData_co2$value_str, levels = c("void","1st half","less than half","more than half","2nd half"))

require(scales)
p_mat_co2 <- ggplot(longData_co2[longData_co2$value>0,], aes(x = Var2/100, y = Var1)) +
  geom_raster(aes(fill=as.factor(value_str)))+
  scale_fill_manual(values = c("#e7298a","#fdbb84","#99d8c9","#7570b3"), name = "Experts' selection")+
  xlab("Normalized incubation time")+
  ylab(expression("Running ID of "*CO[2]~" time series"))+
  # ggtitle("Experts' selection for CO2 incubations")+
  theme_article() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11),
                     legend.position = 'bottom')



# order by starting time
# mymatch4_ord <- mymatch4[order(table_draws$code_row_ch4, (table_draws$start.time_expert_ch4-table_draws$start.time_auto)/(table_draws$end.time_auto-table_draws$start.time_auto)),]
mymatch4_ord <- mymatch4[order((table_draws$start.time_expert_ch4-table_draws$start.time_auto)/(table_draws$end.time_auto-table_draws$start.time_auto)),]

rownames(mymatch4_ord) <- as.character(seq(1,dim(mymatch4)[1]))
longData_ch4<-reshape2::melt(mymatch4_ord)
longData_ch4$value_str <- "void"
longData_ch4$value_str[longData_ch4$value==1] <- "1st half"
longData_ch4$value_str[longData_ch4$value==2] <- "less than half"
longData_ch4$value_str[longData_ch4$value==3] <- "more than half"
longData_ch4$value_str[longData_ch4$value==4] <- "2nd half"

longData_ch4$value_str <- factor(x = longData_ch4$value_str, levels = c("void","1st half","less than half","more than half","2nd half"))

p_mat_ch4 <- ggplot(longData_ch4[longData_ch4$value>0,], aes(x = Var2/100, y = Var1)) +
  geom_raster(aes(fill=as.factor(value_str)))+
  scale_fill_manual(values = c("#e7298a","#fdbb84","#99d8c9","#7570b3"), name = "Experts' selection")+
  xlab("Normalized incubation time")+
  ylab(expression("Running ID of "*CH[4]~" time series"))+
  # ggtitle("Experts' selection for ch4 incubations")+
  theme_article() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11),
                          legend.position = 'bottom')



table_draws$code_row_co2_str <- "void"
table_draws$code_row_co2_str[table_draws$code_row_co2==1] <- "1st half"
table_draws$code_row_co2_str[table_draws$code_row_co2==2] <- "less than half"
table_draws$code_row_co2_str[table_draws$code_row_co2==3] <- "more than half"
table_draws$code_row_co2_str[table_draws$code_row_co2==4] <- "2nd half"
table_draws$code_row_co2_str <- factor(x = table_draws$code_row_co2_str, levels = c("void","1st half","less than half","2nd half","more than half"))

table_draws$code_row_co2_str_glob <- "more than half"
table_draws$code_row_co2_str_glob[table_draws$code_row_co2!=3] <- "less than half"

p_repartition_co2 <- ggplot(table_draws, aes(code_row_co2_str_glob, fill = as.factor(code_row_co2_str)))+geom_bar()+
  scale_fill_manual(values = c("#e7298a","#fdbb84","#7570b3","#99d8c9"), name = "Experts' selection")+
  xlab(expression("Selected span for "*CO[2]))+
  ylab("Time series count")+
  theme_article() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11),
                          legend.position = 'none',
                          legend.title = element_blank())




table_draws$code_row_ch4_str <- "void"
table_draws$code_row_ch4_str[table_draws$code_row_ch4==1] <- "1st half"
table_draws$code_row_ch4_str[table_draws$code_row_ch4==2] <- "less than half"
table_draws$code_row_ch4_str[table_draws$code_row_ch4==3] <- "more than half"
table_draws$code_row_ch4_str[table_draws$code_row_ch4==4] <- "2nd half"
table_draws$code_row_ch4_str <- factor(x = table_draws$code_row_ch4_str, levels = c("void","1st half","less than half","2nd half","more than half"))

table_draws$code_row_ch4_str_glob <- "more than half"
table_draws$code_row_ch4_str_glob[table_draws$code_row_ch4!=3] <- "less than half"

p_repartition_ch4 <- ggplot(table_draws, aes(code_row_ch4_str_glob, fill = as.factor(code_row_ch4_str)))+geom_bar()+
  scale_fill_manual(values = c("#e7298a","#fdbb84","#7570b3","#99d8c9"), name = "Experts' selection")+
  xlab(expression("Selected span for "*CH[4]))+
  ylab("Time series count")+
  theme_article() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11),
                          legend.position = 'none',
                          legend.title = element_blank())


p_duration_co2 <- ggplot(table_draws, aes(code_row_co2_str_glob, duration_expert_co2, fill = as.factor(code_row_co2_str)))+geom_boxplot()+
  scale_fill_manual(values = c("#e7298a","#fdbb84","#7570b3","#99d8c9"), name = "experts'selection")+
  xlab(expression("Selected span for "*CO[2]))+
  ylab("Duration of selection [secs]")+
  theme_article() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11),
                          legend.position = 'none',
                          legend.title = element_blank())

p_duration_ch4 <- ggplot(table_draws, aes(code_row_ch4_str_glob, duration_expert_ch4, fill = as.factor(code_row_ch4_str)))+geom_boxplot()+
  scale_fill_manual(values = c("#e7298a","#fdbb84","#7570b3","#99d8c9"), name = "experts'selection")+
  xlab(expression("Selected span for "*CH[4]))+
  ylab("Duration of selection [secs]")+
  theme_article() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=9),
                          plot.title=element_text(size=11),
                          legend.position = 'none',
                          legend.title = element_blank())

p_mat_exprts_selections <- ggpubr::ggarrange(p_mat_co2, p_mat_ch4, ggpubr::ggarrange(p_repartition_co2, p_repartition_ch4,
                  p_duration_co2, p_duration_ch4,
                  labels = c("c", "d","e","f"), align = "hv", common.legend = F), nrow = 1, widths = c(0.2,0.2,0.6),
                  labels = c("a", "b"), align = "hv", common.legend = T, legend = "bottom")



p_mat_exprts_selections <- ggpubr::ggarrange(p_mat_co2, p_mat_ch4, p_repartition_co2, p_repartition_ch4, nrow = 1, widths = c(0.3,0.3,0.2,0.2),
                                             labels = c("a", "b","c","d"), align = "hv", common.legend = T, legend = "bottom")


p_mat_exprts_selections


ggsave(plot = p_mat_exprts_selections, filename = "Fig2_overview_time_series_selection_experts.jpeg", path = plots_path, 
       width = 10, height = 4., dpi = 300, units = 'in', scale = 1.3)



# p_co2_ch4_fluxes <- ggplot(stats_all, aes(flux_co2, flux_ch4))+geom_point(alpha=0.5)+
#   # xlab("Average CO2 [ppm]")+
#   # ylab("Standard deviation CO2 [ppm]")+
#   scale_x_log10()+
#   scale_y_log10()+
#   theme_article()
# ggMarginal(p_co2_ch4_fluxes, type="density")
# 
# 
# ggplot(stats_all, aes(as.factor(n_bubbles), ebull_ch4))+geom_boxplot(alpha=0.5)+
#   # xlab("Average CO2 [ppm]")+
#   # ylab("Standard deviation CO2 [ppm]")+
#   # scale_x_log10()+
#   scale_y_log10()+
#   theme_article()





# ---- Differences between experts CO2 ----


ind_multiple_assess <- which(duplicated(table_draws$UniqueID))
list_multiple_experts <- unique(table_draws$UniqueID[ind_multiple_assess])
length(list_multiple_experts)

isF <- T
runID <- 0
df_multiple_users_co2 <- NULL
for(i in list_multiple_experts){
  runID = runID+1
  
  
  flux_blind = CO2_flux.auto$best.flux[which(CO2_flux.auto$UniqueID==i)]
  tab_users = CO2_flux.man[which(CO2_flux.man$UniqueID==i),]
  
  n = dim(tab_users)[1]
  
  if (n>=3){
    df_multiple_users_co2 <- rbind(df_multiple_users_co2, 
                                   data.frame(id = i,
                                              runID = as.factor(runID),
                                              var = "CO2",
                                              n = n,
                                              g.fact = CO2_flux.auto$g.fact[which(CO2_flux.auto$UniqueID==i)],
                                              LM.r2 = CO2_flux.auto$LM.r2[which(CO2_flux.auto$UniqueID==i)],
                                              LM.MAE = CO2_flux.auto$LM.MAE[which(CO2_flux.auto$UniqueID==i)],
                                              LM.RMSE = CO2_flux.auto$LM.RMSE[which(CO2_flux.auto$UniqueID==i)],
                                              flux_blind = flux_blind,
                                              flux_users_mean = mean(tab_users$best.flux),
                                              flux_users_sd = sd(tab_users$best.flux),
                                              flux_users_CV = sd(tab_users$best.flux)/mean(tab_users$best.flux)))
  }
}
df_multiple_users_co2$abs_diff <- abs(df_multiple_users_co2$flux_blind-df_multiple_users_co2$flux_users_mean)
df_multiple_users_co2$rel_diff <- (df_multiple_users_co2$flux_blind-df_multiple_users_co2$flux_users_mean)/df_multiple_users_co2$flux_users_mean*100


# ggplot(df_multiple_users_co2, aes(abs(flux_users_CV)))+geom_density()+scale_x_log10()

df_multiple_users_co2$is_gfact_over_2 <- "linear"
df_multiple_users_co2$is_gfact_over_2[df_multiple_users_co2$g.fact>=2] <- "non-linear"

df_multiple_users_co2$r2_over_0.9 <- "LM.r2 < 0.9"
df_multiple_users_co2$r2_over_0.9[df_multiple_users_co2$LM.r2>=.9] <- "LM.r2 >= 0.9"


length(which(df_multiple_users_co2$flux_users_CV<0.1))/dim(df_multiple_users_co2)[1]*100
length(which(df_multiple_users_co2$flux_users_CV>0.5))/dim(df_multiple_users_co2)[1]*100


# plot(log10(df_multiple_users_co2$LM.RMSE), log10(df_multiple_users_co2$LM.MAE))


plt_MAE_sd_co2 <- ggplot(df_multiple_users_co2, aes(LM.MAE, abs(flux_users_sd), colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM MAE FCO2 [\u00b5mol " * m^-2 ~ s^-1* "]"))+
  ylab("SD of experts' FCO2 [\u00b5mol " * m^-2 ~ s^-1* "]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  # scale_colour_manual(values = c("#1b9e77","#fc8d62"))+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_sd_co2


plt_MAE_CV_co2 <- ggplot(df_multiple_users_co2, aes(LM.MAE, abs(flux_users_CV), colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM MAE FCO2 [\u00b5mol " * m^-2 ~ s^-1* "]"))+
  ylab("CV of experts' FCO2 [-]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  # scale_colour_manual(values = c("#1b9e77","#fc8d62"))+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_CV_co2

plt_mean_CV_co2 <- ggplot(df_multiple_users_co2, aes(abs(flux_users_mean), abs(flux_users_CV), colour = log10(LM.MAE)))+
  geom_point()+
  
  xlab(expression(group("|", bar(F[CO2]^expert), "|")~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  
  # xlab(expression(tilde(*F[CO2]^expert)~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  ylab("CV of estimated FCO2 [dimensionless]")+
  scale_x_log10()+scale_y_log10()+theme_article()+
  # scale_colour_manual(values = c("#1b9e77","#fc8d62"))+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

# plt_mean_CV_co2 <- ggMarginal(plt_mean_CV_co2, type="density", groupColour = T)
plt_mean_CV_co2


p_expert_vs_blind_CO2 <- ggplot(df_multiple_users_co2, aes(flux_blind, flux_users_mean, colour = log10(LM.MAE)))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_errorbar(aes(ymin=flux_users_mean-flux_users_sd, 
                    ymax=flux_users_mean+flux_users_sd), linewidth = 1.2, alpha = 0.5)+
  geom_point()+
  xlab(expression(F[CO2]^auto~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(F[CO2]^expert~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  scale_size(range = c(0.1,5))+
  theme_article()+
  # scale_colour_manual(values = c("#1b9e77","#fc8d62"))+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)




p_out_CO2 <- ggpubr::ggarrange(p_expert_vs_blind_CO2, plt_MAE_sd_co2, nrow = 1, 
                                           common.legend = T, legend = "bottom", labels = c("a)","b)"))
p_out_CO2


ggsave(plot = p_out_CO2, filename = "Fig4_expert_vs_blind_CO2.jpeg", path = plots_path, 
       width = 6, height = 3.5, dpi = 300, units = 'in', scale = 1.2)




p_mean_SD_co2 <- ggplot(df_multiple_users_co2, aes(abs(flux_users_mean), abs(flux_users_sd), colour = log10(LM.MAE)))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab(expression(bar(F[CO2]^expert)~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("SD of "*F[CO2]^expert~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

# p_mean_SD_co2_dens <- ggMarginal(p_mean_SD_co2, type="density", groupColour = T)


p_mean_CV_co2 <- ggplot(df_multiple_users_co2, aes(abs(flux_users_mean), abs(flux_users_CV), colour = log10(LM.MAE)))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab(expression(bar(F[CO2]^expert)~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CO2]^expert))+
  theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

# p_mean_CV_co2_dens <- ggMarginal(p_mean_CV_co2, type="density", groupColour = T)

p_mean_SD_CV_co2 <- ggarrange(p_mean_SD_co2, p_mean_CV_co2, ncol = 2, common.legend = T)

# flux_ref = CO2_flux.man$best.flux
# flux_model = CO2_flux.auto$best.flux[match(CO2_flux.man$UniqueID, CO2_flux.auto$UniqueID)]
# df_exceedance_CO2 <- get_df_exceedance(abs(flux_model - flux_ref)/flux_model*100)
# ind_closest_10 <- which.min(abs(df_exceedance_CO2$t-10))
# 
# p_exceed <- ggplot(df_exceedance_CO2, aes(t, p))+geom_path()+geom_point()+
#   geom_hline(yintercept = 0)+
#   theme_bw()+
#   geom_segment(aes(x=-0,xend=10,
#                    y=df_exceedance_CO2$p[ind_closest_10], yend=df_exceedance_CO2$p[ind_closest_10]), color ="red")+
#   geom_segment(aes(x=10,xend=10,
#                    y=-Inf, yend=df_exceedance_CO2$p[ind_closest_10]), color ="red")+
#   # scale_x_log10()+
#   scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
#   ylab("Proportion of timeseries [%]")+
#   xlab("Relative difference [% of Expert flux]")+
#   ggtitle(paste0(round(df_exceedance_CO2$p[ind_closest_10]*10)/10,"% timeseries with < 10% difference"))
# 
# p_exceed
# ggsave(plot = p_exceed, filename = "overview_expert_vs_blind_co2.jpg", path = plots_path, 
#        width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)



# p2 <- ggplot(df_multiple_users_co2, aes(abs_diff, abs(flux_users_CV)))+geom_point()+
#   scale_x_log10()+
#   scale_y_log10()+
#   xlab("abs((Automated - mean(Expert))  [\u00b5mol m-2 s-1]")+
#   ylab("abs(CV(Expert)) [d.l.]")+
#   scale_size(range = c(0.1,5))+
#   scale_colour_viridis_c(direction = 1, option = "A", begin = 0.1, end = 0.9)+
#   theme_article()




# p1 <- p1 + geom_point(aes(size=g.fact, colour = g.fact))
# p2 <- p2 + geom_point(aes(size=g.fact, colour = g.fact))
# 
# p_expert_vs_blind_CO2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
# p_expert_vs_blind_CO2
# ggsave(plot = p_expert_vs_blind_CO2, filename = "expert_vs_blind_CO2_curvature.jpeg", path = plots_path, 
#        width = 8, height = 3.2, dpi = 300, units = 'in', scale = 1.0)




# plot incubations with the highest disagreements between experts
list_ids <- df_multiple_users_co2$id[which(abs(df_multiple_users_co2$flux_users_CV)>2)]

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_co2-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_co2-(table_draws_sel$start.time_auto)


p_disagree_CO2 <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf, fill=userID), alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  facet_wrap(UniqueID~., scales = "free")
p_disagree_CO2
ggsave(plot = p_disagree_CO2, filename = "expert_choice_CO2.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.0)




# plot incubations with the highest curvatures
list_ids <- df_multiple_users_co2$id[which(df_multiple_users_co2$g.fact>2)]

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_co2-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_co2-(table_draws_sel$start.time_auto)

ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf), fill="grey50", alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
  theme_article()+
  facet_wrap(UniqueID~., scales = "free")




# ---- Differences between experts CH4 --> COMPUTING ----

isF <- T
runID <- 0
df_multiple_users_CH4 <- NULL
for(i in list_multiple_experts){
  runID = runID+1
  
  i_flux_blind = which(CH4_flux.auto$UniqueID==i)
  tab_users = CH4_flux.man[which(CH4_flux.man$UniqueID==i),]
  tab_diff_users = CH4_diff_flux.man[which(CH4_diff_flux.man$UniqueID==i),]
  
  n = dim(tab_users)[1]
  
  if (n>=3){
    df_multiple_users_CH4 <- rbind(df_multiple_users_CH4, 
                                   data.frame(id = i,
                                              runID = as.factor(runID),
                                              var = "CH4",
                                              n = n,
                                              g.fact = CH4_flux.auto$g.fact[i_flux_blind],
                                              LM.r2 = CH4_flux.auto$LM.r2[which(CH4_flux.auto$UniqueID==i)],
                                              LM.MAE = CH4_flux.auto$LM.MAE[which(CH4_flux.auto$UniqueID==i)],
                                              
                                              total_blind = CH4_flux.auto$total.flux[i_flux_blind],
                                              total_users_mean = mean(tab_users$total.flux),
                                              total_users_sd = sd(tab_users$total.flux),
                                              total_users_CV = sd(tab_users$total.flux)/mean(tab_users$total.flux),
                                              
                                              ebull_blind = CH4_flux.auto$ebullition.flux[i_flux_blind],
                                              ebull_users_mean = mean(tab_users$ebullition.flux),
                                              ebull_users_sd = sd(tab_users$ebullition.flux),
                                              ebull_users_CV = sd(tab_users$ebullition.flux)/mean(tab_users$ebullition.flux),
                                              
                                              diffusion_blind = CH4_flux.auto$diffusion.flux[i_flux_blind],
                                              diffusion_users_mean = mean(tab_diff_users$best.flux),
                                              diffusion_users_sd = sd(tab_diff_users$best.flux),
                                              diffusion_users_CV = sd(tab_diff_users$best.flux)/mean(tab_diff_users$best.flux)))
  }
}


df_multiple_users_CH4$total_abs_diff <- abs(df_multiple_users_CH4$total_blind-df_multiple_users_CH4$total_users_mean)
df_multiple_users_CH4$total_rel_diff <- (df_multiple_users_CH4$total_blind-df_multiple_users_CH4$total_users_mean)/df_multiple_users_CH4$total_users_mean*100

df_multiple_users_CH4$diffusion_abs_diff <- abs(df_multiple_users_CH4$diffusion_blind-df_multiple_users_CH4$diffusion_users_mean)
df_multiple_users_CH4$diffusion_rel_diff <- (df_multiple_users_CH4$diffusion_blind-df_multiple_users_CH4$diffusion_users_mean)/df_multiple_users_CH4$diffusion_users_mean*100

df_multiple_users_CH4$ebull_abs_diff <- abs(df_multiple_users_CH4$ebull_blind-df_multiple_users_CH4$ebull_users_mean)
df_multiple_users_CH4$ebull_rel_diff <- (df_multiple_users_CH4$ebull_blind-df_multiple_users_CH4$ebull_users_mean)/df_multiple_users_CH4$ebull_users_mean*100


df_multiple_users_CH4$ebull_contrib <- df_multiple_users_CH4$ebull_blind/df_multiple_users_CH4$total_blind
df_multiple_users_CH4$ebull_contrib [which(df_multiple_users_CH4$ebull_contrib >1)] <- NA
df_multiple_users_CH4$Bubbling <- "< 10%"
df_multiple_users_CH4$Bubbling[which(df_multiple_users_CH4$ebull_contrib>0.1)] <- "> 10%"



names(df_multiple_users_CH4)


# ---- Differences between experts total CH4 ----

df_multiple_users_CH4$is_gfact_over_2 <- "linear"
df_multiple_users_CH4$is_gfact_over_2[df_multiple_users_CH4$g.fact>=2] <- "non-linear"

df_multiple_users_CH4$r2_over_0.9 <- "LM.r2 < 0.9"
df_multiple_users_CH4$r2_over_0.9[df_multiple_users_CH4$LM.r2>=.9] <- "LM.r2 >= 0.9"



length(which(df_multiple_users_CH4$total_users_CV<0.1))/dim(df_multiple_users_CH4)[1]*100
length(which(df_multiple_users_CH4$total_users_CV>0.5))/dim(df_multiple_users_CH4)[1]*100
length(which(df_multiple_users_CH4$total_users_CV>1))/dim(df_multiple_users_CH4)[1]*100



p_mean_SD_ch4tot <- ggplot(df_multiple_users_CH4, aes(abs(total_users_mean), abs(total_users_sd), colour = log10(LM.MAE)))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab(expression(bar(F[CH4~tot]^expert)~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("SD of "*F[CH4~tot]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

# p_mean_SD_co2_dens <- ggMarginal(p_mean_SD_co2, type="density", groupColour = T)


p_mean_CV_ch4tot <- ggplot(df_multiple_users_CH4, aes(abs(total_users_mean), abs(total_users_CV), colour = log10(LM.MAE)))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab(expression(bar(F[CH4~tot]^expert)~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CH4~tot]^expert))+
  theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)


p_mean_SD_CV_co2_ch4tot <- ggarrange(p_mean_SD_co2, p_mean_CV_co2,
                                     p_mean_SD_ch4tot, p_mean_CV_ch4tot, nrow = 2, ncol = 2, 
                                     labels = c("a)","b)",
                                                "c)","d)"),
                                     common.legend = T, legend = "bottom")

p_mean_SD_CV_co2_ch4tot

ggsave(plot = p_mean_SD_CV_co2_ch4tot, filename = "FigX_mean_SD_CV_Fco2_Fch4tot.jpeg", path = plots_path, 
       width = 8, height = 5.5, dpi = 300, units = 'in', scale = 1.2)






plt_MAE_sd_CH4tot <- ggplot(df_multiple_users_CH4, aes(LM.MAE, abs(total_users_sd), 
                                                       colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("SD of "*F[CH4~tot]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_sd_CH4tot


plt_MAE_sd_CH4ebull <- ggplot(df_multiple_users_CH4, aes(LM.MAE, abs(ebull_users_sd), 
                                                       colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("SD of "*F[CH4~ebull]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_sd_CH4ebull



plt_MAE_sd_CH4diffusion <- ggplot(df_multiple_users_CH4, aes(LM.MAE, abs(diffusion_users_sd), 
                                                             colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("SD of "*F[CH4~diff]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_sd_CH4diffusion




plt_MAE_CV_CH4tot <- ggplot(df_multiple_users_CH4, aes(LM.MAE, abs(total_users_CV), 
                                                       colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CH4~tot]^expert))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_CV_CH4tot



plt_MAE_CV_CH4ebull <- ggplot(df_multiple_users_CH4, aes(LM.MAE, abs(ebull_users_CV), 
                                                         colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CH4~ebull]^expert))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_CV_CH4ebull



plt_MAE_CV_CH4diffusion <- ggplot(df_multiple_users_CH4, aes(LM.MAE, abs(diffusion_users_CV), 
                                                             colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CH4~diff]^expert))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

plt_MAE_CV_CH4diffusion



plt_mean_CV_CH4tot <- ggplot(df_multiple_users_CH4, aes(abs(total_users_mean), abs(total_users_CV), colour = log10(LM.MAE)))+
  geom_point()+
  xlab(expression(bar(F[CH4~tot]^expert)~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CH4~tot]^expert))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)

# plt_mean_CV_CH4 <- ggMarginal(plt_mean_CV_CH4, type="density", groupColour = T)
plt_mean_CV_CH4tot


p_expert_vs_blind_CH4tot <- ggplot(df_multiple_users_CH4[which(df_multiple_users_CH4$total_users_mean > 0),],
                                   aes(total_blind, total_users_mean, colour = log10(LM.MAE)))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  # geom_point(aes(size = ebull_contrib))+
  geom_errorbar(aes(ymin=total_users_mean-total_users_sd, 
                    ymax=total_users_mean+total_users_sd), linewidth = 1.2, alpha = 0.5)+
  xlab(expression(F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(F[CH4~tot]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  theme_article()+
  scale_colour_viridis_c(direction = -1, option = "B", end = 0.9)


ggplot(df_multiple_users_CH4[which(df_multiple_users_CH4$total_users_mean > 0),],
       aes(total_blind, total_users_mean, colour = ebull_contrib))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  # geom_point(aes(size = ebull_contrib))+
  geom_errorbar(aes(ymin=total_users_mean-total_users_sd, 
                    ymax=total_users_mean+total_users_sd), linewidth = 1.2, alpha = 0.5)+
  xlab(expression(F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(F[CH4~tot]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  theme_article()+
  scale_colour_viridis_c(direction = -1, option = "A", end = 0.8)






ggplot(df_multiple_users_CH4, aes(abs(ebull_contrib), abs(diffusion_users_CV), 
                                  colour = Bubbling))+
  geom_point()+
  # scale_x_log10()+scale_y_log10()+
  theme_article()+
  scale_colour_viridis_d(direction = -1, option = "A", end = 0.8)






df_contrib <- df_multiple_users_CH4[,c("runID","total_users_CV","diffusion_users_CV","ebull_users_CV","ebull_contrib")]
names(df_contrib) <- c("runID","Total flux","Diffusion","Ebullition","ebull_contrib")
df_contrib_gath <- gather(df_contrib[!is.na(df_contrib$ebull_contrib),], var, CV, -runID, -ebull_contrib)

length(which(df_contrib$ebull_contrib>0.5))
length(which(df_contrib$diffusion>1 & df_contrib$ebull_contrib>0.5))

length(which(df_contrib$total>1 & df_contrib$ebull_contrib<0.5))
mean(df_contrib$diffusion[which(df_contrib$ebull_contrib>0.1 & df_contrib$ebull_contrib<0.5)])
sd(df_contrib$diffusion[which(df_contrib$ebull_contrib>0.1 & df_contrib$ebull_contrib<0.5)])

mean(df_contrib$total[which(df_contrib$ebull_contrib>0.1 & df_contrib$ebull_contrib<0.5)])
sd(df_contrib$total[which(df_contrib$ebull_contrib>0.1 & df_contrib$ebull_contrib<0.5)])

mean(df_contrib$diffusion[which(df_contrib$ebull_contrib>0.5)])
sd(df_contrib$diffusion[which(df_contrib$ebull_contrib>0.5)])
mean(df_contrib$total[which(df_contrib$ebull_contrib>0.5)])
sd(df_contrib$total[which(df_contrib$ebull_contrib>0.5)])

mean(df_contrib$ebullition[which(df_contrib$ebull_contrib>0.1)], na.rm = T)
sd(df_contrib$ebullition[which(df_contrib$ebull_contrib>0.1)], na.rm = T)




df_contrib_gath$ebull_ratio_range <- "0"

df_contrib_gath$ebull_ratio_range[which(abs(df_contrib_gath$ebull_contrib)>0 & 
                                           abs(df_contrib_gath$ebull_contrib)<=0.2)] <- "0 - 0.2"

df_contrib_gath$ebull_ratio_range[which(abs(df_contrib_gath$ebull_contrib )>0.2 & 
                                           abs(df_contrib_gath$ebull_contrib)<=0.4)] <- "0.2 - 0.4"

df_contrib_gath$ebull_ratio_range[which(abs(df_contrib_gath$ebull_contrib)>0.4 & 
                                           abs(df_contrib_gath$ebull_contrib)<=0.6)] <- "0.4 - 0.6"

df_contrib_gath$ebull_ratio_range[which(abs(df_contrib_gath$ebull_contrib )>0.6 & 
                                           abs(df_contrib_gath$ebull_contrib)<=0.8)] <- "0.6 - 0.8"

df_contrib_gath$ebull_ratio_range[which(abs(df_contrib_gath$ebull_contrib)>0.8)] <- "0.8 - 1"


df_contrib_gath$ebull_ratio_range <- factor(df_contrib_gath$ebull_ratio_range, levels = c("0","0 - 0.2", "0.2 - 0.4", "0.4 - 0.6",  "0.6 - 0.8",  "0.8 - 1"))


df_bubb_plt <- df_contrib_gath[which(!is.na(df_contrib_gath$ebull_contrib) & !is.na(df_contrib_gath$CV)),]

p_CVs_ch4_bubbling_a <- ggplot(df_bubb_plt[df_bubb_plt$var=="Ebullition",],
                             aes(ebull_ratio_range, abs(CV)))+
  geom_boxplot(outliers = F, alpha=0.5)+
  geom_jitter(width = 0.2, alpha=.2)+
  xlab("")+
  # xlab(expression("Ebullition ratio "*F[CH4~ebull]^auto~~"/"*~~F[CH4~tot]^auto))+
  ylab(expression("CV of "*F[CH4~ebull]^expert))+
  theme_article()+facet_wrap(var~.)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))


p_CVs_ch4_bubbling_b <- ggplot(df_bubb_plt[df_bubb_plt$var=="Diffusion",],
                               aes(ebull_ratio_range, abs(CV)))+
  geom_boxplot(outliers = F, alpha=0.5)+
  geom_jitter(width = 0.2, alpha=.2)+
  # xlab("")+
  xlab(expression("Ebullition ratio "*F[CH4~ebull]^auto~~"/"*~~F[CH4~tot]^auto))+
  ylab(expression("CV of "*F[CH4~diff]^expert))+
  theme_article()+facet_wrap(var~.)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))



p_CVs_ch4_bubbling_c <- ggplot(df_bubb_plt[df_bubb_plt$var=="Total flux",],
                               aes(ebull_ratio_range, abs(CV)))+
  geom_boxplot(outliers = F, alpha=0.5)+
  geom_jitter(width = 0.2, alpha=.2)+
  xlab("")+
  # xlab(expression("Ebullition ratio "*F[CH4~ebull]^auto~~"/"*~~F[CH4~tot]^auto))+
  ylab(expression("CV of "*F[CH4~tot]^expert))+
  theme_article()+facet_wrap(var~.)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))


p_CVs_ch4_bubbling <- ggpubr::ggarrange(p_CVs_ch4_bubbling_a, p_CVs_ch4_bubbling_b, p_CVs_ch4_bubbling_c,
                                        nrow = 1, labels = c("a","b","c"), align = "h")

p_CVs_ch4_bubbling

ggsave(plot = p_CVs_ch4_bubbling, filename = "Fig4_CVs_ch4_bubbling.jpeg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.2)



ggplot(df_multiple_users_CH4[!is.na(df_multiple_users_CH4$ebull_contrib),],
       aes(factor(ceiling(ebull_contrib*10/2)*2/10), abs(total_users_CV)))+
  # scale_x_log10()+
  scale_y_log10()+
  geom_jitter(width = 0.2, alpha=.2)+
  geom_boxplot(outliers = F, alpha=0.2)+
  xlab(expression(F[CH4~ebull]^auto~~"/"*~~F[CH4~tot]^auto))+
  ylab(expression("CV of "*F[CH4~tot]^expert))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()


head(df_multiple_users_CH4)

# ggplot(df_multiple_users_CH4, aes(total_blind))+geom_density()+theme_article()
# ggplot(df_multiple_users_CH4, aes(total_blind))+geom_histogram()+theme_article()



plt_mean_CV_CH4ebull <- ggplot(df_multiple_users_CH4, aes(abs(ebull_users_mean), abs(ebull_users_CV), colour = Bubbling))+
  geom_point()+
  xlab(expression(bar(F[CH4~ebull]^expert)~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CH4~ebull]^expert))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_d(direction = -1, option = "A", end = 0.8)

# plt_mean_CV_CH4ebull <- ggMarginal(plt_mean_CV_CH4ebull, type="density", groupColour = T)
plt_mean_CV_CH4ebull


plt_mean_CV_CH4diffusion <- ggplot(df_multiple_users_CH4, aes(abs(diffusion_users_mean), abs(diffusion_users_CV), colour = Bubbling))+
  geom_point()+
  xlab(expression(bar(F[CH4~diff]^expert)~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression("CV of "*F[CH4~diff]^expert))+
  scale_x_log10()+scale_y_log10()+theme_article()+
  scale_colour_viridis_d(direction = -1, option = "A", end = 0.8)

# plt_mean_CV_CH4diffusion <- ggMarginal(plt_mean_CV_CH4diffusion, type="density", groupColour = T)
plt_mean_CV_CH4diffusion


p_expert_vs_blind_CH4diffusion <- ggplot(df_multiple_users_CH4, 
                                         aes(diffusion_blind, diffusion_users_mean, colour = Bubbling))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbar(aes(ymin=diffusion_users_mean-diffusion_users_sd, 
                    ymax=diffusion_users_mean+diffusion_users_sd), linewidth = 1.2, alpha = .5)+
  xlab(expression(F[CH4~diff]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(F[CH4~diff]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  theme_article()+
  scale_colour_viridis_d(direction = -1, option = "A", end = 0.8)




# df_ebull <- df_multiple_users_CH4[df_multiple_users_CH4$ebull_blind>0 & df_multiple_users_CH4$ebull_users_mean>0,]
p_expert_vs_blind_CH4ebullition <- ggplot(df_multiple_users_CH4, aes(ebull_blind, ebull_users_mean, colour = Bubbling))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbar(aes(ymin=ebull_users_mean-ebull_users_sd, 
                    ymax=ebull_users_mean+ebull_users_sd), linewidth = 1.2, alpha = .5)+
  xlab(expression(F[CH4~ebull]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(F[CH4~ebull]^expert~~"[nmol " * m^-2 ~ s^-1* "]"))+
  theme_article()+
  scale_colour_viridis_d(direction = -1, option = "A", end = 0.8)



p_expert_vs_blind_CH4 <- ggpubr::ggarrange(plt_mean_CV_CH4tot, p_expert_vs_blind_CH4tot, 
                                           plt_mean_CV_CH4diffusion, p_expert_vs_blind_CH4diffusion,
                                           plt_mean_CV_CH4ebull, p_expert_vs_blind_CH4ebullition,
                                           nrow = 3, ncol = 2, labels = c("a)","b)",
                                                                          "c)","d)",
                                                                          "e)","f)"),
                                           common.legend = T, legend = "bottom")
p_expert_vs_blind_CH4
ggsave(plot = p_expert_vs_blind_CH4, filename = "Fig5_expert_vs_blind_CH4.jpeg", path = plots_path, 
       width = 8, height = 7.5, dpi = 300, units = 'in', scale = 1.2)




p_expert_vs_blind_CO2_CH4 <- ggpubr::ggarrange(p_expert_vs_blind_CO2, p_expert_vs_blind_CH4tot, 
                                               p_expert_vs_blind_CH4diffusion, p_expert_vs_blind_CH4ebullition,
                                                 nrow = 2, ncol = 2, labels = c("a)","b)",
                                                                                "c)","d)"),
                                                 common.legend = T, legend = "bottom")
p_expert_vs_blind_CO2_CH4
ggsave(plot = p_expert_vs_blind_CO2_CH4, filename = "Fig_experts_vs_blind.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.2)






p_expert_vs_blind_CO2_CH4 <- ggpubr::ggarrange(p_expert_vs_blind_CO2, p_expert_vs_blind_CH4tot, 
                                               plt_MAE_sd_co2, plt_MAE_sd_CH4tot,
                                               nrow = 2, ncol = 2, labels = c("a)","b)",
                                                                              "c)","d)"),
                                               common.legend = T, legend = "bottom")
p_expert_vs_blind_CO2_CH4
ggsave(plot = p_expert_vs_blind_CO2_CH4, filename = "FigX_experts_vs_blind_CO2_CH4tot.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.2)







p_LM_MAE_SD_experts_CO2_CH4 <- ggpubr::ggarrange(plt_MAE_sd_co2, plt_MAE_sd_CH4tot, 
                                           plt_MAE_sd_CH4diffusion, plt_MAE_sd_CH4ebull,
                                           nrow = 2, ncol = 2, labels = c("a)","b)",
                                                                          "c)","d)"),
                                           common.legend = T, legend = "bottom")
p_LM_MAE_SD_experts_CO2_CH4
ggsave(plot = p_LM_MAE_SD_experts_CO2_CH4, filename = "Fig_LM_MAE_SD_experts.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.2)




p_LM_MAE_CV_experts_CO2_CH4 <- ggpubr::ggarrange(plt_MAE_CV_co2, plt_MAE_CV_CH4tot, 
                                                 plt_MAE_CV_CH4diffusion, plt_MAE_CV_CH4ebull,
                                                 nrow = 2, ncol = 2, labels = c("a)","b)",
                                                                                "c)","d)"),
                                                 common.legend = T, legend = "bottom")
p_LM_MAE_CV_experts_CO2_CH4
ggsave(plot = p_LM_MAE_CV_experts_CO2_CH4, filename = "Fig_LM_MAE_CV_experts.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.2)






p_avg_SD_CV_expert_vs_auto_CO2_CH4 <- ggpubr::ggarrange(p_mean_SD_co2, p_mean_CV_co2,p_expert_vs_blind_CO2,
                                                        p_mean_SD_ch4tot, p_mean_CV_ch4tot,p_expert_vs_blind_CH4tot,
                                               nrow = 2, ncol = 3, labels = c("a","b","c",
                                                                              "d","e","f"),
                                               common.legend = T, legend = "bottom")
p_avg_SD_CV_expert_vs_auto_CO2_CH4
ggsave(plot = p_avg_SD_CV_expert_vs_auto_CO2_CH4, filename = "FigX_avg_SD_CV_expert_vs_auto_CO2_CH4tot.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.2)




p_avg_SD_expert_vs_auto_CO2_CH4 <- ggpubr::ggarrange(p_mean_SD_co2,p_expert_vs_blind_CO2,
                                                        p_mean_SD_ch4tot,p_expert_vs_blind_CH4tot,
                                                        nrow = 2, ncol = 2, labels = c("a","b",
                                                                                       "c","d"),
                                                        common.legend = T, legend = "bottom")
p_avg_SD_expert_vs_auto_CO2_CH4
ggsave(plot = p_avg_SD_expert_vs_auto_CO2_CH4, filename = "FigX_avg_SD_expert_vs_auto_CO2_CH4tot.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.2)







# ---------- relationship between CV CO2 and CV Ch4 -----------------


df_multiple_users_co2$flux_users_CV_ch4 <- df_multiple_users_CH4$total_users_CV[match(df_multiple_users_co2$id, df_multiple_users_CH4$id)]
# df_multiple_users_co2$CV_ch4 <- stats_all$CV_ch4[match(df_multiple_users_co2$id, stats_all$UniqueID)]
df_multiple_users_co2$FCH4ebull_blind <- df_multiple_users_CH4$ebull_blind [match(df_multiple_users_co2$id, df_multiple_users_CH4$id)]
df_multiple_users_co2$FCH4tot_blind <- df_multiple_users_CH4$total_blind [match(df_multiple_users_co2$id, df_multiple_users_CH4$id)]

df_multiple_users_co2$ebull_contrib <- df_multiple_users_co2$FCH4ebull_blind/df_multiple_users_co2$FCH4tot_blind

df_multiple_users_co2$ebull_contrib[which(df_multiple_users_co2$ebull_contrib > 1)] <- NA

p_CVs_co2_ch4 <- ggplot(df_multiple_users_co2, aes(abs(flux_users_CV), abs(flux_users_CV_ch4)))+
                          geom_point(aes(size = FCH4ebull_blind))+
  scale_x_log10()+
  scale_y_log10()+
  # geom_smooth(method = "lm", fill="grey", color = "grey25")+
  xlab(expression("CV of "*F[CO2]^expert))+
  ylab(expression("CV of "*F[CH4~tot]^expert))+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = 1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()
p_CVs_co2_ch4

cor(df_multiple_users_co2$flux_users_CV, df_multiple_users_co2$flux_users_CV_ch4)
summary(lm(formula = flux_users_CV_ch4 ~ flux_users_CV, data = df_multiple_users_co2))


ggplot(df_multiple_users_co2,
       aes(factor(ceiling(abs(flux_users_CV_ch4)*10/2)*2/10), abs(flux_users_CV)))+
  # scale_x_log10()+
  scale_y_log10()+
  geom_jitter()+
  geom_boxplot(outliers = F, alpha=0.2)+
  xlab(expression("CV of "*F[CH4~tot]^expert))+
  ylab(expression("CV of "*F[CH2]^expert))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()


df_multiple_users_co2$CV_ch4_range <- ""

df_multiple_users_co2$CV_ch4_range[which(abs(df_multiple_users_co2$flux_users_CV_ch4)>0 & 
                                           abs(df_multiple_users_co2$flux_users_CV_ch4)<=0.05)] <- "0.0-0.05"

df_multiple_users_co2$CV_ch4_range[which(abs(df_multiple_users_co2$flux_users_CV_ch4)>0.05 & 
                                           abs(df_multiple_users_co2$flux_users_CV_ch4)<=0.1)] <- "0.05-0.1"

df_multiple_users_co2$CV_ch4_range[which(abs(df_multiple_users_co2$flux_users_CV_ch4)>0.1 & 
                                           abs(df_multiple_users_co2$flux_users_CV_ch4)<=0.2)] <- "0.1-0.2"

df_multiple_users_co2$CV_ch4_range[which(abs(df_multiple_users_co2$flux_users_CV_ch4)>0.2 & 
                                           abs(df_multiple_users_co2$flux_users_CV_ch4)<=0.3)] <- "0.2-0.3"

df_multiple_users_co2$CV_ch4_range[which(abs(df_multiple_users_co2$flux_users_CV_ch4)>0.3 & 
                                           abs(df_multiple_users_co2$flux_users_CV_ch4)<=0.4)] <- "0.3-0.4"

df_multiple_users_co2$CV_ch4_range[which(abs(df_multiple_users_co2$flux_users_CV_ch4)>0.4)] <- ">0.4"

df_multiple_users_co2$CV_ch4_range <- factor(df_multiple_users_co2$CV_ch4_range, levels = c("0.0-0.05", "0.05-0.1", "0.1-0.2",  "0.2-0.3",  "0.3-0.4",  ">0.4",""))


df_multiple_users_co2$significantEbull <- "FCH4 diff > FCH4 ebull"
df_multiple_users_co2$significantEbull[which(df_multiple_users_co2$ebull_contrib>.5)] <- "FCH4 diff < FCH4 ebull"

p_CVs_co2_ch4_boxplot <- ggplot(df_multiple_users_co2[which(df_multiple_users_co2$CV_ch4_range != ""),],
       aes(factor(CV_ch4_range), abs(flux_users_CV)))+
  scale_y_log10()+
  geom_jitter(width = 0.2, aes(colour = significantEbull), size=3, alpha=.5)+
  geom_boxplot(outliers = F, alpha=0.6)+
  xlab(expression("CV of "*F[CH4~tot]^expert))+
  ylab(expression("CV of "*F[C02]^expert))+
  scale_colour_viridis_d(direction = 1, option = "A", begin = 0.1, end = 0.75)+
  theme_article()+theme(legend.title = element_blank(), legend.position = c(0.6,0.1))
p_CVs_co2_ch4_boxplot


ggsave(plot = p_CVs_co2_ch4_boxplot, filename = "Fig4_CV_ch4_vs_CV_co2_boxplot.jpeg", path = plots_path, 
       width = 3, height = 2.5, dpi = 300, units = 'in', scale = 1.5)


# examples of high CVs
table(df_multiple_users_co2$significantEbull[which(df_multiple_users_co2$CV_ch4_range==">0.4")])
mean(df_multiple_users_co2$flux_users_CV[which(df_multiple_users_co2$CV_ch4_range==">0.4")])

list_highCVs <- df_multiple_users_co2[which(df_multiple_users_co2$flux_users_CV>.5 & df_multiple_users_co2$CV_ch4_range==">0.4"),]

list_ids <- list_highCVs$id
mydata_sel <- load_this(mylist = list_ids)
# 
# ggplot()+
#   geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
#   theme_article()+
#   xlab("Elapsed time [secs]")+
#   facet_wrap(UniqueID~., scales = "free")


# list_ids <- c("s2-du-p2-1-o-d-09:17", "s1-cu-a2-16-o-d-11:58")
# 
# list_ids <- c("s2-du-p2-1-o-d-09:17", "s1-cu-a2-14-o-d-11:32")
# 
# list_ids <- c("s2-cu-r2-10-o-d-11:10", "s1-cu-a2-16-o-d-11:58")

list_ids <- c("s1-cu-a2-16-o-d-11:58","s3-va-a1-4-o-d-10:04")


mydata_sel <- load_this(mylist = list_ids)
# table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
# table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
# table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)


CVs_co2 <- c(abs(round(list_highCVs$flux_users_CV[which(list_highCVs$id==list_ids[1])]*100)/100),
             abs(round(list_highCVs$flux_users_CV[which(list_highCVs$id==list_ids[2])]*100)/100))
mydata_sel$incubation <- ""
mydata_sel$incubation[which(mydata_sel$UniqueID == list_ids[1])] <- "incubation A"
mydata_sel$incubation[which(mydata_sel$UniqueID == list_ids[2])] <- "incubation B"


CVs_ch4 <- c(abs(round(list_highCVs$flux_users_CV_ch4[which(list_highCVs$id==list_ids[1])]*100)/100),
             abs(round(list_highCVs$flux_users_CV_ch4[which(list_highCVs$id==list_ids[2])]*100)/100))
mydata_sel$incubation <- ""
mydata_sel$incubation[which(mydata_sel$UniqueID == list_ids[1])] <- "incubation A"
mydata_sel$incubation[which(mydata_sel$UniqueID == list_ids[2])] <- "incubation B"

xs <- 0.7*c((max(mydata_sel$Etime[mydata_sel$incubation=="incubation A"])),
           (max(mydata_sel$Etime[mydata_sel$incubation=="incubation B"])))

y_co2 <- c(0.1*(max(mydata_sel$CO2dry_ppm[mydata_sel$incubation=="incubation A"])-min(mydata_sel$CO2dry_ppm[mydata_sel$incubation=="incubation A"]))+min(mydata_sel$CO2dry_ppm[mydata_sel$incubation=="incubation A"]),
           0.1*(max(mydata_sel$CO2dry_ppm[mydata_sel$incubation=="incubation B"])-min(mydata_sel$CO2dry_ppm[mydata_sel$incubation=="incubation B"]))+min(mydata_sel$CO2dry_ppm[mydata_sel$incubation=="incubation B"]))

           
y_ch4 <- c(0.1*(max(mydata_sel$CH4dry_ppb[mydata_sel$incubation=="incubation A"])-min(mydata_sel$CH4dry_ppb[mydata_sel$incubation=="incubation A"]))+min(mydata_sel$CH4dry_ppb[mydata_sel$incubation=="incubation A"]),
           0.1*(max(mydata_sel$CH4dry_ppb[mydata_sel$incubation=="incubation B"])-min(mydata_sel$CH4dry_ppb[mydata_sel$incubation=="incubation B"]))+min(mydata_sel$CH4dry_ppb[mydata_sel$incubation=="incubation B"]))



data_geomtext <- data.frame(incubation = c("incubation A","incubation B"),
                            incubation_short = c("A","B"),
                            xs = xs,
                            y_co2 = y_co2,
                            y_ch4 = y_ch4,
                            label_co2 = paste0("CV_FCO2=",CVs_co2),
                            label_ch4 = paste0("CV_FCH4tot=",CVs_ch4),
                            CV_ch4 = CVs_ch4,
                            CV_co2 = CVs_co2,
                            CV_ch4_range = as.factor(c(">0.4",">0.4")))



p_ex_co2 <- ggplot()+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = incubation))+
  theme_article()+
  geom_text(data = data_geomtext, aes(xs, y_co2, label = label_co2), parse = F, size=3)+
  xlab("Elapsed time [secs]")+
  ylab(expression(""*CO[2~~dry]~"[ppm]"))+
  facet_wrap(incubation~., scales = "free")


p_ex_ch4 <- ggplot()+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = incubation))+
  theme_article()+
  geom_text(data = data_geomtext, aes(xs, y_ch4, label = label_ch4), size=3)+
  xlab("Elapsed time [secs]")+
  ylab(expression(""*CH[4~~dry]~"[ppb]"))+
  facet_wrap(incubation~., scales = "free")



p_CVs_co2_ch4 <- ggplot(df_multiple_users_co2,
                        aes(abs(flux_users_CV_ch4), abs(flux_users_CV)))+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(aes(colour = significantEbull), size=3, alpha=.5)+
  geom_point(data = df_multiple_users_co2[df_multiple_users_co2$id %in% list_ids,], 
             aes(colour = significantEbull), size=5, alpha=.7)+
  
  geom_text(data =data_geomtext, aes(CV_ch4, CV_co2, label = incubation_short), nudge_x=0.02, nudge_y = -0.1, size=4)+
  xlab(expression("CV of "*F[CH4~tot]^expert))+
  ylab(expression("CV of "*F[CO2]^expert))+
  scale_colour_viridis_d(direction = 1, option = "A", begin = 0.1, end = 0.75)+
  theme_article()+theme(legend.title = element_blank(), legend.position = c(0.6,0.1))


p_CVs_co2_ch4_boxplot <- ggplot(df_multiple_users_co2[which(df_multiple_users_co2$CV_ch4_range != ""),],
                                aes(factor(CV_ch4_range), abs(flux_users_CV)))+
  scale_y_log10()+
  geom_jitter(width = 0.2, aes(colour = significantEbull), size=3, alpha=.5)+
  geom_boxplot(outliers = F, alpha=0.6)+
  geom_text(data =data_geomtext, aes(CV_ch4, CV_co2, label = incubation_short), nudge_x=0.5, size=3.5)+
  xlab(expression("CV of "*F[CH4~tot]^expert))+
  ylab(expression("CV of "*F[C02]^expert))+
  scale_colour_viridis_d(direction = 1, option = "A", begin = 0.1, end = 0.75)+
  theme_article()+theme(legend.title = element_blank(), legend.position = c(0.6,0.1))#+ theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
# p_CVs_co2_ch4_boxplot


p_CVs_co2_ch4_boxplot_with_timeseries <- ggarrange(p_CVs_co2_ch4, ggarrange(p_ex_co2, p_ex_ch4, ncol = 1, align = "hv", labels=c("b","c")), ncol = 2, labels = c("a","",""))
p_CVs_co2_ch4_boxplot_with_timeseries

ggsave(plot = p_CVs_co2_ch4_boxplot_with_timeseries, filename = "Fig5_CV_ch4_vs_CV_co2_boxplot.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)



# ---------- distance to blind vs LM.MAE -----------------

# ggplot(df_multiple_users_CH4, aes(abs(total_blind), abs(total_users_mean/total_blind)))+geom_point()+
#   scale_x_log10()+
#   scale_y_log10()+
#   theme_article()
# 
# df_multiple_users_CH4[which(df_multiple_users_CH4$total_users_mean/df_multiple_users_CH4$total_blind>100),]
# df_multiple_users_CH4[which(df_multiple_users_CH4$total_users_mean/df_multiple_users_CH4$total_blind<0.1),]



p_MAE_absdiff_co2 <- ggplot(df_multiple_users_co2, aes(abs(LM.MAE), abs(abs_diff)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CO2]^auto~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(group("|", F[CO2]^auto~~-~~bar(F[CO2]^expert), "|")~~"[\u00b5mol " * m^-2 ~ s^-1* "]"))+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()

p_MAE_absdiff_ch4tot <- ggplot(df_multiple_users_CH4, aes(abs(LM.MAE), abs(total_abs_diff)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(group("|", F[CH4~tot]^auto~~-~~bar(F[CH4~tot]^expert), "|")~~"[nmol " * m^-2 ~ s^-1* "]"))+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()

p_MAE_absdiff_co2_ch4tot <- ggarrange(p_MAE_absdiff_co2, p_MAE_absdiff_ch4tot, align = "hv", labels = c("a","b"))
p_MAE_absdiff_co2_ch4tot

ggsave(plot = p_MAE_absdiff_co2_ch4tot, filename = "Fig4_lmMAE_absdiff_co2_ch4tot.jpeg", path = plots_path, 
       width = 8, height = 3.5, dpi = 300, units = 'in', scale = 1.2)



ggplot(df_multiple_users_CH4, aes(abs(LM.MAE), abs(ebull_blind/total_blind)))+
  geom_point()+
  xlab(expression("LM.MAE of "*F[CH4~tot]^auto~~"[nmol " * m^-2 ~ s^-1* "]"))+
  ylab(expression(group("|", F[CH4~tot]^auto~~-~~bar(F[CH4~tot]^expert), "|")~~"[nmol " * m^-2 ~ s^-1* "]"))+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()



# ---------------------------





flux_ref = CH4_flux.man$total.flux
flux_model = CH4_flux.auto$total.flux[match(CH4_flux.man$UniqueID, CH4_flux.auto$UniqueID)]
df_exceedance_CH4 <- get_df_exceedance(abs(flux_model - flux_ref)/flux_model*100)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")+
  ggtitle(paste0(round(df_exceedance_CH4$p[ind_closest_10]*10)/10,"% timeseries with < 10% difference"))

p_exceed
ggsave(plot = p_exceed, filename = "overview_expert_vs_blind_ch4_total_flux.jpg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)



df_multiple_users_CH4 <- df_multiple_users_CH4[df_multiple_users_CH4$ebull_contrib >=0 & df_multiple_users_CH4$ebull_contrib<=1,]


df_multiple_users_co2$ebull_contrib <- df_multiple_users_CH4$ebull_contrib[match(df_multiple_users_co2$id, df_multiple_users_CH4$id)]



p_CVs_co2_ch4 <- ggplot(df_multiple_users_co2, aes(abs(flux_users_CV), abs(flux_users_CV_ch4), colour = ebull_contrib, size=ebull_contrib))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("CO2 CV(Expert) [d.l.]")+
  ylab("CH4 CV(Expert) [d.l.]")+
  scale_size(range = c(1,5))+
  scale_colour_viridis_c(direction = 1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()

ggsave(plot = p_CVs_co2_ch4, filename = "CV_ch4_vs_CV_co2_ebullcontrib.jpeg", path = plots_path, 
       width = 5, height = 4, dpi = 300, units = 'in', scale = 1.0)





p_tot <- ggplot(df_multiple_users_CH4[!is.na(df_multiple_users_CH4$ebull_contrib),],
                aes(factor(ceiling(ebull_contrib*10/2)*2/10), total_abs_diff))+
  # scale_x_log10()+
  # scale_y_log10()+
  geom_jitter()+
  geom_boxplot(outliers = F, alpha=0.2)+
  xlab("Blind bubbling contribution [d.l.]")+
  ylab("abs((Blind - mean(Expert)) [ch4 flux]")+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  ggtitle("Total CH4 flux")+
  theme_article()


p_diff <- ggplot(df_multiple_users_CH4[!is.na(df_multiple_users_CH4$ebull_contrib),],
                 aes(factor(ceiling(ebull_contrib*10/2)*2/10), diffusion_abs_diff))+
  # scale_x_log10()+
  # scale_y_log10()+
  geom_jitter()+
  geom_boxplot(outliers = F, alpha=0.2)+
  xlab("Blind bubbling contribution [d.l.]")+
  ylab("abs((Blind - mean(Expert)) [ch4 flux]")+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  ggtitle("Diffusion CH4 flux")+
  theme_article()

p_ebull <- ggplot(df_multiple_users_CH4[!is.na(df_multiple_users_CH4$ebull_contrib),],
                  aes(factor(ceiling(ebull_contrib*10/2)*2/10), ebull_abs_diff))+
  # scale_x_log10()+
  # scale_y_log10()+
  geom_jitter()+
  geom_boxplot(outliers = F, alpha=0.2)+
  xlab("Blind bubbling contribution [d.l.]")+
  ylab("abs((Blind - mean(Expert)) [ch4 flux]")+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  ggtitle("Ebullition CH4 flux")+
  theme_article()


ggarrange(p_tot, p_diff, p_ebull, ncol = 3)


p2 <- ggplot(df_multiple_users_CH4[which(df_multiple_users_CH4$total_users_mean > 0),], 
             aes(total_abs_diff, abs(total_users_CV)))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("abs((Automated - mean(Expert)) [nmol m-2 s-1]")+
  ylab("abs(CV(Expert)) [d.l.]")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = 1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()

p_expert_vs_blind_CH4tot <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
p_expert_vs_blind_CH4tot
ggsave(plot = p_expert_vs_blind_CH4tot, filename = "expert_vs_blind_CH4tot.jpeg", path = plots_path, 
       width = 8, height = 3.2, dpi = 300, units = 'in', scale = 1.0)

p1 <- p1 + geom_point(aes(colour = ebull_contrib, size = ebull_contrib))
p2 <- p2 + geom_point(aes(colour = ebull_contrib, size = ebull_contrib))
p_expert_vs_blind_CH4tot <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
p_expert_vs_blind_CH4tot
ggsave(plot = p_expert_vs_blind_CH4tot, filename = "expert_vs_blind_CH4tot_ebullContrib.jpeg", path = plots_path, 
       width = 8, height = 3.2, dpi = 300, units = 'in', scale = 1.0)


# plot incubations with the highest disagreements between experts
list_ids <- tail(df_multiple_users_CH4$id[order(df_multiple_users_CH4$total_users_CV)], 9)

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

p_disagree_CH4tot <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf, fill=userID), 
            # fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  facet_wrap(UniqueID~., scales = "free")

ggsave(plot = p_disagree_CH4tot, filename = "expert_choice_CH4tot.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.0)


# ---- Differences between experts DIFFUSION CH4 ----

flux_ref = CH4_diff_flux.man$best.flux
flux_model = CH4_flux.auto$diffusion.flux[match(CH4_diff_flux.man$UniqueID, CH4_flux.auto$UniqueID)]
df_exceedance_CH4 <- get_df_exceedance(abs(flux_model - flux_ref)/flux_model*100)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")+
  ggtitle(paste0(round(df_exceedance_CH4$p[ind_closest_10]*10)/10,"% timeseries with < 10% difference"))

p_exceed
ggsave(plot = p_exceed, filename = "overview_expert_vs_blind_ch4_diffusion.jpg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)





p2 <- ggplot(df_multiple_users_CH4, 
             aes(diffusion_abs_diff, abs(diffusion_users_CV)))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("abs((Automated - mean(Expert)) [nmol m-2 s-1]")+
  ylab("abs(CV(Expert)) [d.l.]")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()

p_expert_vs_blind_CH4diff <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
ggsave(plot = p_expert_vs_blind_CH4diff, filename = "expert_vs_blind_CH4diff.jpeg", path = plots_path, 
       width = 8, height = 3.2, dpi = 300, units = 'in', scale = 1.0)


p1 <- p1 + geom_point(aes(colour = ebull_contrib, size = ebull_contrib))
p2 <- p2 + geom_point(aes(colour = ebull_contrib, size = ebull_contrib))
p_expert_vs_blind_CH4diff <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
p_expert_vs_blind_CH4diff
ggsave(plot = p_expert_vs_blind_CH4diff, filename = "expert_vs_blind_CH4diff_ebullContrib.jpeg", path = plots_path, 
       width = 8, height = 3.2, dpi = 300, units = 'in', scale = 1.0)



# plot incubations with the highest disagreements between experts
list_ids <- tail(df_multiple_users_CH4$id[order(df_multiple_users_CH4$diffusion_users_CV)], 9)

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

table_draws_sel[order(table_draws_sel$UniqueID),]

p_disagree_CH4diff <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf, fill=userID), 
            # fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  facet_wrap(UniqueID~., scales = "free")
ggsave(plot = p_disagree_CH4diff, filename = "expert_choice_CH4diff.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.0)




# ---- Differences between experts EBULLITION CH4 ----


flux_ref = CH4_flux.man$ebullition.flux[CH4_flux.man$ebullition.flux>0]
flux_model = CH4_flux.auto$ebullition.flux[match(CH4_flux.man$UniqueID[CH4_flux.man$ebullition.flux>0], CH4_flux.auto$UniqueID)]
df_exceedance_CH4 <- get_df_exceedance(abs(flux_model - flux_ref)/flux_model*100)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")+
  ggtitle(paste0(round(df_exceedance_CH4$p[ind_closest_10]*10)/10,"% timeseries with < 10% difference"))

p_exceed
ggsave(plot = p_exceed, filename = "overview_expert_vs_blind_ch4_diffusion.jpg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)



ggplot(df_multiple_users_CH4, aes(ebull_contrib))+geom_density()+theme_article()




p2 <- ggplot(df_ebull, aes(ebull_abs_diff, abs(ebull_users_CV)))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("abs((Automated - mean(Expert)) [nmol m-2 s-1]")+
  ylab("abs(CV(Expert)) [d.l.]")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()


p_expert_vs_blind_CH4ebull <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
ggsave(plot = p_expert_vs_blind_CH4ebull, filename = "expert_vs_blind_CH4ebull.jpeg", path = plots_path, 
       width = 8, height = 3.2, dpi = 300, units = 'in', scale = 1.0)

p1 <- p1 + geom_point(aes(colour = ebull_contrib, size = ebull_contrib))
p2 <- p2 + geom_point(aes(colour = ebull_contrib, size = ebull_contrib))
p_expert_vs_blind_CH4ebull <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
p_expert_vs_blind_CH4ebull
ggsave(plot = p_expert_vs_blind_CH4ebull, filename = "expert_vs_blind_CH4ebull.jpeg_ebullContrib.jpeg", path = plots_path, 
       width = 8, height = 3.2, dpi = 300, units = 'in', scale = 1.0)


# plot incubations with the highest disagreements

list_ids <- tail(df_multiple_users_CH4$id[order(df_multiple_users_CH4$ebull_users_mean)], 9)

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

p_disagree_CH4ebull <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf, fill= userID), 
            # fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  facet_wrap(UniqueID~., scales = "free")
ggsave(plot = p_disagree_CH4ebull, filename = "expert_choice_CH4ebull.jpeg", path = plots_path, 
       width = 8, height = 6, dpi = 300, units = 'in', scale = 1.0)




# examples of very weird choices from dear experts
list_ids <- c("s3-da-r1-1-o-d-08:22", "s4-cu-a1-6-o-d-08:44", "s4-da-r1-5-o-d-10:01")
mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf, fill = username), 
            # fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  facet_grid(username~UniqueID, scales = "free")


# ---- relationship between disagreements on CH4diff and CH4tot or CH4ebull? ----


ggplot(df_multiple_users_CH4, aes(total_rel_diff, diffusion_rel_diff))+
  geom_point()+
  geom_smooth(method = 'lm')+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("Total abs((Blind - mean(Expert))/mean(Expert) [%]")+
  ylab("Diffusion abs((Blind - mean(Expert))/mean(Expert) [%]")








head(CH4_flux.man$uniqAssessID, 50)
head(CH4_diff_flux.man$uniqAssessID, 50)

CH4_diff_flux.man$diffusion.flux <- CH4_diff_flux.man$best.flux
CH4_diff_flux.man$ebullition.flux <- CH4_flux.man$ebullition.flux[match(CH4_diff_flux.man$uniqAssessID, CH4_flux.man$uniqAssessID)]
CH4_diff_flux.man$total.flux <- CH4_flux.man$total.flux[match(CH4_diff_flux.man$uniqAssessID, CH4_flux.man$uniqAssessID)]
CH4_diff_flux.man$total.flux_auto <- CH4_flux.auto$best.flux[match(CH4_diff_flux.man$UniqueID, CH4_flux.auto$UniqueID)]
CH4_diff_flux.man$ebullition.flux_auto <- CH4_flux.auto$ebullition.flux[match(CH4_diff_flux.man$UniqueID, CH4_flux.auto$UniqueID)]

CH4_diff_flux.man$ebullition.flux_bydifference <- CH4_diff_flux.man$total.flux - CH4_diff_flux.man$diffusion.flux


ggplot(CH4_diff_flux.man, aes(total.flux, total.flux_auto))+geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Total CH4 flux based on delta C over incubation")+
  ylab("Total CH4 based on curve fitting")+theme_article()



ggplot(CH4_diff_flux.man, aes(ebullition.flux, ebullition.flux_bydifference))+geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  # scale_x_log10()+
  # scale_y_log10()+
  xlab("CH4 ebullition based on flux.separator")+
  ylab("CH4 ebullition = total - manual diffusion")+theme_article()


ggplot(CH4_diff_flux.man, aes(ebullition.flux, ebullition.flux_auto))+geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Ebullition with flux.separator after manual selection")+
  ylab("Ebullition with flux.separator no manual selection")+theme_article()


head(df_multiple_users_CH4)

ggplot(df_multiple_users_CH4, aes(ebull_abs_diff, diffusion_abs_diff))+geom_point()+
  scale_x_log10()+
  scale_y_log10()



plot()

# ---- Do we really need flux.separator to estimate total.flux? ----

CH4_flux.auto$ebullition.contrib <- CH4_flux.auto$ebullition.flux / CH4_flux.auto$total.flux

flux_ref = CH4_flux.auto$total.flux
flux_model = CH4_flux.auto$best.flux
df_exceedance_CH4 <- get_df_exceedance(abs(flux_model - flux_ref)/flux_model*100)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")+
  ggtitle(paste0(round(df_exceedance_CH4$p[ind_closest_10]*10)/10,"% timeseries with < 10% difference"))

p_exceed

p_goflux_vs_deltaC <- ggplot(CH4_flux.auto[CH4_flux.auto$ebullition.contrib >=0 & CH4_flux.auto$ebullition.contrib <= 1,], 
                             aes(best.flux, total.flux, size = ebullition.contrib))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  scale_size(range = c(0.1,3))+
  xlab("Total CH4 flux with goFlux [nmol m-2 s-1]")+
  ylab("Total CH4 flux using delta C [nmol m-2 s-1]")+
  theme_article()+
  ggtitle(paste0(round(df_exceedance_CH4$p[ind_closest_10]*10)/10,"% timeseries with < 10% difference"))


ggsave(plot = p_goflux_vs_deltaC, filename = "goflux_vs_ch4_total_flux.jpg", path = plots_path,
       width = 6, height = 4, dpi = 300, units = 'in', scale = 1.0)






# ---- An example of how flux.separator works ----

require(goFlux)
require(ggnewscale)

auxfile$obs.length <- auxfile$duration

list_ids <- tail(df_multiple_users_CH4$id[order(df_multiple_users_CH4$ebull_contrib)], 9)

mydata_sel <- load_this(mylist = list_ids)
# print(plot_incubations(mydata_sel))

automaticflux(dataframe = mydata_sel, myauxfile = auxfile, shoulder = 0, 
              gastype = "CH4dry_ppb", fluxSeparation = T, displayPlots = T, method = "trust.it.all")


clickflux(dataframe = mydata_sel, myauxfile = auxfile, shoulder = 0, 
          gastype = "CH4dry_ppb", fluxSeparation = T, displayPlots = T, plot.lim = c(2000, max(mydata_sel$CH4dry_ppb)))




# ---- Does CH4 bubbling cause major disruptions in CO2 flux measurements ? ----

stats_all$ebullition.detected <- stats_all$ebull_ch4>0

ggplot(stats_all, aes(n_flag_co2+1, colour = ebullition.detected))+
  geom_density(alpha=0.5)+
  scale_x_log10()+
  # scale_y_log10()+
  scale_colour_viridis_d(option = "B", begin = 0, end = 0.7)+
  theme_article()+theme(legend.position = c(0.8,0.2))



list_ids <- stats_all$UniqueID[which(stats_all$ebull_ch4>100 & stats_all$ebull_ch4<1000)]


mydata_sel <- load_this(mylist = list_ids)
mydata_sel<- mydata_sel[order(mydata_sel$POSIX.time),]
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_co2-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_co2-(table_draws_sel$start.time_auto)

ggplot()+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
  theme_article()+
  facet_wrap(.~UniqueID, scales = "free")


ggplot()+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  facet_wrap(.~UniqueID, scales = "free")


ggplot()+
  geom_path(data = mydata_sel, aes(CO2dry_ppm, CH4dry_ppb, group = UniqueID))+
  theme_article()+
  facet_wrap(.~UniqueID, scales = "free")


print(plot_incubations(mydata_sel))



list_ids <- df_ebull$id[which(df_ebull$ebull_contrib>.9)]
mydata_sel <- load_this(mylist = list_ids)
mydata_sel<- mydata_sel[order(mydata_sel$POSIX.time),]
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf, fill = username), 
            # fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
  theme_article()+
  facet_wrap(.~UniqueID, scales = "free")


ggplot()+
  geom_path(data = mydata_sel, aes(CO2dry_ppm, CH4dry_ppb, group = UniqueID))+
  theme_article()+
  facet_wrap(.~UniqueID, scales = "free")


mydata_sel$Etime <- as.numeric(mydata_sel$Etime)
print(plot_incubations(dataframe = mydata_sel))





p_mean_sd <- ggplot(stats_all, aes(mean_ch4, sd_ch4, colour = ebullition.detected))+
  geom_point(alpha=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  scale_colour_viridis_d(option = "B", begin = 0, end = 0.7)+
  theme_article()+theme(legend.position = c(0.8,0.2))

ggMarginal(p_mean_sd, type="density", groupColour = T)



list_ids <- stats_all$UniqueID[which(!stats_all$ebullition.detected & stats_all$sd_ch4>1e+4)]
mydata_sel <- load_this(mylist = list_ids)
mydata_sel<- mydata_sel[order(mydata_sel$POSIX.time),]
ggplot()+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  facet_wrap(.~UniqueID, scales = "free")

CH4_flux.auto[which(CH4_flux.auto$UniqueID == "s3-da-a2-8-o-d-08:20"),]
stats_all[stats_all$UniqueID=="s3-da-a2-8-o-d-08:20",]


id = "s3-da-a2-8-o-d-08:20"
mydata_sel <- load_this(mylist = id)
print(plot_incubations(mydata_sel))


# do we see any relationship between bubbling and sd(CH4)?

list_ids_suspicious <- CH4_flux.auto$UniqueID[which(CH4_flux.auto$ebullition.flux>1000)]

stats_all_sel <- stats_all[! stats_all$UniqueID %in% list_ids_suspicious,]
stats_all_sel <- stats_all_sel[stats_all_sel$ebull_ch4>0,]

# do we see larger sd(CH4) when ebullition flux is larger?
ggplot(stats_all_sel, aes(ebull_ch4, sd_ch4))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()


ggplot(stats_all, aes(ebull_ch4, sd_ch4))+
  geom_point()+
  # scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method = "lm")+
  theme_article()

# do we see large sd(CO2) when ebullition flux is larger?
ggplot(stats_all[! stats_all$UniqueID %in% list_ids_suspicious,], aes(sd_co2, ebull_ch4))+geom_point()+theme_article()+
  scale_x_log10()+
  scale_y_log10()

ggplot(stats_all[! stats_all$UniqueID %in% list_ids_suspicious,], aes(ebull_ch4, var_dCdt_co2))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()



# do we see poor fits for CO2 when bubbling is detected?
stats_all$LM.MAE_co2 <- CO2_flux.auto$LM.MAE[match(stats_all$UniqueID, CO2_flux.auto$UniqueID)]
stats_all$LM.SE_co2 <- CO2_flux.auto$LM.SE[match(stats_all$UniqueID, CO2_flux.auto$UniqueID)]
stats_all$LM.r2_co2 <- CO2_flux.auto$LM.r2[match(stats_all$UniqueID, CO2_flux.auto$UniqueID)]



stats_all_joined <- stats_all %>%
  left_join(CO2_flux.auto, by = "UniqueID")



ggplot(stats_all_joined[stats_all_joined$ebull_ch4>0,], aes(ebull_ch4, var_co2))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()















# ---- CO2 fluxes LM vs HM ----


table(CO2_flux.man$model)

prop_best_HM <- table(CO2_flux.man$model)[1]/sum(table(CO2_flux.man$model))*100

CO2_flux.man$diff_abs_LM_HM <- (CO2_flux.man$LM.flux - CO2_flux.man$HM.flux)
CO2_flux.man$diff_rel_LM_HM <- CO2_flux.man$diff_abs_LM_HM/CO2_flux.man$LM.flux*100

ggplot(data = CO2_flux.man)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(LM.flux, HM.flux))+
  # scale_x_log10()+
  # scale_y_log10()+
  xlab("LM automated CO2 flux [mmol/m2/s]")+
  ylab("HM CO2 flux [mmol/m2/s]")+
  theme_bw()

prop_below10 <- length(which(CO2_flux.man$diff_rel_LM_HM<10))/length(CO2_flux.man$UniqueID)*100

message(paste("For both CO2 and CH4, better fit with HM for ",round(prop_best_HM),"% but difference with LM < 10% for ", round(prop_below10),"%"))



# negative fluxes

# examples of very weird choices from dear experts
list_ids <- CO2_flux.man$UniqueID[which(CO2_flux.man$HM.flux< -1)]

mydata_sel <- load_this(mylist = list_ids)
mydata_sel<- mydata_sel[order(mydata_sel$POSIX.time),]
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                        ymin = -Inf, ymax = Inf, fill = username), 
            # fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CO2dry_ppm, group = UniqueID))+
  theme_article()+
  facet_grid(UniqueID~username, scales = "free")






# ---- CO2 fluxes expert vs "blind" ----

sprd_CO2 <- CO2[,c("variable","UniqueID","flux_method","LM.flux","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(LM.flux))
#list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_co2<100)]
#ind_rmv <- match(list_rm, sprd_CO2$UniqueID)
#sprd_CO2 <- sprd_CO2[-ind_rmv,]

#sprd_CO2 <- sprd_CO2[sprd_CO2$Expert> -0.4,]

sprd_CO2$relative_diff <- abs((sprd_CO2$Expert-sprd_CO2$Blind)/sprd_CO2$Expert)*100

df_exceedance_CO2 <- get_df_exceedance(sprd_CO2$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CO2$t-10))


p_xy <- ggplot(data = sprd_CO2)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(Blind, Expert))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Blind automated CO2 flux [mmol/m2/s]")+
  ylab("Expert CO2 flux [mmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CO2 flux",", difference < 10% for ",round(df_exceedance_CO2$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CO2, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CO2$p[ind_closest_10], yend=df_exceedance_CO2$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CO2$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

p_exceed_co2 <- ggarrange(p_xy, p_exceed, nrow = 1)


ggsave(plot = p_exceed_co2, filename = "overview_expert_vs_blind_co2.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)



sprd_CO2$UniqueID[which(sprd_CO2$Blind>0.8)]

sprd_CO2$UniqueID[which(sprd_CO2$Blind < -.4)]



sprd_CO2$absolute_diff <- abs(sprd_CO2$Expert - sprd_CO2$Blind)

table_all_co2 <- table_all[table_all$variable=="CO2",]
sprd_CO2$LM.RMSE <- table_all_co2$LM.RMSE[match(sprd_CO2$UniqueID, table_all_co2$UniqueID)]


p_diff_LM.RMSE <- ggplot(sprd_CO2, aes(LM.RMSE, absolute_diff))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("RMSE Linear model [mmol/m2/s]")+
  ylab("abs(Expert - Automated) [mmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")


ggsave(plot = p_diff_LM.RMSE, filename = "overview_expert_vs_blind_co2_LM_RMSE.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)






# ---- CH4 diffusion expert vs "dydt" ----

sprd_CH4 <- table_ebull[,c("UniqueID","flux_method","diffusion","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(diffusion))
# list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_CH4<100)]
# ind_rmv <- match(list_rm, sprd_CH4$UniqueID)
# sprd_CH4 <- sprd_CH4[-ind_rmv,]

# sprd_CH4 <- sprd_CH4[sprd_CH4$Expert> -0.4,]

sprd_CH4$relative_diff <- abs((sprd_CH4$Expert-sprd_CH4$dydt)/sprd_CH4$Expert)*100

df_exceedance_CH4 <- get_df_exceedance(sprd_CH4$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))


p_xy <- ggplot(data = sprd_CH4)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(dydt, Expert))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("dydt automated CH4 flux [mmol/m2/s]")+
  ylab("Expert CH4 flux [mmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CH4 flux",", difference < 10% for ",round(df_exceedance_CH4$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

p_exceed_CH4 <- ggarrange(p_xy, p_exceed, nrow = 1)


ggsave(plot = p_exceed_CH4, filename = "overview_expert_vs_dydt_CH4.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)



sprd_CH4$UniqueID[which(sprd_CH4$dydt>0.8)]

sprd_CH4$UniqueID[which(sprd_CH4$dydt < -.4)]



sprd_CH4$absolute_diff <- abs(sprd_CH4$Expert - sprd_CH4$dydt)

table_all_CH4 <- table_all[table_all$variable=="CH4",]
sprd_CH4$LM.RMSE <- table_all_CH4$LM.RMSE[match(sprd_CH4$UniqueID, table_all_CH4$UniqueID)]


p_diff_LM.RMSE <- ggplot(sprd_CH4, aes(LM.RMSE, absolute_diff))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("RMSE Linear model [mmol/m2/s]")+
  ylab("abs(Expert - Automated) [mmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")


ggsave(plot = p_diff_LM.RMSE, filename = "overview_expert_vs_dydt_CH4_LM_RMSE.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)




# ---- CH4 ebullition expert vs "dydt" ----

sprd_ebull <- table_ebull[,c("UniqueID","flux_method","ebullition","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(ebullition))
list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_ch4<10)]

ind_rmv <- match(list_rm, sprd_ebull$UniqueID)
sprd_ebull <- sprd_ebull[-ind_rmv,]



sprd_ebull$relative_diff <- abs((sprd_ebull$Expert-sprd_ebull$dydt)/sprd_ebull$Expert)*100

df_exceedance_CH4 <- get_df_exceedance(sprd_ebull$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))


p_xy <- ggplot(data = sprd_ebull)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(dydt, Expert))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("Automated CH4 ebullition [nmol/m2/s]")+
  ylab("Expert CH4 ebullition [nmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CH4 ebullition",", difference < 10% for ",round(df_exceedance_CH4$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

p_exceed_ch4_ebull <- ggarrange(p_xy, p_exceed, nrow = 1)

ggsave(plot = p_exceed_ch4_ebull, filename = "overview_expert_vs_blind_ch4_ebull.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)




# ---- CH4 total expert vs "dydt" ----

sprd_total <- table_ebull[,c("UniqueID","flux_method","total_estimated","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(total_estimated))
list_rm <- table_draws$UniqueID[which(table_draws$duration_expert_ch4<10)]

ind_rmv <- match(list_rm, sprd_total$UniqueID)
sprd_total <- sprd_total[-ind_rmv,]
sprd_total <- sprd_total[which(sprd_total$dydt<1500),]


sprd_total$relative_diff <- abs((sprd_total$Expert-sprd_total$dydt)/sprd_total$Expert)*100

df_exceedance_CH4 <- get_df_exceedance(sprd_total$relative_diff)
ind_closest_10 <- which.min(abs(df_exceedance_CH4$t-10))

p_xy <- ggplot(data = sprd_total)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(dydt, Expert))+
  xlab("Automated CH4 total flux [nmol/m2/s]")+
  ylab("Expert CH4 total flux [nmol/m2/s]")+
  theme_bw()+ggtitle(paste0("CH4 total flux",", difference < 10% for ",round(df_exceedance_CH4$p[ind_closest_10]),"% of the measurements"))

p_exceed <- ggplot(df_exceedance_CH4, aes(t, p))+geom_path()+geom_point()+
  geom_hline(yintercept = 0)+
  theme_bw()+
  geom_segment(aes(x=-0,xend=10,
                   y=df_exceedance_CH4$p[ind_closest_10], yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  geom_segment(aes(x=10,xend=10,
                   y=-Inf, yend=df_exceedance_CH4$p[ind_closest_10]), color ="red")+
  # scale_x_log10()+
  scale_x_continuous(transform = "log10", breaks = c(1,10,100,1000))+
  ylab("Proportion of timeseries [%]")+
  xlab("Relative difference [% of Expert flux]")

p_exceed_ch4_tot <- ggarrange(p_xy, p_exceed, nrow = 1)

ggsave(plot = p_exceed_ch4_tot, filename = "overview_expert_vs_blind_ch4.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)




sprd_total$UniqueID[which(sprd_total$Expert< (-5))]




sprd_total$absolute_diff <- abs(sprd_total$Expert - sprd_total$dydt)

sprd_total$ebullition_expert <- sprd_ebull$Expert[match(sprd_total$UniqueID, sprd_ebull$UniqueID)]
sprd_total$diffusion_expert <- sprd_diffus$Expert[match(sprd_total$UniqueID, sprd_diffus$UniqueID)]
sprd_total$LM.RMSE <- table_ebull$LM.RMSE[match(sprd_total$UniqueID, table_ebull$UniqueID)]
sprd_total$LM.r2 <- table_ebull$LM.r2[match(sprd_total$UniqueID, table_ebull$UniqueID)]



p_diff_LM.RMSE <- ggplot(sprd_total[sprd_total$Expert>0,], aes(LM.RMSE, absolute_diff, size = ebullition_expert))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("RMSE Linear model [nmol/m2/s]")+
  ylab("abs(Expert - Automated) [nmol/m2/s]")+
  scale_size("Ebullition [nmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")

ggsave(plot = p_diff_LM.RMSE, filename = "overview_expert_vs_blind_LM_RMSE.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)





ggplot(sprd_total[sprd_total$Expert>0,], aes(Expert, diffusion_expert, size = ebullition_expert))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("Expert total [nmol/m2/s]")+
  ylab("Expert diffusive [nmol/m2/s]")+
  scale_size("Expert ebullition [nmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")



# 
# 
# p_prop <- ggplot(sprd_ebull[order(sprd_ebull$Expert),], 
#                  aes(seq_along(sprd_ebull$Expert)/dim(sprd_ebull)[1]*100,Expert))+
#   geom_hline(yintercept = 100)+
#   geom_hline(yintercept = 50, alpha=0.2)+
#   # geom_jitter(alpha=0.5)+
#   geom_point(alpha=0.5, width=0.3, aes(colour = sprd_ebull$Expert[order(sprd_ebull$Expert)] > 50))+
#   theme_article()+
#   # ylab("CH4 flux nmol/m2/s")+
#   scale_fill_viridis_d(begin = 0.2, end = 0.9)+
#   scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
#   ylim(c(0,120))+
#   # scale_y_log10()+
#   # facet_wrap(strata~.)+
#   theme(legend.position = 'none')+
#   xlab("Proportion of timeseries [%]")+
#   ylab("Contribution of bubbling to total CH4 flux [%]")+
#   ggtitle(paste0("Bubbling > 50% of F_total for ",round(n_50/dim(sprd_ebull)[1]*100),"% of the measurements"))





table_ch4_ebull <- sprd_total[!is.na(sprd_total$ebullition_expert),]
table_ch4_ebull$p_ebullition <- table_ch4_ebull$ebullition_expert/table_ch4_ebull$Expert *100

n_50 <- length(which(table_ch4_ebull$p_ebullition > 50))

p_prop <- ggplot(table_ch4_ebull[order(table_ch4_ebull$p_ebullition),], 
                 aes(seq_along(table_ch4_ebull$p_ebullition)/dim(table_ch4_ebull)[1]*100,p_ebullition))+
  geom_hline(yintercept = 100)+
  geom_hline(yintercept = 50, alpha=0.2)+
  # geom_jitter(alpha=0.5)+
  geom_point(alpha=0.5, width=0.3, aes(colour = table_ch4_ebull$p_ebullition[order(table_ch4_ebull$p_ebullition)] > 50))+
  theme_article()+
  # ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  ylim(c(0,120))+
  # scale_y_log10()+
  # facet_wrap(strata~.)+
  theme(legend.position = 'right')+
  xlab("Proportion of timeseries [%]")+
  ylab("Contribution of bubbling to total CH4 flux [%]")


p_ebull <- ggplot(table_ch4_ebull, 
                  aes(" ", Expert))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2, aes(colour = table_ch4_ebull$p_ebullition > 50))+
  geom_boxplot(alpha=0.8, width=0.2)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  theme_article()+
  scale_y_log10()+
  ylab("CH4 total flux [nmol/m2/s]")+
  theme(legend.position = 'none')+
  xlab("")#+ggtitle(paste0("Bubbling > 50% of total flux for ",round(n_50/dim(table_ch4_ebull)[1]*100),"% of the measurements"))


p_contrib_ebull_expert <- ggarrange(p_ebull, p_prop, nrow = 1)

ggsave(plot = p_contrib_ebull_expert, filename = "overview_expert_vs_blind_contrib_ebull.jpeg", path = plots_path, 
       width = 8, height = 4, dpi = 300, units = 'in', scale = 1.2)






sprd_total$ebullition_blind <- sprd_ebull$dydt[match(sprd_total$UniqueID, sprd_ebull$UniqueID)]
sprd_total$diffusion_blind <- sprd_diffus$dydt[match(sprd_total$UniqueID, sprd_diffus$UniqueID)]


ggplot(sprd_total[sprd_total$dydt>0,], aes(dydt, diffusion_blind, size = ebullition_blind))+geom_point()+
  # geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()+
  xlab("Expert total [nmol/m2/s]")+
  ylab("Expert diffusive [nmol/m2/s]")+
  scale_size("Expert ebullition [nmol/m2/s]")+
  # ggtitle("Poor linear fits correspond to largest absolute differences")+
  theme(legend.position = "top")




table_ch4_ebull <- sprd_total[!is.na(sprd_total$ebullition_blind),]
table_ch4_ebull$p_ebullition <- table_ch4_ebull$ebullition_blind/table_ch4_ebull$dydt *100

n_50 <- length(which(table_ch4_ebull$p_ebullition > 50))

p_prop <- ggplot(table_ch4_ebull[order(table_ch4_ebull$p_ebullition),], 
                 aes(seq_along(table_ch4_ebull$p_ebullition)/dim(table_ch4_ebull)[1]*100,p_ebullition))+
  geom_hline(yintercept = 100)+
  geom_hline(yintercept = 50, alpha=0.2)+
  # geom_jitter(alpha=0.5)+
  geom_point(alpha=0.5, width=0.3, aes(colour = table_ch4_ebull$p_ebullition[order(table_ch4_ebull$p_ebullition)] > 50))+
  theme_article()+
  # ylab("CH4 flux nmol/m2/s")+
  scale_fill_viridis_d(begin = 0.2, end = 0.9)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  ylim(c(0,120))+
  # scale_y_log10()+
  # facet_wrap(strata~.)+
  theme(legend.position = 'right')+
  xlab("Proportion of timeseries [%]")+
  ylab("Contribution of bubbling to total CH4 flux [%]")


p_ebull <- ggplot(table_ch4_ebull, 
                  aes(" ", Expert))+
  # geom_hline(yintercept = 100)+
  geom_jitter(alpha=0.5, width = 0.2, aes(colour = table_ch4_ebull$p_ebullition > 50))+
  geom_boxplot(alpha=0.8, width=0.2)+
  scale_colour_viridis_d("over 50%",begin = 0.2, end = 0.9, option = "A", direction = -1)+
  theme_article()+
  scale_y_log10()+
  ylab("CH4 total flux [nmol/m2/s]")+
  theme(legend.position = 'none')+
  xlab("")+
  ggtitle(paste0("Bubbling > 50% of total flux for ",round(n_50/dim(table_ch4_ebull)[1]*100),"% of the measurements"))


p_contrib_ebull_blind <- ggarrange(p_ebull, p_prop, nrow = 1)









# ---- Stats in timeseries cuts ----


ggplot(table_draws, aes(draw))+geom_density()

table_draws$diff_t_start_co2 <- table_draws$start.time_expert_co2-table_draws$start.time_auto
table_draws$diff_t_start_ch4 <- table_draws$start.time_expert_ch4-table_draws$start.time_auto
table_draws$diff_t_end_co2 <- table_draws$end.time_expert_co2-table_draws$end.time_auto
table_draws$diff_t_end_ch4 <- table_draws$end.time_expert_ch4-table_draws$end.time_auto


ggplot(table_draws, aes(diff_t_start_co2))+geom_density()

table_draws$duration_expert_co2 <- table_draws$end.time_expert_co2 - table_draws$start.time_expert_co2

ggplot(table_draws, aes(reorder(draw, duration_expert_co2, FUN=min), duration_expert_co2))+geom_point()+theme_article()

table_draws$UniqueID[which(table_draws$duration_expert_co2 < 100)]


# ---- Stats in fluxes differences ----

table_results_sprd <- table_all[,c("variable","UniqueID","flux_method","best.flux","timestamp_processing","username")] %>%
  pivot_wider(names_from = flux_method, values_from = c(best.flux))

ggplot(data = table_results_sprd)+
  # geom_abline(slope = 0,intercept = 0, color = 'lightgrey')+
  # geom_segment(data = data.frame(UniqueID = CO2_flux_res_auto$UniqueID,
  #                                meth1 = CO2_flux_res_auto$best.flux,
  #                                meth2 = CO2_flux_res_manID$best.flux), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_violin(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5)+
  geom_jitter(aes(variable,abs((Blind-Expert)/Expert)*100), alpha = 0.5, width=0.2)+
  ylab("CO2 flux relative difference [%]")+
  theme_bw()+
  # scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






plot(sort(table(paste0(table_results_sprd$variable,table_results_sprd$UniqueID))))

plot(sort(table(table_results_sprd$username)))


# ---- Ebullition ----

ggplot(data = table_ebull[which(table_ebull$flux_method=="Expert" & table_ebull$diffusion>0 & table_ebull$ebullition>0),])+ # [table_ebull$flux_method=="Expert",]
  # geom_abline(slope = 1,intercept = 0, color = 'grey')+
  # geom_segment(data = data.frame(UniqueID = CH4_res_meth1$UniqueID,
  #                                meth1 = CH4_res_meth1$ebullition,
  # meth2 = CH4_res_meth2$ebullition), aes(x=UniqueID, xend=UniqueID, y = meth1, yend = meth2), linewidth=1, alpha = 0.5)+
  geom_point(aes(diffusion, ebullition,
                 colour = log10(total_estimated)), size=4, alpha = 0.5)+
  # ylab("ebullition component [nmol/m2/s]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()+ 
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9, option = "F", direction = -1)


ggplot(table_ebull[which(table_ebull$flux_method=="Expert" & table_ebull$diffusion>0 & table_ebull$ebullition>0),])+
  geom_density(aes(ebullition/diffusion, fill="ebullition/diffusion"))+
  # geom_boxplot(aes(ebullition, fill="ebullition"))+
  theme_article()+coord_flip()


table_ebull$UniqueID[which(table_ebull$total_estimated>1000)]

















# ---- Co2 vs CH4 ----

# draw <- sample(seq_along(auxfile$subsite), nb_draw)
draw <- seq_along(auxfile$subsite)

# draw <- which(auxfile$UniqueID == "s1-cu-a2-3-o-d-06:59")

table_draw <- data.frame(username = username,
                         draw = draw,
                         subsite = auxfile$subsite[draw],
                         UniqueID = auxfile$UniqueID[draw])
auxfile <- auxfile[draw,]
auxfile$username <- username

mydata_all <- NULL
for(k in seq_along(auxfile$UniqueID)){
  message("Loading data for ",auxfile$UniqueID[k])
  gas <- unique(auxfile$gas_analiser[k])
  
  setwd(RData_path)
  if(gas== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gas
  }
  load(file = paste0(auxfile$subsite[k],"_",gs_suffix,".RData"))
  mydata <- mydata[,c("POSIX.time", 
                      "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                      "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
  
  mydata <- mydata[which(mydata$POSIX.time>=auxfile$start.time[k] & 
                           mydata$POSIX.time<=auxfile$start.time[k]+auxfile$duration[k]),]
  mydata$UniqueID <- auxfile$UniqueID[k]
  mydata$Etime <- mydata$POSIX.time - min(mydata$POSIX.time)
  
  table_draw$corr_co2_ch4[k] <- cor(mydata$CO2dry_ppm, mydata$CH4dry_ppb)
  mod_lm <- lm(data = mydata, formula = CH4dry_ppb~CO2dry_ppm)
  table_draw$slope_co2_ch4[k] <- coefficients(mod_lm)[2]
  table_draw$r2_co2_ch4[k] <- summary(mod_lm)$adj.r.squared
  
  mydata_all <- rbind(mydata_all, mydata)
  rm(mydata)
}



# ggplot(mydata_all, aes(Etime, CH4dry_ppb))+geom_point()+geom_smooth(method = 'lm')+facet_wrap(UniqueID~.)+theme_article()
# ggplot(mydata_all, aes(CO2dry_ppm, CH4dry_ppb))+geom_point()+geom_smooth(method = 'lm')+facet_wrap(UniqueID~.)+theme_article()



ggplot(table_draw, aes(corr_co2_ch4))+geom_density()

ggplot(table_draw, aes(slope_co2_ch4, corr_co2_ch4))+geom_point()


table_draw$UniqueID[which(table_draw$slope_co2_ch4>5000)]


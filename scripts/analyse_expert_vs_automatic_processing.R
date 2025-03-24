
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


# ---- Directories ----

# You have to make sure this is pointing to the write folder on your local machine
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 


datapath <- paste0(dropbox_root,"/GHG/RAW data")
loggerspath <- paste0(datapath,"/RAW Data Logger")
RData_path <- paste0("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/RESTORE4Cs/data/Harmonized_GHG")
# plots_path <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/meetings_presentations/2024_PPNW_Girona"
plots_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/plots/")
results_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/results/")


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

source("C:/Projects/myGit/aquaGHG/R/get_dxdy.R")
source("C:/Projects/myGit/aquaGHG/R/get_dCdt_density.R")
source("C:/Projects/myGit/aquaGHG/R/find_first_linear_chunk.R")
source("C:/Projects/myGit/aquaGHG/R/find_bubbles.R")

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
      mydata$Etime <- mydata$POSIX.time - min(mydata$POSIX.time)
      
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





# ---- Loading data ----


# loading auxfile
setwd(dirname(results_path))
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

# ---- Retrieve timestamp processing for table_draw and making unique IDs ----

# this was done in recalculate_fluxes_with_aquaGHG.R



# ---- overview database ----

message("contributions by experts")
table(table_draws$username)

tab_users <- as.data.frame(table(table_draws$username))
tab_users$userID <- paste0("user",seq(1,dim(tab_users)[1]))

tab_users$userID <- factor(x = tab_users$userID, levels = tab_users$userID[order(tab_users$Freq, decreasing = T)])

table_draws$userID <- tab_users$userID[match(table_draws$username, tab_users$Var1)]

p_users <- ggplot(tab_users, aes(userID, Freq))+geom_point()+theme_article()+
  ggtitle(paste0("n = ",dim(tab_users)[1]))+
  xlab("User ID")+ylab("Counts")+coord_flip()

p_users



n <- length(unique(table_draws$UniqueID))
table_n <- data.frame(id = names(table(table_draws$UniqueID)),
                      n = as.numeric(table(table_draws$UniqueID)))


p_overview <- ggplot(table_n, aes(n))+geom_histogram()+theme_article()+xlab("Repetitions")+ylab("Counts")+
  ggtitle(paste0("n = ",n))+
  scale_x_continuous(breaks = breaks_pretty())

p_overview




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

p_overview_dist <- ggpubr::ggarrange(p_users, p_overview, p_dist_co2, p_dist_ch4, ncol = 2, nrow = 2)
p_overview_dist

ggsave(plot = p_overview_dist, filename = "overview_database_ghg_expert.jpeg", path = plots_path, 
       width = 6, height = 6, dpi = 300, units = 'in', scale = 1.0)



stats_all$flux_co2 <- CO2_flux.auto$best.flux[match(stats_all$UniqueID, CO2_flux.auto$UniqueID)]
stats_all$flux_ch4 <- CH4_flux.auto$best.flux[match(stats_all$UniqueID, CH4_flux.auto$UniqueID)]
stats_all$ebull_ch4 <- CH4_flux.auto$ebullition.flux[match(stats_all$UniqueID, CH4_flux.auto$UniqueID)]
stats_all$diffus_ch4 <- CH4_diff_flux.man$best.flux[match(stats_all$UniqueID, CH4_diff_flux.man$UniqueID)]



p_co2_ch4_fluxes <- ggplot(stats_all, aes(flux_co2, flux_ch4))+geom_point(alpha=0.5)+
  # xlab("Average CO2 [ppm]")+
  # ylab("Standard deviation CO2 [ppm]")+
  scale_x_log10()+
  scale_y_log10()+
  theme_article()
ggMarginal(p_co2_ch4_fluxes, type="density")


ggplot(stats_all, aes(as.factor(n_bubbles), ebull_ch4))+geom_boxplot(alpha=0.5)+
  # xlab("Average CO2 [ppm]")+
  # ylab("Standard deviation CO2 [ppm]")+
  # scale_x_log10()+
  scale_y_log10()+
  theme_article()





# ---- Flagged anomalous ----
# threshold_anomalous <- 30 # seconds or n observations
# 
# ggplot(table_draws, aes(duration_expert_co2 ))+
#   geom_histogram(alpha=0.5)+theme_article()+geom_vline(xintercept = threshold_anomalous)
# ggplot(table_draws, aes(duration_expert_ch4))+
#   geom_histogram(alpha=0.5)+theme_article()+geom_vline(xintercept = threshold_anomalous)
# 
# # table_all$uniqAssessID <- paste(table_all$timestamp_processing, table_all$UniqueID, sep="/")
# # table_ebull$uniqAssessID <- paste(table_ebull$timestamp_processing, table_ebull$UniqueID, sep="/")
# 
# tmpstmp_flag_co2 <- unique(table_draws$uniqAssessID[which(table_draws$duration_expert_co2 < threshold_anomalous)])
# tmpstmp_flag_ch4 <- unique(table_draws$uniqAssessID[which(table_draws$duration_expert_ch4 < threshold_anomalous)])
# 
# incubs_flagged_co2 <- unique(table_all$UniqueID[which(table_all$flux_method=="Expert" & table_all$variable=="CO2" & table_all$nb.obs<threshold_anomalous)])
# incubs_flagged_ch4 <- unique(table_ebull$UniqueID[which(table_ebull$nb.obs<threshold_anomalous)])[-1]
# 
# tmp <- rbind(data.frame(variable = "co2",
#                         ID = incubs_flagged_co2),
#              data.frame(variable = "ch4",
#                         ID = incubs_flagged_ch4))
# 
# tmp_flag_both <- tmp[duplicated(tmp$ID), ]
# tmp_flag_NOTboth <- tmp[ ! tmp$ID %in% tmp_flag_both$ID, ]
#  
# # what is special about these flagged incubations?
# stats_all$co2flagged <- "valid"
# stats_all$co2flagged[which(stats_all$UniqueID%in%incubs_flagged_co2)] <- "flagged"
# stats_all$ch4flagged <- "valid"
# stats_all$ch4flagged[which(stats_all$UniqueID%in%incubs_flagged_ch4)] <- "flagged"
# 
# 
# p_mean_sd_co2_f <- ggplot(stats_all, aes(mean_co2, sd_co2, colour = co2flagged))+geom_point(alpha=0.5)+
#   xlab("Average CO2 [ppm]")+
#   ylab("Standard deviation CO2 [ppm]")+
#   scale_x_log10()+
#   scale_y_log10()+
#   theme_article()+
#   scale_colour_viridis_d(end = .9, option = "A")+
#   theme(legend.title = element_blank())+ggtitle(paste0("n flagged = ", length(incubs_flagged_co2)))
# 
# p_mean_sd_ch4_f <- ggplot(stats_all, aes(mean_ch4, sd_ch4, colour = ch4flagged))+geom_point(alpha=0.5)+
#   xlab("Average CH4 [ppb]")+
#   ylab("Standard deviation CH4 [ppb]")+
#   scale_x_log10()+
#   scale_y_log10()+
#   theme_article()+
#   scale_colour_viridis_d(end = .9, option = "A")+
#   theme(legend.title = element_blank())+ggtitle(paste0("n flagged = ", length(incubs_flagged_ch4)))
# 
# # p_dist_co2 <- ggMarginal(p_mean_sd_co2, type="density", groupFill = T)
# p_flagged <- ggpubr::ggarrange(p_mean_sd_co2_f, p_mean_sd_ch4_f, common.legend = T, legend = "right")
# 
# ggsave(plot = p_flagged, filename = "flagged.jpeg", path = plots_path, 
#        width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)



# 
# names(stats_all)
# 
# 
# 
# ggplot(stats_all, aes(ch4flagged, mean_ch4, colour = ch4flagged))+geom_jitter(alpha=0.5)+
#   # scale_y_log10()+
#   scale_colour_viridis_d(begin = 0.1, end = 0.7, option = "A", direction = -1)+
#   theme_article()
# 
# ggplot(stats_all, aes(mean_ch4, sd_ch4, colour = ch4flagged))+geom_jitter(alpha=0.5)+
#   # scale_y_log10()+
#   scale_colour_viridis_d(begin = 0.1, end = 0.7, option = "A", direction = -1)+
#   theme_article()


# # evidence the ones with clear cuts
# 
# table_draws$duration_auto <- table_draws$end.time_auto - table_draws$start.time_auto
# table_draws$duration_ratio <- table_draws$duration_expert_co2/table_draws$duration_auto
# 
# ggplot(table_draws, aes(duration_ratio))+geom_histogram()+theme_article()+geom_vline(xintercept = .75)
# 
# list_much_shorter <- table_draws$UniqueID[which(table_draws$duration_ratio>.10 & table_draws$duration_ratio<.25)]
# 
# mydata_all <- load_this(mylist = list_flag_co2)
# 
# 
# ggplot(mydata_all, aes(Etime, CO2dry_ppm, colour = flag_score_num, group_by = UniqueID))+geom_path()+theme_article()+
#   facet_wrap(UniqueID~., scales = "free_y")+scale_colour_viridis_c(option = "C", direction = -1)




# ---- Differences between experts CO2 ----

# first we remove from the tables the incubations flagged as anomalous
# table_draws <- table_draws[ ! table_draws$uniqAssessID %in% tmpstmp_flag_co2, ]
# table_all_clean <- table_all[ ! table_all$uniqAssessID %in% tmpstmp_flag_co2, ]

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


names(df_multiple_users_co2)


p1 <- ggplot(df_multiple_users_co2, aes(flux_blind, flux_users_mean, colour = n))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  geom_errorbar(aes(ymin=flux_users_mean-flux_users_sd, 
                    ymax=flux_users_mean+flux_users_sd))+
  xlab("Blind LM (CO2)")+
  ylab("Expert LM (CO2)")+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()


p2 <- ggplot(df_multiple_users_co2, aes(abs_diff, abs(flux_users_CV), colour = n))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("abs((Blind - mean(Expert)) [CO2 flux]")+
  ylab("abs(CV(Expert)) [d.l.]")+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()

p_expert_vs_blind_CO2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
p_expert_vs_blind_CO2
ggsave(plot = p_expert_vs_blind_CO2, filename = "expert_vs_blind_CO2.jpeg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)




# plot incubations with the highest disagreements between experts
list_ids <- df_multiple_users_co2$id[which(abs(df_multiple_users_co2$flux_users_CV)>1)]

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_co2-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_co2-(table_draws_sel$start.time_auto)


p_disagree_CO2 <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                              ymin = -Inf, ymax = Inf), fill="grey50", alpha=0.2)+
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




# ---- Differences between experts CH4 ----

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
                                              
                                              total_blind = CH4_flux.auto$total.flux[i_flux_blind],
                                              total_users_mean = mean(tab_users$total.flux),
                                              total_users_sd = sd(tab_users$total.flux),
                                              total_users_CV = sd(tab_users$total.flux)/mean(tab_users$total.flux),
                                              
                                              ebull_blind = CH4_flux.auto$ebullition.flux[i_flux_blind],
                                              ebull_users_mean = mean(tab_users$ebullition.flux),
                                              ebull_users_sd = sd(tab_users$ebullition.flux),
                                              ebull_users_CV = sd(tab_users$ebullition.flux)/mean(tab_users$ebullition.flux),
                                              
                                              diffusion_blind = CH4_flux.auto$ebullition.flux[i_flux_blind],
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

names(df_multiple_users_CH4)

# Total flux

df_multiple_users_CH4$ebull_contrib <- df_multiple_users_CH4$ebull_blind/df_multiple_users_CH4$total_blind


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




p1 <- ggplot(df_multiple_users_CH4[which(df_multiple_users_CH4$total_users_mean > 0),],
             aes(total_blind, total_users_mean, colour = n))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(size = ebull_contrib))+
  geom_errorbar(aes(ymin=total_users_mean-total_users_sd, 
                    ymax=total_users_mean+total_users_sd))+
  xlab("Blind (ch4 total)")+
  ylab("Expert (ch4 total)")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()


p2 <- ggplot(df_multiple_users_CH4[which(df_multiple_users_CH4$total_users_mean > 0),], 
             aes(total_abs_diff, abs(total_users_CV), colour = n, size = ebull_contrib))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("abs((Blind - mean(Expert)) [ch4 flux]")+
  ylab("abs(CV(Expert)) [d.l.]")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()

p_expert_vs_blind_CH4tot <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
ggsave(plot = p_expert_vs_blind_CH4tot, filename = "expert_vs_blind_CH4tot.jpeg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)


# plot incubations with the highest disagreements between experts
list_ids <- df_multiple_users_CH4$id[which(abs(df_multiple_users_CH4$total_users_CV)>0.9)]

# list_ids <- df_multiple_users_CH4$id[which(abs(df_multiple_users_CH4$ebull_blind)>1500)]

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

p_disagree_CH4tot <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                              ymin = -Inf, ymax = Inf), 
            fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  facet_wrap(UniqueID~., scales = "free")

ggsave(plot = p_disagree_CH4tot, filename = "expert_choice_CH4tot.jpeg", path = plots_path, 
       width = 12, height = 10, dpi = 300, units = 'in', scale = 1.0)


# Diffusion
p1 <- ggplot(df_multiple_users_CH4, aes(diffusion_blind, diffusion_users_mean, colour = n))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(size = ebull_contrib))+
  geom_errorbar(aes(ymin=diffusion_users_mean-diffusion_users_sd, 
                    ymax=diffusion_users_mean+diffusion_users_sd))+
  xlab("Blind (ch4 diffusion)")+
  ylab("Expert (ch4 diffusion)")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()


p2 <- ggplot(df_multiple_users_CH4, 
             aes(diffusion_abs_diff, abs(diffusion_users_CV), colour = n, size = ebull_contrib))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("abs((Blind - mean(Expert)) [ch4 flux]")+
  ylab("abs(CV(Expert)) [d.l.]")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()

p_expert_vs_blind_CH4diff <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
ggsave(plot = p_expert_vs_blind_CH4diff, filename = "expert_vs_blind_CH4diff.jpeg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)


# plot incubations with the highest disagreements between experts
list_ids <- df_multiple_users_CH4$id[which(abs(df_multiple_users_CH4$diffusion_users_CV)>1)]

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)

table_draws_sel[order(table_draws_sel$UniqueID),]

p_disagree_CH4diff <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                              ymin = -Inf, ymax = Inf), 
            fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group_by = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  facet_wrap(UniqueID~., scales = "free")
ggsave(plot = p_disagree_CH4diff, filename = "expert_choice_CH4diff.jpeg", path = plots_path, 
       width = 12, height = 10, dpi = 300, units = 'in', scale = 1.0)




# Ebullition
ggplot(df_multiple_users_CH4, aes(ebull_contrib))+geom_density()+theme_article()


df_ebull <- df_multiple_users_CH4[df_multiple_users_CH4$ebull_contrib>0.05,]
p1 <- ggplot(df_ebull, aes(ebull_blind, ebull_users_mean, colour = n))+
  scale_x_log10()+
  scale_y_log10()+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(size = ebull_contrib))+
  geom_errorbar(aes(ymin=ebull_users_mean-ebull_users_sd, 
                    ymax=ebull_users_mean+ebull_users_sd))+
  xlab("Blind (ch4 ebullition)")+
  ylab("Expert (ch4 ebullition)")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()

p2 <- ggplot(df_ebull, aes(ebull_abs_diff, abs(ebull_users_CV), colour = n, size = ebull_contrib))+geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  xlab("abs((Blind - mean(Expert)) [ch4 flux]")+
  ylab("abs(CV(Expert)) [d.l.]")+
  scale_size(range = c(0.1,3))+
  scale_colour_viridis_c(direction = -1, option = "A", begin = 0.1, end = 0.9)+
  theme_article()


p_expert_vs_blind_CH4ebull <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "right")
ggsave(plot = p_expert_vs_blind_CH4ebull, filename = "expert_vs_blind_CH4ebull.jpeg", path = plots_path, 
       width = 6, height = 3.2, dpi = 300, units = 'in', scale = 1.0)



# plot incubations with the highest ebullition contribution

list_ids <- df_ebull$id[which(df_ebull$ebull_contrib>.9)]

mydata_sel <- load_this(mylist = list_ids)
table_draws_sel <- table_draws[which(table_draws$UniqueID%in%list_ids),]
table_draws_sel$Etime_start <- table_draws_sel$start.time_expert_ch4-(table_draws_sel$start.time_auto)
table_draws_sel$Etime_stop <- table_draws_sel$end.time_expert_ch4-(table_draws_sel$start.time_auto)


p_disagree_CH4ebull <- ggplot()+
  geom_rect(data = table_draws_sel, aes(xmin = Etime_start, xmax = Etime_stop, 
                                              ymin = -Inf, ymax = Inf), 
            fill="grey50", 
            alpha=0.2)+
  geom_path(data = mydata_sel, aes(as.numeric(Etime), CH4dry_ppb, group = UniqueID))+
  theme_article()+
  xlab("Elapsed time [secs]")+
  facet_wrap(UniqueID~., scales = "free")
ggsave(plot = p_disagree_CH4ebull, filename = "expert_choice_CH4ebull.jpeg", path = plots_path, 
       width = 12, height = 10, dpi = 300, units = 'in', scale = 1.0)




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





























# ---- CO2 fluxes LM vs HM ----

CO2 <- table_all[table_all$variable=="CO2",]
# get rid of incubations flagged as anomalous
CO2 <- CO2[ ! CO2$UniqueID %in% list_flag_co2, ]


table(CO2$model)

prop_best_HM <- table(table_all$model)[1]/sum(table(table_all$model))*100

table_all$diff_abs_LM_HM <- table_all$LM.flux - table_all$HM.flux
table_all$diff_rel_LM_HM <- table_all$diff_abs_LM_HM/table_all$LM.flux*100

ggplot(data = table_all)+
  geom_abline(slope = 1,intercept = 0, color = 'grey')+
  geom_point(aes(LM.flux, HM.flux))+
  scale_x_log10()+
  scale_y_log10()+
  xlab("LM automated CO2 flux [mmol/m2/s]")+
  ylab("HM CO2 flux [mmol/m2/s]")+
  theme_bw()

prop_below10 <- length(which(table_all$diff_rel_LM_HM<10))/length(table_all$UniqueID)*100

message(paste("For both CO2 and CH4, better fit with HM for ",round(prop_best_HM),"% but difference with LM < 10% for ", round(prop_below10),"%"))


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


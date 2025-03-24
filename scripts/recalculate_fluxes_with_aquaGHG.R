# ---
# Authors: Camille Minaudo
# Project: "RESTORE4Cs"
# date: "Oct 2024"
# https://github.com/camilleminaudo/ghg-flux-expert
# ---

# --- Description of this script
# This script loads a selection of randomly chosen incubations and asks the expert 
# to manually select what looks like safe data to calculate CO2 flux. 
# For CH4, experts are asked to manually select what looks to them like diffusion 
# only (or no ebullition).


# clearing workspace and console
rm(list = ls()) # clear workspace
cat("/014") # clear console


library(lubridate)
library(tidyr)
library(dplyr)
library(pbapply)
library(ggplot2)
library(egg)
library(goFlux)
library(purrr)
library(ggnewscale)
library(stringr)
library(data.table)

files.sources = list.files(path = "C:/Projects/myGit/aquaGHG/R/", full.names = T)
for (f in files.sources){source(f)}



# ---- Directories ----
dropbox_root <- "C:/Users/Camille Minaudo/Dropbox/RESTORE4Cs - Fieldwork/Data" 

datapath <- paste0(dropbox_root,"/GHG/RAW data")
RData_path <- paste0("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/PROJECTS/RESTORE4Cs/data/Harmonized_GHG") # be careful with the data in the Dropbox. On 24/03/2025, message sent to Miguel about it
plots_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/plots/")
results_path <- paste0(dropbox_root,"/GHG/GHG_expert_vs_automated/results/")


# ---- Loading auxfile ----
setwd(dirname(results_path))
old <- read.csv(file = "auxfile.csv")
new <- read.csv(file = "auxfile_20250320.csv")

old$start.time <- as.POSIXct(old$start.time, tz = 'UTC')
new$start.time <- as.POSIXct(new$start.time, tz = 'UTC')

new$start.time_old <- old$start.time[match(new$UniqueID, old$UniqueID)]

new$diff_start <- abs(as.numeric(new$start.time_old) - as.numeric(new$start.time))
list_different <- new$UniqueID[which(new$diff_start>60)]

auxfile <- old
auxfile <- auxfile %>% mutate(start.time = as.POSIXct(round(start.time,"secs")))
auxfile$obs.length <- floor(auxfile$duration)

# example for issue with data stored in RData_path
# id = "s4-cu-p2-2-o-d-10:28"
# data = load_incubation(auxfile_i = auxfile[which(auxfile$UniqueID==id),], RData_path = RData_path)
# 
# ggplot(data, aes(Etime, CO2dry_ppm))+geom_path()+theme_article()+ggtitle(id)


# ---- Loading ghg-expert results ----
load_fs <- function(path, pattern){
  
  list_f <- list.files(path = path, pattern = pattern, full.names = T)
  list_f <- list_f[grep(pattern = "csv", x = list_f)]
  
  isF <- T
  for (f in list_f){
    df_tmp <- read.csv(f)
    if(isF){
      isF <- F
      df_out <- df_tmp
    } else {
      df_out <- rbind(df_out,df_tmp)
    }
  }
  return(df_out)
}

table_all <- load_fs(path = results_path, pattern = "BLIND_vs_EXPERT_co2_ch4_fluxes_")
table_ebull <- load_fs(path = results_path, pattern = "ch4_ebullition")
table_ebull$ebullition[which(table_ebull$ebullition>table_ebull$total_estimated)] <- table_ebull$total_estimated[which(table_ebull$ebullition>table_ebull$total_estimated)]
table_ebull$ebullition[which(table_ebull$ebullition < 0)] <- 0
table_draws <- load_fs(path = results_path, pattern = "table_draw")

draws <- data.frame(draw = seq_along(unique(table_draws$UniqueID)),
                    UniqueID = unique(table_draws$UniqueID))

table_draws$draw_corrected <- draws$draw[match(table_draws$UniqueID, draws$UniqueID)]




# clean names
table_all <- table_all %>%
  mutate(username = str_replace_all(username, fixed("Camille Minaudo"), "Camille"),
         username = str_replace_all(username, fixed("Dmorant"), "DMorant"),
         username = str_replace_all(username, fixed("Jorge Montes"), "Jorge"),
         username = str_replace_all(username, fixed("Katrin Attermeyer"), "Katrin"),
         username = str_replace_all(username, fixed("Miguel Cabrera"), "Miguel"))

table_ebull <- table_ebull %>%
  mutate(username = str_replace_all(username, fixed("Camille Minaudo"), "Camille"),
         username = str_replace_all(username, fixed("Dmorant"), "DMorant"),
         username = str_replace_all(username, fixed("Jorge Montes"), "Jorge"),
         username = str_replace_all(username, fixed("Katrin Attermeyer"), "Katrin"),
         username = str_replace_all(username, fixed("Miguel Cabrera"), "Miguel"))

table_draws <- table_draws %>%
  mutate(username = str_replace_all(username, fixed("Camille Minaudo"), "Camille"),
         username = str_replace_all(username, fixed("Dmorant"), "DMorant"),
         username = str_replace_all(username, fixed("Jorge Montes"), "Jorge"),
         username = str_replace_all(username, fixed("Katrin Attermeyer"), "Katrin"),
         username = str_replace_all(username, fixed("Miguel Cabrera"), "Miguel"))



# retrieve timestamp processing for table_draw and making unique IDs

table_draw_fixed <- NULL
for(id in unique(table_draws$UniqueID)){
  message("Retrieving info for ",id)
  
  table_draw_all_users <- table_draws[which(table_draws$UniqueID==id),]
  tab_search <- table_all[which(table_all$flux_method=="Expert" & table_all$UniqueID==id & table_all$variable=="CO2"),]
  
  # build ID
  tab_search$IDtmp <- ""
  table_draw_all_users$IDtmp <- ""
  for (user in tab_search$username){
    tab_search$IDtmp[which(tab_search$username==user)] <- paste0(user, rev(seq(1, length(which(tab_search$username==user)))))
    table_draw_all_users$IDtmp[which(table_draw_all_users$username==user)] <- paste0(user, rev(seq(1, length(which(table_draw_all_users$username==user)))))
  }
  
  table_draw_all_users$timestamp_processing <- tab_search$timestamp_processing[match(table_draw_all_users$IDtmp, tab_search$IDtmp)]
  table_draw_fixed <- rbind(table_draw_fixed, table_draw_all_users)
}

table_draw_fixed$uniqAssessID <- paste(table_draw_fixed$draw_corrected,table_draw_fixed$IDtmp, sep = "/")
table_draws <- table_draw_fixed[order(table_draw_fixed$draw_corrected),]

table_draws$diff_t_start_co2 <- table_draws$start.time_expert_co2-table_draws$start.time_auto
table_draws$diff_t_start_ch4 <- table_draws$start.time_expert_ch4-table_draws$start.time_auto
table_draws$diff_t_end_co2 <- table_draws$end.time_expert_co2-table_draws$end.time_auto
table_draws$diff_t_end_ch4 <- table_draws$end.time_expert_ch4-table_draws$end.time_auto
table_draws$duration_expert_co2 <- table_draws$end.time_expert_co2 - table_draws$start.time_expert_co2
table_draws$duration_expert_ch4 <- table_draws$end.time_expert_ch4 - table_draws$start.time_expert_ch4


table_all$uniqAssessID <- table_draw_fixed$uniqAssessID[match(paste0(table_all$username,table_all$UniqueID,table_all$timestamp_processing), 
                                                              paste0(table_draw_fixed$username,table_draw_fixed$UniqueID,table_draw_fixed$timestamp_processing))]

table_ebull$uniqAssessID <- table_draw_fixed$uniqAssessID[match(paste0(table_ebull$username,table_ebull$UniqueID,table_ebull$timestamp_processing), 
                                                                paste0(table_draw_fixed$username,table_draw_fixed$UniqueID,table_draw_fixed$timestamp_processing))]

# ---- Re-compute fluxes ----

load_incubation <- function(auxfile_i, RData_path){
  
  message("Loading data for ",auxfile_i$UniqueID)
  gas <- unique(auxfile_i$gas_analiser)
  
  setwd(RData_path)
  if(gas== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gas
  }
  load(file = paste0(auxfile_i$subsite,"_",gs_suffix,".RData"))
  mydata <- mydata[,c("POSIX.time", 
                      "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                      "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
  
  mydata <- mydata[which(mydata$POSIX.time>=auxfile_i$start.time & 
                           mydata$POSIX.time<=auxfile_i$start.time+auxfile_i$duration),]
  
  if(dim(mydata)[1]>0){
    mydata$Etime <- as.numeric(mydata$POSIX.time) - min(as.numeric(mydata$POSIX.time))
    mydata$UniqueID <- auxfile_i$UniqueID
  } else {
    warning(paste0("For ", auxfile_i$UniqueID, ", no measurements were found!"))
  }
  
  return(mydata)
}




load_incubation_subsite <- function(subsite, auxfile, RData_path){
  
  message("Loading data for ",subsite)
  
  auxfile_i <- auxfile[auxfile$subsite==subsite,]
  
  
  gass <- unique(auxfile_i$gas_analiser)
  
  for(gas in gass){
    
    setwd(RData_path)
    if(gas== "LI-COR"){
      gs_suffix <- "LI-7810"
    } else {
      gs_suffix <- gas
    }
    load(file = paste0(subsite,"_",gs_suffix,".RData"))
    mydata <- mydata[,c("POSIX.time", 
                        "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                        "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
    
    mydata_id <- NULL
    for (id in unique(auxfile_i$UniqueID)){
      auxfile_j <- auxfile_i[which(auxfile_i$UniqueID==id),]
      
      mydata.tmp <- mydata[which(mydata$POSIX.time>=auxfile_j$start.time & 
                               mydata$POSIX.time<=auxfile_j$start.time+auxfile_j$duration),]
      mydata.tmp$Etime <- as.numeric(mydata.tmp$POSIX.time) - min(as.numeric(mydata.tmp$POSIX.time))
      mydata.tmp$UniqueID <- auxfile_j$UniqueID
      
      mydata_id <- rbind(mydata_id, mydata.tmp)
      rm(mydata)
    }
    
  }
  
  return(mydata)
}

# discard all incubations where errors in fieldsheets have been detected

table_draws <- table_draws[!table_draws$UniqueID %in% list_different,]


CO2_flux.auto <- CH4_flux.auto <- CO2_flux.man <- CH4_flux.man <- CH4_diff_flux.man <- df.no_measurements <- NULL
list_ids <- unique(table_draws$UniqueID)
for (k in seq_along(list_ids)){
  i = list_ids[k]
  message(paste0("processing ",i))
  
  table_draw_i <- table_draws[which(table_draws$UniqueID==i),]
  auxfile_i <- auxfile[which(auxfile$UniqueID==i),]
  
  if(dim(auxfile_i)[1]==0){
    message(paste0("Could not find corresponding auxfile for ",i))
    df.no_measurements <- rbind(df.no_measurements,
                                data.frame(UniqueID=i,
                                           message = "no auxfile"))
  } else {
    mydata <- load_incubation(auxfile_i, RData_path)
    
    if(dim(mydata)[1]==0){
      df.no_measurements <- rbind(df.no_measurements,
                                  data.frame(UniqueID=i,
                                             message = "no measurements"))
    } else {
      
      # print(plot_incubations(dataframe = mydata))
      
      # message("... calculate CO2 flux with aquaGHG, no manual selection")
      CO2_flux.auto_i <- automaticflux(dataframe = mydata, myauxfile = auxfile_i, shoulder = 0, gastype = "CO2dry_ppm", 
                                       fluxSeparation = FALSE, displayPlots = FALSE, method = "trust.it.all")
      # message("... calculate CH4 flux with aquaGHG, no manual selection")
      CH4_flux.auto_i <- automaticflux(dataframe = mydata, myauxfile = auxfile_i, shoulder = 0, gastype = "CH4dry_ppb", 
                                       fluxSeparation = T, displayPlots = FALSE, method = "trust.it.all")
      
      CO2_flux.auto <- rbind(CO2_flux.auto, CO2_flux.auto_i)
      CH4_flux.auto <- rbind(CH4_flux.auto, CH4_flux.auto_i)
      
      for(j in table_draw_i$uniqAssessID){
        message(paste0("...processing ",j))
        table_draw_j <- table_draw_i[which(table_draw_i$uniqAssessID==j),]
        
        # change auxfile values based on outputs from expert
        auxfile_j_co2 <- auxfile_j_ch4 <- auxfile_i
        auxfile_j_co2$start.time <- as.POSIXct(table_draw_j$start.time_expert_co2, tz = 'UTC')
        auxfile_j_co2$obs.length <- table_draw_j$duration_expert_co2
        
        auxfile_j_ch4$start.time <- as.POSIXct(table_draw_j$start.time_expert_ch4, tz = 'UTC')
        auxfile_j_ch4$obs.length <- table_draw_j$duration_expert_ch4
        
        if (auxfile_j_co2$obs.length > 30){
          # message("... calculate CO2 flux with aquaGHG, using expert selection")
          CO2_flux.man_j <- automaticflux(dataframe = mydata, myauxfile = auxfile_j_co2, shoulder = 0, gastype = "CO2dry_ppm", 
                                          fluxSeparation = FALSE, displayPlots = FALSE, method = "trust.it.all")
          # message("... calculate CH4 flux with aquaGHG, using expert selection")
          CH4_flux.man_j <- automaticflux(dataframe = mydata, myauxfile = auxfile_j_co2, shoulder = 0, gastype = "CH4dry_ppb", 
                                          fluxSeparation = T, displayPlots = FALSE, method = "trust.it.all")
          
          # assigning draw number to ease further data analysis
          CO2_flux.man_j$draw <- CH4_flux.man_j$draw <- table_draw_j$draw_corrected
          
          # assigning uniqAssessID
          CO2_flux.man_j$uniqAssessID <- CH4_flux.man_j$uniqAssessID <- j
          
          CO2_flux.man <- rbind(CO2_flux.man, CO2_flux.man_j)
          CH4_flux.man <- rbind(CH4_flux.man, CH4_flux.man_j)
        }
        if (auxfile_j_ch4$obs.length > 30){
          # message("... calculate CH4 diffusion flux with aquaGHG, using expert selection")
          CH4_diff_flux.man_j <- automaticflux(dataframe = mydata, myauxfile = auxfile_j_ch4, shoulder = 0, gastype = "CH4dry_ppb", 
                                               fluxSeparation = FALSE, displayPlots = FALSE, method = "trust.it.all")
          
          # assigning draw number to ease further data analysis
          CH4_diff_flux.man_j$draw <- table_draw_j$draw
          CH4_diff_flux.man_j$uniqAssessID <- j
          CH4_diff_flux.man <- rbind(CH4_diff_flux.man, CH4_diff_flux.man_j)
        }
        rm(j)
      }
    }
    rm(auxfile_i, table_draw_i)
    }
}


# ggplot(CH4_diff_flux.man, aes(as.factor(draw), best.flux))+geom_boxplot()+theme_article()



# Save files

setwd(results_path)
save(list = c("CO2_flux.auto",
              "CH4_flux.auto",
              "CO2_flux.man",
              "CH4_flux.man",
              "CH4_diff_flux.man",
              "df.no_measurements",
              "table_draws"), file = "results_ghg_experts_recalculated.Rdata")







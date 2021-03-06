###############################################
# Cette script doit aider à charger les données
# Author : Pham Nguyen Hoang
###############################################


# Charger les packages nécessaire. Il faut que tu bien d'installer toutes les pakages nécessaires.
library(sp)
library(maptools)
library(RColorBrewer)
library(lattice)
library(vars)
library(lubridate)
library(rEDM)
library(dplyr) 
library(tidyr)
library(metap)



# Charger les données climatiques
load("climatic_data_Vietnam_Laos_Thailand.rdata")
load("demo_meteo.RData")

# Charger les données de la dengue
main_data <- read.csv("main_data_for_kmeans.csv", row.names = 1)
# Charger les données géographies du Vietnam
load("VNM_adm2.RData", vn <- new.env())
vn <- vn$gadm
vn  <- thinnedSpatialPoly(vn,tolerance=0.05, minarea=0.001)
provinces_vn <- vn@data$VARNAME_2

# Étape trier les données dengues
dengue_data_vn <- main_data[rownames(main_data) %in% provinces_vn,]
dengue_data_vn <- as.data.frame(t(dengue_data_vn))

# Étape trier les données climatiques
meteo_data_vn <- merge(data_vietnam,stations, by.x = "station", by.y = "station", all.x = TRUE, sort = TRUE)
meteo_data_vn_1 <- meteo_data_vn[which(meteo_data_vn$year >=1994),]
meteo_data_vn_1 <- meteo_data_vn_1[with(meteo_data_vn_1, order(year,month)),]
provinces_meteo <- as.character(unique(meteo_data_vn_1$station))
stations_excl <- c("Bavi","Xuanloc","Laocai","Bacninh","Pharang","TSNhat")
provinces_meteo <- provinces_meteo[!(provinces_meteo %in% stations_excl)]
meteo_data_vn_1$altitude <- NULL


#la fonction pour choisir les données climatique
select_variable <- function(v,s,data=meteo_data_vn_1) 
{
  # select the station and the variable (plus year and month):
  out <- subset(data,station==s,c(v,"year","month"))
  # order chronologically:
  out <- out[with(out,order(year,month)),]
  # return output (i.e. only the variable, without year and month):
  return(out[,1])
}


#Choisir les données climatiques
variables <- c("Ta","Tx","Tm","Rf","rH","Sh","aH","latitude","longitude")
clim_ts <- lapply(provinces_meteo,function(v) as.data.frame(sapply(variables,function(x)select_variable(x,v))))
names(clim_ts) <- provinces_meteo

# Exemple pour choisir les données
# Tu peut regarder la struture des données pour bien comprendre l'utilisation des données.

# Choisir les données de la dengue au Provinces "An Giang"
dengue_data_vn$`An Giang`
#Choisir les données de la dengue au Provinces numéro 1
dengue_data_vn[,1]
# Choisir 7 donné climatiques au Province Bac Can
clim_ts$Baccan
# Choisir la facteur climatiques temperature moyenne au Province Bac Can
clim_ts$Baccan$Ta

# Bon Courage

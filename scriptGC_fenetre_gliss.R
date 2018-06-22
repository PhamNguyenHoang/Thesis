################################################################################ 
## Cette Script appliqué la méthode Causalité de Granger par fenêtre glissant ##
##                        Author : Pham Nguyen Hoang                          ##
################################################################################

library(sp)
library(biwavelet)
library(TSclust)
library(TSdist)
library(fpc)
library(rgeos)
library(vars)
library(maptools)
library(RColorBrewer)
library(lattice)
library(vars)
library(lubridate)

load("climatic_data_Vietnam_Laos_Thailand.rdata")
load("demo_meteo.RData")
main_data <- read.csv("main_data_for_kmeans.csv", row.names = 1)
load("VNM_adm2.RData", vn <- new.env())
vn <- vn$gadm
vn  <- thinnedSpatialPoly(vn,tolerance=0.05, minarea=0.001)
## Dengue data
provinces_vn <- vn@data$VARNAME_2

dengue_data_vn <- main_data[rownames(main_data) %in% provinces_vn,]

dengue_data_vn <- as.data.frame(t(dengue_data_vn))
meteo_data_vn <- merge(data_vietnam,stations, by.x = "station", by.y = "station", all.x = TRUE, sort = TRUE)
meteo_data_vn_1 <- meteo_data_vn[which(meteo_data_vn$year >=1994),]

meteo_data_vn_1 <- meteo_data_vn_1[with(meteo_data_vn_1, order(year,month)),]

provinces_meteo <- as.character(unique(meteo_data_vn_1$station))
stations_excl <- c("Bavi","Xuanloc","Laocai","Bacninh","Pharang","TSNhat")
provinces_meteo <- provinces_meteo[!(provinces_meteo %in% stations_excl)]

meteo_data_vn_1$altitude <- NULL

select_variable <- function(v,s,data=meteo_data_vn_1) 
{
  # select the station and the variable (plus year and month):
  out <- subset(data,station==s,c(v,"year","month"))
  # order chronologically:
  out <- out[with(out,order(year,month)),]
  # return output (i.e. only the variable, without year and month):
  return(out[,1])
}

variables <- c("Ta","Tx","Tm","Rf","rH","Sh","aH","latitude","longitude")
clim_ts <- lapply(provinces_meteo,function(v) as.data.frame(sapply(variables,function(x)select_variable(x,v))))
names(clim_ts) <- provinces_meteo

cal_dis_eucl <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))



sig_level = 0.05
tmp = 1:64
excl = c(26,27, 28, 34, 35)
tmp = tmp[!(tmp %in% excl)]

# paramètre indique la longueur de la fenetre
fenetre_long <- 12 * 9 # unité = mois, modifie nombre d'année

# paramètre indique le décalage de la fenetre
fenetre_lag <- 3  # unité = mois

num_rep <- seq(1,144 - fenetre_long + fenetre_lag, by = fenetre_lag)

#matrix stocke les resultats
matrix_result1 <- NULL
matrix_result2 <- NULL
matrix_result3 <- NULL
matrix_result4 <- NULL
matrix_result5 <- NULL
matrix_result6 <- NULL
matrix_result7 <- NULL

matrix_result_1_4 <- NULL
matrix_result_1_5 <- NULL
matrix_result_1_6 <- NULL
matrix_result_1_7 <- NULL
matrix_result_4_5 <- NULL
matrix_result_4_6 <- NULL
matrix_result_4_7 <- NULL
matrix_result_5_6 <- NULL
matrix_result_5_7 <- NULL
matrix_result_6_7 <- NULL

for (t in num_rep)
{
  causality_result_1 <- NULL
  causality_result_2 <- NULL
  causality_result_3 <- NULL
  causality_result_4 <- NULL
  causality_result_5 <- NULL
  causality_result_6 <- NULL
  causality_result_7 <- NULL
  
  causality_result_1_4 <- NULL
  causality_result_1_5 <- NULL
  causality_result_1_6 <- NULL
  causality_result_1_7 <- NULL
  causality_result_4_5 <- NULL
  causality_result_4_6 <- NULL
  causality_result_4_7 <- NULL
  causality_result_5_6 <- NULL
  causality_result_5_7 <- NULL
  causality_result_6_7 <- NULL
  
  for(k in tmp)
  {
    num = k # number of province in the list
    vari = c(1,2,3,4,5,6,7) 
    prov_c <- vn@data$VARNAME_2[as.numeric(num)]
    vari_c <- variables[as.numeric(vari)]
    
    #take a coordinate of provinces 
    x <- mean(vn@polygons[[as.numeric(num)]]@Polygons[[1]]@coords[,1]) 
    y <- mean(vn@polygons[[as.numeric(num)]]@Polygons[[1]]@coords[,2])
    prov_coor <- cbind(x,y)
    prov_coor <- as.data.frame(prov_coor)
    colnames(prov_coor) <- c("longitude","latitude")
    
    
    d <- dengue_data_vn[,which( colnames(dengue_data_vn) == prov_c )]
    
    dis_eu <- NULL
    for ( j in 1:63)
    {
      x1 = cbind(mean(clim_ts[[j]]$longitude),mean(clim_ts[[j]]$latitude))
      
      a <- cal_dis_eucl(prov_coor,x1)
      dis_eu <- rbind(dis_eu,c(a, attributes(clim_ts[j])))
    }
    
    dis_eu <- as.data.frame(dis_eu)
    dis_eu$V1 <- as.numeric(dis_eu$V1)
    climat_result <- dis_eu[ which(dis_eu[,1] == min(dis_eu[,1]) ) , ]
    
    climat_data_1 <- select_variable(vari_c[1],climat_result$names)
    climat_data_2 <- select_variable(vari_c[2],climat_result$names)
    climat_data_3 <- select_variable(vari_c[3],climat_result$names)
    climat_data_4 <- select_variable(vari_c[4],climat_result$names)
    climat_data_5 <- select_variable(vari_c[5],climat_result$names)
    climat_data_6 <- select_variable(vari_c[6],climat_result$names)
    climat_data_7 <- select_variable(vari_c[7],climat_result$names)
    climat_1 <- climat_data_1[49:192]
    climat_2 <- climat_data_2[49:192]
    climat_3 <- climat_data_3[49:192]
    climat_4 <- climat_data_4[49:192]
    climat_5 <- climat_data_5[49:192]
    climat_6 <- climat_data_6[49:192]
    climat_7 <- climat_data_7[49:192]
    #Test 1
    
    da <- d[t : (t+fenetre_long - 1)]
    climat_1a <- climat_1[t : (t+fenetre_long - 1)]
    climat_2a <- climat_2[t : (t+fenetre_long - 1)]
    climat_3a <- climat_3[t : (t+fenetre_long - 1)]
    climat_4a <- climat_4[t : (t+fenetre_long - 1)]
    climat_5a <- climat_5[t : (t+fenetre_long - 1)]
    climat_6a <- climat_6[t : (t+fenetre_long - 1)]
    climat_7a <- climat_7[t : (t+fenetre_long - 1)]
    
    test_case_1a = data.frame(cbind(da,climat_1a))
    test_case_2a = data.frame(cbind(da,climat_2a))
    test_case_3a = data.frame(cbind(da,climat_3a))
    test_case_4a = data.frame(cbind(da,climat_4a))
    test_case_5a = data.frame(cbind(da,climat_5a))
    test_case_6a = data.frame(cbind(da,climat_6a))
    test_case_7a = data.frame(cbind(da,climat_7a))
    
    test_case_1_4a = data.frame(cbind(da,climat_1a,climat_4a))
    test_case_1_5a = data.frame(cbind(da,climat_1a,climat_5a))
    test_case_1_6a = data.frame(cbind(da,climat_1a,climat_6a))
    test_case_1_7a = data.frame(cbind(da,climat_1a,climat_7a))
    test_case_4_5a = data.frame(cbind(da,climat_4a,climat_5a))
    test_case_4_6a = data.frame(cbind(da,climat_4a,climat_6a))
    test_case_4_7a = data.frame(cbind(da,climat_4a,climat_7a))
    test_case_5_6a = data.frame(cbind(da,climat_5a,climat_6a))
    test_case_5_7a = data.frame(cbind(da,climat_5a,climat_7a))
    test_case_6_7a = data.frame(cbind(da,climat_6a,climat_7a))
    
    ############################################################################
    Var_dengue_1a <- VAR(test_case_1a, p = 1, type = "const")
    Var_dengue_2a <- VAR(test_case_2a, p = 1, type = "const")
    Var_dengue_3a <- VAR(test_case_3a, p = 1, type = "const")
    Var_dengue_4a <- VAR(test_case_4a, p = 1, type = "const")
    Var_dengue_5a <- VAR(test_case_5a, p = 1, type = "const")
    Var_dengue_6a <- VAR(test_case_6a, p = 1, type = "const")
    Var_dengue_7a <- VAR(test_case_7a, p = 1, type = "const")
    
    Var_dengue_1_4a <- VAR(test_case_1_4a, p = 1, type = "const")
    Var_dengue_1_5a <- VAR(test_case_1_5a, p = 1, type = "const")
    Var_dengue_1_6a <- VAR(test_case_1_6a, p = 1, type = "const")
    Var_dengue_1_7a <- VAR(test_case_1_7a, p = 1, type = "const")
    Var_dengue_4_5a <- VAR(test_case_4_5a, p = 1, type = "const")    
    Var_dengue_4_6a <- VAR(test_case_4_6a, p = 1, type = "const")
    Var_dengue_4_7a <- VAR(test_case_4_7a, p = 1, type = "const")
    Var_dengue_5_6a <- VAR(test_case_5_6a, p = 1, type = "const")
    Var_dengue_5_7a <- VAR(test_case_5_7a, p = 1, type = "const")
    Var_dengue_6_7a <- VAR(test_case_6_7a, p = 1, type = "const")
    
    ############################################################################
    cause_1a <- causality(Var_dengue_1a, cause = c("climat_1a"))
    cause_2a <- causality(Var_dengue_2a, cause = c("climat_2a"))
    cause_3a <- causality(Var_dengue_3a, cause = c("climat_3a"))
    cause_4a <- causality(Var_dengue_4a, cause = c("climat_4a"))
    cause_5a <- causality(Var_dengue_5a, cause = c("climat_5a"))
    cause_6a <- causality(Var_dengue_6a, cause = c("climat_6a"))
    cause_7a <- causality(Var_dengue_7a, cause = c("climat_7a"))
    
    cause_1_4a <- causality(Var_dengue_1_4a, cause = c("climat_1a","climat_4a"))
    cause_1_5a <- causality(Var_dengue_1_5a, cause = c("climat_1a","climat_5a"))
    cause_1_6a <- causality(Var_dengue_1_6a, cause = c("climat_1a","climat_6a"))
    cause_1_7a <- causality(Var_dengue_1_7a, cause = c("climat_1a","climat_7a"))
    cause_4_5a <- causality(Var_dengue_4_5a, cause = c("climat_4a","climat_5a"))
    cause_4_6a <- causality(Var_dengue_4_6a, cause = c("climat_4a","climat_6a"))
    cause_4_7a <- causality(Var_dengue_4_7a, cause = c("climat_4a","climat_7a"))
    cause_5_6a <- causality(Var_dengue_5_6a, cause = c("climat_5a","climat_6a"))
    cause_5_7a <- causality(Var_dengue_5_7a, cause = c("climat_5a","climat_7a"))
    cause_6_7a <- causality(Var_dengue_6_7a, cause = c("climat_6a","climat_7a"))
    
    reject_1a <- NULL
    reject_2a <- NULL
    reject_3a <- NULL
    reject_4a <- NULL
    reject_5a <- NULL
    reject_6a <- NULL
    reject_7a <- NULL
    
    reject_1_4a <- NULL
    reject_1_5a <- NULL
    reject_1_6a <- NULL
    reject_1_7a <- NULL
    reject_4_5a <- NULL
    reject_4_6a <- NULL
    reject_4_7a <- NULL
    reject_5_6a <- NULL
    reject_5_7a <- NULL
    reject_6_7a <- NULL
    
    ifelse(cause_1a$Granger$p.value <= sig_level,reject_1a <- 1 , reject_1a <- 0)
    ifelse(cause_2a$Granger$p.value <= sig_level,reject_2a <- 1 , reject_2a <- 0)
    ifelse(cause_3a$Granger$p.value <= sig_level,reject_3a <- 1 , reject_3a <- 0)
    ifelse(cause_4a$Granger$p.value <= sig_level,reject_4a <- 1 , reject_4a <- 0)
    ifelse(cause_5a$Granger$p.value <= sig_level,reject_5a <- 1 , reject_5a <- 0)
    ifelse(cause_6a$Granger$p.value <= sig_level,reject_6a <- 1 , reject_6a <- 0)
    ifelse(cause_7a$Granger$p.value <= sig_level,reject_7a <- 1 , reject_7a <- 0)
    
    ifelse(cause_1_4a$Granger$p.value <= sig_level,reject_1_4a <- 1 , reject_1_4a <- 0)
    ifelse(cause_1_5a$Granger$p.value <= sig_level,reject_1_5a <- 1 , reject_1_5a <- 0)
    ifelse(cause_1_6a$Granger$p.value <= sig_level,reject_1_6a <- 1 , reject_1_6a <- 0)
    ifelse(cause_1_7a$Granger$p.value <= sig_level,reject_1_7a <- 1 , reject_1_7a <- 0)
    ifelse(cause_4_5a$Granger$p.value <= sig_level,reject_4_5a <- 1 , reject_4_5a <- 0)
    ifelse(cause_4_6a$Granger$p.value <= sig_level,reject_4_6a <- 1 , reject_4_6a <- 0)
    ifelse(cause_4_7a$Granger$p.value <= sig_level,reject_4_7a <- 1 , reject_4_7a <- 0)
    ifelse(cause_5_6a$Granger$p.value <= sig_level,reject_5_6a <- 1 , reject_5_6a <- 0)
    ifelse(cause_5_7a$Granger$p.value <= sig_level,reject_5_7a <- 1 , reject_5_7a <- 0)
    ifelse(cause_6_7a$Granger$p.value <= sig_level,reject_6_7a <- 1 , reject_6_7a <- 0)
    
    
    #######################################################################  
    causality_result_1 <- rbind(causality_result_1, c(prov_c, cause_1a$Granger$statistic, cause_1a$Granger$p.value, reject_1a ))
    causality_result_2 <- rbind(causality_result_2, c(prov_c, cause_2a$Granger$statistic, cause_2a$Granger$p.value, reject_2a ))
    causality_result_3 <- rbind(causality_result_3, c(prov_c, cause_3a$Granger$statistic, cause_3a$Granger$p.value, reject_3a ))
    causality_result_4 <- rbind(causality_result_4, c(prov_c, cause_4a$Granger$statistic, cause_4a$Granger$p.value, reject_4a ))
    causality_result_5 <- rbind(causality_result_5, c(prov_c, cause_5a$Granger$statistic, cause_5a$Granger$p.value, reject_5a ))
    causality_result_6 <- rbind(causality_result_6, c(prov_c, cause_6a$Granger$statistic, cause_6a$Granger$p.value, reject_6a ))
    causality_result_7 <- rbind(causality_result_7, c(prov_c, cause_7a$Granger$statistic, cause_7a$Granger$p.value, reject_7a ))
    
    causality_result_1_4 <- rbind(causality_result_1_4, c(prov_c, cause_1_4a$Granger$statistic, cause_1_4a$Granger$p.value, reject_1_4a ))
    causality_result_1_5 <- rbind(causality_result_1_5, c(prov_c, cause_1_5a$Granger$statistic, cause_1_5a$Granger$p.value, reject_1_5a ))
    causality_result_1_6 <- rbind(causality_result_1_6, c(prov_c, cause_1_6a$Granger$statistic, cause_1_6a$Granger$p.value, reject_1_6a ))
    causality_result_1_7 <- rbind(causality_result_1_7, c(prov_c, cause_1_7a$Granger$statistic, cause_1_7a$Granger$p.value, reject_1_7a ))
    causality_result_4_5 <- rbind(causality_result_4_5, c(prov_c, cause_4_5a$Granger$statistic, cause_4_5a$Granger$p.value, reject_4_5a ))
    causality_result_4_6 <- rbind(causality_result_4_6, c(prov_c, cause_4_6a$Granger$statistic, cause_4_6a$Granger$p.value, reject_4_6a ))
    causality_result_4_7 <- rbind(causality_result_4_7, c(prov_c, cause_4_7a$Granger$statistic, cause_4_7a$Granger$p.value, reject_4_7a ))
    causality_result_5_6 <- rbind(causality_result_5_6, c(prov_c, cause_5_6a$Granger$statistic, cause_5_6a$Granger$p.value, reject_5_6a ))
    causality_result_5_7 <- rbind(causality_result_5_7, c(prov_c, cause_5_7a$Granger$statistic, cause_5_7a$Granger$p.value, reject_5_7a ))
    causality_result_6_7 <- rbind(causality_result_6_7, c(prov_c, cause_6_7a$Granger$statistic, cause_6_7a$Granger$p.value, reject_6_7a ))
    
    
  }
  
  causality_result_1 <- as.data.frame(causality_result_1)
  causality_result_2 <- as.data.frame(causality_result_2)
  causality_result_3 <- as.data.frame(causality_result_3)
  causality_result_4 <- as.data.frame(causality_result_4)
  causality_result_5 <- as.data.frame(causality_result_5)
  causality_result_6 <- as.data.frame(causality_result_6)
  causality_result_7 <- as.data.frame(causality_result_7)
  
  causality_result_1_4 <- as.data.frame(causality_result_1_4)
  causality_result_1_5 <- as.data.frame(causality_result_1_5)
  causality_result_1_6 <- as.data.frame(causality_result_1_6)
  causality_result_1_7 <- as.data.frame(causality_result_1_7)
  causality_result_4_5 <- as.data.frame(causality_result_4_5)
  causality_result_4_6 <- as.data.frame(causality_result_4_6)
  causality_result_4_7 <- as.data.frame(causality_result_4_7)
  causality_result_5_6 <- as.data.frame(causality_result_5_6)
  causality_result_5_7 <- as.data.frame(causality_result_5_7)
  causality_result_6_7 <- as.data.frame(causality_result_6_7)
  
  causality_result_1$`F-Test` <- NULL
  causality_result_1$V3 <- NULL
  causality_result_2$`F-Test` <- NULL
  causality_result_2$V3 <- NULL
  causality_result_3$`F-Test` <- NULL
  causality_result_3$V3 <- NULL
  causality_result_4$`F-Test` <- NULL
  causality_result_4$V3 <- NULL
  causality_result_5$`F-Test` <- NULL
  causality_result_5$V3 <- NULL
  causality_result_6$`F-Test` <- NULL
  causality_result_6$V3 <- NULL
  causality_result_7$`F-Test` <- NULL
  causality_result_7$V3 <- NULL
  
  causality_result_1_4$`F-Test` <- NULL
  causality_result_1_4$V3 <- NULL
  causality_result_1_5$`F-Test` <- NULL
  causality_result_1_5$V3 <- NULL
  causality_result_1_6$`F-Test` <- NULL
  causality_result_1_6$V3 <- NULL
  causality_result_1_7$`F-Test` <- NULL
  causality_result_1_7$V3 <- NULL
  causality_result_4_5$`F-Test` <- NULL
  causality_result_4_5$V3 <- NULL
  causality_result_4_6$`F-Test` <- NULL
  causality_result_4_6$V3 <- NULL
  causality_result_4_7$`F-Test` <- NULL
  causality_result_4_7$V3 <- NULL
  causality_result_5_6$`F-Test` <- NULL
  causality_result_5_6$V3 <- NULL
  causality_result_5_7$`F-Test` <- NULL
  causality_result_5_7$V3 <- NULL
  causality_result_6_7$`F-Test` <- NULL
  causality_result_6_7$V3 <- NULL
  
  
  colum1 <- NULL
  colum2 <- NULL
  colum3 <- NULL
  colum4 <- NULL
  colum5 <- NULL
  colum6 <- NULL
  colum7 <- NULL
  
  colum_1_4 <- NULL
  colum_1_5 <- NULL
  colum_1_6 <- NULL
  colum_1_7 <- NULL
  colum_4_5 <- NULL
  colum_4_6 <- NULL
  colum_4_7 <- NULL
  colum_5_6 <- NULL
  colum_5_7 <- NULL
  colum_6_7 <- NULL
  
  colum1 <- cbind(colum1,as.numeric(as.character(causality_result_1$V4)))
  colum2 <- cbind(colum2,as.numeric(as.character(causality_result_2$V4)))
  colum3 <- cbind(colum3,as.numeric(as.character(causality_result_3$V4)))
  colum4 <- cbind(colum4,as.numeric(as.character(causality_result_4$V4)))
  colum5 <- cbind(colum5,as.numeric(as.character(causality_result_5$V4)))
  colum6 <- cbind(colum6,as.numeric(as.character(causality_result_6$V4)))
  colum7 <- cbind(colum7,as.numeric(as.character(causality_result_7$V4)))
  
  colum_1_4 <- cbind(colum_1_4,as.numeric(as.character(causality_result_1_4$V4)))
  colum_1_5 <- cbind(colum_1_5,as.numeric(as.character(causality_result_1_5$V4)))
  colum_1_6 <- cbind(colum_1_6,as.numeric(as.character(causality_result_1_6$V4)))
  colum_1_7 <- cbind(colum_1_7,as.numeric(as.character(causality_result_1_7$V4)))
  colum_4_5 <- cbind(colum_4_5,as.numeric(as.character(causality_result_4_5$V4)))
  colum_4_6 <- cbind(colum_4_6,as.numeric(as.character(causality_result_4_6$V4)))
  colum_4_7 <- cbind(colum_4_7,as.numeric(as.character(causality_result_4_7$V4)))
  colum_5_6 <- cbind(colum_5_6,as.numeric(as.character(causality_result_5_6$V4)))
  colum_5_7 <- cbind(colum_5_7,as.numeric(as.character(causality_result_5_7$V4)))
  colum_6_7 <- cbind(colum_6_7,as.numeric(as.character(causality_result_6_7$V4)))
  
#   colum <- rbind(colum,sum(as.numeric(as.character(causality_result_1$V4)),na.rm = TRUE))
#   colum <- rbind(colum,sum(as.numeric(as.character(causality_result_2$V4)),na.rm = TRUE))
#   colum <- rbind(colum,sum(as.numeric(as.character(causality_result_3$V4)),na.rm = TRUE))
#   colum <- rbind(colum,sum(as.numeric(as.character(causality_result_4$V4)),na.rm = TRUE))
#   colum <- rbind(colum,sum(as.numeric(as.character(causality_result_5$V4)),na.rm = TRUE))
#   colum <- rbind(colum,sum(as.numeric(as.character(causality_result_6$V4)),na.rm = TRUE))
#   colum <- rbind(colum,sum(as.numeric(as.character(causality_result_7$V4)),na.rm = TRUE))
  
  matrix_result1 <- cbind(matrix_result1,colum1)
  matrix_result2 <- cbind(matrix_result2,colum2)
  matrix_result3 <- cbind(matrix_result3,colum3)
  matrix_result4 <- cbind(matrix_result4,colum4)
  matrix_result5 <- cbind(matrix_result5,colum5)
  matrix_result6 <- cbind(matrix_result6,colum6)
  matrix_result7 <- cbind(matrix_result7,colum7)
  
  matrix_result_1_4 <- cbind(matrix_result_1_4,colum_1_4)
  matrix_result_1_5 <- cbind(matrix_result_1_5,colum_1_5)
  matrix_result_1_6 <- cbind(matrix_result_1_6,colum_1_6)
  matrix_result_1_7 <- cbind(matrix_result_1_7,colum_1_7)
  matrix_result_4_5 <- cbind(matrix_result_4_5,colum_4_5)
  matrix_result_4_6 <- cbind(matrix_result_4_6,colum_4_6)
  matrix_result_4_7 <- cbind(matrix_result_4_7,colum_4_7)
  matrix_result_5_6 <- cbind(matrix_result_5_6,colum_5_6)
  matrix_result_5_7 <- cbind(matrix_result_5_7,colum_5_7)
  matrix_result_6_7 <- cbind(matrix_result_6_7,colum_6_7)
  
}



# rownames(matrix_result1) <- causality_result_1$V1
# rownames(matrix_result2) <- causality_result_1$V1
# rownames(matrix_result3) <- causality_result_1$V1
# rownames(matrix_result4) <- causality_result_1$V1
# rownames(matrix_result5) <- causality_result_1$V1
# rownames(matrix_result6) <- causality_result_1$V1
# rownames(matrix_result7) <- causality_result_1$V1
# 
# 
# rownames(matrix_result_1_4) <- causality_result_1$V1
# rownames(matrix_result_1_5) <- causality_result_1$V1
# rownames(matrix_result_1_6) <- causality_result_1$V1
# rownames(matrix_result_1_7) <- causality_result_1$V1
# rownames(matrix_result_4_5) <- causality_result_1$V1
# rownames(matrix_result_4_6) <- causality_result_1$V1
# rownames(matrix_result_4_7) <- causality_result_1$V1
# rownames(matrix_result_5_6) <- causality_result_1$V1
# rownames(matrix_result_5_7) <- causality_result_1$V1
# rownames(matrix_result_6_7) <- causality_result_1$V1

sum1 <- colSums(x = matrix_result1)
sum2 <- colSums(x = matrix_result2)
sum3 <- colSums(x = matrix_result3)
sum4 <- colSums(x = matrix_result4)
sum5 <- colSums(x = matrix_result5)
sum6 <- colSums(x = matrix_result6)
sum7 <- colSums(x = matrix_result7)

sum1_4 <- colSums(x = matrix_result_1_4)
sum1_5 <- colSums(x = matrix_result_1_5)
sum1_6 <- colSums(x = matrix_result_1_6)
sum1_7 <- colSums(x = matrix_result_1_7)
sum4_5 <- colSums(x = matrix_result_4_5)
sum4_6 <- colSums(x = matrix_result_4_6)
sum4_7 <- colSums(x = matrix_result_4_7)
sum5_6 <- colSums(x = matrix_result_5_6)
sum5_7 <- colSums(x = matrix_result_5_7)
sum6_7 <- colSums(x = matrix_result_6_7)


result_univar <- NULL
result_multivar <- NULL

result_univar <- rbind(result_univar, sum1)
result_univar <- rbind(result_univar, sum2)
result_univar <- rbind(result_univar, sum3)
result_univar <- rbind(result_univar, sum4)
result_univar <- rbind(result_univar, sum5)
result_univar <- rbind(result_univar, sum6)
result_univar <- rbind(result_univar, sum7)

result_multivar <- rbind(result_multivar, sum1_4)
result_multivar <- rbind(result_multivar, sum1_5)
result_multivar <- rbind(result_multivar, sum1_6)
result_multivar <- rbind(result_multivar, sum1_7)
result_multivar <- rbind(result_multivar, sum4_5)
result_multivar <- rbind(result_multivar, sum4_6)
result_multivar <- rbind(result_multivar, sum4_7)
result_multivar <- rbind(result_multivar, sum5_6)
result_multivar <- rbind(result_multivar, sum5_7)
result_multivar <- rbind(result_multivar, sum6_7)



date_begin_seq1 <- as.Date("1998/01/01")
date_begin_seq2 <- date_begin_seq1
month(date_begin_seq2) <- month(date_begin_seq2) + fenetre_long - 1

time_seq1 <- format(seq (date_begin_seq1, by = paste(fenetre_lag, "month", sep = " "), length.out = ( (144 - fenetre_long)/fenetre_lag) + 1) , format = "%m/%Y")
time_seq2 <- format(seq (date_begin_seq2, by = paste(fenetre_lag, "month", sep = " "), length.out = ( (144 - fenetre_long)/fenetre_lag) + 1) , format = "%m/%Y")
time_seq <- NULL

temp <- seq(1,length(num_rep), by = 1 )

for ( l in  temp)
{
  time_seq[l] <- paste(time_seq1[l], time_seq2[l], sep = " - ")
}


# matplot(t(result_univar),type = "l",
#         main = "Frequences of Provinces has significative between DHF and Absolute Humidity in Vietnam classify by latitude",
#         xlab = "", ylab = "", xaxt = "n" , pch = 15, cex = 1.0, yaxt = "n", col = "red")
# axis(1, at = temp, label = FALSE , cex.axis = 0.6)
# axis(2, at = seq(0,max(result_univar), by = 1) , cex.axis = 0.5, las = 2)
# # abline(h =seq(1,nrow(result_univar), by = 1), col="darkgray", lty="dotted")
# text(x=temp, y=par()$usr[3]-1,
#      labels=time_seq, srt=60, adj=1, xpd=TRUE, cex = 0.6)

colnames(result_univar) <- time_seq
rownames(result_univar) <- c("ta","tx","tm","rf","rh","sh","ah")

colnames(result_multivar) <- time_seq
rownames(result_multivar) <- c("ta + rf","ta + rh","ta + sh","ta + ah","rf + rh","rf + sh","rf + ah", "rh + sh", "rh + ah", "sh + ah")

write.csv2(result_univar,"Result_univar_LenWin9.csv")
write.csv2(result_multivar,"Result_multivar_LenWin9.csv")

#matplot(result_univar)






# 
# ###################################
# # đoạn code hiển thị đồ thị tầng suất bệnh phân bố từng thời điểm ở từng tỉnh thành sắp xếp theo vĩ độ
# rs <- NULL
# ##classify by frequence
# ##rs <- rowSums(matrix_result7)
# 
# ##classify by lattitude
# 
# for (i in tmp)
# {
#   rs <- rbind(rs,mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]))
# }
# 
# 
# matrix_result_temp <- matrix_result7
# 
# matrix_result_temp <- cbind(matrix_result_temp,rs)
# matrix_result_temp <- matrix_result_temp[order(rs),]
# #matrix_result_temp <- subset(matrix_result_temp, select =  - rs)
# 
# 
# 
# for (i in seq(1,nrow(matrix_result_temp), by = 1) )
# {
#   for (j in seq(1, ncol(matrix_result_temp)  , by = 1)     )
#   {
#     ifelse(matrix_result_temp[i,j] == 1, matrix_result_temp[i,j] <-  i,matrix_result_temp[i,j] <- NA  )
#   }
#   
# }
# 
# 
# 
# 
# date_begin_seq1 <- as.Date("1998/01/01")
# date_begin_seq2 <- date_begin_seq1
# month(date_begin_seq2) <- month(date_begin_seq2) + fenetre_long - 1
# 
# time_seq1 <- format(seq (date_begin_seq1, by = paste(fenetre_lag, "month", sep = " "), length.out = ( (144 - fenetre_long)/fenetre_lag) + 1) , format = "%m/%Y")
# time_seq2 <- format(seq (date_begin_seq2, by = paste(fenetre_lag, "month", sep = " "), length.out = ( (144 - fenetre_long)/fenetre_lag) + 1) , format = "%m/%Y")
# time_seq <- NULL
# temp <- seq(1,length(num_rep), by = 1 )
# for ( l in  temp)
# {
#   time_seq[l] <- paste(time_seq1[l], time_seq2[l], sep = " - ")
# }
# 
# 
# 
# 
# 
# 
# 
# matplot(t(matrix_result_temp),type = "p",
#         main = "Frequences of Provinces has significative between DHF and Absolute Humidity in Vietnam classify by latitude",
#         xlab = "", ylab = "", xaxt = "n" , pch = 15, cex = 1.0, yaxt = "n", col = "red")
# axis(1, at = temp, label = FALSE , cex.axis = 0.6)
# axis(2, at = seq(1,nrow(matrix_result_temp), by = 1), label = rownames(matrix_result_temp) , cex.axis = 0.5, las = 2)
# abline(h =seq(1,nrow(matrix_result_temp), by = 1), col="darkgray", lty="dotted")
# text(x=temp, y=par()$usr[3]-1,
#      labels=time_seq, srt=60, adj=1, xpd=TRUE, cex = 0.6)




#######################################################################################
## Cette Script teste les codes de Sugihara sur la méthode Convergence Cross Mapping ##
##                        Author : Pham Nguyen Hoang                                 ##
#######################################################################################


library(sp)
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
library(rEDM)
library(dplyr) 
library(tidyr)

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

mean_anomal <- function(x)
{ 
  # x: time-series values to compute seasonal mean and anomaly
  doy <- rep(1:12,len = 153)
  I_use <- which(!is.na(x))
  doy_sm <- rep(doy[I_use],3) + rep(c(-12,0,12),each=length(I_use)) 
  x_sm <- rep(x[I_use],3)
  xsp <- smooth.spline(doy_sm, y = x_sm, w = NULL, spar = 0.8, cv = NA, all.knots = TRUE,keep.data = TRUE, df.offset = 0)
  
  xbar <- data.frame(doy=doy) %>% left_join(data.frame(doy=xsp$x,xbar=xsp$y),by='doy') 

  out = data.frame(mean=xbar$xbar,anomaly=(x - xbar$xbar)) 
  names(out) <- c('mean','anomaly') 
  return(out)

}



# packageVersion("rEDM")
# packageVersion("tidyr")
# 
# maindata <- read.csv("DataCCM.txt",header = TRUE)
# data.vn <- maindata %>% filter(country == 'Japan') %>% 
#   filter(year >= 1996)
# 
# 
# data.vn.trans <- data.vn %>% select(-country) %>% spread(variable,value) %>% mutate(date = ISOdate(year,month,day)) %>% select(-year,-month,-day) %>% select(date,flu,everything())
# 
# set.seed(599213)
# 
# 
# make_block <- function(data,cols,delays,lib=c(1,NROW(data)))
#   {
#   lib <- matrix(lib,ncol = 2) 
#   data <- as.matrix(data)
#   ncol <- length(cols) 
#   nrow <- dim(data)[1] 
#   
#   block <- array(NA,dim = c(nrow,ncol)) 
#   colnames(block) <- 1:ncol
#   for (i in 1:ncol)
#     { 
#     I <- 1:nrow 
#     I_delay <- intersect(I,I+delays[i]) 
#     block[I_delay-delays[i],i] <- data[I_delay,cols[i]] 
#     if (delays[i] < 0)
#       { 
#       # remove data points that fall at start of lib segments 
#       block[lib[,1] - (0:(delays[i]+1)),i] <- NA 
#       colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t-',abs(delays[i]),sep="")
#       }
#     else if (delays[i] > 0) 
#       { # remove data points that fall at end of lib segments 
#       block[lib[,2] - (0:(delays[i]+1)),i] <- NA
#       colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t+',abs(delays[i]),sep="")
#       }
#     else 
#       { 
#         colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t',sep="")
#       } 
#     }
#   return(block)
# }
# 
# 
# yearday_anom <- function(t,x)
# { 
#   # t: date formatted with POSIXt 
#   # x: time-series values to compute seasonal mean and anomaly
#   doy <- as.numeric(strftime(t, format = "%j")) 
#   I_use <- which(!is.na(x))
#   doy_sm <- rep(doy[I_use],3) + rep(c(-366,0,366),each=length(I_use)) 
#   x_sm <- rep(x[I_use],3)
#   xsp <- smooth.spline(doy_sm, y = x_sm, w = NULL, spar = 0.8, cv = NA, all.knots = TRUE,keep.data = TRUE, df.offset = 0)
#   xbar <- data.frame(t=t,doy=doy) %>% left_join(data.frame(doy=xsp$x,xbar=xsp$y),by='doy') %>% select(xbar)
# 
#   out = data.frame(t=t,mean=xbar,anomaly=(x - xbar)) 
#   names(out) <- c('date','mean','anomaly') 
#   return(out)
# 
# }
# 

# AH.surrs <- do.call(cbind, lapply(1:500, function(i) 
#   {
#   I_na <- is.na(AH.tilde)
#   out <- AH.bar
#   out[I_na] <- NA
#   out[!I_na] <- out[!I_na] + sample(AH.tilde[!I_na], sum(!I_na), replace = FALSE)
#   return(out)
#   } ))
# 
# plot(data.vn.trans$date,data.vn.trans$AH, type = "l", col = "blue",xlab = "date", ylab= "Absolute Humidity")
# lines(data.vn.trans$date,AH.surrs[,2],col = "red")
# 
# legend(x = "bottomleft", 
#        legend = c("Observed", "Seasonal Surrogate"), 
#        col = c("blue", "red"), lwd = 1, lty = c(1,1), bg = "white")  


###Cette fonction sera éliminer les valeurs 0 dans les séries temporelles
make_pred_nozero <- function(time_series,E)
  {
  I_zero_strings <- which(time_series==0)
  I_zero_strings <- Reduce(intersect, lapply((0:E),function(offset) I_zero_strings-offset))
  I_zero_strings <- c(0,I_zero_strings,153)
  N_zero_strings <- length(I_zero_strings)
  lib_nozeroes <- cbind(I_zero_strings[1:(N_zero_strings-1)]+1,I_zero_strings[2:(N_zero_strings)]-1)
  
  lib_out <- lib_nozeroes[which(lib_nozeroes[,2] > lib_nozeroes[,1]),]
  return(lib_out)
  }

tmp = 1:64
excl = c(26,27, 28, 34, 35)
tmp = tmp[!(tmp %in% excl)]
# Determiner la valeur optional de E - Embedded Dimension 

ccm_apply <- function (num_prov, num_clim)
{
   # number of province in the list
  vari = c(1,2,3,4,5,6,7) 
  prov_c <- as.character(vn@data$VARNAME_2[as.numeric(num_prov)])
  vari_c <- variables[as.numeric(vari)]
  
  #take a coordinate of provinces 
  x <- mean(vn@polygons[[as.numeric(num_prov)]]@Polygons[[1]]@coords[,1]) 
  y <- mean(vn@polygons[[as.numeric(num_prov)]]@Polygons[[1]]@coords[,2])
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
  
  climat_data_1 <- select_variable(vari_c[num_clim],climat_result$names)

  climat_1 <- climat_data_1[49:201]

  
  
block_temp <- as.data.frame(cbind(d, climat_1))
lib_ccm <-  c(1,NROW(d))

out.temp <- do.call(rbind, lapply(1:8, function(E_i) {
  pred_ccm <- make_pred_nozero(block_temp$d, E_i)
  ccm(block = block_temp, 
      E = E_i,
      lib = lib_ccm,
      pred = pred_ccm,
      lib_sizes = NROW(block_temp),
      exclusion_radius = 0,
      random_libs = FALSE,
      num_samples = 1,
      tp = -1,
      lib_column = 1,
      target_column = 2)
}))

E_start <- out.temp$E[which.max(out.temp$rho)]

pred_ccm <- make_pred_nozero(block_temp$d,E_start)

data.result.ccm <-   ccm(block = block_temp, 
                         E = E_start,
                         lib = lib_ccm,
                         pred = pred_ccm,
                         lib_sizes = NROW(block_temp),
                         exclusion_radius = 0,
                         random_libs = FALSE,
                         num_samples = 1,
                         tp = 0)

data.result.ccm <- cbind(data.result.ccm, prov_c )

return(data.result.ccm)
  }




ccm_apply_surr <- function (num_prov, num_clim)
{
  # number of province in the list
  vari = c(1,2,3,4,5,6,7) 
  prov_c <- as.character(vn@data$VARNAME_2[as.numeric(num_prov)])
  vari_c <- variables[as.numeric(vari)]
  
  #take a coordinate of provinces 
  x <- mean(vn@polygons[[as.numeric(num_prov)]]@Polygons[[1]]@coords[,1]) 
  y <- mean(vn@polygons[[as.numeric(num_prov)]]@Polygons[[1]]@coords[,2])
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
  
  climat_data_1 <- select_variable(vari_c[num_clim],climat_result$names)

  climat_1 <- climat_data_1[49:201]

  block_temp <- as.data.frame(cbind(d, climat_1))
  lib_ccm <-  c(1,NROW(d))
  

  #  Calculer Null distribution avec surrogate
  out <- mean_anomal(climat_1)
  AH.bar <- out$mean
  AH.tilde <- out$anomaly
  clm.surrs <- do.call(cbind, lapply(1:500, function(i) 
    {
    I_na <- is.na(AH.tilde)
    out <- AH.bar
    out[I_na] <- NA
    out[!I_na] <- out[!I_na] + sample(AH.tilde[!I_na], sum(!I_na), replace = FALSE)
    return(out)
    } ))
  
  data.result.ccm.surr <- data.frame()
  
  for (i_surr in 1:dim(clm.surrs)[2])
    {
    block_temp <- data.frame(d,AH=clm.surrs[,i_surr])
    out.temp <- do.call( rbind, lapply(1:8, function(E_i)
      { 
      pred_ccm <- make_pred_nozero(block_temp$d,E_i) 
      ccm(block=block_temp,
          E=E_i,
          lib=lib_ccm,
          pred=pred_ccm,
          lib_sizes = NROW(block_temp),
          exclusion_radius=0,
          random_libs = FALSE,
          num_sample=1,
          tp = -1,
          lib_column = 1,
          target_column = 2)
      })) 
    
    E_star <- out.temp$E[which.max(out.temp$rho)]
    pred_ccm <- make_pred_nozero(block_temp$d,E_star)
    out.temp <- ccm(block=block_temp,
                    E=E_star,
                    lib=lib_ccm,
                    pred = pred_ccm,
                    lib_sizes = NROW(block_temp),
                    exclusion_radius=0, 
                    random_libs = FALSE, 
                    num_sample=1, 
                    tp = 0)
    
    data.result.ccm.surr <- data.result.ccm.surr %>% bind_rows(out.temp) 

  }
  
  return(data.result.ccm.surr)
}









boxplot(data.result.ccm.surr$rho, ylim=c(0,1), ylab='rho (CCM)', xlab='Denmark', main='Flu cross-map AH')
points(1,data.result.ccm$rho, pch=8, col='red', cex=2, lwd=2)

boxplot(result_surr[[1]]$rho, 
        result_surr[[2]]$rho, 
        result_surr[[3]]$rho,
        result_surr[[4]]$rho,
        main = "testing",
        at = c(1,2,3,4),
        names = c("11","22","33","44"),
        las = 2,
        horizontal = TRUE)




q0 = data.result.ccm$rho
qm = mean(data.result.ccm.surr$rho)
delq = sd(data.result.ccm.surr$rho)
s = (abs(q0 - qm)) / delq

data.result.ccm$rho






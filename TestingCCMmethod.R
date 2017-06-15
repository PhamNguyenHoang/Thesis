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
library(metap)



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

# maindata <- read.csv("DataCCM.txt",header = TRUE)
# data.vn <- maindata %>% filter(country == 'Japan') %>% 
#   filter(year >= 1996)
# 
# data.vn.trans <- data.vn %>% select(-country) %>% spread(variable,value) %>% mutate(date = ISOdate(year,month,day)) %>% select(-year,-month,-day) %>% select(date,flu,everything())

# set.seed(599213)
# 
# 
make_block <- function(data,cols,delays,lib=c(1,NROW(data)))
  {
  lib <- matrix(lib,ncol = 2) 
  data <- as.matrix(data)
  ncol <- length(cols) 
  nrow <- dim(data)[1] 
  
  block <- array(NA,dim = c(nrow,ncol)) 
  colnames(block) <- 1:ncol
  for (i in 1:ncol)
    { 
    I <- 1:nrow 
    I_delay <- intersect(I,I+delays[i]) 
    block[I_delay-delays[i],i] <- data[I_delay,cols[i]] 
    if (delays[i] < 0)
      { 
      # remove data points that fall at start of lib segments 
      block[lib[,1] - (0:(delays[i]+1)),i] <- NA 
      colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t-',abs(delays[i]),sep="")
      }
    else if (delays[i] > 0) 
      { # remove data points that fall at end of lib segments 
      block[lib[,2] - (0:(delays[i]+1)),i] <- NA
      colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t+',abs(delays[i]),sep="")
      }
    else 
      { 
        colnames(block)[i] <- paste(colnames(data)[cols[i]],'_t',sep="")
      } 
    }
  return(block)
}

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


make_sea_surr_data <- function(ts,num_surr,T_period)
{
  n <- length(ts)
  I_season <- suppressWarnings(matrix(1:T_period, nrow=n, ncol=1))
  
  # Calculate seasonal cycle using smooth.spline
  seasonal_F <- smooth.spline(c(I_season - T_period, I_season, I_season + T_period), 
                              c(ts, ts, ts))
  seasonal_cyc <- predict(seasonal_F,I_season)$y
  seasonal_resid <- ts - seasonal_cyc
  
  return(sapply(1:num_surr, function(i) {
    seasonal_cyc + sample(seasonal_resid, n)
  }))
}


tmp = 1:64
excl = c(26,27, 28, 34, 35)
tmp = tmp[!(tmp %in% excl)]

time_seq <- seq(ISOdate(1998,1,1), by = "month", length.out = 153)
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

##############


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

  # block_temp <- as.data.frame(cbind(d, climat_1))
  lib_ccm <-  c(1,NROW(d))
  

  #  Calculer Null distribution avec surrogate
  
  clm.surrs <- make_sea_surr_data(climat_1,num_surr,12)
#   plot(time_seq,climat_1, type = "l", col = "blue",xlab = "date", ylab= "Absolute Humidity")
#   lines(time_seq,tes.surr[,2],col = "red")
  
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

num_surr = 50
result <- lapply(tmp, function(x) ccm_apply(x,7))
result_surr <- lapply(tmp, function(x) ccm_apply_surr(x,7))

    
## Classify by lattitude
rs <- NULL
for(i in tmp)
{
    rs <- rbind(rs,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
}

rs <- as.data.frame(rs)
rs$rank <-  rank(as.numeric(as.character(rs$V2)), ties.method = "first")
rho <- NULL
for(j in seq(1,NROW(rs),by =  1) )
{
      rho <- rbind(rho,result[[as.numeric(j)]]$rho )
}
rs <- cbind(rs,rho)

rs$V1 <- as.character(rs$V1)
rs$V2 <- as.numeric(as.character(rs$V2))

sum <- data.frame(count = numeric(NROW(rs)), pvalue =  numeric(NROW(rs))   )


for (t in seq(1,NROW(rs),by =  1))
{
  sum$count[t] = 1
  for ( k in seq(1, num_surr, by = 1))
  {
    if(result_surr[[t]]$rho[k] > result[[t]]$rho)
      sum$count[t] <- sum$count[t] + 1
      
  }
}

for (t in seq(1,NROW(rs),by =  1))
{
  sum$pvalue[t] <- sum$count[t]/(num_surr+1)
}




# rs <- rbind(rs, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
# rs <- rbind(rs, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
# rs <- rbind(rs, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
# rs <- rbind(rs, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
# rs <- rbind(rs, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))
# 
# rs1 <- rs
# rs1$V2 <- NULL
# rs1$rank <- NULL
# rs1$rho <- as.numeric(rs1$rho)
## Plotting graph

# vn1 <- vn
# vn1@data <- merge(vn1@data, rs1, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)
# 
# spplot(vn1,13,col.regions=colorRampPalette(brewer.pal(9, "YlOrRd"))(10), at = seq(0,1, by = 0.1),
#        main = "Rainfall")

boxplot(result_surr[[1]]$rho, 
        result_surr[[2]]$rho, 
        result_surr[[3]]$rho,
        result_surr[[4]]$rho,
        result_surr[[5]]$rho,
        result_surr[[6]]$rho,
        result_surr[[7]]$rho,
        result_surr[[8]]$rho,
        result_surr[[9]]$rho,
        result_surr[[10]]$rho,
        result_surr[[11]]$rho,
        result_surr[[12]]$rho,
        result_surr[[13]]$rho,
        result_surr[[14]]$rho,
        result_surr[[15]]$rho,
        result_surr[[16]]$rho,
        result_surr[[17]]$rho,
        result_surr[[18]]$rho,
        result_surr[[19]]$rho,
        result_surr[[20]]$rho,
        result_surr[[21]]$rho,
        result_surr[[22]]$rho,
        result_surr[[23]]$rho,
        result_surr[[24]]$rho,
        result_surr[[25]]$rho,
        result_surr[[26]]$rho,
        result_surr[[27]]$rho,
        result_surr[[28]]$rho,
        result_surr[[29]]$rho,
        result_surr[[30]]$rho,
        result_surr[[31]]$rho,
        result_surr[[32]]$rho,
        result_surr[[33]]$rho,
        result_surr[[34]]$rho,
        result_surr[[35]]$rho,
        result_surr[[36]]$rho,
        result_surr[[37]]$rho,
        result_surr[[38]]$rho,
        result_surr[[39]]$rho,
        result_surr[[40]]$rho,
        result_surr[[41]]$rho,
        result_surr[[42]]$rho,
        result_surr[[43]]$rho,
        result_surr[[44]]$rho,
        result_surr[[45]]$rho,
        result_surr[[46]]$rho,
        result_surr[[47]]$rho,
        result_surr[[48]]$rho,
        result_surr[[49]]$rho,
        result_surr[[50]]$rho,
        result_surr[[51]]$rho,
        result_surr[[52]]$rho,
        result_surr[[53]]$rho,
        result_surr[[54]]$rho,
        result_surr[[55]]$rho,
        result_surr[[56]]$rho,
        result_surr[[57]]$rho,
        result_surr[[58]]$rho,
        result_surr[[59]]$rho,
        ylim = c(0,1),
        main = "Absolute Humidity",
        at = rs$rank,
        names = rs$V1,
        las = 2,
        horizontal = TRUE,
        cex = 0.5, cex.axis = 0.5)

for(i in seq(1,59, by = 1))
{
  if(sum$pvalue[i] < 0.05)
  {
    points(rs$rho[i],rs$rank[i], pch = 19, col = "red",cex = 0.8) 
  }
         
    else
    {
      points(rs$rho[i],rs$rank[i], pch = 21, col = "red",cex = 0.8)   
    }
}


######## plotting rho map #########
result1 <- lapply(tmp, function(x) ccm_apply(x,1))
result2 <- lapply(tmp, function(x) ccm_apply(x,2))
result3 <- lapply(tmp, function(x) ccm_apply(x,3))
result4 <- lapply(tmp, function(x) ccm_apply(x,4))
result5 <- lapply(tmp, function(x) ccm_apply(x,5))
result6 <- lapply(tmp, function(x) ccm_apply(x,6))
result7 <- lapply(tmp, function(x) ccm_apply(x,7))

rs1 <- NULL
rs2 <- NULL
rs3 <- NULL
rs4 <- NULL
rs5 <- NULL
rs6 <- NULL
rs7 <- NULL
rho1 <- NULL
rho2 <- NULL
rho3 <- NULL
rho4 <- NULL
rho5 <- NULL
rho6 <- NULL
rho7 <- NULL

for(i in tmp)
{
  rs1 <- rbind(rs1,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
  rs2 <- rbind(rs2,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
  rs3 <- rbind(rs3,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
  rs4 <- rbind(rs4,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
  rs5 <- rbind(rs5,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
  rs6 <- rbind(rs6,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
  rs7 <- rbind(rs7,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
}

rs1 <- as.data.frame(rs1)
rs2 <- as.data.frame(rs2)
rs3 <- as.data.frame(rs3)
rs4 <- as.data.frame(rs4)
rs5 <- as.data.frame(rs5)
rs6 <- as.data.frame(rs6)
rs7 <- as.data.frame(rs7)
rs1$rank <-  rank(as.numeric(as.character(rs1$V2)), ties.method = "first")
rs2$rank <-  rank(as.numeric(as.character(rs2$V2)), ties.method = "first")
rs3$rank <-  rank(as.numeric(as.character(rs3$V2)), ties.method = "first")
rs4$rank <-  rank(as.numeric(as.character(rs4$V2)), ties.method = "first")
rs5$rank <-  rank(as.numeric(as.character(rs5$V2)), ties.method = "first")
rs6$rank <-  rank(as.numeric(as.character(rs6$V2)), ties.method = "first")
rs7$rank <-  rank(as.numeric(as.character(rs7$V2)), ties.method = "first")

for(j in seq(1,NROW(rs1),by =  1) )
{
  if(result1[[as.numeric(j)]]$rho < 0)
    rho1 <- rbind(rho1, 0)
  else 
  {
    rho1 <- rbind(rho1,result1[[as.numeric(j)]]$rho )
  }
  
  if(result2[[as.numeric(j)]]$rho < 0)
    rho2 <- rbind(rho2, 0)
  else 
  {
    rho2 <- rbind(rho2,result2[[as.numeric(j)]]$rho )
  }
  if(result3[[as.numeric(j)]]$rho < 0)
    rho3 <- rbind(rho3, 0)
  else 
  {
    rho3 <- rbind(rho3,result3[[as.numeric(j)]]$rho )
  }
  if(result4[[as.numeric(j)]]$rho < 0)
    rho4 <- rbind(rho4, 0)
  else 
  {
    rho4 <- rbind(rho4,result4[[as.numeric(j)]]$rho )
  }
  if(result5[[as.numeric(j)]]$rho < 0)
    rho5 <- rbind(rho5, 0)
  else 
  {
    rho5 <- rbind(rho5,result5[[as.numeric(j)]]$rho )
  }
  if(result6[[as.numeric(j)]]$rho < 0)
    rho6 <- rbind(rho6, 0)
  else 
  {
    rho6 <- rbind(rho6,result6[[as.numeric(j)]]$rho )
  }
  if(result7[[as.numeric(j)]]$rho < 0)
    rho7 <- rbind(rho7, 0)
  else 
  {
    rho7 <- rbind(rho7,result7[[as.numeric(j)]]$rho )
  }
}
rs1 <- cbind(rs1,rho1)
rs2 <- cbind(rs2,rho2)
rs3 <- cbind(rs3,rho3)
rs4 <- cbind(rs4,rho4)
rs5 <- cbind(rs5,rho5)
rs6 <- cbind(rs6,rho6)
rs7 <- cbind(rs7,rho7)




rs1$V1 <- as.character(rs1$V1)
rs1$V2 <- as.numeric(as.character(rs1$V2))
rs2$V1 <- as.character(rs2$V1)
rs2$V2 <- as.numeric(as.character(rs2$V2))
rs3$V1 <- as.character(rs3$V1)
rs3$V2 <- as.numeric(as.character(rs3$V2))
rs4$V1 <- as.character(rs4$V1)
rs4$V2 <- as.numeric(as.character(rs4$V2))
rs5$V1 <- as.character(rs5$V1)
rs5$V2 <- as.numeric(as.character(rs5$V2))
rs6$V1 <- as.character(rs6$V1)
rs6$V2 <- as.numeric(as.character(rs6$V2))
rs7$V1 <- as.character(rs7$V1)
rs7$V2 <- as.numeric(as.character(rs7$V2))


rs1 <- rbind(rs1, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
rs1 <- rbind(rs1, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
rs1 <- rbind(rs1, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
rs1 <- rbind(rs1, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
rs1 <- rbind(rs1, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))

rs2 <- rbind(rs2, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
rs2 <- rbind(rs2, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
rs2 <- rbind(rs2, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
rs2 <- rbind(rs2, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
rs2 <- rbind(rs2, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))

rs3 <- rbind(rs3, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
rs3 <- rbind(rs3, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
rs3 <- rbind(rs3, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
rs3 <- rbind(rs3, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
rs3 <- rbind(rs3, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))

rs4 <- rbind(rs4, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
rs4 <- rbind(rs4, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
rs4 <- rbind(rs4, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
rs4 <- rbind(rs4, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
rs4 <- rbind(rs4, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))

rs5 <- rbind(rs5, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
rs5 <- rbind(rs5, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
rs5 <- rbind(rs5, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
rs5 <- rbind(rs5, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
rs5 <- rbind(rs5, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))

rs6 <- rbind(rs6, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
rs6 <- rbind(rs6, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
rs6 <- rbind(rs6, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
rs6 <- rbind(rs6, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
rs6 <- rbind(rs6, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))

rs7 <- rbind(rs7, c(vn@data$VARNAME_2[as.numeric(26)], 0, "NA", "NA" ))
rs7 <- rbind(rs7, c(vn@data$VARNAME_2[as.numeric(28)], 0, "NA", "NA" ))
rs7 <- rbind(rs7, c(vn@data$VARNAME_2[as.numeric(34)], 0, "NA", "NA" ))
rs7 <- rbind(rs7, c(vn@data$VARNAME_2[as.numeric(35)], 0, "NA", "NA" ))
rs7 <- rbind(rs7, c(vn@data$VARNAME_2[as.numeric(27)], 0, "NA", "NA" ))

rs1$V2 <- NULL
rs1$rank <- NULL
rs1$rho1 <- as.numeric(rs1$rho1)

rs2$V2 <- NULL
rs2$rank <- NULL
rs2$rho2 <- as.numeric(rs2$rho2)

rs3$V2 <- NULL
rs3$rank <- NULL
rs3$rho3 <- as.numeric(rs3$rho3)

rs4$V2 <- NULL
rs4$rank <- NULL
rs4$rho4 <- as.numeric(rs4$rho4)

rs5$V2 <- NULL
rs5$rank <- NULL
rs5$rho5 <- as.numeric(rs5$rho5)

rs6$V2 <- NULL
rs6$rank <- NULL
rs6$rho6 <- as.numeric(rs6$rho6)

rs7$V2 <- NULL
rs7$rank <- NULL
rs7$rho7 <- as.numeric(rs7$rho7)
## Plotting graph

vn1 <- vn
vn1@data <- merge(vn1@data, rs1, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)
vn1@data <- merge(vn1@data, rs2, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)
vn1@data <- merge(vn1@data, rs3, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)
vn1@data <- merge(vn1@data, rs4, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)
vn1@data <- merge(vn1@data, rs5, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)
vn1@data <- merge(vn1@data, rs6, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)
vn1@data <- merge(vn1@data, rs7, by.x = "VARNAME_2", by.y ="V1",sort = FALSE)

spplot(vn1,c(13,14,15,16,17,18,19),col.regions=colorRampPalette(brewer.pal(9, "YlOrRd"))(10), 
       names.attr = c("Average Temperature","Max Temperature","Min Temperature","Absolute Humidity","Rainfall","Relative Humidity","Hours of Sunshine"),
       at = seq(0,1, by = 0.1),
       colorkey=TRUE, scales = list(draw = FALSE),
       main = "The Result of Convergence Cross Mapping Test")




spplot(vn1,c(13,14,15,16,17,18,19),names.attr = c("Average Temperature","Max Temperature","Min Temperature","Absolute Humidity","Rainfall","Relative Humidity","Hours of Sunshine"),
       col.regions=colors_10, #at = seq(-0.5, 1.5, by = 1), 
       colorkey=FALSE, scales = list(draw = FALSE),
       main = "The Result of Granger Causality Test",
       sub = "Null Hypothesis : Climatic factor for each province do not Granger-cause the DHF in that provinces")



############################################################################################
## Multivariable EDM


multivariableEDM <- function(num_prov, num_clim)
{
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
  
  out.temp <- do.call( rbind, lapply(1:8, function(E_i)
  { 
    pred_mfi <- make_pred_nozero(d,E_i) 
    simplex(d, E=E_i, pred=pred_mfi, tp = 2)
    
  })) 
  
  E_star <- out.temp$E[which.max(out.temp$rho)[1]]
  
  block_in <- NULL
  block_in <- cbind(block_in, d, climat_1)
  
  block_mfi <- make_block(block_in,c(rep(1,E_star),2) , c(tp=2,0:-(E_star- 2),0)) 
  block_mfi <- as.data.frame(block_mfi)
  pred_mfi <- make_pred_nozero(block_mfi$d_t,E_star)
  for(j in 1:NCOL(block_mfi)) 
    block_mfi[,j] <- (block_mfi[,j] - mean(block_mfi[,j], na.rm = TRUE)) / sd(block_mfi[,j], na.rm = TRUE)
  
  out.multivar <- block_lnlp(block = block_mfi, method = 'simplex', num_neighbors=E_star+2, pred = pred_mfi, tp=0, target_column=1, columns=2:(E_star+1))
  
  out.univar <- block_lnlp(block = block_mfi, method = 'simplex', num_neighbors=E_star+2, pred = pred_mfi, tp=0, target_column=1, columns=2:(E_star))
  
  delta_rho <- out.multivar$rho - out.univar$rho
  
  
  return(delta_rho)
}

multivariableEDM2 <- function(num_prov, num_clim1,num_clim2)
{
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
  
  climat_data_1 <- select_variable(vari_c[num_clim1],climat_result$names)
  climat_data_2 <- select_variable(vari_c[num_clim2],climat_result$names)
  climat_1 <- climat_data_1[49:201]
  climat_2 <- climat_data_2[49:201]
  
  out.temp <- do.call( rbind, lapply(1:8, function(E_i)
  { 
    pred_mfi <- make_pred_nozero(d,E_i) 
    simplex(d, E=E_i, pred=pred_mfi, tp = 2)
    
  })) 
  
  E_star <- out.temp$E[which.max(out.temp$rho)[1]]
  
  block_in <- NULL
  block_in <- cbind(block_in, d, climat_1, climat_2)
  
  block_mfi <- make_block(block_in,c(rep(1,E_star),2,3) , c(tp=2,0:-(E_star- 2),0,0)) 
  block_mfi <- as.data.frame(block_mfi)
  pred_mfi <- make_pred_nozero(block_mfi$d_t,E_star)
  for(j in 1:NCOL(block_mfi)) 
    block_mfi[,j] <- (block_mfi[,j] - mean(block_mfi[,j], na.rm = TRUE)) / sd(block_mfi[,j], na.rm = TRUE)
  
  out.multivar <- block_lnlp(block = block_mfi, method = 'simplex', num_neighbors=E_star+2, pred = pred_mfi, tp=0, target_column=1, columns=2:(E_star+2))
  
  out.univar <- block_lnlp(block = block_mfi, method = 'simplex', num_neighbors=E_star+2, pred = pred_mfi, tp=0, target_column=1, columns=2:(E_star))
  
  delta_rho <- out.multivar$rho - out.univar$rho
  
  
  return(delta_rho)
}





tmp1 = 1:64
excl1 = c(26,27, 28, 34, 35, 58)
tmp1 = tmp1[!(tmp1 %in% excl1)]

result_multi1 <- sapply(tmp1, function(x) multivariableEDM(x, 1))
result_multi2 <- sapply(tmp1, function(x) multivariableEDM(x, 2))
result_multi3 <- sapply(tmp1, function(x) multivariableEDM(x, 3))
result_multi4 <- sapply(tmp1, function(x) multivariableEDM(x, 4))
result_multi5 <- sapply(tmp1, function(x) multivariableEDM(x, 5))
result_multi6 <- sapply(tmp1, function(x) multivariableEDM(x, 6))
result_multi7 <- sapply(tmp1, function(x) multivariableEDM(x, 7))

result_multi8 <- sapply(tmp1, function(x) multivariableEDM2(x, 1, 4))
result_multi9 <- sapply(tmp1, function(x) multivariableEDM2(x, 1, 5))
result_multi10 <- sapply(tmp1, function(x) multivariableEDM2(x, 1, 6))
result_multi11 <- sapply(tmp1, function(x) multivariableEDM2(x, 1, 7))
result_multi12 <- sapply(tmp1, function(x) multivariableEDM2(x, 4, 5))
result_multi13 <- sapply(tmp1, function(x) multivariableEDM2(x, 4, 6))
result_multi14 <- sapply(tmp1, function(x) multivariableEDM2(x, 4, 7))

result_multi15 <- sapply(tmp1, function(x) multivariableEDM2(x, 5, 6))
result_multi16 <- sapply(tmp1, function(x) multivariableEDM2(x, 5, 7))
result_multi17 <- sapply(tmp1, function(x) multivariableEDM2(x, 6, 7))


rs <- NULL
for(i in tmp1)
{
  rs <- rbind(rs,c(as.character(vn@data$VARNAME_2[as.numeric(i)]),mean(vn@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2]) ))
}

rs <- as.data.frame(rs)
rs$rank <-  rank(as.numeric(as.character(rs$V2)), ties.method = "first")

result_multi1_n <- NULL
result_multi1_c <- NULL
result_multi1_s <- NULL

result_multi2_n <- NULL
result_multi2_c <- NULL
result_multi2_s <- NULL

result_multi3_n <- NULL
result_multi3_c <- NULL
result_multi3_s <- NULL

result_multi4_n <- NULL
result_multi4_c <- NULL
result_multi4_s <- NULL

result_multi5_n <- NULL
result_multi5_c <- NULL
result_multi5_s <- NULL

result_multi6_n <- NULL
result_multi6_c <- NULL
result_multi6_s <- NULL

result_multi7_n <- NULL
result_multi7_c <- NULL
result_multi7_s <- NULL

result_multi8_n <- NULL
result_multi8_c <- NULL
result_multi8_s <- NULL

result_multi9_n <- NULL
result_multi9_c <- NULL
result_multi9_s <- NULL

result_multi10_n <- NULL
result_multi10_c <- NULL
result_multi10_s <- NULL

result_multi11_n <- NULL
result_multi11_c <- NULL
result_multi11_s <- NULL

result_multi12_n <- NULL
result_multi12_c <- NULL
result_multi12_s <- NULL

result_multi13_n <- NULL
result_multi13_c <- NULL
result_multi13_s <- NULL

result_multi14_n <- NULL
result_multi14_c <- NULL
result_multi14_s <- NULL

result_multi15_n <- NULL
result_multi15_c <- NULL
result_multi15_s <- NULL

result_multi16_n <- NULL
result_multi16_c <- NULL
result_multi16_s <- NULL

result_multi17_n <- NULL
result_multi17_c <- NULL
result_multi17_s <- NULL

result_multi1 <- as.data.frame(cbind(result_multi1, rs$rank))
result_multi1_s <- result_multi1 %>% filter(V2 < 25)
result_multi1_c <- result_multi1 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi1_n <- result_multi1 %>% filter(V2 >= 40)

result_multi2 <- as.data.frame(cbind(result_multi2, rs$rank))
result_multi2_s <- result_multi2 %>% filter(V2 < 25)
result_multi2_c <- result_multi2 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi2_n <- result_multi2 %>% filter(V2 >= 40)

result_multi3 <- as.data.frame(cbind(result_multi3, rs$rank))
result_multi3_s <- result_multi3 %>% filter(V2 < 25)
result_multi3_c <- result_multi3 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi3_n <- result_multi3 %>% filter(V2 >= 40)

result_multi4 <- as.data.frame(cbind(result_multi4, rs$rank))
result_multi4_s <- result_multi4 %>% filter(V2 < 25)
result_multi4_c <- result_multi4 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi4_n <- result_multi4 %>% filter(V2 >= 40)

result_multi5 <- as.data.frame(cbind(result_multi5, rs$rank))
result_multi5_s <- result_multi5 %>% filter(V2 < 25)
result_multi5_c <- result_multi5 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi5_n <- result_multi5 %>% filter(V2 >= 40)

result_multi6 <- as.data.frame(cbind(result_multi6, rs$rank))
result_multi6_s <- result_multi6 %>% filter(V2 < 25)
result_multi6_c <- result_multi6 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi6_n <- result_multi6 %>% filter(V2 >= 40)

result_multi7 <- as.data.frame(cbind(result_multi7, rs$rank))
result_multi7_s <- result_multi7 %>% filter(V2 < 25)
result_multi7_c <- result_multi7 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi7_n <- result_multi7 %>% filter(V2 >= 40)

result_multi8 <- as.data.frame(cbind(result_multi8, rs$rank))
result_multi8_s <- result_multi8 %>% filter(V2 < 25)
result_multi8_c <- result_multi8 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi8_n <- result_multi8 %>% filter(V2 >= 40)

result_multi9 <- as.data.frame(cbind(result_multi9, rs$rank))
result_multi9_s <- result_multi9 %>% filter(V2 < 25)
result_multi9_c <- result_multi9 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi9_n <- result_multi9 %>% filter(V2 >= 40)

result_multi10 <- as.data.frame(cbind(result_multi10, rs$rank))
result_multi10_s <- result_multi10 %>% filter(V2 < 25)
result_multi10_c <- result_multi10 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi10_n <- result_multi10 %>% filter(V2 >= 40)

result_multi11 <- as.data.frame(cbind(result_multi11, rs$rank))
result_multi11_s <- result_multi11 %>% filter(V2 < 25)
result_multi11_c <- result_multi11 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi11_n <- result_multi11 %>% filter(V2 >= 40)

result_multi12 <- as.data.frame(cbind(result_multi12, rs$rank))
result_multi12_s <- result_multi12 %>% filter(V2 < 25)
result_multi12_c <- result_multi12 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi12_n <- result_multi12 %>% filter(V2 >= 40)

result_multi13 <- as.data.frame(cbind(result_multi13, rs$rank))
result_multi13_s <- result_multi13 %>% filter(V2 < 25)
result_multi13_c <- result_multi13 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi13_n <- result_multi13 %>% filter(V2 >= 40)

result_multi14 <- as.data.frame(cbind(result_multi14, rs$rank))
result_multi14_s <- result_multi14 %>% filter(V2 < 25)
result_multi14_c <- result_multi14 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi14_n <- result_multi14 %>% filter(V2 >= 40)

result_multi15 <- as.data.frame(cbind(result_multi15, rs$rank))
result_multi15_s <- result_multi15 %>% filter(V2 < 25)
result_multi15_c <- result_multi15 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi15_n <- result_multi15 %>% filter(V2 >= 40)

result_multi16 <- as.data.frame(cbind(result_multi16, rs$rank))
result_multi16_s <- result_multi16 %>% filter(V2 < 25)
result_multi16_c <- result_multi16 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi16_n <- result_multi16 %>% filter(V2 >= 40)

result_multi17 <- as.data.frame(cbind(result_multi17, rs$rank))
result_multi17_s <- result_multi17 %>% filter(V2 < 25)
result_multi17_c <- result_multi17 %>% filter(25 <= V2) %>% filter(V2 < 40)
result_multi17_n <- result_multi17 %>% filter(V2 >= 40)

opar <- par(no.readonly=TRUE)
par(mfrow=c(1,4))
boxplot(result_multi1[,1], 
        result_multi2[,1], 
        result_multi3[,1],
        result_multi4[,1],
        result_multi5[,1],
        result_multi6[,1],
        result_multi7[,1],
        result_multi8[,1],
        result_multi9[,1],
        result_multi10[,1],
        result_multi11[,1],
        result_multi12[,1],
        result_multi13[,1],
        result_multi14[,1],
        result_multi15[,1],
        result_multi16[,1],
        result_multi17[,1],
        ylim = c(-0.6, 0.6),
        main = "Global of Vietnam",
        at = 1:17,
        names = c("ta","tx","tm","rf","rh","ah","sh","ta-rf","ta-rh","ta-sh","ta-ah","rf-rh","rf-sh","rf-ah","rh-sh","rh-ah","sh-ah"),
        las = 1,
        horizontal = FALSE,
        cex = 0.5, cex.axis = 0.5)

boxplot(result_multi1_n[,1], 
        result_multi2_n[,1], 
        result_multi3_n[,1],
        result_multi4_n[,1],
        result_multi5_n[,1],
        result_multi6_n[,1],
        result_multi7_n[,1],
        result_multi8_n[,1],
        result_multi9_n[,1],
        result_multi10_n[,1],
        result_multi11_n[,1],
        result_multi12_n[,1],
        result_multi13_n[,1],
        result_multi14_n[,1],
        result_multi15_n[,1],
        result_multi16_n[,1],
        result_multi17_n[,1],
        ylim = c(-0.6, 0.6),
        main = "North of Vietnam",
        at = 1:17,
        names = c("ta","tx","tm","rf","rh","ah","sh","ta-rf","ta-rh","ta-sh","ta-ah","rf-rh","rf-sh","rf-ah","rh-sh","rh-ah","sh-ah"),
        las = 1,
        horizontal = FALSE,
        cex = 0.5, cex.axis = 0.5)

boxplot(result_multi1_c[,1], 
        result_multi2_c[,1], 
        result_multi3_c[,1],
        result_multi4_c[,1],
        result_multi5_c[,1],
        result_multi6_c[,1],
        result_multi7_c[,1],
        result_multi8_c[,1],
        result_multi9_c[,1],
        result_multi10_c[,1],
        result_multi11_c[,1],
        result_multi12_c[,1],
        result_multi13_c[,1],
        result_multi14_c[,1],
        result_multi15_c[,1],
        result_multi16_c[,1],
        result_multi17_c[,1],
        ylim = c(-0.6, 0.6),
        main = "Centre of Vietnam",
        at = 1:17,
        names = c("ta","tx","tm","rf","rh","ah","sh","ta-rf","ta-rh","ta-sh","ta-ah","rf-rh","rf-sh","rf-ah","rh-sh","rh-ah","sh-ah"),
        las = 1,
        horizontal = FALSE,
        cex = 0.5, cex.axis = 0.5)

boxplot(result_multi1_s[,1], 
        result_multi2_s[,1], 
        result_multi3_s[,1],
        result_multi4_s[,1],
        result_multi5_s[,1],
        result_multi6_s[,1],
        result_multi7_s[,1],
        result_multi8_s[,1],
        result_multi9_s[,1],
        result_multi10_s[,1],
        result_multi11_s[,1],
        result_multi12_s[,1],
        result_multi13_s[,1],
        result_multi14_s[,1],
        result_multi15_s[,1],
        result_multi16_s[,1],
        result_multi17_s[,1],
        ylim = c(-0.6, 0.6),
        main = "South of Vietnam",
        at = 1:17,
        names = c("ta","tx","tm","rf","rh","ah","sh","ta-rf","ta-rh","ta-sh","ta-ah","rf-rh","rf-sh","rf-ah","rh-sh","rh-ah","sh-ah"),
        las = 1,
        horizontal = FALSE,
        cex = 0.5, cex.axis = 0.5)


par(opar)








boxplot(result_multi1[,1], 
        result_multi2[,1], 
        result_multi3[,1],
        result_multi4[,1],
        result_multi5[,1],
        result_multi6[,1],
        result_multi7[,1],
        result_multi8[,1],
        result_multi9[,1],
        result_multi10[,1],
        result_multi11[,1],
        result_multi12[,1],
        result_multi13[,1],
        result_multi14[,1],
        result_multi15[,1],
        result_multi16[,1],
        result_multi17[,1],
        ylim = c(-0.6, 0.6),
        main = "Multivariable EDM - Global of Vietnam",
        at = 1:17,
        names = c("ta","tx","tm","rf","rh","ah","sh","ta-rf","ta-rh","ta-sh","ta-ah","rf-rh","rf-sh","rf-ah","rh-sh","rh-ah","sh-ah"),
        las = 1,
        horizontal = FALSE,
        cex = 0.5, cex.axis = 0.5)



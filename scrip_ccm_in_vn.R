#####################################################################################
## Cette Script appliqué la méthode Convergence Cross Mapping on dengue data in VN ##
##                        Author : Pham Nguyen Hoang                               ##
#####################################################################################

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
library(multispatialCCM)

load("climatic_data_Vietnam_Laos_Thailand.rdata")
load("demo_meteo.RData")
load("climatic_data.RData")
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

tmp = 1:64
excl = c(26,27, 28, 34, 35)
tmp = tmp[!(tmp %in% excl)]
seq_dengue <- NULL
seq_climat_1 <- NULL
seq_climat_2 <- NULL
seq_climat_3 <- NULL
seq_climat_4 <- NULL
seq_climat_5 <- NULL
seq_climat_6 <- NULL
seq_climat_7 <- NULL

result_1 <- NULL
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
  climat_1 <- climat_data_1[49:201]
  climat_2 <- climat_data_2[49:201]
  climat_3 <- climat_data_3[49:201]
  climat_4 <- climat_data_4[49:201]
  climat_5 <- climat_data_5[49:201]
  climat_6 <- climat_data_6[49:201]
  climat_7 <- climat_data_7[49:201]
#   d <- as.data.frame(d)
#   climat_1 <- as.data.frame(climat_1)
#   climat_2 <- as.data.frame(climat_2)
#   climat_3 <- as.data.frame(climat_3)
#   climat_4 <- as.data.frame(climat_4)
#   climat_5 <- as.data.frame(climat_5)
#   climat_6 <- as.data.frame(climat_6)
#   climat_7 <- as.data.frame(climat_7)
#   
#   
#   d <- rbind(d,NA)
#   climat_1 <- rbind(climat_1,NA)
#   climat_2 <- rbind(climat_2,NA)
#   climat_3 <- rbind(climat_3,NA)
#   climat_4 <- rbind(climat_4,NA)
#   climat_5 <- rbind(climat_5,NA)
#   climat_6 <- rbind(climat_6,NA)
#   climat_7 <- rbind(climat_7,NA)
#   
#   seq_dengue <- rbind(seq_dengue,d)
#   seq_climat_1 <- rbind(seq_climat_1,climat_1)
#   seq_climat_2 <- rbind(seq_climat_2,climat_2)
#   seq_climat_3 <- rbind(seq_climat_3,climat_3)
#   seq_climat_4 <- rbind(seq_climat_4,climat_4)
#   seq_climat_5 <- rbind(seq_climat_5,climat_5)
#   seq_climat_6 <- rbind(seq_climat_6,climat_6)
#   seq_climat_7 <- rbind(seq_climat_7,climat_7)

#}

# Accm <- seq_dengue$d
# Bccm <- seq_climat_1$climat_1
  
Accm <- d
Bccm <- climat_1
  
# maxE<-5 #Maximum E to test
#Matrix for storing output
# Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("A", "B")
# for(E in 2:maxE) {
#   #Uses defaults of looking forward one prediction step (predstep)
#   #And using time lag intervals of one time step (tau)
#   Emat[E-1,"A"]<-SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho
#   Emat[E-1,"B"]<-SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
# }
# 
# matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
#         xlab="E", ylab="rho", lwd=2)
# legend("bottomleft", c("A", "B"), lty=1:2, col=1:2, lwd=2, bty="n")
#Results will vary depending on simulation.
#Using the seed we provide,
#maximum E for A should be 2, and maximum E for B should be 3.
#For the analyses in the paper, we use E=2 for all simulations.
E_A<-2
E_B<-2

#increasing time distance
#See manuscript and R code for details
# signal_A_out<-SSR_check_signal(A=Accm, E=E_A, tau=1,
#                                predsteplist=1:10)
# signal_B_out<-SSR_check_signal(A=Bccm, E=E_B, tau=1,
#                                predsteplist=1:10)
#Run the CCM test
#E_A and E_B are the embedding dimensions for A and B.
#tau is the length of time steps used (default is 1)
#iterations is the number of bootsrap iterations (default 100)
# Does A "cause" B?
#Note - increase iterations to 100 for consistant results
CCM_boot_A<-CCM_boot(Accm, Bccm, E_A, tau=1, iterations=10)
# Does B "cause" A?
CCM_boot_B<-CCM_boot(Bccm, Accm, E_B, tau=1, iterations=10)
#Test for significant causal signal
#See R function for details
(CCM_significance_test<-ccmtest(CCM_boot_A,
                                CCM_boot_B))

result_1 <- rbind(result_1,c(prov_c, CCM_significance_test))


}
#Plot results
plotxlimits<-range(c(CCM_boot_A$Lobs, CCM_boot_B$Lobs))
#Plot "A causes B"
plot(CCM_boot_A$Lobs, CCM_boot_A$rho, type="l", col=1, lwd=2,
     xlim=c(plotxlimits[1], plotxlimits[2]), ylim=c(0,1),
     xlab="L", ylab="rho")
#Add +/- 1 standard error
matlines(CCM_boot_A$Lobs,
         cbind(CCM_boot_A$rho-CCM_boot_A$sdevrho,
               CCM_boot_A$rho+CCM_boot_A$sdevrho),
         lty=3, col=1)
#Plot "B causes A"
lines(CCM_boot_B$Lobs, CCM_boot_B$rho, type="l", col=2, lty=2, lwd=2)
#Add +/- 1 standard error
matlines(CCM_boot_B$Lobs,
         cbind(CCM_boot_B$rho-CCM_boot_B$sdevrho,
               CCM_boot_B$rho+CCM_boot_B$sdevrho),
         lty=3, col=2)
legend("topleft",
       c("A causes B", "B causes A"),
       lty=c(1,2), col=c(1,2), lwd=2, bty="n")






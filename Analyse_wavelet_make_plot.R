
library(sp)
library(biwavelet)
library(TSclust)
library(TSdist)
library(fpc)
library(rgeos)
library(vars)


load("climatic_data_Vietnam_Laos_Thailand.rdata")
load("demo_meteo.RData")
main_data <- read.csv("main_data_for_kmeans.csv", row.names = 1)
load("VNM_adm2.RData", vn <- new.env())
vn <- vn$gadm

provinces_vn <- vn@data$VARNAME_2

dengue_data_vn <- main_data[rownames(main_data) %in% provinces_vn,]

dengue_data_vn <- as.data.frame(t(dengue_data_vn))
### Meteo data
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


## calculate a euclien distance between two matrix
cal_dis_eucl <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, 
       xlab="Fréquence (Hz)", ylab="Force", main = "Représentation de la transformé Fourier du signal",
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}


num <- 41
vari <- 1


prov_c <- vn@data$VARNAME_2[as.numeric(num)]
vari_c <- variables[as.numeric(vari)]

cat(paste("Vous avez choisir le province ", prov_c, " et le variable ", vari_c,sep = ""))
#### Calculate a distance between selecting station and other station for finding a closer stations
x <- mean(vn@polygons[[as.numeric(num)]]@Polygons[[1]]@coords[,1]) 
y <- mean(vn@polygons[[as.numeric(num)]]@Polygons[[1]]@coords[,2])
prov_coor <- cbind(x,y)
prov_coor <- as.data.frame(prov_coor)
colnames(prov_coor) <- c("longitude","latitude")
prov_coor <-  as.matrix(prov_coor)
dis_eu <- NULL
for ( j in 1:63)
{
  x1 = cbind(mean(clim_ts[[j]]$longitude),mean(clim_ts[[j]]$latitude))
  
  a <- cal_dis_eucl(prov_coor,x1)
  dis_eu <- rbind(dis_eu,c(a, attributes(clim_ts[j])))
}
dis_eu <- as.data.frame(dis_eu)
dis_eu$V1 <- as.numeric(dis_eu$V1)
prov_result <- dis_eu[ which(dis_eu[,1] == min(dis_eu[,1]) ) , ]
cat(paste("\nLa province la plus proche de ",prov_c, " dans la donn?? climatique est ", prov_result$names, sep = "") )

year = seq(1994,2011, by = (17/204))
time_serie <-  year[49:201]
dengue_data_vn_traiter <- dengue_data_vn[,which( colnames(dengue_data_vn) == prov_c )  ]
variable_data <- select_variable(vari_c,prov_result$names)
variable_data_traiter <- variable_data[49:201]
dengue_data_vn_traiter_n <- (dengue_data_vn_traiter - mean(dengue_data_vn_traiter)) / max(abs( (dengue_data_vn_traiter - mean(dengue_data_vn_traiter) )))
variable_data_traiter_n <- (variable_data_traiter - mean(variable_data_traiter)) / max(abs((variable_data_traiter - mean(variable_data_traiter))))

dengue_data_frame = cbind(time_serie,dengue_data_vn_traiter_n)
variable_data_frame = cbind(time_serie,variable_data_traiter_n)
province_coherence = wtc(dengue_data_frame,variable_data_frame,quiet = TRUE)


# plot( c(1998,2011), c(-1,1), type = "n", main = "Analyzed Signals")
# lines(time_serie,dengue_data_vn_traiter_n,col="red")
# lines(time_serie,variable_data_traiter_n,col="blue")
# legend(2008,0.95,c(prov_c,vari_c),lty = c(1,1), lwd=c(2.5,2.5), col = c("red", "blue") )
# 
# dengue_data_vn_ff <- fft(dengue_data_vn_traiter_n)
# variable_data_ff <- fft(variable_data_traiter_n)


# # plot time series and their fourier transform
# # plot dengue times series 
# layout(matrix(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4), 2,10, byrow = TRUE))
# 
# plot(x = time_serie, y = dengue_data_vn_traiter_n,type = "l", xlim = c(1998,2011) , ylim = c(-1,1),
#      xlab = "Times", ylab = "Dengue infection rate", main = "Taux d'infection du province Cantho")
# 
# plot.frequency.spectrum(dengue_data_vn_ff,xlimits = c(0,20))
# 
# plot(x = time_serie, y = variable_data_traiter_n,type = "l", xlim = c(1998,2011) , ylim = c(-1,1),
#      xlab = "Times", ylab = "Température moyenne", main = "Température moyenne du province Cantho")
# 
# plot.frequency.spectrum(variable_data_ff,xlimits = c(0,20))
# 

#################################################################
#plot wavelet coherence. 

layout(matrix(c(1,1,1,1,2,2,2,2),4 ,2, byrow = FALSE))
plot( c(1998,2011), c(-1,1), type = "n", main = "Les signaux observés", xlab = "Temps", ylab= "Taux d'inflection", sub = "(a)")

lines(time_serie,variable_data_traiter_n,col="blue")
#legend("topleft",c(prov_c,vari_c),lty = c(1,1), lwd=c(1,1), col = c("red", "blue") )

#par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
plot(province_coherence,plot.cb = FALSE, plot.phase = TRUE, main = " Cohérence des ondelettes", xlab = "Temps", ylab = "Période" ,sub = "(b)")



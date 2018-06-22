################################################################
# Script appliquer le k-means sur les séries temporelles.
# Author : Pham Nguyen Hoang
# Date : 20/5/2018
################################################################

library(sp)
library(TSclust)
library(TSdist)
library(fpc)
library(rgeos)
#library(vars)
library(maptools)
library(RColorBrewer)
library(lattice)
#library(vars)
#library(lubridate)
library(rEDM)
library(dplyr) 
library(tidyr)
library(metap)
library(kml)


load("ALL_sans_Indo.RData")
main_data <- read.csv("main_data_for_kmeans.csv", row.names = 1)
cal_dis_eucl <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
#nombre de répétition pour le test surrogate
num_surr = 10

####### 

# main_data_2 <- t(main_data)
# main_data_2 <- as.data.frame(main_data_2)
# dengue_dis <- diss(main_data_2,"DWT")
# hist(dengue_dis,100,col= "red",border="black",main=NULL)
# plot(hclust(dengue_dis,"ward.D2"),ylab="heights",main=NULL,xlab=NA,sub=NA,cex=.5)
# 
# dengue_dis2 <- as.matrix(dengue_dis)
# dengue_dis2 <- sweep(dengue_dis2,1,rowMeans(dengue_dis2))
# dengue_dis2 <- sweep(dengue_dis2,2,colMeans(dengue_dis2))
# dengue_dis2 <- (dengue_dis2 + mean(dengue_dis2))/2
# dengue_dis_pca <- prcomp(dengue_dis2,scale=T)
# plot(dengue_dis_pca,main=NULL,col="red")
# 
# dengue_dis_pca <- as.data.frame(predict(dengue_dis_pca))
# with(dengue_dis_pca,plot(PC1,PC2,pch=19,col="red"))
# 
# dengue_dis_km <- kmeans(dengue_dis_pca,3,100,25)
# colors <- adjustcolor(c("red","blue","green"),.5)
# with(dengue_dis_pca,plot(PC1,PC2,pch=19,col=colors[dengue_dis_km$cluster]))
# 
# DKA <-kmeansruns(dengue_dis_pca,krange=1:20,critout=TRUE,runs=2,criterion="asw")
# 
# plot(1:20,DKA$crit,type="b",xlab="number of clusters",
#      ylab="average silhouette width",pch=19,col= "red")
# 
# DKA_2 <-kmeansruns(dengue_dis_pca,krange=1:20,critout=TRUE,runs=2,criterion="ch")
## Calculer le somme cummulaire
main_data_sc <- main_data
for(i in seq(1,267,by = 1))
{
a <- as.numeric(main_data_sc[i,1:153])
main_data_sc[i,1:153] <- cumsum(a)
}

main_data_2_sc <- t(main_data_sc)
main_data_2_sc <- as.data.frame(main_data_2_sc)

## Appliquer le méthode clustering
dengue_dis <- diss(main_data_2_sc,"DWT")
#hist(dengue_dis,100,col= "red",border="black",main=NULL)
#plot(hclust(dengue_dis,"ward.D2"),ylab="heights",main=NULL,xlab=NA,sub=NA,cex=.5)
dengue_dis2 <- as.matrix(dengue_dis)
dengue_dis2 <- sweep(dengue_dis2,1,rowMeans(dengue_dis2))
dengue_dis2 <- sweep(dengue_dis2,2,colMeans(dengue_dis2))
# add matrix mean and divide by 2:
dengue_dis2 <- (dengue_dis2 + mean(dengue_dis2))/2
dengue_dis_pca <- prcomp(dengue_dis2,scale=T)
#plot(dengue_dis_pca,main=NULL,col="red")

dengue_dis_pca <- as.data.frame(predict(dengue_dis_pca))
with(dengue_dis_pca,plot(PC1,PC2,pch=19,col="red"))
dengue_dis_km <- kmeans(dengue_dis_pca,3,100,25)
colors <- adjustcolor(c("red","blue","green"),.5)
with(dengue_dis_pca,plot(PC1,PC2,pch=19,col=colors[dengue_dis_km$cluster]))

DKA_sc <-kmeansruns(dengue_dis_pca,krange=1:20,critout=TRUE,runs=2,criterion="asw")
#plot(1:20,DKA$crit,type="b",xlab="number of clusters",
#     ylab="average silhouette width",pch=19,col= "red")

#DKA_sc2 <-kmeansruns(dengue_dis_pca,krange=1:20,critout=TRUE,runs=2,criterion="ch")
#plot(1:20,DKA_2$crit,type="b",xlab="number of clusters", ylab="average silhouette width",pch=19,col= "red")





all_country_sans_indo_1 <- all_country

cluster_1_sc <- DKA_sc$cluster
k3 = max(cluster_1_sc)

colors_10 <- adjustcolor(c("red","blue","green","yellow","cyan","pink","orange"
                           ,"black","brown","gray", "darkviolet","deeppink","gold","forestgreen","darksalmon","darkolivegreen1"
                           ,"darkmagenta","deepskyblue","blueviolet","chocolate"),.5)

all_country_sans_indo_1@data <- merge(all_country_sans_indo_1@data
                                      , cluster_1_sc, by.x = "VARNAME_1", by.y = "row.names", sort = FALSE, all.x = TRUE)

spplot(all_country_sans_indo_1,215, col.regions= colors_10[1:k3]
       ,at = seq(0.5,0.5+k3,by =1), main = paste("Donnée-sommes cumulaire avec k = "
                                                 ,k3, "Validation Method : Silhouette"))

result_k_means <- NULL
for(i in 1:301)
{
x <- mean(all_country_sans_indo_1@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,1]) 
y <- mean(all_country_sans_indo_1@polygons[[as.numeric(i)]]@Polygons[[1]]@coords[,2])
z <- all_country_sans_indo_1$y[i]
a <- cbind(x,y,z)
result_k_means <- rbind(result_k_means,a)

}
result_k_means <- as.data.frame(result_k_means)
colnames(result_k_means) <- c("x","y","k")
rownames(result_k_means) <- all_country_sans_indo_1@data$VARNAME_1
result_k_means_cut <- result_k_means[1:267,]

result_k_means_cut <- result_k_means_cut[order(result_k_means_cut$k),]


result_k_means_split <- split(result_k_means_cut,result_k_means_cut$k)
mat1 <- result_k_means_split$`1`
mat2 <- result_k_means_split$`2`
mat3 <- result_k_means_split$`3`

mat_dis_1 <- NULL
for ( i in seq(1, nrow(mat1), by = 1))
{
  x1 <- cbind( mat1[i,]$x, mat1[i,]$y)
  temp <- NULL
  for (j in seq(1, nrow(mat1), by = 1))
  {
    x2 <- cbind( mat1[j,]$x, mat1[j,]$y)
    a <- cal_dis_eucl(x1,x2)
    temp <- rbind(temp,a)
  }
  mat_dis_1 <- cbind(mat_dis_1, temp)
}

mat_dis_2 <- NULL
for ( i in seq(1, nrow(mat2), by = 1))
{
  x1 <- cbind( mat2[i,]$x, mat2[i,]$y)
  temp <- NULL
  for (j in seq(1, nrow(mat2), by = 1))
  {
    x2 <- cbind( mat2[j,]$x, mat2[j,]$y)
    a <- cal_dis_eucl(x1,x2)
    temp <- rbind(temp,a)
  }
  mat_dis_2 <- cbind(mat_dis_2, temp)
}

mat_dis_3 <- NULL
for ( i in seq(1, nrow(mat3), by = 1))
{
  x1 <- cbind( mat3[i,]$x, mat3[i,]$y)
  temp <- NULL
  for (j in seq(1, nrow(mat3), by = 1))
  {
    x2 <- cbind( mat3[j,]$x, mat3[j,]$y)
    a <- cal_dis_eucl(x1,x2)
    temp <- rbind(temp,a)
  }
  mat_dis_3 <- cbind(mat_dis_3, temp)
}

mat_vec <- NULL


for(i in seq(1,nrow(mat_dis_1), by = 1))
{
  for (j in seq(i, nrow(mat_dis_1), by = 1))
  {
    mat_vec <- cbind(mat_vec,mat_dis_1[i,j])
  }
}

for(i in seq(1,nrow(mat_dis_2), by = 1))
{
  for (j in seq(i, nrow(mat_dis_2), by = 1))
  {
    mat_vec <- cbind(mat_vec,mat_dis_2[i,j])
  }
}

for(i in seq(1,nrow(mat_dis_3), by = 1))
{
  for (j in seq(i, nrow(mat_dis_3), by = 1))
  {
    mat_vec <- cbind(mat_vec,mat_dis_3[i,j])
  }
}


mat_vec <- mat_vec[which(mat_vec[1,] != 0)]

#le moyenne original
moyen_origine = mean(mat_vec)



#######################
#Calculer avec les test surrogates
######################
result_surrogate <- NULL
#num_surr = 1000
#numprov = 267

for (t in seq(1,num_surr,by = 1))
{
  k_means_sample <- result_k_means_cut
  xx <- sample(1:3, nrow(k_means_sample), replace = T)
  k_means_sample$k <- xx
  
  k_means_sample <- k_means_sample[order(k_means_sample$k),]
  k_means_sample <- split(k_means_sample,k_means_sample$k)
  mat1 <- k_means_sample$`1`
  mat2 <- k_means_sample$`2`
  mat3 <- k_means_sample$`3`
  
  mat_dis_1 <- NULL
  for ( i in seq(1, nrow(mat1), by = 1))
  {
    x1 <- cbind( mat1[i,]$x, mat1[i,]$y)
    temp <- NULL
    for (j in seq(1, nrow(mat1), by = 1))
    {
      x2 <- cbind( mat1[j,]$x, mat1[j,]$y)
      a <- cal_dis_eucl(x1,x2)
      temp <- rbind(temp,a)
    }
    mat_dis_1 <- cbind(mat_dis_1, temp)
  }
  
  mat_dis_2 <- NULL
  for ( i in seq(1, nrow(mat2), by = 1))
  {
    x1 <- cbind( mat2[i,]$x, mat2[i,]$y)
    temp <- NULL
    for (j in seq(1, nrow(mat2), by = 1))
    {
      x2 <- cbind( mat2[j,]$x, mat2[j,]$y)
      a <- cal_dis_eucl(x1,x2)
      temp <- rbind(temp,a)
    }
    mat_dis_2 <- cbind(mat_dis_2, temp)
  }
  
  mat_dis_3 <- NULL
  for ( i in seq(1, nrow(mat3), by = 1))
  {
    x1 <- cbind( mat3[i,]$x, mat3[i,]$y)
    temp <- NULL
    for (j in seq(1, nrow(mat3), by = 1))
    {
      x2 <- cbind( mat3[j,]$x, mat3[j,]$y)
      a <- cal_dis_eucl(x1,x2)
      temp <- rbind(temp,a)
    }
    mat_dis_3 <- cbind(mat_dis_3, temp)
  }
  
  mat_vec <- NULL
  
  
  
  for(i in seq(1,nrow(mat_dis_1), by = 1))
  {
    for (j in seq(i, nrow(mat_dis_1), by = 1))
    {
      mat_vec <- cbind(mat_vec,mat_dis_1[i,j])
    }
  }
  
  for(i in seq(1,nrow(mat_dis_2), by = 1))
  {
    for (j in seq(i, nrow(mat_dis_2), by = 1))
    {
      mat_vec <- cbind(mat_vec,mat_dis_2[i,j])
    }
  }
  
  for(i in seq(1,nrow(mat_dis_3), by = 1))
  {
    for (j in seq(i, nrow(mat_dis_3), by = 1))
    {
      mat_vec <- cbind(mat_vec,mat_dis_3[i,j])
    }
  }
  
  
  mat_vec <- mat_vec[which(mat_vec[1,] != 0)]
  
  m1 = mean(mat_vec)
  
  result_surrogate <- cbind(result_surrogate, m1)
}





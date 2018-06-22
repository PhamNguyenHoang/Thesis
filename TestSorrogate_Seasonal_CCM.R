#Testing seasonal & Surrogate

  
  maindata <- read.csv("DataCCM.txt",header = TRUE)
  data.vn <- maindata %>% filter(country == 'Vietnam') %>% 
    filter(year >= 1996)
  
  data.vn.trans <- data.vn %>% select(-country) %>% spread(variable,value) %>% mutate(date = ISOdate(year,month,day)) %>% select(-year,-month,-day) %>% select(date,flu,everything())
  
  set.seed(599213)
  
  
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
  
  
  yearday_anom <- function(t,x)
  { 
    # t: date formatted with POSIXt 
    # x: time-series values to compute seasonal mean and anomaly
    doy <- as.numeric(strftime(t, format = "%j")) 
    I_use <- which(!is.na(x))
    doy_sm <- rep(doy[I_use],3) + rep(c(-12,0,12),each=length(I_use)) 
    x_sm <- rep(x[I_use],3)
    xsp <- smooth.spline(doy_sm, y = x_sm, w = NULL, spar = 1, cv = NA, all.knots = TRUE,keep.data = TRUE, df.offset = 0)
    xbar <- data.frame(t=t,doy=doy) %>% left_join(data.frame(doy=xsp$x,xbar=xsp$y),by='doy') %>% select(xbar)
  
    out = data.frame(t=t,mean=xbar,anomaly=(x - xbar)) 
    names(out) <- c('date','mean','anomaly') 
    return(out)
  
  }
  
  out <- yearday_anom(data.vn.trans$date,data.vn.trans$AH) 
  AH.bar <- out$mean 
  AH.tilde <- out$anomaly
  AH.surrs <- do.call(cbind, lapply(1:500, function(i) 
    {
    I_na <- is.na(AH.tilde)
    out <- AH.bar
    out[I_na] <- NA
    out[!I_na] <- out[!I_na] + sample(AH.tilde[!I_na], sum(!I_na), replace = FALSE)
    return(out)
    } ))
  
  
  aaaa <- surrogate(climat_1, method = "aaft")
  plot(aaaa, type = "time" )
  lines(1:153,climat_1,  col = "blue")
  
  plot(data.vn.trans$date,data.vn.trans$AH, type = "l", col = "blue",xlab = "date", ylab= "Absolute Humidity")
  lines(data.vn.trans$date,AH.surrs[,11],col = "red")
  
  legend(x = "bottomleft", 
         legend = c("Observed", "Seasonal Surrogate"), 
         col = c("blue", "red"), lwd = 1, lty = c(1,1), bg = "white")  
  
  

  
  
  time.tes <- data.vn.trans$date[1:153]
  ah.tes <- data.vn.trans$AH[1:153]
  
  plot(time.tes, climat_1, type = "l", col = "blue",xlab = "date", ylab= "Absolute Humidity")
  lines(time.tes,AH.surrs[,4],col = "red")
  
  
  
  make_surrogate_data <- function(ts, method = c("random_shuffle", "ebisuzaki", "seasonal"), 
                                  num_surr = 100, T_period = 1)
  {  
    method <- match.arg(method)
    if(method == "random_shuffle")
    {
      return(sapply(1:num_surr, function(i) {
        sample(ts, size = length(ts))
      }))
    }
    else if(method == "ebisuzaki")
    {
      if(any(!is.finite(ts)))
        stop("input time series contained invalid values")
      
      n <- length(ts)
      n2 <- floor(n/2)
      
      mu <- mean(ts)
      sigma <- sd(ts)
      a <- fft(ts)
      amplitudes <- abs(a)
      amplitudes[1] <- 0
      
      return(sapply(1:num_surr, function(i) {
        if(n %% 2 == 0) # even length
        {
          thetas <- 2*pi*runif(n2-1)
          angles <- c(0, thetas, 0, -rev(thetas))
          recf <- amplitudes * exp(complex(imaginary = angles))
          recf[n2] <- complex(real = sqrt(2) * amplitudes[n2] * cos(runif(1)*2*pi))
        }
        else # odd length
        {
          thetas <- 2*pi*runif(n2)
          angles <- c(0, thetas, -rev(thetas))
          recf <- amplitudes * exp(complex(imaginary = angles))
        }
        temp <- Re(fft(recf, inverse = T) / n)
        
        # adjust variance of the surrogate time series to match the original            
        return(temp / sd(temp) * sigma)
      }))
    }
    else
    {
      if(any(!is.finite(ts)))
        stop("input time series contained invalid values")
      
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
  }
  
  
  test_nonlinearity <- function(ts, method = "ebisuzaki", num_surr = 200, T_period = 1, E = 1, ...)
  {
    compute_stats <- function(ts, ...)
    {
      results <- s_map(ts, stats_only = TRUE, silent = TRUE, ...)
      delta_rho <- max(results$rho) - results$rho[results$theta == 0]
      delta_mae <- results$mae[results$theta == 0] - min(results$mae)
      return(c(delta_rho = delta_rho, delta_mae = delta_mae))
    }
    
    actual_stats <- compute_stats(ts, ...)
    delta_rho <- actual_stats["delta_rho"]
    delta_mae <- actual_stats["delta_mae"]
    names(delta_rho) <- NULL
    names(delta_mae) <- NULL
    surrogate_data <- make_surrogate_data(ts, method, num_surr, T_period)
    null_stats <- data.frame(t(apply(surrogate_data, 2, compute_stats, ...)))
    
    return(data.frame(delta_rho = delta_rho, 
                      delta_mae = delta_mae, 
                      num_surr = num_surr, 
                      E = E, 
                      delta_rho_p_value = (sum(null_stats$delta_rho > delta_rho)+1) / (num_surr+1), 
                      delta_mae_p_value = (sum(null_stats$delta_mae > delta_mae)+1) / (num_surr+1)))
  }
  
  
  test11 <- test_nonlinearity(climat_1, method = "seasonal", num_surr = 100, T_period = 12, E = 1)
  
  
  compute_stats <- function(ts, ...)
  {
    results <- s_map(ts, stats_only = TRUE, silent = TRUE, ...)
    delta_rho <- max(results$rho) - results$rho[results$theta == 0]
    delta_mae <- results$mae[results$theta == 0] - min(results$mae)
    return(c(delta_rho = delta_rho, delta_mae = delta_mae))
  }
  
  
  
  ###Test exploration 
  
  df.long <- read.csv("DataCCM.txt",header = TRUE)
  
  df.long <- df.long %>% filter(country == 'Denmark') %>% filter(year >= 1996)
  df.in <- df.long %>% select(-country) %>% spread(variable,value) %>% mutate(date = ISOdate(year,month,day)) %>% select(-year,-month,-day) %>% select(date,flu,everything())
  
  
  df.in.SE <- df.in %>% mutate(flu = flu / (mean(flu,na.rm=TRUE)*365/7))
  
  idx <- is.finite(df.in.SE$flu) & is.finite(df.in.SE$AH) 
  df.in.SE <- df.in.SE[idx,]
  block.SE <- make_block(df.in.SE[,c('flu','AH')],c(rep(1,E.SE),2),c(tp=tp.SE,0:-(E.SE-2),0) ) %>% as.data.frame()
  block.norm.SE <- block.SE 
  for(j in 1:NCOL(block.norm.SE)) 
    block.norm.SE[,j] <- (block.norm.SE[,j] - mean(block.norm.SE[,j], na.rm = TRUE)) / sd(block.norm.SE[,j], na.rm = TRUE) 
  
  block.norm.SE <- data.frame(time = df.in.SE$date, block.norm.SE) 
  
  dAH.norm <- dAH / sd(block.SE$AH_t,na.rm=TRUE)
  
  pred.SE <- c(1,NROW(block.SE)) 
  lib.SE <- pred.SE + NROW(block.SE)
  
  block.temp1 <- block.norm.SE %>% mutate(AH_t = AH_t + dAH.norm/2) %>% bind_rows(block.norm.SE)
  
  out <- block_lnlp(block=block.temp1, 
                    lib = lib.SE, 
                    pred = pred.SE,
                    method = 's-map',
                    tp = 0,
                    num_neighbors = 0,
                    columns = 1+(1:E.SE),
                    target_column = 1,
                    stats_only = FALSE,
                    first_column_time = TRUE,
                    exclusion_radius = 0,
                    theta = theta.SE)[[1]]$model_output[1:NROW(block.SE),]
  flu.pred.pos <- out$pred
  
  
  
  pred.SE <- c(1,NROW(block.SE)) 
  
  lib.SE <- pred.SE + NROW(block.SE)
  
  block.temp2 <- block.norm.SE %>% mutate(AH_t = AH_t - dAH.norm/2) %>% bind_rows(block.norm.SE)
  
  out <- block_lnlp(block=block.temp2,
                    lib = lib.SE,
                    pred = pred.SE,
                    method = 's-map',
                    tp = 0,
                    num_neighbors = 0,
                    columns = 1+(1:E.SE),
                    target_column = 1,
                    stats_only = FALSE, 
                    first_column_time = TRUE,
                    exclusion_radius = 0,
                    theta = theta.SE)[[1]]$model_output[1:NROW(block.SE),]
  flu.pred.neg <- out$pred
  dFlu.pred <- flu.pred.pos - flu.pred.neg
  dFlu.pred <- dFlu.pred * sd(block.SE$d_t+2,na.rm=TRUE)
  
  
  
  
  
  
  
  

library(dplyr)
library(sars)
library(ggplot2)
library(gridExtra)
library(foreach)
library(doParallel)
library(cluster)
library(mgcv)


##Load Data

wwd <- getwd()
ldf <- readRDS("C:\\Users\\Tom\\Desktop\\ISAR_Extrap\\Data\\Datasets.rds")
setwd(wwd)

if (any(sapply(ldf, nrow) < 8)) stop("ewsrdgfserg")

#get predictors for the subsetted dataset (i.e. the one which the predictions are based on)

get_preds <- function(d, th = 0.5){
  
  #subset the data based on th
  d2 <- data.frame("a" = d$a, "s" = d$s)
  d3 <- d2[order(d2$a),] # orders dataset based on area
  y <- d3$a[nrow(d2)]
  maxSp <- d3$s[nrow(d3)]
  yDeng <- th * y
  x <- filter(d3, a < yDeng) 
  
  #get preds
  Amax <- max(x$a)
  Amin <- min(x$a)
  Smax <- max(x$s) + 1
  Smin <- min(x$s) + 1
  Ascale <- Amax / Amin
  Sscale <- Smax / Smin
  Noisl <- nrow(x)
  
  res <- c(Amin, Amax, Ascale, Smin, Smax, Sscale, Noisl)
  return(res)
}

#checks normality, homogeneity, neg values and P value of slope
pow_check <- function(d, th = 0.5){
  
  d2 <- data.frame("a" = d$a, "s" = d$s)
  d3 <- d2[order(d2$a),] # orders dataset based on area
  y <- d3$a[nrow(d2)]
  maxSp <- d3$s[nrow(d3)]
  yDeng <- th * y
  nD <- filter(d3, a < yDeng) 
 
  s <- sar_power(nD, normaTest = "shapiro", 
                 homoTest = "cor.fitted", homoCor = "pearson",
                 grid_start = "partial")
  nt <- s$normaTest[[2]]$p.value
  ht <- s$homoTest[[2]]$p.value
  nv <- ifelse(s$neg_check, 1, 0) #1 for negative; 0 for all ok
  p <- as.vector(summary(s)$Parameters[2,4])
  res <- round(c(nt, ht, nv, p), 6)
  return(res)
}

###main extrapolation function
extrap <- function(d, th = 0.5, CI = FALSE, n = 100, NT = "shapiro", All = FALSE,
                   grid_start = "partial", ...){
  if (CI && All) stop ("All and CI cannot both be TRUE")
  d2 <- data.frame("a" = d$a, "s" = d$s)
  d3 <- d2[order(d2$a),] # orders dataset based on area
  y <- d3$a[nrow(d2)]
  maxSp <- d3$s[nrow(d3)]
  yDeng <- th * y
  nD <- filter(d3, a < yDeng) 
  
  #check new dataset has > 5 data points
  xx <- nrow(nD)
  if (xx < 6) stop("less than 6 islands")

  #power model
  s <- sar_power(nD, grid_start = grid_start)
  predPow <- sar_pred(s, y)$Prediction
  
  #CONFIDENCE INTERVALS:power
  if (CI){
    predCI <- vector(length = n)
    k <- 1
    while (k < (n + 1)){
      dum <- sample(seq(1, nrow(nD), 1), nrow(nD), replace = TRUE)
      dDF <- nD[dum,]
      dum2 <- tryCatch(sar_power(dDF), error = function(e) NA)
      if (length(dum2) == 1) next
      predCI[k] <- sar_pred(dum2, y)$Prediction
      k <- k + 1
    }
    CIs_pow <-  quantile(predCI, c(0.025, 0.975))
  }#eo CI

  ##multi SAR predictions
  
  #if wanting to run without linear model, add obj = x3 into
  #the sar_average function
  # x3 <- c("power",
  #         "powerR","epm1","epm2","p1","p2","loga","koba",
  #         "monod","negexpo","chapman","weibull3","asymp",
  #         "ratio","gompertz","weibull4","betap","logistic", "heleg")
  
 sm <- tryCatch(suppressMessages(sar_average(data = nD, crit = "AICc", 
                                              normaTest = NT, verb = FALSE, 
                                              homoTest = "cor.fitted", 
                                              homoCor = "pearson",
                                              grid_start = grid_start)), error = function(e) NA)
  

  #return all NAs if sm fails
  if (length(sm) == 1){
    if (CI) return(rep(NA, 14))
    if (All) return(rep(NA, 70))
    return(rep(NA, 10))
  }
  predMulti <- sar_pred(sm, y)$Prediction
  
  #best model (& get predictions of best model)
  bm <- which(sm$details$delta_ics == 0) %>% names()
  bmi <- sm$details$fits[bm]
  bmp <- sar_pred(bmi[[1]], y)$Prediction
  

  #confidence intervals:multi
  if (CI){
  x3 <- c("power",
            "powerR","epm1","epm2","p1","p2","loga","koba",
            "monod","negexpo","chapman","weibull3","asymp",
            "ratio","gompertz","weibull4","betap","logistic", "heleg", "linear")
    
    allN <- sm$details$mod_names %>% names()
    whN <- which(x3 %in% allN)
    x3S <- x3[whN]
    if (length(x3S) == 0) stop("error with model names in CI code: A")
    if (length(x3S) != length(sm$details$mod_names)) stop("error with model names in CI code: B")
    #x3S just equals allN so do a check for sense
    if (!all(allN ==x3S)) stop("error with model names in CI code: C")
    
    predCI <- vector(length = n)
    k <- 1
    k2 <- 1
    while (k < (n + 1)){
      k2 <- k2 + 1
      #in case k never reaches n
      if (k2 > 250) {
        if (k >=25){
          ppn <- paste0(sum(d2$a), ".csv")
          write.csv(nD, file = ppn)
          break
        } else{
        return(rep(NA, 14))
        }
      }
      dum <- sample(seq(1, nrow(nD), 1), nrow(nD), replace = TRUE)
      dDF <- nD[dum,]
      #it is fitting the models that were fitted in sm (x3S) so that is why normatTest and cor test were set to 
      #none, to ensure these models are all fitted. But does still exclude neg fitted values
      # dum2 <- tryCatch(suppressMessages(sar_average(data = dDF, obj = x3S, crit = "AICc", normaTest = "none",
      #         homoTest = "none",
      #      verb = FALSE, ...)), error = function(e) NA)
      dum2 <- tryCatch(suppressMessages(sar_average(data = dDF, obj = x3S, crit = "AICc", normaTest = "none",
                                                    homoTest = "none",
                                                    verb = FALSE, neg_check = TRUE)), error = function(e) NA)
      
      
      if (length(dum2) == 1 || length(dum2$details$mod_names) != length(x3S)) next
      pd <-  sar_pred(dum2, y)$Prediction
      if (is.na(pd) || is.infinite(pd) || pd < 0) next
      predCI[k] <- pd
      k <- k + 1
    }
    CIs_mul <- quantile(predCI, c(0.025, 0.975))
  }#eo CI
  
  #weight of power (if in it)
  if ("Power" %in% sm$details$mod_names){
    weiP <- sm$details$weights_ics["power"]
  } else {
    weiP <- NA
  }
  
    #LEE metric of Dengler
  leeP <- log10(predPow) - log10(maxSp) 
  leeM <- log10(predMulti) - log10(maxSp) 
  leeB <- log10(bmp) - log10(maxSp)

  
  ##if All == TRUE, get extrap predictions for all twenty ISAR models
  if (All){
    x3 <- c("power",
      "powerR","epm1","epm2","p1","p2","loga","koba",
      "monod","negexpo","chapman","weibull3","asymp",
      "ratio","gompertz","weibull4","betap","logistic", "heleg", "linear")
    
     #make vector of all sar model functions
    funcs <- sapply(x3, function(x) paste0("sar_", x))
    #run all sar model functions with dataset nD
    allFits0 <- suppressWarnings(lapply(funcs, function(x){
      
      if(x != "sar_linear"){
        do.call(x, args = list("data" = nD))
      } else{
        do.call(x, args = list("data" = nD))
      }
      
    }
      
      ))

    #work out any bad models and remove;
    #this removes bad models based on the multi-model SAR curve construction
    #(i.e. any models that fail validation checks)
    badMods <- x3[which(!x3 %in% names(sm$details$weights_ics))]
    if (length(badMods) >= 1){
      if (!all(badMods %in% names(allFits0))) stop("aaaaaaaa")
      wb <- which(names(allFits0) %in% badMods)
      allFits0[wb] <- NA
      allFits <- allFits0
    } else{
      allFits <- allFits0
    }
    
    #get the predictions for each 
    allPred <- lapply(allFits, function(x) {
      if (length(x) == 1) return(NA)
      xx <- sar_pred(x, y)$Prediction
  #    if (xx == 0) return(NA)#if singular gradient at parameter estimates, then returns 0
      xx
  })
    allVec <- unlist(allPred)
    #get model AICc weight values
    del <- vector(length = length(x3))
    for (i in seq_along(x3)){
      if (x3[i] %in% names(sm$details$weights_ics)){
        w <- which(names(sm$details$weights_ics) == x3[i])
        del[i] <- sm$details$weights_ics[w]
      } else {
        del[i] <- NA
      }
    }
    #do a check to make sure the model with max weight is the best model from above
    w2 <- which.max(del)
    if (bm != x3[w2]) stop("model name order wrong")
    if (weiP != del[1] && !is.na(weiP)) stop("Power weights do not match")
    
    #get LEE metric for each model
    leeA <- vector(length = length(x3))
    for (i in seq_along(x3)){
      if (is.na(del[i])){
        leeA[i] <- NA
      } else {
        leeA[i] <- log10(allVec[i]) - log10(maxSp)
      }
    }
    #some models can give negative Extrapolated prediction values. If this is the case, the lee
    #value is NaN (not NA). If this happens skip the error check as will error if any NaNs
    if (!any(is.nan(leeA))) if (length(sm$details$ics) != length(na.omit(leeA))) stop("model numbers wrong")
    res <- c(c(maxSp, predPow, predMulti, weiP, xx, leeP, leeM, bm, bmp, leeB), allVec, del, leeA)
    return(res)
  }#eo if All
  
  #results
  if (CI){
    res <- c(maxSp, predPow, predMulti, weiP, xx, leeP, leeM, bm, bmp, leeB, CIs_pow, CIs_mul)
  } else {
    res <- c(maxSp, predPow, predMulti, weiP, xx, leeP, leeM, bm, bmp, leeB)
  }
  return(res)
}
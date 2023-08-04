###################################################
#######RUN CODE###################
###################################################

source("R code\\source_code_V3.R")

############################################################
############Run main extrapolation function#################
############################################################

#####For when wanting to re-run with out linear model, add 'obj = x6' into the extra function below
#x6 = c("power", "powerR","epm1","epm2","p1","p2","loga","koba","mmf","monod","negexpo","chapman",
 #       "weibull3","asymp","ratio","gompertz","weibull4","betap","heleg")

cores = 10
cl = makeCluster(cores); on.exit(stopCluster(cl))
registerDoParallel(cl)
i = 1 #Dummy line for RStudio warnings


##run across all dataset
resM = foreach(i=seq(from=1, to=length(ldf), by=1))  %dopar% {
  library(dplyr)
  library(sars)
 extrap(d = ldf[[i]], th = 0.5, NT = "shapiro", 
                      neg_check = TRUE, grid_start = "partial")
}
  
  
resM <- matrix(unlist(resM), nrow = length(ldf), ncol = 10, byrow = T)  
  
colnames(resM) <- c("Obs", "Pow", "Multi", "Power_weight", "rows", "Lee_Pow", "Lee_Mult","Best_mod",
                    "Best", "LEE_Best")
resM <- as.data.frame(resM)

#save(resM, file = "resM_noLin.R")
#load(file.choose())

if (any(is.na(resM[,8]))){
  whna <- which(is.na(resM[, 8]))
  resM <- resM[-whna, ]
  cat(length(whna))
}

#check if any predictions are the same
cc <- c()
for (i in 1:nrow(resM)){
  
  if (resM[i,]$Pow == resM[i,]$Multi){
    cc[i] <- i
  } else{
    cc[i] <- 0
  }} 
if (any(cc != 0)) stop("frety")
  
#need to convert certain cols to numeric (as are factors)
resM[ ,c(1:7, 9:10)] <- apply(resM[ ,c(1:7, 9:10)], 2, function(x) as.numeric(as.character(x)))

#check which prediction is close
bp <- vector(length = nrow(resM))#just power vs multi
bp2 <- vector(length = nrow(resM))#power, multi, best

for (i in 1:nrow(resM)){
  e <- which.min(abs(unlist(resM[i,2:3]) - unlist(resM[i,1])))
  bp[i] <- ifelse(e == 1, "1", "0") #1 = power, 0 = multi
  e2 <- (abs(unlist(resM[i,c(2:3, 9)]) - unlist(resM[i,1])))
  dum <- which.min(e2)
  if (dum == 1){
    if (e2[1] == e2[3]){ #because if they are identical which.min picks e2[1] which is the power
      bp2[i] <- 3
    } else {
      bp2[i] <- 1
    }
  } else {
    bp2[i] <- dum #1 = power, 2 multi, 3 best
  }
}#eo for

resM$Best_pred <- as.numeric(as.character(bp))#1 = power
resM$Best_pred2 <- as.numeric(as.character(bp2))

table(resM$Best_pred)
table(resM$Best_pred2)

table(resM$Best_mod)

#mean weight of power model
NN <- filter(resM, (!is.na(Power_weight)))
mean(NN$Power_weight)

write.csv(resM, file = "resM_0.5.csv")

##summary stats
median(resM$Lee_Pow)
quantile(resM$Lee_Pow, probs = c(0.025, 0.975))
median(abs(resM$Lee_Pow))
quantile(abs(resM$Lee_Pow), probs = c(0.025, 0.975))

median(resM$Lee_Mult)
quantile(resM$Lee_Mult, probs = c(0.025, 0.975))
median(abs(resM$Lee_Mult))
quantile(abs(resM$Lee_Mult), probs = c(0.025, 0.975))

#over or under
length(which(resM$Lee_Pow > 0)) / 119
length(which(resM$Lee_Mult > 0)) / 119

###################################################
####RUN AGAIN BUT WITH CONFIDENCE INTERVALS########
#####################################################
##main parallel for loop
resM_CI = foreach(i=seq(from=1, to=length(ldf), by=1))  %dopar% {
  library(dplyr)
  library(sars)
 caro <- extrap(d = ldf[[i]], th = 0.5, NT = "shapiro",
         CI = TRUE, n = 100, neg_check = TRUE,
         grid_start = "partial")
  pp = paste0(i, ".csv")
  write.csv(4, file = pp)
  caro
}

#save(resM_CI, file = "resM_CI.R")
#load(file.choose())

resMCI <-  t(matrix(nrow = 14, unlist(resM_CI))) #parallel returns a list; so turn into a df
resMCI <- as.data.frame(resMCI)

colnames(resMCI) <- c("Obs", "Pow", "Multi", "Power_weight", "rows", "Lee_Pow", "Lee_Mult","Best_mod",
                    "Best", "LEE_Best", "CIP[1]", "CIP[2]", "CIM[1]", "CIM[2]")

if (any(is.na(resMCI[,8]))){
  whna2 <- which(is.na(resMCI[, 8]))
  resMCI <- resMCI[-whna2, ]
  cat(length(whna2))
}

#need to convert certain cols to numeric (as are factors)
resMCI[ ,c(1:7, 9:14)] <- apply(resMCI[ ,c(1:7, 9:14)], 2, 
                                function(x) as.numeric(as.character(x)))
resMCI[ ,11:14] <- apply(resMCI[ ,11:14], 2, round, digits = 2)

write.csv(resMCI, file = "resM_0.5_CIsC.csv")

#calculate median CI width
pCr <- resMCI$'CIP[2]' - resMCI$'CIP[1]'
median(pCr)
mCr <- resMCI$'CIM[2]' - resMCI$'CIM[1]'
median(mCr)

############################################
####RUN AGAIN BUT WITH ALL TWENTY MODELS
########################################

#-Inf cases in model LEE values are when there is the singular gradient warning
#as it returns 0 species predicted and so can't be log-transformed

#NaNs in the Lee values mean a model had a negative extrapolated predicted
#value; when written to csv it is converted to NA

#remember: best model column is best model to subset of data, not best
#extrapolation prediction

resA = foreach(i=seq(from=1, to=length(ldf), by=1))  %dopar% {
  library(dplyr)
  library(sars)
extrap(d = ldf[[i]], th = 0.5, 
                      NT = "shapiro", All = TRUE, neg_check = TRUE)
}


resA <- matrix(unlist(resA), nrow = length(ldf), 
               ncol = 70, byrow = T)  

colnames(resA) <- c("Obs", "Pow", "Multi", "Power_weight", "rows", "Lee_Pow", 
                    "Lee_Mult","Best_mod",
                    "Best", "LEE_Best", "power","powerR",   "epm1", 
                    "epm2", "p1", "p2","loga" , "koba" ,   
                    "monod" ,   "negexpo",  "chapman" , 
                    "weibull3" ,"asymp"  ,  "ratio" ,   "gompertz",
                    "weibull4", "betap", "logistic",   "heleg",    "linear", 
                    "power_weight",  "powerR_weight",   "epm1_weight", 
                    "epm2_weight", "p1_weight",       "p2_weight",
                    "loga_weight" ,    "koba_weight" ,   
                    "monod_weight" ,   "negexpo_weight",  
                    "chapman_weight" , "weibull3_weight" ,
                    "asymp_weight"  ,  "ratio_weight" ,   "gompertz_weight",
                    "weibull4_weight", "betap_weight", "logistic_weight",
                    "heleg_weight", "linear_weight",
                    "power_LEE",  "powerR_LEE",   "epm1_LEE", 
                    "epm2_LEE", "p1_LEE",   "p2_LEE","loga_LEE" , 
                    "koba_LEE" ,  "monod_LEE" ,   "negexpo_LEE",  
                    "chapman_LEE" , "weibull3_LEE" ,
                    "asymp_LEE"  ,  "ratio_LEE" ,   "gompertz_LEE",
                    "weibull4_LEE", "betap_LEE",  "logistic_LEE",
                    "heleg_LEE",    "linear_LEE")

resA <- as.data.frame(resA)

#save(resA, file = "resA.R")
#load(file.choose())

anyNA(resA[,8])#do a check

if (any(is.na(resA[,8]))){
  whna <- which(is.na(resA[, 8]))
  resA <- resA[-whna, ]
  cat(length(whna))
}

#need to convert certain cols to numeric (as are factors)
resA[ ,c(1:7, 9:70)] <- apply(resA[ ,c(1:7, 9:70)], 2, 
                              function(x) as.numeric(as.character(x)))

#check which prediction is closest
bp <- vector(length = nrow(resA))#all models vs multi
bp2 <- vector(length = nrow(resA))#all individual models (not multi)

for (i in 1:nrow(resA)){
  e <- which.min(abs(unlist(resA[i,c(3,11:30)]) - unlist(resA[i,1])))
  if (length(e) > 1) stop(paste("Multiple models giving same predictions: A", i))
  bp[i] <- names(e)
  e2 <- which.min(abs(unlist(resA[i,11:30]) - unlist(resA[i,1])))
  if (names(e2) == "Multi") stop(paste("multi should not be in this model set", i))
  if (length(e2) > 1) stop(paste("Multiple models giving same predictions: B", i))
  bp2[i] <- names(e2)
}#eo for

resA$Best_pred <- bp
resA$Best_pred2 <- bp2


table(resA$Best_pred)
table(resA$Best_pred2)

write.csv(resA, file = "resA2C.csv")

##regressions for each model: LEE ~ Weight
#function to return regression model and vector of coefficient, P value and R2

#NB THE FUNCTION DOES NOT WORK WITH THE CSV VERSION OF RESA AS THIS CONVERTS THE
#NAN TO NA SO YOU NEED TO RUN THE RESA CODE ABOVE TO CREATE IT AND THEN RUN
#MODREG ETC

modReg <- function(d1, d2, type = "gam"){#d2 = LEE, d1 weight, type = either lm or gam
  if (any(is.nan(d2))){ #if any Nan in the LEE values (see above) turn to NA 
    wn <- which(is.nan(d2))
    d2[wn] <- NA
    d1[wn] <- NA    
  }
  wn2 <- which(is.na(d1))
  wn3 <- which(is.na(d2))
  if (!identical(wn3,wn2)) stop("Exit")
  d11 <- na.omit(d1)
  d22 <- na.omit(d2)
  if (any(d22 == "-Inf")){ #if any -Inf in the LEE values (see above) turn to NA and remove
    wi <- which(d22 == "-Inf")
    d22[wi] <- NA
    d11[wi] <- NA
    wn4 <- which(is.na(d11))
    wn5 <- which(is.na(d22))
    if (!identical(wn4,wn5)) stop("Exit")
    d11 <- na.omit(d11)
    d22 <- na.omit(d22)
  }

  if (length(d11) != length(d22)) stop("ds dont match")
  if (type == "lm"){
    r <- lm(abs(d22) ~ d11)
    r1 <- r$coefficients[2]
    s <- summary(r)
    r2 <- s$coefficients[8]
    r3 <-  s$r.squared
    res <- list(r, c("Coef" = r1,"P" = r2,"R2" = r3))
  } else if (type == "gam"){
    r <- gam(abs(d22) ~ s(d11))
    s <- summary(r)
    r1 <- s$edf
    r2 <- s$s.table[4]
    r3 <- s$r.sq
    res <- list(r, c("edf" = r1,"P" = r2,"R2" = r3))
  }
  return(res)
}

#create new version of dataset with just columns we need
resA2 <- resA[ ,31:70]

#loop modReg across all models and output a table
out <- matrix(nrow = 20, ncol = 3)
k <- 21 #start for the LEE values

for (i in 1:20){
  m <- modReg(d1 = resA2[ ,i], d2 = resA2[ ,k], 
              type = "gam")
  out[i, ] <- m[[2]]
  k <- k + 1
}

colnames(out) <- c("edf", "P", "R2")
out <- as.data.frame(out)
out2 <- apply(out, 2, round, digits = 4)
nn <- gsub("_weight", "", colnames(resA2)[1:20])
rownames(out2) <- nn
write.csv(round(out2, 2), file ="out2_gamC.csv")

#P value corrected: 0.05 / 20 = 0.0025

#Note, you get very similar results fitting the GAMS
#using REML; only difference is that the logistic model
#is no longer a significant fit. also very similar if you
#use raw LEE values as the response rather than abs()

#get mean AICc weight for each model
meanA <- apply(resA2[, 1:20], 2, function(x) mean(x,na.rm  = TRUE))
wcs(round(meanA, 2))

#############################################################
##GAM analyses
############################################################

#filename orders in ldf and pm2 all match up

#get predictors
predMat <- matrix(nrow = length(ldf), ncol = 7)
for (i in 1:length(ldf)){
  predMat[i, ] <- get_preds(ldf[[i]])
}
predMat <- as.data.frame(predMat)
colnames(predMat) <- c("Amin", "Amax", "Ascale", "Smin", "Smax", "Sscale", "Noisl")

#join with lat, taxa data
pm2 <- read.csv("Data\\preds.csv")[,2:3]
colnames(pm2) <- c("Lat", "Taxon")
predMat <- cbind(predMat, pm2)
#remove from predMat the dataset(s) that were removed from resM
predMat <- predMat[-whna,]

#all preds need log-transforming
predMat[ ,1:7] <- apply(predMat[ ,1:7], 2, log)

regMat <- cbind(resM, predMat)

#Fit LM, check vifs.
#Use gam to fit, so that same functions used for both
tes <- gam(Lee_Pow ~ Amin + Ascale + Smin  + Sscale +
            Noisl + Lat + Taxon, data = regMat)
car::vif(tes)
summary(tes)
#cant have Amax/Smax and the scale equivalents, so have removed the max

#check distribution of preds
pm2 <- select(predMat, c(Amin, Ascale, Smin,
                         Sscale, Noisl, Lat))
apply(pm2, 2, hist)

##Using GAMs
#Note - we used the default setting, which fits using
#the default GCV.Cp method. However, many people recommend
#the use of REML instead. Fitting using REML instead 
#results in a very similar full model.
g <- gam(Lee_Pow ~ s(Amin) + s(Ascale) + s(Smin) +
           s(Sscale) + s(Noisl) + s(Lat) + 
           Taxon, data = regMat)

# gAlt <- gam(Lee_Pow ~ s(Amin) + s(Ascale) + s(Smin) +
#            s(Sscale) + s(Noisl) + s(Lat) + 
#            Taxon, data = regMat, method = "REML")

##Compare the lm with the GAM (using default GAM fitting method),
#using AIC, following Zuur (2013 - beginner's guide to 
#GAMS with R), p.137
AIC(g); AIC(tes)

#Model selection using dredge, and comparing alternative
#GAM model structures (all with default fitting method) using 
#AICc following Zuur (2013), p. 109
options(na.action = "na.fail")
dd <- MuMIn::dredge(g)
filter(dd, delta <=2)
MuMIn::sw(dd)

#Taking the four variables present in the models with
#delta AICc < 2
g2 <- gam(Lee_Pow ~ s(Ascale) + s(Lat) +
            s(Smin) + s(Sscale), data = regMat)
summary(g2)

#Again, very similar when fitted using REML
# g2Alt <- gam(Lee_Pow ~ s(Ascale) + s(Lat) +
#             s(Smin)  + s(Sscale), data = regMat,
#           method = "REML")
# summary(g2Alt)

#check model residuals
r <- g2$residuals
hist(r)
plot(r, fitted(g2))
gam.check(g2)

#saveRDS(list(g,g2), file = "gam_results.rds")

##Make plot for paper
#the select = 1,2,3,4 matches up with the xlab and titles
#get intercept to shift plot
int <- g2$coefficients[1]

jpeg(file = "GAM_fig.jpeg", width = 25, height = 25, res = 600, units = "cm")
par(mfrow = c(2, 2))
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(g2, select = 1, main = expression(italic('A')['scale']), cex.main = 1.8, 
     cex.axis = 1.3, cex.lab = 1.5, xlab = expression(italic('A')['scale']), shift = int,
     ylab = "LEE-POW (partial effect)")
plot(g2, select = 2, main = expression("Latitude"), cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.5, shift = int,
     ylab = "LEE-POW (partial effect)")
plot(g2, select = 3, main = expression(italic('S')['min']), cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.5,
     xlab = expression(italic('S')['min']), shift = int, ylab = "LEE-POW (partial effect)")
plot(g2, select = 4, main = expression(italic('S')['scale']), cex.main = 1.8, cex.axis = 1.3, cex.lab = 1.5,
     xlab = expression(italic('S')['scale']), shift = int, ylab = "LEE-POW (partial effect)")
dev.off()

#Note - in certain places, e.g.,
#https://stats.stackexchange.com/questions/471267/plotting-gams-on-response-scale-with-multiple-smooth-and-linear-terms
#it appears that a better approach to plotting 
#the individual GAM smoother terms than shifting
#by the intercept, when there are multiple smoother terms in
#the model, is to use the model to predict values and plot those.
#In this case, it makes near identical plots.
# new_data1 <- with(regMat, expand.grid(Ascale = seq(min(Ascale),
#                                                  max(Ascale),
#                                              length = 200),
#                                      Lat = median(Lat),
#                                      Smin = median(Smin),
#                                      Sscale = median(Sscale)))
# pred1 <- predict(g2, new_data1, se.fit = TRUE)
# pred1 <- cbind(pred1, new_data1)
# pred1$LI <- pred1$fit - (2*pred1$se.fit)
# pred1$UI <- pred1$fit + (2*pred1$se.fit)
# 
# new_data2 <- with(regMat, expand.grid(Lat = seq(min(Lat),
#                                                    max(Lat),
#                                                    length = 200),
#                                       Ascale = median(Ascale),
#                                       Smin = median(Smin),
#                                       Sscale = median(Sscale)))
# pred2 <- predict(g2, new_data2, se.fit = TRUE)
# pred2 <- cbind(pred2, new_data2)
# pred2$LI <- pred2$fit - (2*pred2$se.fit)
# pred2$UI <- pred2$fit + (2*pred2$se.fit)
# 
# new_data3 <- with(regMat, expand.grid(Smin = seq(min(Smin),
#                                                 max(Smin),
#                                                 length = 200),
#                                       Ascale = median(Ascale),
#                                       Lat = median(Lat),
#                                       Sscale = median(Sscale)))
# pred3 <- predict(g2, new_data3, se.fit = TRUE)
# pred3 <- cbind(pred3, new_data3)
# pred3$LI <- pred3$fit - (2*pred3$se.fit)
# pred3$UI <- pred3$fit + (2*pred3$se.fit)
# 
# new_data4 <- with(regMat, expand.grid(Sscale = seq(min(Sscale),
#                                                 max(Sscale),
#                                                 length = 200),
#                                       Ascale = median(Ascale),
#                                       Smin = median(Smin),
#                                       Lat = median(Lat)))
# pred4 <- predict(g2, new_data4, se.fit = TRUE)
# pred4 <- cbind(pred4, new_data4)
# pred4$LI <- pred4$fit - (2*pred4$se.fit)
# pred4$UI <- pred4$fit + (2*pred4$se.fit)
# 
# jpeg(file = "GAM_fig_Alt.jpeg", width = 25, height = 25, res = 600, units = "cm")
# par(mfrow = c(2,2))
# plot(pred1$Ascale, pred1$fit, ylab = "LEE-POW (partial effect)",
#      type = "l", ylim = c(-0.2, 0.3), yaxt = "n")
# axis(2, at = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3),
#      labels = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3))
# lines(pred1$Ascale, pred1$LI, lty = 3)
# lines(pred1$Ascale, pred1$UI, lty = 3)
# plot(pred2$Lat, pred2$fit, ylab = "LEE-POW (partial effect)",
#      type = "l", ylim = c(-0.2, 0.3), yaxt = "n")
# axis(2, at = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3),
#      labels = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3))
# lines(pred2$Lat, pred2$LI, lty = 3)
# lines(pred2$Lat, pred2$UI, lty = 3)
# plot(pred3$Smin, pred3$fit, ylab = "LEE-POW (partial effect)",
#      type = "l", ylim = c(-0.2, 0.3), yaxt = "n")
# axis(2, at = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3),
#      labels = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3))
# lines(pred3$Smin, pred3$LI, lty = 3)
# lines(pred3$Smin, pred3$UI, lty = 3)
# plot(pred4$Sscale, pred4$fit, ylab = "LEE-POW (partial effect)",
#      type = "l", ylim = c(-0.2, 0.3), yaxt = "n")
# axis(2, at = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3),
#      labels = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3))
# lines(pred4$Sscale, pred4$LI, lty = 3)
# lines(pred4$Sscale, pred4$UI, lty = 3)
# dev.off()
##########################################################
#check which datasets the power model passes checks##########
##############################################################

PC <- t(vapply(ldf, pow_check, FUN.VALUE = numeric(4)))
PC <- as.data.frame(PC)
colnames(PC) <- c("Norm", "Homo", "Neg", "P")
GB <- apply(PC, 1, function(x){
  ifelse((x[1] >= 0.05 && x[2] >= 0.05 && x[3] == 0 && x[4] < 0.05), 1, 0)#1 = good; 0 = bad
})
PC$GB <- GB

table(PC$GB)

write.csv(PC, "power_checks.csv")

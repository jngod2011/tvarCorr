## Download PISY tree ring width data and analyse

## set working directory to cloned GitHub folder
#setwd(paste(getwd()))

## clear global environment
rm(list=ls())

## load libraries
library(dplR)
library(dlm)
library(car) # for delta method
library(ggplot2)
library(zoo)

################################################
## Download PISY tree ring width data from 
################################################
## API URL for PISY data: https://www.ncdc.noaa.gov/paleo-search/study/search.json?headersOnly=true&dataPublisher=NOAA&dataTypeId=18&locations=Continent>Europe>Northern Europe>British Isles>Ireland|Continent>Europe>Northern Europe>British Isles>United Kingdom&keywords=earth science>paleoclimate>tree-ring>width&species=PISY

## download PISY tree ring width data, merge to single dataframe and save to data folder
## notes: Data downloaded 5th December; uncomment to download again

# location.names <- c("brit015"="Inverey", "brit016"="Shieldaig", "brit017"="Ballochbuite", "brit018"="Loch Maree",
#                     "brit019"="Franchise Wood", "brit020"="Glen Derry", "brit021"="Drimmie Schottland", "brit024"="Glen Affric", "brit026"="Coulin", "brit048w"="Plockton", "brit049w"="Mallaig")
# dat <- list()
# #file.names <- list.files("PISY_1089559831_2017-12-05/width_ring_raw_data")
# file.names <- c("brit015.rwl", "brit018.rwl", "brit019.rwl", "brit026.rwl", "brit020.rwl")

# for(i in 1:length(file.names)){
#      name <- strsplit(file.names[i], split=".rwl")[[1]]
#      dat[[name]] <- read.rwl(paste("PISY_1089559831_2017-12-05/width_ring_raw_data", file.names[i], sep="/"))
#      colnames(dat[[name]]) <- paste(paste(substring(strsplit(location.names[name], split=" ")[[1]], 1 ,1), collapse=""), colnames(dat[[name]]), sep="")
#      ## sample 15 random trees from each location 
#      set.seed(1)
#      dat[[name]] <- dat[[name]][,sample(colnames(dat[[name]]),15)]
#      dat[[name]]$year <- rownames(dat[[name]])
# }
# ## merge all dataframes in the list to a single dataframe and save
# dat.tree.rings <- Reduce(function(...) merge(..., all=TRUE), dat)
#write.csv2(dat.tree.rings, file="/data/tree_ring_data/PISY/PISY_tree_rings.csv")

#######################################################
## load PISY tree ring data and prepare for time-varying correlation analysis
#######################################################
## load downloaded PISY tree ring data into R
dat.tree.rings <- read.csv2(file=paste(getwd(), "/data/tree_ring_data/PISY/PISY_tree_rings.csv", sep="")) ## loads data from cloned gitHub repository
rownames(dat.tree.rings) <- dat.tree.rings$year
dat.tree.rings <- dat.tree.rings[,-1] # remove X column
ring.years <- dat.tree.rings$year
dat.tree.rings <- dat.tree.rings[, colnames(dat.tree.rings) != "year"]
tree.IDs <- colnames(dat.tree.rings)
tree.pairs <- unique(t(apply(expand.grid(tree.IDs, tree.IDs), 1, sort)))
length(tree.pairs)

## plot all tree ring time series
# dev.new()
# pdf("PISY_tree_ring_plots.pdf", height = 10, width = 10*1.618)
# par(mfrow=c(8,8), oma=c(0,0,0,0), mar=c(4,4,1,1))
# for(i in 1:length(tree.IDs)){
#      id <- tree.IDs[i]
#      plot(x=as.numeric(rownames(dat.tree.rings)), y=dat.tree.rings[,id], type="l", main=id, ylab="inc width measurement", xlab="Year", ylim=c(0,10))
# }
# dev.off()

## create list of tree pairwise comparison dataframes
## note: removed rows where NAs are in both columns
tree.pair.ts.list <- list()
for (i in 1:nrow(tree.pairs)){
  pair <- tree.pairs[i,]
  pair.name <- paste(pair[1], pair[2], sep=" + ")
  print(c(i, pair.name)) # to trace progress of the for loop
  tmp <- cbind(dat.tree.rings[,pair[1]], dat.tree.rings[,pair[2]]) # temporary pariwise dataframe 
  colnames(tmp) <- c(pair[1], pair[2])
  tree.pair.ts.list[[pair.name]] <- ts(tmp, start=ring.years[1], end=ring.years[length(ring.years)]) # convrt tmp dataframe to timeseries
  tmp.index <- index(tree.pair.ts.list[[pair.name]])[rowSums(is.na(tree.pair.ts.list[[pair.name]])) != ncol(tree.pair.ts.list[[pair.name]])] ## find index of rows where NA's are only in one of the two columns --> to remove long periods of where no data is available for either tree
  tree.pair.ts.list[[pair.name]] <- ts(tree.pair.ts.list[[pair.name]][rowSums(is.na(tree.pair.ts.list[[pair.name]])) != ncol(tree.pair.ts.list[[pair.name]]),], start=tmp.index[1], end=tmp.index[length(tmp.index)]) ## only select rows where at least one tree has an entry 
}
length(tree.pair.ts.list)

## take first difference of the pariwise tree time series
diff.tree.pair.ts.list <- mapply(FUN=function(x){diff(x)}, x=tree.pair.ts.list, SIMPLIFY = FALSE)

##############################
## find and remove pairwise comparisons of trees with themselves and of trees that do not have an overlapping growth period
## find self-correlation pairs
selfcorrelation.pair.names <- c(NULL)
for(i in 1:length(diff.tree.pair.ts.list)){
  pair.name <- names(diff.tree.pair.ts.list)[i]
  if(unlist(strsplit(pair.name,split=" + ", fixed=TRUE))[1] == unlist(strsplit(pair.name,split=" + ", fixed=TRUE))[2]){ 
    print(c(i, pair.name))
    selfcorrelation.pair.names <- append(selfcorrelation.pair.names, pair.name)
  }
}

## find tree pairs that do not have overlapping growth periods
no.shared.timeperiod.pair.names <- c(NULL)
tmp.error <- c(NULL) 
for(i in 1:length(diff.tree.pair.ts.list)){
  pair.name <- names(diff.tree.pair.ts.list)[i]
  try(tmp.error <- na.contiguous(diff.tree.pair.ts.list[[pair.name]])) #  no overlapping time-period of growth among trees
  if(is.null(tmp.error) == TRUE){
    print(c(i, pair.name))
    no.shared.timeperiod.pair.names <- append(no.shared.timeperiod.pair.names, pair.name)
  }
  tmp.error <- NULL
}

## remove pairs of self correlations and no shared time periods from analysis
diff.tree.pair.ts.list <- diff.tree.pair.ts.list[names(diff.tree.pair.ts.list)[!names(diff.tree.pair.ts.list) %in% c(selfcorrelation.pair.names, no.shared.timeperiod.pair.names)]]
length(diff.tree.pair.ts.list)
#####################################

#######################################################
## time-varying bivariat random walk plus noise state space model
## notes: seperate single Verror (observation error) for each variable (here tree time series); 
#######################################################
tvar.corr.tree.RWN.build <- function(parm, y.ts, ...){
  ## TRANSFORM to (-10, 10)
  parm <- -10 + 20 / (1 + exp(-parm))
  
  mod <- dlm(FF = 1, V = 1, GG = 1, W = 1, m0 = 0, C0 = 1e+7) ## RWN DLM with 2 observations
  mod <- dlmSum(mod, mod) # combine 2 indepedent RWN DLMs with 2 observation to bivariate RWN with 2 observations per variable
  
  ##-----------------------------
  # note: no time-variance for observations and no correlation among observations estimated to keep model simpler
  mod$V <- diag(exp(parm[c("V11", "V21")])) # seperate single Verror among variables
  
  ##-----------------------------
  #### build remaining dlm to make it time-varying at process stage 
  ## extract start and end DOY of each time series reading
  time <- 1:nrow(y.ts)
  
  ## implement JW and WX
  w <- list(NA)
  for (i in 1:length(matrix(NA,2,2)[lower.tri(matrix(NA,2,2),TRUE)])){
    w[[i]] <- i
  }
  JW <- matrix(NA,2,2)
  JW[lower.tri(JW,TRUE)] <- unlist(w)
  JW[upper.tri(JW)] <- t(JW)[upper.tri(t(JW))]
  
  mod$JW <- JW
  
  alpha.var1 <- exp(parm["a0"] + parm["a1"] * time / max(time)) 
  alpha.var2 <- exp(parm["b0"] + parm["b1"] * time / max(time)) 
  alpha.corr <- tanh(parm["z0"] + parm["delta"] * time / max(time))
  alpha.covar <- alpha.corr * sqrt(alpha.var1) * sqrt(alpha.var2)
  
  ## set WX columns
  WX <- cbind(alpha.var1, alpha.covar, alpha.var2)
  colnames(WX) <- c("alpha.var1", "alpha.covar", "alpha.var2")
  
  ## combine VX and WX into one matrix
  mod$X <- WX
  return(mod)
}

##################################################
## fit RWN model to tree pairwise comparisons
##################################################

## set starting values and apply transformation to [-10,10]
start.parm.RWN <- c(rep(log(1e-01), 2), log(0.1), 0, log(0.1), 0, 0, 0) # order of start: var1.obs var2.obs, , a0, a1, b0, b1, z0, delta1
# start.parm.RWN <- c(rep(log(1e-04), 2), log(0.1), 0, log(0.1), 0, 0, 0) # order of start: var1.obs var2.obs, , a0, a1, b0, b1, z0, delta1
names(start.parm.RWN) <- c("V11", "V21", "a0", "a1", "b0", "b1", "z0", "delta")
start.parm.RWN <- -log(20/(start.parm.RWN + 10) - 1) # working on the transformed scale between -10 and 10 (from preliminary analysis BFGS hangs in particular for Vnull of simulated data)

## fit model to all pairwise comparisons
fit.list <- list()
for (i in 1:length(diff.tree.pair.ts.list)){
  pair.name <- names(diff.tree.pair.ts.list)[[i]]
  print(c(i, pair.name))
  
  pair.ts <- diff.tree.pair.ts.list[[pair.name]] ## first differenced time series
  fit.list[[pair.name]] <- try(dlmMLE(pair.ts, y.ts=pair.ts, parm=start.parm.RWN, build = tvar.corr.tree.RWN.build, hessian=TRUE, control = c(trace=TRUE, maxit=1e4), method="BFGS"))
  
  ## code to save fit.list periodically in case of crash
  #if(i %in% seq(1, length(tree.pair.ts.list), by = 10) | i==length(tree.pair.ts.list)){
  #    saveRDS(fit.list, file=paste(getwd(), "/R/tree_ring/PISY/first_diff_PISY_Verror_start_1e-01.rds", sep=""))
  #}
}

###########
## check fits for errors and unsuccessful convergences
###########
## load saved fit.list
## fit.list <- readRDS(paste(getwd(), "/R/tree_ring/PISY/fits_first_diff_PISY_Verror_start_1e-01.rds", sep=""))
## print errors
for(i in 1:length(fit.list)){
  if(class(fit.list[[i]]) == "try-error" ){ #} | fit.list[[i]]$convergence !=0 ){
    print(names(fit.list[i]))
  }
}

## print unsuccessfull convergences
for(i in 1:length(fit.list)){
  if(fit.list[[i]]$convergence !=0 ){
    print(names(fit.list[i]))
  }
}

## model and estimations
mod.list <- list()
filt.list <- list()
smooth.list <- list()
for(i in 1:length(fit.list)){
  pair.name <- names(fit.list)[[i]]
  # if(class(fit.list[[pair.name]]) == "try-error" ){ next} # go straight to next itteration if the model didn't fit for a particular pairwise comparison
  mod.list[[pair.name]] <- try(tvar.corr.tree.RWN.build(fit.list[[pair.name]]$par, y.ts=diff.tree.pair.ts.list[[pair.name]]))
  filt.list[[pair.name]] <- dlmFilter(diff.tree.pair.ts.list[[pair.name]], mod=mod.list[[pair.name]])
  smooth.list[[pair.name]] <- dlmSmooth(filt.list[[pair.name]])
}

results <- list("fit"=fit.list, "mod"=mod.list, "filt"=filt.list, "smooth"=smooth.list)
# saveRDS(results, file=paste(getwd(), "/R/tree_ring/PISY/results_fits_mod_filt_smooth_first_diff_PISY_Verror_1e-01.rds", sep=""))
# results <- readRDS(paste(getwd(), "/R/tree_ring/PISY/results_fits_mod_filt_smooth_first_diff_PISY_Verror_1e-01.rds", sep=""))

###############################################################
## extract and plot all estimated time varying correlations
###############################################################
# fit.list <- readRDS("row_removal_if_NA_in_all_columns/results_first_diff_PISY_tree_rings_single_Verror_RWN_fits_and_estimates_tree_rings_pairwise_comparisons.rds")[["fit"]]
fit.list <- readRDS(paste(getwd(), "/R/tree_ring/PISY/results_fits_mod_filt_smooth_first_diff_PISY_Verror_1e-01.rds", sep=""))[["fit"]]
#########################################

selfcorrelation.pair.names <- c(NULL)
for(i in 1:length(fit.list)){
  pair.name <- names(fit.list)[i]
  if(unlist(strsplit(pair.name,split=" + ", fixed=TRUE))[1] == unlist(strsplit(pair.name,split=" + ", fixed=TRUE))[2]){ 
    print(c(i, pair.name))
    selfcorrelation.pair.names <- append(selfcorrelation.pair.names, pair.name)
  }
}

## find tree pairs that do not have overlapping growth periods
no.shared.timeperiod.pair.names <- c(NULL)
tmp.error <- c(NULL) 
for(i in 1:length(fit.list)){
  pair.name <- names(fit.list)[i]
  #try(tmp.error <- na.contiguous(tree.pair.ts.list[[pair.name]])) #  no overlapping time-period of growth among trees
  try(tmp.error <- na.contiguous(diff.tree.pair.ts.list[[pair.name]])) #  no overlapping time-period of growth among trees
  if(is.null(tmp.error) == TRUE){
    print(c(i, pair.name))
    no.shared.timeperiod.pair.names <- append(no.shared.timeperiod.pair.names, pair.name)
  }
  tmp.error <- NULL
}
length(no.shared.timeperiod.pair.names)

########## remove self correlations and no shared time periods from fit.list
#fit.list <- fit.list[names(fit.list)[!names(fit.list) %in% c(selfcorrelation.pair.names, no.shared.timeperiod.pair.names)]]
#length(fit.list)

## check for eigenvalues
pair.name <- names(fit.list)[597]
pars <- c("z0", "delta")
eigen(fit.list[[pair.name]]$hessian[pars,pars])$values



## apply deltaMethod to estimate correlation variance and standard error
xpred <- par.est <- vcov.hat <- corr.est.uncert.list <- X.star <- eta.hat.est <- se.eta.hat.est <- uncert.list <- corr.hat.est <- upr.95CI.corr.est <- lwr.95CI.corr.est <- delta.df <- error.delta <-  list() # create buckets
error.delta <- c(NULL)

for(i in 1:length(fit.list)){
  #for(i in 959:length(fit.list)){
  #for(i in 597:600){
  pair.name <- names(fit.list)[i]
  
  ## this part is obsolete since previously removed
  # if self correlation pairs or no shared growth period then skip this dataset
  # if(pair.name %in% selfcorrelation.pair.names | pair.name %in% no.shared.timeperiod.pair.names) {next} # do not apply deltaMethod to self_correlations or pairwise comparisons that share not common time period
  
  print(c(i, pair.name))
  pars <- c("z0", "delta")
  
  ## set xpred, if no overlapping time period for correlation  go to next iteration
  xpred[[pair.name]] <- seq(0, 1, length=nrow(na.omit(diff.tree.pair.ts.list[[pair.name]]))) # generate x prediction values only over correlation range without NAs
  
  par.est[[pair.name]] <- fit.list[[pair.name]]$par
  
  
  
  try(vcov.hat[[pair.name]] <- solve(fit.list[[pair.name]]$hessian)[pars,pars])# estimated var-covariance matrix
  if(is.null(vcov.hat[[pair.name]])){error.delta[[pair.name]] <- i} ## extract error pair.names and index number in fit.list
  if(class(vcov.hat[[pair.name]]) == "NULL"){  #test if error in solving hessian, paste error message if so and record correlation pair.name
    delta.df[[pair.name]] <- paste("Error solving hessian")
    #error.delta <- append(error.delta, pair.name)
  } else {
    # use deltaMethod() to approximate variance of non-linear function, due to bounding parameters between [-10,10]
    delta.approx <- sapply(xpred[[pair.name]], function(z){
      deltaMethod(object = fit.list[[pair.name]]$par[pars], g = "(-10 + 20 / (1 + exp(-z0))) + (-10 + 20 / (1 + exp(-delta))) * x", vcov. = vcov.hat[[pair.name]], constant = c(x = z))
    })
    
    delta.df[[pair.name]] <- as.data.frame(t(delta.approx))
  }
  # if(i == seq(1, length(fit.list), by=10) | length(fit.list)){
  #      delta.results <- list("delta.df"=delta.df, "error.delta"=error.delta, "vcov.hat"=vcov.hat, "xpred"=xpred)
  #      saveRDS(delta.results, file="row_removal_if_NA_in_all_columns/test_fits/verror_start_value_1e-04/results_delta_method_first_diff_PISY_tree_ring_correlation_Verror_start_1e-04.rds")
  # }
}


## load deltaMethod results

#results.delta <- readRDS(file="row_removal_if_NA_in_all_columns/first_diff_PISY_tree_ring_correlation_delta_method_results.rds")
results.delta <- readRDS(file="row_removal_if_NA_in_all_columns/test_fits/verror_start_value_1e-04/results_delta_method_first_diff_PISY_tree_ring_correlation_Verror_start_1e-04.rds")
delta.df <- results.delta[["delta.df"]]
error.delta <- results.delta[["error.delta"]]
vcov.hat <-  results.delta[["vcov.hat"]]
xpred <- results.delta[["xpred"]]

delta.df <- delta.df[names(delta.df)[!names(delta.df) %in% names(error.delta)]]

## check variance matrix and eigenvalues
# test.vcov <- eigen.list <- list()
# for(i in 1:length(error.delta)){
#      diff.pair.name <- error.delta[i]
#      pars <- c("z0", "delta")
#      #try(test.vcov[[diff.pair.name]] <- solve(fit.list[[diff.pair.name]]$hessian)[pars,pars])
#      eigen.list[[diff.pair.name]] <- eigen(fit.list[[diff.pair.name]]$hessian)
# }

## select all covariance matrices with negative eigenvalues which does not full fill the properties of semi-positive definite covariance matrices
## based on theorem: The real symmetric matrix covariance matrix is positive definite if and only if its eigenvalues are positive. It is positive-semi-definite if and only if its eigenvalues are nonnegative (>=0).
## & remove analyses with non-positive eigenvalues vcov
# non.positive.eigenvalues.vcov <- vcov.hat[(mapply(mapply(vcov.hat, FUN=function(x){eigen(x)[1]}), FUN=function(z){any(c(z[1]<=0, z[2]<=0))}))]
# delta.df <- delta.df[ ! names(delta.df) %in% names(non.positive.eigenvalues.vcov)] ## with non positive eigenvalues removed

#correct.names <- names(delta.df)[!names(delta.df) %in%  names(non.positive.eigenvalues.vcov)]


corr.hat <- corr.lwr <- corr.upr <- list() 
for(i in 1:length(delta.df)){
  pair.name <- names(delta.df)[i]
  if(delta.df[[pair.name]] == "Error solving hessian"){ next
  } else {
    #print(i)
    corr.hat[[pair.name]] <- tanh(unlist(delta.df[[pair.name]]$Estimate))
    corr.lwr[[pair.name]] <- tanh(unlist(delta.df[[pair.name]][, "2.5 %"]))
    corr.upr[[pair.name]] <- tanh(unlist(delta.df[[pair.name]][, "97.5 %"]))
  }
}


###########################################################################
## plot pairwise comparisons by region and location on DOY scale
## create data frame for plotting correlations of pairwise RWN fitting
tree1 <- list() ## bucket for name of tree sample 1 of each pairwise RWN fitting
tree2 <- list() ## bucket for name of tree sample 2 of each pairwise RWN fitting
location1 <- list() ## bucket for location of tree sample 1 of each pairwise RWN fitting
location2 <- list() ## bucket for locationof tree sample 2 of each pairwise RWN fitting
region1 <- list() ## bucket for region of tree 1 of each pairwise RWN fitting
region2 <- list() ## bucket for region of tree 2 of each pairwise RWN fitting
# bucket for extracted estimated correlation parameters of pairwise RWN fits and correlation start, end DOY and DOY for the correlations
# z0 <- list()
# delta <- list()
corr.start <- list()
corr.end <- list()
corr.DOY <- list()

for(i in 1:length(delta.df)){
  pair.name <- names(delta.df)[i] 
  
  tree1[[pair.name]] <- unlist(strsplit(pair.name, split=" + ", fixed=TRUE))[1]
  tree2[[pair.name]] <-  unlist(strsplit(pair.name, split=" + ", fixed=TRUE))[2]
  
  location1[[pair.name]] <- if(grepl("I", tree1[[pair.name]])){"Inverey"} else if(grepl("LM", tree1[[pair.name]])){"Loch Maree"} else if(grepl("FW", tree1[[pair.name]])){"Franchise Wood"} else if(grepl("C", tree1[[pair.name]])){"Coulin"} else if(grepl("GD", tree1[[pair.name]])){"Glen Derry"}
  location2[[pair.name]] <- if(grepl("I", tree2[[pair.name]])){"Inverey"} else if(grepl("LM", tree2[[pair.name]])){"Loch Maree"} else if(grepl("FW", tree2[[pair.name]])){"Franchise Wood"} else if(grepl("C", tree2[[pair.name]])){"Coulin"} else if(grepl("GD", tree2[[pair.name]])){"Glen Derry"}
  
  region1[[pair.name]] <- if(grepl("I", tree1[[pair.name]])){"NorthEast"} else if(grepl("LM", tree1[[pair.name]])){"NorthWest"} else if(grepl("FW", tree1[[pair.name]])){"South"} else if(grepl("C", tree1[[pair.name]])){"NorthWest"} else if(grepl("GD", tree1[[pair.name]])){"NorthEast"}
  region2[[pair.name]] <- if(grepl("I", tree2[[pair.name]])){"NorthEast"} else if(grepl("LM", tree2[[pair.name]])){"NorthWest"} else if(grepl("FW", tree2[[pair.name]])){"South"} else if(grepl("C", tree2[[pair.name]])){"NorthWest"} else if(grepl("GD", tree2[[pair.name]])){"NorthEast"}
  
  corr.start[[pair.name]] <-  start(na.contiguous(diff.tree.pair.ts.list[[pair.name]]))[1] # na.contiguous pair ts
  corr.end[[pair.name]] <-  end(na.contiguous(diff.tree.pair.ts.list[[pair.name]]))[1] # na.contiguous pair ts
  corr.DOY[[pair.name]] <- time(na.contiguous(diff.tree.pair.ts.list[[pair.name]])) # na.contiguous pair ts
}

## convert to factors
region1 <- factor(unlist(region1))
region2 <- factor(unlist(region2))
location1 <- factor(unlist(location1), levels=c("Loch Maree", "Coulin", "Inverey", "Glen Derry", "Franchise Wood"))
location2 <- factor(unlist(location2), levels=c("Loch Maree", "Coulin", "Inverey", "Glen Derry", "Franchise Wood"))

## sort regions and locations to get all the plots into one corner of the correlation matrix
regions.df <- data.frame(region1, region2)
regions.switched <- t(apply(regions.df, 1, sort))

# switch order of locations if order of region has been switched in the previous step
locations.df <- data.frame(factor(unlist(location1)), factor(unlist(location2)))
test.order <- cbind(regions.df[,1], regions.switched[,1], locations.df)
locations.switched <- t(apply(as.matrix(test.order), 1, FUN = function(x){if(x[1] != x[2]){rev(x[3:4])} else{x[3:4]}})) ## locations are switched if regions were switched

prop.conf.interv <- xpred.DOY <- time <-  list()
for(i in 1:length(corr.hat)){
  pair.name <- names(corr.hat)[i]
  #pair.name <-correct.names[i] 
  ## 95% confidence interval width as a proportion of the correlation width of 2
  prop.conf.interv[[pair.name]] <- (corr.upr[[pair.name]] - corr.lwr[[pair.name]]) / 2
  xpred.DOY[[pair.name]] <- as.vector(time(na.contiguous(diff.tree.pair.ts.list[[pair.name]]))) ##DOY only over range of correlation without NAs in either time series
  time[[pair.name]] <- 1:nrow(na.contiguous(diff.tree.pair.ts.list[[pair.name]]))
}

weighted.est.corr.df <- data.frame("pair.name"=names(corr.hat), "tree1"=unlist(tree1), "tree2"=unlist(tree2), "region1"=regions.switched[,1], "region2"=regions.switched[,2], "location1"=locations.switched[,1], "location2"=locations.switched[,2], "corr.start"=unlist(corr.start), "corr.end"=unlist(corr.end))

weighted.est.corr.long <- weighted.est.corr.df[rep(1:length(xpred.DOY), times=as.vector(sapply(xpred.DOY, FUN=function(x){length(x)}))-1), c("pair.name", "tree1", "tree2", "region1", "region2", "location1", "location2", "corr.start", "corr.end")]
weighted.est.corr.long$DOY <- unlist(mapply(corr.DOY, FUN=function(x){x[-1]}))
weighted.est.corr.long$corr <- unlist(mapply(corr.hat, FUN=function(x){x[-1]}))
weighted.est.corr.long$xstart <- unlist(mapply(xpred.DOY, FUN=function(x){x[-length(x)]}))
weighted.est.corr.long$xend <- unlist(mapply(xpred.DOY, FUN=function(x){x[-1]}))
weighted.est.corr.long$ystart <- unlist(mapply(corr.hat, FUN=function(x){x[-length(x)]}))
weighted.est.corr.long$yend <- unlist(mapply(corr.hat, FUN=function(x){x[-1]}))
weighted.est.corr.long$prop.line.width <- unlist(mapply(prop.conf.interv, FUN=function(x){x[-1]}))


## correlation matrix plot by region and location
theme_set(theme_bw())
weighted.tree.ring.corr.plot <- ggplot(weighted.est.corr.long, aes(x=xstart, y=ystart, xend=xend, yend=yend, size=1-prop.line.width)) +#, colour=grey(prop.line.width))) + #
  geom_segment() +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_grid(region1 + location1 ~ region2 + location2) +
  scale_size_continuous(name="Line width inverse porportional to 95%CI width", range = c(0.005, 0.5), limits=c(0,1)) + # 
  scale_colour_grey(guide=FALSE) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Correlation",limits=c(-1,1)) +
  #guides(colour=guide_colorbar()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=10),
        axis.text.y = element_text(size=10),
        legend.position ="bottom") 

weighted.tree.ring.corr.plot 

# pdf("/Users/Staeudchen/Documents/Work_and_Volunteering/Ireland/GMIT/GMIT_Reserach_Assistant_State_Space_Modeling_2016/article_draft/LaTex_draft/figures/QUPE_width_and_col_weighted_tree_ring_single_Verror_timevar_correlation_plot.pdf", height = 10, width = 10*1.618)
# pdf("PISY_width_and_col_weighted_tree_ring_single_Verror_timevar_correlation_plot.pdf", height = 10, width = 10*1.618)
# weighted.tree.ring.corr.plot
# dev.off()

# add data for braod trend calculation via smoothing
weighted.est.corr.long$zcorr <- atanh(weighted.est.corr.long$corr)
weighted.est.corr.long$time <- unlist(mapply(time, FUN=function(x){x[-1]}))

library(mgcv)
## include weights on the smooth using inverse of proportion
#gam.fit <- gam(zcorr ~ s(time), data = weighted.est.corr.long)
gam.fit <- gam(zcorr ~ s(DOY,location1), data = weighted.est.corr.long)
#gam.fit <- gam(zcorr ~ s(time*(1/prop.line.width)), data = weighted.est.corr.long)

pred.df <- data.frame(DOY = weighted.est.corr.long$DOY)
rownames(pred.df) <- rownames(weighted.est.corr.long)
# time.range <- min(weighted.est.corr.long$time):max(weighted.est.corr.long$time)
# pred.df <- data.frame(time = time.range)
#pred.df <- data.frame(time = unlist(time))

## on (-1,1) scale
pred.df$pred <- tanh(predict(gam.fit, newdata = pred.df))

theme_set(theme_bw())
weighted.and.trend.tree.corr.plot <- ggplot(data=weighted.est.corr.long) +
  geom_segment(aes(x=xstart, y=ystart, xend=xend, yend=yend, size=1-prop.line.width, colour=grey(prop.line.width))) +
  geom_smooth(aes(x = DOY, y = corr, weight=1-prop.line.width), colour = "red", size=0.5) +
  #geom_line(data = pred.df, aes(x=DOY, y = pred), colour = "red") +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_grid(region1 + location1 ~ region2 + location2) +
  scale_size_continuous(name="Estimated correlation line width inverse porportional to 95%CI width", range = c(0.005, 0.5), limits=c(0,1)) + # 
  scale_colour_grey(guide=FALSE) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Correlation",limits=c(-1,1)) +
  #guides(colour=guide_colorbar()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x  = element_text(angle=90, vjust=0.5, size=10),
        axis.text.y = element_text(size=10),
        legend.position ="bottom") 

# pdf("/Users/Staeudchen/Documents/Work_and_Volunteering/Ireland/GMIT/GMIT_Reserach_Assistant_State_Space_Modeling_2016/article_draft/LaTex_draft/figures/PISY_weighted_and_broad_trend_tree_ring_corr_plot.pdf", height = 10, width = 10*1.618)
# pdf("PISY_weighted_and_broad_trend_tree_ring_corr_plot.pdf", height = 10, width = 10*1.618)
# weighted.and.trend.tree.corr.plot
# dev.off()

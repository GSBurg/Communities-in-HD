#### Scripts for the paper: Communities in high definition: spatial and environmental factors shape the microdistribution of aquatic invertebrates ####
## by Burgazzi et al. 2020
#  https://doi.org/10.1111/fwb.13599

# The original dataset used in the paper has been collected in Parma, Enza and Nure streams (in Northern Italy)
# Each stream has been sampled twice (one sampling in summer and one sampling in winter), resulting in a total of 300 macroinvertebrate samples
# For each campaign macroinvertebrate samples have been collected according to 50-points random sampling grids
# The macroinvertebrate dataset (macro) is a 300x78 community matrix
# The covariate dataset (env) is composed by 4 environmental numerical variables (flow velocity, water depth, substrate size and benthic organic matter) and by spatial coordinates within each grid


## Import the datasets===============================================
env <- read.table("env.txt", head=T) #dataset of environmental variables, spatial coordinates and invertebrate metrics
macro <- read.table("macro.txt", head=T) #community matrix
xy <- env[, c(3,4)] #create a dataset of coordinates by subsetting the environmental dataset
env$campaign <- as.factor(paste(env$river, env$season, sep="_"))


## GAM models========================================================
library(mgcv)

# the following model is for taxon richness, for the other models just replace "ntaxa" in the formula with the name of the other macroinvertebrate metrics
# the smoothing function is applied to coordinates
# the "by" option allow to produce a different smooth for each factor level

metric.gam <- gam(log(ntaxa) ~ campaign + BOM + vel + depth + sub + 
                   s(x, y, by = campaign), 
                 data=env, method = "REML", family="gaussian")

summary(metric.gam)
anova(metric.gam)


## HMSC models ======================================================
# the following scripts are modified from Tikhonov et al. (2020) "Joint species distribution modelling with the r-package Hmsc"
# https://doi.org/10.1111/2041-210X.13345

## load libraries
set.seed(1)
source("load_libraries.r")
source("evaluateModelFit.r")

## assign the datasets
Y <- macro
X <- env[,5:8] #create a dataset for environmental covariates
grid_coords <- env[,3:4]

# a different HMSC model is done for each sampling campaign
# subsets for each sampling campaign is created from the global datasets
# in the following example is done for the first campaign (Parma stream during summer)
# the models for the other campaigns can be done by changing the row numbers in "camp"

camp <- 1:50 #this is parma summer campaign. Change these row numbers for the other campaigns
Y_camp <- Y[camp,] 
Y_camp <- Y_camp[ , colSums(Y_camp) > 0 ]
X_camp <- as.data.frame(scale(X[camp,]))
grid_coords_camp <- grid_coords[camp,]

## exclusion of rare taxa
abu_taxa_camp <- apply(Y_camp, 2, sum)
abu_perc_camp <- (abu_taxa_camp/sum(abu_taxa_camp))*100
Y_camp_red <- Y_camp[,which(abu_perc_camp>0.5)]

## study design and random effect structure
studyDesign_camp <- data.frame(replica=paste("r",c(camp), sep=""))
sRL_camp <- grid_coords_camp
rownames(sRL_camp) <- studyDesign_camp[,1]
rL_camp <- HmscRandomLevel(sData = sRL_camp) #sData argument for spatial data
rL_camp <- setPriors(rL_camp, nfMin = 5, nfMax = 10)

## construct and fit the model
thin = 100
samples = 1000
nChains = 4
Y_camp_red <- as.matrix(Y_camp_red) 
Hmsc_model_camp_red <- Hmsc(Y=Y_camp_red, XFormula=~vel+depth+sub+BOM, XData=X_camp, XScale=FALSE,
                          distr = "lognormal poisson", studyDesign = studyDesign_camp, ranLevels = list(replica=rL_camp))

# be careful, time consuming
Hmsc_model_camp_red <-  sampleMcmc(Hmsc_model_camp_red, samples = samples, thin=thin, adaptNf=rep(ceiling(0.4*samples*thin),1), transient = ceiling(0.5*samples*thin), nChains = nChains)

## compute mixing statistics
mpost_camp_red = convertToCodaObject(Hmsc_model_camp_red)

es.beta_camp_red = effectiveSize(mpost_camp_red$Beta)
ge.beta_camp_red = gelman.diag(mpost_camp_red$Beta,multivariate=FALSE)$psrf

es.gamma_camp_red = effectiveSize(mpost_camp_red$Gamma)
ge.gamma_camp_red = gelman.diag(mpost_camp_red$Gamma,multivariate=FALSE)$psrf

es.V_camp_red = effectiveSize(mpost_camp_red$V)
ge.V_camp_red = gelman.diag(mpost_camp_red$V,multivariate=FALSE)$psrf

es.omega_camp_red = effectiveSize(mpost_camp_red$Omega[[1]])
ge.omega_camp_red = gelman.diag(mpost_camp_red$Omega[[1]],multivariate=FALSE)$psrf

mixing_camp_red = list(es.beta=es.beta_camp_red, ge.beta=ge.beta_camp_red,
                     es.gamma=es.gamma_camp_red, ge.gamma=ge.gamma_camp_red,
                     es.V=es.V_camp_red, ge.V=ge.V_camp_red,
                     es.omega=es.omega_camp_red, ge.omega=ge.omega_camp_red)
mixing_camp_red

## compute model fit
m <- Hmsc_model_camp_red
predY_camp_red = computePredictedValues(Hmsc_model_camp_red, expected=FALSE)
MF_camp_red = evaluateModelFit(hM=m, predY=predY_camp_red)
MF_camp_red$SR2$A
MF_camp_red$O.AUC$A

## compute model fit with 5 folds cross-validation
m <- Hmsc_model_camp_red
partition_camp_CV_red=createPartition(hM=m, nfolds=5, column="replica")
# be careful, time consuming
predY_camp_CV_red = computePredictedValues(m, expected=FALSE, partition=partition_camp_CV_red)
MF_camp_CV_red = evaluateModelFit(hM=m, predY=predY_camp_CV_red)
MF_camp_CV_red$SR2$A
MF_camp_CV_red$O.AUC$A

## plot variance partitioning
group = c(1,2,3,4,5)
groupnames = c("int","vel","depth","sub","BOM")
VP_camp_red = computeVariancePartitioning(Hmsc_model_camp_red, group=group,groupnames=groupnames)

mean_camp <- apply(VP_camp_red$vals[-1,], 1, mean)
leg <- c(paste("vel (mean =", round(mean_camp[1], digits = 3)*100, "%)"),
         paste("depth (mean =", round(mean_camp[2], digits = 3)*100, "%)"),
         paste("sub (mean =", round(mean_camp[3], digits = 3)*100, "%)"),
         paste("BOM (mean =", round(mean_camp[4], digits = 3)*100, "%)"),
         paste("space (mean =", round(mean_camp[5], digits = 3)*100, "%)")
         )

par(mar=c(4,9,7,1))
barplot(VP_camp_red$vals[-1,], horiz=TRUE, las=1, col=c("navy","green3","gray50","red2","beige") , 
        border="black", xlab="", legend=leg,
        args.legend=list(xjust=1.25, yjust=0, ncol=1, cex=1.09, bty="n")
        )

## plot association networks
OmegaCor_camp_red = computeAssociations(Hmsc_model_camp_red)
supportLevel = 0.90

plotOrder_camp_red = corrMatOrder(OmegaCor_camp_red[[1]]$mean,order="AOE")
toPlot_camp_red = ((OmegaCor_camp_red[[1]]$support>supportLevel) + (OmegaCor_camp_red[[1]]$support<(1-supportLevel))>0)*OmegaCor_camp_red[[1]]$mean
corrplot(toPlot_camp_red[plotOrder_camp_red,plotOrder_camp_red], diag = F, type = "lower", tl.col="black", tl.cex=0.7, cex.main=1.3, tl.srt = 45, mar=c(0.1,0.1,1,0.1), col = colorRampPalette(c("blue", "white", "red"))(200))

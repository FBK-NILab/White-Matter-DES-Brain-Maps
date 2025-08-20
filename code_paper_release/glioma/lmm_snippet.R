# http://dpmartin42.github.io/posts/Piecewise-growth
# https://www.lexjansen.com/pharmasug-cn/2015/ST/PharmaSUG-China-2015-ST08.pdf
# https://joshuawiley.com/MonashHonoursStatistics/LMM_Comparison.html#effect-sizes
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13434


setwd("C:/Users/ludo2/OneDrive/Desktop/Post_doc_FBK/DES_SUBCORTICAL_PAPER/LMM")
library(lme4)
library(tidyverse)
library(performance)
library(partR2)

datavar <- read.delim("SEMANTIC_lmm_weighted.csv",sep = ',')

# We are going to implement a piecewise linear regression model centered at t1 (post one week). We will model t0-t1 separately from t1-t2
# To do this, we will add two (dummy) colums to our df. Important: t1 will be zero in both columns

dummy_t0=rep(0,length(datavar$OCCASION))
dummy_t2=rep(0,length(datavar$OCCASION))
dummy_t0[datavar$OCCASION==0]=1
dummy_t2[datavar$OCCASION==2]=1
datavar$t0=dummy_t0
datavar$t2=dummy_t2

# WM_overlap:

# pairwise comparisons

# SURGERY_TYPE, AGE, SEX, EDUCATION, LESION_CM3, MGMT, IDH, GRADE, LOBE, SIDE

plgModel_1 <- lmer(BEHAV_OUT ~ 1 + WM_OVERLAP*t0 + WM_OVERLAP*t2 + (1 | SUB_ID),
                   REML=FALSE,data = datavar)
plgModel_2 <- lmer(BEHAV_OUT ~ 1 + SIDE*t0 + SIDE*t2 + (1 | SUB_ID),
                   REML=FALSE,data = datavar)
plgModel_3 <- lmer(BEHAV_OUT ~ 1 + SIDE*t0 + SIDE*t2 + WM_OVERLAP*t0 + WM_OVERLAP*t2 + (1 | SUB_ID),
                   REML=FALSE,data = datavar)

aa <- model_performance(plgModel_3)
bb <- model_performance(plgModel_2)
cc <- model_performance(plgModel_1)

# effect size for other covariate
(aa$R2_marginal-bb$R2_marginal)/(1-aa$R2_marginal)

# effect size for wm_overlap
(aa$R2_marginal-cc$R2_marginal)/(1-aa$R2_marginal)

# Check assumptions

plot(plgModel_2)
rfx = data.frame(ranef(plgModel_2)) # edit model
shapiro.test(rfx$condval)
shapiro.test(datavar$BEHAV_OUT - fitted(plgModel_2)) # edit model


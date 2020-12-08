# 11/20/2020
# use other functions to get results
library(survival)
library(magrittr)
library(ggplot2)

#### continuous outcome simulations ####
source("~/Box/Penn/research/trials/code/public/contsims_function.R")
gam04pow08_rand = twoarmcont(NSIMS = 5000, n = seq(190, 420, by = 10),
                             gam = 0.4, pow = 0.8, rand = T) #390 210
gam04pow09_rand = twoarmcont(NSIMS = 5000, n = seq(200, 550, by = 10),
                             gam = 0.4, pow = 0.9, rand = T) #530 270

gam05pow08_rand = twoarmcont(NSIMS = 5000, n = seq(130, 270, by = 10),
                             gam = 0.5, pow = 0.8, rand = T) #260 130
gam05pow09_rand = twoarmcont(NSIMS = 5000, n = seq(150, 350, by = 10),
                             gam = 0.5, pow = 0.9, rand = T) #340 170

gam06pow08_rand = twoarmcont(NSIMS = 5000, n = seq(80, 200, by = 10),
                              gam = 0.6, pow = 0.8, rand = T) #180 90
gam06pow09_rand = twoarmcont(NSIMS = 5000, n = seq(110, 250, by = 10),
                              gam = 0.6, pow = 0.9, rand = T) # 240 130

gam075pow08_rand = twoarmcont(NSIMS = 5000, n = seq(50, 120, by = 1),
                              gam = 0.75, pow = 0.8, rand = T) #112 61
gam075pow09_rand = twoarmcont(NSIMS = 5000, n = seq(75, 160, by = 1),
                              gam = 0.75, pow = 0.9, rand = T) #151 78


#### survival outcome simulations ####
source("~/Box/Penn/research/trials/code/public/survsims_function.R")
set.seed(12345)
gam01pow08_rand = powercomp_surv(NSIMS = 1000, n = seq(200, 600, by= 10), 
                                 gam = 0.1, pow = 0.8, rand = T) #510 210
gam01pow09_rand = powercomp_surv(NSIMS = 5000, n = seq(280, 800, by= 10), 
                                 gam = 0.1, pow = 0.9, rand = T) #690 290

gam015pow08_rand = powercomp_surv(NSIMS = 5000, n = seq(50, 240, by= 10), 
                                 gam = 0.15, pow = 0.8, rand = T) #230 100
gam015pow09_rand = powercomp_surv(NSIMS = 5000, n = seq(100, 350, by= 10), 
                                 gam = 0.15, pow = 0.9, rand = T) #320 140

gam02pow08_rand = powercomp_surv(NSIMS = 5000, n = seq(40, 130, by= 1), 
                          gam = 0.2, pow = 0.8, rand = T) #128 58
gam02pow09_rand = powercomp_surv(NSIMS = 5000, n = seq(75, 200, by= 1), 
                                 gam = 0.2, pow = 0.9, rand = T) #173 78

#### gbm surv simulations ####
gbmsurvdat = read.csv("~/Box/Penn/research/trials/data/surv_full_02062019.csv")
gbmsurvdat$cens <- ifelse(gbmsurvdat$Survival < 36, 0, 1) #censoring > 3-year
gbmsurvdat$surv3 <- ifelse(gbmsurvdat$Survival < 36, gbmsurvdat$Survival, 36)
source('~/Box/Penn/research/trials/code/public/gbmsurvsims_function.R')

set.seed(12345)
gam165pow08_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(60, 134, by= 1), dat = gbmsurvdat, 
                                        gam = 1.65, pow = 0.8,  rand = T) # 109 64
gam165pow09_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(60, 134, by= 1), dat = gbmsurvdat, 
                                        gam = 1.65, pow = 0.9,  rand = T) # doesn't reach, 87

gam17pow08_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(50, 120, by= 1), dat = gbmsurvdat, 
                                        gam = 1.7, pow = 0.8,  rand = T) # 94 58
gam17pow09_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(80, 134, by= 1), dat = gbmsurvdat, 
                                        gam = 1.7, pow = 0.9,  rand = T) # 126 80

gam18pow08_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(45, 85, by= 1), dat = gbmsurvdat, 
                                        gam = 1.8, pow = 0.8,  rand = T) # 79 48
gam18pow09_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(45, 115, by= 1), dat = gbmsurvdat, 
                                        gam = 1.8, pow = 0.9,  rand = T) # 102 63

gam19pow08_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(35, 80, by= 1), dat = gbmsurvdat, 
                                        gam = 1.9, pow = 0.8,  rand = T) # 67 40 
gam19pow09_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(45, 95, by= 1), dat = gbmsurvdat, 
                                        gam = 1.9, pow = 0.9,  rand = T) # 87 53

gam20pow08_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(35, 65, by= 1), dat = gbmsurvdat, 
                                             gam = 2, pow = 0.8,  rand = T) # 57 35
gam20pow09_gbmsurv = powercomp_gbm_surv(NSIMS = 5000, n = seq(35, 85, by= 1), dat = gbmsurvdat, 
                                        gam = 2, pow = 0.9,  rand = T) # 76 46

#### adni surv simulations ####
source("~/Box/Penn/research/trials/code/adnisurv_powercomp_function.R")
dat = read.csv("~/Box/Penn/research/trials/data/adni_mci_surv.csv")


# 80% power
gam175_pow08_adnisurv = powercomp_surv_adni(NSIMS = 10000, n = seq(260, 300, by = 1), 
                                   gam = 1.75, pow = 0.8, data = dat) #296 262
gam18_pow08_adnisurv = powercomp_surv_adni(NSIMS = 5000, n = seq(230, 260, by = 1), 
                                  gam = 1.8, pow = 0.8, data = dat) #260 233
gam19_pow08_adnisurv = powercomp_surv_adni(NSIMS = 5000, n = seq(160, 200, by = 1), 
                                  gam = 1.9, pow = 0.8, data = dat) #193 171
gam20_pow08_adnisurv = powercomp_surv_adni(NSIMS = 5000, n = seq(130, 160, by = 1), 
                                  gam = 2.0, pow = 0.8, data = dat) #145 135

#90% power
gam175_pow09_adnisurv = powercomp_surv_adni(NSIMS = 1000, n = seq(260, 330, by = 10), 
                                   gam = 1.75, pow = 0.9, data = dat) # doesn't reach
gam18_pow09_adnisurv = powercomp_surv_adni(NSIMS = 10000, n = seq(300, 330, by = 1), 
                                  gam = 1.8, pow = 0.9, data = dat) # doesn't reach
gam19_pow09_adnisurv = powercomp_surv_adni(NSIMS = 5000, n = seq(210, 270, by = 1), 
                                  gam = 1.9, pow = 0.9, data = dat) #256  264 228
gam20_pow09_adnisurv = powercomp_surv_adni(NSIMS = 5000, n = seq(160, 210, by = 1), 
                                  gam = 2.0, pow = 0.9, data = dat) #193 175


#### adni cont simulations ####
dat = read.csv("~/Box/Penn/research/trials/data/adnimem_2year.csv")
source('~/Box/Penn/research/trials/code/public/adnimem_function.R')

# pow = 0.8
gam18_pow08 = powercomp_adnicont(NSIMS = 5000, n = seq(195, 250, by = 1), 
                                 gam = 0.18, pow = 0.8, data = dat) # 238 197
gam19_pow08 = powercomp_adnicont(NSIMS = 1000, n = seq(150, 283, by = 10), 
                                 gam = 0.2, pow = 0.8, data = dat) # 
gam02pow08 = powercomp_adnicont(NSIMS = 5000, n = seq(150, 200, by = 1), 
                                gam = 0.2, pow = 0.8, data = dat) #189 164
gam025pow08 = powercomp_adnicont(NSIMS = 5000, n = seq(100, 150, by = 1), 
                                 gam = 0.25, pow = 0.8, data = dat) #125 107
gam03pow08 = powercomp_adnicont(NSIMS = 5000, n = seq(75, 100, by = 1), 
                                gam = 0.3, pow = 0.8, data = dat) #88 76
gam05pow08 = powercomp_adnicont(NSIMS = 5000, n = seq(30, 40, by = 1), 
                                gam = 0.5, pow = 0.8, data = dat) #36 31

# pow = 0.9
gam18pow09_adnicont = powercomp_adnicont(NSIMS = 5000, n = seq(200, 283, by = 1), 
                                gam = 0.18, pow = 0.9, data = dat, boot = F) #
gam02pow09 = powercomp_adnicont(NSIMS = 5000, n = seq(200, 260, by = 1), 
                                  gam = 0.2, pow = 0.9, data = dat, boot = F) #258 212
gam025pow09 = powercomp_adnicont(NSIMS = 5000, n = seq(130, 200, by = 1), 
                                 gam = 0.25, pow = 0.9, data = dat, boot = F) #171 142
gam03pow09 = powercomp_adnicont(NSIMS = 5000, n = seq(70, 130, by = 1), 
                                gam = 0.3, pow = 0.9, data = dat) #118 100
gam05pow09 = powercomp_adnicont(NSIMS = 5000, n = seq(40, 50, by = 1), 
                                gam = 0.5, pow = 0.9, data = dat) # 47 40

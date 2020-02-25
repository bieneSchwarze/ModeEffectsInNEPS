##############################################################
##############################################################
##                                                          ##
## Analyze differences in competences measured in different ##
##   modes PBA, CBA, and WBA.                               ##
##                                                          ##
## Differences are studied with respect to                  ##
## (i)   self-selection into modes,                         ##
## (ii)  differential test functioning, and                 ##
## (iii) mode-specific prediction error.                    ## 
##                                                          ##
## authors: Sabine & Timo                                   ##
## date: 2020-02-23                                         ## 
##                                                          ##
##############################################################
##############################################################
# Note: The assumption is that students are randomly assigned to groups 
#       who were initally invited for group testing and to groups who 
#       were initially invited for online testing.



# empty work space
rm(list=ls())

# curren version 
version <- "20200223"



# ----------------------------------------------------------------------------
# Load libraries
# ----------------------------------------------------------------------------

library(mice)
library(readstata13)
library(haven)
library(dplyr)
library(doBy)
library(TAM)
library(mirt)
library(MASS)
library(statmod)
library(ggplot2)
library(rms)
library(msm)

# load additional functions for DTF analyses
source("Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/syntax/dtf.R")

# load prepared data
load(paste0("Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/results/_dataPrepared-", version, ".RData"))



# ------------------------------------------------------
# Conduct analyses for self-selection into modes, 
#    differential test functioning, and 
#    mode-specific prediction error.
# ------------------------------------------------------

modMList <- vector(mode="list", length=imp1$m)
modMainEff_PBA_List <- vector(mode="list", length=imp1$m)
modMainEff_CBA_List <- vector(mode="list", length=imp1$m)
modMainEff_WBA_List <- vector(mode="list", length=imp1$m)
pseudoMR2 <- rep(NA,imp1$m)
modG1List <- vector(mode="list", length=imp1$m)
modG2List <- vector(mode="list", length=imp1$m)
modG3List <- vector(mode="list", length=imp1$m)
modG4List <- vector(mode="list", length=imp1$m)
pseudoR2 <- matrix(NA, nrow=4, ncol=imp1$m)
psychList <- vector(mode="list", length=imp1$m)
predList_sc <- list(pred10 = vector(mode="list", length=imp1$m),
                    pred11 = vector(mode="list", length=imp1$m),
                    pred20 = vector(mode="list", length=imp1$m),
                    pred21 = vector(mode="list", length=imp1$m),
                    pred30 = vector(mode="list", length=imp1$m),
                    pred31 = vector(mode="list", length=imp1$m))
r2_sc <- list(pred10 = rep(NA,imp1$m),
              pred11 = rep(NA,imp1$m),
              pred20 = rep(NA,imp1$m),
              pred21 = rep(NA,imp1$m),
              pred30 = rep(NA,imp1$m),
              pred31 = rep(NA,imp1$m))

# run over all imputed samples
for(i in 1:imp1$m){
  
  cat("It:", i, "\n")
  datIm <- complete(imp1,i)
  datIm$birthYear <- as.factor(ifelse(datIm$birthYear <= 1989, 0, ifelse(datIm$birthYear<=1990, 1, 2))) # otherwise: models cannot be fitted

  # -----------------------------
  # 1. Self-selection into modes? (only possible for PBA, CBA, and WBA. Not for switchers.)
  # -----------------------------  
  
  datM <- datIm
  datM$mode <- as.factor(ifelse(datIm$WBA_E %in% 1, 0, ifelse(datIm$PBA_E %in% 1, 1, ifelse(datIm$CBA_E %in% 1, 2, NA)))) # mode: WBA: 0, PBA: 1, CBA: 2
  datM$part <- ifelse(datIm$PBA_R %in% 1 | datIm$CBA_R %in% 1 | datIm$WBA_R %in% 1, 1, 0)
  namM <- colnames(datM)
  namM <- namM[!(namM %in% c("ID_t", "ID_i", "w1CompPart", "w2PartCAWI", "w3PartCATI", "w4PartCAWI",
                             "PBA_E", "PBA_R", "CBA_E", "CBA_R", "WBA_E", "WBA_R", "WBA_switch_E", "WBA_switch_R",
                             "GPA_W4", "GPA_W6", "part", "mode", "uniFunding", "kidsHH_W1", "kidsHH_W3",
                             "studydropout", "helpless", "selfconcept"))] 
            # "funding of insitution" cannot be used: in private institutions no CBA 
            # "kids in HH" in Waves 1 and 3 are only used for imputation & deriving propensity scores
  form <- paste("part ~",paste(paste(namM,"*mode", sep=""), collapse="+"), sep="")
  modM <- lrm(formula=as.formula(form), dat=datM) 
  pseudoMR2[i] <- modM$stats["R2"]
  modMList[[i]] <- cbind(modM$coefficients, sqrt(diag(modM$var)))  
  
  # Get main effects for WBA, PBA and CBA
  # WBA
  modMainEff_WBA_List[[i]] <-  modMList[[i]][c(2,5:30),]
  # PBA and WBA
  # (i) Add main and interaction effects properly 
  p1_pba <- modMList[[i]][c(2,5:30),1] + modMList[[i]][seq(from=31, to=83, by=2),1] # PBA
  p1_cba <- modMList[[i]][c(2,5:30),1] + modMList[[i]][seq(from=32, to=84, by=2),1] # cBA
  # (ii) Make dummies in data set & give them names the "deltamethod" function in R can read, i.e. x1 to x81
  datM_temp <- datM[,namM]
  datM_dum <- cbind(datM_temp[, c(1:10)], model.matrix(~datM_temp$hhSize_W5)[,-1],model.matrix(~datM_temp$DegreeAbroad_noDegree)[,-1],
                    datM_temp[, c(13,14)], model.matrix(~datM_temp$birthYear)[,-1], datM_temp[, c(16:21)],
                    model.matrix(~datM_temp$uniType-1)[,-1], datM_temp[,"mas1_sc1u"]) 
  colnames(datM_dum)[11] <- "hhSize_W5_1"
  colnames(datM_dum)[12] <- "hhSize_W5_2"
  colnames(datM_dum)[13] <- "DegreeAbroad_noDegree1_1"
  colnames(datM_dum)[14] <- "DegreeAbroad_noDegree1_2"
  colnames(datM_dum)[17] <- "birthYear_1"
  colnames(datM_dum)[18] <- "birthYear_2"
  colnames(datM_dum)[25] <- "uniType_1"
  colnames(datM_dum)[26] <- "uniType_2"
  colnames(datM_dum)[27] <- "mas1_sc1u"  
  namXX <- cbind(colnames(datM_dum), paste("x",1:ncol(datM_dum), sep=""))
  colnames(datM_dum) <- namXX[,2] 
  datM_dum_mode <- datM_dum
  for(k in 1:ncol(datM_dum)){
    namK <- paste("x",ncol(datM_dum)+k,sep="")
    nn <- as.numeric(as.character(datM_dum[,k]))*ifelse(as.numeric(as.character(datM$mode))==1,1,0)
    datM_dum_mode <- cbind(datM_dum_mode,nn)
    colnames(datM_dum_mode)[ncol(datM_dum_mode)] <- namK
    namXX <- rbind(namXX,c(paste(namXX[k,1], " * mode=1", sep=""), namK))
  }
  rm(k)
  for(k in 1:ncol(datM_dum)){
    namK <- paste("x",nrow(namXX)+1,sep="")
    nn <- as.numeric(as.character(datM_dum[,k]))*ifelse(as.numeric(as.character(datM$mode))==2,1,0)
    datM_dum_mode <- cbind(datM_dum_mode,nn)
    colnames(datM_dum_mode)[ncol(datM_dum_mode)] <- namK
    namXX <- rbind(namXX,c(paste(namXX[k,1], " * mode=2", sep=""), namK))
  }  
  datM_dum_mode <- cbind(datM[,"part"], datM_dum_mode)
  colnames(datM_dum_mode)[1] <- "y"
  colnames(datM_dum_mode) <- make.unique(colnames(datM_dum_mode))
  for(k in 1:ncol(datM_dum_mode)){
    datM_dum_mode[,k] <- as.numeric(as.character(datM_dum_mode[,k])) # remove factor type,otherwise there will be an error since x11 may also be a result of factor value 1 for x1 alternatively it is the names of x11
  }
  # Estimate model with dummies anew
  form <- paste(colnames(datM_dum_mode)[1], paste("~",paste(colnames(datM_dum_mode)[-1], collapse="+"), sep=""),sep="")
  modM <- lrm(formula=as.formula(form), dat=datM_dum_mode)   
  p2_pba <- rep(NA, length(p1_pba))
  p2_cba <- rep(NA, length(p1_cba))
  for(k in 1:length(p2_pba)){
    form <- paste(paste(paste("~ x",k,sep=""), " + x", sep=""), ncol(datM_dum)+k, sep="")
    p2_pba[k] <- deltamethod(as.formula(form), coef(modM), vcov(modM))
    form <- paste(paste(paste("~ x",k,sep=""), " + x", sep=""), 2*ncol(datM_dum)+k, sep="")
    p2_cba[k] <- deltamethod(as.formula(form), coef(modM), vcov(modM))
  } 
  modMainEff_PBA_List[[i]] <- cbind(p1_pba,p2_pba)
  modMainEff_CBA_List[[i]] <- cbind(p1_cba,p2_cba)  
  
  # -----------------------------
  # 2. Derive Propensity Scores
  # -----------------------------
  
  # A. Participation of students who were initially invited to group testing PBA, derive p1
  datG1 <- datIm[datIm$PBA_E %in% 1,]
  datG1 <- datG1[,-which(colnames(datG1) %in% c("ID_i", "PBA_E", "CBA_E", "CBA_R", "WBA_E", "WBA_R", "WBA_switch_E", "WBA_switch_R"))] 
  nam <- colnames(datG1)[-which(colnames(datG1) %in% c("ID_t", "PBA_R", "GPA_W4", "GPA_W6", "kidsHH_W1", "kidsHH_W3", "studydropout", "helpless", "selfconcept"))]
  form <- paste("PBA_R ~",paste(nam, collapse="+"), sep="")
  modG1 <-  lrm(formula=as.formula(form), dat=datG1) 
  pseudoR2[1,i] <- modG1$stats["R2"]
  modG1List[[i]] <- cbind(modG1$coefficients, sqrt(diag(modG1$var)))
  snam <- names(anova(modG1)[,"P"])
  seff <- snam[anova(modG1)[,"P"]<0.05] # select significant effects
  seff <- seff[-which(seff %in% "TOTAL")]
  form <- paste("PBA_R ~",paste(seff, collapse="+"), sep="")  
  modG1_m <- lrm(formula=as.formula(form), dat=datG1)  
  p1 <- predict(object=modG1_m, type="fitted.ind"); # hist(p1) 
  datG1$adj <- 1/p1
  TN_gr1 <- datG1[datG1$PBA_R %in% 1,]; # range(TN_gr1$adj)
  TN_gr1 <- TN_gr1[, !(colnames(TN_gr1) %in% "PBA_R")]
  datG1$adjGr4 <- 1/(1-p1) 

  # B. Participation of students who were initially invited to group testing CBA, derive p2  
  datG2 <- datIm[datIm$CBA_E %in% 1,]
  datG2 <- datG2[,-which(colnames(datG2) %in% c("ID_i", "CBA_E", "PBA_E", "PBA_R", "WBA_E", "WBA_R", "WBA_switch_E", "WBA_switch_R"))] 
  nam <- colnames(datG2)[-which(colnames(datG2) %in% c("ID_t", "CBA_R", "uniFunding", "GPA_W4", "GPA_W6", "studydropout", "helpless", "selfconcept"))] # no CBA at private institutions
  form <- paste("CBA_R ~",paste(nam, collapse="+"), sep="")
  modG2 <-  lrm(formula=as.formula(form), dat=datG2) 
  pseudoR2[2,i] <- modG2$stats["R2"]
  modG2List[[i]] <- cbind(modG2$coefficients, sqrt(diag(modG2$var)))
  snam <- names(anova(modG2)[,"P"])
  seff <- snam[anova(modG2)[,"P"]<0.05] # select significant effects
  seff <- seff[-which(seff %in% "TOTAL")]
  form <- paste("CBA_R ~",paste(seff, collapse="+"), sep="")  
  modG2_m <- lrm(formula=as.formula(form), dat=datG2)  
  p2 <- predict(object=modG2_m, type="fitted.ind"); # hist(p2) 
  datG2$adj <- 1/p2
  TN_gr2 <- datG2[datG2$CBA_R %in% 1,]; # range(TN_gr2$adj)
  TN_gr2 <- TN_gr2[, !(colnames(TN_gr2) %in% "CBA_R")]
  datG2$adjGr4 <- 1/(1-p2) 
  
  # C. Participation of students who were invited to online testing from the beginning, derive p^on_TN
  datG3 <- datIm[datIm$WBA_E %in% 1,]
  datG3 <- datG3[,-which(colnames(datG3) %in% c("ID_i", "PBA_E", "PBA_R", "CBA_E", "CBA_R", "WBA_E", "WBA_switch_E", "WBA_switch_R"))] 
  nam <- colnames(datG3)[-which(colnames(datG3) %in% c("ID_t", "WBA_R", "GPA_W4", "GPA_W6", "kidsHH_W1", "kidsHH_W3", "studydropout", "helpless", "selfconcept"))] 
  form <- paste("WBA_R ~",paste(nam, collapse="+"), sep="")
  modG3 <-  lrm(formula=as.formula(form), dat=datG3) 
  pseudoR2[3,i] <- modG3$stats["R2"]
  modG3List[[i]] <- cbind(modG3$coefficients, sqrt(diag(modG3$var)))
  snam <- names(anova(modG3)[,"P"])
  seff <- snam[anova(modG3)[,"P"]<0.05] # select significant effects
  seff <- seff[-which(seff %in% "TOTAL")]
  form <- paste("WBA_R ~",paste(seff, collapse="+"), sep="")  
  modG3_m <- lrm(formula=as.formula(form), dat=datG3)  
  p3 <- predict(object=modG3_m, type="fitted.ind"); # hist(p3) 
  datG3$adj <- 1/p3
  TN_gr3 <- datG3[datG3$WBA_R %in% 1,]; # range(TN_gr3$adj)
  TN_gr3 <- TN_gr3[, !(colnames(TN_gr3) %in% "WBA_R")]

  # D. Participation of students who did not participate in group testing (but were invited) and were asked to switch to online testing, derive p^sw_TN
  datG4 <- datIm[datIm$WBA_switch_E %in% 1,]
  datG4$groupModus <- ifelse(datG4$ID_t %in% datG1$ID_t, 1, 0) # PBA: 1, CBA: 0
  datG4 <- datG4[,-which(colnames(datG4) %in% c("ID_i", "PBA_E", "PBA_R", "CBA_E", "CBA_R", "WBA_E", "WBA_R", "WBA_switch_E"))] 
  nam <- colnames(datG4)[-which(colnames(datG4) %in% c("ID_t", "WBA_switch_R", "GPA_W4", "GPA_W6", "kidsHH_W1", "kidsHH_W3", "studydropout", "helpless", "selfconcept"))] 
  form <- paste("WBA_switch_R ~",paste(nam, collapse="+"), sep="")
  modG4 <-  lrm(formula=as.formula(form), dat=datG4) 
  pseudoR2[4,i] <- modG4$stats["R2"]
  modG4List[[i]] <- cbind(modG4$coefficients, sqrt(diag(modG4$var)))
  snam <- names(anova(modG4)[,"P"])
  seff <- snam[anova(modG4)[,"P"]<0.05] # select significant effects
  seff <- seff[-which(seff %in% "TOTAL")]
  form <- paste("WBA_switch_R ~",paste(seff, collapse="+"), sep="")  
  modG4_m <- lrm(formula=as.formula(form), dat=datG4)  
  datG4$p4 <- predict(object=modG4_m, type="fitted.ind"); # hist(p4) 
  datG4 <- merge(datG4, rbind(datG1[, c("ID_t", "adjGr4")], datG2[, c("ID_t", "adjGr4")]), by="ID_t", all.x=TRUE)
  datG4$adj <- datG4$adjGr4 * 1/datG4$p4
  TN_gr4 <- datG4[datG4$WBA_switch_R %in% 1, !(colnames(datG4) %in% c("groupModus", "p4", "adjGr4", "WBA_switch_R"))]
  
  # Bind it
  TN_gr1$group <- 1 # PBA
  TN_gr2$group <- 2 # CBA
  TN_gr3$group <- 3 # WBA
  TN_gr4$group <- 4 # WBA (switched)
  DAT <- rbind(TN_gr1, TN_gr2, TN_gr3, TN_gr4) # N=8254 participants, TODO SZ: REVISE N=8215

  # --------------------------------
  # 3. Differential Test Functioning 
  # --------------------------------
  
  # add competence data
  idMiss <- comp[comp$misSCI | comp$misSCIwle,]$ID_t # Remove those students that have either missing values in the science test or no wle for science
  comp <- comp[!(comp$ID_t %in% idMiss),] # Entries for N=8254 students. Fits. TODO SZ: REVISE N=8215
  DAT <- full_join(DAT, comp[, c("ID_t", "scs3_sc1u", citems$sci)], by = "ID_t") 
  rm(idMiss)
  
  # 0.5 scoring of PCMs
  pars <- mirt(data = DAT[, citems$sci],
               model = 1, pars = "values",
               itemtype = "Rasch")
  f <- pars$name == "a1" &
       pars$item %in% c("scs3033s_c", "scs3112s_c", "scs3131s_c", "scs3623s_c", "scs3021s_c",
                        "scs3132s_c", "scs3022s_c", "scs3133s_c", "scs3012s_c", "scs3643s_c",
                        "scs3642s_c", "scs3642s_c", "scs3061s_c", "scs3031s_c")
  pars$value[f] <- 0.5
  rm(f)
  
  # rescale adjustment factor to match sample size
  DAT$adj2 <- DAT$adj * nrow(Dat) / sum(DAT$adj)

  # fit models for science in each group
  modsSCI <- list(gr1 = mirt(data = DAT[DAT$group == 1, citems$sci],
                             model = 1,
                             itemtype = "Rasch", 
                             pars = pars,
                             survey.weights = DAT$adj2[DAT$group == 1],
                             SE = TRUE, SE.type = "Oakes", 
                             verbose = FALSE),
                  gr2 = mirt(data = DAT[DAT$group == 2, citems$sci],
                             model = 1,
                             itemtype = "Rasch", 
                             pars = pars,
                             survey.weights = DAT$adj2[DAT$group == 2],
                             SE = TRUE, SE.type = "Oakes",
                             verbose = FALSE),
                  gr3 = mirt(data = DAT[DAT$group == 3, citems$sci],
                             model = 1,
                             itemtype = "Rasch", 
                             pars = pars,
                             survey.weights = DAT$adj2[DAT$group == 3],
                             SE = TRUE, SE.type = "Oakes",
                             verbose = FALSE),
                  gr4 = mirt(data = DAT[DAT$group == 4, citems$sci],
                             model = 1,
                             itemtype = "Rasch", 
                             pars = pars,
                             survey.weights = DAT$adj2[DAT$group == 4],
                             SE = TRUE, SE.type = "Oakes",
                             verbose = FALSE))

  # extract item parameters from fitted models
  xsi <- lapply(modsSCI, function(y) {
                b <- t(sapply(coef(y), function(x) {
                        s <- rep(NA, 5)
                        names(s) <- paste0("d", 0:4)
                        v <- x["par", colnames(x) %in% paste0("d", c("", 0:4))]
                        if (length(v) == 1) v <- c(0, v)
                        s[seq_len(length(v))] <- v
                        s
                    })) * -1
                b <- b[seq_len(nrow(b) - 1), ]
                b
  })
  se <- lapply(modsSCI, function(y) {
                b <- t(sapply(coef(y, printSE = TRUE), function(x) {
                  s <- rep(0, 5)
                  names(s) <- paste0("d", 0:4)
                  v <- x["SE", colnames(x) %in% paste0("d", c("", 0:4))]
                  if (length(v) == 1) v <- c(0, v)
                  else v[1] <- 0
                  s[seq_len(length(v))] <- v
                  s
                }))
                b <- b[seq_len(nrow(b) - 1), ]
                b
  })
  paramSCI <- list(xsi = xsi, se = se)
  rm(xsi, se)
  
  # extract variances from fitted models
  varSCI <- t(sapply(modsSCI, function(y) {
              coef(y, printSE = TRUE)$GroupPars[, "COV_11"]
            })) 
  
  # calculate reliabilities
  relSCI <- t(sapply(modsSCI, function(y) {
              c(empirical_rxx(fscores(y, method = "WLE", full.scores.SE = TRUE)),
                marginal_rxx(y))
              }))
  
  # item fit
  Q <- matrix(1, length(citems$sci), 1)
  Q[apply(DAT[, citems$sci], 2, max, na.rm = TRUE) > 1, ] <- 0.5
  fitSCI <- list(gr1 = tam.fit(tam.mml(resp = DAT[DAT$group == 1, citems$sci], 
                                      irtmodel = "PCM", Q = Q, 
                                      verbose = FALSE), progress = FALSE)$itemfit[, c("Infit", "Infit_t")],
                 gr2 = tam.fit(tam.mml(resp = DAT[DAT$group == 2, citems$sci], 
                                       irtmodel = "PCM", Q = Q, 
                                       verbose = FALSE), progress = FALSE)$itemfit[, c("Infit", "Infit_t")],
                 gr3 = tam.fit(tam.mml(resp = DAT[DAT$group == 3, citems$sci], 
                                       irtmodel = "PCM", Q = Q, 
                                       verbose = FALSE), progress = FALSE)$itemfit[, c("Infit", "Infit_t")],
                 gr4 = tam.fit(tam.mml(resp = DAT[DAT$group == 4, citems$sci], 
                                       irtmodel = "PCM", Q = Q, 
                                       verbose = FALSE), progress = FALSE)$itemfit[, c("Infit", "Infit_t")])
  rm(Q)
  
  # estimate DTF (absolute scale, expected item test score)
  dtfSCI <- IRT.dtf(xsi.list = paramSCI$xsi, se.list = paramSCI$se, n = 21, R = 100) #   to increase precision of SE

  # estimate DIF (logit scale)
  # note: works fine, but we don't report the results, so skip it and speed up the estimation
  #group <- unclass(DAT$group)
  #difSCI <- IRT.dif(resp = DAT[, citems$sci], group = group, wgt = DAT$adj)
  #difSCI <- list(grp12= difSCI[, c("b12", "se12")],
  #               grp13= difSCI[, c("b13", "se13")],
  #               grp14= difSCI[, c("b14", "se14")],
  #               grp23= difSCI[, c("b23", "se23")],
  #               grp24= difSCI[, c("b24", "se24")],
  #               grp34= difSCI[, c("b34", "se34")])
  
  # expected test score
  # note: a little awkward, but the mirt function does not respect the half scoring
  theta <- seq(-6, 6, by = .1)
  w <- pars$value[pars$name == "a1"]
  etsSCI <- 0
  for (j in seq_len(length(citems$sci))) {
    extr1 <- extract.item(modsSCI$gr1, j)
    extr2 <- extract.item(modsSCI$gr2, j)
    extr3 <- extract.item(modsSCI$gr3, j)
    extr4 <- extract.item(modsSCI$gr4, j)
    iscores <- cbind(expected.item(extr1, theta) * w[j],
                     expected.item(extr2, theta) * w[j],
                     expected.item(extr3, theta) * w[j],
                     expected.item(extr4, theta) * w[j])
    etsSCI <- etsSCI + iscores
  }
  rm(theta, w, extr1, extr2, extr3, extr4, iscores, j)

  # description of output:
  #
  # - dtfSCI:          results for differential test functioning (DTF) 
  #   + sDTF          signed DTF (= differences in expected test scores)
  #   + uDTF          unsigned DTF (= absolute differences in expected test scores)
  #   + sDTF21        signed DTF at different competence levels for groups 1 and 2
  #   + uDTF21        unsigned DTF at different competence levels for groups 1 and 2
  #   + sDTF31 etc.   same as for sDTF21 but for different groups
  #   + uDTF31 etc.   same as for uDTF21 but for different groups
  #
  # - difSCI:          results for differential item functioning (DIF) 
  #                      contains differences in item parameters between two groups (b) and
  #                      respective standard errors (se)
  #                   main effect refers to the latent difference between groups
  #
  # - varSCI: latent variances in each group
  #
  # - relSCI: reliabilities in each group
  #
  # - etsSCI: expected test scores for different proficiencies in each group
  #
  # - fitSCI: infit for all items in each group
  #
  # note: b = parameter, se = standard error
  #       for reliabilities no SEs are given (-> no significance tests)  

  # collect results
   psychList[[i]] <- list(dtfSCI = dtfSCI, # differential test functioning for science
                          difSCI = NULL,#difSCI, # differential item functioning for science
                          varSCI = varSCI, # (latent) variances for science
                          relSCI = relSCI, # reliabilities for science
                          etsSCI = etsSCI, # expected test score for science
                          fitSCI = fitSCI) # item fit for science
   rm(dtfSCI, varSCI, relSCI, etsSCI, fitSCI)
                          
  # ------------------------------------------------
  # 4. Predicting outcomes in Wave 6 for the distinct modes
  # ------------------------------------------------

   # estimate science WLE
   mod <- mirt(data = DAT[, citems$sci],
               model = 1,
               itemtype = "Rasch", 
               pars = pars,
               survey.weights = DAT$adj2,
               SE = FALSE, 
               verbose = FALSE)
   DAT$sciwle <- c(fscores(mod, method = "WLE", full.scores.SE = FALSE))
   rm(mod)

   # Use science wle for prediction
   f <- DAT$naturalScience == 1 # select natural science students
   DAT$female_ec <- DAT$female - 0.5 # effect coding
   DAT$teacherEdu_ec <- ifelse(DAT$teacherEdu == 0, -0.5, 0.5)
   modPr10 <- lm(scale(GPA_W6) ~ as.factor(group) + scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   modPr11 <- lm(scale(GPA_W6) ~ as.factor(group) * scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   modPr20 <- lm(scale(helpless) ~ as.factor(group) + scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   modPr21 <- lm(scale(helpless) ~ as.factor(group) * scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   modPr30 <- lm(scale(selfconcept) ~ as.factor(group) + scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   modPr31 <- lm(scale(selfconcept) ~ as.factor(group) * scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   modPr40 <- lm(scale(studydropout) ~ as.factor(group) + scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   modPr41 <- lm(scale(studydropout) ~ as.factor(group) * scale(sciwle) + female_ec + teacherEdu_ec, 
                 data = DAT[f, ], weights = adj)
   predList_sc$pred10[[i]] <- coef(summary(modPr10))[,1:2] 
   predList_sc$pred11[[i]] <- coef(summary(modPr11))[,1:2] 
   predList_sc$pred20[[i]] <- coef(summary(modPr20))[,1:2] 
   predList_sc$pred21[[i]] <- coef(summary(modPr21))[,1:2] 
   predList_sc$pred30[[i]] <- coef(summary(modPr30))[,1:2] 
   predList_sc$pred31[[i]] <- coef(summary(modPr31))[,1:2] 
   predList_sc$pred40[[i]] <- coef(summary(modPr40))[,1:2] 
   predList_sc$pred41[[i]] <- coef(summary(modPr41))[,1:2] 
   r2_sc$pred10[i] <- summary(modPr10)$r.squared
   r2_sc$pred11[i] <- summary(modPr11)$r.squared
   r2_sc$pred20[i] <- summary(modPr20)$r.squared
   r2_sc$pred21[i] <- summary(modPr21)$r.squared
   r2_sc$pred30[i] <- summary(modPr30)$r.squared
   r2_sc$pred31[i] <- summary(modPr31)$r.squared
   r2_sc$pred40[i] <- summary(modPr40)$r.squared
   r2_sc$pred41[i] <- summary(modPr41)$r.squared
   rm(f, pars, modPr10, modPr11, modPr20, modPr21, modPr30, modPr31, modPr40, modPr41)
}
rm(i, DAT, datG1, datG2, datG3, datG4, datIm, datM, datM_dum, datM_dum_mode, datM_temp,
   modG1, modG1_m, modG2, modG2_m, modG3, modG3_m, modG4, modG4_m, modM, modsSCI, namXX,
   paramSCI, TN_gr1, TN_gr2, TN_gr3, TN_gr4, form, k, nam, namK, namM, nn, p1, p1_cba,
   p1_pba, p2, p2_cba, p2_pba, p3, seff, snam)



# --------------------------------------------------------------------------
# Pool results for self-selection
# --------------------------------------------------------------------------

# function to pool individual results
getPooled <- function(RES){ # RES is matrix with m rows and 2 cols: [est. theta, est. std.Err for theta]   
  estC <- mean(RES[,1])
  uC <- mean(RES[,2]^2)
  bC <- 1/(imp1$m-1)*sum((RES[,1]-estC)^2)
  varC <- uC + (1+1/imp1$m)*bC
  return(cbind(estC, sqrt(varC)))
}
alpha <- 0.05

# Pool self-selection model
M <- do.call(cbind,modMList)
res_modM <- matrix(NA, ncol=2, nrow=nrow(M))
rownames(res_modM) <- rownames(modMList[[1]])
colnames(res_modM) <- c("beta","stdBeta")
for(cc in 1:nrow(modMList[[1]])){
  theta <- M[cc,seq(from=1, to=ncol(M), by=2)]
  stE <- M[cc,seq(from=2, to=ncol(M), by=2)]
  res_modM[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modM[,1] - qt(1-alpha/2, nrow(Dat))*res_modM[,2]
ci_up  <- res_modM[,1] + qt(1-alpha/2, nrow(Dat))*res_modM[,2]
res_modM <- round(cbind(res_modM, ci_low, ci_up),2)
rm(M, cc, theta, stE, ci_low, ci_up)

# Pool self-selection model: WBA
M <- do.call(cbind,modMainEff_WBA_List)
res_modM <- matrix(NA, ncol=2, nrow=nrow(M))
rownames(res_modM) <- rownames(modMainEff_WBA_List[[1]])
colnames(res_modM) <- c("beta","stdBeta")
for(cc in 1:nrow(modMainEff_WBA_List[[1]])){
  theta <- M[cc,seq(from=1, to=ncol(M), by=2)]
  stE <- M[cc,seq(from=2, to=ncol(M), by=2)]
  res_modM[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modM[,1] - qt(1-alpha/2, nrow(Dat))*res_modM[,2]
ci_up  <- res_modM[,1] + qt(1-alpha/2, nrow(Dat))*res_modM[,2]
res_modM_WBA <- round(cbind(res_modM, ci_low, ci_up),2)
rm(M, cc, theta, stE, ci_low, ci_up)

# Pool self-selection model: PBA
M <- do.call(cbind,modMainEff_PBA_List)
res_modM <- matrix(NA, ncol=2, nrow=nrow(M))
rownames(res_modM) <- rownames(modMainEff_PBA_List[[1]])
colnames(res_modM) <- c("beta","stdBeta")
for(cc in 1:nrow(modMainEff_PBA_List[[1]])){
  theta <- M[cc,seq(from=1, to=ncol(M), by=2)]
  stE <- M[cc,seq(from=2, to=ncol(M), by=2)]
  res_modM[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modM[,1] - qt(1-alpha/2, nrow(Dat))*res_modM[,2]
ci_up  <- res_modM[,1] + qt(1-alpha/2, nrow(Dat))*res_modM[,2]
res_modM_PBA <- round(cbind(res_modM, ci_low, ci_up),2)
rm(M, cc, theta, stE, ci_low, ci_up)

# Pool self-selection model: CBA
M <- do.call(cbind,modMainEff_CBA_List)
res_modM <- matrix(NA, ncol=2, nrow=nrow(M))
rownames(res_modM) <- rownames(modMainEff_CBA_List[[1]])
colnames(res_modM) <- c("beta","stdBeta")
for(cc in 1:nrow(modMainEff_CBA_List[[1]])){
  theta <- M[cc,seq(from=1, to=ncol(M), by=2)]
  stE <- M[cc,seq(from=2, to=ncol(M), by=2)]
  res_modM[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modM[,1] - qt(1-alpha/2, nrow(Dat))*res_modM[,2]
ci_up  <- res_modM[,1] + qt(1-alpha/2, nrow(Dat))*res_modM[,2]
res_modM_CBA <- round(cbind(res_modM, ci_low, ci_up),2)
rm(M, cc, theta, stE, ci_low, ci_up)

# Pool non-response model for group testing: PBA
A <- do.call(cbind,modG1List)
res_modA <- matrix(NA, ncol=2, nrow=nrow(A))
rownames(res_modA) <- rownames(modG1List[[1]])
colnames(res_modA) <- c("beta","stdBeta")
for(cc in 1:nrow(modG1List[[1]])){
  theta <- A[cc,seq(from=1, to=ncol(A), by=2)]
  stE <- A[cc,seq(from=2, to=ncol(A), by=2)]
  res_modA[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modA[,1] - qt(1-alpha/2, sum(Dat$PBA_E %in% 1))*res_modA[,2]
ci_up  <- res_modA[,1] + qt(1-alpha/2, sum(Dat$PBA_E %in% 1))*res_modA[,2]
res_modA <- cbind(res_modA, ci_low, ci_up)
rm(A, cc, theta, stE, ci_low, ci_up)

# Pool non-response model for group testing: CBA
B <- do.call(cbind,modG2List)
res_modB <- matrix(NA, ncol=2, nrow=nrow(B))
rownames(res_modB) <- rownames(modG2List[[1]])
colnames(res_modB) <- c("beta","stdBeta")
for(cc in 1:nrow(modG2List[[1]])){
  theta <- B[cc,seq(from=1, to=ncol(B), by=2)]
  stE <- B[cc,seq(from=2, to=ncol(B), by=2)]
  res_modB[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modB[,1] - qt(1-alpha/2, sum(Dat$CBA_E %in% 1))*res_modB[,2]
ci_up  <- res_modB[,1] + qt(1-alpha/2, sum(Dat$CBA_E %in% 1))*res_modB[,2]
res_modB <- cbind(res_modB, ci_low, ci_up)
rm(B, cc, theta, stE, ci_low, ci_up)

# Pool non-response model for WBA
C <- do.call(cbind,modG3List)
res_modC <- matrix(NA, ncol=2, nrow=nrow(C))
rownames(res_modC) <- rownames(modG3List[[1]])
colnames(res_modC) <- c("beta","stdBeta")
for(cc in 1:nrow(modG3List[[1]])){
  theta <- C[cc,seq(from=1, to=ncol(C), by=2)]
  stE <- C[cc,seq(from=2, to=ncol(C), by=2)]
  res_modC[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modC[,1] - qt(1-alpha/2, sum(Dat$WBA_E %in% 1))*res_modC[,2]
ci_up  <- res_modC[,1] + qt(1-alpha/2, sum(Dat$WBA_E %in% 1))*res_modC[,2]
res_modC <- cbind(res_modC, ci_low, ci_up)
rm(C, cc, theta, stE, ci_low, ci_up)

# Pool non-response model for WBA (switcher)
D <- do.call(cbind,modG4List)
res_modD <- matrix(NA, ncol=2, nrow=nrow(D))
rownames(res_modD) <- rownames(modG4List[[1]])
colnames(res_modD) <- c("beta","stdBeta")
for(cc in 1:nrow(modG4List[[1]])){
  theta <- D[cc,seq(from=1, to=ncol(D), by=2)]
  stE <- D[cc,seq(from=2, to=ncol(D), by=2)]
  res_modD[cc,] <- getPooled(cbind(theta, stE))
}
ci_low <- res_modD[,1] - qt(1-alpha/2, sum(Dat$WBA_switch_E %in% 1))*res_modD[,2]
ci_up  <- res_modD[,1] + qt(1-alpha/2, sum(Dat$WBA_switch_E %in% 1))*res_modD[,2]
res_modD <- cbind(res_modD, ci_low, ci_up)
rm(D, cc, theta, stE, ci_low, ci_up)



# ------------------------------------------------------------------------------
# Pool results of DIF / DTF -> insert everything into object "psychRes" 
# ------------------------------------------------------------------------------

# result object
psychRes <- vector(length=length(names(psychList[[1]])), mode="list")
names(psychRes) <- names(psychList[[1]])
psychRes[[1]] <- vector(length=length(names(psychList[[1]]$dtfSC)), mode="list")
names(psychRes[[1]]) <- names(psychList[[1]]$dtfSC)
psychRes[[2]] <- vector(length=length(names(psychList[[1]]$difSC)), mode="list")
names(psychRes[[2]]) <- names(psychList[[1]]$difSC)
psychRes[[3]] <- vector(length=length(names(psychList[[1]]$varSCI)), mode="list")
names(psychRes[[3]]) <- names(psychList[[1]]$varSCI)
psychRes[[4]] <- vector(length=length(names(psychList[[1]]$relSCI)), mode="list")
names(psychRes[[4]]) <- names(psychList[[1]]$relSCI)
psychRes[[5]] <- vector(length=length(names(psychList[[1]]$etsSCI)), mode="list")
names(psychRes[[5]]) <- names(psychList[[1]]$etsSCI)
psychRes[[6]] <- vector(length=length(names(psychList[[1]]$fitSCI)), mode="list")
names(psychRes[[6]]) <- names(psychList[[1]]$fitSCI)
alpha <- 0.05
lg1 <- sum(Dat$PBA_E %in% 1) # group sizes
lg2 <- sum(Dat$CBA_E %in% 1)
lg3 <- sum(Dat$WBA_E %in% 1)
lg4 <- sum(Dat$WBA_switch_E %in% 1)
nn <- c(lg1+lg2+lg3+lg4,lg1+lg2+lg3+lg4,
        lg1+lg2, lg3+lg1, lg4+lg1, lg3+lg2, lg4+lg2, lg4+lg3, 
        lg1+lg2, lg3+lg1, lg4+lg1, lg3+lg2, lg4+lg2, lg4+lg3)

# pool dtfSC and difSC
for(i in 1:2){
  for(j in seq_len(length(psychList[[1]][[i]]))){ #sDTF, uDTF, sDTF21, ..., sDTF43, uDTF21, ... , uDTF43 or sDIF, uDIF, sDIF21, ..., sDIF43, uDIF21, ... , uDIF43
    #cat("i: ",i," - j:",j,"\n")
    resTemp <- NULL
    for(m in 1:imp1$m){ # collect results from all imputed data sets
      if (m == 1) resTemp <- psychList[[m]][[i]][[j]]
      else resTemp <- cbind(resTemp, psychList[[m]][[i]][[j]])
    }
    res_temp_mat <- matrix(NA, ncol=2, nrow=nrow(resTemp))
    rownames(res_temp_mat) <- rownames(psychList[[m]][[i]][[j]])
    colnames(res_temp_mat) <- c("est","stdBeta")
    for(cc in 1:nrow(res_temp_mat)){ # pool it
      theta <-  resTemp[cc,seq(from=1, to=ncol(resTemp), by=2)]
      stE <-  resTemp[cc,seq(from=2, to=ncol(resTemp), by=2)]
      res_temp_mat[cc,] <- getPooled(cbind(theta, stE))
    }
    # compute confidence intervals
    ci_low <- res_temp_mat[,1] - qt(1-alpha/2, nn[j])*res_temp_mat[,2] # get the sample size of the considered samples right (group sizes)
    ci_up  <- res_temp_mat[,1] + qt(1-alpha/2, nn[j])*res_temp_mat[,2]
    psychRes[[i]][[j]] <- cbind(res_temp_mat[,1], ci_low, ci_up)
  }
}
psychRes$dtfSCI$uDTFp <- psychRes$dtfSCI$uDTF / 36 * 100
rm(i, j, ci_low, ci_up, theta, stE, res_temp_mat, cc, resTemp, m, nn)

# pool variances
resTemp <- NULL
for(m in 1:imp1$m){ # collect results from all imputed data sets
  if (m == 1) resTemp <- psychList[[m]][[3]]
  else resTemp <- cbind(resTemp, psychList[[m]][[3]])
}
res_temp_mat <- matrix(NA, ncol=2, nrow=nrow(resTemp))
rownames(res_temp_mat) <- rownames(psychList[[m]][[3]])
colnames(res_temp_mat) <- c("est","stdBeta")
for(cc in 1:nrow(res_temp_mat)){ # pool it
  theta <-  resTemp[cc,seq(from=1, to=ncol(resTemp), by=2)]
  stE <-  resTemp[cc,seq(from=2, to=ncol(resTemp), by=2)]
  res_temp_mat[cc,] <- getPooled(cbind(theta, stE))
}
# compute confidence intervals
ci_low <- res_temp_mat[,1] - qt(1-alpha/2, c(lg1, lg2, lg3, lg4))*res_temp_mat[,2] # get the sample size of the considered samples right (group sizes)
ci_up  <- res_temp_mat[,1] + qt(1-alpha/2, c(lg1, lg2, lg3, lg4))*res_temp_mat[,2]
psychRes[[3]] <- cbind(res_temp_mat[,1], ci_low, ci_up)
rm(ci_low, ci_up, cc, theta, stE, res_temp_mat, resTemp, m)

# add ranges for relSC to "psychRes" object
allV <- NULL
for(m in 1:imp1$m){  # collect results from all imputed data sets
  allV <- rbind(allV, c(psychList[[m]][[4]][, 1], psychList[[m]][[4]][, 2]))
}
psychRes[[4]] <- apply(allV,2,range)
psychRes[[4]] <- rbind(psychRes[[4]], apply(allV, 2, mean))
rownames(psychRes[[4]]) <- c("min", "max", "mean")
colnames(psychRes[[4]]) <- c(paste0("ERgr", 1:4), paste0("MRgr", 1:4))
rm(allV, m)

# pool expected test scores
allV <- 0
for(m in 1:imp1$m){  # collect results from all imputed data sets
  allV <- allV + psychList[[m]]$etsSCI
}
psychRes[[5]] <- allV / imp1$m
rm(allV, m)

# pool infit
psychRes[[6]] <- list()
for (i in seq_len(length(psychList[[1]]$fitSCI))) {
  psychRes[[6]][[i]] <- 0
  for(m in 1:imp1$m){  # collect results from all imputed data sets
    psychRes[[6]][[i]] <- psychRes[[6]][[i]] + psychList[[m]]$fitSCI[[i]]
  }
  psychRes[[6]][[i]] <- psychRes[[6]][[i]] / imp1$m
}
rm(i, m, lg1, lg2, lg3, lg4, alpha)



# ------------------------------------------------------------------------------
# Pool results of prediction error
# ------------------------------------------------------------------------------

# parameter estimates
alpha <- .05
res_modP <- list()
for (i in names(predList_sc)) {
  P_sc <- do.call(cbind,predList_sc[[i]])
  res_modP[[i]] <- matrix(NA, ncol=2, nrow=nrow(P_sc))
  rownames(res_modP[[i]]) <- rownames(predList_sc[[i]][[1]])
  colnames(res_modP[[i]]) <- c("beta","stdBeta")
  for(cc in 1:nrow(predList_sc[[i]][[1]])){
    theta <- P_sc[cc,seq(from=1, to=ncol(P_sc), by=2)]
    stE <- P_sc[cc,seq(from=2, to=ncol(P_sc), by=2)]
    res_modP[[i]][cc,] <- getPooled(cbind(theta, stE))
  }
  ci_low <- res_modP[[i]][,1] - qt(1-alpha/2, sum(Dat$naturalScience %in% 1))*res_modP[[i]][,2]
  ci_up  <- res_modP[[i]][,1] + qt(1-alpha/2, sum(Dat$naturalScience %in% 1))*res_modP[[i]][,2]
  res_modP[[i]] <- cbind(res_modP[[i]], ci_low, ci_up)
}
rm(ci_low, ci_up, cc, theta, stE, i, P_sc, alpha)



# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

obj <- c("modMList", "modMainEff_PBA_List", "modMainEff_CBA_List", 
         "modMainEff_WBA_List", "modG1List", "modG2List",
         "modG3List", "modG4List", "pseudoMR2",
         "psychList",
         "predList_sc", "r2_sc",
         "res_modM","res_modM_PBA", "res_modM_CBA", "res_modM_WBA",
         "res_modA", "res_modB", "res_modC", "res_modD",
         "psychRes", "res_modP")
save(list = obj, file = paste0("Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/results/_final-", version, ".RData"))
rm(obj)



# ------------------------------------------------------------------------------
# Plot results of pooled self-selection analysis
# ------------------------------------------------------------------------------

# load results
load(paste0("Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/results/_final-", version, ".RData"))

# only interaction effects
# res_modMa_PBA <- rbind(res_modM[3,], res_modM[seq(from=31, to=nrow(res_modM), by=2),])
# rownames(res_modMa_PBA)
# namM <- c("PBAorCBA_vs_WBA", "NaturalScience_vs_Non", "TeacherEdu_vs_Non", "EnjoyStud", "ProbGrad",
#               "Big5:extravers", "Big5:agreeabl", "Big5:conscien", "Big5:neurot", "Big5:open", "SelfEsteem",
#               "2PerInHH_vs_1Per", "3+PerInHH_vs_1Per", "Abitur_vs_NoAbi", "AbiturNonGerm_vs_NoAbi","KidsinHH_vs_Non", 
#               "AverageGrade2ndSem", "BirthYear1989and1990_vs_older1989", "BirthYearYounger1990_vs_older1989",
#               "Female_vs_Male", "BirthInGerm_vs_NonGerm", "MotherTongueGerm_vs_NonGerm", "EastGerm_vs_WestGerm",
#               "EduYearsMother", "EduYearsFather", "UniApplScie_vs_GeneralUniv", "OtherType_vs_GeneralUniv", "MathCompScoreWave1")
# modelPBAFrame <- data.frame(Predictor = namM,
#                           B = res_modMa_PBA[,1],
#                           CI_low= res_modMa_PBA[,3],
#                           CI_up = res_modMa_PBA[,4],
#                           Mode="PBA_vs_WBA")
# modelPBAFrame$Predictor <- as.factor(modelPBAFrame$Predictor)
# res_modMa_CBA <- rbind(res_modM[4,], res_modM[seq(from=32, to=nrow(res_modM), by=2),])
# rownames(res_modMa_CBA)
# modelCBAFrame <- data.frame(Predictor = namM,
#                             B = res_modMa_CBA[,1],
#                             CI_low= res_modMa_CBA[,3],
#                             CI_up = res_modMa_CBA[,4],
#                             Mode="CBA_vs_WBA")
# modelCBAFrame$Predictor <- as.factor(modelCBAFrame$Predictor)
# allModelFrame <- data.frame(rbind(modelPBAFrame, modelCBAFrame))
# 
# levNam <- rev(c("BirthYear1989and1990_vs_older1989", "BirthYearYounger1990_vs_older1989",
#                 "Female_vs_Male", "BirthInGerm_vs_NonGerm", "MotherTongueGerm_vs_NonGerm",
#                 "KidsinHH_vs_Non", "EduYearsMother", "EduYearsFather",
#                 "2PerInHH_vs_1Per", "3+PerInHH_vs_1Per",
#                 "Abitur_vs_NoAbi", "AbiturNonGerm_vs_NoAbi",
#                 "TeacherEdu_vs_Non", "NaturalScience_vs_Non", "EastGerm_vs_WestGerm",
#                 "UniApplScie_vs_GeneralUniv", "OtherType_vs_GeneralUniv", "MathCompScoreWave1",
#                 "EnjoyStud", "ProbGrad", "AverageGrade2ndSem", 
#                 "Big5:extravers", "Big5:agreeabl", "Big5:conscien", "Big5:neurot", "Big5:open", "SelfEsteem",
#                 "PBAorCBA_vs_WBA"))
# allModelFrame$Predictor <- factor(allModelFrame$Predictor, levels=c(levNam))
# 
# zp <- ggplot(allModelFrame, aes(colour = Mode))
# zp <- zp + geom_hline(yintercept = 0, colour = gray(1/2), lty=2)
# zp <- zp + geom_linerange(aes(x=Predictor, ymin= CI_low, ymax=CI_up),
#                             lwd=1, position=position_dodge(width = 1/2))
# zp <- zp + geom_point(aes(x=Predictor, y=B), 
#                         shape=21, fill="WHITE",
#                         position=position_dodge(width = 1/2))
# zp <- zp + coord_flip() + theme_bw()
# print(zp)

# mode-specific main effects
rownames(res_modM_PBA)
namM <- c("NaturalScience_vs_Non", "TeacherEdu_vs_Non", "EnjoyStud", "ProbGrad",
          "Big5:extravers", "Big5:agreeabl", "Big5:conscien", "Big5:neurot", "Big5:open", "SelfEsteem",
          "2PerInHH_vs_1Per", "3+PerInHH_vs_1Per", "Abitur_vs_NoAbi", "AbiturNonGerm_vs_NoAbi","KidsinHH_vs_Non", 
          "AverageGrade2ndSem", "BirthYear1989and1990_vs_older1989", "BirthYearYounger1990_vs_older1989",
          "Female_vs_Male", "BirthInGerm_vs_NonGerm", "MotherTongueGerm_vs_NonGerm", "EastGerm_vs_WestGerm",
          "EduYearsMother", "EduYearsFather", "UniApplScie_vs_GeneralUniv", "OtherType_vs_GeneralUniv", "MathCompScoreWave1")
modelPBAFrame <- data.frame(Predictor = namM,
                            B = res_modM_PBA[,1],
                            CI_low= res_modM_PBA[,3],
                            CI_up = res_modM_PBA[,4],
                            Mode="PBA")
modelPBAFrame$Predictor <- as.factor(modelPBAFrame$Predictor)

rownames(res_modM_CBA)
modelCBAFrame <- data.frame(Predictor = namM,
                            B = res_modM_CBA[,1],
                            CI_low= res_modM_CBA[,3],
                            CI_up = res_modM_CBA[,4],
                            Mode="CBA")
modelCBAFrame$Predictor <- as.factor(modelCBAFrame$Predictor)

rownames(res_modM_WBA)
modelWBAFrame <- data.frame(Predictor = namM,
                            B = res_modM_WBA[,1],
                            CI_low= res_modM_WBA[,3],
                            CI_up = res_modM_WBA[,4],
                            Mode="WBA")
modelWBAFrame$Predictor <- as.factor(modelWBAFrame$Predictor)

allModelFrame <- data.frame(rbind(modelPBAFrame, modelCBAFrame, modelWBAFrame))

levNam <- rev(c("BirthYear1989and1990_vs_older1989", "BirthYearYounger1990_vs_older1989",
                "Female_vs_Male", "BirthInGerm_vs_NonGerm", "MotherTongueGerm_vs_NonGerm",
                "KidsinHH_vs_Non", "EduYearsMother", "EduYearsFather",
                "2PerInHH_vs_1Per", "3+PerInHH_vs_1Per",
                "Abitur_vs_NoAbi", "AbiturNonGerm_vs_NoAbi",
                "TeacherEdu_vs_Non", "NaturalScience_vs_Non", "EastGerm_vs_WestGerm",
                "UniApplScie_vs_GeneralUniv", "OtherType_vs_GeneralUniv", "MathCompScoreWave1",
                "EnjoyStud", "ProbGrad", "AverageGrade2ndSem", 
                "Big5:extravers", "Big5:agreeabl", "Big5:conscien", "Big5:neurot", "Big5:open", "SelfEsteem"))

allModelFrame$Predictor <- factor(allModelFrame$Predictor, levels=c(levNam))
levels(allModelFrame$Predictor) <- list("Year of birth: before 1989 vs. 1989/90" = "BirthYear1989and1990_vs_older1989",
                                        "Year of birth: before 1989 vs. after 1990" = "BirthYearYounger1990_vs_older1989",
                                        "Gender: male vs. female" = "Female_vs_Male",
                                        "Country of birth: other vs. Germany" = "BirthInGerm_vs_NonGerm",
                                        "Mother tongue: non-German vs. German" = "MotherTongueGerm_vs_NonGerm",
                                        "Children in houshold: no vs. yes" = "KidsinHH_vs_Non",
                                        "Number of years in education of mother" = "EduYearsMother", 
                                        "Number of years in education of father" = "EduYearsFather",
                                        "Houshold size: 1 vs. 2 persons" = "2PerInHH_vs_1Per",
                                        "Houshold size: 1 vs. 3 or more persons" = "3+PerInHH_vs_1Per",
                                        "University admission certificate: none vs. German" = "Abitur_vs_NoAbi",
                                        "University admission certificate: none vs. non-German" = "AbiturNonGerm_vs_NoAbi",
                                        "Field of study in teacher education: no vs. yes" = "TeacherEdu_vs_Non",
                                        "Field of study in natural sciences: no vs. yes" = "NaturalScience_vs_Non",
                                        "Region of Germany: West vs. East" = "EastGerm_vs_WestGerm",
                                        "University type: general vs. applied sciences" = "UniApplScie_vs_GeneralUniv",
                                        "University type: general vs. other" = "OtherType_vs_GeneralUniv",
                                        "Mathematical literacy score in wave 1" = "MathCompScoreWave1",
                                        "Study enjoyment" = "EnjoyStud",
                                        "Perceived probability of graduation" = "ProbGrad",
                                        "Average grade after the 2nd semester" = "AverageGrade2ndSem",
                                        "Big Five: Extraversion" = "Big5:extravers", 
                                        "Big Five: Agreeableness" = "Big5:agreeabl", 
                                        "Big Five: Conscientiousness" = "Big5:conscien", 
                                        "Big Five: Neuroticisim" = "Big5:neurot", 
                                        "Big Five: Openness" = "Big5:open",
                                        "Self-esteem" = "SelfEsteem")

zp <- ggplot(allModelFrame, aes(colour = Mode))
zp <- zp + geom_hline(yintercept = 0, colour = gray(1/2), lty=2)
zp <- zp + geom_linerange(aes(x=Predictor, ymin= CI_low, ymax=CI_up),
                          lwd=1, position=position_dodge(width = 1/2))
zp <- zp + geom_point(aes(x=Predictor, y=B), 
                      shape=21, fill="WHITE",
                      position=position_dodge(width = 1/2))
zp <- zp + coord_flip() + theme_bw() + xlab("")
print(zp)
ggsave("../plots/selfselectionbymode.tif", device = "tiff", dpi = 300,
       width = 20, height = 20, units = "cm")



# ------------------------------------------------------------------------------
# Plot results of pooled DTF analyses
# ------------------------------------------------------------------------------

# load results
load(paste0("Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/results/_final-", version, ".RData"))

# plot expected test score
tiff("../plots/expectedtestscore.tif", width = 480*3, height = 900, res = 300,
     pointsize = 10)
plot(seq(-6, 6, by = .1), psychRes[[5]][, 1],
     main = "Test characteristic curves",
     xlab = "Latent proficiency",
     ylab = "Expected test score",
     type = "l", lwd = 2, col = "black")
points(seq(-6, 6, by = .1), psychRes[[5]][, 2], type = "l", 
       lty = 5, lwd = 1, col = "black")
points(seq(-6, 6, by = .1), psychRes[[5]][, 3], type = "l", 
       lty = 1, lwd = 2, col = "gray42")
points(seq(-6, 6, by = .1), psychRes[[5]][, 4], type = "l", 
       lty = 3, lwd = 1, col = "gray42")
legend(1, 25, c("PBT", "CBT", "WBT", "WBT-switch"),
       lty = c(1, 5, 1, 3), col = c("black", "black", "gray42", "gray42"),
       lwd = c(2, 1, 2, 1), bty = "n")
dev.off()

# plot DTF
tiff("../plots/DTFa.tif", width = 480*3, height = 910, res = 300,
     pointsize = 8)
plot(rownames(psychRes$dtfSCI$sDTF31),
     psychRes$dtfSCI$sDTF31[, 1] * -1,
     main = "Average difference in test scores between PBT and WBT",
     xlab = "Latent proficiency",
     ylab = "sDTF",
     type = "l", lwd = 2,
     ylim = c(-1, 3),
     axes = FALSE)
polygon(x = c(rownames(psychRes$dtfSCI$sDTF31), rev(rownames(psychRes$dtfSCI$sDTF31))),
        y = c(psychRes$dtfSCI$sDTF31[, 2], rev(psychRes$dtfSCI$sDTF31[, 3])) * -1, 
        col = adjustcolor("black", alpha.f = .1), border = NA)
points(rownames(psychRes$dtfSCI$sDTF31),
       rep(0, nrow(psychRes$dtfSCI$sDTF31)), 
       type = "l", lty = 3,
       col = "gray42")
axis(1, at = seq(-8, 8, 2))
axis(2, at = seq(-1, 3, 1))
dev.off()

tiff("../plots/DTFb.tif", width = 480*3, height = 910, res = 300,
     pointsize = 8)
plot(rownames(psychRes$dtfSCI$sDTF21),
     psychRes$dtfSCI$sDTF21[, 1] * -1,
     main = "Average difference in test scores between PBT and CBT",
     xlab = "Latent proficiency",
     ylab = "sDTF",
     type = "l", lwd = 2,
     ylim = c(-1, 3),
     axes = FALSE)
polygon(x = c(rownames(psychRes$dtfSCI$sDTF21), rev(rownames(psychRes$dtfSCI$sDTF21))),
        y = c(psychRes$dtfSCI$sDTF21[, 2], rev(psychRes$dtfSCI$sDTF21[, 3])) * -1, 
        col = adjustcolor("black", alpha.f = .1), border = NA)
points(rownames(psychRes$dtfSCI$sDTF21),
       rep(0, nrow(psychRes$dtfSCI$sDTF21)), 
       type = "l", lty = 3,
       col = "gray42")
axis(1, at = seq(-8, 8, 2))
axis(2, at = seq(-1, 3, 1))
dev.off()

tiff("../plots/DTFc.tif", width = 480*3, height = 910, res = 300,
     pointsize = 8)
plot(rownames(psychRes$dtfSCI$sDTF32),
     psychRes$dtfSCI$sDTF32[, 1] * -1,
     main = "Average difference in test scores between CBT and WBT",
     xlab = "Latent proficiency",
     ylab = "sDTF",
     type = "l", lwd = 2,
     ylim = c(-1, 3),
     axes = FALSE)
polygon(x = c(rownames(psychRes$dtfSCI$sDTF32), rev(rownames(psychRes$dtfSCI$sDTF32))),
        y = c(psychRes$dtfSCI$sDTF32[, 2], rev(psychRes$dtfSCI$sDTF32[, 3])) * -1, 
        col = adjustcolor("black", alpha.f = .1), border = NA)
points(rownames(psychRes$dtfSCI$sDTF32),
       rep(0, nrow(psychRes$dtfSCI$sDTF32)), 
       type = "l", lty = 3,
       col = "gray42")
axis(1, at = seq(-8, 8, 2))
axis(2, at = seq(-1, 3, 1))
dev.off()



# ------------------------------------------------------------------------------
# Results for self-selection analyses
# ------------------------------------------------------------------------------

# load results
load(paste0("Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/results/_final-", version, ".RData"))

# Nagelkerke's R2, for models describing self-selection
round(mean(pseudoMR2),2) 

# Nagelkerke's R2, for non-response models used to derive the propensity scores
round(apply(pseudoR2,1,mean),2)  




# ------------------------------------------------------------------------------
# Results for differential test functioning
# ------------------------------------------------------------------------------

# latent variance
round(psychRes$varSCI, 2)

# reliabilities
round(psychRes$relSCI, 2)

# infit
lapply(psychRes$fitSCI, round, digits = 2)
lapply(psychRes$fitSCI, psych::describe)

# sDTF
round(psychRes$dtfSCI$sDTF, 2)

# uDTF
round(psychRes$dtfSCI$uDTF, 2)

# percentage uDTF
round(psychRes$dtfSCI$uDTFp, 2)



# ------------------------------------------------------------------------------
# Results for analyses of prediction error
# ------------------------------------------------------------------------------

# load results
load(paste0("Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/results/_final-", version, ".RData"))

# GPA
round(res_modP$pred10, 2)
round(res_modP$pred11, 2)

# helplessness
round(res_modP$pred20, 2) 
round(res_modP$pred21, 2)

# self-concept
round(res_modP$pred30, 2)
round(res_modP$pred31, 2)

# study dropout
round(res_modP$pred40, 2)
round(res_modP$pred41, 2)

# R2
round(sapply(r2_sc, mean), 2)

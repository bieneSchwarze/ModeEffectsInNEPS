############################################################################
############################################################################
##                                                                        ##
## Prepare data for studying mode effects in competence tests SC5 Wave 5. ##
## The studied modes are PBA, CBA, and WBA.                               ##
## This includes multiple imputation of missing values.                   ##
##                                                                        ##
## authors: Sabine & Timo                                                 ##
## date: 2020-02-23                                                       ##
##                                                                        ##
############################################################################
############################################################################



# empty work space
rm(list=ls())



# ----------------------------------------------------------------------------
# Load libraries
# ----------------------------------------------------------------------------

library(mice)
library(readstata13)
library(haven)
library(dplyr)
library(doBy)
library(TAM)
library(MASS)
library(statmod)
library(ggplot2)
library(MBESS)
library(lavaan)
library(rms)



# ----------------------------------------------------------------------------
# Read & edit design and imputation data
# ----------------------------------------------------------------------------

setwd("Z:\\Projects\\p000139_Methoden_Survey_Psych\\B57Project\\data")

# Methods data containing information about all modes, 
#   assignments to modes and whether the students attended in the tests.
meth <- read.dta13("MethodsCompetenceTestsWave5SC5.dta") # This data set should be available on request from the LIfBi FDZ.
table(meth$PBA_E); table(meth$PBA_R)
table(meth$CBA_E); table(meth$CBA_R)
table(meth$WBA_E); table(meth$WBA_R)
table(meth$WBA_switch_E); table(meth$WBA_switch_R)

# Load SUF data.
setwd("Z:\\SUF\\On-site\\SC5\\SC5_O_12-0-0\\SPSS\\en")
basic <- read_sav("SC5_Basics_O_12-0-0.sav")
targetCATI <- read_sav("SC5_pTargetCATI_O_12-0-0.sav")
targetCAWI <- read_sav("SC5_pTargetCAWI_O_12-0-0.sav")
institution <- read_sav("SC5_xInstitution_O_12-0-0.sav")
cohortProfile <- read_sav("SC5_CohortProfile_O_12-0-0.sav")
comp <- read_sav("SC5_xTargetCompetencies_O_12-0-0.sav")
xx <- comp[!(comp$ID_t %in% meth$ID_t),]; nrow(xx) # N=49 cases not in MethodsData -> These N=49 are all final dropouts before the Wave 3 competence tests.

# Get participation status from CohortProfile.
w1_part <- cohortProfile[cohortProfile$wave %in% 1, c("ID_t", "ID_i", "tx80522")] # Competence data available in Wave 1.
table(w1_part$tx80522)
attr(w1_part$tx80522, "labels") <- NULL
w2_part <- cohortProfile[cohortProfile$wave %in% 2, c("ID_t", "tx80220")] # Participated in Wave 2 (CAWI): yes: 1, temp. dropout: 2, fin. dropout: 3
table(w2_part$tx80220)
attr(w2_part$tx80220, "labels") <- NULL
w3_part <- cohortProfile[cohortProfile$wave %in% 3, c("ID_t", "tx80220")] # Participated in Wave 3 (CATI): yes: 1, temp. dropout: 2, fin. dropout: 3
table(w3_part$tx80220)
attr(w3_part$tx80220, "labels") <- NULL
w4_part <- cohortProfile[cohortProfile$wave %in% 4, c("ID_t", "tx80220")] # Participated in Wave 4 (CAWI): yes: 1, temp. dropout: 2, fin. dropout: 3
table(w4_part$tx80220)
attr(w4_part$tx80220, "labels") <- NULL
w5_part <- cohortProfile[cohortProfile$wave %in% 5, c("ID_t", "tx80220")] # Participated in Wave 5 (CATI): yes: 1, temp. dropout: 2, fin. dropout: 3
table(w5_part$tx80220)
attr(w5_part$tx80220, "labels") <- NULL
# Risk set for study is defined as such not dropping out until Wave 4. 
# Explanation: Competence test Wave 5 was conducted in parallel to CATI Wave 5. Thus, students who dropped out during CATI Wave 5 could still participate in competence test Wave 5.
idsRisk <- w4_part[w4_part$tx80220 %in% c(1,2),"ID_t"]
length(idsRisk$ID_t) # N=17626
partSet <- data.frame(ID_t=idsRisk, w1CompPart=w1_part[w1_part$ID_t %in% idsRisk$ID_t,3], w2PartCAWI=w2_part[w2_part$ID_t %in% idsRisk$ID_t,2], 
                      w3PartCATI=w3_part[w3_part$ID_t %in% idsRisk$ID_t,2],  w4PartCAWI=w4_part[w4_part$ID_t %in% idsRisk$ID_t,2], 
                      ID_i=w1_part[w1_part$ID_t %in% idsRisk$ID_t,2])
colnames(partSet) <- c("ID_t", "w1CompPart", "w2PartCAWI", "w3PartCATI", "w4PartCAWI", "ID_i")
partSet$w2PartCAWI <- recodeVar(partSet$w2PartCAWI, unique(partSet$w2PartCAWI), c(1,0))
partSet$w3PartCATI <- recodeVar(partSet$w3PartCATI, unique(partSet$w3PartCATI), c(0,1))
partSet$w4PartCAWI <- recodeVar(partSet$w4PartCAWI, unique(partSet$w4PartCAWI), c(1,0))
partSet$ID_i[partSet$ID_i<0] <- NA
table(is.na(partSet$ID_i))

# Anybody in methods data who already dropped out until Wave 4?
table(!(meth$ID_t %in% partSet$ID_t)) # No. 

# For each student in the panel cohort, is there an entry in the Wave 5 competence methods data set? 
table(!(idsRisk$ID_t %in%  meth$ID_t)) # No. For N=153 students we have no info on the Wave 5 competence test.
# -> We have to remove these students from our analysis.
partSet <- partSet[partSet$ID_t %in% meth$ID_t,] # Remain: N=17473.
partSet <- merge(partSet, meth, by="ID_t")

# Check whether competence data available 
w5_comp <- cohortProfile[cohortProfile$wave %in% 5, c("ID_t", "tx80522")] # Competence data available in Wave 5.
w5_comp_ids <- w5_comp[w5_comp$tx80522 %in% 1, "ID_t"]
table(w5_comp_ids$ID_t %in% comp$ID_t) # Okay.

# Get variables from Basic data set.
varsBasic <- c("ID_t", "t70000y", "t700001", "tx29005", "tx29003", "t751001_g1", "t731301_g3", "t731351_g3")
basic <- basic[, varsBasic]
basic$t700001 <- recodeVar(basic$t700001, unique(basic$t700001), c(1,0))
attr(basic$t700001, "labels") <- NULL
basic$tx29005 <- recodeVar(basic$tx29005, unique(basic$tx29005), c(1,0))
attr(basic$tx29005, "labels") <- NULL
basic$tx29003 <- recodeVar(basic$tx29003, unique(basic$tx29003), c(1,0, NA))
attr(basic$tx29003, "labels") <- NULL
basic$t751001_g1 <- recodeVar(basic$t751001_g1, unique(basic$t751001_g1), c(0,1, NA))
attr(basic$t751001_g1, "labels") <- NULL
colnames(basic) <- c("ID_t", "birthYear", "female", "birthInGermany", "mothTongueGerman", "eastGermany", "yearsEduMother", "yearsEduFather") 

# Get data from institution data set.
varsInst <- c("ID_i", "tg92301_O", "tg92401_O")
institution <- institution[, varsInst] 
institution$tg92301_O <- recodeVar(institution$tg92301_O, attributes(institution$tg92301_O)$labels, c(1,1,0)) # public uni & clerical: 1, private & clerical: 0
attr(institution$tg92301_O, "labels") <- NULL
table(institution$tg92301_O, exclude=NULL)
institution$tg92401_O <- recodeVar(institution$tg92401_O, attributes(institution$tg92401_O)$labels, c(0,2,2,2,1,2,2,2)) # uni: 0, applied science: 1, remaining: 2
attr(institution$tg92401_O, "labels") <- NULL
table(institution$tg92401_O, exclude=NULL)
colnames(institution) <- c("ID_i", "uniFunding", "uniType")
institution <- institution[order(institution$ID_i),]
institution <- institution[!duplicated(institution$ID_i),] # Take the first entry for each institution. (Multiple rows belong the the distinct study fields. Irrelevant for us.)

# Get data from CATI data (1st Wave).
varsCATI_W1 <- c("ID_t", "tg04001_g1R", "tg24201_g1", "t30142a", "t300400", "t743022", "t743023")
targetCATI_W1 <- targetCATI[targetCATI$wave %in% 1, varsCATI_W1] 
targetCATI_W1$tg04001_g1R <- recodeVar(targetCATI_W1$tg04001_g1R, attributes(targetCATI_W1$tg04001_g1R)$labels, 
                                       c(rep(NA,4), rep(0,27), rep(1,9), rep(0,24), rep(2,4)))
table(targetCATI_W1$tg04001_g1R, exclude=NULL)
attr(targetCATI_W1$tg04001_g1R, "labels") <- NULL
targetCATI_W1$tg24201_g1 <- recodeVar(targetCATI_W1$tg24201_g1, attributes(targetCATI_W1$tg24201_g1)$labels, c(NA,NA,1,0))
table(targetCATI_W1$tg24201_g1, exclude=NULL)
attr(targetCATI_W1$tg24201_g1, "labels") <- NULL
targetCATI_W1$t30142a <- recodeVar(targetCATI_W1$t30142a, attributes(targetCATI_W1$t30142a)$labels, c(NA,NA,NA,1:5))
table(targetCATI_W1$t30142a, exclude=NULL)
attr(targetCATI_W1$t30142a, "labels") <- NULL
targetCATI_W1$t300400 <- recodeVar(targetCATI_W1$t300400, attributes(targetCATI_W1$t300400)$labels, c(NA,NA,NA,1:5))
table(targetCATI_W1$t300400, exclude=NULL)
attr(targetCATI_W1$t300400, "labels") <- NULL
targetCATI_W1$t743022 <- recodeVar(targetCATI_W1$t743022, attributes(targetCATI_W1$t743022)$labels, c(rep(NA,5),0,1))
table(targetCATI_W1$t743022, exclude=NULL)
attr(targetCATI_W1$t743022, "labels") <- NULL
targetCATI_W1$t743023 <- recodeVar(targetCATI_W1$t743023, attributes(targetCATI_W1$t743023)$labels, c(rep(NA,5),0,1))
table(targetCATI_W1$t743023, exclude=NULL)
attr(targetCATI_W1$t743023, "labels") <- NULL
targetCATI_W1$kidsHH_W1 <- ifelse(targetCATI_W1$t743022 %in% 1 | targetCATI_W1$t743023 %in% 1, 1, 
                               ifelse(is.na(targetCATI_W1$t743022) | is.na(targetCATI_W1$t743022), NA,0))
table(targetCATI_W1$kidsHH_W1, exclude=NULL)
targetCATI_W1 <- targetCATI_W1[,!(colnames(targetCATI_W1) %in% c("t743022","t743023"))]
colnames(targetCATI_W1) <- c("ID_t", "naturalScience", "teacherEdu", "enjoyStud", "probGrad", "kidsHH_W1")

# Get data from CATI data (3rd Wave).
varsCATI_W3 <- c("ID_t", "t66800a_g1", "t66800b_g1", "t66800c_g1", "t66800d_g1", "t66800e_g1", "t66003a_g1", "t743022", "t743023") 
targetCATI_W3 <- targetCATI[targetCATI$wave %in% 3, varsCATI_W3] 
table(targetCATI_W3$t66800a_g1, exclude=NULL)
table(targetCATI_W3$t66800b_g1, exclude=NULL)
table(targetCATI_W3$t66800c_g1, exclude=NULL)
table(targetCATI_W3$t66800d_g1, exclude=NULL)
table(targetCATI_W3$t66800e_g1, exclude=NULL)
table(targetCATI_W3$t66003a_g1, exclude=NULL)
targetCATI_W3$t743022 <- recodeVar(targetCATI_W3$t743022, attributes(targetCATI_W3$t743022)$labels, c(rep(NA,5),0,1))
table(targetCATI_W3$t743022, exclude=NULL)
attr(targetCATI_W3$t743022, "labels") <- NULL
targetCATI_W3$t743023 <- recodeVar(targetCATI_W3$t743023, attributes(targetCATI_W3$t743023)$labels, c(rep(NA,5),0,1))
table(targetCATI_W3$t743023, exclude=NULL)
attr(targetCATI_W3$t743023, "labels") <- NULL
targetCATI_W3$kidsHH_W3 <- ifelse(targetCATI_W3$t743022 %in% 1 | targetCATI_W3$t743023 %in% 1, 1, 
                                  ifelse(is.na(targetCATI_W3$t743022) | is.na(targetCATI_W3$t743022), NA,0))
table(targetCATI_W3$kidsHH_W3, exclude=NULL)
targetCATI_W3 <- targetCATI_W3[,!(colnames(targetCATI_W3) %in% c("t743022","t743023"))]
colnames(targetCATI_W3) <- c("ID_t", "extraversion", "agreeableness", "conscientiousness", "neuroticism", "openess", "selfEsteem", "kidsHH_W3")

# Get data from CATI data (5th Wave).
varsCATI_W5 <- c("ID_t", "t743022","t743023", "t741001", "tf11105") 
targetCATI_W5 <- targetCATI[targetCATI$wave %in% 5, varsCATI_W5] 
targetCATI_W5$t743022 <- recodeVar(targetCATI_W5$t743022, attributes(targetCATI_W5$t743022)$labels, c(rep(NA,5),0,1))
table(targetCATI_W5$t743022, exclude=NULL)
attr(targetCATI_W5$t743022, "labels") <- NULL
targetCATI_W5$t743023 <- recodeVar(targetCATI_W5$t743023, attributes(targetCATI_W5$t743023)$labels, c(rep(NA,5),0,1))
table(targetCATI_W5$t743023, exclude=NULL)
attr(targetCATI_W5$t743023, "labels") <- NULL
targetCATI_W5$kidsHH_W5 <- ifelse(targetCATI_W5$t743022 %in% 1 | targetCATI_W5$t743023 %in% 1, 1, 
                                  ifelse(is.na(targetCATI_W5$t743022) | is.na(targetCATI_W5$t743022), NA,0))
table(targetCATI_W5$kidsHH_W5, exclude=NULL)
targetCATI_W5 <- targetCATI_W5[,!(colnames(targetCATI_W5) %in% c("t743022","t743023"))]
attr(targetCATI_W5$tf11105, "labels") <- NULL
colnames(targetCATI_W5) <- c("ID_t", "hhSize_W5", "DegreeAbroad_noDegree", "kidsHH_W5") # Degree: no: 0, in Germany: 1, abroad: 2

# Combine CATI data.
Dat <- merge(partSet, targetCATI_W1, by="ID_t", all.x=TRUE)
Dat <- merge(Dat, targetCATI_W3, by="ID_t", all.x=TRUE)
Dat <- merge(Dat, targetCATI_W5, by="ID_t", all.x=TRUE)

# Got data from CAWI data from 2nd Wave.
varsCAWI_W2 <- c("ID_t", "tg52020") 
targetCAWI_W2 <- targetCAWI[targetCAWI$wave %in% 2, varsCAWI_W2] 
colnames(targetCAWI_W2) <- c("ID_t", "GPA_W2")

# Got data from CAWI data from 4th Wave.
varsCAWI_W4 <- c("ID_t", "tg52020") 
targetCAWI_W4 <- targetCAWI[targetCAWI$wave %in% 4, varsCAWI_W4] 
colnames(targetCAWI_W4) <- c("ID_t", "GPA_W4")

# Got data from CAWI data from 6th Wave.
varsCAWI_W6 <- c("ID_t", "tg52020", # GPA
                 "tg53221", "tg53222", "tg53223", # intention to quit
                 "tg53224", "tg53225",
                 "t66010d", "t66010b", "t66010a", # helplessness
                 "t66007a", "t66007b", "t66007d", "t66007e")  # self-concept
targetCAWI_W6 <- targetCAWI[targetCAWI$wave %in% 6, varsCAWI_W6] 
colnames(targetCAWI_W6) <- c("ID_t", "GPA_W6",
                             paste0("studydropout", 1:5),
                             paste0("helpless", 1:3),
                             paste0("selfconcept", 1:4))
targetCAWI_W6$studydropout5 <- 5- targetCAWI_W6$studydropout5
targetCAWI_W6$studydropout <- rowMeans(targetCAWI_W6[, paste0("studydropout", 1:5)])
targetCAWI_W6$helpless <- rowMeans(targetCAWI_W6[, paste0("helpless", 1:3)])
targetCAWI_W6$selfconcept <- rowMeans(targetCAWI_W6[, paste0("selfconcept", 1:4)])
ci.reliability(targetCAWI_W6[, paste0("studydropout", 1:5)], 
               type ="categorical", interval.type = "none")$est
ci.reliability(targetCAWI_W6[, paste0("helpless", 1:3)], 
               type ="categorical", interval.type = "none")$est
ci.reliability(targetCAWI_W6[, paste0("selfconcept", 1:4)], 
               type ="categorical", interval.type = "none")$est
targetCAWI_W6[, c(paste0("studydropout", 1:5),
                  paste0("helpless", 1:3),
                  paste0("selfconcept", 1:4))] <- NULL

# Combine CAWI data.
Dat_o <- merge(idsRisk, targetCAWI_W2, by="ID_t", all.x=TRUE)
Dat_o <- merge(Dat_o, targetCAWI_W4, by="ID_t", all.x=TRUE)
Dat_o <- merge(Dat_o, targetCAWI_W6, by="ID_t", all.x=TRUE)

# Combine BASIC, CATI, CAWI, and institution data.
Dat <- merge(Dat, Dat_o, by="ID_t")
Dat <- merge(Dat, basic, by="ID_t", all.x=TRUE)
table(Dat$ID_i %in% institution$ID_i)
instIDs <- unique(Dat$ID_i)
instIDs <- instIDs[!(is.na(instIDs))]
length(instIDs) # N=225
institution <- institution[institution$ID_i %in% instIDs,]
setdiff(instIDs, unique(institution$ID_i)) # No info for two of the institutions.
Dat <- merge(Dat, institution, by="ID_i", all.x=TRUE)
Dat <- Dat[order(Dat$ID_t),]    
   


# ----------------------------------------------------------------------------
# Read & edit competence data
# ----------------------------------------------------------------------------

# Variable names for scientific literacy test.
citems <- list(sci = c("scs36310_c", "scs36320_c", "scs36220_c", "scs3623s_c", "scs30510_c", "scs30520_c",
                       "scs31210_c", "scs31220_c", "scs31240_c", "scs30920_c", "scs30930_c", "scs30940_c",
                       "scs3021s_c", "scs3022s_c", "scs36020_c", "scs3643s_c", "scs3642s_c", "scs3031s_c",
                       "scs3033s_c", "scs3112s_c", "scs3131s_c", "scs3132s_c", "scs3133s_c", "scs3012s_c",
                       "scs30130_c", "scs3061s_c", "scs30630_c", "scs30640_c", "scs30810_c"))

# Collapse science items with less than 200 respondents (see Pohl & Carstenen, 2011).
comp$scs3061s_c <- recodeVar(comp$scs3061s_c, 0:4, c(0, 0, 1, 2, 3))
comp$scs3012s_c <- recodeVar(comp$scs3012s_c, 0:4, c(0, 0, 1, 2, 3))
comp$scs3133s_c <- recodeVar(comp$scs3133s_c, 0:4, c(0, 0, 1, 2, 3))
comp$scs3132s_c <- recodeVar(comp$scs3132s_c, 0:4, c(0, 0, 1, 2, 3))

# Recode mode based on tx-codes
comp$mode <- recodeVar(comp$tx80211_w5, 330:337, c(rep("PBA", 4), "CBA", "CBA", "WBA", "WBA"))
attr(comp$mode, "labels") <- NULL
table(comp$mode, useNA = "always")

# Select relevant variables.
comp <- dplyr::select(comp, ID_t, mode, mas1_sc1u, scs3_sc1u, one_of(citems$sci))

# Check response distributions.
sapply(citems$sci, function(x) { table(comp[[x]], useNA = "always") }) # science

# Missing indicators.
comp$inComp <- 1 # case is in competence file
comp$misSCI <- rowSums(is.na(comp[, citems$sci])) == length(citems$sci) # missing science
table(comp$misSCI)
comp$misSCIwle <- ifelse(is.na(comp$scs3_sc1u),TRUE,FALSE) # mark those students who have missing in the science wle 

# Add previous test score to Dat.
Dat <- merge(Dat, comp[,c("ID_t", "mas1_sc1u")], by="ID_t", all.x = TRUE)


  
# ----------------------------------------------------------------------------
# Check response status against available competence data
# ----------------------------------------------------------------------------

comp$ID_t <- as.integer(comp$ID_t)
table(Dat$ID_t %in% comp$ID_t) # In total N=5711 of the Wave 5 panel members have no competence value measured.
table(comp$ID_t %in% Dat$ID_t) # N=49 students in competence data set are final dropouts until Wave 4.

# Check PBA participants.
idsPBA_R <- Dat[Dat$PBA_R %in% 1, "ID_t"] 
idsPBA_c <- comp[comp$mode %in% "PBA", ]$ID_t
table(idsPBA_R %in% idsPBA_c)
table(idsPBA_c %in% idsPBA_R) # Okay.

# Check CBA participants.
idsCBA_R <- Dat[Dat$CBA_R %in% 1, "ID_t"] 
idsCBA_c <- comp[comp$mode %in% "CBA", ]$ID_t
table(idsCBA_R %in% idsCBA_c)
table(idsCBA_c %in% idsCBA_R) # Okay.

# Check WBA participants.
idsWBA_R <- Dat[Dat$WBA_R == 1 | Dat$WBA_switch_R == 1, "ID_t"] 
idsWBA_c <- comp[comp$mode %in% "WBA", ]$ID_t
table(idsWBA_R %in% idsWBA_c)
table(idsWBA_c %in% idsWBA_R) # Okay.

# Check whether all participants have competence scores in science.
# If not they are non-participants.
part <- apply(Dat[, c("PBA_R", "CBA_R", "WBA_R", "WBA_switch_R")],1,sum)
idsPart <- Dat[part==1,"ID_t"]
length(idsPart)

idsMisSCI <- comp[comp$misSCI,]$ID_t
length(idsMisSCI)
part <- apply(Dat[, c("PBA_R", "CBA_R", "WBA_R", "WBA_switch_R")],1,sum)
idsPart <- Dat[part==1,"ID_t"]
length(idsPart)
table(idsMisSCI %in% idsPart) # N=161 students who have partic==Yes in Dat, have to science score. -> Set their participation status to zero.
Dat[Dat$ID_t %in%  idsPart[idsPart %in% idsMisSCI], c("PBA_R", "CBA_R", "WBA_R", "WBA_switch_R")] <- 0

idsMisSCIwle <- comp[comp$misSCIwle,]$ID_t
length(idsMisSCIwle)
part <- apply(Dat[, c("PBA_R", "CBA_R", "WBA_R", "WBA_switch_R")],1,sum)
idsPart <- Dat[part==1,"ID_t"]
length(idsPart)
table(idsMisSCIwle %in% idsPart) # There are still N=164 cases in the participants data with missing wle in science. 
# idsMisSCIwle[idsMisSCIwle %in% idsPart] # TODO TG: why?
Dat[Dat$ID_t %in%  idsPart[idsPart %in% idsMisSCIwle], c("PBA_R", "CBA_R", "WBA_R", "WBA_switch_R")] <- 0

# Assign proper data type.
Dat$uniFunding <- as.factor(Dat$uniFunding)
Dat$uniType <- as.factor(Dat$uniType)
Dat$teacherEdu <- as.factor(Dat$teacherEdu)
Dat$enjoyStud <- as.numeric(Dat$enjoyStud)
Dat$naturalScience <- as.factor(Dat$naturalScience)
Dat$probGrad <- as.numeric(Dat$probGrad)
Dat$kidsHH_W1 <- as.factor(Dat$kidsHH_W1)
Dat$kidsHH_W3 <- as.factor(Dat$kidsHH_W3)
Dat$kidsHH_W5 <- as.factor(Dat$kidsHH_W5)
Dat$extraversion <- as.numeric(Dat$extraversion)
Dat$conscientiousness <- as.numeric(Dat$conscientiousness)
Dat$neuroticism <- as.numeric(Dat$neuroticism)
Dat$openess <- as.numeric(Dat$openess)
Dat$agreeableness <- as.numeric(Dat$agreeableness)
Dat$selfEsteem <- as.numeric(Dat$selfEsteem)
Dat$DegreeAbroad_noDegree <- as.factor(Dat$DegreeAbroad_noDegree)
Dat$GPA_W2 <- as.numeric(Dat$GPA_W2)
Dat$GPA_W2[Dat$GPA_W2 > 5] <- NA
Dat$GPA_W4 <- as.numeric(Dat$GPA_W4)
Dat$GPA_W4[Dat$GPA_W4 > 5] <- NA
Dat$GPA_W6 <- as.numeric(Dat$GPA_W6)
Dat$GPA_W6[Dat$GPA_W6 == 0 | Dat$GPA_W6 > 5] <- NA
Dat$studydropout <- as.numeric(Dat$studydropout)
Dat$helpless <- as.numeric(Dat$helpless)
Dat$selfconcept <- as.numeric(Dat$selfconcept)
Dat$mas1_sc1u <- as.numeric(Dat$mas1_sc1u)
Dat$hhSize_W5 <- ifelse(is.na(Dat$hhSize_W5), NA, ifelse(Dat$hhSize_W5 %in% 1, 0, ifelse(Dat$hhSize_W5 %in% 2, 1, 2)))
Dat$hhSize_W5 <- as.factor(Dat$hhSize_W5)



# -------------------------------------------------------------------------------
# Descriptives
# -------------------------------------------------------------------------------

propMiss <- function(vv){
  tab <- prop.table(table(is.na(vv)))
  return(tab[names(tab)=="TRUE"])
}

# A. Participation of students who were initially invited to group testing PBA, derive p^gr_TN
table(Dat$PBA_E); table(Dat$PBA_R)

# B. Participation of students who were initially invited to group testing CBA, derive p^gr_TN  
table(Dat$CBA_E); table(Dat$CBA_R)

# C. Participation of students who were invited to online testing from the beginning, derive p^on_TN
table(Dat$WBA_E); table(Dat$WBA_R)

# D. Participation of students who did not participate in group testing (but were invited) and were asked to switch to online testing, derive p^sw_TN
table(Dat$WBA_switch_E); table(Dat$WBA_switch_R)

# Descriptives for model variables
getStats <- function(DD){
  for(ii in 1:ncol(DD)){
    v <- DD[,ii]
    if(is.factor(v)){
      u <- unique(v[!is.na(v)])
      if(length(u)==2){ # Coding: factor has two levels (0 and 1)
        v <- as.numeric(as.character(v))
        mE <- round(mean(v, na.rm = TRUE),2)
        miP <- round(propMiss(v)*100,2)
        cat(colnames(DD)[ii],": Prop: ",mE, " -- MissP: ", miP, "\n")        
      } else { # Coding: factor has three levels (0, 1, 2) -> This is special to our case, max. 3 levels.
        miP <- round(propMiss(as.numeric(as.character(v)))*100,2)
        v1 <- ifelse(as.numeric(as.character(v)) %in% 1,1,0)
        v2 <- ifelse(as.numeric(as.character(v)) %in% 2,1,0)
        mE1 <- round(mean(v1, na.rm = TRUE),2)
        mE2 <- round(mean(v2, na.rm = TRUE),2)
        cat(colnames(DD)[ii],": Prop Fact 1: ",mE1, ": Prop Fact 2: ",mE2, " -- MissP: ", miP, "\n")  
      }
    } else {
      mE <- round(mean(v, na.rm = TRUE),2)
      std <- round(sqrt(var(v, na.rm=TRUE)),2)
      miP <- round(propMiss(v)*100,2)
      cat(colnames(DD)[ii],": Mean: ",mE, "  -- StdE: ", std, " -- MissP: ", miP, "\n")
    }
  }  
} 

getStats(Dat[Dat$PBA_E %in% 1,])
D1 <- Dat[Dat$PBA_E %in% 1,]
D1$by <- ifelse(D1$birthYear < 1989, 0, ifelse(D1$birthYear <=1990, 1, 2))       
round(table(D1$by, exclude = NULL)/nrow(D1),2)

getStats(Dat[Dat$CBA_E %in% 1,])
D2 <- Dat[Dat$CBA_E %in% 1,]
D2$by <- ifelse(D2$birthYear < 1989, 0, ifelse(D2$birthYear <=1990, 1, 2))       
round(table(D2$by, exclude = NULL)/nrow(D2),2)

getStats(Dat[Dat$WBA_E %in% 1,])
D3 <- Dat[Dat$WBA_E %in% 1,]
D3$by <- ifelse(D3$birthYear < 1989, 0, ifelse(D3$birthYear <=1990, 1, 2))       
round(table(D3$by, exclude = NULL)/nrow(D3),2)

getStats(Dat[Dat$WBA_switch_E %in% 1,])
D4 <- Dat[Dat$WBA_switch_E %in% 1,]
D4$by <- ifelse(D4$birthYear < 1989, 0, ifelse(D4$birthYear <=1990, 1, 2))       
round(table(D4$by, exclude = NULL)/nrow(D4),2)

getStats(Dat[Dat$naturalScience %in% 1,])



# ----------------------------------------------------------------------------
# Missing data handling and imputation
# ----------------------------------------------------------------------------

HD <- md.pattern(Dat, plot=FALSE)
round(HD[nrow(HD),]/nrow(Dat)*100,2) # proportion of missing values

predM <- mice::make.predictorMatrix(data=Dat)
impM <- mice::make.method(data=Dat)
predM1 <- predM
predM1[,"ID_t"] <- 0
predM1[,"ID_i"] <- 0
predM1["ID_t",] <- 0
predM1["ID_i",] <- 0
impM1 <- impM 
imp1 <- mice::mice(Dat, predictorMatrix = predM1, method = impM1, print = FALSE,
                   maxit=30, m=20, seed=145) # use: m=20 & maxit=30

# -> GPA 4th semester has so many missings that I get a warning during imputation. 
# Thus, for further analyses I recommend to let this variable aside.


# ----------------------------------------------------------------------------
# Save workspace
# ----------------------------------------------------------------------------

save(Dat, comp, citems, imp1, file = "Z:/Projects/p000139_Methoden_Survey_Psych/B57Project/results/_dataPrepared-20200223.RData")



# ------------------------------------------------------------------------------
# Cohort description -----------------------------------------------------------
# ------------------------------------------------------------------------------

# standardized mean difference: https://www.rdatagen.net/post/2023-09-26-nice-looking-table-1-with-standardized-mean-difference/

#
# Packages ---------------------------------------------------------------------
#


# install.packages("table1")
library("table1")                                                               # table 1

# install.packages("stringr")                
library("stringr")                                                              # to create the vairable for pracid


#
# Set directions ---------------------------------------------------------------
#


data_dir <- "C:/Users/lg704/OneDrive - University of Exeter/Studies/SGLT2i_ethnicity_wJohn/data/"


#
# Load data --------------------------------------------------------------------
#

load(paste0(data_dir,"final.dataset.sglt2.dpp4.val.Rda"))

load(paste0(data_dir,"final.hb.Rdata"))
load(paste0(data_dir,"final.dc.Rdata"))
load(paste0(data_dir,"final.wt.Rdata"))

final.hb_analysis <- final.hb
final.wt_analysis <- final.wt
final.dc_analysis <- final.dc

rm(final.hb)
rm(final.wt)
rm(final.dc)



#
# Merge to include needed variables --------------------------------------------
#


final.dataset$patiddrugclass       <- paste0(final.dataset$patid, final.dataset$drugclass)
final.hb_analysis$patiddrugclass   <- paste0(final.hb_analysis$patid, final.hb_analysis$drugclass)
final.wt_analysis$patiddrugclass   <- paste0(final.wt_analysis$patid, final.wt_analysis$drugclass)
final.dc_analysis$patiddrugclass   <- paste0(final.dc_analysis$patid, final.dc_analysis$drugclass)


added_variables <- c("patiddrugclass","drugsubstances", "MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1",
                      "sex", "prehdl", "pretriglyceride", "prealbuminblood", "prebilirubin", "t2dmduration")


final.hb <- merge(final.hb_analysis, final.dataset[ , added_variables], by = "patiddrugclass")
final.wt <- merge(final.wt_analysis, final.dataset[ , added_variables], by = "patiddrugclass")
final.dc <- merge(final.dc_analysis, final.dataset[ , added_variables], by = "patiddrugclass")


#
# Data preparation -------------------------------------------------------------
#

final.hb$prealt <- exp(final.hb$prealtlog)
final.dc$prealt <- exp(final.dc$prealtlog)
final.wt$prealt <- exp(final.wt$prealtlog)


final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Alogliptin & Linagliptin",      "Alogliptin")
final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Alogliptin & Sitagliptin",      "Alogliptin")
final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Canagliflozin & Dapagliflozin", "Canagliflozin")
final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Canagliflozin & Empagliflozin", "Canagliflozin")
final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Dapagliflozin & Empagliflozin", "Dapagliflozin")
final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Linagliptin & Saxagliptin",      "Linagliptin")
final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Linagliptin & Sitagliptin",      "Linagliptin")
final.hb$drugsubstances <- replace(final.hb$drugsubstances, final.hb$drugsubstances == "Sitagliptin & Vildagliptin",     "Sitagliptin")


final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Alogliptin & Linagliptin",      "Alogliptin")
final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Alogliptin & Sitagliptin",      "Alogliptin")
final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Canagliflozin & Dapagliflozin", "Canagliflozin")
final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Canagliflozin & Empagliflozin", "Canagliflozin")
final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Dapagliflozin & Empagliflozin", "Dapagliflozin")
final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Linagliptin & Saxagliptin",      "Linagliptin")
final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Linagliptin & Sitagliptin",      "Linagliptin")
final.dc$drugsubstances <- replace(final.dc$drugsubstances, final.dc$drugsubstances == "Sitagliptin & Vildagliptin",     "Sitagliptin")


final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Alogliptin & Linagliptin",      "Alogliptin")
final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Alogliptin & Sitagliptin",      "Alogliptin")
final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Canagliflozin & Dapagliflozin", "Canagliflozin")
final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Canagliflozin & Empagliflozin", "Canagliflozin")
final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Dapagliflozin & Empagliflozin", "Dapagliflozin")
final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Linagliptin & Saxagliptin",      "Linagliptin")
final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Linagliptin & Sitagliptin",      "Linagliptin")
final.wt$drugsubstances <- replace(final.wt$drugsubstances, final.wt$drugsubstances == "Sitagliptin & Vildagliptin",     "Sitagliptin")

final.hb$drugline_chr   <- as.character(final.hb$drugline)
final.hb$drugline_chr   <- ifelse(final.hb$drugline_chr == "4" | final.hb$drugline_chr == "5", "4+", final.hb$drugline_chr)

final.dc$drugline_chr   <- as.character(final.dc$drugline)
final.dc$drugline_chr   <- ifelse(final.dc$drugline_chr == "4" | final.dc$drugline_chr == "5", "4+", final.dc$drugline_chr)

final.wt$drugline_chr   <- as.character(final.wt$drugline)
final.wt$drugline_chr   <- ifelse(final.wt$drugline_chr == "4" | final.wt$drugline_chr == "5", "4+", final.wt$drugline_chr)


final.dc$stopdrug_6m_3mFU_chr   <- as.character(final.dc$stopdrug_6m_3mFU)

length(unique(str_sub(final.hb$patid, - 3, - 1)))                               # number of practice 

    
final.hb$MFN    <- as.character(final.hb$MFN)
final.hb$SU     <- as.character(final.hb$SU)
final.hb$DPP4   <- as.character(final.hb$DPP4)
final.hb$SGLT2  <- as.character(final.hb$SGLT2)
final.hb$DPP4   <- as.character(final.hb$DPP4)
final.hb$TZD    <- as.character(final.hb$TZD)
final.hb$GLP1   <- as.character(final.hb$GLP1)

final.dc$MFN    <- as.character(final.dc$MFN)
final.dc$SU     <- as.character(final.dc$SU)
final.dc$DPP4   <- as.character(final.dc$DPP4)
final.dc$SGLT2  <- as.character(final.dc$SGLT2)
final.dc$DPP4   <- as.character(final.dc$DPP4)
final.dc$TZD    <- as.character(final.dc$TZD)
final.dc$GLP1   <- as.character(final.dc$GLP1)

final.wt$MFN    <- as.character(final.wt$MFN)
final.wt$SU     <- as.character(final.wt$SU)
final.wt$DPP4   <- as.character(final.wt$DPP4)
final.wt$SGLT2  <- as.character(final.wt$SGLT2)
final.wt$DPP4   <- as.character(final.wt$DPP4)
final.wt$TZD    <- as.character(final.wt$TZD)
final.wt$GLP1   <- as.character(final.wt$GLP1)


#
# Table 1 setup ----------------------------------------------------------------
#

variables     <- c("agetx", "t2dmduration", "sex", "ethnicity","drugsubstances", "ncurrtx", "drugline_chr", "prehba1cmmol", "prebmi", 
                   "egfr_ckdepi", "prealt", "prehdl", "pretriglyceride", "prealtlog", "prealt", "prealbuminblood", "prebilirubin",
                   "MFN", "SU", "DPP4", "SGLT2", "TZD", "GLP1")

variables.hb  <- c(variables, "posthba1c_final", "hba1cmonth")

variables.wt  <- c(variables, "postweightfinal")

variables.dc  <- c(variables, "stopdrug_6m_3mFU_chr")



#
# Cohort description of final.hb -----------------------------------------------
#


formula_data_description  <- as.formula(paste("~ ", paste0(c(variables.hb), collapse =  " + "), "| drugclass "))
table_data_description    <- table1(formula_data_description, droplevels = TRUE, render.continuous = "Mean (SD)", data = final.hb)
print(table_data_description)

dim(final.hb)
table(final.hb$drugclass)
length(unique(final.hb$prac))                                                   # 3 last digits of patid indicate practice ID
length(unique(final.hb$patid))
table(final.hb$ethnicity)
round(table(final.hb$ethnicity)/dim(final.hb)[1]*100, 2)

#
# Cohort description of final.wt -----------------------------------------------
#


formula_data_description  <- as.formula(paste("~ ", paste0(c(variables.wt), collapse =  " + "), "| drugclass "))
table_data_description    <- table1(formula_data_description, droplevels = TRUE, render.continuous = "Mean (SD)", data = final.wt)
print(table_data_description)


#
# Cohort description of final.dc -----------------------------------------------
#


formula_data_description  <- as.formula(paste("~ ", paste0(c(variables.dc), collapse =  " + "), "| drugclass "))
table_data_description    <- table1(formula_data_description, droplevels = TRUE, render.continuous = "Mean (SD)", data = final.dc)
print(table_data_description)







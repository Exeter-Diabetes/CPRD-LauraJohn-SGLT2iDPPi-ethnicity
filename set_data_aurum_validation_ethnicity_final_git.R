# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

# Validation of a treatment selection algorithm for optimal choice of SGLT2 and 
# DPP4 inhibitor therapies in people with type 2 diabetes across major UK 
# ethnicity groups 

# Laura M Güdemann, Katherine G Young, Pedro Cardoso, Bilal A Mateen,
# Rury R Holman, Naveed Sattar, Ewan R Pearson, Andrew T Hattersley, 
# Angus G Jones, Beverley M Shields, John M Dennis, 
# on behalf of the MASTERMIND consortium 

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

rm(list=ls())


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Load packages ----------------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

require(tidyverse)
library(rms)
library(tableone)
#library(plyr)
library(ggthemes)
library(forestplot)
library(broom)
library(ggplotify)
library(patchwork)
library(scales)


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Set paths --------------------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --


output_dir <- "C:/Users/lg704/OneDrive - University of Exeter/Studies/SGLT2i_ethnicity_wJohn/study_code/output_Johns_code/" 
data_dir   <- "C:/Users/lg704/OneDrive - University of Exeter/Studies/SGLT2i_ethnicity_wJohn/data/"


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Load data --------------------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

# Dataset defined in data preparation 
load(paste0(data_dir,"final.dataset.sglt2.dpp4.val.Rda"))

# Dataset with additional 12 month data (sensitivity analysis)

load(paste0(data_dir,"final.dataset.sglt2.dpp4.val_25.Rda"))

# GOLD SGLT2-DPP4 model (Dennis et al. 2022, Lancet Digital Health)
load(paste0(data_dir,"m1_hba1cmodel.Rdata"))

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Some data preparation for later sensitivity analysis -------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

# table(final.dataset$ethnicity.backup)

final.dataset$deprivation5 <- as.factor(((as.numeric(final.dataset$deprivation) + 1) %/% 2 ))

# merge new 12 month outcome data into main dataset for sensivitiy analysis 

dim(f.d)
dim(final.dataset)
table(final.dataset$pated%in%f.d$pated)

final.dataset_merged <- merge(final.dataset, f.d, by = "pated")

# check and tidy up so code can run
dim(final.dataset_merged)
dim(final.dataset)

final.dataset <- final.dataset_merged
rm(final.dataset_merged)

# Check if cohorts for 6-month and 12-month are identical so code can be alterted easily for the senstivitiy analysis

table(is.na(final.dataset$posthba1cfinal_12m))
table(is.na(final.dataset$posthba1cfinal))
table(is.na(final.dataset$posthba1cfinal) == is.na(final.dataset$posthba1cfinal_12m)) 

table(is.na(final.dataset$hba1cmonth))
table(is.na(final.dataset$hba1cmonth12))
table(is.na(final.dataset$hba1cmonth) == is.na(final.dataset$hba1cmonth12)) 

str(final.dataset$posthba1cfinal_12m)
str(final.dataset$posthba1cfinal)
str(final.dataset$hba1cmonth)
str(final.dataset$hba1cmonth12)

# --> great! Calculation of this sensitivity analysis is possible with the same study cohort as main analysis 
# (below in the code)




# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Some global settings ---------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

# Global settings
  
  # model.name <- "unadj"
  # OR
  model.name <- "full"
  # OR
  # model.name <- "simple"
  
  
  pdfwidth <- 14
  pdfheight <- 10
  pngwidth <- 3200
  pngheight <- 2400
  pngres <- 200
  row_names <- matrix(c("Predicted SGLT2i glycaemic benefit",paste0(intToUtf8(8805),"5 mmol/mol"), "3-5 mmol/mol","0-3 mmol/mol",
                        "Predicted DPP4i glycaemic benefit","0-3 mmol/mol", paste0(intToUtf8(8805),"3 mmol/mol")))
  tick.hb <- c(-10,-7.5,-5,-2.5,0,2.5,5,7.5, 10)
  tick.hb.resp <- c(-18,-15,-12,-9,-6,-3,0,3)
  tick.dc <- c(0,5,10,15,20,25,30,35,40)
  tick.wt <- c(-4,-3,-2,-1,0)
  B <- 1000

# Add dummy legend
  dummy <- final.dataset %>% 
    sample_n(1000) %>% 
    mutate(drugclass = ifelse(drugclass=="SGLT2","SGLT2-inhibitor","DPP4-inhibitor")) %>%
    ggplot(aes(x = prehba1c, y = prebmi, group=drugclass)) + 
    scale_color_manual(values=c("#4118de","#f1a340")) +
    geom_point(aes(colour=drugclass), size = 1.5) + theme_bw() +
    theme(legend.text = element_text(colour="black", size=rel(1))) + 
    theme(legend.title=element_blank())  + 
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal"
    ) +
    guides(colour = guide_legend(override.aes = list(size=4)))
  
  # Create user-defined function, which extracts legends from ggplots #https://statisticsglobe.com/add-common-legend-to-combined-ggplot2-plots-in-r/
  extract_legend <- function(my_ggp) {
    step1 <- ggplot_gtable(ggplot_build(my_ggp))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
  }
  # Apply user-defined function to extract legend
  shared_legend <- extract_legend(dummy)  

  
  
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Some helpful functions -------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
  
  
# Define HTE treatment diff plot function
hist_plot <- function(data,sx,sy,y) {
  #label for hist
  annotation <- data.frame(
    x = c(sx,sy),
    y = c(y),
    label = c("Favours SGLT2i", "Favours DPP4i")
  )
  #define data
  dat <- data %>% dplyr::select(hba1c_diff) %>% mutate(above=ifelse(hba1c_diff> 0, "Favours DPP4i", "Favours SGLT2i")) 
  c_low <- quantile(dat$hba1c_diff,.001)
  c_upp <- quantile(dat$hba1c_diff,.999)
  c_lowr  <- 2*round(c_low/2)
  c_uppr <- 2*round(c_upp/2)
  c_low <- min(dat$hba1c_diff,.001)
  c_upp <- quantile(dat$hba1c_diff,.999)
  c_lowr  <- 2*round(c_low/2)
  c_uppr <- 2*round(c_upp/2)
  
  #plot
  ggplot(data=dat, aes(x=hba1c_diff,fill=above)) +
    geom_histogram(position="identity", alpha=0.5,color="black",breaks=seq(-15,10,by=1)) +
    geom_vline(aes(xintercept=0), linetype="dashed")+
    labs(title="",x="Predicted HbA1c difference (mmol/mol)", y = "Number of people") +
    #scale_x_continuous(limits=c(c_low,c_upp),breaks=c(seq(c_lowr,c_uppr,by=2))) +
    scale_fill_manual(values=c("#998ec3","#f1a340"))+
    theme_classic() + theme(legend.position = "none")
  #theme(legend.position = c(0.87, 0.97)) + theme(legend.title = element_blank())
  #geom_text(data=annotation, aes(x=x, y=y, label=label),color=c("#f1a340","#998ec3"),size=4, fontface="bold" )
}   

# Function to fit a series of models and output the coefficient(s) of interest with CIs and p-value
hte.model.coefs <- function(x,nmodels) {
  mnumber = c(1:nmodels)
  models <- as.list(1:nmodels)
  nobs <- vector()
  coef <- vector()
  lower <- vector()
  upper <- vector()
  pvalue <- vector()
  data <- x #!!!
  
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(f[[i]]),data=data)
    nobs <- append(nobs,nobs(models[[i]]))
    coef <- append(coef,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower <- append(lower,confint_all[2,1])
    upper <- append(upper,confint_all[2,2])
    pvalue <- append(pvalue,summary(models[[i]])$coefficients[2,4])
  }
  
  datasetname = c(deparse(substitute(x)),deparse(substitute(x)),deparse(substitute(x)))
  x <- data.frame(datasetname,modelname,cbind(nobs,coef,lower,upper,pvalue))
  rownames(x) <- c()
  return(x)
}  

# Function to output HTE by subgroup
hte_plot <- function(data,pred,obs,obslowerci,obsupperci,ymin.ymax) {
  
  #ymin <- min(data$lci); ymax <- max(data$uci);yminr  <- 2*round(ymin/2);  ymaxr <- 2*round(ymax/2)
  ymin  <- -14;  ymax <- 14
  
  ggplot(data=data,aes_string(x=pred,y=obs)) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey60") + geom_hline(yintercept=0, linetype="dashed", color = "grey60") +
    geom_abline(intercept=0,slope=1, color="red", lwd=0.75) + ggtitle("") +
    geom_point(alpha=1) + theme_classic() +
    geom_errorbar(aes_string(ymin=obslowerci, ymax=obsupperci), colour="black", width=.1) +
    ylab("Observed HbA1c difference (mmol/mol)*") + xlab("Predicted HbA1c difference (mmol/mol)") +
    scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(ymin,ymax,by=2))) +
    scale_y_continuous(limits=c(ymin,ymax),oob=rescale_none,breaks=c(seq(ymin,ymax,by=2))) + 
    # scale_x_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    # scale_y_continuous(limits=c(ymin,ymax),breaks=c(seq(yminr,ymaxr,by=2))) +
    theme_base(base_size = 8)  +
    theme(plot.background = element_blank()) + theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8)) +
    theme(text = element_text(size = 10),
          # axis.ticks.x = element_blank(),
          # axis.ticks.y = element_blank(),
          axis.ticks.x = element_line(colour =  "grey50"),
          axis.ticks.y = element_line(colour =  "grey50"),
          axis.line = element_line(colour =  "grey50" ),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
}

# Model and extract HTE coefficients for subgroups
hte_model <- function(data) {
  final <- data
  #Define subsets of interest to calc HTE 
  #overall
  overall <- final
  #sglt2.best
  sglt2.best <- final %>% dplyr::filter(bestdrug=="SGLT2")
  #dpp4.best
  dpp4.best <- final %>% dplyr::filter(bestdrug=="DPP4")
  
  #match trials
  #df1 <- final %>% dplyr::filter(hba1c_diff<= -10) 
  df2 <- final %>% dplyr::filter(hba1c_diff<= -5) 
  df3 <- final %>% dplyr::filter(hba1c_diff<= -3 & hba1c_diff > -5) 
  df4 <- final %>% dplyr::filter(hba1c_diff> -3 & hba1c_diff <= 0) 
  df5 <- final %>% dplyr::filter(hba1c_diff> 0 & hba1c_diff < 3) 
  df6 <- final %>% dplyr::filter(hba1c_diff>= 3) 
  #df7 <- final %>% dplyr::filter(hba1c_diff>= 5) 
  
  #Run HTE models and extract coefs for each subset
  dflist <- list(overall,sglt2.best,dpp4.best,df2,df3,df4,df5,df6)
  res.list <- lapply(dflist, function(df) { 
    hte.model.coefs(df,3) #!!!
  })
  
  res.list <- bind_rows(res.list, .id="column_label")
  res.lab <- c(rep("overall",3),rep("sglt2.best",3),rep("dpp4.best",3),rep("sglt.best5",3),rep("sglt.best3",3),rep("sglt.best0-3",3),
               rep("dpp4.best0-3",3),rep("dpp4.best3",3))
  res.list <- data.frame(res.lab,res.list) %>% dplyr::select(-column_label,-datasetname)
  res.list %>% dplyr::filter(modelname==model.name) #final adjusted list
  
  #Calibration plot HTE
  
  #Unadjusted obs vs pred 
  #Define tenths
  final <- final %>% mutate(hba1c_diff.q = ntile(hba1c_diff, 10))    
  
  #define dataset with predicted values
  t1 <- final %>% 
    group_by(hba1c_diff.q) %>%
    dplyr::summarise(N = length(hba1c_diff),
                     hba1c_diff.pred = mean(hba1c_diff))
  
  #check some patients actually prescribed both drugs in each tenth
  # ddply(final, c("hba1c_diff.q","drugclass"), dplyr::summarise,
  #       N    = length(posthba1c_final),
  #       posthba1c_final.m = mean(posthba1c_final),
  #       se = sd(posthba1c_final)/sqrt((length(posthba1c_final))))
  
  #obs vs pred, by decile of predicted treatment difference
  #For Formula 1-3
  mnumber = c(1:10)
  models  <- as.list(1:10)
  
  hba1c_diff.obs.unadj <- vector()
  lower.unadj <- vector()
  upper.unadj <- vector()
  hba1c_diff.obs.sim <- vector()
  lower.sim <- vector()
  upper.sim <- vector() 
  hba1c_diff.obs.adj <- vector()
  lower.adj <- vector()
  upper.adj <- vector() 
  
  #Unadj
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula1),data=final,subset=hba1c_diff.q==i)
    hba1c_diff.obs.unadj <- append(hba1c_diff.obs.unadj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.unadj <- append(lower.unadj,confint_all[2,1])
    upper.unadj <- append(upper.unadj,confint_all[2,2])
  }
  #Simple 
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula2),data=final,subset=hba1c_diff.q==i)
    hba1c_diff.obs.sim <- append(hba1c_diff.obs.sim,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.sim <- append(lower.sim,confint_all[2,1])
    upper.sim <- append(upper.sim,confint_all[2,2])
  }
  #Full
  for(i in mnumber) {
    models[[i]] <- lm(as.formula(formula3),data=final,subset=hba1c_diff.q==i)
    hba1c_diff.obs.adj <- append(hba1c_diff.obs.adj,models[[i]]$coefficients[2])
    confint_all <- confint(models[[i]], levels=0.95)
    lower.adj <- append(lower.adj,confint_all[2,1])
    upper.adj <- append(upper.adj,confint_all[2,2])
  }
  
  #Final data.frame  
  t1 <- data.frame(t1,cbind(hba1c_diff.obs.unadj,lower.unadj,upper.unadj,
                            hba1c_diff.obs.sim,lower.sim,upper.sim,
                            hba1c_diff.obs.adj,lower.adj,upper.adj))
  
  #unadj
  # plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.unadj,lci=lower.unadj,uci=upper.unadj)
  # hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")
  # #simple adj
  # plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.sim,lci=lower.sim,uci=upper.sim)
  # hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci") 
  #splie adj
  plotdata <- t1 %>% dplyr::mutate(obs=hba1c_diff.obs.adj,lci=lower.adj,uci=upper.adj)
  hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")
  
  #outputs
  hist <- hist_plot(final,-2.5,2.3,1100)
  hte <- hte_plot(plotdata,"hba1c_diff.pred","obs","lci","uci")  
  res.list <- res.list %>% filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall")  %>%
    mutate(order=c(1,5,2,3,4,6,7)) %>% 
    arrange(order) 
  return(list(hist,hte,res.list))
}

#Forest plot for calibration comparison within ethnicity
fp_plot <- function(coef,cim,cip) {
  fp <-
    forestplot(row_names,
               mean = coef,
               lower= cim,
               upper = cip,
               hrzl_lines = gpar(col="#444444"),#lineheight=unit(2,'cm'),
               is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,2)),
               #title="c) TZD - Oedema",
               xticks = tick.hb,
               zero = 0,
               #boxsize=0.1,
               # graphwidth = unit(2,"inches"),
               # lineheight = unit(0.7,"inches"),
               ci.vertices=TRUE,
               col=fpColors(box=c("#66c2a5","#fc8d62","#8da0cb"), lines=c("#66c2a5","#fc8d62","#8da0cb"), zero = "gray50"),
               #lty.ci = c(1,2,3,4),
               xlab="Observed HbA1c difference (mmol/mol)* [negative favours SGLT2i]",cex=1,
               new_page = TRUE,
               fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
               boxsize = .2, # We set the box size to better visualize the type
               #line.margin = .1, # We need to add this to avoid crowding
               txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),ticks  = gpar(cex = 1)),
               #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
               # legend_args = fpLegend(
               #   pos = list("topright"),
               #   title = "Group",
               #   r = unit(.1, "snpc"),
               #   gp = gpar(col = "#CCCCCC", lwd = 1.5))
               legend = c("Uncalibrated","Recalibrated"),#,
               legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
               #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
               #xlog = TRUE
    )
  return(fp)
}

# Forest plot for calibration comparison by ethnicity
fp_plot.eth <- function(coef,cim,cip) {
  fp <-
    forestplot(row_names,
               mean = coef,
               lower= cim,
               upper = cip,
               hrzl_lines = gpar(col="#444444"),#lineheight=unit(2,'cm'),
               is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,2)),
               #title="c) TZD - Oedema",
               xticks = tick.hb,
               zero = 0,
               #boxsize=0.1,
               # graphwidth = unit(2,"inches"),
               # lineheight = unit(0.7,"inches"),
               ci.vertices=TRUE,
               col=fpColors(box=c("#66c2a5","#fc8d62","#8da0cb","#752302"), lines=c("#66c2a5","#fc8d62","#8da0cb","#752302"), zero = "gray50"),
               #lty.ci = c(1,2,3,4),
               xlab="Observed HbA1c difference (mmol/mol)* [negative favours SGLT2i]",cex=1,
               new_page = TRUE,
               fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI, fpDrawPointCI),
               boxsize = .1, # We set the box size to better visualize the type
               #line.margin = .1, # We need to add this to avoid crowding
               txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),ticks  = gpar(cex = 1)),
               #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
               # legend_args = fpLegend(
               #   pos = list("topright"),
               #   title = "Group",
               #   r = unit(.1, "snpc"),
               #   gp = gpar(col = "#CCCCCC", lwd = 1.5))
               legend = c("White (77.8%)", "South Asian (13.7%)", "Black (4.2%)", "Mixed/Other (2.4%)"),#,
               legend_args = fpLegend(pos = list(x=.850, y=0.95))#, 
               #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
               #xlog = TRUE
    )
  return(fp)
}



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Define each of the HbA1c, weight, discontinuation cohorts --------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --


## Overall cohort --------------------------------------------------------------

# Set drugline as factor

final.dataset <- final.dataset %>% mutate(drugclass=as.factor(drugclass),
                                          drugline=as.factor(drugline))

# Collapse ethnicity
final.dataset <- final.dataset %>% mutate(ethnicity.backup = ethnicity,
                                ethnicity=fct_collapse(ethnicity,mixed.other=c("Mixed","Other")))

## HbA1c cohort ----------------------------------------------------------------

final.hb <- final.dataset %>%
  #select(patid, pated, posthba1cfinal, ethnicity, ethnicity16, ethnicity.backup, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx, hba1cmonth, deprivation5) %>% # only for sensitivity analysis
  select(patid, pated, posthba1cfinal, posthba1cfinal_12m, ethnicity, ethnicity16, ethnicity.backup, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx, hba1cmonth, hba1cmonth12) %>%
  dplyr::rename("prealtlog"="prealt",
                "posthba1c_final"="posthba1cfinal",
                "prehba1cmmol"="prehba1c",
                "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog)) %>% 
  drop_na() #%>% 

final.hb.backup <- final.hb

# Proportion of patients by ethnicity 
prop.tab.hb    <- round(prop.table(table(final.hb$ethnicity))*100,1)
prop.tab.hb_sa <- round(prop.table(table(final.hb$ethnicity.backup))*100,1)


## Weight cohort ---------------------------------------------------------------

final.wt <- final.dataset %>%
  select(patid, pated, preweight, postweightfinal, ethnicity, ethnicity16, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx) %>%
  dplyr::rename("prealtlog"="prealt",
                "prehba1cmmol"="prehba1c",
                "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog),
         wtchange=postweightfinal-preweight) %>% 
  drop_na() #%>% 

## Discontinuation cohort ------------------------------------------------------

final.dc <- final.dataset %>%
  select(patid, pated, stopdrug_6m_3mFU, ethnicity, ethnicity16, ethnicity.backup, prehba1c, drugclass, drugline, ncurrtx, preegfr, prealt, prebmi, agetx) %>%
  dplyr::rename("prealtlog"="prealt",
                "prehba1cmmol"="prehba1c",
                "egfr_ckdepi"="preegfr") %>%
  mutate(prealtlog = log(prealtlog)) %>% 
  drop_na() #%>% 




# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Define model formula for the validation -------------------------------------- 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Formulas for main analysis --------------------------------------------------

formula1 <- "posthba1c_final~drugclass"
##formula2 <- "posthba1c_final~factor(drugclass)+prehba1cmmol+ncurrtx+drugline+rcs(hba1cmonth,3)+egfr_ckdepi+prealtlog" # not used
formula2 <- "posthba1c_final~drugclass+ncurrtx+drugline"
formula3 <- "posthba1c_final~drugclass+rcs(prehba1cmmol,3)+ncurrtx+drugline+rcs(hba1cmonth,3)+rcs(egfr_ckdepi,3)+rcs(prealtlog,3)+rcs(agetx,3)+rcs(prebmi,3)" # used in paper

## Formulas for senstivity analysis --------------------------------------------

#formula1 <- "posthba1c_final~drugclass" # same as before 
#formula2 <- "posthba1c_final~drugclass+ncurrtx+drugline+deprivation5"
#formula3 <- "posthba1c_final~drugclass+rcs(prehba1cmmol,3)+ncurrtx+drugline+deprivation5+rcs(hba1cmonth,3)+rcs(egfr_ckdepi,3)+rcs(prealtlog,3)+rcs(agetx,3)+rcs(prebmi,3)" # used in paper 


modelname <- c("unadj","simple","full")
f <- as.list(c(formula1,formula2,formula3))


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Predict HbA1c outcomes from the model ----------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

# Predict outcomes
  final.hb <- final.hb %>% mutate(drug=drugclass) %>% mutate(drugclass="DPP4")
  final.hb$DPP4.pred.lm <- predict(m1,final.hb)
  final.hb <- final.hb %>% mutate(drugclass="SGLT2")
  final.hb$SGLT2.pred.lm <- predict(m1,final.hb)                                
  final.hb <- final.hb %>% mutate(drugclass=drug) %>%
    mutate(hba1c_diff = SGLT2.pred.lm-DPP4.pred.lm,
           bestdrug=ifelse(hba1c_diff<=0,"SGLT2","DPP4"),
           drugclass==drug)

  final.hb <- final.hb %>% 
    mutate(DPP4.pred.lm.uncal=DPP4.pred.lm,
           SGLT2.pred.lm.uncal=SGLT2.pred.lm,
           hba1c_diff.uncal=hba1c_diff,
           bestdrug.uncal=bestdrug) %>%
    select(-DPP4.pred.lm,-SGLT2.pred.lm,-hba1c_diff,-bestdrug)
  
# brief summary of predicted treatment difference
  describe(final.hb$hba1c_diff)                                                 # no hba1c_diff for SA!!!
  hist(final.hb$hba1c_diff,breaks=50); abline(v = 0, col="black", lwd=3, lty=2)
  table(final.hb$bestdrug)

# #Run HTE models and extract coefs
#   hte.output.overall <- hte_model(final.hb)
# 
# #outputs
#   hist.overall <- hte.output.overall[1]
#   hte.overall <- hte.output.overall[2]  
#   res.list.overall <- data.frame(hte.output.overall[3]) %>% filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall") 
#   hist.overall  
#   hte.overall
#   res.list.overall
# 
#   val <- res.list.overall %>% filter(res.lab != "sglt2.best" & res.lab != "dpp4.best")
# 
# #Plot
# 
#   x <- rep(NA,4)
#   coef = rbind(x, data.frame(cbind(val[1:5,4])))
#   cim = rbind(x, data.frame(cbind(val[1:5,5])))
#   cip = rbind(x, data.frame(cbind(val[1:5,6])))
#   
#   insertrow <- function(existingDF, newrow, r) {
#     existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
#     existingDF[r,] <- newrow
#     existingDF
#   }
#   
#   coef <- data.matrix(insertrow(coef, NA, 5))
#   cim <- data.matrix(insertrow(cim, NA, 5))
#   cip <- data.matrix(insertrow(cip, NA, 5))
#   
#   
#   row_names.n <- matrix(c("Predicted SGLT2i benefit",paste0(intToUtf8(8805),">5 mmol/mol (n=42,057)"), "3-5 mmol/mol (n=19,808)","0-3 mmol/mol (n=21,171)",
#                         "Predicted DPP4i benefit","0-3 mmol/mol (n=10,244)",paste0(intToUtf8(8805),"3 mmol/mol (n=4,207)")))
#   fp <-
#     forestplot(row_names.n,
#                mean = coef,
#                lower= cim,
#                upper = cip,
#                hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
#                is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,2)),
#                #title="c) TZD - Oedema",
#                xticks = tick.hb,
#                zero = 0,
#                #boxsize=0.1,
#                # graphwidth = unit(2,"inches"),
#                # lineheight = unit(0.7,"inches"),
#                ci.vertices=TRUE,
#                col=fpColors(box=c("#66c2a5"), lines=c("#66c2a5"), zero = "gray50"),
#                #lty.ci = c(1,2,3,4),
#                xlab="Average HbA1c difference (mmol/mol; negative favours SGLT2i)",cex=1,
#                new_page = TRUE,
#                #fn.ci_norm = c(fpDrawNormalCI),
#                boxsize = .1, # We set the box size to better visualize the type
#                #line.margin = .1, # We need to add this to avoid crowding
#                txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),ticks  = gpar(cex = 1)),
#                #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
#                # legend_args = fpLegend(
#                #   pos = list("topright"),
#                #   title = "Group",
#                #   r = unit(.1, "snpc"),
#                #   gp = gpar(col = "#CCCCCC", lwd = 1.5))
#                #legend = c("White (77.8%)", "Asian (13.7%)", "Black (4.2%)", "Mixed/Other (2.4%)"),#,
#                #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
#                #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
#                #xlog = TRUE
#     )
#   
#   fp.overall <- fp
#   fp.overall
# 
#   final.hb <- final.hb %>% 
#     mutate(DPP4.pred.lm.uncal=DPP4.pred.lm,
#            SGLT2.pred.lm.uncal=SGLT2.pred.lm,
#            hba1c_diff.uncal=hba1c_diff,
#            bestdrug.uncal=bestdrug) %>%
#     select(-DPP4.pred.lm,-SGLT2.pred.lm,-hba1c_diff,-bestdrug)
  

  
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Perform the model update for the HbA1c model ---------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
  
# global settings
p.value     <- 0.05
ncolx       <- length(m1$coefficients)-1
sample_frac <- 1

set.seed(8731)

# Predict outcome on therapy received
final.hb$pred <- predict(m1,final.hb)

# Testing function
closedtest <- function(cohort, dataset, observed, predicted, p.value){
  
  #Original model
  
  #Residuals
  resid <- observed-predicted
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #original log-likelihood
  n <- length(resid)
  logLik.original <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #Update intercept
  
  #Model with updated intercept
  m <- lm(observed-predicted~1,data=dataset)
  
  #Extract coefficient
  m1.intercept <- cbind(m$coefficients[1],confint(m)[1],confint(m)[2])
  
  #Residuals (actual - predicted)
  resid <- residuals(m)
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #log-likelihood
  n <- length(resid)
  logLik.intercept <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #Update slope & intercept
  
  #Model with updated slope &  intercept
  m <- lm(observed~predicted,data=dataset)
  
  #Extract coefficient
  m2.intercept <- cbind(m$coefficients[1],confint(m)[1,1],confint(m)[1,2])
  m2.slope <- cbind(m$coefficients[2],confint(m)[2,1],confint(m)[2,2])
  
  #Residuals (actual - predicted)
  resid <- residuals(m)
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #log-likelihood
  n <- length(resid)
  logLik.recal <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #test significance
  
  #1. Test recal in the large against the original model (no extra coeffs estimated)
  #2. If 2. is significant, test full recal against the recal in the large model using p df + 1 (1 extra coef estimated)
  #3. If 3. is significant, select 3. as final model, if not select 2. If neither, select 1 
  
  #ll diff
  dev_intercept <- -2*logLik.original + 2*logLik.intercept
  dev_recal <- -2*logLik.intercept + 2*logLik.recal
  
  #Diff in ll
  ncolx <- ncolx
  test1 <- (1-pchisq(dev_intercept, ncolx)) < p.value
  test2 <- (1-pchisq(dev_recal, ncolx+1)) < p.value
  
  #p.value
  p1 <- (1-pchisq(dev_intercept, ncolx))
  p2 <- (1-pchisq(dev_recal, ncolx+1))
  
  #Which model is chosen
  test_intercept <- 1 * (!test1)
  test_recal <- 2 * ((!test1)&(!test2))
  
  index_test <- (test_intercept + test_recal)
  
  res <- data.frame(cohort=c(cohort,cohort,cohort),
                    n=c(nrow(dataset),nrow(dataset),nrow(dataset)),
                    model=c("Original","Updated intercept","Recalibrated"),
                    loglikelihood=c(logLik.original,logLik.intercept,logLik.recal),
                    intercept=c(NA,
                                m1.intercept[1],
                                m1.intercept[2]),
                    intercept.w.ci=c(NA,
                                     paste0(round(m1.intercept[1],2), " (",paste0(round(m1.intercept[2],2),", ",paste0(round(m1.intercept[3],2)),")")),
                                     paste0(round(m2.intercept[1],2), " (",paste0(round(m2.intercept[2],2),", ",paste0(round(m2.intercept[3],2)),")"))),
                    
                    slope=c(NA,NA,m2.slope[1]),
                    slope.w.ci=c(NA,
                                 NA,
                                 paste0(round(m2.slope[1],2), " (",paste0(round(m2.slope[2],2),", ",paste0(round(m2.slope[3],2)),")"))),
                    p.value=c(NA,
                              round(p1,5),
                              round(p2,5)),
                    model.selected=c(ifelse(test1==FALSE & test2==FALSE,"Yes","No"),
                                     ifelse(test1==TRUE & test2==FALSE,"Yes","No"),
                                     ifelse(test2==TRUE,"Yes","No"))
  )
  return(res)
}



## Application of testing procedure --------------------------------------------

### DPP4i.White subgroup -------------------------------------------------------

# Setup
cohort <- "DPP4i.White"

dataset <- final.hb %>% 
  filter(ethnicity=="White" & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
DPP4i.White <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.White

### DPP4i.Asian subgroup -------------------------------------------------------

# Setup
cohort  <- "DPP4i.Asian"

dataset <- final.hb %>% 
  filter(ethnicity=="South Asian" & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
DPP4i.Asian <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.Asian


### DPP4i.Black subgroup -------------------------------------------------------

# Setup
cohort  <- "DPP4i.Black"

dataset <- final.hb %>% 
  filter(ethnicity=="Black" & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
DPP4i.Black <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.Black

### DPP4i.Mixed.other subgroup -------------------------------------------------

# Setup
cohort  <- "DPP4i.Mixed.Other"

dataset <- final.hb %>% 
  filter((ethnicity=="mixed.other") & drugclass=="DPP4") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
DPP4i.Mixed.Other <- closedtest(cohort,dataset,observed,predicted,p.value)
DPP4i.Mixed.Other


### SGLT2i.White subgroup ------------------------------------------------------

# Setup
cohort  <- "SGLT2i.White"

dataset <- final.hb %>% 
  filter(ethnicity=="White" & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
SGLT2i.White <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2i.White


### SGLT2i.Asian subgroup ------------------------------------------------------

# Setup
cohort  <- "SGLT2i.Asian"

dataset <- final.hb %>% 
  filter(ethnicity=="South Asian" & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
SGLT2.Asian <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2.Asian


### SGLT2.Black subgroup -------------------------------------------------------

# Setup
cohort  <- "SGLT2i.Black"

dataset <- final.hb %>% 
  filter(ethnicity=="Black" & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
SGLT2i.Black <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2i.Black


### SGLT2i.Mixed.other subgroup ------------------------------------------------

# Setup
cohort  <- "SGLT2i.Mixed.Other"

dataset <- final.hb %>% 
  filter((ethnicity=="mixed.other") & drugclass=="SGLT2") %>%
  sample_frac(sample_frac)

observed  <- dataset$posthba1c_final
predicted <- dataset$pred

# Test
SGLT2i.Mixed.Other <- closedtest(cohort,dataset,observed,predicted,p.value)
SGLT2i.Mixed.Other


## Summarize final results of the testing procedure ----------------------------

#final.hb results
closedtest.final.hb <- rbind(DPP4i.White,SGLT2i.White,DPP4i.Asian,SGLT2.Asian,DPP4i.Black,SGLT2i.Black,DPP4i.Mixed.Other,SGLT2i.Mixed.Other)
closedtest.final.hb 

ctfm <- closedtest.final.hb %>% filter(model.selected=="Yes")

ctfm.save <- closedtest.final.hb %>% filter(model.selected=="Yes") %>% 
  select("Therapy/Ethnicity arm"=cohort,"Number of people"=n,"Model selected" = model, "Updated intercept (if required)" = intercept.w.ci)
write.csv(ctfm,file=paste0(output_dir,"closedloopmodelupdate.csv"))

final.hb <- final.hb %>% mutate(DPP4.pred.lm.recal=ifelse(ethnicity == "White",DPP4.pred.lm.uncal+ctfm$intercept[1],
                                                    ifelse(ethnicity == "South Asian",DPP4.pred.lm.uncal+ctfm$intercept[3],
                                                           ifelse(ethnicity == "Black",DPP4.pred.lm.uncal+ctfm$intercept[5],
                                                                  ifelse(ethnicity == "mixed.other",DPP4.pred.lm.uncal+ctfm$intercept[7],NA
                                                                  )))),
                          SGLT2.pred.lm.recal=ifelse(ethnicity == "White",SGLT2.pred.lm.uncal+ctfm$intercept[2],SGLT2.pred.lm.uncal),
                          hba1c_diff.recal = SGLT2.pred.lm.recal-DPP4.pred.lm.recal,
                          bestdrug.recal=ifelse(hba1c_diff.recal<=0,"SGLT2","DPP4"))

table(final.hb$bestdrug.uncal,final.hb$bestdrug.recal,final.hb$ethnicity)                                                  





# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Evaluate model performance for each ethnicity group --------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

## Calculation of RMSE of prediction error -------------------------------------

# to get a clearer picture of model performance (per ethnicity group)

# recalibrated model predicted outcomes for treatments initiated 
final.hb$model_predicted_outcomes                                <- final.hb$DPP4.pred.lm.recal
final.hb$model_predicted_outcomes[final.hb$drugclass == "SGLT2"] <- final.hb$SGLT2.pred.lm.recal[final.hb$drugclass == "SGLT2"]

final.hb %>%
  group_by(ethnicity) %>%
  summarise(
    rmse = sqrt(mean((posthba1c_final - model_predicted_outcomes)^2, na.rm = TRUE))
  )


## Plot observed vs predicted HbA1c values -------------------------------------

# data subsets 

final.hb_White      <- final.hb %>% filter(ethnicity == "White")
final.hb_SouthAsian <- final.hb %>% filter(ethnicity == "South Asian")
final.hb_Black      <- final.hb %>% filter(ethnicity == "Black")
final.hb_MixedOther <- final.hb %>% filter(ethnicity == "mixed.other")


# plots for ethnicity groups 

lw_White      <- loess(posthba1c_final ~ model_predicted_outcomes, data = final.hb_White)
lw_SouthAsian <- loess(posthba1c_final ~ model_predicted_outcomes, data = final.hb_SouthAsian)
lw_Black      <- loess(posthba1c_final ~ model_predicted_outcomes, data = final.hb_Black)
lw_MixedOther <- loess(posthba1c_final ~ model_predicted_outcomes, data = final.hb_MixedOther)


# save PDF plot
grDevices::cairo_pdf(paste0(output_dir,"observed_vs_predicted_HbA1c.pdf"),width=8,height=8)

par(mfrow = c(2, 2))

# plot for White ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_White, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)",
     main = "White")
lines(final.hb_White$model_predicted_outcomes, lw_White$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)

# plot for South Asian ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_SouthAsian, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)",
     main = "South Asian")
lines(final.hb_SouthAsian$model_predicted_outcomes, lw_SouthAsian$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)

# plot for Black ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_Black, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)",
     main = "Black")
lines(final.hb_Black$model_predicted_outcomes, lw_Black$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)

# plot for Mixed/Other ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_MixedOther, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)", 
     main = "Mixed/Other")
lines(final.hb_MixedOther$model_predicted_outcomes, lw_MixedOther$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)


dev.off()


# Save png plot
png(paste0(output_dir,"observed_vs_predicted_HbA1c.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
par(mfrow = c(2, 2))

# plot for White ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_White, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)",
     main = "White")
lines(final.hb_White$model_predicted_outcomes, lw_White$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)

# plot for South Asian ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_SouthAsian, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)",
     main = "South Asian")
lines(final.hb_SouthAsian$model_predicted_outcomes, lw_SouthAsian$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)

# plot for Black ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_Black, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)",
     main = "Black")
lines(final.hb_Black$model_predicted_outcomes, lw_Black$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)

# plot for Mixed/Other ethnicity
plot(posthba1c_final ~ model_predicted_outcomes, data = final.hb_MixedOther, pch = 19, cex = 0.1,
     ylim = c(30, 130), xlim = c(30, 130),
     ylab = "Observed 6 month HbA1c (mmol/mol)",
     xlab = "Observed 6 month HbA1c (mmol/mol)", 
     main = "Mixed/Other")
lines(final.hb_MixedOther$model_predicted_outcomes, lw_MixedOther$fitted, col = "yellow", lwd = 3)
abline(a = 0, b = 1, col = "darkgreen", lty = 2, lwd = 2)

dev.off()

par(mfrow = c(1, 1))

final.hb$model_predicted_outcomes <- NULL


## Asian -----------------------------------------------------------------------


#final.hb$posthba1cfinal <- final.hb$posthba1cfinal_12m                          # only for sensitivity analysis
#final.hb$hba1cmonth     <- final.hb$hba1cmonth12                                # only for sensitivity analysis


final <- final.hb %>% filter(ethnicity == "South Asian")

# GOLD model (same as code above)
final$hba1c_diff  <- final$hba1c_diff.uncal                                     # GOLD model
final$bestdrug    <- final$bestdrug.uncal                                       # GOLD model
overall.uncal     <- hte_model(final)

# recalibrated intercept
final$hba1c_diff  <- final$hba1c_diff.recal                                     # recalibrated
final$bestdrug    <- final$bestdrug.recal 
overall.recal     <- hte_model(final)

# outputs

# Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))

hist.asian.both  <- wrap_elements(p1) + wrap_elements(p2) 
hist.asian.uncal <- wrap_elements(p1) 
hist.asian.recal <- wrap_elements(p2) 
hist.asian.both

# Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))
#p3 <- grid2grob(print(overall.recallp[2]))

cal.asian.both  <- wrap_elements(p1) + wrap_elements(p2) 
cal.asian.uncal <- wrap_elements(p1) 
cal.asian.recal <- wrap_elements(p2) 
cal.asian.both

# HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x   <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))         # ,val[11:15,4]))))
cim  = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))         # ,val[11:15,5]))))
cip  = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))         # ,val[11:15,6]))))

coef <- rbind(coef[1:4,],c(NA), coef[5:6,])
cim  <- rbind(cim[1:4,],c(NA), cim[5:6,])
cip  <- rbind(cip[1:4,],c(NA), cip[5:6,])

fp.asian.recal  <- fp_plot(coef,cim,cip)
fp.asian.recal 
val.asian.recal <- val



## White -----------------------------------------------------------------------

final <- final.hb %>% filter(ethnicity == "White")

# uncal model (same as code above)
final$hba1c_diff <- final$hba1c_diff.uncal                                      # uncal model
final$bestdrug   <- final$bestdrug.uncal                                        # GOLD model
overall.uncal    <- hte_model(final)

# recalibrated intercept
final$hba1c_diff <- final$hba1c_diff.recal                                      # recalibrated
final$bestdrug   <- final$bestdrug.recal 
overall.recal    <- hte_model(final)

# outputs

# Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))

hist.white.both  <- wrap_elements(p1) + wrap_elements(p2) 
hist.white.uncal <- wrap_elements(p1) 
hist.white.recal <- wrap_elements(p2) 
hist.white.both

# Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))

cal.white.both  <- wrap_elements(p1) + wrap_elements(p2) 
cal.white.uncal <- wrap_elements(p1)
cal.white.recal <- wrap_elements(p2) 
cal.white.both

# HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x   <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))
cim  = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))
cip  = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim  <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip  <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.white.recal  <- fp_plot(coef,cim,cip)
fp.white.recal 
val.white.recal <- val


## Black -----------------------------------------------------------------------

final <- final.hb %>% filter(ethnicity == "Black")

# uncal model (same as code above)
final$hba1c_diff <- final$hba1c_diff.uncal                                      # uncal model
final$bestdrug   <- final$bestdrug.uncal                                        # GOLD model
overall.uncal    <- hte_model(final)

# recalibrated intercept
final$hba1c_diff <- final$hba1c_diff.recal                                      # recalibrated
final$bestdrug   <- final$bestdrug.recal 
overall.recal    <- hte_model(final)

# outputs

# Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))

hist.black.both  <- wrap_elements(p1) + wrap_elements(p2) 
hist.black.uncal <- wrap_elements(p1) 
hist.black.recal <- wrap_elements(p2) 
hist.black.both

# Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))

cal.black.both  <- wrap_elements(p1) + wrap_elements(p2) 
cal.black.uncal <- wrap_elements(p1) 
cal.black.recal <- wrap_elements(p2) 
cal.black.both

# HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x   <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))#,val[11:15,4]))))
cim  = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))#,val[11:15,5]))))
cip  = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))#,val[11:15,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim  <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip  <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.black.recal  <- fp_plot(coef,cim,cip)
fp.black.recal 
val.black.recal <- val


# Mixed.Other ------------------------------------------------------------------

final <- final.hb %>% filter(ethnicity == "mixed.other")

# uncal model (same as code above)
final$hba1c_diff <- final$hba1c_diff.uncal                                      # uncal model
final$bestdrug   <- final$bestdrug.uncal                                        # GOLD model
overall.uncal    <- hte_model(final)

# recalibrated intercept
final$hba1c_diff <- final$hba1c_diff.recal                                      # recalibrated
final$bestdrug   <- final$bestdrug.recal 
overall.recal    <- hte_model(final)

# outputs

# Treatment diff histograms
p1 <- grid2grob(print(overall.uncal[1]))
p2 <- grid2grob(print(overall.recal[1]))

hist.mixed.other.both  <- wrap_elements(p1) + wrap_elements(p2) 
hist.mixed.other.uncal <- wrap_elements(p1)
hist.mixed.other.recal <- wrap_elements(p2) 
hist.mixed.other.both

# Calibration plots
p1 <- grid2grob(print(overall.uncal[2]))
p2 <- grid2grob(print(overall.recal[2]))

cal.mixed.other.both  <- wrap_elements(p1) + wrap_elements(p2) 
cal.mixed.other.uncal <- wrap_elements(p1) 
cal.mixed.other.recal <- wrap_elements(p2) 
cal.mixed.other.both

# HTE subgroups
res.list.overall.uncal <- data.frame(overall.uncal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Uncalibrated")
res.list.overall.recal <- data.frame(overall.recal[3]) %>% 
  filter(modelname==model.name & res.lab != "sglt.best10" & res.lab != "dpp4.best5"  & res.lab != "overall" & res.lab != "sglt2.best" & res.lab != "dpp4.best") %>%
  mutate(analysis="Recalibrated intercept")
res.list.overall.comparison <- rbind(res.list.overall.uncal,res.list.overall.recal)#,res.list.overall.recallp)
res.list.overall.comparison

val <- res.list.overall.comparison

x   <- rep(NA,2)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4]))))
cim  = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5]))))
cip  = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6]))))

coef <- rbind(coef[1:4,],c(NA),coef[5:6,])
cim  <- rbind(cim[1:4,],c(NA),cim[5:6,])
cip  <- rbind(cip[1:4,],c(NA),cip[5:6,])

fp.mixed.other.recal  <- fp_plot(coef,cim,cip)
fp.mixed.other.recal 
val.mixed.other.recal <- val


# Mixed (sensitivity analysis) -------------------------------------------------

#final <- final.hb %>% filter(ethnicity.backup == "Mixed")


# uncal model (same as code above)
#final$hba1c_diff <- final$hba1c_diff.uncal                                      # uncal model
#final$bestdrug   <- final$bestdrug.uncal                                        # GOLD model
#overall.uncal    <- hte_model(final)

# recalibrated intercept
#final$hba1c_diff <- final$hba1c_diff.recal                                      # recalibrated
#final$bestdrug   <- final$bestdrug.recal 
#overall.recal    <- hte_model(final)

# outputs

# Treatment diff histograms
#p1 <- grid2grob(print(overall.uncal[1]))
#p2 <- grid2grob(print(overall.recal[1]))

#hist.mixed.both  <- wrap_elements(p1) + wrap_elements(p2) 
#hist.mixed.uncal <- wrap_elements(p1)
#hist.mixed.recal <- wrap_elements(p2) 
#hist.mixed.both

# Calibration plots
#p1 <- grid2grob(print(overall.uncal[2]))
#p2 <- grid2grob(print(overall.recal[2]))

#cal.mixed.both  <- wrap_elements(p1) + wrap_elements(p2) 
#cal.mixed.uncal <- wrap_elements(p1) 
#cal.mixed.recal <- wrap_elements(p2) 
#cal.mixed.both


# Other (senstivity analysis) --------------------------------------------

#final <- final.hb %>% filter(ethnicity.backup == "Other")

# uncal model (same as code above)
#final$hba1c_diff <- final$hba1c_diff.uncal                                      # uncal model
#final$bestdrug   <- final$bestdrug.uncal                                        # GOLD model
#overall.uncal    <- hte_model(final)

# recalibrated intercept
#final$hba1c_diff <- final$hba1c_diff.recal                                      # recalibrated
#final$bestdrug   <- final$bestdrug.recal 
#overall.recal    <- hte_model(final)

# outputs

# Treatment diff histograms
#p1 <- grid2grob(print(overall.uncal[1]))
#p2 <- grid2grob(print(overall.recal[1]))

#hist.other.both  <- wrap_elements(p1) + wrap_elements(p2) 
#hist.other.uncal <- wrap_elements(p1)
#hist.other.recal <- wrap_elements(p2) 
#hist.other.both

# Calibration plots
#p1 <- grid2grob(print(overall.uncal[2]))
#p2 <- grid2grob(print(overall.recal[2]))

#cal.other.both  <- wrap_elements(p1) + wrap_elements(p2) 
#cal.other.uncal <- wrap_elements(p1) 
#cal.other.recal <- wrap_elements(p2) 
#cal.other.both



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Output fp.eth_recal plot (not included in paper) -----------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --


# Uncalibrated

val.uncal <- rbind(val.white.recal[1:5,],val.asian.recal[1:5,],val.black.recal[1:5,],val.mixed.other.recal[1:5,])

val <- val.uncal
x   <- rep(NA,4)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4],val[11:15,4],val[16:20,4]))))
cim  = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5],val[11:15,5],val[16:20,5]))))
cip  = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6],val[11:15,6],val[16:20,6]))))

coef <- rbind(coef[1:4,],c(NA,NA,NA,NA),coef[5:6,])
cim  <- rbind(cim[1:4,],c(NA,NA,NA,NA),cim[5:6,])
cip  <- rbind(cip[1:4,],c(NA,NA,NA,NA),cip[5:6,])

fp.eth.uncal <-fp_plot.eth(coef,cim,cip)

# Intercept
val.int <- rbind(val.white.recal[6:10,],val.asian.recal[6:10,],val.black.recal[6:10,],val.mixed.other.recal[6:10,])

val <- val.int
x   <- rep(NA,4)
coef = data.matrix(rbind(x, data.frame(cbind(val[1:5,4],val[6:10,4],val[11:15,4],val[16:20,4]))))
cim  = data.matrix(rbind(x, data.frame(cbind(val[1:5,5],val[6:10,5],val[11:15,5],val[16:20,5]))))
cip  = data.matrix(rbind(x, data.frame(cbind(val[1:5,6],val[6:10,6],val[11:15,6],val[16:20,6]))))

coef <- rbind(coef[1:4,],c(NA,NA,NA,NA),coef[5:6,])
cim  <- rbind(cim[1:4,],c(NA,NA,NA,NA),cim[5:6,])
cip  <- rbind(cip[1:4,],c(NA,NA,NA,NA),cip[5:6,])

fp.eth.int <-fp_plot.eth(coef,cim,cip)

#Plot together
p1 <- grid2grob(print(fp.eth.uncal))
p2 <- grid2grob(print(fp.eth.int))

fp.eth <- wrap_elements(p1) + wrap_elements(p2) + #/ wrap_elements(p3) +
  plot_annotation("Left=uncalibrated, Right=recalibrated intercept")


# grDevices::cairo_pdf(paste0(output_dir,"fp.eth_recal_both.pdf"),width=20,height=8)
# fp.eth 
# dev.off()
# 
# png(paste0(output_dir,"fp.eth_recal_both.png"),width=pngwidth,height=1250,res=pngres,restoreConsole=TRUE)
# fp.eth 
# dev.off()
# 
# grDevices::cairo_pdf(paste0(output_dir,"fp.eth_uncal.pdf"),width=10,height=8)
# wrap_elements(p1)
# dev.off()
# 
# png(paste0(output_dir,"fp.eth_uncal.png"),width=1800,height=1250,res=pngres,restoreConsole=TRUE)
# wrap_elements(p1)
# dev.off()

grDevices::cairo_pdf(paste0(output_dir,"fp.eth_recal.pdf"),width=10,height=8)
wrap_elements(p2)
dev.off()

png(paste0(output_dir,"fp.eth_recal.png"),width=1800,height=1250,res=pngres,restoreConsole=TRUE)
wrap_elements(p2)
dev.off()




# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Save histogram plot and calibration plot (hist.eth_recal & cal.eth_recal) ----
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --

thm <- theme(plot.title = element_text(face = 1, size = 12))


## Save hist.eth_recal ---------------------------------------------------------

# recalibrated
hist.white.recal.p <- wrap_elements(hist.white.recal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
hist.asian.recal.p <- wrap_elements(hist.asian.recal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
hist.black.recal.p <- wrap_elements(hist.black.recal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
hist.mixed.other.recal.p <- wrap_elements(hist.mixed.other.recal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))

hist.eth <- 
  hist.white.recal.p +
  hist.asian.recal.p + 
  hist.black.recal.p + 
  hist.mixed.other.recal.p +
  plot_layout(ncol = 2) 
hist.eth

grDevices::cairo_pdf(paste0(output_dir,"hist.eth_recal.pdf"),width=8,height=8)
hist.eth
dev.off()

png(paste0(output_dir,"hist.eth_recal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
hist.eth
dev.off()

table(final.hb$ethnicity,final.hb$bestdrug.recal)
round(prop.table(table(final.hb$ethnicity,final.hb$bestdrug.recal),1)*100,1)

# #uncalibrated
# hist.white.uncal.p <- wrap_elements(hist.white.uncal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
# hist.asian.uncal.p <- wrap_elements(hist.asian.uncal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
# hist.black.uncal.p <- wrap_elements(hist.black.uncal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
# hist.mixed.other.uncal.p <- wrap_elements(hist.mixed.other.uncal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))
# 
# hist.eth <- 
#   hist.white.uncal.p +
#   hist.asian.uncal.p + 
#   hist.black.uncal.p + 
#   hist.mixed.other.uncal.p +
#   plot_layout(ncol = 2) 
# 
# grDevices::cairo_pdf(paste0(output_dir,"hist.eth_uncal.pdf"),width=8,height=8)
# hist.eth
# dev.off()
# 
# png(paste0(output_dir,"hist.eth_uncal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
# hist.eth
# dev.off()



## Save cal.eth_recal ----------------------------------------------------------

# Recalibrated
cal.white.recal.p       <- wrap_elements(cal.white.recal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
cal.asian.recal.p       <- wrap_elements(cal.asian.recal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
cal.black.recal.p       <- wrap_elements(cal.black.recal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
cal.mixed.other.recal.p <- wrap_elements(cal.mixed.other.recal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))

cal.eth <- 
  cal.white.recal.p +
  cal.asian.recal.p + 
  cal.black.recal.p + 
  cal.mixed.other.recal.p +
  plot_layout(ncol = 2) 

grDevices::cairo_pdf(paste0(output_dir,"cal.eth_recal.pdf"),width=8,height=8)
cal.eth
dev.off()

png(paste0(output_dir,"cal.eth_recal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
cal.eth
dev.off()

# #Uncalibrated
# cal.white.uncal.p <- wrap_elements(cal.white.uncal+plot_annotation(title = paste0("White (",prop.tab.hb[1],"%)"), theme = thm))
# cal.asian.uncal.p <- wrap_elements(cal.asian.uncal+plot_annotation(title = paste0("South Asian (",prop.tab.hb[2],"%)"), theme = thm))
# cal.black.uncal.p <- wrap_elements(cal.black.uncal+plot_annotation(title = paste0("Black (",prop.tab.hb[3],"%)"), theme = thm))
# cal.mixed.other.uncal.p <- wrap_elements(cal.mixed.other.uncal+plot_annotation(title = paste0("Mixed or Other (",prop.tab.hb[4],"%)"), theme = thm))
# 
# cal.eth <- 
#   cal.white.uncal.p +
#   cal.asian.uncal.p + 
#   cal.black.uncal.p + 
#   cal.mixed.other.uncal.p +
#   plot_layout(ncol = 2) 
# 
# grDevices::cairo_pdf(paste0(output_dir,"cal.eth_uncal.pdf"),width=8,height=8)
# cal.eth
# dev.off()
# 
# png(paste0(output_dir,"cal.eth_uncal.png"),width=1800,height=1800,res=pngres,restoreConsole=TRUE)
# cal.eth
# dev.off()


## Save cal.eth_recal_SA (sensitivity analysis) --------------------------------

#cal.mixed.recal.p <- wrap_elements(cal.mixed.recal + plot_annotation(title = paste0("Mixed (",prop.tab.hb_sa[5],"%)"), theme = thm))
#cal.other.recal.p <- wrap_elements(cal.other.recal + plot_annotation(title = paste0("Other (",prop.tab.hb_sa[4],"%)"), theme = thm))

#cal.eth_SA <- 
#  cal.mixed.recal.p +
#  cal.other.recal.p + 
#  plot_layout(nrow = 1) 

#grDevices::cairo_pdf(paste0(output_dir, "cal.eth_recal_SA.pdf"), width = 8, height = 4.8)
#cal.eth_SA
#dev.off()

#png(paste0(output_dir,"cal.eth_recal_SA.png"), width = 1800, height = 1000, res = pngres, restoreConsole = TRUE)
#cal.eth_SA
#dev.off()



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Output HbA1c response (change from baseline) plots ---------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --


  # Variable prep
  final.hb<- final.hb %>% mutate(hba1c.breaks = cut(hba1c_diff.recal, breaks=c(min(hba1c_diff.recal)-0.01,-5,-3,0,3,max(hba1c_diff.recal)+0.01)),
                                 hba1c.change=posthba1c_final-prehba1cmmol)
  
  # Unadjusted
  hb.res.unadjusted <- final.hb %>% 
    mutate(hba1c.change=posthba1c_final-prehba1cmmol) %>%
    dplyr::group_by(ethnicity,hba1c.breaks,drugclass) %>% 
    dplyr::summarise(n=length(hba1c.change),
                     hba1c.resp = mean(hba1c.change),
                     lci = mean(hba1c.change) - (1.96*(sd(hba1c.change)/sqrt(length(hba1c.change)))),
                     uci = mean(hba1c.change) + (1.96*(sd(hba1c.change)/sqrt(length(hba1c.change))))) %>%
    select(ethnicity,hba1c.breaks,drugclass,n,hba1c.resp,lci,uci)
  
  # Long to wide for plotting
  hb.res.unadjusted <- hb.res.unadjusted %>% pivot_wider(
    names_from = drugclass,
    values_from = c(n,hba1c.resp,lci,uci)
  )

  # hb.res.unadjusted$diff <- hb.res.unadjusted$hba1c.change_SGLT2- hb.res.unadjusted$hba1c.change_DPP4

  #Plot observed and 95% CI by ethnicity and HbA1c defined subgroup
  white <- hb.res.unadjusted %>% filter(ethnicity=="White") 
  asian <- hb.res.unadjusted %>% filter(ethnicity=="South Asian")
  black <- hb.res.unadjusted %>% filter(ethnicity=="Black") 
  mixed <- hb.res.unadjusted %>% filter(ethnicity=="mixed.other") 
  
  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "South Asian",
                "Black",
                "Mixed or Other")
  
  hb.plot <- list()
  
  
  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,6],plotdata[,5]))
    cim = data.matrix(cbind(plotdata[,8],plotdata[,7]))
    cip = data.matrix(cbind(plotdata[,10],plotdata[,9]))
    
    hb.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.hb.resp,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Average glycaemic response (mmol/mol)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  
  #Plot together
  p.hb.white <- grid2grob(print(hb.plot[[1]]))
  p.hb.asian <- grid2grob(print(hb.plot[[2]]))
  p.hb.black <- grid2grob(print(hb.plot[[3]]))
  p.hb.mixed.other <- grid2grob(print(hb.plot[[4]]))
  
  hb.eth.obs <- (wrap_elements(p.hb.white) + wrap_elements(p.hb.asian)) /
    (wrap_elements(p.hb.black) + wrap_elements(p.hb.mixed.other))/
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  hb.eth.obs
  
  #grDevices::cairo_pdf(paste0(output_dir,"hb_eth_recal_observed_unadjusted.pdf"),width=pdfwidth,height=pdfheight)
  grDevices::cairo_pdf(paste0(output_dir,"hb_eth_recal_observed_unadjusted.pdf"), width = 14, height  = 14)
  hb.eth.obs
  dev.off()
  
  png(paste0(output_dir,"hb_eth_recal_observed_unadjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  hb.eth.obs
  dev.off()
  
  write.csv(hb.res.unadjusted,file=paste0(output_dir,"hb_eth_recal_observed_unadjusted.csv"))
  
  
  
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Output Discontinuation plots -------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
  

# Set up 

  # Set dummy hba1cmonth at 6 months
  final.dc <- final.dc %>% dplyr::mutate(hba1cmonth=6)
  
  # Predict HbA1c outcome uncalibrated
  final.dc <- final.dc %>% mutate(drug=drugclass,
                                  drugclass="DPP4")
  final.dc$DPP4.pred.lm.uncal <- predict(m1,final.dc)
  final.dc <- final.dc %>% mutate(drugclass="SGLT2")
  final.dc$SGLT2.pred.lm.uncal <- predict(m1,final.dc)
  final.dc <- final.dc %>% 
    mutate(hba1c_diff.uncal = SGLT2.pred.lm.uncal-DPP4.pred.lm.uncal,
           bestdrug.uncal=ifelse(hba1c_diff.uncal<=0,"SGLT2","DPP4"),
           drugclass=drug)
  head(final.dc)
 
  # Recalibrate HbA1c outcome to AURUM
  final.dc <- final.dc %>% mutate(DPP4.pred.lm.recal=ifelse(ethnicity == "White",DPP4.pred.lm.uncal+ctfm$intercept[1],
                                                            ifelse(ethnicity == "South Asian",DPP4.pred.lm.uncal+ctfm$intercept[3],
                                                                   ifelse(ethnicity == "Black",DPP4.pred.lm.uncal+ctfm$intercept[5],
                                                                          ifelse(ethnicity == "mixed.other",DPP4.pred.lm.uncal+ctfm$intercept[7],NA
                                                                          )))),
                                  SGLT2.pred.lm.recal=ifelse(ethnicity == "White",SGLT2.pred.lm.uncal+ctfm$intercept[2],SGLT2.pred.lm.uncal),
                                  hba1c_diff.recal = SGLT2.pred.lm.recal-DPP4.pred.lm.recal,
                                  bestdrug.recal=ifelse(hba1c_diff.recal<=0,"SGLT2","DPP4"))
  
  # Generate other variables for modelling
  
  # Define hba1c.breaks based on recalibrated HbA1c outcome
  final.dc <- final.dc %>% mutate(hba1c.breaks = cut(hba1c_diff.recal, breaks=c(min(hba1c_diff.recal,na.rm=T)-0.01,-5,-3,0,3,max(hba1c_diff.recal,na.rm=T)+0.01)))
  describe(final.dc$hba1c.breaks)
  
  # Define datadist for modelling
  ddist <- datadist(final.dc); options(datadist='ddist') 
  
## Observed discontinuation, by ethnicity and HbA1c benefit subgroup -----------
  
  # Harrell method
  dc.res.unadjusted <- final.dc %>% 
    dplyr::group_by(ethnicity,hba1c.breaks,drugclass) %>% 
    dplyr::summarise(n=length(stopdrug_6m_3mFU),
                     n.discontinued = sum(stopdrug_6m_3mFU),
                     #prop = sum(stopdrug_6m_3mFU) / length(stopdrug_6m_3mFU),
                     prop.conf = binconf(x=sum(stopdrug_6m_3mFU), n = length(stopdrug_6m_3mFU)),
                     prop=prop.conf[1],
                     lci=prop.conf[2],
                     uci=prop.conf[3]) %>%
    select(ethnicity,hba1c.breaks,drugclass,n,n.discontinued,prop,lci,uci)
  
  # Long to wide for plotting
  dc.res.unadjusted <- dc.res.unadjusted %>% pivot_wider(
    names_from = drugclass,
    values_from = c(n,n.discontinued,prop,lci,uci)
  )
  
  # Plot observed and 95% CI by ethnicity and HbA1c defined subgroup
  
  white <- dc.res.unadjusted %>% filter(ethnicity=="White") 
  asian <- dc.res.unadjusted %>% filter(ethnicity=="South Asian")
  black <- dc.res.unadjusted %>% filter(ethnicity=="Black") 
  mixed <- dc.res.unadjusted %>% filter(ethnicity=="mixed.other") 

  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "South Asian",
                "Black",
                "Mixed or Other")
  
  dc.plot <- list()
  

  for(i in 1:4) 
  { 
    #Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    #plot
    coef = data.matrix(cbind(plotdata[,8]*100,plotdata[,7]*100))
    cim = data.matrix(cbind(plotdata[,10]*100,plotdata[,9]*100))
    cip = data.matrix(cbind(plotdata[,12]*100,plotdata[,11]*100))
    dc.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.dc,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Proportion discontinuing therapy (%)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  # Plot together
  p.dc.white <- grid2grob(print(dc.plot[[1]]))
  p.dc.asian <- grid2grob(print(dc.plot[[2]]))
  p.dc.black <- grid2grob(print(dc.plot[[3]]))
  p.dc.mixed.other <- grid2grob(print(dc.plot[[4]]))
  
  dc.eth.obs <- (wrap_elements(p.dc.white) + wrap_elements(p.dc.asian)) /
    (wrap_elements(p.dc.black) + wrap_elements(p.dc.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  dc.eth.obs
  
  #grDevices::cairo_pdf(paste0(output_dir,"dc_eth_recal_observed_unadjusted.pdf"),width=pdfwidth,height=pdfheight)
  grDevices::cairo_pdf(paste0(output_dir,"dc_eth_recal_observed_unadjusted.pdf"), width = 14, height  = 14)
  dc.eth.obs
  dev.off()
  
  png(paste0(output_dir,"dc_eth_recal_observed_unadjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  dc.eth.obs
  dev.off()
  
  write.csv(dc.res.unadjusted,file=paste0(output_dir,"dc_eth_recal_observed_unadjusted.csv"))


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --
# Output Weight plots ----------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --  
  

  # Set up 
  
  # Set dummy hba1cmonth at 6 months
  final.wt <- final.wt %>% dplyr::mutate(hba1cmonth=6)
  
  # Predict HbA1c outcome uncalibrated
  final.wt <- final.wt %>% mutate(drug=drugclass,
                                  drugclass="DPP4")
  final.wt$DPP4.pred.lm.uncal  <- predict(m1,final.wt)
  final.wt                     <- final.wt %>% mutate(drugclass="SGLT2")
  final.wt$SGLT2.pred.lm.uncal <- predict(m1,final.wt)
  final.wt                     <- final.wt %>% 
    mutate(hba1c_diff.uncal = SGLT2.pred.lm.uncal-DPP4.pred.lm.uncal,
           bestdrug.uncal=ifelse(hba1c_diff.uncal<=0,"SGLT2","DPP4"),
           drugclass=drug)

  # Recalibrate HbA1c outcome to AURUM
  final.wt <- final.wt %>% mutate(DPP4.pred.lm.recal=ifelse(ethnicity == "White",DPP4.pred.lm.uncal+ctfm$intercept[1],
                                                            ifelse(ethnicity == "South Asian",DPP4.pred.lm.uncal+ctfm$intercept[3],
                                                                   ifelse(ethnicity == "Black",DPP4.pred.lm.uncal+ctfm$intercept[5],
                                                                          ifelse(ethnicity == "mixed.other",DPP4.pred.lm.uncal+ctfm$intercept[7],NA
                                                                          )))),
                                  SGLT2.pred.lm.recal=ifelse(ethnicity == "White",SGLT2.pred.lm.uncal+ctfm$intercept[2],SGLT2.pred.lm.uncal),
                                  hba1c_diff.recal = SGLT2.pred.lm.recal-DPP4.pred.lm.recal,
                                  bestdrug.recal=ifelse(hba1c_diff.recal<=0,"SGLT2","DPP4"))
  
# Generate other variables for modelling
  
  # Define hba1c.breaks based on recalibrated HbA1c outcome
  final.wt <- final.wt %>% mutate(hba1c.breaks = cut(hba1c_diff.recal, breaks=c(min(hba1c_diff.recal,na.rm=T)-0.01,-5,-3,0,3,max(hba1c_diff.recal,na.rm=T)+0.01)))
  describe(final.wt$hba1c.breaks)
  
  # Set HbA1c diff to recalibrated HbA1c outcome
  final.wt <- final.wt %>% dplyr::rename(hba1c_diff="hba1c_diff.recal",
                                         bestdrug="bestdrug.recal")
  
  # Define datadist for modelling
  ddist <- datadist(final.wt); options(datadist='ddist') 
  
  
## Observed weight, by ethnicity and HbA1c benefit subgroup --------------------
  
  # Harrell method
  wt.res.unadjusted <- final.wt %>% 
    dplyr::group_by(ethnicity,hba1c.breaks,drugclass) %>% 
    dplyr::summarise(n=length(wtchange),
                     wt.change = mean(wtchange),
                     lci = mean(wtchange) - (1.96*(sd(wtchange)/sqrt(length(wtchange)))),
                     uci = mean(wtchange) + (1.96*(sd(wtchange)/sqrt(length(wtchange))))) %>%
    select(ethnicity,hba1c.breaks,drugclass,n,wt.change,lci,uci)
  
  # Long to wide for plotting
  wt.res.unadjusted <- wt.res.unadjusted %>% pivot_wider(
    names_from = drugclass,
    values_from = c(n,wt.change,lci,uci)
  )
  
  # Plot observed and 95% CI by ethnicity and HbA1c defined subgroup
  
  white <- wt.res.unadjusted %>% filter(ethnicity=="White") 
  asian <- wt.res.unadjusted %>% filter(ethnicity=="South Asian")
  black <- wt.res.unadjusted %>% filter(ethnicity=="Black") 
  mixed <- wt.res.unadjusted %>% filter(ethnicity=="mixed.other") 
  
  L <- list(white,asian,black,mixed)
  
  names(L) <- c("White",
                "South Asian",
                "Black",
                "Mixed or Other")
  
  wt.plot <- list()
  
  for(i in 1:4) 
  { 
    # Subgroups by predicted treatment difference
    plotdata <- L[[i]] %>% ungroup() %>% as.data.frame()  %>%
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 1) %>% 
      add_row(!!! setNames(list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),names(.)), .before = 5) 
    
    # plot
    coef = data.matrix(cbind(plotdata[,6],plotdata[,5]))
    cim = data.matrix(cbind(plotdata[,8],plotdata[,7]))
    cip = data.matrix(cbind(plotdata[,10],plotdata[,9]))
    wt.plot[[i]] <-
      forestplot(row_names,
                 mean = coef,
                 lower= cim,
                 upper = cip,
                 hrzl_lines = gpar(col="#444444"),lineheight=unit(2,'cm'),
                 is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
                 xticks = tick.wt,
                 zero = 0,
                 #boxsize=0.1,
                 # graphwidth = unit(2,"inches"),
                 # lineheight = unit(0.7,"inches"),
                 ci.vertices=TRUE,
                 col=fpColors(box=c("#f1a340","#4118de"), lines=c("#f1a340","#4118de"), zero = "gray50"),
                 lty.ci = c(1,2),           ,
                 xlab="Average weight change (kg)",cex=1,
                 title = names(L[i]),
                 new_page = TRUE,
                 fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
                 boxsize = .2, # We set the box size to better visualize the type
                 #line.margin = .1, # We need to add this to avoid crowding
                 txt_gp = fpTxtGp(legend  = gpar(cex = 1),xlab  = gpar(cex = 1),summary=gpar(cex = 1), ticks  = gpar(cex = 1), title=gpar(cex=1.3))
                 #txt_gp = fpTxtGp(label= gpar(cex = 0.7),ticks  = gpar(cex = 0.7),xlab  = gpar(cex = 0.7)),
                 #legend = c("SGLT2i","DPP4i"),
                 #legend_args = fpLegend(pos = list(x=.70, y=0.95))#, 
                 #gp=gpar(col="#CCCCCC", fill="#F9F9F9"))#,
                 #xlog = TRUE
      )
    
  }  
  
  
  # Plot together
  p.wt.white <- grid2grob(print(wt.plot[[1]]))
  p.wt.asian <- grid2grob(print(wt.plot[[2]]))
  p.wt.black <- grid2grob(print(wt.plot[[3]]))
  p.wt.mixed.other <- grid2grob(print(wt.plot[[4]]))
  
  wt.eth.obs <- (wrap_elements(p.wt.white) + wrap_elements(p.wt.asian)) /
    (wrap_elements(p.wt.black) + wrap_elements(p.wt.mixed.other)) /
    shared_legend +
    plot_layout(height=c(10,10,1))
  
  wt.eth.obs
  
  # grDevices::cairo_pdf(paste0(output_dir,"wt_eth_recal_observed_unadjusted.pdf"),width=pdfwidth,height=pdfheight)
  grDevices::cairo_pdf(paste0(output_dir,"wt_eth_recal_observed_unadjusted.pdf"), width = 14, height  = 14)
  wt.eth.obs
  dev.off()
  
  png(paste0(output_dir,"wt_eth_recal_observed_unadjusted.png"),width=pngwidth,height=pngheight,res=pngres,restoreConsole=TRUE)
  wt.eth.obs
  dev.off()
  
  write.csv(wt.res.unadjusted,file=paste0(output_dir,"wt_eth_recal_observed_unadjusted.csv"))
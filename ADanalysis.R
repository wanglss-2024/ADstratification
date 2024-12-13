#------------ Import libraries--------------#
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(PheNorm)
library(MAP)
library(randomForest)
library(glmnet)
library(dplyr)
library(lubridate)

################################################
### Prepare the data
################################################
setwd("~/AD_stratification/analysis")
library(pROC)
library(dplyr)
library(wordcloud2)
library(scales)
library(readxl)
library(fmsb)


load("data/AD_dat_2_0.RData")
# N = 10000
# n = 300
Y.all = Y.all[,colnames(Y.all) %in% features.once]
p = ncol(Y.all); p
q = 5 # also try 10, 15, 20

set.seed(0)
X.all$AD = X.all$AD.x
ind.L = which((!is.na(X.all$AD)) & (X.all$AD !="Not")); length(ind.L)
ind.U.all = which(is.na(X.all$AD))
ind.U = sample(which(is.na(X.all$AD)), N, replace=F)
ind.L.train = sample(ind.L, n, replace=F)
ind.L.test = ind.L[!ind.L %in% ind.L.train]; length(ind.L.test)

dat.L.train = prepare_data(Y.all[ind.L.train,], X.all[ind.L.train,]); nrow(dat.L.train)
dat.L.test = prepare_data(Y.all[ind.L.test,], X.all[ind.L.test,]); nrow(dat.L.test)
dat.L.train$AD = ifelse(dat.L.train$AD == "Possible", 1, 0)
dat.L.test$AD = ifelse(dat.L.test$AD == "Possible", 1, 0)
dat.U.all = prepare_data(Y.all, X.all); nrow(dat.U.all)
dat.U.all$AD = 3
dat.U = prepare_data(Y.all[ind.U,], X.all[ind.U,]); nrow(dat.U)
dat.U$AD = 3
dat.all.train = rbind(dat.L.train, dat.U)

dat.L.train = dat.L.train[rowSums(dat.L.train$Abundance)>0,]
dat.L.test = dat.L.test[rowSums(dat.L.test$Abundance)>0,]
dat.U.all = dat.U.all[rowSums(dat.U$Abundance)>0,]
n = nrow(dat.L.train); n
N = nrow(dat.U); N

########################################
### Fit the model
########################################
lb = list(B = matrix(-Inf, nrow=3, ncol=q),
          L = matrix(0, nrow=q, ncol=q),
          M = matrix(-Inf, nrow=N, ncol=q),
          S = matrix(0, nrow=N, ncol=q))
# Unsupervised fixed V
pln.fixedv.unsup = PLNfixedVunsup(Abundance ~  Intercept + log(utils+1), data  = dat.U, control=PLN_param(V=V[,1:q], #eigen(cov(Y.all[ind.U,]))$vectors[,1:q],
                                                                                                          rank=q,
                                                                                                          config_optim = list("algorithm" = "MMA",
                                                                                                                              "lower_bounds" = lb,
                                                                                                                              "ftol_out" = 1e-6,
                                                                                                                              "maxit_out" = 10,
                                                                                                                              "L0" = diag(0.5, q),
                                                                                                                              "tau0" = matrix(0.5, ncol=2, nrow=N))))
########################################
### KM curve
########################################
# labeled data
library(survival)
library(survminer)
library(patchwork)

pred.U = predict(pln.fixedv.unsup, dat.U.all, type="posterior")[,2]
AD.U = ifelse(pred.U >= mean(pred.U), "Group 1", "Group 0")
dat.U.pred = dat.U.all %>% mutate(AD = AD.U) %>% select(-Abundance, -Offset)
dat.U.0 = X.all[is.na(rowSums(Y.all)),]
dat.U.0$AD = "Group 1"
colnames(dat.U.0) = colnames(dat.U.pred)
dat.all.pred = rbind(dat.U.pred, dat.U.0) 

dat.U.pred$AD = ifelse(dat.U.pred$AD=="Group 0", 0, 1)
df.mean = dat.U.pred %>% group_by(AD) %>% summarise(across(white:codect_2yr, ~ mean(.x, na.rm = TRUE)))

CS.nh = coxph(Surv(time = event_code_month, event_code == "NursingHome") ~ AD + white + female + age + log(codect_2yr+1) + 
                     chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+ 
                     cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                     metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                     dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data=dat.all.pred, x=TRUE)
cox.nh = survfit(CS.nh, df.mean)
CS.nh = coxph(Surv(time = event_code_month, event_code == "NursingHome") ~ AD, data=dat.all.pred, x=TRUE)
g1p = ggsurvplot(cox.nh, 
                  data = df.mean, 
                  xlim=c(-2, 140),
                 break.time.by = 25,
                 size = 1,                 # change line size
                 conf.int = TRUE,          # Add confidence interval
                 pval = F,              # Add p-value
                 xlab = "Months after AD diagnosis",
                 ylab="Adjusted\nSurvival\nProbability",
                 title = "Nursing Home Admission",
                 legend.labs = c("Group 1", "Group 0"),
                 risk.table = TRUE,        # Add risk table
                 risk.table.col = "strata",# Risk table color by groups
                 risk.table.height = 0.25, # Useful to change when you have multiple groups
                 surv.median.line = "hv",
                 ggtheme = theme_bw()      # Change ggplot2 theme)
)$plot + theme(axis.title.y = element_text(angle = 0, vjust = 1.05, hjust = 1.1))
kmfit.nh.tab = survfit(Surv(time = nh_code_month, nh_code) ~ factor(AD), data=dat.all.pred)
g1t = ggsurvtable(kmfit.nh.tab, risk.table.type = "nrisk_cumevents", data=dat.all.pred, color = "strata", y.text=F,
                  legend="none", ggtheme = theme_bw(), break.time.by =25, xlim=c(-2, 140),fontsize = 3)$risk.table +xlab("")+ylab("")+theme(legend.position = "none")
g1.all = g1p + g1t + plot_layout(ncol = 1, heights=c(4,1)); g1.all


kmfit.death = survfit(Surv(time = dat.U.pred$death_month, event = dat.U.pred$death==1) ~ factor(AD), data=dat.all.pred)
CS.death = coxph(Surv(time = dat.U.pred$death_month, event = dat.U.pred$death==1) ~ strata(AD) + white + female + age + log(codect_2yr+1) + 
                chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+ 
                cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data=dat.all.pred, x=TRUE)
cox.death = survfit(CS.death, df.mean)
CS.death = coxph(Surv(time = dat.U.pred$death_month, event = dat.U.pred$death==1) ~ AD, data=dat.all.pred, x=TRUE)

g2p = ggsurvplot(cox.death,
                 data = df.mean, 
                 xlim=c(-2, 140),
                 break.time.by = 25,
                 size = 1,                 # change line size
                 conf.int = TRUE,          # Add confidence interval
                 pval = F,              # Add p-value
                 xlab = "Months after AD diagnosis",
                 ylab="Adjusted\nSurvival\nProbability",
                 title = "Mortality",
                 legend.labs = c("Group 1", "Group 0"),
                 risk.table = TRUE,        # Add risk table
                 risk.table.col = "strata",# Risk table color by groups
                 risk.table.height = 0.25, # Useful to change when you have multiple groups
                 surv.median.line = "hv",
                 ggtheme = theme_bw()      # Change ggplot2 theme)
)$plot + theme(axis.title.y = element_text(angle = 0, vjust = 1.05, hjust = 1.1))
kmfit.death.tab = survfit(Surv(time = death_month, death) ~ factor(AD), data=dat.all.pred)
g2t = ggsurvtable(kmfit.death.tab, risk.table.type = "nrisk_cumevents", data=dat.all.pred, color = "strata", y.text=F,
                  legend="none", ggtheme = theme_bw(), break.time.by =25, xlim=c(-2, 140),fontsize = 3)$risk.table +xlab("")+ylab("")+theme(legend.position = "none")
g2.all = g2p + g2t + plot_layout(ncol = 1, heights=c(4,1)); g2.all


########################################
### Word Cloud and Correlation Plots
########################################
library(ggwordcloud)
library(heatmaply)
library(reshape2)

ONCE.cod.des = read.csv(file="features/ONCE_AD_cod.csv") %>%
  select(c(Variable, Description))
ONCE.nlp.des = read.csv(file="features/ONCE_AD_nlp.csv") %>% select(c(cui, term))

pred.U = predict(pln.fixedv.unsup, dat.U.all, type="posterior")[,2]
AD.U = ifelse(pred.U >= mean(pred.U), "Group 1", "Group 0")
dat.U.pred = dat.U.all %>% mutate(AD = AD.U)

# Difference between groups
res = matrix(NA, nrow=nrow(dat.U.pred$Abundance), ncol=ncol(dat.U.pred$Abundance))
for (j in 1:ncol(dat.U.pred$Abundance)){
  res[,j] = lm(log(dat.U.pred$Abundance[,j]+1)~log(dat.U.pred$utils+1))$residuals
}
res = data.frame(res)
res$AD = dat.U.pred$AD
res.group = res %>% group_by(AD) %>% summarise_all("mean")
res.dif = data.frame(feature = colnames(dat.U.all$Abundance), 
                     importance=t(res.group[2,-1] - res.group[1,-1])* apply(log(dat.U.pred$Abundance+1), MARGIN=2, FUN=sd)) %>% 
  left_join(ONCE.cod.des, by=c("feature"="Variable")) %>% left_join(ONCE.nlp.des, by=c("feature"="cui")) %>%
  mutate(Term = ifelse(is.na(Description), term, Description))

pos = res.dif %>% arrange(desc(importance)) %>% head(20)
neg = res.dif %>% arrange(importance) %>% head(30) %>% filter(feature != "C0168634")
pos$color = "#00BFC4"
neg$color = "#F8766D"

all = rbind(pos%>% head(15), 
            neg %>%head(15))

id.samp = sample(1:nrow(all), nrow(all), replace=F)
ggplot(all[id.samp,], aes(label = Term, size=percentile(abs(importance))+10, color = color)) +
  geom_text_wordcloud_area() +
  scale_size_area(max_size = 30) +
  theme_minimal() #+

p.select.neg = 15
p.select.pos = 15
cor.mat = cor(res[,-ncol(res)])
rownames(cor.mat) = colnames(cor.mat) =  colnames(dat.U.all$Abundance)
cor.mat.pos = cor.mat[pos$feature[1:p.select.pos], pos$feature[1:p.select.pos]]
cor.mat.neg = cor.mat[neg$feature[1:p.select.neg], neg$feature[1:p.select.neg]]

rownames(cor.mat.neg) = paste0(c(neg$feature[1:p.select.neg]), " (",
                               c(neg$Term[1:p.select.neg]), ")")
rownames(cor.mat.pos) = paste0(c(pos$feature[1:p.select.pos]), " (",
                               c(pos$Term[1:p.select.pos]), ")")

heatmaply_cor(x = cor.mat.neg, xlab = "Features", 
              ylab = "Features", k_col = 1, k_row = 1) #740 * 500
heatmaply_cor(x = cor.mat.pos, xlab = "Features", 
              ylab = "Features", k_col = 1, k_row = 1) #990 * 500


cor.mat.all = cor.mat[c(pos$feature[1:p.select.pos],neg$feature[1:p.select.neg]), 
                      c(pos$feature[1:p.select.pos], neg$feature[1:p.select.neg])]

rownames(cor.mat.all) = c(paste0(c(pos$feature[1:p.select.pos]), " (",
                                 c(pos$Term[1:p.select.pos]), ")"),
                          paste0(c(neg$feature[1:p.select.neg]), " (",
                                 c(neg$Term[1:p.select.neg]), ")"))
heatmaply_cor(x = cor.mat.all, xlab = "Features", 
              ylab = "Features", k_col = 1, k_row = 1) #1250*680


########################################
### Comorbidities
########################################
library(tidyr)
pval_sig = Vectorize(function(p) {
  if ((p <= 0.05)& (p >0.01)) {return ("*")}
  else if ((p <= 0.01)& (p >0.001)) {return ("**")}
  else if ((p <= 0.001)& (p >0.0001)) {return("***")}
  else if (p <= 0.0001) {return("****")}
  else {return("")}
})
prop = dat.all.pred %>% group_by(AD) %>% summarise(n = n(), across(chf_2yr:depre_2yr, ~sum(.))) %>% gather(features, count, chf_2yr:depre_2yr)
prop = prop %>% mutate(prop = count / n, 
                       lower = prop-1.96*sqrt(prop*(1-prop)/n), 
                       upper = prop+1.96*sqrt(prop*(1-prop)/n))

pval = prop %>% group_by(features) %>% summarise(pval = prop.test(x=count, n=n)$p.value) %>% mutate(pval.sig = pval_sig(pval))
features.df = read.csv("mapping_comorbidity.csv")
prop = prop %>% left_join(pval, by="features") %>% 
  mutate(AD = ifelse(AD=="Group 1", "Slow decline", "Fast decline"),
         features = features) %>% 
  left_join(features.df, by=c("features"="ac"))
base = ggplot(prop, aes(x=reorder(nm,prop), y=prop)) +
  geom_point(aes(color=as.factor(AD), group=as.factor(AD)), position=position_dodge(width = 0.5)) +
  geom_errorbar(width=.1, position=position_dodge(width = 0.5), aes(color=as.factor(AD), group=as.factor(AD), ymin=lower, ymax=upper)) +
  geom_text(aes(x=nm, y=0.95, label=pval.sig),color="red") +
  ylim(0, 1) +
  xlab("Comorbidity") + ylab("% with comorbidity") +
  labs(color = "Race") +
  coord_flip() +
  theme_bw() + 
  ggtitle("Baseline Comorbidities by Group") +
  scale_color_manual(name = "Group", values=c("#F8766D","#00BFC4"))
base


########################################
### Demographics Table
########################################
library(comorbidity)
library(table1)

# Granular race/ethinicty
demo.imp = read.csv("data/UPMC demographics imputed race ethnicity separate.csv")
d = demo.imp %>% select(c("PATIENT_STUDY_ID", "race", "ethnicity")) %>% mutate(patient_num = as.character(PATIENT_STUDY_ID))

dat.all.pred = dat.all.pred %>% left_join(d, by = "patient_num")

dat.all.pred$race_eth = case_when((dat.all.pred$race == "White") & (dat.all.pred$ethnicity == "Non Hispanic") ~ "Non-Hispanic White",
                                  dat.all.pred$ethnicity == "Hispanic" ~ "Hispanic",
                                  dat.all.pred$race == "Asian" ~ "Asian",
                                  dat.all.pred$race == "Black" ~ "Black",
                                  dat.all.pred$race == "American Indian or Alaska Native" ~ "American Indian or Alaska Native",
                                  dat.all.pred$race == "Native Hawaiian or Other Pacific Islander" ~ "Native Hawaiian or Other Pacific Islander")

weights = c(7, 5, -1, 4, 2, 0, 0, 7, 6, 3, 0, 0, 0, 5, 11, 0, 0, 9, 12, 4, 0, 3, -4, 6, 5, -2, -2, 0, -7, 0, -3)
dat.all.pred$elix_ind = as.vector(weights %*% t(as.matrix(dat.all.pred %>% select(chf_2yr:depre_2yr))))
dat.all.pred$ad_med = case_when(is.na(dat.all.pred$med_month) ~ "no",
                                    dat.all.pred$med_month < 0 ~ "before diagnosis",
                                    dat.all.pred$med_month >= 0 ~ "after diagnosis")

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.0001)))
}
table1(~ age + factor(female) + factor(race_eth) + elix_ind + 
         factor(hypunc_2yr) + factor(carit_2yr) + factor(depre_2yr) + factor(fed_2yr) +
         factor(ad_med) + event_code_month
       | factor(AD), data=dat.all.pred, overall = T)
table1(~ age + factor(female) + factor(race_eth) + elix_ind + 
         factor(hypunc_2yr) + factor(carit_2yr) + factor(depre_2yr) + factor(fed_2yr) +
         factor(ad_med) + event_code_month
         | factor(AD), data=dat.all.pred, overall = F, extra.col=list(`P-value`=pvalue))


### Time of follow-up
# To get median time of followup from reverse KM estimator
fit.censor.upmc = survfit(Surv(dat.all.pred$event_month, dat.all.pred$event_code=="No") ~ 1)
quantile(fit.censor.upmc, probs = c(0.25, 0.5, 0.75))$quantile


########################################
### Table for prevalence of features
########################################
feat.ind = (dat.U.pred$Abundance > 0) %>% as.data.frame()
feat.ind$AD = dat.U.pred$AD
feat.prev = feat.ind %>% group_by(AD) %>% summarise(across(C0002395:`RXNORM:77492`, ~ mean(.x, na.rm = TRUE))) %>% as.data.frame()
feat.prev.df = data.frame(feature = colnames(feat.prev)[-1], 
                          all = colMeans(feat.ind %>% select(-AD)),
                          fast = unlist(c(feat.prev[1,-1])), 
                          slow = unlist(c(feat.prev[2,-1]))) %>%
  left_join(res.dif %>% select(feature, importance), by="feature") %>%
  left_join(ONCE.cod.des, by=c("feature"="Variable")) %>% left_join(ONCE.nlp.des, by=c("feature"="cui")) %>%
  mutate(Term = ifelse(is.na(Description), term, Description)) %>%
  arrange(desc(importance)) %>%
  select(feature, Term, importance, all, fast, slow)
  
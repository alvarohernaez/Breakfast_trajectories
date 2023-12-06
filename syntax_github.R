rm(list=ls())

library(compareGroups)
library(lme4)
library(car)
library(Hmisc)
library(foreign)
library(nlme)
library(haven)
library(lmtest)
library(tidyr)


### GUAPAS ###
##############

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.00001,signif(x,1),
                   ifelse(abs(x)<0.0001,signif(x,1),
                          ifelse(abs(x)<0.001,signif(x,1),
                                 ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                        ifelse(abs(x)<1,sprintf("%.2f",round(x,2)),
                                               ifelse(abs(x)<10,sprintf("%.2f",round(x,2)),
                                                      ifelse(abs(x)<100,sprintf("%.1f",round(x,1)),
                                                             ifelse(abs(x)>=100,round(x,0),round(x,0)))))))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

ic_guapa2<-function(x,y,z)
{
  ic<-paste(x," (",y," to ",z,")",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",
                      ifelse(abs(x)<0.01,sprintf("%.3f",round(x,3)),
                             ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                    ifelse(abs(x)<1,sprintf("%.3f",round(x,3)),guapa(x))))))
  return(pval)
}

pval_guapa2<-function(x)
{
  pval<-ifelse(x<0.00001," < 0.00001",
               ifelse(x<0.001," < 0.001",
                      ifelse(abs(x)<0.01,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),
                             ifelse(abs(x)<0.1,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),
                                    ifelse(abs(x)<1,paste(" = ",sprintf("%.3f",round(x,3)),sep=""),guapa(x))))))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

beta_se_ic_guapa2 <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa2(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa3 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),3)
  ic95a<-round(exp(x-(z*y)),3)
  ic95b<-round(exp(x+(z*y)),3)
  ic_ok<-ic_guapa2(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

z<-qnorm(1-0.05/2)

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }


setwd("...")
dir.create("./Data")
dir.create("./Splines_medialong")
dir.create("./Splines_medialong/women")
dir.create("./Splines_medialong/men")
dir.create("./Backup_medialong")
dir.create("./Backup_medialong/women")
dir.create("./Backup_medialong/men")
dir.create("./Splines_repeat")
dir.create("./Forest_repeat")
dir.create("./Breakfast_items")


### CALCULATION OF LONGITUDINAL AVERAGES AND DIFFERENCES ###

load("./Data/database.RData")
dat$ajus<-1
dat$ethnic<-0
dat$kcalp_long<-apply(dat[,c("kcal_brk_p00","kcal_brk_p02","kcal_brk_p03")],
                      1, function(row) mean(row, na.rm = TRUE))
dat$kcalday_long<-apply(dat[,c("kcal_day00","kcal_day02","kcal_day03")],
                        1, function(row) mean(row, na.rm = TRUE))
dat$score_long<-apply(dat[,c("brk_score00","brk_score02","brk_score03")],
                      1, function(row) mean(row, na.rm = TRUE))

dat$egfr_v00<-CKDEpi.creat(dat$creat_v00, dat$sexo, (dat$edad+0), dat$ethnic)
dat$egfr_v66<-CKDEpi.creat(dat$creat_v66, dat$sexo, (dat$edad+0.5), dat$ethnic)
dat$egfr_v01<-CKDEpi.creat(dat$creat_v01, dat$sexo, (dat$edad+1), dat$ethnic)
dat$egfr_v03<-CKDEpi.creat(dat$creat_v03, dat$sexo, (dat$edad+3), dat$ethnic)


### LONG DATABASES ###

dat<-rename.vars(dat,from=c("kcal_brk_p00","kcal_brk_p02","kcal_brk_p03",
                            "brk_score00","brk_score02","brk_score03"),
                 to=c("kcalp00","kcalp02","kcalp03",
                      "score00","score02","score03"))

dat$kcalp_cat<-with(dat,ifelse(kcalp_long<20,1,
                               ifelse(kcalp_long<=30,0,
                                      ifelse(kcalp_long>30,2,NA))))
quant01<-as.numeric(quantile(dat$score_long,prob=c(0.333333,0.666666),na.rm=TRUE,type=6))[1]
quant02<-as.numeric(quantile(dat$score_long,prob=c(0.333333,0.666666),na.rm=TRUE,type=6))[2]
# Tertile 1: 59.33064
# Tertile 2: 69.37053
#dat$score_cat<-car::recode(dat$score_long,"lo:quant01=1;quant01:quant02=2;quant02:hi=0")
#dat$score_cat<-with(dat,ifelse(score_cat==2,NA,score_cat))
quant<-as.numeric(quantile(dat$score_long,prob=c(0.5),na.rm=TRUE,type=6))[1]
dat$score_cat<-car::recode(dat$score_long,"lo:quant=1;quant:hi=0")


dat$kcalp_basal<-dat$kcalp00
dat$score_basal<-dat$score00
dat$glucosa_basal<-dat$glucosa_v00
dat$hb1ac_basal<-dat$hb1ac_v00
dat$ldl_basal<-dat$ldl_v00
dat$hdl_basal<-dat$hdl_v00
dat$trigli_basal<-dat$trigli_v00
dat$pas_basal<-dat$pas_v00
dat$pad_basal<-dat$pad_v00
dat$imc_basal<-dat$imc_v00
dat$imc00<-dat$imc_v00
dat$cintura_basal<-dat$cintura_v00
dat$egfr_basal<-dat$egfr_v00
dat$grup_int<-dat$intervencion
dat$diabetes00<-dat$diab_s1
dat$hipercolest00<-dat$colest_s1
dat$hta00<-dat$hta_s1
dat$tabaco00<-dat$tabacos1

dat$seg00<-0
dat$seg66<-(dat$fecha_v66-dat$fecha_v00)/365.25
dat$seg66<-with(dat,ifelse(is.na(seg66),0.5,seg66))
dat$seg01<-(dat$fecha_v01-dat$fecha_v00)/365.25
dat$seg01<-with(dat,ifelse(is.na(seg01),1,seg01))
dat$seg02<-(dat$fecha_v02-dat$fecha_v00)/365.25
dat$seg02<-with(dat,ifelse(is.na(seg02),2,seg02))
dat$seg03<-(dat$fecha_v03-dat$fecha_v00)/365.25
dat$seg03<-with(dat,ifelse(is.na(seg03),3,seg03))

dat$glucosa_v02<-NA
dat$hb1ac_v02<-NA
dat$ldl_v02<-NA
dat$hdl_v02<-NA
dat$trigli_v02<-NA
dat$egfr_v02<-NA
dat$kcalp66<-NA
dat$kcalp01<-NA
dat$score66<-NA
dat$score01<-NA

dat$grup_int<-with(dat,ifelse(grup_int=="A",0,1))


glucosa_vars<-c("id","ajus","kcalp_cat","kcalp_long","kcalp_basal",
                "score_cat","score_long","score_basal",
                "brk_mix","p17_long","af_long","kcalday_long",
                "edad","sexo","grup_int","escolar","diabetes00","hipercolest00","hta00","tabaco00","imc00",
                "glucosa_basal","hb1ac_basal","ldl_basal","hdl_basal","trigli_basal",
                "pas_basal","pad_basal","imc_basal","cintura_basal","egfr_basal",
                "glucosa_v00","glucosa_v66","glucosa_v01","glucosa_v02","glucosa_v03")
glucosa_wide<-dat[,glucosa_vars]
glucosa_long<-reshape(glucosa_wide, 
                      varying = c("glucosa_v00","glucosa_v66","glucosa_v01","glucosa_v02","glucosa_v03"), 
                      v.names = "glucosa",
                      timevar = "seg", 
                      times = c("0","0.5","1","2","3"), 
                      direction = "long")

hb1ac_vars<-c("id","hb1ac_v00","hb1ac_v66","hb1ac_v01","hb1ac_v02","hb1ac_v03")
hb1ac_wide<-dat[,hb1ac_vars]
hb1ac_long<-reshape(hb1ac_wide, 
                    varying = c("hb1ac_v00","hb1ac_v66","hb1ac_v01","hb1ac_v02","hb1ac_v03"), 
                    v.names = "hb1ac",
                    timevar = "seg", 
                    times = c("0","0.5","1","2","3"), 
                    direction = "long")

coltot_vars<-c("id","coltot_v00","coltot_v66","coltot_v01","coltot_v02","coltot_v03")
coltot_wide<-dat[,coltot_vars]
coltot_long<-reshape(coltot_wide, 
                     varying = c("coltot_v00","coltot_v66","coltot_v01","coltot_v02","coltot_v03"), 
                     v.names = "coltot",
                     timevar = "seg", 
                     times = c("0","0.5","1","2","3"), 
                     direction = "long")

ldl_vars<-c("id","ldl_v00","ldl_v66","ldl_v01","ldl_v02","ldl_v03")
ldl_wide<-dat[,ldl_vars]
ldl_long<-reshape(ldl_wide, 
                  varying = c("ldl_v00","ldl_v66","ldl_v01","ldl_v02","ldl_v03"), 
                  v.names = "ldl",
                  timevar = "seg", 
                  times = c("0","0.5","1","2","3"), 
                  direction = "long")

hdl_vars<-c("id","hdl_v00","hdl_v66","hdl_v01","hdl_v02","hdl_v03")
hdl_wide<-dat[,hdl_vars]
hdl_long<-reshape(hdl_wide, 
                  varying = c("hdl_v00","hdl_v66","hdl_v01","hdl_v02","hdl_v03"), 
                  v.names = "hdl",
                  timevar = "seg", 
                  times = c("0","0.5","1","2","3"), 
                  direction = "long")

trigli_vars<-c("id","trigli_v00","trigli_v66","trigli_v01","trigli_v02","trigli_v03")
trigli_wide<-dat[,trigli_vars]
trigli_long<-reshape(trigli_wide, 
                     varying = c("trigli_v00","trigli_v66","trigli_v01","trigli_v02","trigli_v03"), 
                     v.names = "trigli",
                     timevar = "seg", 
                     times = c("0","0.5","1","2","3"), 
                     direction = "long")

egfr_vars<-c("id","egfr_v00","egfr_v66","egfr_v01","egfr_v02","egfr_v03")
egfr_wide<-dat[,egfr_vars]
egfr_long<-reshape(egfr_wide, 
                   varying = c("egfr_v00","egfr_v66","egfr_v01","egfr_v02","egfr_v03"), 
                   v.names = "egfr",
                   timevar = "seg", 
                   times = c("0","0.5","1","2","3"), 
                   direction = "long")

imc_vars<-c("id","imc_v00","imc_v66","imc_v01","imc_v02","imc_v03")
imc_wide<-dat[,imc_vars]
imc_long<-reshape(imc_wide, 
                  varying = c("imc_v00","imc_v66","imc_v01","imc_v02","imc_v03"), 
                  v.names = "imc",
                  timevar = "seg", 
                  times = c("0","0.5","1","2","3"), 
                  direction = "long")

cintura_vars<-c("id","cintura_v00","cintura_v66","cintura_v01","cintura_v02","cintura_v03")
cintura_wide<-dat[,cintura_vars]
cintura_long<-reshape(cintura_wide, 
                      varying = c("cintura_v00","cintura_v66","cintura_v01","cintura_v02","cintura_v03"), 
                      v.names = "cintura",
                      timevar = "seg", 
                      times = c("0","0.5","1","2","3"), 
                      direction = "long")

pas_vars<-c("id","pas_v00","pas_v66","pas_v01","pas_v02","pas_v03")
pas_wide<-dat[,pas_vars]
pas_long<-reshape(pas_wide, 
                  varying = c("pas_v00","pas_v66","pas_v01","pas_v02","pas_v03"), 
                  v.names = "pas",
                  timevar = "seg", 
                  times = c("0","0.5","1","2","3"), 
                  direction = "long")

pad_vars<-c("id","pad_v00","pad_v66","pad_v01","pad_v02","pad_v03")
pad_wide<-dat[,pad_vars]
pad_long<-reshape(pad_wide, 
                  varying = c("pad_v00","pad_v66","pad_v01","pad_v02","pad_v03"), 
                  v.names = "pad",
                  timevar = "seg", 
                  times = c("0","0.5","1","2","3"), 
                  direction = "long")

kcalp_vars<-c("id","kcalp00","kcalp66","kcalp01","kcalp02","kcalp03")
kcalp_wide<-dat[,kcalp_vars]
kcalp_long<-reshape(kcalp_wide, 
                    varying = c("kcalp00","kcalp66","kcalp01","kcalp02","kcalp03"), 
                    v.names = "kcalp",
                    timevar = "seg", 
                    times = c("0","0.5","1","2","3"), 
                    direction = "long")

score_vars<-c("id","score00","score66","score01","score02","score03")
score_wide<-dat[,score_vars]
score_long<-reshape(score_wide, 
                    varying = c("score00","score66","score01","score02","score03"), 
                    v.names = "score",
                    timevar = "seg", 
                    times = c("0","0.5","1","2","3"), 
                    direction = "long")

segyears_vars<-c("id","seg00","seg66","seg01","seg02","seg03")
segyears_wide<-dat[,segyears_vars]
segyears_long<-reshape(segyears_wide, 
                       varying = c("seg00","seg66","seg01","seg02","seg03"), 
                       v.names = "segyears",
                       timevar = "seg", 
                       times = c("0","0.5","1","2","3"), 
                       direction = "long")


dat_long<-as.data.frame(cbind(glucosa_long,hb1ac_long[,3],ldl_long[,3],hdl_long[,3],
                              trigli_long[,3],pas_long[,3],pad_long[,3],imc_long[,3],cintura_long[,3],
                              egfr_long[,3],kcalp_long[,3],score_long[,3],segyears_long[,3]))
colnames(dat_long)<-c("id","ajus","kcalp_cat","kcalp_long","kcalp_basal",
                      "score_cat","score_long","score_basal",
                      "brk_mix","p17_long","af_long","kcalday_long",
                      "edad","sexo","grup_int","escolar","diabetes00","hipercolest00","hta00","tabaco00","imc00",
                      "glucosa_basal","hb1ac_basal","ldl_basal","hdl_basal","trigli_basal",
                      "pas_basal","pad_basal","imc_basal","cintura_basal","egfr_basal",
                      "seg","glucosa","hb1ac","ldl","hdl","trigli","pas","pad",
                      "imc","cintura","egfr","kcalp","score","segyears")
dat_long$seg<-as.numeric(dat_long$seg)*12
dat_long$segyears<-as.numeric(dat_long$segyears)*12

save(dat_long,file="./Data/dat_long.RData")


### WIDE DATABASE ###

dat<-dat[,c("id","ajus","kcalp00","kcalp_long","score00","score_long",
            "edad","sexo","grup_int","escolar","diabetes00","hipercolest00","hta00","tabaco00","imc00","p17_long","af_long","kcalday_long",
            "glucosa_v00","glucosa_long",
            "hb1ac_v00","hb1ac_long",
            "ldl_v00","ldl_long",
            "hdl_v00","hdl_long",
            "trigli_v00","trigli_long",
            "pas_v00","pas_long",
            "pad_v00","pad_long",
            "imc_v00","imc_long",
            "cintura_v00","cintura_long",
            "egfr_v00","egfr_long")]
save(dat,file="./Data/dat.RData")


### MIXED LINEAR MODELS ###
###########################

#https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/

### ANALYSES ADJUSTED FOR DIABETES / HCHOL / HTA ###

z<-qnorm(1-0.05/2)
se <- function(x) sqrt(var(x) / length(x))

vars00<-c("Glucose, mg/dL\n(predicted values, 95% CI)",
          "Hb1Ac, %\n(predicted values, 95% CI)",
          "LDL cholesterol, mg/dL\n(predicted values, 95% CI)",
          "HDL cholesterol, mg/dL\n(predicted values, 95% CI)",
          "Triglycerides, mg/dL\n(predicted values, 95% CI)",
          "Systolic blood pressure, mmHg\n(predicted values, 95% CI)",
          "Diastolic blood pressure, mmHg\n(predicted values, 95% CI)")
vars01<-c("glucosa","hb1ac","ldl","hdl","trigli","pas","pad","imc","cintura","creat")
vars02<-c("glucosa_basal","hb1ac_basal","ldl_basal","hdl_basal","trigli_basal","pas_basal","pad_basal")
vars03<-c("diabetes00","diabetes00","hipercolest00","hipercolest00","hipercolest00","hta00","hta00")

vars04<-c("0","1","2","3","4","5","6","7","8","9",
          "10","11","12","13","14","15","16","17","18","19",
          "20","21","22","23","24","25","26","27","28","29",
          "30","31","32","33","34","35","36")
#vars04<-c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9",
#          "1","1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9",
#          "2","2.1","2.2","2.3","2.4","2.5","2.6","2.7","2.8","2.9",
#          "3","3.1","3.2","3.3","3.4","3.5","3.6","3.7","3.8","3.9",
#          "4","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","4.9","5")

tab<-NULL
tab2<-NULL
for(i in 1:length(vars01))
  
{
  load("./Data/dat_long.RData")
  dat_long<-subset2(dat_long,"!is.na(dat_long[,vars02[i]])
                    & !is.na(dat_long$kcalp_long) & !is.na(dat_long$score_long) & !is.na(dat_long$escolar) & !is.na(dat_long$grup_int)")
  xxx<-dat_long[,c("id",vars01[i],"kcalp_cat","score_long","segyears","edad","sexo","grup_int","escolar","tabaco00","kcalday_long",vars03[i])]
  names(xxx)<-c("id","variable","group","ajus01","seg","edad","sexo","grup_int","escolar","tabaco","kcal","cov")
  xxx$id2<-as.factor(xxx$id)
  
  mod_lin<-lme(variable~seg*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal+cov, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_par<-lme(variable~bs(seg,df=4)+group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal+cov, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_gam<-lme(variable~bs(seg,df=4)*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal+cov, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  
  pval_group<-pval_guapa(summary(mod_gam)$tTable[6,5])
  pval_time_group<-pval_guapa(lrtest(mod_par,mod_gam)[2,5])
  pval_lrtest<-pval_guapa(lrtest(mod_lin,mod_gam)[2,5])
  
  gam_predict<-expand.grid(group=factor(c(0,1,2)),
                           seg=seq(0,36,by=1),
                           id=unique(dat_long$id))
  xxx<-as.data.frame(unique(cbind(dat_long$id,dat_long$kcalp_cat,dat_long$score_long,dat_long$edad,dat_long$sexo,dat_long$grup_int,dat_long$escolar,
                                  dat_long$tabaco,dat_long$kcalday_long,dat_long[,vars03[i]])))
  names(xxx)<-c("id","group_ok","ajus01","edad","sexo","grup_int","escolar","tabaco","kcal","cov")
  gam_predict<-merge2(gam_predict,xxx,by.id=c("id"),all.x=TRUE,sort=FALSE)
  gam_predict$group<-as.numeric(gam_predict$group)-1
  gam_predict$group_ok<-as.numeric(gam_predict$group_ok)
  gam_predict$id2<-as.factor(gam_predict$id)
  
  gam_predict<-gam_predict %>% 
    mutate(fit=predict(mod_gam,gam_predict,level=0,type="response"))
  
  gam_predict<-gam_predict[gam_predict$group==gam_predict$group_ok,]
  gam_predict0<-gam_predict[gam_predict$group==0,c("id","seg","fit")]
  gam_predict1<-gam_predict[gam_predict$group==1,c("id","seg","fit")]
  gam_predict2<-gam_predict[gam_predict$group==2,c("id","seg","fit")]
  
  group0<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict0[gam_predict0$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group0<-as.data.frame(rbind(group0,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group0)<-c("seg","group0_fit","group0_lo","group0_hi")
  
  group1<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict1[gam_predict1$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group1<-as.data.frame(rbind(group1,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group1)<-c("seg","group1_fit","group1_lo","group1_hi")
  
  group2<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict2[gam_predict2$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group2<-as.data.frame(rbind(group2,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group2)<-c("seg","group2_fit","group2_lo","group2_hi")
  
  plot.data<-merge2(group0,group1,by.id=c("seg"),all.x=TRUE,sort=FALSE)
  plot.data<-merge2(plot.data,group2,by.id=c("seg"),all.x=TRUE,sort=FALSE)
  #leg<-paste("Time*group interaction p-value",pval_time_group,
  #           "\nNon-linearity p-value",pval_lrtest,sep="")
  maxx1<-with(plot.data,ifelse(max(group0_hi,na.rm=TRUE)>max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),
                               ifelse(max(group0_hi,na.rm=TRUE)<max(group1_hi,na.rm=TRUE),max(group1_hi,na.rm=TRUE),
                                      ifelse(max(group0_hi,na.rm=TRUE)==max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),NA))))
  maxx2<-with(plot.data,ifelse(max(group0_hi,na.rm=TRUE)>max(group2_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),
                               ifelse(max(group0_hi,na.rm=TRUE)<max(group2_hi,na.rm=TRUE),max(group2_hi,na.rm=TRUE),
                                      ifelse(max(group0_hi,na.rm=TRUE)==max(group2_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),NA))))
  maxx<-ifelse(maxx1>maxx2,maxx1,
               ifelse(maxx1<maxx2,maxx2,
                      ifelse(maxx1==maxx2,maxx1,NA)))
  
  figure<-ggplot(data=plot.data, aes_string(x='seg', y='group0_fit')) + 
    geom_ribbon(aes_string(ymin='group0_lo', ymax='group0_hi'), alpha=0.25, fill="#2b2b2b") +
    geom_line(aes_string(x='seg', y='group0_fit'), color='#2b2b2b') + 
    geom_ribbon(aes_string(ymin='group1_lo', ymax='group1_hi'), alpha=0.25, fill="#b023fc") +
    geom_line(aes_string(x='seg', y='group1_fit'), color='#b023fc') +
    geom_ribbon(aes_string(ymin='group2_lo', ymax='group2_hi'), alpha=0.25, fill="#3093fc") +
    geom_line(aes_string(x='seg', y='group2_fit'), color='#3093fc') +
    theme_bw() +
    scale_x_continuous(expand=c(0,0), breaks=c(0,6,12,18,24,30,36)) +
    labs(x=c("Follow-up time (months)"),y=vars00[i]) +
    #annotate("text", x=max(plot.data$seg)*0.95, y=maxx, label=leg, vjust=1, hjust=1, size=4.5) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm"))
  
  namefile<-paste("./Splines_repeat/",vars01[i],"_kcalp.jpg",sep="")
  ggsave(filename=namefile, units="px", width=8000, height=8000, dpi=1200, bg="transparent")
  par(las=1,cex=1.2,bty="n",lheight=0.9)
  figure
  dev.off()
  
  gam_predictx<-spread(gam_predict[gam_predict$seg==0 | gam_predict$seg==6 | gam_predict$seg==12 | gam_predict$seg==18 | gam_predict$seg==24 | gam_predict$seg==30 | gam_predict$seg==36,
                                   c("id","seg","fit","group")], seg, fit)
  names(gam_predictx)<-c("id","group","v00","v06","v12","v18","v24","v30","v36")
  
  mod00<-lm(v00~as.factor(group),data=gam_predictx)
  mod06<-lm(v06~as.factor(group),data=gam_predictx)
  mod12<-lm(v12~as.factor(group),data=gam_predictx)
  mod18<-lm(v18~as.factor(group),data=gam_predictx)
  mod24<-lm(v24~as.factor(group),data=gam_predictx)
  mod30<-lm(v30~as.factor(group),data=gam_predictx)
  mod36<-lm(v36~as.factor(group),data=gam_predictx)
  
  class_g1<-c("<20 vs. 20-29.9")
  class_g2<-c(">=30 vs. 20-29.9")
  beta00_g1<-intervals(mod00)[2,1]
  lo00_g1<-intervals(mod00)[2,2]
  hi00_g1<-intervals(mod00)[2,3]
  coef00_g1<-ic_guapa(guapa(beta00_g1),guapa(lo00_g1),guapa(hi00_g1))
  pval00_g1<-pval_guapa(intervals(mod00)[2,4])
  beta06_g1<-intervals(mod06)[2,1]
  lo06_g1<-intervals(mod06)[2,2]
  hi06_g1<-intervals(mod06)[2,3]
  coef06_g1<-ic_guapa(guapa(beta06_g1),guapa(lo06_g1),guapa(hi06_g1))
  pval06_g1<-pval_guapa(intervals(mod06)[2,4])
  beta12_g1<-intervals(mod12)[2,1]
  lo12_g1<-intervals(mod12)[2,2]
  hi12_g1<-intervals(mod12)[2,3]
  coef12_g1<-ic_guapa(guapa(beta12_g1),guapa(lo12_g1),guapa(hi12_g1))
  pval12_g1<-pval_guapa(intervals(mod12)[2,4])
  beta18_g1<-intervals(mod18)[2,1]
  lo18_g1<-intervals(mod18)[2,2]
  hi18_g1<-intervals(mod18)[2,3]
  coef18_g1<-ic_guapa(guapa(beta18_g1),guapa(lo18_g1),guapa(hi18_g1))
  pval18_g1<-pval_guapa(intervals(mod18)[2,4])
  beta24_g1<-intervals(mod24)[2,1]
  lo24_g1<-intervals(mod24)[2,2]
  hi24_g1<-intervals(mod24)[2,3]
  coef24_g1<-ic_guapa(guapa(beta24_g1),guapa(lo24_g1),guapa(hi24_g1))
  pval24_g1<-pval_guapa(intervals(mod24)[2,4])
  beta30_g1<-intervals(mod30)[2,1]
  lo30_g1<-intervals(mod30)[2,2]
  hi30_g1<-intervals(mod30)[2,3]
  coef30_g1<-ic_guapa(guapa(beta30_g1),guapa(lo30_g1),guapa(hi30_g1))
  pval30_g1<-pval_guapa(intervals(mod30)[2,4])
  beta36_g1<-intervals(mod36)[2,1]
  lo36_g1<-intervals(mod36)[2,2]
  hi36_g1<-intervals(mod36)[2,3]
  coef36_g1<-ic_guapa(guapa(beta36_g1),guapa(lo36_g1),guapa(hi36_g1))
  pval36_g1<-pval_guapa(intervals(mod36)[2,4])
  beta00_g2<-intervals(mod00)[3,1]
  lo00_g2<-intervals(mod00)[3,2]
  hi00_g2<-intervals(mod00)[3,3]
  coef00_g2<-ic_guapa(guapa(beta00_g2),guapa(lo00_g2),guapa(hi00_g2))
  pval00_g2<-pval_guapa(intervals(mod00)[3,4])
  beta06_g2<-intervals(mod06)[3,1]
  lo06_g2<-intervals(mod06)[3,2]
  hi06_g2<-intervals(mod06)[3,3]
  coef06_g2<-ic_guapa(guapa(beta06_g2),guapa(lo06_g2),guapa(hi06_g2))
  pval06_g2<-pval_guapa(intervals(mod06)[3,4])
  beta12_g2<-intervals(mod12)[3,1]
  lo12_g2<-intervals(mod12)[3,2]
  hi12_g2<-intervals(mod12)[3,3]
  coef12_g2<-ic_guapa(guapa(beta12_g2),guapa(lo12_g2),guapa(hi12_g2))
  pval12_g2<-pval_guapa(intervals(mod12)[3,4])
  beta18_g2<-intervals(mod18)[3,1]
  lo18_g2<-intervals(mod18)[3,2]
  hi18_g2<-intervals(mod18)[3,3]
  coef18_g2<-ic_guapa(guapa(beta18_g2),guapa(lo18_g2),guapa(hi18_g2))
  pval18_g2<-pval_guapa(intervals(mod18)[3,4])
  beta24_g2<-intervals(mod24)[3,1]
  lo24_g2<-intervals(mod24)[3,2]
  hi24_g2<-intervals(mod24)[3,3]
  coef24_g2<-ic_guapa(guapa(beta24_g2),guapa(lo24_g2),guapa(hi24_g2))
  pval24_g2<-pval_guapa(intervals(mod24)[3,4])
  beta30_g2<-intervals(mod30)[3,1]
  lo30_g2<-intervals(mod30)[3,2]
  hi30_g2<-intervals(mod30)[3,3]
  coef30_g2<-ic_guapa(guapa(beta30_g2),guapa(lo30_g2),guapa(hi30_g2))
  pval30_g2<-pval_guapa(intervals(mod30)[3,4])
  beta36_g2<-intervals(mod36)[3,1]
  lo36_g2<-intervals(mod36)[3,2]
  hi36_g2<-intervals(mod36)[3,3]
  coef36_g2<-ic_guapa(guapa(beta36_g2),guapa(lo36_g2),guapa(hi36_g2))
  pval36_g2<-pval_guapa(intervals(mod36)[3,4])
  tab<-rbind(tab,cbind(beta00_g1,lo00_g1,hi00_g1,beta06_g1,lo06_g1,hi06_g1,beta12_g1,lo12_g1,hi12_g1,beta18_g1,lo18_g1,hi18_g1,
                       beta24_g1,lo24_g1,hi24_g1,beta30_g1,lo30_g1,hi30_g1,beta36_g1,lo36_g1,hi36_g1))
  tab<-rbind(tab,cbind(beta00_g2,lo00_g2,hi00_g2,beta06_g2,lo06_g2,hi06_g2,beta12_g2,lo12_g2,hi12_g2,beta18_g2,lo18_g2,hi18_g2,
                       beta24_g2,lo24_g2,hi24_g2,beta30_g2,lo30_g2,hi30_g2,beta36_g2,lo36_g2,hi36_g2))
  tab2<-rbind(tab2,cbind(class_g1,pval_group,pval_time_group,pval_lrtest,
                         coef00_g1,pval00_g1,coef06_g1,pval06_g1,coef12_g1,pval12_g1,coef18_g1,pval18_g1,coef24_g1,pval24_g1,coef30_g1,pval30_g1,coef36_g1,pval36_g1))
  tab2<-rbind(tab2,cbind(class_g2,pval_group,pval_time_group,pval_lrtest,
                         coef00_g2,pval00_g2,coef06_g2,pval06_g2,coef12_g2,pval12_g2,coef18_g2,pval18_g2,coef24_g2,pval24_g2,coef30_g2,pval30_g2,coef36_g2,pval36_g2))
  
  
  load("./Data/dat_long.RData")
  dat_long<-subset2(dat_long,"!is.na(dat_long[,vars02[i]]) & !is.na(dat_long$score_cat) & !is.na(dat_long$kcalp_long)
                    & !is.na(dat_long$score_long) & !is.na(dat_long$escolar) & !is.na(dat_long$grup_int)")
  score_group<-paste("<",guapa(as.numeric(quantile(dat_long$score_long,prob=c(0.33333),na.rm=TRUE,type=6))[1]),
                     " vs. >=",guapa(as.numeric(quantile(dat_long$score_long,prob=c(0.66666),na.rm=TRUE,type=6))[1]),sep="")
  xxx<-dat_long[,c("id",vars01[i],"score_cat","kcalp_long","segyears","edad","sexo","grup_int","escolar","tabaco00","kcalday_long",vars03[i])]
  names(xxx)<-c("id","variable","group","ajus01","seg","edad","sexo","grup_int","escolar","tabaco","kcal","cov")
  xxx$id2<-as.factor(xxx$id)
  
  mod_lin<-lme(variable~seg*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal+cov, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_par<-lme(variable~bs(seg,df=4)+group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal+cov, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_gam<-lme(variable~bs(seg,df=4)*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal+cov, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  
  pval_group<-pval_guapa(summary(mod_gam)$tTable[6,5])
  pval_time_group<-pval_guapa(lrtest(mod_par,mod_gam)[2,5])
  pval_lrtest<-pval_guapa(lrtest(mod_lin,mod_gam)[2,5])
  
  gam_predict<-expand.grid(group=factor(c(0,1)),
                           seg=seq(0,36,by=1),
                           id=unique(dat_long$id))
  xxx<-as.data.frame(unique(cbind(dat_long$id,dat_long$score_cat,dat_long$kcalp_long,dat_long$edad,dat_long$sexo,dat_long$grup_int,dat_long$escolar,
                                  dat_long$tabaco,dat_long$kcalday_long,dat_long[,vars03[i]])))
  names(xxx)<-c("id","group_ok","ajus01","edad","sexo","grup_int","escolar","tabaco","kcal","cov")
  gam_predict<-merge2(gam_predict,xxx,by.id=c("id"),all.x=TRUE,sort=FALSE)
  gam_predict$group<-as.numeric(gam_predict$group)-1
  gam_predict$group_ok<-as.numeric(gam_predict$group_ok)
  gam_predict$id2<-as.factor(gam_predict$id)
  
  gam_predict<-gam_predict %>% 
    mutate(fit=predict(mod_gam,gam_predict,level=0,type="response"))
  
  gam_predict<-gam_predict[gam_predict$group==gam_predict$group_ok,]
  gam_predict0<-gam_predict[gam_predict$group==0,c("id","seg","fit")]
  gam_predict1<-gam_predict[gam_predict$group==1,c("id","seg","fit")]
  
  group0<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict0[gam_predict0$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group0<-as.data.frame(rbind(group0,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group0)<-c("seg","group0_fit","group0_lo","group0_hi")
  
  group1<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict1[gam_predict1$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group1<-as.data.frame(rbind(group1,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group1)<-c("seg","group1_fit","group1_lo","group1_hi")
  
  plot.data<-merge2(group0,group1,by.id=c("seg"),all.x=TRUE,sort=FALSE)
  #leg<-paste("Time*group interaction p-value",pval_time_group,
  #           "\nNon-linearity p-value",pval_lrtest,sep="")
  maxx<-with(plot.data,ifelse(max(group0_hi,na.rm=TRUE)>max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),
                              ifelse(max(group0_hi,na.rm=TRUE)<max(group1_hi,na.rm=TRUE),max(group1_hi,na.rm=TRUE),
                                     ifelse(max(group0_hi,na.rm=TRUE)==max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),NA))))
  
  figure<-ggplot(data=plot.data, aes_string(x='seg', y='group0_fit')) + 
    geom_ribbon(aes_string(ymin='group0_lo', ymax='group0_hi'), alpha=0.2, fill="#2b2b2b") +
    geom_line(aes_string(x='seg', y='group0_fit'), color='#2b2b2b') + 
    geom_ribbon(aes_string(ymin='group1_lo', ymax='group1_hi'), alpha=0.2, fill="#d1024f") +
    geom_line(aes_string(x='seg', y='group1_fit'), color='#d1024f') +
    theme_bw() +
    scale_x_continuous(expand=c(0,0), breaks=c(0,6,12,18,24,30,36)) +
    labs(x=c("Follow-up time (months)"),y=vars00[i]) +
    #annotate("text", x=max(plot.data$seg)*0.95, y=maxx, label=leg, vjust=1, hjust=1, size=4.5) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm"))
  
  namefile<-paste("./Splines_repeat/",vars01[i],"_score.jpg",sep="")
  ggsave(filename=namefile, units="px", width=8000, height=8000, dpi=1200, bg="transparent")
  par(las=1,cex=1.2,bty="n",lheight=0.9)
  figure
  dev.off()
  
  gam_predictx<-spread(gam_predict[gam_predict$seg==0 | gam_predict$seg==6 | gam_predict$seg==12 | gam_predict$seg==18 | gam_predict$seg==24 | gam_predict$seg==30 | gam_predict$seg==36,
                                   c("id","seg","fit","group")], seg, fit)
  names(gam_predictx)<-c("id","group","v00","v06","v12","v18","v24","v30","v36")
  
  mod00<-lm(v00~as.factor(group),data=gam_predictx)
  mod06<-lm(v06~as.factor(group),data=gam_predictx)
  mod12<-lm(v12~as.factor(group),data=gam_predictx)
  mod18<-lm(v18~as.factor(group),data=gam_predictx)
  mod24<-lm(v24~as.factor(group),data=gam_predictx)
  mod30<-lm(v30~as.factor(group),data=gam_predictx)
  mod36<-lm(v36~as.factor(group),data=gam_predictx)
  
  beta00<-intervals(mod00)[2,1]
  lo00<-intervals(mod00)[2,2]
  hi00<-intervals(mod00)[2,3]
  coef00<-ic_guapa(guapa(beta00),guapa(lo00),guapa(hi00))
  pval00<-pval_guapa(intervals(mod00)[2,4])
  beta06<-intervals(mod06)[2,1]
  lo06<-intervals(mod06)[2,2]
  hi06<-intervals(mod06)[2,3]
  coef06<-ic_guapa(guapa(beta06),guapa(lo06),guapa(hi06))
  pval06<-pval_guapa(intervals(mod06)[2,4])
  beta12<-intervals(mod12)[2,1]
  lo12<-intervals(mod12)[2,2]
  hi12<-intervals(mod12)[2,3]
  coef12<-ic_guapa(guapa(beta12),guapa(lo12),guapa(hi12))
  pval12<-pval_guapa(intervals(mod12)[2,4])
  beta18<-intervals(mod18)[2,1]
  lo18<-intervals(mod18)[2,2]
  hi18<-intervals(mod18)[2,3]
  coef18<-ic_guapa(guapa(beta18),guapa(lo18),guapa(hi18))
  pval18<-pval_guapa(intervals(mod18)[2,4])
  beta24<-intervals(mod24)[2,1]
  lo24<-intervals(mod24)[2,2]
  hi24<-intervals(mod24)[2,3]
  coef24<-ic_guapa(guapa(beta24),guapa(lo24),guapa(hi24))
  pval24<-pval_guapa(intervals(mod24)[2,4])
  beta30<-intervals(mod30)[2,1]
  lo30<-intervals(mod30)[2,2]
  hi30<-intervals(mod30)[2,3]
  coef30<-ic_guapa(guapa(beta30),guapa(lo30),guapa(hi30))
  pval30<-pval_guapa(intervals(mod30)[2,4])
  beta36<-intervals(mod36)[2,1]
  lo36<-intervals(mod36)[2,2]
  hi36<-intervals(mod36)[2,3]
  coef36<-ic_guapa(guapa(beta36),guapa(lo36),guapa(hi36))
  pval36<-pval_guapa(intervals(mod36)[2,4])
  tab<-rbind(tab,cbind(beta00,lo00,hi00,beta06,lo06,hi06,beta12,lo12,hi12,beta18,lo18,hi18,beta24,lo24,hi24,beta30,lo30,hi30,beta36,lo36,hi36))
  tab2<-rbind(tab2,cbind(score_group,pval_group,pval_time_group,pval_lrtest,
                         coef00,pval00,coef06,pval06,coef12,pval12,coef18,pval18,coef24,pval24,coef30,pval30,coef36,pval36))
  
  
  }

rownames(tab)<-paste(c("glucosa","glucosa","glucosa",
                       "hb1ac","hb1ac","hb1ac",
                       "ldl","ldl","ldl",
                       "hdl","hdl","hdl",
                       "trigli","trigli","trigli",
                       "pas","pas","pas",
                       "pad","pad","pad")
                     ,c("kcalp_lo","kcalp_hi","score"),sep="_")
rownames(tab2)<-paste(c("glucosa","glucosa","glucosa",
                        "hb1ac","hb1ac","hb1ac",
                        "ldl","ldl","ldl",
                        "hdl","hdl","hdl",
                        "trigli","trigli","trigli",
                        "pas","pas","pas",
                        "pad","pad","pad")
                      ,c("kcalp_lo","kcalp_hi","score"),sep="_")
write.table(tab,file="./Forest_repeat/forest01.csv",sep=";",col.names=NA)
write.table(tab2,file="./Forest_repeat/coef01.csv",sep=";",col.names=NA)


### ANALYSES NOT ADJUSTED FOR DIABETES / HCHOL / HTA ###

vars00<-c("Body mass index, kg/m2\n(predicted values, 95% CI)",
          "Waist circumference, cm\n(predicted values, 95% CI)",
          "Glomerular filtration rate, mL/min/1.73 m2\n(predicted values, 95% CI)")
vars01<-c("imc","cintura","egfr")
vars02<-c("imc_basal","cintura_basal","egfr_basal")

vars04<-c("0","1","2","3","4","5","6","7","8","9",
          "10","11","12","13","14","15","16","17","18","19",
          "20","21","22","23","24","25","26","27","28","29",
          "30","31","32","33","34","35","36")

tab<-NULL
tab2<-NULL
for(i in 1:length(vars01))
  
{
  load("./Data/dat_long.RData")
  dat_long<-subset2(dat_long,"!is.na(dat_long[,vars02[i]])
                    & !is.na(dat_long$kcalp_long) & !is.na(dat_long$score_long) & !is.na(dat_long$escolar) & !is.na(dat_long$grup_int)")
  xxx<-dat_long[,c("id",vars01[i],"kcalp_cat","score_long","segyears","edad","sexo","grup_int","escolar","tabaco00","kcalday_long")]
  names(xxx)<-c("id","variable","group","ajus01","seg","edad","sexo","grup_int","escolar","tabaco","kcal")
  xxx$id2<-as.factor(xxx$id)
  
  mod_lin<-lme(variable~seg*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_par<-lme(variable~bs(seg,df=4)+group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_gam<-lme(variable~bs(seg,df=4)*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  
  pval_group<-pval_guapa(summary(mod_gam)$tTable[6,5])
  pval_time_group<-pval_guapa(lrtest(mod_par,mod_gam)[2,5])
  pval_lrtest<-pval_guapa(lrtest(mod_lin,mod_gam)[2,5])
  
  gam_predict<-expand.grid(group=factor(c(0,1,2)),
                           seg=seq(0,36,by=1),
                           id=unique(dat_long$id))
  xxx<-as.data.frame(unique(cbind(dat_long$id,dat_long$kcalp_cat,dat_long$score_long,dat_long$edad,dat_long$sexo,dat_long$grup_int,dat_long$escolar,
                                  dat_long$tabaco,dat_long$kcalday_long)))
  names(xxx)<-c("id","group_ok","ajus01","edad","sexo","grup_int","escolar","tabaco","kcal")
  gam_predict<-merge2(gam_predict,xxx,by.id=c("id"),all.x=TRUE,sort=FALSE)
  gam_predict$group<-as.numeric(gam_predict$group)-1
  gam_predict$group_ok<-as.numeric(gam_predict$group_ok)
  gam_predict$id2<-as.factor(gam_predict$id)
  
  gam_predict<-gam_predict %>% 
    mutate(fit=predict(mod_gam,gam_predict,level=0,type="response"))
  
  gam_predict<-gam_predict[gam_predict$group==gam_predict$group_ok,]
  gam_predict0<-gam_predict[gam_predict$group==0,c("id","seg","fit")]
  gam_predict1<-gam_predict[gam_predict$group==1,c("id","seg","fit")]
  gam_predict2<-gam_predict[gam_predict$group==2,c("id","seg","fit")]
  
  group0<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict0[gam_predict0$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group0<-as.data.frame(rbind(group0,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group0)<-c("seg","group0_fit","group0_lo","group0_hi")
  
  group1<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict1[gam_predict1$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group1<-as.data.frame(rbind(group1,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group1)<-c("seg","group1_fit","group1_lo","group1_hi")
  
  group2<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict2[gam_predict2$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group2<-as.data.frame(rbind(group2,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group2)<-c("seg","group2_fit","group2_lo","group2_hi")
  
  plot.data<-merge2(group0,group1,by.id=c("seg"),all.x=TRUE,sort=FALSE)
  plot.data<-merge2(plot.data,group2,by.id=c("seg"),all.x=TRUE,sort=FALSE)
  #leg<-paste("Time*group interaction p-value",pval_time_group,
  #           "\nNon-linearity p-value",pval_lrtest,sep="")
  maxx1<-with(plot.data,ifelse(max(group0_hi,na.rm=TRUE)>max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),
                               ifelse(max(group0_hi,na.rm=TRUE)<max(group1_hi,na.rm=TRUE),max(group1_hi,na.rm=TRUE),
                                      ifelse(max(group0_hi,na.rm=TRUE)==max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),NA))))
  maxx2<-with(plot.data,ifelse(max(group0_hi,na.rm=TRUE)>max(group2_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),
                               ifelse(max(group0_hi,na.rm=TRUE)<max(group2_hi,na.rm=TRUE),max(group2_hi,na.rm=TRUE),
                                      ifelse(max(group0_hi,na.rm=TRUE)==max(group2_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),NA))))
  maxx<-ifelse(maxx1>maxx2,maxx1,
               ifelse(maxx1<maxx2,maxx2,
                      ifelse(maxx1==maxx2,maxx1,NA)))
  
  figure<-ggplot(data=plot.data, aes_string(x='seg', y='group0_fit')) + 
    geom_ribbon(aes_string(ymin='group0_lo', ymax='group0_hi'), alpha=0.25, fill="#2b2b2b") +
    geom_line(aes_string(x='seg', y='group0_fit'), color='#2b2b2b') + 
    geom_ribbon(aes_string(ymin='group1_lo', ymax='group1_hi'), alpha=0.25, fill="#b023fc") +
    geom_line(aes_string(x='seg', y='group1_fit'), color='#b023fc') +
    geom_ribbon(aes_string(ymin='group2_lo', ymax='group2_hi'), alpha=0.25, fill="#3093fc") +
    geom_line(aes_string(x='seg', y='group2_fit'), color='#3093fc') +
    theme_bw() +
    scale_x_continuous(expand=c(0,0), breaks=c(0,6,12,18,24,30,36)) +
    labs(x=c("Follow-up time (months)"),y=vars00[i]) +
    #annotate("text", x=max(plot.data$seg)*0.95, y=maxx, label=leg, vjust=1, hjust=1, size=4.5) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm"))
  
  namefile<-paste("./Splines_repeat/",vars01[i],"_kcalp.jpg",sep="")
  ggsave(filename=namefile, units="px", width=8000, height=8000, dpi=1200, bg="transparent")
  par(las=1,cex=1.2,bty="n",lheight=0.9)
  figure
  dev.off()
  
  gam_predictx<-spread(gam_predict[gam_predict$seg==0 | gam_predict$seg==6 | gam_predict$seg==12 | gam_predict$seg==18 | gam_predict$seg==24 | gam_predict$seg==30 | gam_predict$seg==36,
                                   c("id","seg","fit","group")], seg, fit)
  names(gam_predictx)<-c("id","group","v00","v06","v12","v18","v24","v30","v36")
  
  mod00<-lm(v00~as.factor(group),data=gam_predictx)
  mod06<-lm(v06~as.factor(group),data=gam_predictx)
  mod12<-lm(v12~as.factor(group),data=gam_predictx)
  mod18<-lm(v18~as.factor(group),data=gam_predictx)
  mod24<-lm(v24~as.factor(group),data=gam_predictx)
  mod30<-lm(v30~as.factor(group),data=gam_predictx)
  mod36<-lm(v36~as.factor(group),data=gam_predictx)
  
  class_g1<-c("<20 vs. 20-29.9")
  class_g2<-c(">=30 vs. 20-29.9")
  beta00_g1<-intervals(mod00)[2,1]
  lo00_g1<-intervals(mod00)[2,2]
  hi00_g1<-intervals(mod00)[2,3]
  coef00_g1<-ic_guapa(guapa(beta00_g1),guapa(lo00_g1),guapa(hi00_g1))
  pval00_g1<-pval_guapa(intervals(mod00)[2,4])
  beta06_g1<-intervals(mod06)[2,1]
  lo06_g1<-intervals(mod06)[2,2]
  hi06_g1<-intervals(mod06)[2,3]
  coef06_g1<-ic_guapa(guapa(beta06_g1),guapa(lo06_g1),guapa(hi06_g1))
  pval06_g1<-pval_guapa(intervals(mod06)[2,4])
  beta12_g1<-intervals(mod12)[2,1]
  lo12_g1<-intervals(mod12)[2,2]
  hi12_g1<-intervals(mod12)[2,3]
  coef12_g1<-ic_guapa(guapa(beta12_g1),guapa(lo12_g1),guapa(hi12_g1))
  pval12_g1<-pval_guapa(intervals(mod12)[2,4])
  beta18_g1<-intervals(mod18)[2,1]
  lo18_g1<-intervals(mod18)[2,2]
  hi18_g1<-intervals(mod18)[2,3]
  coef18_g1<-ic_guapa(guapa(beta18_g1),guapa(lo18_g1),guapa(hi18_g1))
  pval18_g1<-pval_guapa(intervals(mod18)[2,4])
  beta24_g1<-intervals(mod24)[2,1]
  lo24_g1<-intervals(mod24)[2,2]
  hi24_g1<-intervals(mod24)[2,3]
  coef24_g1<-ic_guapa(guapa(beta24_g1),guapa(lo24_g1),guapa(hi24_g1))
  pval24_g1<-pval_guapa(intervals(mod24)[2,4])
  beta30_g1<-intervals(mod30)[2,1]
  lo30_g1<-intervals(mod30)[2,2]
  hi30_g1<-intervals(mod30)[2,3]
  coef30_g1<-ic_guapa(guapa(beta30_g1),guapa(lo30_g1),guapa(hi30_g1))
  pval30_g1<-pval_guapa(intervals(mod30)[2,4])
  beta36_g1<-intervals(mod36)[2,1]
  lo36_g1<-intervals(mod36)[2,2]
  hi36_g1<-intervals(mod36)[2,3]
  coef36_g1<-ic_guapa(guapa(beta36_g1),guapa(lo36_g1),guapa(hi36_g1))
  pval36_g1<-pval_guapa(intervals(mod36)[2,4])
  beta00_g2<-intervals(mod00)[3,1]
  lo00_g2<-intervals(mod00)[3,2]
  hi00_g2<-intervals(mod00)[3,3]
  coef00_g2<-ic_guapa(guapa(beta00_g2),guapa(lo00_g2),guapa(hi00_g2))
  pval00_g2<-pval_guapa(intervals(mod00)[3,4])
  beta06_g2<-intervals(mod06)[3,1]
  lo06_g2<-intervals(mod06)[3,2]
  hi06_g2<-intervals(mod06)[3,3]
  coef06_g2<-ic_guapa(guapa(beta06_g2),guapa(lo06_g2),guapa(hi06_g2))
  pval06_g2<-pval_guapa(intervals(mod06)[3,4])
  beta12_g2<-intervals(mod12)[3,1]
  lo12_g2<-intervals(mod12)[3,2]
  hi12_g2<-intervals(mod12)[3,3]
  coef12_g2<-ic_guapa(guapa(beta12_g2),guapa(lo12_g2),guapa(hi12_g2))
  pval12_g2<-pval_guapa(intervals(mod12)[3,4])
  beta18_g2<-intervals(mod18)[3,1]
  lo18_g2<-intervals(mod18)[3,2]
  hi18_g2<-intervals(mod18)[3,3]
  coef18_g2<-ic_guapa(guapa(beta18_g2),guapa(lo18_g2),guapa(hi18_g2))
  pval18_g2<-pval_guapa(intervals(mod18)[3,4])
  beta24_g2<-intervals(mod24)[3,1]
  lo24_g2<-intervals(mod24)[3,2]
  hi24_g2<-intervals(mod24)[3,3]
  coef24_g2<-ic_guapa(guapa(beta24_g2),guapa(lo24_g2),guapa(hi24_g2))
  pval24_g2<-pval_guapa(intervals(mod24)[3,4])
  beta30_g2<-intervals(mod30)[3,1]
  lo30_g2<-intervals(mod30)[3,2]
  hi30_g2<-intervals(mod30)[3,3]
  coef30_g2<-ic_guapa(guapa(beta30_g2),guapa(lo30_g2),guapa(hi30_g2))
  pval30_g2<-pval_guapa(intervals(mod30)[3,4])
  beta36_g2<-intervals(mod36)[3,1]
  lo36_g2<-intervals(mod36)[3,2]
  hi36_g2<-intervals(mod36)[3,3]
  coef36_g2<-ic_guapa(guapa(beta36_g2),guapa(lo36_g2),guapa(hi36_g2))
  pval36_g2<-pval_guapa(intervals(mod36)[3,4])
  tab<-rbind(tab,cbind(beta00_g1,lo00_g1,hi00_g1,beta06_g1,lo06_g1,hi06_g1,beta12_g1,lo12_g1,hi12_g1,beta18_g1,lo18_g1,hi18_g1,
                       beta24_g1,lo24_g1,hi24_g1,beta30_g1,lo30_g1,hi30_g1,beta36_g1,lo36_g1,hi36_g1))
  tab<-rbind(tab,cbind(beta00_g2,lo00_g2,hi00_g2,beta06_g2,lo06_g2,hi06_g2,beta12_g2,lo12_g2,hi12_g2,beta18_g2,lo18_g2,hi18_g2,
                       beta24_g2,lo24_g2,hi24_g2,beta30_g2,lo30_g2,hi30_g2,beta36_g2,lo36_g2,hi36_g2))
  tab2<-rbind(tab2,cbind(class_g1,pval_group,pval_time_group,pval_lrtest,
                         coef00_g1,pval00_g1,coef06_g1,pval06_g1,coef12_g1,pval12_g1,coef18_g1,pval18_g1,coef24_g1,pval24_g1,coef30_g1,pval30_g1,coef36_g1,pval36_g1))
  tab2<-rbind(tab2,cbind(class_g2,pval_group,pval_time_group,pval_lrtest,
                         coef00_g2,pval00_g2,coef06_g2,pval06_g2,coef12_g2,pval12_g2,coef18_g2,pval18_g2,coef24_g2,pval24_g2,coef30_g2,pval30_g2,coef36_g2,pval36_g2))
  
  
  load("./Data/dat_long.RData")
  dat_long<-subset2(dat_long,"!is.na(dat_long[,vars02[i]]) & !is.na(dat_long$score_cat) & !is.na(dat_long$kcalp_long)
                    & !is.na(dat_long$score_long) & !is.na(dat_long$escolar) & !is.na(dat_long$grup_int)")
  score_group<-paste("<",guapa(as.numeric(quantile(dat_long$score_long,prob=c(0.33333),na.rm=TRUE,type=6))[1]),
                     " vs. >=",guapa(as.numeric(quantile(dat_long$score_long,prob=c(0.66666),na.rm=TRUE,type=6))[1]),sep="")
  xxx<-dat_long[,c("id",vars01[i],"score_cat","kcalp_long","segyears","edad","sexo","grup_int","escolar","tabaco00","kcalday_long")]
  names(xxx)<-c("id","variable","group","ajus01","seg","edad","sexo","grup_int","escolar","tabaco","kcal")
  xxx$id2<-as.factor(xxx$id)
  
  mod_lin<-lme(variable~seg*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_par<-lme(variable~bs(seg,df=4)+group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  mod_gam<-lme(variable~bs(seg,df=4)*group+ajus01+edad+sexo+grup_int+as.factor(escolar)+tabaco+kcal, 
               random = ~1 | id2, 
               correlation = corCAR1(form = ~seg | id2), 
               control = lmeControl(opt="optim"), 
               data=xxx,
               method='REML',
               na.action=na.exclude)
  
  pval_group<-pval_guapa(summary(mod_gam)$tTable[6,5])
  pval_time_group<-pval_guapa(lrtest(mod_par,mod_gam)[2,5])
  pval_lrtest<-pval_guapa(lrtest(mod_lin,mod_gam)[2,5])
  
  gam_predict<-expand.grid(group=factor(c(0,1)),
                           seg=seq(0,36,by=1),
                           id=unique(dat_long$id))
  xxx<-as.data.frame(unique(cbind(dat_long$id,dat_long$score_cat,dat_long$kcalp_long,dat_long$edad,dat_long$sexo,dat_long$grup_int,dat_long$escolar,
                                  dat_long$tabaco,dat_long$kcalday_long)))
  names(xxx)<-c("id","group_ok","ajus01","edad","sexo","grup_int","escolar","tabaco","kcal")
  gam_predict<-merge2(gam_predict,xxx,by.id=c("id"),all.x=TRUE,sort=FALSE)
  gam_predict$group<-as.numeric(gam_predict$group)-1
  gam_predict$group_ok<-as.numeric(gam_predict$group_ok)
  gam_predict$id2<-as.factor(gam_predict$id)
  
  gam_predict<-gam_predict %>% 
    mutate(fit=predict(mod_gam,gam_predict,level=0,type="response"))
  
  gam_predict<-gam_predict[gam_predict$group==gam_predict$group_ok,]
  gam_predict0<-gam_predict[gam_predict$group==0,c("id","seg","fit")]
  gam_predict1<-gam_predict[gam_predict$group==1,c("id","seg","fit")]
  
  group0<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict0[gam_predict0$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group0<-as.data.frame(rbind(group0,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group0)<-c("seg","group0_fit","group0_lo","group0_hi")
  
  group1<-NULL
  for(j in 1:length(vars04))
  {
    xxx<-na.omit(gam_predict1[gam_predict1$seg==vars04[j],])
    meanx<-mean(xxx$fit)
    stdex<-se(xxx$fit)
    group1<-as.data.frame(rbind(group1,cbind(as.numeric(vars04[j]),meanx,meanx-z*stdex,meanx+z*stdex)))
  }
  names(group1)<-c("seg","group1_fit","group1_lo","group1_hi")
  
  plot.data<-merge2(group0,group1,by.id=c("seg"),all.x=TRUE,sort=FALSE)
  #leg<-paste("Time*group interaction p-value",pval_time_group,
  #           "\nNon-linearity p-value",pval_lrtest,sep="")
  maxx<-with(plot.data,ifelse(max(group0_hi,na.rm=TRUE)>max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),
                              ifelse(max(group0_hi,na.rm=TRUE)<max(group1_hi,na.rm=TRUE),max(group1_hi,na.rm=TRUE),
                                     ifelse(max(group0_hi,na.rm=TRUE)==max(group1_hi,na.rm=TRUE),max(group0_hi,na.rm=TRUE),NA))))
  
  figure<-ggplot(data=plot.data, aes_string(x='seg', y='group0_fit')) + 
    geom_ribbon(aes_string(ymin='group0_lo', ymax='group0_hi'), alpha=0.2, fill="#2b2b2b") +
    geom_line(aes_string(x='seg', y='group0_fit'), color='#2b2b2b') + 
    geom_ribbon(aes_string(ymin='group1_lo', ymax='group1_hi'), alpha=0.2, fill="#d1024f") +
    geom_line(aes_string(x='seg', y='group1_fit'), color='#d1024f') +
    theme_bw() +
    scale_x_continuous(expand=c(0,0), breaks=c(0,6,12,18,24,30,36)) +
    labs(x=c("Follow-up time (months)"),y=vars00[i]) +
    #annotate("text", x=max(plot.data$seg)*0.95, y=maxx, label=leg, vjust=1, hjust=1, size=4.5) +
    theme(axis.title.x = element_text(vjust=0.5, size=18, face="bold"), 
          axis.title.y = element_text(vjust=0.5, size=18, face="bold"),
          axis.text.x = element_text(size=16, colour = 'black'),
          axis.text.y = element_text(size=16, colour = 'black'),
          axis.ticks.x = element_line(colour = 'black'),
          axis.ticks.y = element_line(colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm"))
  
  namefile<-paste("./Splines_repeat/",vars01[i],"_score.jpg",sep="")
  ggsave(filename=namefile, units="px", width=8000, height=8000, dpi=1200, bg="transparent")
  par(las=1,cex=1.2,bty="n",lheight=0.9)
  figure
  dev.off()
  
  gam_predictx<-spread(gam_predict[gam_predict$seg==0 | gam_predict$seg==6 | gam_predict$seg==12 | gam_predict$seg==18 | gam_predict$seg==24 | gam_predict$seg==30 | gam_predict$seg==36,
                                   c("id","seg","fit","group")], seg, fit)
  names(gam_predictx)<-c("id","group","v00","v06","v12","v18","v24","v30","v36")
  
  mod00<-lm(v00~as.factor(group),data=gam_predictx)
  mod06<-lm(v06~as.factor(group),data=gam_predictx)
  mod12<-lm(v12~as.factor(group),data=gam_predictx)
  mod18<-lm(v18~as.factor(group),data=gam_predictx)
  mod24<-lm(v24~as.factor(group),data=gam_predictx)
  mod30<-lm(v30~as.factor(group),data=gam_predictx)
  mod36<-lm(v36~as.factor(group),data=gam_predictx)
  
  beta00<-intervals(mod00)[2,1]
  lo00<-intervals(mod00)[2,2]
  hi00<-intervals(mod00)[2,3]
  coef00<-ic_guapa(guapa(beta00),guapa(lo00),guapa(hi00))
  pval00<-pval_guapa(intervals(mod00)[2,4])
  beta06<-intervals(mod06)[2,1]
  lo06<-intervals(mod06)[2,2]
  hi06<-intervals(mod06)[2,3]
  coef06<-ic_guapa(guapa(beta06),guapa(lo06),guapa(hi06))
  pval06<-pval_guapa(intervals(mod06)[2,4])
  beta12<-intervals(mod12)[2,1]
  lo12<-intervals(mod12)[2,2]
  hi12<-intervals(mod12)[2,3]
  coef12<-ic_guapa(guapa(beta12),guapa(lo12),guapa(hi12))
  pval12<-pval_guapa(intervals(mod12)[2,4])
  beta18<-intervals(mod18)[2,1]
  lo18<-intervals(mod18)[2,2]
  hi18<-intervals(mod18)[2,3]
  coef18<-ic_guapa(guapa(beta18),guapa(lo18),guapa(hi18))
  pval18<-pval_guapa(intervals(mod18)[2,4])
  beta24<-intervals(mod24)[2,1]
  lo24<-intervals(mod24)[2,2]
  hi24<-intervals(mod24)[2,3]
  coef24<-ic_guapa(guapa(beta24),guapa(lo24),guapa(hi24))
  pval24<-pval_guapa(intervals(mod24)[2,4])
  beta30<-intervals(mod30)[2,1]
  lo30<-intervals(mod30)[2,2]
  hi30<-intervals(mod30)[2,3]
  coef30<-ic_guapa(guapa(beta30),guapa(lo30),guapa(hi30))
  pval30<-pval_guapa(intervals(mod30)[2,4])
  beta36<-intervals(mod36)[2,1]
  lo36<-intervals(mod36)[2,2]
  hi36<-intervals(mod36)[2,3]
  coef36<-ic_guapa(guapa(beta36),guapa(lo36),guapa(hi36))
  pval36<-pval_guapa(intervals(mod36)[2,4])
  tab<-rbind(tab,cbind(beta00,lo00,hi00,beta06,lo06,hi06,beta12,lo12,hi12,beta18,lo18,hi18,beta24,lo24,hi24,beta30,lo30,hi30,beta36,lo36,hi36))
  tab2<-rbind(tab2,cbind(score_group,pval_group,pval_time_group,pval_lrtest,
                         coef00,pval00,coef06,pval06,coef12,pval12,coef18,pval18,coef24,pval24,coef30,pval30,coef36,pval36))
  
  
}

rownames(tab)<-paste(c("imc","imc","imc",
                       "cintura","cintura","cintura",
                       "egfr","egfr","egfr")
                     ,c("kcalp_lo","kcalp_hi","score"),sep="_")
rownames(tab2)<-paste(c("imc","imc","imc",
                        "cintura","cintura","cintura",
                        "egfr","egfr","egfr")
                      ,c("kcalp_lo","kcalp_hi","score"),sep="_")
write.table(tab,file="./Forest_repeat/forest02.csv",sep=";",col.names=NA)
write.table(tab2,file="./Forest_repeat/coef02.csv",sep=";",col.names=NA)


coef01<-read.csv2("./Forest_repeat/coef01.csv",header=TRUE,sep=";",dec=".")
coef02<-read.csv2("./Forest_repeat/coef02.csv",header=TRUE,sep=";",dec=".")
coef<-rbind(coef01,coef02)
write.table(coef,file="./predvalues_coef.csv",sep=";",col.names=TRUE,row.names=FALSE)

forest01<-read.csv2("./Forest_repeat/forest01.csv",header=TRUE,sep=";",dec=".")
forest02<-read.csv2("./Forest_repeat/forest02.csv",header=TRUE,sep=";",dec=".")
forest<-rbind(forest01,forest02)
write.table(forest,file="./predvalues_forest.csv",sep=";",col.names=TRUE,row.names=FALSE)



### FORESTPLOTS OF DIFFERENCES IN PREDICTED DIFFERENCES ###
###########################################################

### FORMAT FILES ###

forest<-read.csv2("./predvalues_forest.csv",header=TRUE,sep=";",dec=".")
forest$det<-sub("_.*", "", forest$X)
forest$group<-sub(".*_", "", forest$X)

temp_lo<-forest[forest$group=="lo",]
temp_hi<-forest[forest$group=="hi",]
temp_score<-forest[forest$group=="score",]

vars01<-c("glucosa","hb1ac","ldl","hdl","trigli","pas","pad",
          "imc","cintura","egfr")
vars02<-c("glucose, mg/dL\n(difference in mean predicted values, 95% CI)",
          "Hb1Ac, %\n(difference in mean predicted values, 95% CI)",
          "LDL cholesterol, mg/dL\n(difference in mean predicted values, 95% CI)",
          "HDL cholesterol, mg/dL\n(difference in mean predicted values, 95% CI)",
          "triglycerides, mg/dL\n(difference in mean predicted values, 95% CI)",
          "systolic blood pressure, mmHg\n(difference in mean predicted values, 95% CI)",
          "diastolic blood pressure, mmHg\n(difference in mean predicted values, 95% CI)",
          "body mass index, kg/m2\n(difference in mean predicted values, 95% CI)",
          "waist circumference, cm\n(difference in mean predicted values, 95% CI)",
          "glomerular filtration rate, mL/min/1.73 m2\n(difference in mean predicted values, 95% CI)")
vars03<-c(1,2,4,5,6,7,8,9,10,11)

tab<-NULL
tab2<-NULL
tab3<-NULL
for(i in 1:length(vars01))
  
{
  xxx<-temp_lo[temp_lo$det==vars01[i],]
  names(xxx)<-NULL
  xxx<-rbind(cbind(as.matrix(xxx[2:4])),cbind(as.matrix(xxx[5:7])),cbind(as.matrix(xxx[8:10])),
             cbind(as.matrix(xxx[11:13])),cbind(as.matrix(xxx[14:16])),cbind(as.matrix(xxx[17:19])),
             cbind(as.matrix(xxx[20:22])))
  xxx<-as.data.frame(xxx)
  names(xxx)<-c("beta","lo","hi")
  xxx$det<-vars01[i]
  xxx$age<-as.numeric(c(0,6,12,18,24,30,36))
  xxx$coef<-ic_guapa2(guapa(xxx$beta),guapa(xxx$lo),guapa(xxx$hi))
  xxx$out<-paste("Low energy intake in breakfast (<20%) vs. reference (20-30%):\n",
                 vars02[i],sep="")
  xxx$out2<-vars03[i]
  tab<-rbind(tab,cbind(xxx))
  
  xxx<-temp_hi[temp_hi$det==vars01[i],]
  names(xxx)<-NULL
  xxx<-rbind(cbind(as.matrix(xxx[2:4])),cbind(as.matrix(xxx[5:7])),cbind(as.matrix(xxx[8:10])),
             cbind(as.matrix(xxx[11:13])),cbind(as.matrix(xxx[14:16])),cbind(as.matrix(xxx[17:19])),
             cbind(as.matrix(xxx[20:22])))
  xxx<-as.data.frame(xxx)
  names(xxx)<-c("beta","lo","hi")
  xxx$det<-vars01[i]
  xxx$age<-as.numeric(c(0,6,12,18,24,30,36))
  xxx$coef<-ic_guapa2(guapa(xxx$beta),guapa(xxx$lo),guapa(xxx$hi))
  xxx$out<-paste("High energy intake in breakfast (>30%) vs. reference (20-30%):\n",
                 vars02[i],sep="")
  xxx$out2<-vars03[i]
  tab2<-rbind(tab2,cbind(xxx))
  
  xxx<-temp_score[temp_score$det==vars01[i],]
  names(xxx)<-NULL
  xxx<-rbind(cbind(as.matrix(xxx[2:4])),cbind(as.matrix(xxx[5:7])),cbind(as.matrix(xxx[8:10])),
             cbind(as.matrix(xxx[11:13])),cbind(as.matrix(xxx[14:16])),cbind(as.matrix(xxx[17:19])),
             cbind(as.matrix(xxx[20:22])))
  xxx<-as.data.frame(xxx)
  names(xxx)<-c("beta","lo","hi")
  xxx$det<-vars01[i]
  xxx$age<-as.numeric(c(0,6,12,18,24,30,36))
  xxx$coef<-ic_guapa2(guapa(xxx$beta),guapa(xxx$lo),guapa(xxx$hi))
  xxx$out<-paste("Low breakfast quality (< median) vs. reference (> median):\n",
                 vars02[i],sep="")
  xxx$out2<-vars03[i]
  tab3<-rbind(tab3,cbind(xxx))
}

write.table(tab,file="./Forest_repeat/forest_lo.csv",sep=";",col.names=NA)
write.table(tab2,file="./Forest_repeat/forest_hi.csv",sep=";",col.names=NA)
write.table(tab3,file="./Forest_repeat/forest_score.csv",sep=";",col.names=NA)


### FOREST PLOTS ###

vars01<-c(1,2,4,5,6,7,8,9,10,11)
vars02<-c("glucosa","hb1ac","ldl",
          "hdl","trigli","pas","pad",
          "imc","cintura","egfr")
vars03<-c(2.5, 3.25, 3,
          -5.5, 2, 4, 3.5,
          2, 2, 3.5)

for(i in 1:length(vars01))
  
{
  forest_lo<-read.csv2("./Forest_repeat/forest_lo.csv",header=TRUE,sep=";",dec=".")
  dat<-forest_lo[forest_lo$out2==vars01[i],]
  label<-dat$out[1]
  
  figure<-ggplot() +
    geom_point(data=dat, aes(x=age, y=beta), size = 5, shape = 15, color="#b023fc") +
    geom_errorbar(data=dat, aes(x=age, y=beta, ymin=lo, ymax=hi), width=1.25, cex=1, color="#b023fc") +
    geom_text(data=dat, size=6.5, aes(y=max(hi)*vars03[i], x=age, label=coef, hjust='inward')) +
    coord_flip() +
    geom_hline(yintercept=0, linetype=2) +
    scale_x_continuous(breaks = dat$age) +
    xlab("Follow-up (months)") +
    ylab(label) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=16),
          axis.ticks.y=element_line(),
          axis.ticks.x=element_line(),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm"))
  
  namefile<-paste("./Forest_repeat/",vars02[i],"_f_lo.jpg",sep="")
  ggsave(filename=namefile, units="px", width=8000, height=8000, dpi=1200, bg="transparent")
  par(las=1,cex=1.2,bty="n",lheight=0.9)
  figure
  dev.off()
}


vars02<-c("glucosa","hb1ac","ldl",
          "hdl","trigli","pas","pad",
          "imc","cintura","egfr")
vars04<-c(2.5, 3.5, 3,
          -2.25, 2, 4, 3.25,
          2, 2, 3.75)

for(i in 1:length(vars01))
  
{
  forest_hi<-read.csv2("./Forest_repeat/forest_hi.csv",header=TRUE,sep=";",dec=".")
  dat<-forest_hi[forest_hi$out2==vars01[i],]
  label<-dat$out[1]
  
  figure<-ggplot() +
    geom_point(data=dat, aes(x=age, y=beta), size = 5, shape = 15, color="#3093fc") +
    geom_errorbar(data=dat, aes(x=age, y=beta, ymin=lo, ymax=hi), width=1.25, cex=1, color="#3093fc") +
    geom_text(data=dat, size=6.5, aes(y=max(hi)*vars04[i], x=age, label=coef, hjust='inward')) +
    coord_flip() +
    geom_hline(yintercept=0, linetype=2) +
    scale_x_continuous(breaks = dat$age) +
    xlab("Follow-up (months)") +
    ylab(label) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=16),
          axis.ticks.y=element_line(),
          axis.ticks.x=element_line(),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm"))
  
  namefile<-paste("./Forest_repeat/",vars02[i],"_f_hi.jpg",sep="")
  ggsave(filename=namefile, units="px", width=8000, height=8000, dpi=1200, bg="transparent")
  par(las=1,cex=1.2,bty="n",lheight=0.9)
  figure
  dev.off()
}


vars02<-c("glucosa","hb1ac","ldl",
          "hdl","trigli","pas","pad",
          "imc","cintura","egfr")
vars05<-c(2.25, 3.25, 2.75,
          5.5, 2.1, 2.25, 2.25,
          4.25, 2, -12)

for(i in 1:length(vars01))
  
{
  forest_score<-read.csv2("./Forest_repeat/forest_score.csv",header=TRUE,sep=";",dec=".")
  dat<-forest_score[forest_score$out2==vars01[i],]
  label<-dat$out[1]
  
  figure<-ggplot() +
    geom_point(data=dat, aes(x=age, y=beta), size = 5, shape = 15, color="#d1024f") +
    geom_errorbar(data=dat, aes(x=age, y=beta, ymin=lo, ymax=hi), width=1.25, cex=1, color="#d1024f") +
    geom_text(data=dat, size=6.5, aes(y=max(hi)*vars05[i], x=age, label=coef, hjust='inward')) +
    coord_flip() +
    geom_hline(yintercept=0, linetype=2) +
    scale_x_continuous(breaks = dat$age) +
    xlab("Follow-up (months)") +
    ylab(label) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=16),
          axis.ticks.y=element_line(),
          axis.ticks.x=element_line(),
          axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=14),
          plot.margin = unit(c(0.5,0.5,0.25,0.25), "cm"))
  
  namefile<-paste("./Forest_repeat/",vars02[i],"_f_score.jpg",sep="")
  ggsave(filename=namefile, units="px", width=8000, height=8000, dpi=1200, bg="transparent")
  par(las=1,cex=1.2,bty="n",lheight=0.9)
  figure
  dev.off()
}


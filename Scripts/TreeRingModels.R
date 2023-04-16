## ---------------------------
##
## Script name: E-W tree-ring analysis
##
## Author: Dr. Joan Dudney
##
## Date Created: 2022-03-09
##
## Copyright (c) Joan Dudney, 2022
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: this code runs the tree-ring 
## models and analyses checks; also includes
## supplemental tables 1-4
## ---------------------------


library(sjPlot)
library(ggeffects)
library(patchwork)
library(tidyverse)
library(lme4)
library(plotrix)
library(ggpubr)
library(mgcv)
library(nlme)
library(marginaleffects)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                    DATA CLEANING
#                The first part of this code reads in data and creates
#                new weather datasets to test different weather windows
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# READING IN AND CLEANING DATA  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## cleaning tree-ring cores
cores <- read_csv("Data/wbp_rwi_v4.csv") %>% 
  rename("X1" = "...1")

plotnames <- read_csv("Data/plot_names_normalized_v2.csv") 
treelook <- read_csv("Data/wbp_tree_lookup.csv")[-1]

## taking out duplicated cores
dup.cores <-c("FS12120A","FS1510A","FS1530A","FSD2B04B","FSD2B16B","FSD4B07B",
              "FSD4B16A","FSD5C18A","SE04B12B","SE04B28A","SE18A21A","SE24A11A",
              "SE24A27A","SE4309A","SE4316A","SE4319B","SE4329A","SE67A28B","SE67C12B",
              "YO03A06B","YO03A11B","YO03A13A","YO03A16A","YO03A24B","YO03A30B","YO25A29B")

mm <- match(dup.cores, treelook$COREID_ALT)
dup.cores[!is.na(mm)]<- as.character(treelook$COREID_2[na.omit(mm)])


plots <- plotnames%>%
  select(plot_id_needle, plot_label, plot_label_2)%>%
  rename(plot=plot_label_2)

cores_labels <- treelook %>% 
  select(COREID_2, COREID, plot_label_2, plot_id_needle, Tree_Num, lat,long, DBH_cm, Ht_m, Canopy_Pos) %>%   
  filter(!COREID_2 %in% dup.cores) %>% 
  filter(COREID_2!="site_01_06") %>% ## removing wpbr tree
  rename_all(tolower)


## combining datasets
cores_cleaner <- cores %>% 
  pivot_longer(-X1, names_to = "COREID_2") %>% 
  rename_all(tolower)%>%
  rename(year=x1) %>% 
  left_join(cores_labels) %>% 
  filter(coreid_2!="site_01_06") 


## prism data
prism <- read_csv("Data/prism_plots_1900.csv")
length(unique(prism$plot))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CREATING DIFFERENT WEATHER WINDOWS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## annual Jan 1-Dec 31
prism_annual <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, year) %>% 
  summarize(vpdall=mean(vpdmax, na.rm=T), tmaxall = mean(tmax, na.rm=T),
            pptall = sum(ppt, na.rm=T))%>%
  ungroup()

## annual water year (Oct 1-Sept 30)
prism_annual_growing <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  summarize(vpd.ann=mean(vpdmax, na.rm=T), tmax.ann = mean(tmax, na.rm=T),
            ppt.ann = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup()

# non-summer water year (Oct 1-June 30)
prism_growing <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  summarize(vpd=mean(vpdmax, na.rm=T), tmax = mean(tmax, na.rm=T),
            ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup()

## summer months (April-Sept)
prism_seasonSummerAp_S <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(month%in%c(4:9)) %>% 
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  summarize(vpdSummer=mean(vpdmax, na.rm=T), tmaxSummer = mean(tmax, na.rm=T),
            pptSummer = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup()

## shorter summer months (June-Sept)
prism_seasonSummerJune_S <- prism %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(month%in%c(6:9)) %>% 
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  summarize(vpdJ_S=mean(vpdmax, na.rm=T), tmaxJ_S = mean(tmax, na.rm=T),
            pptJ_S = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup()

## combining datasets
treenum <- cores_cleaner %>% 
  distinct(plot_id_needle, tree_num) %>% 
  mutate(tree_id=seq(1:771))


# full dataframe
clim_core1 <- cores_cleaner %>%
  filter(year>1900) %>%
  filter(!is.na(value)) %>% 
  left_join(prism_annual) %>% 
  left_join(prism_growing) %>% 
  #left_join(drought_years) %>% 
  left_join(prism_annual_growing) %>% 
  left_join(prism_seasonSummerAp_S) %>% 
  left_join(prism_seasonSummerJune_S) %>% 
  left_join(treenum) %>% 
  mutate(thresh_thirds = ifelse(tmax<4.993222, "0-33",
                                ifelse(tmax>=4.993222&tmax<6.401111 , "34-66", "67-100"))) %>% 
  mutate(thresh_quants = ifelse(tmax<4.590556, "0-25",
                                ifelse(tmax>=4.590556&tmax<9 , "26-75", "76-100"))) %>% 
  mutate()

length(unique(clim_core1$plot_id_needle))
hist(clim_core1$value)

## looking for outliers
summary(clim_core1$value)
Q1 <- as.numeric(quantile(clim_core1$value,  probs = c(.25), na.rm=T))
Q3 <-  as.numeric(quantile(clim_core1$value,  probs = c(.75), na.rm=T))

IQR <-  Q3 - Q1
Q1 - (1.5*IQR) 
(3*IQR)-Q1 #0.2425516

Q3 + (1.5*IQR)
Q3 + (3*IQR) #2.202722


## removing outliers
clim_core <- clim_core1 %>% 
  filter(value<2.202722) %>% 
  filter(value>0.2425516)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                      GAM Analysis
##              This code runs a model selection using AIC
##              to select the most parsimonious weather variables
##              
##              Second, the code runs a model selection to determine whether
##              adding nonlinear climate terms (smooths) improves the models
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## creating gam dataframe
gam_dat <- clim_core%>% 
  mutate(plot=factor(plot_id_needle), treeid=as.factor(tree_id)) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Selecting the most parsimonious climate variables 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## MODEL SELECTION

## Four models with annual climate, water year, summer_june (i.e., June-September)
## and summer_april (April-September)

model_ann <- gam(value~pptall +tmaxall, data=gam_dat)

model_growing_ann <- gam(value~ppt.ann +tmax.ann, data=gam_dat)

model_wateryr <- gam(value~ppt+tmax, data=gam_dat)

model_june <- gam(value~pptJ_S+tmaxJ_S, data=gam_dat)

model_april <- gam(value~pptSummer+tmaxSummer, data=gam_dat)

## Results: water year is the most parsimonious model, followed by annual variables
## AIC values reported in S2
AIC(model_ann, model_growing_ann, model_wateryr, model_june, model_april)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Selecting the most parsimonious model
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Scaling data
gam_dat_scale <- clim_core%>% 
  mutate(plot=factor(plot_id_needle), treeid=as.factor(tree_id)) %>% 
  mutate_at(scale, .vars = vars(tmax, dbh_cm, ht_m, ppt))%>%
  as.data.frame(.)

# Test1: comparing with and without dbh and height
modgam <- gam(value ~  ppt+tmax,
              data = gam_dat_scale)

modgamsimp = gam(value ~ dbh_cm + ht_m + ppt+tmax,
                 data = gam_dat_scale)

summary(modgamsimp)

## Result1: adding height and dbh lowers AIC
AIC(modgam, modgamsimp)


## Test 2: adding in plot random effect
modgam_plot = gam(value ~ dbh_cm + ht_m + ppt+tmax + s(plot, bs = "re"),
                  data = gam_dat_scale)

## Result 2: adding a plot random effect lowers AIC
AIC(modgamsimp, modgam_plot)


## Test 3: adding in a tree random effect
modgam_plot_tree = gam(value ~ dbh_cm + ht_m + ppt+tmax + s(plot, bs = "re")+ s(tree_id, bs = "re"),
                       data = gam_dat_scale)

## Result 3: adding a tree level random effect does not improve the model
AIC(modgamsimp, modgam_plot, modgam_plot_tree)


## Test 4: adding in nonlinear temperature term
modgam_tmax = gam(value ~ dbh_cm + ht_m + ppt+ s(tmax) + s(plot, bs = "re"),
                  data = gam_dat_scale)

## Result 4: adding in a smoothed term for temperature lowers the AIC
AIC(modgamsimp, modgam_plot, modgam_tmax)


## Test 5: adding in precipitation nonlinear term
modgam_all = gam(value ~ dbh_cm + ht_m + s(ppt)+ s(tmax) + s(plot, bs = "re"),
                 data = gam_dat_scale)

summary(modgam_all)


# Result 6: adding in precipitation slightly lowers AIC
AIC(modgamsimp, modgam_plot, modgam_tmax, modgam_all)


## Summary: adding in smoothed terms lowers AIC by a lot
AIC(modgam_plot) - AIC(modgam_all)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Verifying that GLMMs show similar results to GAMMs   
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Does adding plot random effects and quadratic terms improve fit?
modelglmm <- lmer(value~ppt+tmax+(1|plot_id_needle), data=clim_core)

model_wateryear <- lmer(value~ppt+I(ppt^2)+tmax+I(tmax^2)+dbh_cm+ ht_m+(1|plot_id_needle), data=clim_core)

## Yes, quadratic terms improve fit
AIC(modelglmm, model_wateryear)


## Supplemental table S2
tab_model(model_wateryear, 
          show.aic = T, show.se=T, show.ci = F, digits = 4)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MODEL CHECKS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## checking residuals
qqnorm(resid(model_wateryear))
qqline(resid(model_wateryear)) 

## heteroscedasticity and nonlinearities?
plot(model_wateryear)

## Leverage and Cook's distance
lev<-hat(model.matrix(model_wateryear))

## Plot leverage against standardized residuals
plot(resid(model_wateryear)~lev,las=1,ylab="Standardised residuals",xlab="Leverage")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SPATIAL AUTOCORRELATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Dataset is too large, it is maxing out the memory, so splitting into smaller sub-samples;
## Given this is spatial not temporal autocorrelation, subtests should suffice

small <- clim_core %>%
  filter(year>2000) %>%
  na.omit()

sp_model<-gls(value~ppt+I(ppt^2)+tmax+I(tmax^2), data=small)
plot(Variogram(sp_model))

small2 <- clim_core %>%
  filter(year<2000&year>1980)%>%
  na.omit()

sp_model4<-gls(value~ppt+I(ppt^2)+tmax+I(tmax^2), data=small2)
plot(Variogram(sp_model4))

small3 <- clim_core %>%
  filter(year<1960&year>1930) %>%
  na.omit()

sp_model7<-gls(value~ppt+I(ppt^2)+tmax+I(tmax^2), data=small3)
plot(Variogram(sp_model7))

## all variograms are not showing strong spatial autocorrelation at any distance;
## even without plot random effects

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TEMPORAL AUTOCORRELATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## not a strong year effect = dentrending was effective
ggplot(clim_core, aes(x=year, y=value))+
  geom_point()+
  geom_smooth(method="lm")

summary(lm(value~year, data=clim_core))

## significant because n is very high
acf((residuals(modgam_all_orig)))

# ## checking to see if corAR1 improves things
cleandat <- clim_core%>%
  select(value, ppt, tmax, year, tree_num, plot_id_needle) %>%
  mutate(plot_id_needle = factor(plot_id_needle), year=factor(year))

lme_fit1 <- lme(value~ppt+I(ppt^2)+tmax+I(tmax^2),
                random=~1|plot_id_needle,
                data=na.omit(cleandat))

lme_fit2 <- lme(value~ppt+I(ppt^2)+tmax+I(tmax^2),
                random=~1|plot_id_needle,
                data=na.omit(cleandat),
                correlation = corAR1())

tab_model(lme_fit1, lme_fit2,show.aic = T, show.se=T, show.ci = F, digits = 3)

## does not improve the model
## acf(resid(lme_fit2, type = "normalized"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Using ggpredict to estimate the marginal effects of temp and precip
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modgam_all_orig = mgcv::gam(value ~ dbh_cm + ht_m +s(ppt)+ s(tmax) + s(plot, bs = "re"),
                            data = gam_dat)

summary(modgam_all_orig)

predgam=ggpredict(modgam_all_orig , terms=c("ppt[all]"))

ppt_gam = plot(predgam) + ggtitle("GAMM Precipitation")+
  #geom_jitter(width = .02, height = .02, alpha=.5, color="grey")+
  ylab("Predicted (RWI)")+
  scale_x_continuous(name="Precipitation (mm)", breaks = scales::pretty_breaks(n = 10))

predgam_temp=ggpredict(modgam_all_orig, terms=c("tmax[all]"))

temp_gam=plot(predgam_temp) + ggtitle("GAMM Temperature")+
  #geom_jitter(width = .01, height = .01, alpha=.5, color="grey")+
  ylab("Predicted (RWI)")+
  scale_x_continuous(name="Max temperature (°C)", breaks = scales::pretty_breaks(n = 10))
temp_gam

temp_gam + ppt_gam

m1 <- modgam_all

# default ACF function:
acf(resid(m1), main="acf(resid(m1))")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                    MONTE CARLO SIMULATION
#            This code runs the MC sim for estimating the E-W threshold
#            
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(mvtnorm)
library(clubSandwich)
library(locfit)

##randomly sample from coefficient estimates and variance-covariance matrix

V_CR1 = vcov(modgam_all_orig)
V_CR1=as.matrix(V_CR1)

coef_vector = modgam_all_orig$coefficients

draw = rmvnorm(n = 1000, mean = coef_vector, sigma = V_CR1)

orig_dat <-gam_dat

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ESTIMATING THE PRECIPITATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df_sim <- data.frame()

for (i in 1:1000){

  ## run the monte carlo simulation
  d <-  draw[i,]
  gam_mod <-  modgam_all_orig
  gam_mod$coefficients <-  d
  orig_dat$vals_mcs <-  predict(gam_mod, newdata = orig_dat)

  temp_plot <- ggplot(na.omit(orig_dat), aes(x=ppt, y=vals_mcs))+
    geom_smooth(method='locfit', method.args = list(deg=1))

  ## extracting threshold point
  gb <- ggplot_build(temp_plot)

  exact_x_maximum <- gb$data[[1]]$x[which(diff(sign(diff(gb$data[[1]]$y)))==-2)+1]

  max <- max(exact_x_maximum)

  thresh=data.frame(thresh=max)

  df_sim=rbind(df_sim, thresh)
}

df_sim <- df_sim %>% 
  filter(thresh!="-Inf")
mean(df_sim$thresh)

quantile(df_sim$thresh, probs = c(.05, .95))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ESTIMATING THE TEMPERATURE THRESHOLD
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df_sim <- data.frame()

V_CR1 = vcov(modgam_all_orig)
V_CR1=as.matrix(V_CR1)

coef_vector = modgam_all_orig$coefficients

draw = rmvnorm(n = 1000, mean = coef_vector, sigma = V_CR1)

orig_dat <-gam_dat

for (i in 1:1000){
  
  ## run the monte carlo simulation
  d <-  draw[i,]
  gam_mod <-  modgam_all_orig
  gam_mod$coefficients <-  d
  orig_dat$vals_mcs <-  predict(gam_mod, newdata = orig_dat)
  
  temp_plot <- ggplot(na.omit(orig_dat), aes(x=tmax, y=vals_mcs))+
    geom_smooth(method="loess", span=1, se=F)
  
  ## extracting threshold point
  gb <- ggplot_build(temp_plot)
  
  exact_x_maximum <- gb$data[[1]]$x[which(diff(sign(diff(gb$data[[1]]$y)))==-2)+1]
  
  thresh=data.frame(thresh=exact_x_maximum)
  
  df_sim=rbind(df_sim, thresh)
}

#write_csv(df_sim, "Data/MCsimresults_theshold_temp_Jan16.csv")
df_sim <- df_sim %>% 
  filter(thresh!="-Inf")
mean(df_sim$thresh)
quantile(df_sim$thresh, probs = c(.05, .95))

## reading in data
mc_sim_res <- read_csv("Data/MCsimresults_theshold_temp_Jan16.csv")
precip_res <- read_csv("Data/precipitationthreshold.csv")
quantile(precip_res$thresh, probs = c(.05, .95))
mean(precip_res$thresh)


mean(mc_sim_res$thresh)
mean(mc_sim_res$thresh)
quantile(mc_sim_res$thresh, probs = c(.05, .95))
sd(mc_sim_res$thresh)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                 SUPPLEMENTAL TABLE 3
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newcoredat <- clim_core

## lower quartile
quant0_25 <- newcoredat %>% 
  filter(thresh_quants == "0-25")

mod0_25 = lmer(value~tmax+ppt+(1|plot_id_needle), data=quant0_25)
mod0_25_1 = lmer(value~tmax+ppt+(1|plot_id_needle/tree_id), data=quant0_25)

tab_model(mod0_25,mod0_25_1)
AIC(mod0_25, mod0_25_1)

## mid quantiles
quant26_75 <- newcoredat %>% 
  filter(thresh_quants == "26-75")

mod26_75 = lmer(value~tmax+ppt+(1|plot_id_needle), data=quant26_75)
mod26_75_1 = lmer(value~tmax+ppt+(1|plot_id_needle/tree_id), data=quant26_75)
tab_model(mod26_75, mod26_75_1)
AIC(mod26_75, mod26_75_1)

## upper quartile
quant76_100 <- newcoredat %>% 
  filter(thresh_quants == "76-100")

mod76_100 = lmer(value~tmax+ppt+(1|plot_id_needle), data=quant76_100)
mod76_100_1 = lmer(value~tmax+ppt+(1|plot_id_needle/tree_id), data=quant76_100)
tab_model(mod76_100, mod76_100_1)
AIC(mod76_100, mod76_100_1)

## Table S3
tab_model(mod0_25, mod26_75, mod76_100, show.ci = F, show.est = T, show.se=T, digits=4)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                 SUPPLEMENTAL TABLE 4
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## lower tercile
terc0_33 <- newcoredat %>% 
  filter(thresh_thirds == "0-33")

mod0_33 = lmer(value~tmax+ppt+(1|plot_id_needle), data=terc0_33)
mod0_33_1 = lmer(value~tmax+ppt+(1|plot_id_needle/tree_id), data=terc0_33)

AIC(mod0_33, mod0_33_1)

## mid tercile
terc34_66 <- newcoredat %>% 
  filter(thresh_thirds == "34-66")

mod34_66 = lmer(value~tmax+ppt+(1|plot_id_needle), data=terc34_66)
mod34_66_1 = lmer(value~tmax+ppt+(1|plot_id_needle/tree_id), data=terc34_66)

AIC(mod34_66, mod34_66_1)

## upper tercile
terc67_100 <- newcoredat %>% 
  filter(thresh_thirds == "67-100")

mod67_100 = lmer(value~tmax+ppt+(1|plot_id_needle), data=terc67_100)
tab_model(mod67_100)

## Table S4
tab_model(mod0_33, mod34_66, mod67_100, show.ci = F, show.est = T, show.se=T, digits=4)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                           SUPPLEMENTAL TABLE S5
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alldroughts <- clim_core %>%
  na.omit() %>%
  mutate(extreme=ifelse(year%in%c(2012:2015),"2012-2015",
                        ifelse(year%in%c(1976:1977), "1976-1977",
                               ifelse(year%in%c(1959:1961),"1959-1961","Remaining")))) %>%
  mutate(drought1SD_wet=ifelse(year%in%c(2012:2015),"2012-2015",
                               ifelse(year%in%c(1976:1977), "1976-1977",
                                      ifelse(year%in%c(1959:1961),"1959-1961",
                                             ifelse(year%in%c(1901, 1904, 1906, 1907, 1909, 1911, 1938, 
                                                              1952, 1956, 1967, 1969, 1978, 1980, 1982,
                                                              1983, 1986, 1995, 1997, 1998, 2006, 2011, 2017), "wet", "Remaining"))))) %>%
  mutate(thresh_thirds = ifelse(tmax<4.993222, "0-33",
                                ifelse(tmax>=4.993222&tmax<6.846000 , "34-75", "76-100")))


drought_diffthresh <- alldroughts %>%
  mutate(threshtmax=ifelse(tmax<8.4, "Below", "Above")) %>% 
  mutate(drought = ifelse(extreme == "Remaining", "ave", "drought")) %>% 
  group_by(threshtmax, plot_id_needle, tree_num, drought) %>%
  summarize(meanrwi=mean(value, na.rm=T)) %>%
  ungroup()

droughtdiff_mod <- drought_diffthresh %>% 
  pivot_wider(names_from = drought, values_from = meanrwi) %>% 
  group_by(threshtmax,plot_id_needle, tree_num) %>% 
  summarise(meandrought = mean(drought, na.rm=T), meanave=mean(ave, na.rm=T)) %>% 
  mutate(diffrwi = meandrought-meanave)

means <- droughtdiff_mod %>% 
  group_by(threshtmax) %>% 
  summarise(m=mean(diffrwi, na.rm=T), se=std.error(diffrwi))

modrwi <- lmer(diffrwi ~ threshtmax+(1|plot_id_needle), data=na.omit(droughtdiff_mod))
modrwi1 <- lm(diffrwi ~ threshtmax, data=na.omit(droughtdiff_mod))
AIC(modrwi, modrwi1)

tab_model(modrwi, show.ci = F, show.se = T, digits=4)

pred1 <- ggpredict(modrwi)
plot(pred1)

droughtdiff_fig <- droughtdiff_mod %>% 
  select(threshtmax, diffrwi) %>% 
  mutate(treeid = seq(1:length(diffrwi))) %>% 
  pivot_wider(names_from = threshtmax, values_from = diffrwi)

plotmeans <- droughtdiff_fig %>% 
  group_by(plot_id_needle) %>% 
  summarise(Above=mean(Above), Below=mean(Below))

plotmeans %>% 
  ggplot(aes(x=x) ) +
  xlim(-.5,.5)+
  ylim(-9.5,9.5)+
  xlab("Difference in length (cm)")+
  ylab("Density")+
  # Top
  geom_density( aes(x = Above, y = ..density..), fill="#953b17", alpha=.8 ) +
  geom_label( aes(x=2.5, y=-.8, label="Below E-W Threshold"), size=4.5, color="#1A2445") +
  # Bottom
  geom_density( aes(x = Below, y = -..density..), fill= "#1A2445",alpha=.8) +
  geom_label( aes(x=2.5, y=.8, label="Above E-W Threshold"), size=4.5, color="#953b17") +
  geom_vline(xintercept = 0, linetype="dashed")+
  coord_flip()



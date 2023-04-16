library(sjPlot)
library(ggeffects)
library(patchwork)
library(tidyverse)
library(lme4)
library(ggpubr)
library(ggtext) 
library(lubridate)


theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select


##=================================================================================================================
##                    READING IN AND CLEANING DATA  
##==================================================================================================================

needles <- read_csv("Data/needles_newnames_Sept15.csv")
prism_growing <- read_csv("Data/prism_plots_1900.csv")
isotopes <- read_csv("Data/isotopes_newnames_Sept15.csv")
plotnames <- read_csv("Data/plot_names_normalized_v2.csv")
sampledate <- read_csv("Data/plot_sampledate.csv")
demdat <- read_csv("Data/plot_demdat.csv")
treedat <- read_csv("Data/wbp_tree_lookup.csv")[c(11,14, 20,21,23)]
names(treedat) = tolower(names(treedat))

## combining and cleaning data
cleandate <- sampledate %>% 
  mutate(newdate=mdy(date)) %>% 
  mutate(day=yday(newdate)) %>% 
  mutate(month=month(newdate)) %>% 
  select(-c(date, newdate))


needledat <- needles %>% 
  #filter(year%in%c(2011:2017)) %>%
  mutate(drought=ifelse(year%in%c(2012:2015),"Drought","No drought")) %>% 
  left_join(plotnames) %>% 
  left_join(demdat) %>% 
  left_join(treedat)

prism_new <- prism_growing %>%
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
#filter(year%in%c(2011:2017))

prism_seasonSummerJune_S <- prism_growing %>%
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
  summarize(vpds=mean(vpdmax, na.rm=T), tmaxs = mean(tmax, na.rm=T),
            ppts = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup()

prism_season_ann <- prism_growing %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  summarize(vpdann=mean(vpdmax, na.rm=T), tmaxann = mean(tmax, na.rm=T),
            pptann = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup()


prism_zero <- prism_growing %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(growing!=1900) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing,month) %>% 
  summarize(zero=ifelse(tmax<=0, 1, 0),
            one = ifelse(tmax<=1, 1, 0),
            two = ifelse(tmax<=5, 1, 0))%>%
  group_by(plot_id_needle, growing) %>% 
  summarise(num_zerobelow = sum(zero),
            num_onebelow=sum(one),
            num_fivebelow=sum(two)) %>% 
  rename(year=growing) %>% 
  ungroup()


needclim <- needledat %>% 
  left_join(prism_new) %>%
  left_join(prism_season_ann) %>% 
  left_join(prism_seasonSummerJune_S) %>% 
  left_join(prism_zero) %>% 
  select(!c(plot_label, plot_number,plot_id_tree.ring, plot_id_genetics,plot_id_climate,plot_label_2)) %>% 
  #rename(plotid=plot_id_needle) %>% 
  unite(col='tree_id', c('tree_num', 'plot_id_needle'), sep='-', remove=F)


##ISOTOPES
## changing treeid to tree_num
isotopes1 <- isotopes %>% 
  mutate(tree_num = as.numeric(str_replace_all(treeid, "T", "")))


## combining and cleaning data
isodat <- isotopes1 %>% 
  #filter(year%in%c(2011:2016)) %>%
  mutate(drought=ifelse(year%in%c(2012:2015),"Drought","No drought")) %>% 
  left_join(plotnames) %>% 
  left_join(demdat) %>% 
  left_join(treedat) %>% 
  left_join(prism_zero) %>% 
  mutate(plot=plot_id_needle) %>% 
  select(!c(plot_label, plot_number,plot_id_tree.ring, plot_id_genetics,plot_id_climate,lat, long, plot_label_2, 
            drought, elevation, aspect, slope, dbh_cm, ht_m, canopy_pos)) %>% 
  unite(col='tree_id', c('tree_num', 'plot_id_iso'), sep='-', remove=F) %>% 
  dplyr::select(-treeid) %>% 
  left_join(cleandate)


isoclim <- isodat %>% 
  left_join(prism_new) %>% ##making tree_id
  left_join(prism_seasonSummerJune_S) 


## TREE RINGS

rwi <-  read_csv("Data/RWI_sept4.csv") %>%
  filter(year>2010) %>%
  mutate(drought=ifelse(year%in%c(2012:2015),"Drought","No drought")) %>% 
  select(plot_id_needle, year, value, tree_num, dbh_cm, ht_m) %>% 
  unite(col='tree_id', c('tree_num', 'plot_id_needle'), sep='-', remove=F) 




##=================================================================================================================
##                      Compiling all data
##==================================================================================================================

newiso <- read_csv("Data/iso_clim_wue_k.csv") %>% 
  select(-k) %>% 
  select(plot_id_needle, wue)

kvalues <- read_csv("Data/k_values_notcorrect.csv")

length(unique(needclim$plotid))
unique(isodat$plot_id_needle)

alliso <- isodat %>% 
  left_join(needclim) %>% 
  left_join(kvalues) %>% 
  mutate(thresh85=ifelse(tmax<4.8, "q15",
                         ifelse(tmax>=4.8&tmax<8.6, "q50", "q85"))) %>% 
  mutate(thresh=ifelse(tmax<7.1, "Below", 
                       ifelse(tmax>=7.1&tmax<9.5, "Between", "Above"))) %>% 
  filter(year>2010&year<2018)


checkdatiso = alliso %>% 
  filter_all(any_vars(is.na(.)))

length(unique(alliso$plot_id_iso))
unique(alliso$plot_id_iso)

##mean clim values
meanclim <- alliso %>% 
  group_by(plot_id_iso) %>% 
  summarize(meantmax=mean(tmax, na.rm=T), meantmaxs=mean(tmaxs, na.rm=T),
            meanppt=mean(ppt, na.rm=T), meanppts=mean(ppts, na.rm=T)) %>% 
  mutate_at(vars(meantmax, meantmaxs, meanppt, meanppts), round, 1)


## combined data with tree rings

allcombdat <- alliso %>% 
  left_join(rwi) %>% 
  filter(year>2010&year<2018) %>% 
  left_join(meanclim)

checkdat <- allcombdat %>% 
  filter_all(any_vars(is.na(.)))

length(unique(allcombdat$plot_id_needle))



##=================================================================================================================
##                      
##    
##                            Creating models for isotopes, needles, and rwi between 2011-2017
##    
##==================================================================================================================


## dataframes



##=================================================================================================================
##              ## ISOTOPE MODELS   
##==================================================================================================================
## data
alliso_scale <-allcombdat %>% 
  mutate_at(vars(tmax, ppt,ppts, vpds, tmaxs, dbh_cm, ht_m, elevation, perC,  perN, perC, n), scale)


modiso=lmer(delta_c~ tmax+ppt+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso1=lmer(delta_c~ tmax+I(tmax^2)+ppt+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso2=lmer(delta_c~ tmax+I(tmax^2)+(1|plot_id_needle/tree_id), data=alliso_scale)

tab_model(modiso, modiso1,modiso2,  show.aic = T)

prediso <- ggpredict(modiso2)
plot(prediso)
plot(prediso$tmax)

## vpd or tmax?
modiso=lmer(delta_c~ tmax+vpds+(1|plot_id_needle/tree_id), data=alliso) ## tmax and vpds is best
modiso1=lmer(delta_c~ vpd+vpds+(1|plot_id_needle/tree_id), data=alliso)
modiso2=lmer(delta_c~ tmax+tmaxs+(1|plot_id_needle/tree_id), data=alliso)
modiso3 = lmer(delta_c~ num_zerobelow+num_onebelow + num_fivebelow + (1|plot_id_needle/tree_id), data=alliso)

tab_model(modiso3, show.ci=F, show.se=T)

tab_model(modiso, modiso1,modiso2,modiso3,  show.aic = T)


## what about precip
modiso=lmer(delta_c~ tmax+vpds+ppts+(1|plot_id_needle/tree_id), data=alliso) ## best fit
modiso1=lmer(delta_c~ tmax+vpds+ppts+ppt+(1|plot_id_needle/tree_id), data=alliso)
modiso2=lmer(delta_c~ tmax+vpds+ppt+(1|plot_id_needle/tree_id), data=alliso)


## what about nonlinears?
modiso=lmer(delta_c~ tmax+vpds+ppts+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso1=lmer(delta_c~ tmax+I(tmax^2)+vpds+ppts+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso2=lmer(delta_c~ tmax+vpds+ppts+I(ppts^2)+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso3=lmer(delta_c~ tmax+vpds+I(vpds^2)+ppts+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso4=lmer(delta_c~ tmax+I(tmax^2)+tmaxs+ppt+ppts+(1|plot_id_needle/tree_id), data=alliso_scale)

tab_model(modiso, modiso1,modiso2,modiso3, modiso4, show.aic = T)

modiso1=lmer(delta_c~ tmax+I(tmax^2)+tmaxs+ppt+ppts+(1|plot_id_needle/tree_id), data=alliso_scale)
tab_model(modiso1,modiso2, show.aic = T, show.se=T, show.ci = F)


prediso <- ggpredict(modiso4)
plot(prediso)
plot(prediso$tmax)
tab_model(modiso1, show.aic = T, show.se=T, show.ci = F)

## climate window 
modiso=lmer(delta_c~ ppt+tmax+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso1=lmer(delta_c~ ppts+tmaxs+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso1.5 = lmer(delta_c~ tmax + ppts+ppt+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso1.7 = lmer(delta_c~ tmax + ppts+tmaxs+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso2=lmer(delta_c~ ppt+I(ppt^2)+tmax+I(tmax^2)+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso3=lmer(delta_c~ ppts+I(ppts^2)+tmaxs+I(tmaxs^2)+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso4=lmer(delta_c~ ppts*tmaxs+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso5=lmer(delta_c~ ppts+I(ppts^2)+tmaxs+I(tmaxs^2)+(1|plot_id_needle/tree_id), data=alliso_scale)

tab_model(modiso, modiso1,modiso1.5, modiso1.7, modiso2, modiso3, modiso4, modiso5, show.ci=F, show.aic = T, show.se = T)


modiso1.5 = lmer(delta_c~ tmax + ppts+ppt+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso1.7 = lmer(delta_c~ tmax +ppt+ ppts+tmaxs+(1|plot_id_needle/tree_id), data=alliso_scale) ## best fit by 10 AIC
modiso1.8 = lmer(delta_c~ tmax + ppts+I(ppts^2)+tmaxs+(1|plot_id_needle/tree_id), data=alliso_scale)
modiso1.9 = lmer(delta_c~ tmax * tmaxs+ppts+tmaxs+(1|plot_id_needle/tree_id), data=alliso_scale)

cor.test(alliso_scale$tmax, alliso_scale$tmaxs)
cor.test(alliso_scale$ppt, alliso_scale$ppts)


tab_model(modiso1.5, modiso1.7, modiso1.9,  show.ci=F, show.aic = T, show.se = T)
tab_model(modiso1.7, show.ci=F, show.aic = T, show.se = T)


prediso <- ggpredict(modiso5, terms=c( "tmaxs","ppts"))
prediso <- ggpredict(modiso1.7)
plot(prediso)

alliso %>% 
  ggplot(aes(tmaxs, delta_c, fill=drought))+
  geom_point()+
  geom_smooth(method="lm")

alliso %>% 
  ggplot(aes(tmax, tmaxs,color=plot_id_iso,fill=drought))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(~drought)


allcombdat %>% 
  ggplot(aes(mean_length, delta_c,color=plot_id_iso,fill=drought))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(~meanppt, scales="free")

## Thresholds

above_iso <- alliso_scale %>% 
  filter(thresh=="Above")

below_iso <- alliso_scale %>% 
  filter(thresh=="Below")




droughtyears=alliso_scale %>% 
  filter(drought=="Drought")

nondrought=alliso_scale %>% 
  filter(drought!="Drought")

aboveiso = lmer(delta_c~ tmaxs + tmax+ppt+ppts+(1|plot_id_needle/tree_id), data=above_iso)
aboveiso1 = lmer(delta_c~ tmaxs + ppt+ppts+(1|plot_id_needle/tree_id), data=above_iso)
aboveiso2 = lmer(delta_c~  tmax+ppt+ppts+(1|plot_id_needle/tree_id), data=above_iso)
aboveiso3 = lmer(delta_c~  tmaxs+tmax+ppts+(1|plot_id_needle/tree_id), data=above_iso)
aboveiso4 = lmer(delta_c~ tmax+I(tmax^2)+tmaxs+ppt+ppts+(1|plot_id_needle/tree_id), data=above_iso)

tab_model(aboveiso, aboveiso1,  aboveiso2,aboveiso3, show.ci=F, show.aic = T, show.se = T)

belowiso = lmer(delta_c~ tmaxs + tmax+ppt+ppts+(1|plot_id_needle/tree_id), data=below_iso)

aboveiso = lmer(delta_c~ tmax+tmaxs+ppt+ppts+(1|plot_id_needle/tree_id), data=above_iso)
belowiso = lmer(delta_c~ tmax+tmaxs+ppt+ppts+(1|plot_id_needle/tree_id), data=below_iso)

aboveiso1 = lmer(delta_c~ tmax+tmaxs+ppt+(1|plot_id_needle/tree_id), data=above_iso)
belowiso1 = lmer(delta_c~ tmax+tmaxs+ppt+(1|plot_id_needle/tree_id), data=below_iso)


aboveiso = lmer(delta_c~tmax+I(tmax^2)+ppt+I(ppt^2)+(1|plot_id_needle/tree_id), data=above_iso)
belowiso = lmer(delta_c~ tmax+I(tmax^2)+ppt+I(ppt^2)+(1|plot_id_needle/tree_id), data=below_iso)

gam_dat_iso <- alliso%>% 
  mutate(plot_id_iso=factor(plot_id_iso), tree_id=as.factor(tree_id)) 

above_iso <- gam_dat_iso %>% 
  filter(thresh=="Above")

below_iso <- gam_dat_iso %>% 
  filter(thresh=="Below")

aboveiso = gam(delta_c~s(tmax)+s(ppt)+s(plot_id_iso, bs="re")+s(tree_id, bs="re"), data=above_iso)
belowiso = gam(delta_c~s(tmax)+s(ppt)+s(plot_id_iso, bs="re")+s(tree_id, bs="re"),  data=below_iso)

prediso1 <- ggpredict(aboveiso)
plot(prediso1$tmax)

prediso2 <- ggpredict(belowiso)
plot(prediso2$tmax)


tab_model(aboveiso,belowiso, show.ci=F, show.aic = T, show.se = T)

above_iso <- alliso_scale %>% 
  filter(thresh=="Above")

below_iso <- alliso_scale %>% 
  filter(thresh=="Below")

alliso %>% 
  ggplot(aes(x=thresh, y=delta_c, fill=drought))+
  geom_boxplot()

mod <- lm(mean_length ~ drought+tmax+ppt+delta_c+perN, data = alliso %>% filter(thresh == "Above"))
summary(mod)

mod <- lm(mean_length ~ thresh * ppt, data = alliso %>% filter(thresh != "Between"))
summary(mod)


##=================================================================================================================
##                      
##    this model
##    
##    
##==================================================================================================================
##best model doesn't include drought

aboveiso = lmer(delta_c~ppt*tmax+ppts+tmaxs+perN+perC+(1|plot_id_needle/tree_id), data=above_iso) ##best model
aboveiso1 = lmer(delta_c~ppt*tmax+perN+perC+(1|plot_id_needle/tree_id), data=above_iso)
aboveiso2 = lmer(delta_c~ppt+tmax+perN+perC+(1|plot_id_needle/tree_id), data=above_iso)
belowiso = lmer(delta_c~ppt*tmax+ppts+tmaxs+perN+perC+(1|plot_id_needle/tree_id), data=below_iso)
belowiso1 = lmer(delta_c~ppt*tmax+perN+perC+(1|plot_id_needle/tree_id), data=below_iso)
belowiso2 = lmer(delta_c~ppt+tmax+ppts+tmaxs+perN+perC+(1|plot_id_needle/tree_id), data=below_iso)

tab_model(aboveiso, aboveiso1, aboveiso2,belowiso,belowiso1,belowiso2,show.ci=F, show.aic = T, show.se = T)
## these models include drought (which largely captures what happened during the winter, so removing those variables)

aboveiso = lmer(delta_c~ppts+tmaxs+drought+perN+perC+(1|plot_id_needle/tree_id), data=above_iso)
belowiso = lmer(delta_c~ppts+tmaxs+drought+perN+perC+(1|plot_id_needle/tree_id), data=below_iso)

tab_model(aboveiso,  belowiso,show.ci=F, show.aic = T, show.se = T)

prediso1 <- ggpredict(aboveiso)
prediso1 <- ggpredict(aboveiso,terms=c("tmax", "ppt"))
plot(prediso1)

prediso2 <- ggpredict(belowiso)
prediso2 <- ggpredict(belowiso,terms=c("tmax", "ppt"))
plot(prediso2)

aboveiso = lmer(delta_c~tmax+drought+(1|plot_id_needle/tree_id), data=above_iso)
belowiso = lmer(delta_c~tmax+drought+(1|plot_id_needle/tree_id), data=below_iso)

tab_model(aboveiso,  belowiso,show.ci=F, show.aic = T, show.se = T)

prediso1 <- ggpredict(aboveiso)
plot(prediso1$drought)

prediso2 <- ggpredict(belowiso)
plot(prediso2$drought)


above_iso %>% 
  ggplot(aes(x=ppts, y=delta_c, fill=plot_id_iso, color=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")

above_iso %>% 
  ggplot(aes(x=vpds, y=delta_c, fill=plot_id_iso, color=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")

above_iso %>% 
  ggplot(aes(x=ppt, y=delta_c, fill=plot_id_iso,color=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")

above_iso %>% 
  ggplot(aes(x=ppt, y=delta_c))+
  geom_point()+
  geom_smooth(method="lm")

below_iso %>% 
  ggplot(aes(x=ppt, y=delta_c))+
  geom_point()+
  geom_smooth(method="lm")

below_iso %>% 
  ggplot(aes(x=ppts, y=delta_c))+
  geom_point()+
  geom_smooth(method="lm")

below_iso %>% 
  ggplot(aes(x=ppt, y=delta_c, fill=plot_id_iso,color=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(~meanppt)






aboveiso = lmer(delta_c~ tmaxs*ppts+tmax*ppt+ppt+(1|plot_id_needle/tree_id), data=above_iso)
belowiso = lmer(delta_c~ tmaxs*ppts+tmax*ppt+ppt+(1|plot_id_needle/tree_id), data=below_iso)

nodroughtmod = lmer(delta_c~ tmax+tmaxs+ppt+ppts+perN+perC+mean_length+(1|plot_id_needle/tree_id), data=nondrought)
droughtmod = lmer(delta_c~ tmax+tmaxs+ppt+ppts+perN+perC+mean_length+(1|plot_id_needle/tree_id), data=droughtyears)

tab_model(nodroughtmod,droughtmod, show.ci=F, show.aic = T, show.se = T)

prediso <- ggpredict(nodroughtmod)
plot(prediso$ppt)
plot(prediso$mean_length)

prediso1 <- ggpredict(droughtmod)
plot(prediso1$ppt)
plot(prediso1$mean_length)


##=================================================================================================================
##                      
##    # below I think are the models that show the threshold the best
##    
##    
##==================================================================================================================

#

aboveiso = lmer(delta_c~ ppt+tmax+I(tmax^2)+perN+perC+(1|plot_id_needle/tree_id), data=above_iso)
belowiso = lmer(delta_c~ ppt+tmax+I(tmax^2)+perN+perC+(1|plot_id_needle/tree_id), data=below_iso)

aboveiso1 = lmer(delta_c~ ppt+tmax+perN+perC+drought+(1|plot_id_needle/tree_id), data=above_iso)
belowiso1 = lmer(delta_c~ ppt+tmax+perN+perC+drought+(1|plot_id_needle/tree_id), data=below_iso)


tab_model(aboveiso,aboveiso1, belowiso, belowiso1,show.ci=F, show.aic = T, show.se = T)


allcombdat %>% 
  ggplot(aes(x=factor(month), y=delta_c))+
  geom_boxplot()

allcombdat %>% 
  ggplot(aes(x=day, y=delta_c))+
  geom_point()+
  geom_smooth(span=1)




allcombdat %>% 
  ggplot(aes(x=tmax, y=delta_c, color=drought, fill=drought))+
  geom_point()+
  geom_smooth(span=1)+
  facet_grid(~meantmax, scales = "free")




allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=ppts, y=delta_c,color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=ppts, y=delta_c))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=ppt, y=delta_c))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=tmaxs, y=delta_c))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=tmaxs, y=delta_c,color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=tmax, y=delta_c,color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=ppts, y=delta_c,color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=ppt, y=delta_c,color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")


allcombdat %>% 
  ggplot(aes(x=ppt, y=delta_c,color=plot_id_iso, fill=drought))+
  geom_point()+
  geom_smooth(method="lm")

allcombdat %>% 
  ggplot(aes(x=ppt, y=delta_c))+
  geom_point()+
  geom_smooth()



allcombdat %>% 
  filter(thresh!="Between") %>% 
  ggplot(aes(x=tmax, y=delta_c))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~thresh, scales = "free")

allcombdat %>% 
  ggplot(aes(x=tmax, y=delta_c, color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(span=1)+
  geom_smooth(aes(x=tmax, y=delta_c, fill=tmax))

alliso %>% 
  ggplot(aes(x=tmaxs, y=delta_c, color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(~plot_id_iso, scales = "free")

alliso %>% 
  ggplot(aes(x=tmax, y=delta_c, color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(~plot_id_iso, scales = "free")

alliso %>% 
  ggplot(aes(x=tmaxs, y=vpds, color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(~month, scales = "free")


alliso %>% 
  ggplot(aes(x=ppts, y=vpds, color=plot_id_iso, fill=drought))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~.)


alliso %>% 
  ggplot(aes(x=tmaxs, y=delta_c, color=plot_id_iso, fill=plot_id_iso))+
  geom_point()+
  geom_smooth(method="lm")

prediso <- ggpredict(aboveiso)
plot(prediso$ppt)
plot(prediso$tmax)

prediso1 <- ggpredict(belowiso)
plot(prediso1$tmax)
plot(prediso1$ppt)

prediso <- ggpredict(aboveiso, terms=c("ppts[all]"))
plot(prediso)
prediso1 <- ggpredict(belowiso, terms=c("ppts[all]"))
plot(prediso1)

plot(prediso$tmax)
plot(prediso$ppt)

prediso <- ggpredict(belowiso, terms=c("ppt", "ppts"))
plot(prediso)
plot(prediso$ppts)
plot(prediso$ppt)


## adding complexity

modisocomp=lmer(delta_c~ ppts+I(ppts^2)+tmaxs+I(tmaxs^2)+factor(year)+(1|plot_id_needle/tree_id), data=alliso_scale) 
modisocomp1=lmer(delta_c~ ppts+I(ppts^2)+tmaxs+I(tmaxs^2)+perN+factor(year)+(1|plot_id_needle/tree_id), data=alliso_scale) ## best fit
modisocomp2=lmer(delta_c~ ppts+I(ppts^2)+tmaxs+I(tmaxs^2)+ perN+perC+factor(year)+(1|plot_id_needle/tree_id), data=alliso_scale)
modisocomp3=lmer(delta_c~ ppts+I(ppts^2)+tmaxs+I(tmaxs^2)+perN+perC+elevation+factor(year)+(1|plot_id_needle/tree_id), data=alliso_scale)

tab_model(modisocomp, modisocomp1, modisocomp2, modisocomp3, show.ci=F, show.aic = T, show.se = T)




##=================================================================================================================
##                      NEEDLE MODELS
##==================================================================================================================
## data
alliso_scale_need <- allcombdat %>% 
  mutate_at(vars(tmax, ppt, ppts, tmaxs, dbh_cm, ht_m, elevation, perC, delta_c, perN, perC, n), scale)

##determine which climate window to use -- nonlinear and winter temps

modneed <- lmer(mean_length~ drought+(1|plot_id_needle/tree_id), data=alliso_scale_need)
modneed0.5 <- lmer(mean_length~ ppt+tmax+(1|plot_id_needle/tree_id), data=alliso_scale_need)
modneed1 <- lmer(mean_length~ ppts+tmaxs+(1|plot_id_needle/tree_id), data=alliso_scale_need)
modneed2 <- lmer(mean_length~ ppt+I(ppt^2)+tmax+I(tmax^2)+(1|plot_id_needle/tree_id), data=allcombdat) ## best fit model #AIC 1832
modneed3 <- lmer(mean_length~ ppt+tmax+I(tmax^2)+(1|plot_id_needle/tree_id), data=allcombdat)

tab_model( modneed2, modneed3, show.ci=F, show.se=T, show.aic=T)

modneed <- lmer(mean_length~ ppt+tmax*ppt+(1|plot_id_needle/tree_id), data=allcombdat)
tab_model(modneed)

pred <- ggpredict(modneed, terms=c("tmax", "ppt"))
pred <- ggpredict(modneed2, terms="tmax[all]")
plot = plot(pred) 
plot

gb <- ggplot_build(plot)

exact_x_tmax <- gb$data[[1]]$x[which(diff(sign(diff(gb$data[[1]]$y)))==-2)+1]
exact_x_tmax

plot2=allcombdat %>% 
  ggplot(aes(x=tmax, y=mean_length, fill=drought))+
  geom_smooth()

# adding in complexity

modneedcomp <- lmer(mean_length~ ppt+I(ppt^2)+tmax+I(tmax^2)+perN+elevation+perC+delta_c+(1|plot_id_needle/tree_id), data=alliso_scale_need)
modneedcomp1=lmer(mean_length~ ppt+I(ppt^2)+tmax+I(tmax^2)+dbh_cm+ht_m +perN+elevation+perC+delta_c+(1|plot_id_needle/tree_id), data=alliso_scale_need)
modneedcomp2=lmer(mean_length~ ppt+I(ppt^2)+tmax+I(tmax^2)+dbh_cm+ht_m +perN+elevation+perC+delta_c*ppt+(1|plot_id_needle/tree_id), data=alliso_scale_need) ## best fit model
modneedcomp3=lmer(mean_length~ ppt+I(ppt^2)+tmax+I(tmax^2)+dbh_cm+ht_m +perN+elevation+perC+delta_c*tmax+(1|plot_id_needle/tree_id), data=alliso_scale_need)

tab_model(modneedcomp, modneedcomp1,modneedcomp2, modneedcomp3,show.ci=F, show.se=T, show.aic=T)


## marginal effects
predneed<- ggpredict(modneedcomp2)
plot(predneed)

predneed1 <- ggpredict(modneed2)
plot(predneed1$tmax)
plot(predneed1$drought)



##=================================================================================================================
##                      
##==================================================================================================================

modn=lmer(n~ ppt+tmax+(1|plot_id_iso/tree_id), data=allcombdat)
modn1=lmer(n~ ppt+I(ppt^2)+tmax+I(tmax^2)+(1|plot_id_iso/tree_id), data=allcombdat)
modn2=lmer(n~ ppt+tmax+I(tmax^2)+(1|plot_id_iso/tree_id), data=allcombdat)
tab_model(modn, modn2, modn1,show.ci = F, show.se = T, show.aic = T)


modn3=lmer(n~ drought*thresh+(1|plot_id_iso/tree_id), data=filter(allcombdat, thresh!="Between"))

tab_model(modn3)
predn3<- ggpredict(modn3, terms=c("drought", "thresh"))
plot(predn3)


threshn <- allcombdat %>% 
  filter(thresh!="Between")

modn3=lmer(n~ ppts+tmaxs+ppt+tmax+thresh+(1|plot_id_iso/tree_id), data=threshn)
modn3=lmer(n~ ppts*tmaxs+ppt+tmax+(1|plot_id_iso/tree_id), data=threshn)

tab_model(modn3,show.ci = F, show.se = T, show.aic = T)

predn3<- ggpredict(modn2, terms=c("tmax[all]"))
plot(predn3)



threshn %>% 
  ggplot(aes(x=factor(year), y=n))+
  geom_boxplot()

threshn %>% 
  ggplot(aes(x=factor(year), y=ppts))+
  geom_boxplot()

threshn %>% 
  ggplot(aes(x=factor(year), y=tmaxs, fill=thresh))+
  geom_boxplot()

threshn %>% 
  ggplot(aes(x=factor(year), y=ppts, fill=thresh))+
  geom_boxplot()

threshn %>% 
  ggplot(aes(x=thresh, y=n, fill=thresh))+
  geom_boxplot()+
  stat_compare_means()+
  facet_grid(~year)

threshn %>% 
  ggplot(aes(x=thresh, y=delta_c, fill=drought))+
  geom_boxplot()+
  stat_compare_means()

threshn %>% 
  ggplot(aes(x=ppts, y=n, fill=thresh))+
  geom_point()+
  geom_smooth()+
  facet_grid(~thresh)

library(mgcv)
gam_dat <- allcombdat%>% 
  mutate(plot=factor(plot_id_needle), treeid=as.factor(tree_id)) 

modgam_n = gam(n ~ dbh_cm + ht_m + s(ppt)+s(tmax) +thresh+ s(plot, bs = "re") + s(treeid, bs="re"),
                  data = gam_dat)

predn<- ggpredict(modgam_n)
plot(predn)

aboven = lmer(n~ ppts+tmaxs+(1|plot_id_needle/tree_id), data=above_iso)
belown = lmer(n~ ppts+tmaxs+(1|plot_id_needle/tree_id), data=below_iso)

tab_model(aboven, belown, show.ci=F, show.aic = T, show.se = T)

aboven1 = lmer(n~ ppts+tmaxs+ (1|plot_id_needle/tree_id), data=above_iso)
belown1 = lmer(n~ ppts+tmaxs+(1|plot_id_needle/tree_id), data=below_iso)

tab_model(aboven1, belown1, show.ci=F, show.aic = T, show.se = T)


predn1<- ggpredict(belown1)
plot(predn1)

predn2<- ggpredict(aboven1)
plot(predn2)


modpern=lmer(perN~ ppt+tmax+n+elevation+perC+year+(1|plot_id_iso/tree_id), data=allcombdat)
modpern0.5=lmer(perN~ ppts+tmaxs+n+elevation+perC+year+(1|plot_id_iso/tree_id), data=allcombdat)
modpern1=lmer(perN~ ppt+I(ppt^2)+tmax+I(tmax^2)+n+elevation+perC+year+(1|plot_id_iso/tree_id), data=allcombdat)


modneed=lmer(mean_length~ ppt+tmax+perN+n+elevation+perC+delta_c+(1|plot_id_iso/tree_num), data=alliso1)
modneed0.5=lmer(mean_length~ ppts+tmax+perN+n+elevation+perC+delta_c+(1|plot_id_iso/tree_num), data=alliso1)
modneed1=lmer(mean_length~ ppt+I(ppt^2)+tmax+I(tmax^2)+perN+n+elevation+perC+delta_c*ppt+(1|plot_id_iso/tree_num), data=alliso1) #best fit model
modneed1.5=lmer(mean_length~ ppt+I(ppt^2)+tmax+I(tmax^2)+perN+n+elevation+perC+delta_c*ppts+(1|plot_id_iso/tree_num), data=alliso1) 

tab_model(modneed,modneed0.5, modneed1,modneed1.5, show.ci = F, show.se = T, show.aic = T)


modrwi=lmer(value ~  ppt+tmax+perN+n+elevation+perC+delta_c+mean_length+(1|plot_id_iso/tree_id), data=all_dat)
modrwi0.3=lmer(value ~  ppt+tmax+perN+n+elevation+perC+delta_c*ppt+mean_length+(1|plot_id_iso/tree_id), data=all_dat)
modrwi0.5=lmer(value ~  ppts+tmax+perN+n+elevation+perC+delta_c*ppts+mean_length+(1|plot_id_iso/tree_id), data=all_dat)
modrwi0.7=lmer(value ~  ppt+tmax+perN+n+elevation+perC+delta_c*ppts+mean_length+(1|plot_id_iso/tree_id), data=all_dat)
modrwi1=lmer(value ~  ppt+I(ppt^2)+tmax+I(tmax^2)+perN+n+elevation+perC+delta_c*ppt+mean_length+(1|plot_id_iso/tree_id), data=all_dat)
modrwi1.1=lmer(value ~  ppt+I(ppt^2)+tmax+I(tmax^2)+perN+n+elevation+perC+delta_c*ppts+mean_length+(1|plot_id_iso/tree_id), data=all_dat)

tab_model(modrwi,modrwi0.3,modrwi1,modrwi1.1, show.ci = F, show.se = T, show.aic = T)



alliso %>%
  ggplot(aes(y=mean_length, x=tmaxs, fill=drought,color=drought))+
  geom_point()+
  geom_smooth()



alliso %>%
  ggplot(aes(x=tmaxs, y=delta_c, fill=drought,color=drought))+
  geom_point()+
  geom_smooth()


above_iso <- alliso %>% 
  filter(thresh=="Above")

belowiso <- alliso %>% 
  filter(thresh=="Below")

range(above_iso$n)

range(belowiso$n)

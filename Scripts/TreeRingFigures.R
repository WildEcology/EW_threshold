## ---------------------------
##
## Script name: Tree-ring Figures
##
## Author: Dr. Joan Dudney
##
## Date Created: 2022-03-09
##
## Copyright (c) Joan Dudney, 2022
## Email: jdudney@berkeley.edu
##
## ---------------------------
##
## Notes: this code creates figures 3 & 4
## as well as the majority of supplemental figures
##
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


theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  READING IN AND CLEANING DATA   
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## cleaning tree-ring cores
cores <- read_csv("Data/wbp_rwi_v4.csv")
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
  pivot_longer(-"...1", names_to = "COREID_2") %>% 
  rename_all(tolower)%>%
  rename(year="...1") %>% 
  left_join(cores_labels) %>% 
  filter(coreid_2!="site_01_06") 

## prism data
prism <- read_csv("Data/prism_plots_1900.csv")

# without-summer water year (Oct 1-June 30)
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

## combining datasets
treenum <- cores_cleaner %>% 
  distinct(plot_id_needle, tree_num) %>% 
  mutate(tree_id=seq(1:771))

# full dataframe
clim_core1 <- cores_cleaner %>%
  filter(year>1900) %>%
  filter(!is.na(value)) %>% 
  left_join(prism_growing) %>% 
  left_join(treenum) %>% 
  mutate(thresh_thirds = ifelse(tmax<4.993222, "0-33",
                                ifelse(tmax>=4.993222&tmax<6.443555 , "34-66", "67-100"))) %>% 
  mutate(thresh_quants = ifelse(tmax<4.590556, "0-25",
                                ifelse(tmax>=4.590556&tmax<9 , "26-75", "76-100"))) %>% 
  mutate()

## removing outliers
clim_core <- clim_core1 %>% 
  filter(value<2.202722) %>% 
  filter(value>0.2425516)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  GAMM MODELS AND MARGINAL EFFECTS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## creating gam dataframe
gam_dat <- clim_core%>% 
  mutate(plot=factor(plot_id_needle), treeid=as.factor(tree_id)) 


modgam_all_orig = mgcv::gam(value ~ dbh_cm + ht_m +s(ppt)+ s(tmax) + s(plot, bs = "re"),
                            data = gam_dat)

summary(modgam_all_orig)


## creating marginal effects figures
predgam=ggpredict(modgam_all_orig , terms=c("ppt[all]"))

ppt_gam = plot(predgam) + ggtitle("GAMM Precipitation")+
  ylab("Predicted (RWI)")+
  scale_x_continuous(name="Precipitation (mm)", breaks = scales::pretty_breaks(n = 10))

predgam_temp=ggpredict(modgam_all_orig, terms=c("tmax[all]"))

temp_gam=plot(predgam_temp) + ggtitle("GAMM Temperature")+
  ylab("Predicted (RWI)")+
  scale_x_continuous(name="Max temperature (°C)", breaks = scales::pretty_breaks(n = 10))
temp_gam

temp_gam + ppt_gam


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                    FIGURE 3
#                             Panels a, b, c, d
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## PANEL A
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## reading in MC sim data
mc_sim_res <- read_csv("Data/MCsimresults_theshold_temp_Jan16.csv")
mean(mc_sim_res$thresh)
quantile(mc_sim_res$thresh, c(0.05, .95))

## marginal effect of temperature
temp_gam=plot(predgam_temp) + ggtitle("GAMM Temperature")+
  ylab("Predicted (RWI)")+
  scale_x_continuous(name="Max temperature (°C)", breaks = scales::pretty_breaks(n = 10))
temp_gam

# final plot
predplot <- temp_gam+aes(fill="#225F75", color="#225F75")+
  geom_jitter(height=.005, alpha=.05)+
  annotate("text", y = 0.9, x =7.9, size=c(6),label = "E-W threshold", angle = 90)+
  theme_bw(base_size = 20)+
  scale_fill_manual(values=c("#225F75"))+
  scale_color_manual(values=c("#225F75"))+
  guides(fill="none", color="none")+
  geom_vline(xintercept = 8.4,  size=1.1)+
  geom_vline(xintercept = 7.12, linetype="dashed", size=.9)+
  geom_vline(xintercept = 9.51, linetype="dashed", size=.9)+
  #annotate("text", y = 170, x =6.75,  size=c(5),label = "E-W threshold", angle = 90)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("")+
  ylab("Predicted (RWI)")+
  scale_x_continuous(name="Max. temperature (°C)", breaks = scales::pretty_breaks(n = 10))

predplot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## PANEL B
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## defining the palette
more_colors <- c("#084f63", "#82bead", "#d6c594", "#c0831a","#953b17", "#642B30", "black")


## calculate drought based on mean ppt for each quantile

## estimating different precip quantiles
drought_quants_year <- prism_growing %>% 
  filter(year<2018) %>% 
  group_by(year) %>% 
  summarize(meanppt=mean(ppt), meantmax=mean(tmax)) %>% 
  summarize(quants = quantile(meanppt,  
                              probs = c(.05, .1, 0.15, .20, .25, .3, 
                                                  .35, .40, .50, .60, .8, .9)), 
                              probs = c(.05, .1, 0.15, .20, .25, .3, 
                                                  .35, .40, .50, .60, .8, .9)) 

## creating a data set with precip quantiles
drought_quants_yearcalc <- prism_growing %>% 
  filter(year<2018) %>% 
  group_by(year) %>% 
  summarize(meanppt=mean(ppt), meantmax=mean(tmax)) %>% 
  mutate(quant5 = ifelse(meanppt<516.8695, "drought", "non_drought"),
         quant10 = ifelse(meanppt<595.5522, "drought", "non_drought"),
         quant15 = ifelse(meanppt<622.2433, "drought", "non_drought"),
         quant20 = ifelse(meanppt<691.1462, "drought", "non_drought"),
         quant80 = ifelse(meanppt<1255.5237 ,"allother", "verywet"),
         quant90 = ifelse(meanppt<1415.5746 , "allother", "verywet"))


## combining with tree-ring data
newdat_quants <- clim_core %>% 
  filter(year<2018) %>% 
  left_join(drought_quants_yearcalc)

## estimating the change in growth (corresponding to the 20th quantile)
diffdat <- newdat_quants %>% 
  select(-meanppt, -meantmax) %>% 
  mutate(id=1:length(tree_num)) %>% 
  select(id, value, quant20, thresh_thirds, plot_id_needle, tree_num, year) %>% 
  pivot_wider(values_from = value, names_from = quant20) %>% 
  group_by(thresh_thirds, tree_num, plot_id_needle) %>%
  summarize(mdrought=mean(drought, na.rm=T), meannon=mean(non_drought, na.rm=T),
            diff=mdrought-meannon) %>% 
  group_by(thresh_thirds) %>% 
  summarize(meand=mean(diff, na.rm=T), std=std.error(diff, na.rm=T))


## final figure
diff_fig <- ggplot(diffdat, aes(x=thresh_thirds, y=meand, fill="blue", color="diff")) +
  geom_bar(stat="identity", width = .3) +
  geom_errorbar(aes(ymin = meand-std, ymax = meand+std), width = 0.2, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  scale_fill_manual(values=c("#82bead"))+
  scale_color_manual(values=c("#82bead"))+
  guides(fill="none", color="none")+
  xlab("Temperature terciles (°C)")+
  ylab("Δ RWI (dry yrs. - wet yrs.)")+
  #ylim(-.04,.04)+
  scale_x_discrete(labels =c("Low temps (0-33th)" = "0-33" , 
                             "Mid temps (34-66th)" = "34-66",
                             "High temps (56-100th)" = "67-100"))+
  theme(plot.title = element_text(hjust = 0.5))

diff_fig

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## PANEL C
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## final Panel C
ppt <- ggplot(newdat_quants, aes(x=ppt, y=value,color="color", fill="color"))+
  geom_smooth(method="lm")+
  facet_grid(.~thresh_thirds, scales = "free")+
  ylab("RWI")+
  xlab("Precipitation (mm)")+
  scale_fill_manual(values=c("#252F7E"))+
  scale_color_manual(values=c("#252F7E"))+
  guides(fill="none", color="none")+
  theme(strip.background = element_rect(fill="white"))+ 
  scale_x_continuous(breaks = c(1000, 2000))

ppt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## PANEL D
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax <- ggplot(newdat_quants, aes(x=tmax, y=value, color="color", fill="color"))+
  geom_smooth(method="lm")+
  facet_grid(.~thresh_thirds, scales = "free")+
  ylab("RWI")+
  xlab("Max. temperature (°C)")+
  scale_fill_manual(values=c("#401552"))+
  scale_color_manual(values=c("#401552"))+
  guides(fill="none", color="none")+
  theme(strip.background = element_rect(fill="white"))+
  scale_x_continuous(breaks = c(1,3.5,5.5,7,11.5,16))

tmax


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  FULL FIGURE 3
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


theme_set(
  theme_bw(base_size = 18)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)


## COMBINED FIGURE
predplot/diff_fig/ppt/tmax +  plot_layout(heights = c(2,2,1,1)) +
  plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                    FIGURE 4
#                     Data cleaning and creating panels a and b
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PANEL A
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## estimating the temperature quantiles
quantile(clim_core$tmax, c(.33,.75))

## classifying the most severe droughts
alldroughts<- clim_core %>%
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



## color palette
more_colors <- c("#084f63", "#82bead", "#d6c594", "#c0831a","#953b17", "#642B30", "black")

# Panel A
drought_quantiles_fig <- ggplot(alldroughts, aes(x=tmax, y=value, fill=extreme, color=extreme))+
  geom_smooth(method="lm",se=T, alpha=.5)+ facet_grid(.~thresh_thirds, scales = "free")+
  scale_fill_manual(values=more_colors, name="Extreme droughts")+
  scale_color_manual(values=more_colors,name="Extreme droughts")+
  theme(strip.background = element_rect(fill="white"))+
  ylab("RWI")+
  xlab("Max. temperature (°C)")+
  scale_x_continuous(breaks = c(1,3,5,6,7,10,13,16))+
  theme(legend.text=element_text(size=18),
        legend.title = element_text(size=20), legend.background = element_blank())

drought_quantiles_fig


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PANEL B
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## data
newdat = alldroughts %>% 
  mutate(drought=ifelse(extreme=="Remaining", "Nodrought", "Drought")) %>% 
  mutate(thresh8.4=ifelse(tmax>8.4, "Above", "Below"))

## estimating the change in growth during the extreme drought
drought_diff_dat <- alldroughts %>%
  group_by(thresh_thirds, plot_id_needle, tree_num, extreme) %>%
  summarize(meanrwi=mean(value, na.rm=T), meanppt=mean(ppt,na.rm=T), meantmax=mean(tmax,na.rm=T)) %>%
  ungroup()

## change in the mean growth across quantiles
diffdat_droughtmean <- drought_diff_dat %>%
  select(-meanppt, -meantmax) %>%
  pivot_wider(values_from = meanrwi, names_from = extreme) %>%
  select(-plot_id_needle, -tree_num) %>%
  group_by(thresh_thirds) %>%
  #summarize(m.val=mean(meanrwi, na.rm=T), std.val=std.error(meanrwi,na.rm=T)) %>% 
  summarise_all(mean, na.rm = TRUE) %>%
  rename("d76" = "1976-1977",
         "d12" = "2012-2015",
         "d59" = "1959-1961") %>%
  mutate(diff1976 = d76 -Remaining,
         diff2012 = d12-Remaining,
         diff1959 = d59-Remaining) %>%
  select(-c(d76,d12,d59, Remaining)) %>%
  pivot_longer(-thresh_thirds)

## estimating the standard error across quantiles
diffdat_droughtstd <- drought_diff_dat %>%
  select(-meanppt, -meantmax) %>%
  pivot_wider(values_from = meanrwi, names_from = extreme) %>%
  select(-plot_id_needle, -tree_num) %>%
  group_by(thresh_thirds) %>%
  #summarize(m.val=mean(meanrwi, na.rm=T), std.val=std.error(meanrwi,na.rm=T)) %>% 
  summarise_all(std.error, na.rm = TRUE) %>%
  rename("d76" = "1976-1977",
         "d12" = "2012-2015",
         "d59" = "1959-1961") %>%
  mutate(diff1976 = d76 -Remaining,
         diff2012 = d12-Remaining,
         diff1959 = d59-Remaining) %>%
  select(-c(d76,d12,d59, Remaining)) %>%
  pivot_longer(-thresh_thirds) %>% 
  rename(se = "value")

## combining data
exdrought <- diffdat_droughtmean %>% 
  left_join(diffdat_droughtstd)


## panel b — differences during extreme drought
dplot <- ggplot(exdrought, aes(x=thresh_thirds, y=value, fill=name)) +
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  geom_bar(position = "dodge", stat = "identity",  width = .4, alpha=.7) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se),color="darkgrey", width=.2,
                position=position_dodge(.4))+
  xlab("Temperature quantiles (°C)")+
  ylab("Δ RWI")+
  scale_fill_manual(values=more_colors, name="Extreme droughts",
                    labels =c("diff1959" = "1959-1961" ,
                              "diff1976" = "1976-1977",
                              "diff2012" = "2012-2015"))+
  scale_color_manual(values=more_colors,name="Extreme droughts",
                     labels =c("diff1959" = "1959-1961" ,
                               "diff1976" = "1976-1977",
                               "diff2012" = "2012-2015"))+
  scale_x_discrete(labels =c("Low temps (0-33th)" = "0-33" ,
                             "Mid temps (34-66th)" = "34-66",
                             "High temps (56-100th)" = "67-100"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.text=element_text(size=18),
        legend.title = element_text(size=20), legend.background = element_blank())


dplot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COMBINED FIGURE 4
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

fig4 <- drought_quantiles_fig / dplot

fig4 + plot_layout(widths = c(2,1.5)) +
  plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                    FIGURE S2
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## figures for supplement
precip_drought_fig <- ggplot(alldroughts, aes(x = reorder(extreme, +ppt), y=ppt, fill="grey"))+
  geom_boxplot(aes(fill = reorder(extreme, +ppt)),outlier.shape = NA, width=.5)+
  #ylim(0, 250)+
  scale_fill_grey() +
  ylab("Precipitation (mm)")+
  xlab("Extreme Droughts")+
  guides(fill="none")

precip_drought_fig

temp_drought_fig <- ggplot(alldroughts, aes(x = reorder(extreme, -tmax), y=tmax, fill=extreme))+
  geom_boxplot(aes(fill = reorder(extreme, -tmax)),outlier.shape = NA, width=.5)+
  ylim(0, 12)+
  scale_fill_grey() +
  ylab("Temperature (°C)")+
  xlab("Extreme Droughts")+
  guides(fill="none")

drought_comb_fig <- precip_drought_fig+temp_drought_fig 

drought_comb_fig + 
  plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                             SUPPLEMENTAL FIGURE S3
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PANEL A
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rwifig=ggplot(alldroughts, aes(x=tmax, y=value, color="drought", fill="drought"))+
  geom_point(alpha=.01)+
  geom_smooth()+
  ylab("RWI")+
  xlab("Max. temperature (°C)")+
  scale_fill_manual(values="#82bead")+
  scale_color_manual(values="#82bead")+
  guides(fill=F, color=F)

rwifig

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PANEL B
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rwifig_2= ggplot(alldroughts, aes(x=tmax, y=value, color="black"))+
  geom_smooth()+
  ylab("RWI")+
  xlab("Max. temperature (°C)")+
  scale_fill_manual(values="#82bead")+
  scale_color_manual(values="#82bead")+
  guides(fill=F, color=F)

rwifig_2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PANEL C
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Palette
col=c( "#c0831a", "#953b17", "#82bead" ,"darkblue")

raw_drought=ggplot(filter(alldroughts, drought1SD_wet!="Remaining"), aes(x=tmax, y=value, color=drought1SD_wet))+
  geom_point(alpha=.05)+
  geom_smooth(method="loess", span=1)+
  scale_fill_manual(values=col, name="Drought")+
  scale_color_manual(values=col,name="Drought")+
  ylab("RWI")+
  xlab("Max. temperature (°C)")+
  scale_x_continuous(breaks = c(1,3,5,7,9,11,13,15))

raw_drought

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COMBINED FIGURE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

(rwifig+rwifig_2 ) / raw_drought + plot_layout(heights = c(1.5,2)) +
  plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   SUPPLEMENTAL FIGURE S5
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meansrwi <- alldroughts %>% 
  group_by(plot_id_needle, year, extreme) %>% 
  summarize(meant=round(mean(tmax, na.rm=T), digits = 2), meanppt=round(mean(ppt, na.rm=T), digits=2),meanrwi=round(mean(value, na.rm=T), digits=2))


figs5 <- ggplot(meansrwi, aes(x=meanppt, y=meanrwi, color=extreme, fill=extreme))+
  geom_point()+
  geom_smooth()+
  #geom_smooth(method="lm")+
  geom_vline(xintercept = 1534.0, linetype="dashed")+
  ylab("Mean RWI (plot-level)")+
  guides(color = guide_legend(title = "Extreme drought"), fill = guide_legend(title = "Extreme drought"))+
  xlab("Mean precipitation (mm)")+
  scale_fill_manual(values=c("#084f63", "#82bead","#642B30", "grey"))+
  scale_color_manual(values=c("#084f63", "#82bead","#642B30", "grey"))

figs5



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                    FIGURE S6
#           
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## MC estimated threshold:1534.0 mm [1194.4; 2243.1 mm]
exact_x_ppt = 1534.0

predplot_ppt <- ppt_gam +aes(fill="#225F75", color="#225F75")+
  geom_jitter(height = .005, alpha=.05)+
  theme_bw(base_size = 15)+
  scale_fill_manual(values=c("#225F75"))+
  scale_color_manual(values=c("#225F75"))+
  guides(fill="none", color="none")+
  geom_vline(xintercept = exact_x_ppt, linetype="dashed", size=1.1)+
  geom_vline(xintercept = 1194.4, linetype="dashed", size=.5)+
  geom_vline(xintercept = 2243.1, linetype="dashed", size=.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("")+
  ylab("Predicted (RWI)")+
  scale_x_continuous(name="Precipitation (mm)", breaks = scales::pretty_breaks(n = 10))

predplot_ppt


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                         FIGURE S7:  WET YEAR CONDITIONS
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PANEL A
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# estimating differences in growth
diffdat90 <- newdat_quants %>% 
  select(-meanppt, -meantmax) %>% 
  mutate(id=1:length(tree_num)) %>% 
  select(id, value, quant90, thresh_thirds, plot_id_needle, tree_num, year) %>% 
  pivot_wider(values_from = value, names_from = quant90) %>% 
  group_by(thresh_thirds, tree_num, plot_id_needle) %>%
  summarize(maverage=mean(allother, na.rm=T), meanwet=mean(verywet, na.rm=T),
            diff=meanwet-maverage) %>% 
  group_by(thresh_thirds) %>% 
  summarize(meand=mean(diff, na.rm=T), std=std.error(diff, na.rm=T))

diff_fig90 <- ggplot(diffdat90, aes(x=thresh_thirds, y=meand, fill="blue", color="diff")) +
  geom_bar(stat="identity", width = .5) +
  geom_errorbar(aes(ymin = meand-std, ymax = meand+std), width = 0.2, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  scale_fill_manual(values=c("#82bead"))+
  scale_color_manual(values=c("#82bead"))+
  guides(fill="none", color="none")+
  xlab("Temperature terciles (°C)")+
  ylab("Δ RWI (very wet - remaining)")+
  ggtitle("90th percentile")+
  #ylim(-.04,.04)+
  scale_x_discrete(labels =c("Low temps (0-33th)" = "0-33" , 
                             "Mid temps (34-66th)" = "34-66",
                             "High temps (56-100th)" = "67-100"))+
  theme(plot.title = element_text(hjust = 0.5))

# estimating differences in growth
diffdat80 <- newdat_quants %>% 
  select(-meanppt, -meantmax) %>% 
  mutate(id=1:length(tree_num)) %>% 
  select(id, value, quant80, thresh_thirds, plot_id_needle, tree_num, year) %>% 
  pivot_wider(values_from = value, names_from = quant80) %>% 
  group_by(thresh_thirds, tree_num, plot_id_needle) %>%
  summarize(maverage=mean(allother, na.rm=T), meanwet=mean(verywet, na.rm=T),
            diff=meanwet-maverage) %>% 
  group_by(thresh_thirds) %>% 
  summarize(meand=mean(diff, na.rm=T), std=std.error(diff, na.rm=T))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PANEL B
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_fig80 <- ggplot(diffdat80, aes(x=thresh_thirds, y=meand, fill="blue", color="diff")) +
  geom_bar(stat="identity", width = .5) +
  geom_errorbar(aes(ymin = meand-std, ymax = meand+std), width = 0.2, color="black") +
  geom_hline(yintercept = 0, linetype="dashed", color="black")+
  scale_fill_manual(values=c("#82bead"))+
  scale_color_manual(values=c("#82bead"))+
  guides(fill="none", color="none")+
  xlab("Temperature terciles (°C)")+
  ylab("Δ RWI (wet - remaining)")+
  ggtitle("80th percentile")+
  #ylim(-.04,.04)+
  scale_x_discrete(labels =c("Low temps (0-33th)" = "0-33" , 
                             "Mid temps (34-66th)" = "34-66",
                             "High temps (56-100th)" = "67-100"))+
  theme(plot.title = element_text(hjust = 0.5))

diff_fig80

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COMBINED FIGURE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_fig80 + diff_fig90 + 
  plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   SUPPLEMENTAL FIGURE S8
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# estimating mean RWI across years
meansrwi <- alldroughts %>% 
  group_by(plot_id_needle, year, extreme) %>% 
  summarize(meant=round(mean(tmax, na.rm=T), digits = 2), meanppt=round(mean(ppt, na.rm=T), digits=2),meanrwi=round(mean(value, na.rm=T), digits=2))

# estimating mean RWI
mean(meansrwi$meanrwi)


figs8 <- ggplot(meansrwi, aes(x=year, y=meanrwi, color=extreme, fill=extreme))+
  geom_boxplot(aes(x=factor(year), y= meanrwi), outlier.color = "NA")+
  geom_hline(yintercept = 0.9885326, linetype="dashed")+
  ylab("Mean RWI (plot-level)")+
  guides(color = guide_legend(title = "Extreme drought"), fill = guide_legend(title = "Extreme drought"))+
  xlab("Year")+
  scale_fill_manual(values=c("#084f63", "#82bead","#642B30", "grey"))+
  scale_color_manual(values=c("#084f63", "#82bead","#642B30", "grey"))+
  scale_x_discrete(labels = abbreviate)+
  theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 1, size=8))

figs8




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   SUPPLEMENTAL FIGURE S10
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


theme_set(
  theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

figs10 <- clim_core %>% 
  ggplot(aes(x=tmax, y=value))+
  geom_smooth(method="loess", span=1)+
  facet_wrap(~plot_id_needle, ncol = 4, scales="free")+
  xlab("Max temperature (°C)")+
  ylab("RWI")

figs10

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   SUPPLEMENTAL FIGURE S11
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figs11 <- clim_core %>% 
  ggplot(aes(x=ppt, y=value))+
  geom_smooth(method="loess", span=1)+
  facet_wrap(~plot_id_needle, ncol = 4, scales="free")+
  xlab("Precipitation (mm)")+
  ylab("RWI")

figs11


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                   SUPPLEMENTAL FIGURE S9
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rwifig <- clim_core %>% 
  filter(year>2009) %>%
  mutate(drought=ifelse(year%in%c(2012:2015),"Drought","No Drought")) %>% 
  ggplot(aes(x=factor(year), y=value, fill=drought, color="black"))+
  geom_jitter(width = .2)+
  geom_boxplot()+
  labs(y="RWI",
       x="Year")+
  #guides(color=F, fill=F)+
  theme(legend.title = element_blank(),
        legend.position = c(.83,.15),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank())+
  scale_color_manual(values=c("grey"))+
  scale_fill_manual(values=c("#e9c46a", "#084f63"))+
  guides(color=F)

rwifig









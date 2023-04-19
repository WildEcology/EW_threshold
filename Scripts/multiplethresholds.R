
## ---------------------------
##
## Script name: Multiple co-occurring thresholds
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-11-19
##
## Copyright (c) Joan Dudney, 2021
## Email: dudney@ucsb.edu
##
## ---------------------------
##
##  Notes: The first section of code creates Figure 5 — multiple thresholds 
##  in stable carbon and nitrogen isotopes and needle length. It also estimates
##  whether the E-W threshold (8.4°C) co-occurs with thresholds in stable carbon 
##  and nitrogen isotopes and needle length.
##  
##  The second section of the code creates Table S5.
##
## ---------------------------


library(sjPlot)
library(ggeffects)
library(patchwork)
library(tidyverse)
library(lme4)
library(ggpubr)
library(ggtext)


## set a theme for the figures
theme_set(
  theme_bw(base_size = 17)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# READING IN AND CLEANING DATA  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## reading in data files
needles <- read_csv("Data/needles_newnames_Sept15.csv")
prism_growing <- read_csv("Data/prism_plots_1900.csv")
isotopes <- read_csv("Data/isotopes_newnames_Sept15.csv")
plotnames <- read_csv("Data/plot_names_normalized_v2.csv")
demdat <- read_csv("Data/plot_demdat.csv")
treedat <- read_csv("Data/wbp_tree_lookup.csv")[c(11,14, 20,21,23)]
names(treedat) = tolower(names(treedat))

## combining and cleaning data
needledat <- needles %>% 
  filter(year%in%c(2011:2017)) %>%
  mutate(drought=ifelse(year%in%c(2012:2015),"Drought","No Drought")) %>% 
  left_join(plotnames) %>% 
  left_join(demdat) %>% 
  left_join(treedat)

## estimating climate windows
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
  ungroup() %>% 
  filter(year%in%c(2011:2017))

meanclim <- prism_growing %>%
  mutate(plot=replace(plot, plot == "K_NA", 'K')) %>%
  mutate(plot=replace(plot, plot == "RC_NA", 'RC')) %>%
  mutate(plot=replace(plot, plot == "JM_NA", 'JM')) %>%
  rename(plot_id_needle=plot) %>% 
  mutate(month=as.numeric(month),
         growing=ifelse(month%in%c(10:12), year+1, year)) %>% 
  filter(!month%in%c(7:9)) %>% ## don't include July, August, September
  filter(growing>2010) %>% 
  pivot_wider(names_from=type, values_from=value) %>% 
  group_by(plot_id_needle, growing) %>% 
  summarize(vpd=mean(vpdmax, na.rm=T), tmax = mean(tmax, na.rm=T),
            ppt = sum(ppt, na.rm=T))%>%
  rename(year=growing) %>% 
  ungroup() %>% 
  group_by(plot_id_needle) %>%
  summarize(m.tmax=mean(tmax, na.rm=T), m.ppt=mean(ppt, na.rm=T)) 

## Combining cleaned data frames
needclim <- needledat %>% 
  left_join(prism_new) %>% 
  left_join(meanclim) %>% 
  rename(plotid=plot_id_needle)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##
##                THE FOLLOWING CODE INCLUDES ALL PANELS IN FIGURE 5
##                It starts with needle length, followed by stable carbon isotopes 
##                and stable nitrogen isotopes; the code first calculates the relative change
##                during drought which is then used to create the density plots
##
##
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NEEDLE LENGTH FIGURES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Estimating relative change in needle length during drought above and below the 
## E-W threshold

needles_diff <- needclim %>%
  select(-drought) %>% 
  mutate(drought=ifelse(year%in%c(2012:2015),"drought","nodrought")) %>%
  group_by(plotid,tree_num, drought, m.tmax) %>%
  summarize(m.length=mean(mean_length, na.rm=T)) %>%
  pivot_wider(names_from=drought, values_from=m.length) %>%
  mutate(diff=drought-nodrought) %>% 
  mutate(thresh=ifelse(m.tmax<8.4, "Below", "Above"))

needles_densfig <- needles_diff %>% 
  select(diff, thresh, tree_num, plotid) %>% 
  pivot_wider(names_from = thresh, values_from = diff)

## density figure
needfigdens <- needles_densfig %>% 
  ggplot(aes(x=x) ) +
  xlim(-2.8,2.9)+
  ylim(-1.5,1.5)+
  xlab("Difference in length (cm)")+
  ylab("Density")+
  # Upper
  geom_density(aes(x = Above, y = ..density..), fill="#953b17", alpha=.8 ) +
  geom_label(aes(x = 2.5, y = -0.8, label="Below E-W Threshold"), size=4.5, color="#1A2445") +
  
  # Lower
  geom_density(aes(x = Below, y = -..density..), fill= "#1A2445",alpha = 0.8) +
  geom_label(aes(x = 2.5, y = 0.8, label="Above E-W Threshold"), size = 4.5, color="#953b17") +
  geom_vline(xintercept = 0, linetype="dashed")+
  coord_flip()

needfigdens

## raw data figure
needle <- ggplot(needclim, aes(x=tmax, y=mean_length, fill=drought, color=drought))+
  geom_point(alpha=.2)+
  geom_smooth(method="loess", span=1)+
  ylab("Needle length (cm)")+
  xlab("Maximum temperature (°C)")+
  scale_color_manual(values=c( "#C99B55","#084f63"))+
  scale_fill_manual(values=c("#e9c46a", "#084f63"))+
  theme(legend.title = element_blank(),
        legend.position = c(.83,.15),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank())+
  geom_vline(xintercept = 8.4,  size=1)+
  geom_vline(xintercept = 7.12, linetype="dashed", size=.5)+
  geom_vline(xintercept = 9.51, linetype="dashed", size=.5)+
  annotate("text", y = 5.9, x =7.95,  size=c(5),label = "E-W threshold", angle = 90)

needle


needle + needfigdens


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STABLE CARBON ISOTOPES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## combining and cleaning data
isodat <- isotopes %>% 
  mutate(tree_num = as.numeric(str_replace_all(treeid, "T", ""))) %>% 
  filter(year%in%c(2011:2017)) %>%
  mutate(drought=ifelse(year%in%c(2012:2015),"Drought","No Drought")) %>% 
  left_join(plotnames) %>% 
  left_join(demdat) %>% 
  left_join(treedat)

isoclim <- isodat %>% 
  left_join(prism_new) %>% 
  left_join(meanclim)

## Estimating relative change in carbon isotopes during drought above and below the 
## E-W threshold
iso_diff <- isoclim %>%
  select(-drought) %>% 
  mutate(drought=ifelse(year%in%c(2012:2015),"drought","nodrought")) %>%
  group_by(plot_id_needle,treeid, drought, m.tmax) %>%
  summarize(m.iso=mean(delta_c, na.rm=T)) %>%
  pivot_wider(names_from=drought, values_from=m.iso) %>%
  mutate(diff=drought-nodrought) %>% 
  mutate(thresh=ifelse(m.tmax<8.4, "Below", "Above"))

ggplot(iso_diff, aes(x=m.tmax, y=diff))+
  geom_point()+
  geom_smooth()

iso_densfig <- iso_diff %>% 
  select(diff, thresh, plot_id_needle,treeid) %>% 
  pivot_wider(names_from = thresh, values_from = diff)

## density figure
isofigdens <- iso_densfig %>% 
  ggplot(aes(x=x))+
  xlim(-2,2)+
  ylim(-1.5,1.5)+
  labs(x="Difference in *&delta;*<sup>13</sup>C (&permil;)",
       y="Density")+
  # Top
  geom_density( aes(x = Above, y = ..density..), fill="#953b17", alpha=.8 ) +
  geom_label( aes(x=1.6, y=-.8, label="Below E-W Threshold"), size=4.5, color="#1A2445") +
  # Bottom
  geom_density( aes(x = Below, y = -..density..), fill= "#1A2445",alpha=.8) +
  geom_label( aes(x=1.6, y=.8, label="Above E-W Threshold"), size=4.5, color="#953b17") +
  geom_vline(xintercept = 0, linetype="dashed")+
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())+
  coord_flip()

isofigdens

## carbon isotope figure
cplotdrought = ggplot(isoclim, aes(x=tmax, y=c, fill=drought, color=drought))+
  geom_point(alpha=.2)+
  geom_smooth(method="loess", span=1)+
  labs(y="Needle *&delta;*<sup>13</sup>C (&permil;)",
       x="Maximum temperature (°C)")+
  #guides(color=F, fill=F)+
  theme(legend.title = element_blank(),
        legend.position = c(.83,.15),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank())+
  scale_color_manual(values=c("#C99B55","#084f63"))+
  scale_fill_manual(values=c("#e9c46a", "#084f63"))+
  geom_vline(xintercept = 8.4,  size=1)+
  geom_vline(xintercept = 7.12, linetype="dashed", size=.5)+
  geom_vline(xintercept = 9.51, linetype="dashed", size=.5)+
  annotate("text", y = -24.2, x =7.95,  size=c(5),label = "E-W threshold", angle = 90)

cplotdrought



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NITROGEN ISOTOPES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## Estimating relative change in nitrogen isotopes during drought above and below the 
## E-W threshold

n_diff <- isoclim %>%
  select(-drought) %>% 
  mutate(drought=ifelse(year%in%c(2012:2015),"drought","nodrought")) %>%
  group_by(plot_id_iso,plot_id_needle, tree_num, m.tmax, drought) %>% 
  summarize(meann=mean(n, na.omit=T)) %>% 
  pivot_wider(names_from = drought, values_from = meann) %>% 
  mutate(diff=drought-nodrought)%>% 
  mutate(thresh=ifelse(m.tmax<8.4, "Below", "Above"))


n_densdat <- n_diff %>% 
  select(diff, thresh, plot_id_needle,tree_num) %>% 
  pivot_wider(names_from = thresh, values_from = diff)

## density figure
nfigdens <- n_densdat %>% 
  ggplot(aes(x=x) ) +
  xlim(-1.2,1.2)+
  ylim(-3.5,3.5)+
  labs(x="Difference in *&delta;*<sup>15</sup>N (&permil;)",
       y="Density")+
  # Top
  geom_density( aes(x = Above, y = ..density..), fill="#953b17", alpha=.8 ) +
  geom_label( aes(x=1, y=-1.9, label="Below E-W Threshold"),size=4.5, color="#1A2445") +
  # Bottom
  geom_density( aes(x = Below, y = -..density..), fill= "#1A2445",alpha=.8) +
  geom_label( aes(x=1, y=1.9, label="Above E-W Threshold"), size=4.5, color="#953b17") +
  geom_vline(xintercept = 0, linetype="dashed")+
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_markdown())+
  coord_flip()

nfigdens

## nitrogen isotopes figure
n_tmax <- isoclim %>%
  ggplot(aes(x=tmax, y=n, fill=drought, color=drought))+
  geom_point(alpha=.2)+
  geom_smooth(method="loess", span=1)+
  scale_color_manual(values=c("#C99B55","#084f63"))+
  scale_fill_manual(values=c("#e9c46a", "#084f63"))+
  geom_vline(xintercept = 8.4,  size=1)+
  geom_vline(xintercept = 7.12, linetype="dashed", size=.5)+
  geom_vline(xintercept = 9.51, linetype="dashed", size=.5)+
  annotate("text", y = 3, x =7.95,  size=c(5),label = "E-W threshold", angle = 90)+
  #guides(fill=F, color=F)+
  labs(y="Needle *&delta;*<sup>15</sup>N (&permil;)",
       x="Maximum temperature (°C)")+
  theme(legend.title = element_blank(),
        legend.position = c(.83,.15),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank())

n_tmax 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   COMBINING ALL PANELS TO CREATE FIGURE 5
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allfigs <- (needle + needfigdens)/(cplotdrought+isofigdens)/(n_tmax+nfigdens)
allfigs + plot_layout(widths  = c(2,1)) + plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                            TABLE S5: ESTIMATING CO-OCCURRING THRESHOLDS
#     This code estimates whether the difference in needle length, stable carbon
#     and stable nitrogen isotopes are significantly correlated with the 8.4°C threshold estimated
#     using the tree-ring data
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# NEEDLE LENGTH MODELS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Using AIC to determine the best fit model
modneed <- lm(diff~thresh, data=needles_diff) 
modneed1 <- lmer(diff~thresh+ (1|plotid), data=needles_diff)## lower AIC; better model fit

AIC(modneed, modneed1)
summary(modneed1)

predneed=ggpredict(modneed)

# figure of model results
plot(predneed)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STABLE CARBON ISOTOPES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Using AIC to determine the best fit model
modc <- lm(diff~thresh, data=iso_diff)
modc1 <- lmer(diff~thresh + (1|plot_id_needle), data=iso_diff) ## lower AIC; better model fit
AIC(modc, modc1)
summary(modc1)


predc=ggpredict(modc1)
plot(predc)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# STABLE NITROGEN ISOTOPES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Using AIC to determine the best fit model
modn <- lm(diff~thresh, data=n_diff) ## lower AIC; better model fit
modn1 <- lmer(diff~thresh + (1|plot_id_needle), data=n_diff)
AIC(modn, modn1)
summary(modn1)

predn=ggpredict(modn)
plot(predn)


tab_model(modneed1, modc1, modn,show.ci = F, show.se = T, digits = 3)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#                                 SUPPLEMENTAL S9
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

need_box <- ggplot(needclim, aes(x=factor(year), y=mean_length, fill=drought, color="black"))+
  geom_jitter(width = .2)+
  geom_boxplot()+
  ylab("Needle length (cm)")+
  xlab("Year")+
  scale_color_manual(values=c( "grey"))+
  scale_fill_manual(values=c("#e9c46a", "#084f63"))+
  theme(legend.title = element_blank(),
        legend.position = c(.83,.15),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        legend.background = element_blank())+
  guides(color=F)

need_box

iso_box <- ggplot(isoclim, aes(x=factor(year), y=c, fill=drought, color="black"))+
  geom_jitter(width = .2)+
  geom_boxplot()+
  labs(y="Needle *&delta;*<sup>13</sup>C (&permil;)",
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

iso_box


n_box <- ggplot(isoclim, aes(x=factor(year), y=n, fill=drought, color="black"))+
  geom_jitter(width = .2)+
  geom_boxplot()+
  labs(y="Needle *&delta;*<sup>15</sup>N (&permil;)",
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

n_box

figs <- need_box / iso_box / n_box / rwifig
figs + plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')




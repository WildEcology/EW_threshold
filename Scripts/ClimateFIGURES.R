
## ---------------------------
##
## Script name: Climate data and figures
##
## Author: Dr. Joan Dudney
##
## Date Created: 2021-10-19
##
## Copyright (c) Joan Dudney, 2021
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: This code cleans climate data and produces 
##        Figure 2 panels b,c
##
## ---------------------------

library(sjPlot)
library(ggeffects)
library(patchwork)
library(tidyverse)
library(lme4)

theme_set(
  theme_bw(base_size = 20)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select


##=================================================================================================================
##  reading in and summarizing data                   
##==================================================================================================================

##annual climate data
clim_ann <- read_csv("Data/prism_growingseason.csv")
climdat<- clim_ann %>% 
  filter(year<2018) %>% 
  rename(plot=plot_id_needle)

range(climdat$tmax, na.rm=T)
mean(climdat$tmax, na.rm=T)

range(climdat$ppt, na.rm=T)
mean(climdat$ppt, na.rm=T)

## testing for duplicated data
test_clim <- climdat %>% 
  filter(duplicated(tmax))
  

##Calculating the mean values across the time period
meanmax_temps <- data.frame(mean=mean(climdat$tmax))



## Temperature data
tempsdat=climdat %>%
  select(year, plot, tmax) %>%
  group_by(year) %>% 
  summarize(mtmax=mean(tmax)) %>% 
  pivot_longer(-year)

## Precipitation data
pptdat=climdat %>%
  distinct(year, plot, ppt) %>%
  group_by(year) %>% 
  summarize(ppt=mean(ppt)) %>% 
  pivot_longer(-year)

meanppt=data.frame(mean=mean(pptdat$value)) 

##=================================================================================================================
##  FIGURES                 
##==================================================================================================================

max(tempsdat$value, na.rm=T)
range(tempsdat$value, na.rm=T)

##temp figure
temp_fig=ggplot(tempsdat, aes(x=year, y=value, color=name))+
  scale_color_manual(values = "#412234",
                     labels="",
                     name="")+
  geom_point(alpha=.9)+
  geom_line(alpha=.9)+
  geom_hline(data=meanmax_temps, aes(yintercept = mean), 
             colour = "#412234", linetype=2)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("Max. temperature (°C)")+
  guides(color="none")+
  scale_x_continuous(name="Year", breaks = scales::pretty_breaks(n = 8))+
  theme(legend.title = element_blank(),
        legend.position = c(.17,.5),
        legend.background=element_blank())+
  geom_rect(aes(xmin=2011.5, xmax=2015.5, ymin=3.493918, ymax=max(value)),
            linetype=0, alpha=.01,fill="#CCCCCC")+
  geom_rect(aes(xmin=1958.5, xmax=1961.5, ymin=3.493918, ymax=max(value)),
            linetype=0, alpha=.01,fill="#CCCCCC")+
  geom_rect(aes(xmin=1975.5, xmax=1977.5, ymin=3.493918, ymax=max(value)),
            linetype=0, alpha=.01,fill="#CCCCCC")

temp_fig


##extreme droughts
# 1959-1961 
# 1976-1977
# 2012-2015


##precip figure
precip_fig=ggplot(pptdat, aes(x=year, y=value, color=name))+
  scale_color_manual(values = "#252F7E",
                     labels="",
                     name="")+
  geom_bar(stat="identity", width=.35, alpha=.8)+
  #geom_line(alpha=.5)+
  geom_hline(data=meanppt, aes(yintercept = mean), 
             colour = "#252F7E",alpha=.8, linetype=2)+
  theme(strip.background = element_blank(),
        strip.placement = "outside")+
  ylab("Mean precipitation (mm)")+
  guides(color="none")+
  scale_x_continuous(name="Year", breaks = scales::pretty_breaks(n = 8))+
  theme(legend.title = element_blank(),
        legend.position = c(.17,.5),
        legend.background=element_blank())+
  geom_rect(aes(xmin=2011.5, xmax=2015.5, ymin=0, ymax=max(value)),
            linetype=0, alpha=.01,fill="#CCCCCC")+
  geom_rect(aes(xmin=1958.5, xmax=1961.5, ymin=0, ymax=max(value)),
            linetype=0, alpha=.01,fill="#CCCCCC")+
  geom_rect(aes(xmin=1975.5, xmax=1977.5, ymin=0, ymax=max(value)),
            linetype=0, alpha=.01,fill="#CCCCCC")


precip_fig


temp_fig / precip_fig
#plot_annotation(tag_levels = "a", tag_prefix = '(',tag_suffix = ')')


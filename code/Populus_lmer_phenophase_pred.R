# this script is for Populus phenology modeling among species and phenophases, and future predictions

library(tidyverse)
library(sf)
library(cowplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(cAIC4)

# read data
pop_npn.df3=read_csv("data/pop_NPN_modeling_data.csv")

summary(pop_npn.df3)


# remove outliers
quant = pop_npn.df3[,c(7,10,12)] %>% group_by(Species) %>% 
  summarize(quant25=quantile(doy, probs=c(.25), na.rm = F),
            quant75=quantile(doy, probs=c(.75), na.rm = F),
            IQR=IQR(doy)) %>% 
  mutate(lower=quant25 - 1.5*IQR, upper=quant75 + 1.5*IQR)

pop_npn.df3 = left_join(pop_npn.df3, quant[, c(1:2, 5:6)])
pop_npn.df4 = pop_npn.df3 %>% filter(doy<upper, doy>lower)

# remove trichocarpa due to small sample size
pop_npn.df4=pop_npn.df4 %>% filter(Species != "trichocarpa")

colnames(pop_npn.df4)
cor(pop_npn.df4[, c(12, 3, 5, 45, 47)])


# model selection
m1 <- lmer(doy ~ scale(lat) + scale(elevation)+
             scale(Tave_sp_anm)*scale(PPT_sp_anm) +
             (1+scale(Tave_sp_anm)|spp:Phenophase_ID), 
           data = pop_npn.df4 , REML=T, control = lmerControl(optimizer="bobyqa"))
summary(m1)

m2 <- lmer(doy ~ scale(lat) + scale(elevation)+
             scale(Tmin_sp_anm)*scale(PPT_sp_anm) +
             (1+scale(Tmin_sp_anm)|spp:Phenophase_ID), 
           data = pop_npn.df4 , REML=T, control = lmerControl(optimizer="bobyqa"))
summary(m2)

m3 <- lmer(doy ~ scale(lat) + scale(elevation)+
             scale(Tmax_sp_anm)*scale(PPT_sp_anm) +
             (1+scale(Tmax_sp_anm)|spp:Phenophase_ID), 
           data = pop_npn.df4 , REML=T, control = lmerControl(optimizer="bobyqa"))
summary(m3)


AIC(m1, m2, m3)


m3a <- lmer(doy ~ scale(lat) + scale(elevation)+
              scale(Tmax_sp_anm)+scale(PPT_sp_anm) +
              (1+scale(Tmax_sp_anm)|Phenophase_ID:spp) , 
            data = pop_npn.df4 , REML=T, control = lmerControl(optimizer="bobyqa"))
summary(m3a)
coef(m3a)
r.squaredGLMM(m3a)

cf_npn=coef(m3a)$`Phenophase_ID:spp`
cf_npn$Phenophase_ID=substr(rownames(cf_npn), 1, 3)
cf_npn$spp=substr(rownames(cf_npn), 5, 30)
cf_npn$cf_us=cf_npn$`scale(Tmax_sp_anm)`/sd(pop_npn.df4$Tmax_sp_anm)


# Figure 5
pop_npn.df4$spp1=paste0("P. ", pop_npn.df4$Species, sep="")
pop_npn.df4$spp1=factor(pop_npn.df4$spp1, levels=c("P. deltoides","P. tremuloides","P. grandidentata","P. balsamifera"))

p1=ggplot(pop_npn.df4, aes(x=spp1, y=doy, fill=factor(Phenophase_ID)))+
  geom_boxplot()+
  scale_fill_manual(name="phenophase", values=c("#659794","#F5C98E"), labels=c("first leaf out","first open flower"))+
  labs(y="day of year", x="species")+
  theme_bw(16)+
  theme(axis.text.y = element_text(face = "italic"))+
  coord_flip()

cf_npn$spp1=substr(cf_npn$spp, 9, 25)
cf_npn$spp1=paste0("P. ", cf_npn$spp1, sep="")
cf_npn$spp1=factor(cf_npn$spp1, levels=c("P. deltoides","P. tremuloides","P. grandidentata","P. balsamifera"))

ylab=expression("phenological" ~ "sensitivity" ~ "(day" ~ paste("\u00b0C" ^ "-1", ")"))

ggplot(cf_npn, aes(y=cf_us, x=spp1, shape=factor(Phenophase_ID)))+
  geom_point(position=position_dodge(width=0.3), size=2)+
  scale_shape_manual(name="phenophase", values=c(24, 19), labels=c("first leaf","first open flower"))+
  theme_bw(14)+ 
  labs(y=ylab, x="species")+
  coord_flip()

p2=ggplot(cf_npn, aes(y=cf_us, x=spp1, col=factor(Phenophase_ID)))+
  geom_point(size=3)+
  scale_color_manual(name="phenophase", values=c("#659794","#F5C98E"), labels=c("first leaf out","first open flower"))+
  theme_bw(16)+ 
  theme(axis.text.y = element_text(face = "italic"))+
  labs(y=ylab, x="")+
  coord_flip()

plots=cowplot::plot_grid(p1, p2, labels=c("(a)","(b)"), label_size=16, ncol=2, nrow=1, 
                         rel_widths=c(1, 1))

ggsave(filename = "figures/NPN_populus_dates_cf.tif", plots, device="tiff",
       dpi=300, height= 4, width=15, compression="lzw")



# future prediction -----------------------------------------------
npn.LT=read_csv(file="data_processed/NPN_populus_longtermavg.csv")

# prepare new data
# spatial and random effect groups
npn_sp=pop_npn.df4 %>% dplyr::select(lat, lon, elevation, spp, Phenophase_ID) %>% distinct()


# 30-y normal 1981-2010
npn_normals=read_csv(file="data_raw/npn_locations_Normal_1981_2010SY.csv")
colnames(npn_normals)

npn_normals1 = npn_normals[, c(1:5, 7, 19)] %>% 
  distinct() %>%
  rename(spp=ID2, lat=Latitude, lon=Longitude, elevation=Elevation) 

npn_normals1 = left_join(npn_normals1, npn.LT) %>% 
  mutate(Tmax_sp_anm = Tmax_sp-LT_Tmax_sp,
         PPT_sp_anm = (PPT_sp-LT_PPT_sp)/LT_PPT_sp)

npn.nd = left_join(npn_sp, npn_normals1[, c(2:5, 18:19)])

summary(npn.nd)


# 30-y normal 1981-2010 estimation
npn.nd$fit.norm=predict(m3a, newdata=npn.nd)


# GCM data from ClimateNA
npn_8GCMs=read_csv(file="data_raw/npn_locations_4 GCMsSY_8ens_30y.csv")
colnames(npn_8GCMs)

npn_8GCMs1 = npn_8GCMs[, c(1:6, 8, 20)] %>% 
  distinct() %>%
  rename(year=Year, spp=ID2, lat=Latitude, lon=Longitude, elevation=Elevation) 

npn_8GCMs1 = left_join(npn_8GCMs1, npn.LT) %>% 
  mutate(Tmax_sp_anm = Tmax_sp-LT_Tmax_sp,
         PPT_sp_anm = (PPT_sp-LT_PPT_sp)/LT_PPT_sp)

unique(npn_8GCMs1$year)

npnpd1 = left_join(npn_sp, npn_8GCMs1[npn_8GCMs1$year=="8GCMs_ensemble_ssp245_2041-2070.gcm", c(3:6, 19:20)])
npn.nd$ssp245_4170=predict(m3a, newdata=npnpd1)

npnpd2 = left_join(npn_sp, npn_8GCMs1[npn_8GCMs1$year=="8GCMs_ensemble_ssp245_2071-2100.gcm", c(3:6, 19:20)])
npn.nd$ssp245_7100=predict(m3a, newdata=npnpd2)

npnpd3 = left_join(npn_sp, npn_8GCMs1[npn_8GCMs1$year=="8GCMs_ensemble_ssp585_2041-2070.gcm", c(3:6, 19:20)])
npn.nd$ssp585_4170=predict(m3a, newdata=npnpd3)

npnpd4 = left_join(npn_sp, npn_8GCMs1[npn_8GCMs1$year=="8GCMs_ensemble_ssp585_2071-2100.gcm", c(3:6, 19:20)])
npn.nd$ssp585_7100=predict(m3a, newdata=npnpd4)

summary(npn.nd)

npn.nd=read_csv(file="data_processed/pop_npn_predictions_20221111.csv")
npn.nd$spp1=substr(npn.nd$spp, 9, 25)
npn.nd$spp1=paste0("P. ", npn.nd$spp1, sep="")

npn.nd1=npn.nd %>% pivot_longer(cols=8:12, names_to = "scenario")


# calculate the time gap between phenophase
npn.p = npn.nd1[, c(1:3,5, 8:10)] %>% 
  pivot_wider(names_from = c(scenario, Phenophase_ID), values_from = value)

npn.p$gap_fit.norm=npn.p$fit.norm_483-npn.p$fit.norm_501
npn.p$gap_ssp245_4170=npn.p$ssp245_4170_483-npn.p$ssp245_4170_501
npn.p$gap_ssp245_7100=npn.p$ssp245_7100_483-npn.p$ssp245_7100_501
npn.p$gap_ssp585_4170=npn.p$ssp585_4170_483-npn.p$ssp585_4170_501
npn.p$gap_ssp585_7100=npn.p$ssp585_7100_483-npn.p$ssp585_7100_501


npn.p1 = npn.p[, c(1:4, 15:19)] %>% pivot_longer(cols=5:9, names_to="scenario")
order1 = npn.p1 %>% group_by(spp1) %>% summarise(value.m=median(value, na.rm=T)) %>% arrange(value.m)
npn.p1$spp1=factor(npn.p1$spp1, levels=order1$spp1)

# Figure S6
tiff(filename = "figures/prediction_gap_NPN_spp4_boxplot.tif", 
     res=300, height= 1000, width=2400, compression="lzw")
ggplot(npn.p1, aes(x=spp1, y=value, color=scenario), fill="white")+
  geom_boxplot()+
  scale_color_manual(values=c("dark green","orange","dark orange3","red","darkred"),
                     labels=c("1981-2010","2041-2070 SSP2-4.5","2071-2100 SSP2-4.5","2041-2070 SSP5-8.5","2071-2100 SSP5-8.5"))+
  labs(x="species", y="temporal gaps (days)")+
  theme_bw(14)+
  coord_flip()+
  theme(axis.text.y = element_text(face = "italic"))
dev.off()


npn.p$dg_ssp245_4170=npn.p$gap_ssp245_4170-npn.p$gap_fit.norm
npn.p$dg_ssp245_7100=npn.p$gap_ssp245_7100-npn.p$gap_fit.norm
npn.p$dg_ssp585_4170=npn.p$gap_ssp585_4170-npn.p$gap_fit.norm
npn.p$dg_ssp585_7100=npn.p$gap_ssp585_7100-npn.p$gap_fit.norm

npn.p %>% group_by(spp1) %>% summarize(dg.mean=mean(dg_ssp585_7100, na.rm=T), 
                                       dg.min=min(dg_ssp585_7100, na.rm=T), 
                                       dg.max=max(dg_ssp585_7100, na.rm=T))



npn.p2 = npn.p[, c(1:4, 20:23)] %>% pivot_longer(cols=5:8, names_to="scenario")
order1 = npn.p2 %>% group_by(spp1) %>% summarise(value.m=median(value, na.rm=T)) %>% arrange(value.m)
npn.p2$spp1=factor(npn.p2$spp1, levels=order1$spp1)

# Figure S7
tiff(filename = "figures/prediction_change_gap_NPN_spp4_boxplot.tif", 
     res=300, height= 1000, width=2400, compression="lzw")
ggplot(npn.p2, aes(x=spp1, y=value, color=scenario), fill="white")+
  geom_boxplot()+
  scale_color_manual(values=c("orange","dark orange3","red","darkred"),
                     labels=c("2041-2070 SSP2-4.5","2071-2100 SSP2-4.5","2041-2070 SSP5-8.5","2071-2100 SSP5-8.5"))+
  labs(x="species", y="change of temporal gaps (days)")+
  theme_bw(14)+
  coord_flip()+
  theme(axis.text.y = element_text(face = "italic"))
dev.off()


# calculate the proportion of change in gaps
npn.p$pg_ssp245_4170=npn.p$dg_ssp245_4170/npn.p$gap_fit.norm
npn.p$pg_ssp245_7100=npn.p$dg_ssp245_7100/npn.p$gap_fit.norm
npn.p$pg_ssp585_4170=npn.p$dg_ssp585_4170/npn.p$gap_fit.norm
npn.p$pg_ssp585_7100=npn.p$dg_ssp585_7100/npn.p$gap_fit.norm


npn.p3 = npn.p[, c(1:4, 24:27)] %>% pivot_longer(cols=5:8, names_to="scenario")

order1 = npn.p3 %>% group_by(spp1) %>% summarise(value.m=median(value, na.rm=T)) %>% arrange(value.m)
npn.p3$spp1=factor(npn.p3$spp1, levels=order1$spp1)

# Figure 6
tiff(filename = "figures/proportion_change_gap_NPN_spp4_boxplot.tif", 
     res=300, height= 1000, width=2400, compression="lzw")
ggplot(npn.p3, aes(x=spp1, y=value, color=scenario), fill="white")+
  geom_boxplot()+
  scale_color_manual(values=c("orange","dark orange3","red","darkred"),
                     labels=c("2041-2070 SSP2-4.5","2071-2100 SSP2-4.5","2041-2070 SSP5-8.5","2071-2100 SSP5-8.5"))+
  labs(x="species", y="proportional change of temporal gaps")+
  theme_bw(14)+
  coord_flip()+
  theme(axis.text.y = element_text(face = "italic"))
dev.off()



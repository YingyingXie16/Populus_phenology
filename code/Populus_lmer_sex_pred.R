# this script is for Populus phenology modeling among species, sex, and ecoregions, and future predictions

library(tidyverse)
library(sf)
library(cowplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(cAIC4)


# read data
pop_fl1 = read_csv(file="data/populus_modeling_data_20221031.csv")
summary(pop_fl1)

# remove species fremontii, which is special!
pop_fl2=pop_fl1 %>% filter(species !="Populus fremontii")


# by ecoregion
pop_fl3=pop_fl2 %>% filter(NA_L2CODE !=0) %>% na.omit()


# select region with at least 3 observations for both sex
t= as.data.frame(table(pop_fl3[, c("species","sex","NA_L2CODE")]))
t=t[t$Freq>0,]
t = mutate(t, spp_ecor=paste(species, NA_L2CODE, sep="_"))
t = filter(t, Freq>=3)
tn = t %>% group_by(species, NA_L2CODE) %>% summarise(n=n(), nFreq=sum(Freq))
t = left_join(t, tn)
t1 = filter(t, n==2)

pop_fl3 = pop_fl3 %>% mutate(spp_ecor=paste(species, NA_L2CODE, sep="_")) %>% 
  filter(spp_ecor %in% t1$spp_ecor)

colnames(pop_fl3)

# model selection
pm1 <- lmer(doy ~ scale(lat) + scale(elevation) +
              scale(Tmin_wt_anm)*scale(PPT_wt_anm) +
              (1+scale(Tmin_wt_anm)|species:NA_L2CODE:sex), 
            data = pop_fl3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(pm1)


pm2 <- lmer(doy ~ scale(lat) + scale(elevation) +
              scale(Tmax_wt_anm)*scale(PPT_wt_anm) +
              (1+scale(Tmax_wt_anm)|species:NA_L2CODE:sex), 
            data = pop_fl3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(pm2)


pm3 <- lmer(doy ~ scale(lat) + scale(elevation) +
              scale(Tmin_sp_anm)*scale(PPT_sp_anm) +
              (1+scale(Tmin_sp_anm)|species:NA_L2CODE:sex), 
            data = pop_fl3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(pm3)


pm4 <- lmer(doy ~ scale(lat) + scale(elevation) +
              scale(Tmax_sp_anm)*scale(PPT_sp_anm) +
              (1+scale(Tmax_sp_anm)|species:NA_L2CODE:sex), 
            data = pop_fl3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(pm4)

AICc(pm1, pm2, pm3, pm4)


pm4a <- lmer(doy ~ scale(lat) + scale(elevation) +
               scale(Tmax_sp_anm)*scale(PPT_sp_anm) +
               (1+scale(Tmax_sp_anm)|species:NA_L2CODE:sex), 
             data = pop_fl3, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(pm4a)

pm4b <- lmer(doy ~ scale(lat) + scale(elevation) +
               scale(Tmax_sp_anm)+scale(PPT_sp_anm) +
               (1+scale(Tmax_sp_anm)|species:NA_L2CODE:sex), 
             data = pop_fl3, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(pm4b)


AICc(pm4a, pm4b)


pme <- lmer(doy ~ scale(lat) + scale(elevation) +
              scale(Tmax_sp_anm)+scale(PPT_sp_anm) +
              (1+scale(Tmax_sp_anm)|species:NA_L2CODE:sex),  
            data = pop_fl3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(pme)
coef(pme)
r.squaredGLMM(pme)

# extract coefficients
cf2=coef(pme)$`species:NA_L2CODE:sex`
cf2$sex=sub(".*:", "", rownames(cf2))
cf2$species=sub("\\:.*", "", rownames(cf2))
cf2$NA_L2CODE=gsub(".*:(.*)\\:.*", "\\1", rownames(cf2))
cf2


save(pme, cf2, file="results/modeling_sppregionsex_20221031.Rdata")


# future prediction ------------------------------------------------
# read ClimateNA future projection data
dir1="data directory"
df_8GCMs=read_csv( file=paste0(dir1, "/county_average_ensemble_8GCMs_4var_format.csv"))
df_8GCMs1=filter(df_8GCMs, stat=="mean")
df_normals=read_csv(file=paste0(dir1, "/county_average_ensemble_1981_2010.csv"))
df_normals1=filter(df_normals, stat=="mean")
df_normals1$sc=paste0(df_normals1$state,"_", df_normals1$county,sep="")

# calculate climatic anomalies for future climate projection data
# 1901-2019 average seasonal temp and ppt
df1 = read_csv("data_raw/UScounty_CHELSA_winter_spring_tmax_tmin_ppt_1901_2019.csv")
df2 = read_csv("data_raw/CAcounty_CHELSA_winter_spring_tmax_tmin_ppt_1901_2019.csv")

unique(df1$sc)
df1$datatype="CHELSAv2.1"
df1$datatype[df1$year <2017]="CHELSAcruts"
df2$datatype="CHELSAv2.1"
df2$datatype[df2$year <2017]="CHELSAcruts"

df=rbind.data.frame(df1, df2)
colnames(df)

# calculate long term climate variables
LT=df %>% group_by(sc) %>% summarise(LT_Tmax_sp=mean(tmax_sp, na.rm=T),
                                     LT_Tmin_wt=mean(tmin_wt, na.rm=T),
                                     LT_PPT_wt=mean(ppt_wt, na.rm=T),
                                     LT_PPT_sp=mean(ppt_sp, na.rm=T)) 
summary(LT)


# spatial and random effect groups
df_sp=pop_fl3 %>% dplyr::select(lat, lon, elevation, sc, species, sex, NA_L2CODE) %>% distinct()


# 30-y normal 1981-2010
df_normals1 = left_join(df_normals1, LT) %>% 
  mutate(Tmax_sp_anm = normal_1981_2010_Tmax_sp-LT_Tmax_sp,
         PPT_sp_anm = (normal_1981_2010_PPT_sp-LT_PPT_sp)/LT_PPT_sp)
df.nd = left_join(df_sp, df_normals1[, c(9, 14:15)])

summary(df.nd)


# 30-y normal 1981-2010 estimation
df.nd$fit.norm=predict(pme, newdata=df.nd)


# 30-y future prediction
df_8GCMs1 = left_join(df_8GCMs1, LT) %>% 
  mutate(Tmax_sp_anm = Tmax_sp-LT_Tmax_sp,
         PPT_sp_anm = (PPT_sp-LT_PPT_sp)/LT_PPT_sp)


# 8GCM-ensemble-ssp245-2041-2070
dfpd1 = left_join(df.nd[, c(1:7)], filter(df_8GCMs1, scenario=="ssp245_2041_2070")) %>% distinct()
df.nd$ssp245_4170=predict(pme, newdata=dfpd1)

# 8GCM-ensemble-ssp585-2041-2070
dfpd2 = left_join(df.nd[, c(1:7)], filter(df_8GCMs1, scenario=="ssp585_2041_2070")) %>% distinct()
df.nd$ssp585_4170=predict(pme, newdata=dfpd2)

# 8GCM-ensemble-ssp245-2071-2100
dfpd3 = left_join(df.nd[, c(1:7)], filter(df_8GCMs1, scenario=="ssp245_2071_2100")) %>% distinct()
df.nd$ssp245_7100=predict(pme, newdata=dfpd3)

# 8GCM-ensemble-ssp585-2071-2100
dfpd4 = left_join(df.nd[, c(1:7)], filter(df_8GCMs1, scenario=="ssp585_2071_2100")) %>% distinct()
df.nd$ssp585_7100=predict(pme, newdata=dfpd4)

df.nd=na.omit(df.nd)

# long format predictions
df.pred=df.nd[, c(1:7, 10:14)] %>% pivot_longer(cols=c(8:12), names_to = "scenario", values_to = "pred")


scn5=c("fit.norm","ssp245_4170","ssp245_7100","ssp585_4170","ssp585_7100")

ggplot(df.pred, aes(x=species, y=pred, fill=sex))+
  geom_boxplot(position = position_dodge2(preserve = "single"))+
  facet_wrap(~scenario, ncol=1)+
  theme_bw()

# change in phenology
df.nd$diff_spp245_4170=df.nd$ssp245_4170-df.nd$fit.norm
df.nd$diff_spp245_7100=df.nd$ssp245_7100-df.nd$fit.norm
df.nd$diff_spp585_4170=df.nd$ssp585_4170-df.nd$fit.norm
df.nd$diff_spp585_7100=df.nd$ssp585_7100-df.nd$fit.norm


df.peco=df.pred %>% group_by(scenario, species, sex, NA_L2CODE) %>% 
  summarise(n=n(), pred.m=mean(pred), pred.md=median(pred), pred.sd=sd(pred))


# calculate the time gap between sex
df.peco2 = df.peco[, c(1:4, 6)] %>% 
  pivot_wider(names_from = c(scenario, sex), values_from = pred.m)

df.peco2$gap_fit.norm=df.peco2$fit.norm_F-df.peco2$fit.norm_M
df.peco2$gap_ssp245_4170=df.peco2$ssp245_4170_F-df.peco2$ssp245_4170_M
df.peco2$gap_ssp245_7100=df.peco2$ssp245_7100_F-df.peco2$ssp245_7100_M
df.peco2$gap_ssp585_4170=df.peco2$ssp585_4170_F-df.peco2$ssp585_4170_M
df.peco2$gap_ssp585_7100=df.peco2$ssp585_7100_F-df.peco2$ssp585_7100_M


df.peco3 = df.peco2[, c(1,2, 13:17)] %>% pivot_longer(cols=3:7, names_to="scenario")
order1 = df.peco3 %>% group_by(spp1) %>% summarise(value.m=median(value, na.rm=T)) %>% arrange(value.m)
df.peco3$spp1=factor(df.peco3$spp1, levels=order1$spp1)

# Figure S3
tiff(filename = "figures/prediction_gap_spp7_boxplot.tif", 
     res=300, height= 1500, width=2500, compression="lzw")
ggplot(df.peco3, aes(x=spp1, y=value, color=scenario), fill="white")+
  geom_boxplot()+
  scale_color_manual(values=c("dark green","orange","dark orange3","red","darkred"),
                     labels=c("1981-2010","2041-2070 SSP2-4.5","2071-2100 SSP2-4.5","2041-2070 SSP5-8.5","2071-2100 SSP5-8.5"))+
  labs(x="species", y="temporal gaps (days)")+
  theme_bw(14)+
  coord_flip()+
  theme(axis.text.y = element_text(face = "italic"))
dev.off()



# calculate change of gap
df.peco2$dg_ssp245_4170=df.peco2$gap_ssp245_4170-df.peco2$gap_fit.norm
df.peco2$dg_ssp245_7100=df.peco2$gap_ssp245_7100-df.peco2$gap_fit.norm
df.peco2$dg_ssp585_4170=df.peco2$gap_ssp585_4170-df.peco2$gap_fit.norm
df.peco2$dg_ssp585_7100=df.peco2$gap_ssp585_7100-df.peco2$gap_fit.norm

# calculate the proportion of change in gaps
df.peco2$pg_ssp245_4170=df.peco2$dg_ssp245_4170/df.peco2$gap_fit.norm
df.peco2$pg_ssp245_7100=df.peco2$dg_ssp245_7100/df.peco2$gap_fit.norm
df.peco2$pg_ssp585_4170=df.peco2$dg_ssp585_4170/df.peco2$gap_fit.norm
df.peco2$pg_ssp585_7100=df.peco2$dg_ssp585_7100/df.peco2$gap_fit.norm


df.peco4 = df.peco2[, c(1,26, 18:21)] %>% pivot_longer(cols=3:6, names_to="scenario")
quantile(df.peco4$value, seq(0,1, 0.05))
df.peco4 %>% filter(scenario=="dg_ssp585_7100") %>% summary()

order1 = df.peco4 %>% group_by(spp1) %>% summarise(value.m=median(value, na.rm=T)) %>% arrange(value.m)
df.peco4$spp1=factor(df.peco4$spp1, levels=order1$spp1)

# Figure S4
tiff(filename = "figures/prediction_change_gap_spp7_boxplot.tif", 
     res=300, height= 1500, width=2500, compression="lzw")
ggplot(df.peco4, aes(x=spp1, y=value, color=scenario), fill="white")+
  geom_boxplot()+
  scale_color_manual(values=c("orange","dark orange3","red","darkred"),
                     labels=c("2041-2070 SSP2-4.5","2071-2100 SSP2-4.5","2041-2070 SSP5-8.5","2071-2100 SSP5-8.5"))+
  labs(x="species", y="change of temporal gaps (days)")+
  theme_bw(14)+
  coord_flip()+
  theme(axis.text.y = element_text(face = "italic"))
dev.off()


df.peco5 = df.peco2[, c(1,2, 22:25)] %>% pivot_longer(cols=3:6, names_to="scenario")
df.peco5$spp1=substr(df.peco5$species, 9, 25)
df.peco5$spp1=paste0("P. ", df.peco5$spp1, sep="")
summary(df.peco5)

order1 = df.peco5 %>% group_by(spp1) %>% summarise(value.m=median(value)) %>% arrange(value.m)
df.peco5$spp1=factor(df.peco5$spp1, levels=order1$spp1)

# Figure 4
tiff(filename = "figures/proportion_change_gap_spp7_boxplot.tif", 
     res=300, height= 1500, width=2500, compression="lzw")
ggplot(df.peco5 , aes(x=spp1, y=value, col=scenario), fill="white")+
  geom_boxplot()+
  #stat_boxplot(width = 0.15) +
  #scale_fill_manual(name="scenario", values=c("orange","brown","red","darkred"),
  #                  labels=c("2041-2070 SSP2-4.5","2071-2100 SSP2-4.5","2041-2070 SSP5-8.5","2071-2100 SSP5-8.5"))+
  scale_color_manual(values=c("orange","dark orange3","red","darkred"),
                     labels=c("2041-2070 SSP2-4.5","2071-2100 SSP2-4.5","2041-2070 SSP5-8.5","2071-2100 SSP5-8.5"))+
  geom_hline(yintercept = 0, lty=2)+
  labs(x="species", y="proportional changes in gaps")+
  theme_bw(14)+
  coord_flip()+
  theme(axis.text.y = element_text(face = "italic"))
dev.off()




# Populus fremontii ------------------------------------------------
pop_fr = pop_fl1 %>% filter(species == "Populus fremontii")

# by ecoregion
t=as.data.frame(table(pop_fr[, c("sex","NA_L2CODE")]))

t=t[t$Freq>3,]

t = mutate(t, ecor_sex=paste(NA_L2CODE, sex, sep="_"))

pop_fr2 = pop_fr %>% mutate(ecor_sex=paste(NA_L2CODE, sex, sep="_")) %>% 
  filter(ecor_sex %in% t$ecor_sex)

pop_fr3 = filter(pop_fr2, NA_L2CODE %in% c("8.3", "10.1", "10.2", "11.1"))

# model selection
frem1 <- lmer(doy ~ scale(lat)+scale(elevation)+
                scale(Tmin_wt_anm)*scale(PPT_wt_anm) +
                (1+scale(Tmin_wt_anm)|NA_L2CODE:sex), 
              data = pop_fr3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(frem1)


frem2 <- lmer(doy ~ scale(lat)+scale(elevation)+
                scale(Tmax_wt_anm)*scale(PPT_wt_anm) +
                (1+scale(Tmax_wt_anm)|NA_L2CODE:sex), 
              data = pop_fr3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(frem2)


frem3 <- lmer(doy ~ scale(lat)+scale(elevation)+
                scale(Tmin_sp_anm)*scale(PPT_sp_anm) +
                (1+scale(Tmin_sp_anm)|NA_L2CODE:sex), 
              data = pop_fr3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(frem3)


frem4 <- lmer(doy ~ scale(lat)+scale(elevation)+
                scale(Tmax_sp_anm)*scale(PPT_sp_anm) +
                (1+scale(Tmax_sp_anm)|NA_L2CODE:sex),  
              data = pop_fr3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(frem4)


AICc(frem1, frem2, frem3, frem4)

frem1a <- lmer(doy ~ scale(elevation)+
                 scale(Tmin_wt_anm)+scale(PPT_wt_anm) +
                 (1|NA_L2CODE:sex), 
               data = pop_fr3, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frem1a)


frem1b <- lmer(doy ~ scale(lat)+scale(elevation)+
                 scale(PPT_wt_anm) +
                 (1+scale(Tmin_wt_anm)|NA_L2CODE:sex), 
               data = pop_fr3, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frem1b)

frem1c <- lmer(doy ~ scale(elevation)+
                 scale(PPT_wt_anm) +
                 (1+scale(PPT_wt_anm)|NA_L2CODE:sex), 
               data = pop_fr3, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frem1c)

frem1d <- lmer(doy ~ scale(elevation)+
                 scale(PPT_wt_anm) +
                 (1|NA_L2CODE:sex), 
               data = pop_fr3, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frem1d)

AICc(frem1a, frem1b, frem1c, frem1d)

# best model
frem <- lmer(doy ~ scale(elevation)+
               scale(PPT_wt_anm) +
               (1|sex:NA_L2CODE),  
             data = pop_fr3, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(frem)
coef(frem)
r.squaredGLMM(frem)
# R2m       R2c
# [1,] 0.3910712 0.4720987

cffr1=coef(frem)$`sex:NA_L2CODE`
cffr1$NA_L2CODE=substr(rownames(cffr1),3, 6)
cffr1$sex=substr(rownames(cffr1),1, 1)
cffr1

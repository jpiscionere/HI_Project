
H0=100
DEG_TO_RAD=0.01745328888
RAD_TO_DEG=57.2957914331

library('ggplot2')

library("astrolibR")

library(stringi)

data_groups_full=read.csv("../Data/2mass_groups.csv",skip=4)

data=read.csv("../Data/HOPCAT2005TabSep.csv")

df=data.frame(data)
df_HIPASS=na.omit(df)

coords=glactc(gl = data_groups_full$Galactic.longitude,gb = data_groups_full$Galactic.latitude,j=2,year = 2000,degree = TRUE)

data_groups_full$RA=coords$ra
data_groups_full$Dec=coords$dec
df_2MASS_South<-subset(data_groups_full,Dec < 0)

coords=glactc(gl = df_HIPASS$HicatExtl,gb = df_HIPASS$HicatExtb,j=2,year = 2000,degree = TRUE)
df_HIPASS$RA=coords$ra
df_HIPASS$Dec=coords$dec

length(df_HIPASS$Morphology)
length(which(df_HIPASS$Bj_Mag_AUTO_Calibrated != "XXXXX"))
length(which(df_HIPASS$Morphology != "XXXXX"))

df_HIPASS$Bj_Mag_AUTO_Calibrated[which(df_HIPASS$Bj_Mag_AUTO_Calibrated == "XXXXX")]=NA
df_HIPASS=na.omit(df_HIPASS)
length(df_HIPASS$Bj_Mag_AUTO_Calibrated)
typeof(df_HIPASS$Bj_Mag_AUTO_Calibrated)
df_HIPASS$Morphology=stri_replace_all_fixed(df_HIPASS$Morphology, " ", "")

galaxy_output_df=data.frame(ID=df_HIPASS$ID,
                  RA=df_HIPASS$RA,
                  Dec=df_HIPASS$Dec,
                  Velocity_vel_mon=df_HIPASS$vel_mom,
                  Bj_mag=df_HIPASS$Bj_Mag_AUTO_Calibrated,
                  R_mag=df_HIPASS$Red_Mag_AUTO_Calibrated,
                  I_mag=df_HIPASS$I_Mag_AUTO_Calibrated,
                  S_int=df_HIPASS$Sint,
                  Morphology=df_HIPASS$Morphology)

write.table(galaxy_output_df,"HIPASS_clean.txt",sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)

summary(df_2MASS_South)

max(df_2MASS_South$Projected.second.turnaround.radius.anticipated.by.corrected.intrinsic.luminosity)
min(df_2MASS_South$Distance.D.Vgp.100.Mpc)

asin(DEG_TO_RAD*5.225/(2*DEG_TO_RAD*0.6)) * 2 * RAD_TO_DEG

Maximum_Dec_Separation=asin(DEG_TO_RAD*max(df_2MASS_South$Projected.second.turnaround.radius.anticipated.by.corrected.intrinsic.luminosity) / 
                            (2*DEG_TO_RAD*min(df_2MASS_South$Distance.D.Vgp.100.Mpc)))*2.*RAD_TO_DEG

group_output_df=data.frame(RA=df_2MASS_South$RA,
                          Dec=df_2MASS_South$Dec,
                          Velocity=df_2MASS_South$Velocity.in.CMB.frame..cosmological.adjustment.described.in.CF2.data.paper,
                          Distance=df_2MASS_South$Distance.D.Vgp.100.Mpc,
                          Radius_Turnaround=df_2MASS_South$Projected.second.turnaround.radius.anticipated.by.corrected.intrinsic.luminosity,
                          N_group=df_2MASS_South$Number.of.group..nest..members.in.2MRS.K.11.75)

write.table(group_output_df,"2MASS_Groups_clean.txt",sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

max(group_output_df$Radius_Turnaround)

length(which(df_2MASS_South$Bi.weight.projected.virial.radius.Rij > 0))



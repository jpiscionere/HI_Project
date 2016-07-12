
H0=100
DEG_TO_RAD=0.01745328888
RAD_TO_DEG=57.2957914331

library('ggplot2')

data_groups_full=read.csv("../Data/2mass_groups.csv",skip=4)


summary(data_groups_full)

library("astrolibR")

data=read.csv("../Data/HOPCAT2005TabSep.csv")

df=data.frame(data)

df_HIPASS=na.omit(df)

summary(df_HIPASS)

coords=glactc(gl = data_groups_full$Galactic.longitude,gb = data_groups_full$Galactic.latitude,j=2,year = 2000,degree = TRUE)

data_groups_full$RA=coords$ra
data_groups_full$Dec=coords$dec

df_2MASS_South<-subset(data_groups_full,Dec < 0)

df_2MASS_South$Distance=df_2MASS_South$Velocity.in.CMB.frame..cosmological.adjustment.described.in.CF2.data.paper/H0

df_2MASS_South$X=df_2MASS_South$Distance * sin((90-df_2MASS_South$Dec) * DEG_TO_RAD)*cos(df_2MASS_South$RA * DEG_TO_RAD) 
df_2MASS_South$Y=df_2MASS_South$Distance * sin((90-df_2MASS_South$Dec) * DEG_TO_RAD)*sin(df_2MASS_South$RA * DEG_TO_RAD) 
df_2MASS_South$Z=df_2MASS_South$Distance * cos((90-df_2MASS_South$Dec) * DEG_TO_RAD)

df_HIPASS$Distance=df_HIPASS$vel_mom/H0
length(df_HIPASS$Distance)

coords=glactc(gl = df_HIPASS$HicatExtl,gb = df_HIPASS$HicatExtb,j=2,year = 2000,degree = TRUE)
qplot(coords$ra,coords$dec)
df_HIPASS$RA=coords$ra
df_HIPASS$Dec=coords$dec

df_HIPASS$M_HI=2.365*10^5*(df_HIPASS$vel_mom/H0)^2*df_HIPASS$Sint

df1=data.frame(distance=df_2MASS_South$Distance,survey=rep("2MASS",times=length(df_2MASS_South$Distance)))
df2=data.frame(distance=df_HIPASS$Distance,survey=rep("HIPASS",times=length(df_HIPASS$Distance)))
df.1=rbind(df1,df2)
ggplot(df.1,aes(distance,fill=survey,colour=survey)) + geom_density(alpha = 0.1)  + xlab("Distance (Mpc/h)") + 
    theme_bw() + scale_fill_discrete(guide_legend(title = "Survey")) + scale_color_discrete(guide_legend(title = "Survey"))

df_HIPASS$X=df_HIPASS$Distance * sin((90-df_HIPASS$Dec) * DEG_TO_RAD)*cos(df_HIPASS$RA * DEG_TO_RAD) 
df_HIPASS$Y=df_HIPASS$Distance * sin((90-df_HIPASS$Dec) * DEG_TO_RAD)*sin(df_HIPASS$RA * DEG_TO_RAD) 
df_HIPASS$Z=df_HIPASS$Distance * cos((90-df_HIPASS$Dec) * DEG_TO_RAD)

length(which(df_2MASS_South$Distance < max(df_HIPASS$Distance)))
summary(df_HIPASS$Distance)
length(which(df_2MASS_South$Distance < 43.470 & df_2MASS_South$Number.of.group..nest..members.in.2MRS.K.11.75 >1))

belong_to_group=c(1:length(df_HIPASS$RA))*0
belong_to_galaxy=c(1:length(df_2MASS_South$RA)) * 0
group_HI_mass=c(1:length(df_2MASS_South$RA))*0

df_2MASS_South$Radius_SQR=df_2MASS_South$Projected.second.turnaround.radius.anticipated.by.corrected.intrinsic.luminosity^2

distance_SQR=0
Galaxies_in_Group=matrix(c(1:length(df_2MASS_South$RA),15),ncol=15,nrow=length(df_2MASS_South$RA) + 1,byrow = TRUE)
Galaxies_in_Group=Galaxies_in_Group*0
df_HIPASS$Flux_r=mag2flux(mag=df_HIPASS$Red_Mag_AUTO_Calibrated,zero_pt = 21.10,ABwave = FALSE)


df_2df_2MASS_SouthASS_SouthASS_SouthASS_SouthMASS_SouthMASS_SouthMASS_SouthSouth

max(df_HIPASS$Distance)
max(df_2MASS_South$Distance)
length(which(df_2MASS_South$Distance < max(df_HIPASS$Distance)))
max(sqrt(df_2MASS_South$Radius_SQR))
df_2MASS_South<-subset(df_2MASS_South,df_2MASS_South$Distance <= max(df_HIPASS$Distance) + 
                      max(sqrt(df_2MASS_South$Radius_SQR)))
belong_to_group=c(1:length(df_HIPASS$RA))*0
belong_to_galaxy=c(1:length(df_2MASS_South$RA)) * 0
group_HI_mass=c(1:length(df_2MASS_South$RA))*0
which(is.na.data.frame(df_2MASS_South$Distance))
which(is.na.data.frame(df_HIPASS))

belong_to_group=c(1:length(df_HIPASS$RA))*0
belong_to_galaxy=c(1:length(df_2MASS_South$RA)) * 0
group_HI_mass_2=c(1:length(df_2MASS_South$RA))*0
for(i in c(1:length(df_2MASS_South$Dec))){
    for(j in c(1:length(df_HIPASS$Dec))) {
        if((df_2MASS_South$Distance[i] - df_HIPASS$Distance[j])*(df_2MASS_South$Distance[i] - df_HIPASS$Distance[j]) < 3){
            distance_SQR=(df_2MASS_South$X[i] - df_HIPASS$X[j])*(df_2MASS_South$X[i] - df_HIPASS$X[j]) +
                    (df_2MASS_South$Y[i] - df_HIPASS$Y[j])*(df_2MASS_South$Y[i] - df_HIPASS$Y[j]) + 
                    (df_2MASS_South$Z[i] - df_HIPASS$Z[j])*(df_2MASS_South$Z[i] - df_HIPASS$Z[j])
                belong_to_galaxy[i]=belong_to_galaxy[i] + 1
                belong_to_group[j]=belong_to_group[j] + 1
                group_HI_mass_2[i]=group_HI_mass_2[i] + df_HIPASS$M_HI[j]
          }  
        }
    }
    
}

group_HI_mass=group_HI_mass*0
Galaxies_in_Group=matrix(c(1:length(df_2MASS_South$RA),15),ncol=15,nrow=length(df_2MASS_South$RA) + 1,byrow = TRUE)
Galaxies_in_Group=Galaxies_in_Group*0
for(i in c(1:length(df_2MASS_South$Dec))){
    for(j in c(1:length(df_HIPASS$Dec))) {
        if((df_2MASS_South$Distance[i] - df_HIPASS$Distance[j])*(df_2MASS_South$Distance[i] - df_HIPASS$Distance[j]) < 5.225){
            distance_SQR=(df_2MASS_South$X[i] - df_HIPASS$X[j])*(df_2MASS_South$X[i] - df_HIPASS$X[j]) +
                    (df_2MASS_South$Y[i] - df_HIPASS$Y[j])*(df_2MASS_South$Y[i] - df_HIPASS$Y[j]) + 
                    (df_2MASS_South$Z[i] - df_HIPASS$Z[j])*(df_2MASS_South$Z[i] - df_HIPASS$Z[j])
            if(distance_SQR < df_2MASS_South$Radius_SQR[i]){
                Galaxies_in_Group[i,1] = Galaxies_in_Group[i,1] + 1
                Galaxies_in_Group[i,Galaxies_in_Group[i,1] + 1] = j
                group_HI_mass[i]=group_HI_mass[i] + df_HIPASS$M_HI[j]
                belong_to_group[j]=belong_to_group[j] + 1
          }  
        }
    }
}

j
i
length(df_2MASS_South$Dec)
length(which(belong_to_galaxy > 0))
length(which(belong_to_group > 0 ))
length(which(Galaxies_in_Group[,1] > 0))

df_2MASS_South$Bi.weight.projected.virial.radius.Rij[which(df_2MASS_South$Number.of.group..nest..members.in.2MRS.K.11.75 >4)]

qplot(log10(which(group_HI_mass > 0))) 

i
length(df_2MASS_South$RA)
which(Galaxies_in_Group == 1581)

Galaxies_in_Group[21155,]

library("scatterplot3d")

scatterplot3d(df_HIPASS$X,df_HIPASS$Y,df_HIPASS$Z)

df_1=data.frame(X=df_2MASS_South$X,Y=df_2MASS_South$Y,Z=df_2MASS_South$Z,R=df_2MASS_South$Distance,Survey="2MASS")
df_2=data.frame(X=df_HIPASS$X,Y=df_HIPASS$Y,Z=df_HIPASS$Z,R=df_HIPASS$Distance,Survey="HIPASS")
df_full=rbind(df_1,df_2)

library(car)
library(rgl)

scatterplot3d(x=df_full$X,y=df_full$Y,z=df_full$Z)

qplot(df_full$X,df_full$Y,color=df_full$Survey,alpha=0.5)
qplot(df_full$X,df_full$Z,color=df_full$Survey,alpha=0.5)
qplot(df_full$Y,df_full$Z,color=df_full$Survey,alpha=0.5)

data=read.csv("2mass_groups.csv",skip=4)
df=data.frame(data)

summary(data)

length(which(df$Number.of.group..nest..members.in.2MRS.K.11.75 >2 &
            df_2MASS_South,df_2MASS_South$Distance <= max(df_HIPASS$Distance) + 
                      max(sqrt(df_2MASS_South$Radius_SQR))) )



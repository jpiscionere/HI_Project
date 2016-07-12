
data=read.csv("/Users/jap/HOPCAT2005TabSep.csv")

summary(data)

library('ggplot2')

df=data.frame(data)

length(df$Morphology)

#df[df == "XXXXX"] <- NA

length(df$Morphology)

summary(df)

length(df$Morphology)

which(df$Morphology=="S0")

df$S0 <- ifelse(grepl('S0', df$Morphology), "S0", "Other")

df$S0[which(df$Morphology=="XXXXX")]="Unclassified"

table(df$S0)

df$Bj_Mag=as.numeric(as.character((df$Bj_Mag_AUTO_Calibrated)))

range(na.omit(df$Bj_Mag))

ggplot(df,aes(Sint,color=S0,fill=S0)) + geom_density(alpha=0.1) + xlab("HI Velocity") + theme_bw() + scale_x_log10() 

ggplot(df,aes(x=Bj_Mag,y=Sint,color=S0,fill=S0)) + 
    geom_point(alpha=0.5) + xlab("M_B Mag") + ylab("Integrated HI radio flux") + 
    theme_bw() + scale_y_log10() + xlim(18,11)

df$RA=as.numeric(as.character(df$SExRA))

df$Dec=as.numeric(as.character(df$SExDec))

ggplot(df,aes(RA,Dec)) + geom_point(aes(colour = factor(S0)))

ggplot(df,aes(df$HicatExtl,df$HicatExtb)) + geom_point(aes(colour = factor(S0)))

df_short=data.frame(RA=df$RA,Dec=df$Dec,Morphology=df$Morphology,S0=df$S0,Bj_Mag=df$Bj_Mag,Sint=df$Sint,ID=df$ID)

summary(df_short)

df_short=na.omit(df_short)

summary(df_short)

range(df_short$Dec)

ngroup<-c(1:length(df_short$Morphology))*0

length(df_short$Morphology)
length(df_short$RA)
length(which(df_short$S0=="S0"))



for(i in 1:length(df_short$Morphology)){
        for(j in 1:length(df_short$Morphology)){
            if(df_short$RA[i] - df_short$RA[j] < 2){
                if(df_short$Dec[i] - df_short$Dec[j] < 2){
                    distance2=(df_short$Dec[i]-df_short$Dec[j])^2 + (df_short$RA[i] - df_short$RA[j])^2
                    if(distance2 < 1){
                        ngroup[i]=ngroup[i]+1
                        }
                }
                
            }
        }
}

ngroup

df_short$Ngroup=ngroup -1

ggplot(df_short,aes(Ngroup,color=S0,fill=S0)) + geom_bar(aes(fill = S0),position="dodge") + 
    stat_bin(aes(group=S0, y=..density..)) + 
    xlab("Number of Galaxies Within 5 Degrees") + theme_bw()

df_short$Morphology_Filtered <- ifelse(grepl('S0', df_short$Morphology), "S0", "Other")



df_short$Morphology_Filtered[grepl('pec', df_short$Morphology)=="TRUE"]="Peculiar"

df_short$Morphology_Filtered[c(1:20)]

df_short$Morphology_Filtered[grepl('SB', df_short$Morphology)=="TRUE"]="SB"

df_short$Morphology[c(1:20)]

df_short$Morphology_Filtered[grepl('SA', df_short$Morphology)=="TRUE"]="SA"

df_short$Morphology_Filtered[c(1:20)]

table(df_short$Morphology_Filtered)

df_short$Morphology_Filtered["Sa"]

length(df_short$Morphology_Filtered)

ggplot(df_short,aes(Ngroup,color=Morphology_Filtered,fill=Morphology_Filtered)) + 
    stat_bin(aes(group=Morphology_Filtered, y=..density..),position=position_dodge(width=1),bins=6) + 
    xlab("Number of Galaxies Within 5 Degrees") + theme_bw()

C=299792.458
H0=70

get_HI_mass<-function(velocity,flux){
    #Inputs: D=Distance in Mpc
    #        F=Column density in Jy km s^-1
    D=velocity/H0
    M_HI=2.365*10^5*D^2*F #Msun
    return(M_HI)
}

range(df$Sint)

range(df$vel_mom)

M_HI=2.365*10^5*(df$vel_mom/H0)^2*df$Sint

hist(log10(M_HI))

plot(as.numeric(df$Red_Mag_AUTO_Calibrated),log10(M_HI))

df$HI_Mass=M_HI

ggplot(df,aes(HI_Mass,color=S0,fill=S0)) + geom_density(alpha=0.1) + xlab("MHI/Msun") + theme_bw() + scale_x_log10() 

df_Classified<-subset(df,Morphology!="XXXXX")

length(df_Classified$Morphology)

hist(log10(df_Classified$HI_Mass),breaks=30,xlab=c("Log HI Mass"),main=c("Log HI Mass"))

M_HI_ned=2.365*10^5*(df_Classified$Velocity_Ned/H0)^2*df_Classified$Sint

length(M_HI_ned)
length(df_Classified$HI_Mass)

hist(log10(df_Classified$HI_Mass),col="blue",breaks=30)
hist(log10(na.omit(M_HI_ned)),breaks=30,col="red",add=T)

M_HI_6df=2.365*10^5*(as.numeric(df_Classified$Velocity_6dF)/H0)^2*df_Classified$Sint

hist(log10(M_HI_6df),breaks=30)

ggplot(df_Classified,aes(HI_Mass,color=S0,fill=S0)) + geom_density(alpha=0.1) + xlab("MHI/Msun") + theme_bw() + scale_x_log10() 

ggplot(df_Classified,aes(M_HI_ned,color=S0,fill=S0)) + geom_density(alpha=0.1) + xlab("MHI/Msun") + theme_bw() + scale_x_log10() 

df_Classified$Velocity_LG=df_Classified$vel_mom + 300*sin(df_Classified$HicatExtl) + cos(df_Classified$HicatExtb)

df_Classified$HI_Mass_LG=2.365*10^5*(df_Classified$Velocity_LG/H0)^2*df_Classified$Sint

ggplot(df_Classified,aes(HI_Mass_LG,color=S0,fill=S0)) + geom_density(alpha=0.1) + xlab("MHI/Msun") + theme_bw() + scale_x_log10() 

length(which(df_Classified$S0=="S0" & log10(df_Classified$HI_Mass) < 9 ))/length(which(df_Classified$S0=="S0"))

length(which(df_Classified$S0=="Other" & log10(df_Classified$HI_Mass) < 9 )) / length(which(df_Classified$S0=="Other"))

length(which(df_Classified$S0=="S0" & log10(df_Classified$HI_Mass) > 9 ))/length(which(df_Classified$S0=="S0"))
length(which(df_Classified$S0=="Other" & log10(df_Classified$HI_Mass) > 9 )) / length(which(df_Classified$S0=="Other"))

ggplot(subset(df_Classified,log10(df_Classified$HI_Mass) < 9),aes(HicatExtl,HicatExtb)) + 
geom_point(aes(colour = factor(S0))) + ggtitle("log(HI_Mass) < 9 Galaxies")

df_Classified$hipass_name[which(df_Classified$S0=="S0" & 
                                log10(df_Classified$HI_Mass) < 9)]

df_Classified$Weird<-FALSE

df_Classified$Weird[which(df_Classified$S0=="S0" & 
                                log10(df_Classified$HI_Mass) < 9 )] <- TRUE

length(which(df_Classified$Weird=="TRUE"))

ggplot(subset(df_Classified,df_Classified$S0=="S0"),aes(Sint,vel_mom)) + 
geom_point(aes(colour = factor(Weird))) + ggtitle("Where do the weird galaxies live?") + scale_x_log10()

ggplot(df_Classified,aes(vel_mom/H0,color=Weird,fill=Weird)) + geom_density(alpha=0.1) +xlab("Distance") + theme_bw() + scale_x_log10() 

ggplot(df_Classified,aes(Velocity_LG/H0,color=Weird,fill=Weird)) + geom_density(alpha=0.1) +xlab("Distance") + theme_bw() + scale_x_log10() 

df_Classified$Weird<-FALSE

df_Classified$Weird[which(df_Classified$S0=="S0" & 
                                log10(df_Classified$HI_Mass) < 9 &
                               df_Classified$HicatExtl > 200 &
                               df_Classified$HicatExtl < 300 )] <- TRUE

df_Classified$Galaxy_Name[which(df_Classified$Weird=="TRUE")]

redshifts<-c(0.003815,0.005524,0.003031,0.003552,0.004195,0.005856,0.004995)

names=df_Classified$Galaxy_Name[which(df_Classified$Weird=="TRUE")]

weird_morphs<-c(1:length(names))*0

for(i in c(1:length(names))){
    weird_morphs = df$Morphology[which(df$Galaxy_Name == names[i])] ;
    print(weird_morphs[1])}

velocities=df_Classified$vel_mom[which(df_Classified$Weird=="TRUE")]

length(redshifts)

length(velocities)

C*redshifts/velocities




installed.packages("ANN")

data<- query<- cbind(1:10, 1:10,1:10)
get.knn(data, k=5)

library("FNN")

data<- query<- cbind(1:10, 1:10,1:10)
get.knn(data, k=5)

DEG_TO_RAD=0.01745328888
RAD_TO_DEG=57.2957914331

X_s=sin((90-df_Classified$Dec) * DEG_TO_RAD)*cos(df_Classified$RA * DEG_TO_RAD) 
Y_s=sin((90-df_Classified$Dec) * DEG_TO_RAD)*sin(df_Classified$RA * DEG_TO_RAD) 
Z_s=cos((90-df_Classified$Dec) * DEG_TO_RAD) 

cartesian=data.frame(X=X_s,Y=Y_s,Z=Z_s)

length(X_s)

ngroup=c(1:length(df_Classified$Dec)) * 0
for(i in c(1:length(df_Classified$Dec))){
        for(j in c(1:length(df_Classified$Dec))){
              #  if(abs(df_Classified$Dec[i] - df_Classified$Dec[j]) < 2){
                    cos_Theta=X_s[i] * X_s[j] + Y_s[i] * Y_s[j] + Z_s[i] * Z_s[j]
                    #ngroup[j]=acos(cos_Theta)
                     if(cos_Theta > 0.99984769502 ){
                           ngroup[i]=ngroup[i]+1
                       }
             #   }
                
        }
}

ngroup=ngroup-1

ngroup

hist(ngroup*RAD_TO_DEG)

summary(df_Classified)

table(df_Classified$Morphology_Filtered)

df_Classified$Ngroup=ngroup

ggplot(df_Classified,aes(Ngroup,color=Morphology_Filtered,fill=Morphology_Filtered)) + 
    stat_bin(aes(group=Morphology_Filtered, y=..density..),position=position_dodge(width=1),bins=6) + 
    xlab("Number of Galaxies Within 1 Degree") + theme_bw()

library("FNN")

distances<-matrix(, nrow = length(df_Classified$Dec), ncol = 2)

cosTheta<-c(1:length(df_Classified$Dec))*0

for(i in c(1:length(df_Classified$Dec))){
    for(j in c(2:length(df_Classified$Dec))) {
            cosTheta[j]=X_s[i] * X_s[j] + Y_s[i] * Y_s[j] + Z_s[i] * Z_s[j]
    }
    cosTheta_sorted<-sort(cosTheta, method = "shell", index.return = TRUE,decreasing = TRUE)
    distances[i,1]<-cosTheta_sorted$x[7]
    distances[i,2]<-cosTheta_sorted$ix[7]
}

summary(cosTheta_sorted)

cosTheta_sorted$x[2]

summary(distances)

acos(distances[,1])*RAD_TO_DEG

hist(acos(distances[,1])*RAD_TO_DEG)

df_Classified$SevenNN<-acos(distances[,1])*RAD_TO_DEG

ggplot(df_Classified,aes(SevenNN,color=S0,fill=S0)) + geom_density(alpha=0.1) + xlab("Degrees to 7th Nearest Neighbor") + theme_bw() 

ggplot(df_Classified,aes(SevenNN)) + geom_histogram(bins=50) + xlab("Degrees to 7th Nearest Neighbor") + theme_bw()  



ggplot(df_Classified,aes(SevenNN,..density..,fill=S0)) + geom_histogram( ) + xlab("Degrees to 7th Nearest Neighbor") + theme_bw() + xlim(0,10) 

summary(df_Classified$SevenNN)

df_Classified$RA[which(df_Classified$SevenNN > 10 & df_Classified$S0=='S0')]

df_Classified$Distances<-df_Classified$vel_mom/H0

summary(df_Classified$Distances)

xyz_distances<-df_Classified$Distances*0

rp_distances<-df_Classified$Distances*0

rp_distances_id<-df_Classified$Distances*0

i=0
for(i in c(1:length(df_Classified$Dec))){
    for(j in c(1:length(df_Classified$Dec))) {
            xyz_distances[j]=(df_Classified$Distances[i]*X_s[i] - df_Classified$Distances[j]*X_s[j])^2 + 
                    (df_Classified$Distances[i]*Y_s[i] - df_Classified$Distances[j]*Y_s[j])^2 + 
                    (df_Classified$Distances[i]*Z_s[i] - df_Classified$Distances[j]*Z_s[j])^2            
    }
    rp_sqr_sorted<-sort(xyz_distances, method = "shell", index.return = TRUE)
    rp_distances[i]<-rp_sqr_sorted$x[7]
    rp_distances_id[i]<-rp_sqr_sorted$ix[7]
}

df_Classified$SevenNN_MPC=sqrt(rp_distances)

summary(df_Classified$SevenNN_MPC)

hist(sqrt(rp_sqr_sorted$x))

ggplot(df_Classified,aes(SevenNN_MPC)) + geom_histogram(bins=50) + xlab("Mpc/h to 7th Nearest Neighbor") + theme_bw()  

summary(df_Classified$SevenNN_MPC)

PI*RAD_TO_DEG

length(df$hipass_name)/180^2

7/(length(df$hipass_name)/180^2)

(7/(length(df$hipass_name)/180^2))^(1/3)

ggplot(df_Classified,aes(SevenNN)) + geom_histogram(bins=50,fill="magenta") + xlab("Degrees to 7th Nearest Neighbor") + 
theme_bw()  + geom_vline(xintercept=3.74588192201299,lwd=1.5,lty=2 )  + annotate("text", x = 7, y = 600, label = "Expected",size=8)

ggsave("hipass_degrees_histo.png")

ggplot(df_Classified,aes(x=df_Classified$Redshift,y=df_Classified$Bj_Mag)) + geom_point()

qplot(df_Classified$Redshift)

Volume_of_Survey=4/3*PI*(C*0.004/H0)^3 * 0.5

Volume_of_Survey

length(df_Classified$Dec)/Volume_of_Survey

7/(length(df_Classified$Dec)/Volume_of_Survey)

(7/(length(df_Classified$Dec)/Volume_of_Survey))^(1/3)

ggplot(df_Classified,aes(SevenNN_MPC)) + geom_histogram(bins=50,fill="dodgerblue") + xlab("Mpc/h to 7th Nearest Neighbor") + 
theme_bw() + geom_vline(xintercept=2.84464294973747,lwd=1.5,lty=2) + annotate("text", x = 15, y = 750, label = "Expected",size=8)

plot=ggplot(df_Classified,aes(SevenNN_MPC)) + geom_histogram(bins=50) + xlab("Mpc/h to 7th Nearest Neighbor") + 
theme_bw() + geom_vline(xintercept=2.84464294973747 )

ggsave("histogram_HIPASS.png")

ggplot(df_Classified,aes(df_Classified$HicatExtl,df_Classified$HicatExtb)) + geom_point()

ggplot(df_Classified,aes(df_Classified$RA,df_Classified$Dec)) + geom_point()

library(grid)


r=1

data_groups=read.csv("Downloads/EDDtable25May2016010159.txt")

summary(data_groups)

qplot(data_groups$Galactic.latitude,data_groups$Galactic.longitude)

library("astrolibR")

coords=glactc(gl = data_groups$Galactic.longitude,gb = data_groups$Galactic.latitude,j=2,
              year = 2000,ra=data_groups$RA,dec=data_groups$Dec,degree=TRUE)

summary(coords)

data_groups$RA=coords$ra

data_groups$Dec=coords$dec

qplot(data_groups$RA,data_groups$Dec)

ncol(data_groups)

data_groups_full=read.csv("2MASS_group_catalogue.txt")

coords=glactc(gl = data_groups_full$Galactic.longitude,gb = data_groups_full$Galactic.latitude,j=2,year = 2000,degree = TRUE)

length(coords$ra)

data_groups_full$RA=coords$ra
data_groups_full$Dec=coords$dec

qplot(data_groups_full$RA,data_groups_full$Dec)

df_2MASS_South<-subset(data_groups_full,Dec < 0)

df_2MASS_South$Distance=df_2MASS_South$Velocity.in.CMB.frame..cosmological.adjustment.described.in.CF2.data.paper/H0

df_2MASS_South$X=df_2MASS_South$Distance * sin((90-df_2MASS_South$Dec) * DEG_TO_RAD)*cos(df_2MASS_South$RA * DEG_TO_RAD) 
df_2MASS_South$Y=df_2MASS_South$Distance * sin((90-df_2MASS_South$Dec) * DEG_TO_RAD)*sin(df_2MASS_South$RA * DEG_TO_RAD) 
df_2MASS_South$Z=df_2MASS_South$Distance * cos((90-df_2MASS_South$Dec) * DEG_TO_RAD)

df_HIPASS=na.omit(df)

df_HIPASS

df_HIPASS$Distance=df_HIPASS$vel_mom/H0

df_HIPASS$dec_str[2]

df_HIPASS$X=df_HIPASS$Distance * sin((90-df_HIPASS$Dec) * DEG_TO_RAD)*cos(df_HIPASS$RA * DEG_TO_RAD) 
df_HIPASS$Y=df_HIPASS$Distance * sin((90-df_HIPASS$Dec) * DEG_TO_RAD)*sin(df_HIPASS$RA * DEG_TO_RAD) 
df_HIPASS$Z=df_HIPASS$Distance * cos((90-df_HIPASS$Dec) * DEG_TO_RAD)

question=is.na.data.frame(df_HIPASS)

summary(question)

belong_to_group=c(1:length(df_HIPASS$RA))*0
belong_to_galaxy=c(1:length(df_2MASS_South$RA)) * 0

df_2MASS_South$Radius_SQR=df_2MASS_South$Projected.second.turnaround.radius.anticipated.by.corrected.intrinsic.luminosity^2

distance_SQR=0

#for(i in c(1:length(df_2MASS_South$Dec))){
for(i in c(1:2) ){
    for(j in c(1:length(df_HIPASS$Dec))) {
            distance_SQR=(df_2MASS_South$X[i] - df_HIPASS$X[j])^2 +
                    (df_2MASS_South$Y[i] - df_HIPASS$Y[j])^2 + 
                    (df_2MASS_South$Z[i] - df_HIPASS$Z[j])^2
          if(distance_SQR < df_2MASS_South$Radius_SQR[i]){
              belong_to_galaxy[i]=belong_to_galaxy[i] + 1
              belong_to_group[j]=belong_to_group[j] + 1
          }  
          
    }
}

distance_SQR



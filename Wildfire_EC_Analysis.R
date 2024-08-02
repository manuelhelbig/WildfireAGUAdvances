################################################################################
##### Download & analyse AmeriFlux data from wildfire disturbance sites ########
################################################################################
# written by: Manuel Helbig (helbig@gfz-potsdam.de)
# last update: August 2, 2024
# code used for the publication:
# Helbig et al. (2024): Boreal Forest Fire Causes Daytime Surface Warming During
# Summer to Exceed Surface Cooling During Winter in North America. AGU Advances.
################################################################################
rm(list=ls())
library(amerifluxr)
library(robustbase)
library(bigleaf)
library(devtools)
library(rapportools)
library(forecast)
library(ggplot2)
library(openair)
library(purrr)
library(tidyverse)
library(patchwork)
library(LakeMetabolizer)
#install_github("laubblatt/cleaRskyQuantileRegression",force = TRUE)
library(cleaRskyQuantileRegression)
################################################################################
setwd(".../DATA_REPOSITORY") # set working directory
HEIGHTS <- read.csv("HeightsWildfire.csv")
################################################################################
# download AmeriFlux site data from AmeriFlux database
################################################################################
# read flux dataset (US-RPF)
base <- amf_read_base(file = ".../AMF_US-RPF_BASE-BADM_6-5.zip",
                      unzip = TRUE,
                      parse_timestamp = TRUE)
# latitude and longitude of site
lat = 65.1198
lon = -147.4290
# extract time vector
base$Time <- format(as.POSIXct(base$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "US/Alaska"), format = "%H:%M")
base$Date <- format(as.Date(base$TIMESTAMP,"%Y-%m-%d",tz = "US/Alaska"), format = "%Y-%m-%d")
base$Time <- as.ITime(base$Time, tz = "US/Alaska")
base$Date <- as.IDate(base$Date, tz = "US/Alaska")

# read  radiation fluxes
base$SW_IN=base$SW_IN_PI_F_1_1_1
base$SW_IN[which(is.na(base$SW_IN))]=base$SW_IN_PI_F_1_1_2[which(is.na(base$SW_IN))]
base$SW_OUT=base$SW_OUT_PI_F_1_1_1
base$SW_OUT[which(is.na(base$SW_OUT))]=base$SW_OUT_PI_F_1_1_2[which(is.na(base$SW_OUT))]
base$RN=base$NETRAD_1_1_1
base$RN[which(is.na(base$RN))]=base$NETRAD_1_1_2[which(is.na(base$RN))]

# get latent and sensible heat flux
base$LE=base$LE_PI_F
base$H=base$H_PI_F
base$EF = base$LE_PI_F/(base$LE_PI_F+base$H_PI_F) # evaporative fraction
base$EF[base$EF<0 | base$EF>1]=NA

# air temperature and humidity
base$TA=base$TA_1_1_1
base$TA[which(is.na(base$TA))]=base$TA_1_2_1[which(is.na(base$TA))]
base$RH=base$RH_1_1_1
base$RH[which(is.na(base$RH))]=base$RH_1_2_1[which(is.na(base$RH))]
base$VPD=rH.to.VPD(base$RH/100,base$TA,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),constant=bigleaf.constants())

# ground heat flux
base$G=(base$G_1_1_1+base$G_1_1_2)/2

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='US-RPF1']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='US-RPF1']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='US-RPF1']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*base$H)/(1.293*1004*(base$TA+273.15)*base$USTAR^3)

ind=which(base$YEAR==2013 & base$MONTH == 7 & base$DAY == 1 & base$HOUR == 0 & base$MINUTE == 15)
# account for growth of vegetation
d0 = HEIGHTS$hc[HEIGHTS$SITE=='US-RPF2']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='US-RPF2']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='US-RPF2']
zeta[ind:length(zeta)] = (-0.4*(z-d0)*9.81*base$H[ind:length(zeta)] )/(1.293*1004*(base$TA[ind:length(zeta)] +273.15)*base$USTAR[ind:length(zeta)]^3)

# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA

# get aerodynamic conductance
df <- data.frame(Tair=base$TA,pressure=base$PA_PI_F,wind=base$WS_1_1_1,ustar=base$USTAR,H=base$H)
Ga = aerodynamic.conductance(df,Rb_model="Thom_1972")

# get aerodynamic temperature
surf = surface.conditions(Tair=base$TA,pressure=base$PA_PI_F,LE=base$LE,H=base$H,VPD=base$VPD,Ga=Ga$Ga_h)
base$TSURF=surf$Tsurf
ind = tsoutliers(base$TSURF)
base$TSURF[ind$index]=NA

# get 10 m air temperature (see Novick & Katul 2020: eq 4)
base$T10 = base$TSURF-(base$H/(0.4*1.293*1004*base$USTAR))*(log(10/z0h)-theta)
d0 = HEIGHTS$hc[HEIGHTS$SITE=='US-RPF1']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='US-RPF1']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='US-RPF1']
ind=which(base$YEAR==2013 & base$MONTH == 7 & base$DAY == 1 & base$HOUR == 0 & base$MINUTE == 15)
base$T10[1:ind-1] = base$TSURF[1:ind-1]-(base$H[1:ind-1]/(0.4*1.293*1004*base$USTAR[1:ind-1]))*(log(10/z0h)-theta[1:ind-1])

uYEAR = unique(base$YEAR)
s=0
# create daily time series
RPF <- data.frame(matrix(ncol = 18, nrow = length(uYEAR)*365))
colnames(RPF) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEAR)-1
for (y in 1:length(uYEAR))
{
  for (n in 1:365)
  {
    s=s+1
    RPF$DOY[s]=n
    RPF$YEAR[s]=sYEAR+y
    # select afternoon hours (between 12pm and 4pm)
    ind = which(base$DOY == n & base$YEAR == sYEAR+y & base$HOUR >= 12 & base$HOUR <= 16)
    # select all daytime hours (between 6am and 6pm)
    ind1 = which(base$DOY == n & base$YEAR == sYEAR+y & base$HOUR >= 6 & base$HOUR <= 18);
    # select hours before noon (between 6am and 12pm)
    ind2 = which(base$DOY == n & base$YEAR == sYEAR+y & base$HOUR >= 6 & base$HOUR <= 12);
    # select all measurements
    all = which(base$DOY == n & base$YEAR == sYEAR+y)
    
    mnt = base$MONTH[ind[1]]
    year = base$YEAR[ind[1]]
    
    RPF$BR[s]=sum(base$H[ind1])/sum(base$LE[ind1]) # Bowen ratio (daytime)
    RPF$Ga[s]=median(Ga$Ga_h[ind], na.rm=T) # aerodynamic conductance (afternoon)
    RPF$LE[s]=median(base$LE[ind1], na.rm=T) # latent heat flux (daytime)
    RPF$VPD[s]=median(base$VPD[ind], na.rm=T) # vapor pressure deficit (afternoon)
    RPF$WS[s]=median(base$WS_1_1_1[ind], na.rm=T) #wind speed (afternoon)
    
    # maximum air temperature in the afternoon
    RPF$TA[s]=max(base$TA[ind], na.rm=T) # (afternoon)
    # mean aerodynamic temperature in the afternoon
    RPF$TS[s]=mean(base$TSURF[ind], na.rm=T) # (afternoon)
    # temperature gradient in the afternoon
    RPF$dT[s]=mean(base$TSURF[ind]-base$TA[ind], na.rm=T) # (afternoon)
    RPF$dT10[s]=mean(base$TSURF[ind]-base$T10[ind], na.rm=T) # (afternoon)
    
    # afternoon albedo
    if (is.na(sum(base$SW_OUT_1_1_1[ind])/sum(base$SW_IN_1_1_1[ind]))) {
      RPF$ALB[s]=sum(base$SW_OUT_1_1_2[ind])/sum(base$SW_IN_1_1_2[ind])
    }
    else 
    {
      RPF$ALB[s]=sum(base$SW_OUT_1_1_1[ind])/sum(base$SW_IN_1_1_1[ind])
    }
    
    # afternoon evaporative fraction
    RPF$EF[s]=sum(base$LE[ind])/sum(base$LE[ind]+base$H[ind])
  }
} 

############## filter time series for absolute limits and outliers #############
RPF$BR[abs(RPF$BR)>6]=NA 
RPF$ALB[abs(RPF$ALB)>1 | RPF$ALB<0.01]=NA
RPF$EF[RPF$EF<0]=NA
RPF$EF[RPF$EF>1]=NA
ind = tsoutliers(RPF$dT)
RPF$dT[ind$index]=NA
ind = tsoutliers(RPF$dT10)
RPF$dT10[ind$index]=NA
ind = tsoutliers(RPF$TS)
RPF$TS[ind$index]=NA

#####################################
### read US-Fcr ###
#####################################
# read flux dataset
baseFcr <- amf_read_base(file = ".../AMF_US-Fcr_BASE-BADM_2-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude of site
lat = 65.3968
lon = -148.9348
# extract time vector
baseFcr$Time <- format(as.POSIXct(baseFcr$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "US/Alaska"), format = "%H:%M")
baseFcr$Date <- format(as.Date(baseFcr$TIMESTAMP,"%Y-%m-%d",tz = "US/Alaska"), format = "%Y-%m-%d")
baseFcr$Time <- as.ITime(baseFcr$Time, tz = "US/Alaska")
baseFcr$Date <- as.IDate(baseFcr$Date, tz = "US/Alaska")

# read  radiation fluxes
baseFcr$SW_IN=baseFcr$SW_IN_PI_F_1_1_1
baseFcr$SW_OUT=baseFcr$SW_OUT_PI_F_1_1_1
baseFcr$RN=baseFcr$NETRAD_PI_F_1_1_1
baseFcr$LST=baseFcr$T_CANOPY_1_1_1

# get latent and sensible heat flux
baseFcr$EF = baseFcr$LE_PI_F/(baseFcr$LE_PI_F+baseFcr$H_PI_F)
baseFcr$EF[baseFcr$EF<0 | baseFcr$EF>1]=NA 

# air temperature and humidity
baseFcr$TA=baseFcr$TA_1_1_1
baseFcr$RH=baseFcr$RH_1_1_1
baseFcr$VPD=rH.to.VPD(baseFcr$RH/100,baseFcr$TA,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),constant=bigleaf.constants())

# ground heat flux
baseFcr$G=(baseFcr$G_1_1_1+baseFcr$G_1_1_2)/2

# get aerodynamic conductance
dfFcr <- data.frame(Tair=baseFcr$TA,pressure=baseFcr$PA,wind=baseFcr$WS_1_1_1,ustar=baseFcr$USTAR,H=baseFcr$H_PI_F)
GaFcr = aerodynamic.conductance(dfFcr,Rb_model="Thom_1972")
# get aerodynamic temperature
surfFcr = surface.conditions(Tair=baseFcr$TA,pressure=baseFcr$PA,LE=baseFcr$LE_PI_F,H=baseFcr$H_PI_F,VPD=baseFcr$VPD,Ga=GaFcr$Ga_h)
baseFcr$TSURF=surfFcr$Tsurf
ind = tsoutliers(baseFcr$TSURF)
baseFcr$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='US-FCR']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='US-FCR']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='US-FCR']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseFcr$H )/(1.293*1004*(baseFcr$TA +273.15)*baseFcr$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseFcr$T10 = baseFcr$TSURF-(baseFcr$H/(0.4*1.293*1004*baseFcr$USTAR))*(log(10/z0h)-theta)

uYEARFcr = unique(baseFcr$YEAR)
s=0
# create daily time series
FCR <- data.frame(matrix(ncol = 18, nrow = length(uYEARFcr)*365))
colnames(FCR) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARFcr)-1

for (y in 1:length(uYEARFcr))
{
  for (n in 1:365)
  {
    s=s+1
    FCR$DOY[s]=n
    FCR$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseFcr$DOY == n & baseFcr$YEAR == sYEAR+y & baseFcr$HOUR >= 12 & baseFcr$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseFcr$DOY == n & baseFcr$YEAR == sYEAR+y & baseFcr$HOUR >= 6 & baseFcr$HOUR <= 18);
    ind2 = which(baseFcr$DOY == n & baseFcr$YEAR == sYEAR+y & baseFcr$HOUR >= 6 & baseFcr$HOUR <= 12);
    all = which(baseFcr$DOY == n & baseFcr$YEAR == sYEAR+y)
    mnt = baseFcr$MONTH[ind[1]]
    year = baseFcr$YEAR[ind[1]]
    
    FCR$BR[s]=sum(baseFcr$H[ind1])/sum(baseFcr$LE[ind1]) # Bowen ratio
    FCR$Ga[s]=median(GaFcr$Ga_h[ind], na.rm=T) # aerodynamic conductance
    FCR$LE[s]=median(baseFcr$LE[ind1], na.rm=T)
    FCR$VPD[s]=median(baseFcr$VPD[ind], na.rm=T)
    FCR$WS[s]=median(baseFcr$WS_1_1_1[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    FCR$TA[s]=max(baseFcr$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    FCR$TS[s]=mean(baseFcr$TSURF[ind], na.rm=T)
    # temperature gradient
    FCR$dT[s]=mean(baseFcr$TSURF[ind]-baseFcr$TA[ind], na.rm=T)
    FCR$dT10[s]=mean(baseFcr$TSURF[ind]-baseFcr$T10[ind], na.rm=T)

    # afternoon albedo
    FCR$ALB[s]=sum(baseFcr$SW_OUT_1_1_1[ind])/sum(baseFcr$SW_IN_1_1_1[ind])
    # afternoon evaporative fraction
    FCR$EF[s]=sum(baseFcr$LE[ind])/sum(baseFcr$LE[ind]+baseFcr$H[ind])
  }
}

############## filter time series ##############
FCR$BR[abs(FCR$BR)>6]=NA
FCR$ALB[abs(FCR$ALB)>1 | FCR$ALB<0.01]=NA
FCR$EF[FCR$EF<0]=NA
FCR$EF[FCR$EF>1]=NA
ind = tsoutliers(FCR$dT)
FCR$dT[ind$index]=NA
ind = tsoutliers(FCR$dT10)
FCR$dT10[ind$index]=NA
ind = tsoutliers(FCR$TS)
FCR$TS[ind$index]=NA
#####################################
### read US-Bn3 ###
#####################################
# read flux dataset
baseBn3 <- amf_read_base(file = ".../AMF_US-Bn3_BASE-BADM_1-1.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude of site
lat = 63.914
lon = -145.744
# extract time vector
baseBn3$Time <- format(as.POSIXct(baseBn3$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "US/Alaska"), format = "%H:%M")
baseBn3$Date <- format(as.Date(baseBn3$TIMESTAMP,"%Y-%m-%d",tz = "US/Alaska"), format = "%Y-%m-%d")
baseBn3$Time <- as.ITime(baseBn3$Time, tz = "US/Alaska")
baseBn3$Date <- as.IDate(baseBn3$Date, tz = "US/Alaska")

# displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='US-BN3']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='US-BN3']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='US-BN3']

# adjust wind speed to match height of eddy covariance flux measurements
baseBn3$WS = baseBn3$WS*(log((7.8-d0)/(z0))/log((9-d0)/(z0)))

# get latent and sensible heat flux
baseBn3$EF = baseBn3$LE/(baseBn3$LE+baseBn3$H)
baseBn3$EF[baseBn3$EF<0 | baseBn3$EF>1]=NA 

# read radiation fluxes
baseBn3$RN=baseBn3$NETRAD

# air temperature and humidity
baseBn3$RH[baseBn3$RH>100]=100
baseBn3$VPD=rH.to.VPD(baseBn3$RH/100,baseBn3$TA,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),constant=bigleaf.constants())

# get aerodynamic conductance
dfBn3 <- data.frame(Tair=baseBn3$TA,pressure=baseBn3$PA,wind=baseBn3$WS,ustar=baseBn3$USTAR,H=baseBn3$H)
GaBn3 = aerodynamic.conductance(dfBn3,Rb_model="Thom_1972")
# get aerodynamic temperature
surfBn3 = surface.conditions(Tair=baseBn3$TA,pressure=baseBn3$PA,LE=baseBn3$LE,H=baseBn3$H,VPD=baseBn3$VPD,Ga=GaBn3$Ga_h)
baseBn3$TSURF=surfBn3$Tsurf
ind = tsoutliers(baseBn3$TSURF)
baseBn3$TSURF[ind$index]=NA

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseBn3$H )/(1.293*1004*(baseBn3$TA +273.15)*baseBn3$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA

# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseBn3$T10 = baseBn3$TSURF-(baseBn3$H/(0.4*1.293*1004*baseBn3$USTAR))*(log(10/z0h)-theta)

uYEARBn3 = unique(baseBn3$YEAR)
s=0
# create daily time series
BN3 <- data.frame(matrix(ncol = 18, nrow = length(uYEARBn3)*365))
colnames(BN3) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")

sYEAR = min(uYEARBn3)-1

for (y in 1:length(uYEARBn3))
{
  for (n in 1:365)
  {
    s=s+1
    BN3$DOY[s]=n
    BN3$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseBn3$DOY == n & baseBn3$YEAR == sYEAR+y & baseBn3$HOUR >= 12 & baseBn3$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseBn3$DOY == n & baseBn3$YEAR == sYEAR+y & baseBn3$HOUR >= 6 & baseBn3$HOUR <= 18);
    ind2 = which(baseBn3$DOY == n & baseBn3$YEAR == sYEAR+y & baseBn3$HOUR >= 6 & baseBn3$HOUR <= 12);
    all = which(baseBn3$DOY == n & baseBn3$YEAR == sYEAR+y)
    mnt = baseBn3$MONTH[ind[1]]
    year = baseBn3$YEAR[ind[1]]
    
    BN3$BR[s]=sum(baseBn3$H[ind1])/sum(baseBn3$LE[ind1]) # Bowen ratio
    BN3$Ga[s]=median(GaBn3$Ga_h[ind], na.rm=T) # aerodynamic conductance
    BN3$LE[s]=median(baseBn3$LE[ind1], na.rm=T)
    BN3$VPD[s]=median(baseBn3$VPD[ind], na.rm=T)
    BN3$WS[s]=median(baseBn3$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    BN3$TA[s]=max(baseBn3$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    BN3$TS[s]=mean(baseBn3$TSURF[ind], na.rm=T)
    # temperature gradient
    BN3$dT[s]=mean(baseBn3$TSURF[ind]-baseBn3$TA[ind], na.rm=T)
    BN3$dT10[s]=mean(baseBn3$TSURF[ind]-baseBn3$T10[ind], na.rm=T)
    
    # albedo
    BN3$ALB[s]=sum(baseBn3$SW_OUT[ind])/sum(baseBn3$SW_IN[ind])
    # afternoon evaporative fraction
    BN3$EF[s]=sum(baseBn3$LE[ind])/sum(baseBn3$LE[ind]+baseBn3$H[ind])
  }
}

BN3$BR[abs(BN3$BR)>6]=NA
BN3$ALB[abs(BN3$ALB)>1 | BN3$ALB<0.05]=NA
BN3$EF[BN3$EF<0]=NA
BN3$EF[BN3$EF>1]=NA
ind = tsoutliers(BN3$dT)
BN3$dT[ind$index]=NA
ind = tsoutliers(BN3$dT10)
BN3$dT10[ind$index]=NA
ind = tsoutliers(BN3$TS)
BN3$TS[ind$index]=NA
#####################################
### read US-Bn2 ###
#####################################
# read flux dataset
baseBn2 <- amf_read_base(file = ".../AMF_US-Bn2_BASE-BADM_1-1.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude of site
lat = 63.9198
lon = -145.3782
# extract time vector
baseBn2$Time <- format(as.POSIXct(baseBn2$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "US/Alaska"), format = "%H:%M")
baseBn2$Date <- format(as.Date(baseBn2$TIMESTAMP,"%Y-%m-%d",tz = "US/Alaska"), format = "%Y-%m-%d")
baseBn2$Time <- as.ITime(baseBn2$Time, tz = "US/Alaska")
baseBn2$Date <- as.IDate(baseBn2$Date, tz = "US/Alaska")

# displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='US-BN2']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='US-BN2']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='US-BN2']

# adjust wind speed to match height of eddy covariance flux measurements
baseBn2$WS = baseBn2$WS*(log((10-d0)/(z0))/log((12-d0)/(z0)))

# get latent and sensible heat flux
baseBn2$EF = baseBn2$LE/(baseBn2$LE+baseBn2$H)
baseBn2$EF[baseBn2$EF<0 | baseBn2$EF>1]=NA 

# read  radiation fluxes
baseBn2$RN=baseBn2$NETRAD

# air temperature and humidity
baseBn2$RH[baseBn2$RH>100]=100
baseBn2$VPD=rH.to.VPD(baseBn2$RH/100,baseBn2$TA,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),constant=bigleaf.constants())

# get aerodynamic conductance
dfBn2 <- data.frame(Tair=baseBn2$TA,pressure=baseBn2$PA,wind=baseBn2$WS,ustar=baseBn2$USTAR,H=baseBn2$H)
GaBn2 = aerodynamic.conductance(dfBn2,Rb_model="Thom_1972")
# get aerodynamic temperature
surfBn2 = surface.conditions(Tair=baseBn2$TA,pressure=baseBn2$PA,LE=baseBn2$LE,H=baseBn2$H,VPD=baseBn2$VPD,Ga=GaBn2$Ga_h)
baseBn2$TSURF=surfBn2$Tsurf
ind = tsoutliers(baseBn2$TSURF)
baseBn2$TSURF[ind$index]=NA

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseBn2$H )/(1.293*1004*(baseBn2$TA +273.15)*baseBn2$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseBn2$T10 = baseBn2$TSURF-(baseBn2$H/(0.4*1.293*1004*baseBn2$USTAR))*(log(10/z0h)-theta)

uYEARBn2 = unique(baseBn2$YEAR)
s=0
# create daily time series
BN2 <- data.frame(matrix(ncol = 18, nrow = length(uYEARBn2)*365))
colnames(BN2) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARBn2)-1

for (y in 1:length(uYEARBn2))
{
  for (n in 1:365)
  {
    s=s+1
    BN2$DOY[s]=n
    BN2$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseBn2$DOY == n & baseBn2$YEAR == sYEAR+y & baseBn2$HOUR >= 12 & baseBn2$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseBn2$DOY == n & baseBn2$YEAR == sYEAR+y & baseBn2$HOUR >= 6 & baseBn2$HOUR <= 18);
    ind2 = which(baseBn2$DOY == n & baseBn2$YEAR == sYEAR+y & baseBn2$HOUR >= 6 & baseBn2$HOUR <= 12);
    all = which(baseBn2$DOY == n & baseBn2$YEAR == sYEAR+y)
    mnt = baseBn2$MONTH[ind[1]]
    year = baseBn2$YEAR[ind[1]]
    
    BN2$BR[s]=sum(baseBn2$H[ind1])/sum(baseBn2$LE[ind1]) # Bowen ratio
    BN2$Ga[s]=median(GaBn2$Ga_h[ind], na.rm=T) # aerodynamic conductance
    BN2$LE[s]=median(baseBn2$LE[ind1], na.rm=T)
    BN2$VPD[s]=median(baseBn2$VPD[ind], na.rm=T)
    BN2$WS[s]=median(baseBn2$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    BN2$TA[s]=max(baseBn2$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    BN2$TS[s]=mean(baseBn2$TSURF[ind], na.rm=T)
    #BN2$LST[s]=mean(baseBn2$LST[ind], na.rm=T)
    # temperature gradient
    BN2$dT[s]=mean(baseBn2$TSURF[ind]-baseBn2$TA[ind], na.rm=T)
    BN2$dT10[s]=mean(baseBn2$TSURF[ind]-baseBn2$T10[ind], na.rm=T)
    
    # albedo
    BN2$ALB[s]=sum(baseBn2$SW_OUT[ind])/sum(baseBn2$SW_IN[ind])
    # afternoon evaporative fraction
    BN2$EF[s]=sum(baseBn2$LE[ind])/sum(baseBn2$LE[ind]+baseBn2$H[ind])
  }
}

BN2$ALB=NA
BN2$BR[abs(BN2$BR)>6]=NA
BN2$ALB[abs(BN2$ALB)>1 | BN2$ALB<0.05]=NA
BN2$EF[BN2$EF<0]=NA
BN2$EF[BN2$EF>1]=NA
ind = tsoutliers(BN2$dT)
BN2$dT[ind$index]=NA
ind = tsoutliers(BN2$dT10)
BN2$dT10[ind$index]=NA
ind = tsoutliers(BN2$TS)
BN2$TS[ind$index]=NA
#####################################
### read CA-NS5 ###
#####################################
# read flux dataset
baseNS5 <- amf_read_base(file = ".../AMF_CA-NS5_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude of site
lat = 55.8631
lon = -98.4850
# extract time vector
baseNS5$Time <- format(as.POSIXct(baseNS5$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseNS5$Date <- format(as.Date(baseNS5$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseNS5$Time <- as.ITime(baseNS5$Time, tz = "Canada/Central")
baseNS5$Date <- as.IDate(baseNS5$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseNS5$EF = baseNS5$LE/(baseNS5$LE+baseNS5$H)
baseNS5$EF[baseNS5$EF<0 | baseNS5$EF>1]=NA 

# read  radiation fluxes
baseNS5$RN=baseNS5$NETRAD
baseNS5$PA=101.3

# air temperature and humidity
baseNS5$TA=baseNS5$TA_1_1_1
baseNS5$RH[baseNS5$RH>100]=100
baseNS5$VPD=baseNS5$VPD_PI/10

# get aerodynamic conductance
dfNS5 <- data.frame(Tair=baseNS5$TA,pressure=baseNS5$PA,wind=baseNS5$WS,ustar=baseNS5$USTAR,H=baseNS5$H)
GaNS5 = aerodynamic.conductance(dfNS5,Rb_model="Thom_1972")

# get aerodynamic temperature
surfNS5 = surface.conditions(Tair=baseNS5$TA,pressure=baseNS5$PA,LE=baseNS5$LE,H=baseNS5$H,VPD=baseNS5$VPD,Ga=GaNS5$Ga_h)
baseNS5$TSURF=surfNS5$Tsurf
ind = tsoutliers(baseNS5$TSURF)
baseNS5$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS5']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS5']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-NS5']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseNS5$H )/(1.293*1004*(baseNS5$TA +273.15)*baseNS5$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseNS5$T10 = baseNS5$TSURF-(baseNS5$H/(0.4*1.293*1004*baseNS5$USTAR))*(log(10/z0h)-theta)

uYEARNS5 = unique(baseNS5$YEAR)
s=0
NS5 <- data.frame(matrix(ncol = 18, nrow = length(uYEARNS5)*365))
colnames(NS5) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARNS5)-1

for (y in 1:length(uYEARNS5))
{
  for (n in 1:365)
  {
    s=s+1
    NS5$DOY[s]=n
    NS5$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseNS5$DOY == n & baseNS5$YEAR == sYEAR+y & baseNS5$HOUR >= 12 & baseNS5$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseNS5$DOY == n & baseNS5$YEAR == sYEAR+y & baseNS5$HOUR >= 6 & baseNS5$HOUR <= 18);
    ind2 = which(baseNS5$DOY == n & baseNS5$YEAR == sYEAR+y & baseNS5$HOUR >= 6 & baseNS5$HOUR <= 12);
    all = which(baseNS5$DOY == n & baseNS5$YEAR == sYEAR+y)
    mnt = baseNS5$MONTH[ind[1]]
    year = baseNS5$YEAR[ind[1]]
    
    NS5$BR[s]=sum(baseNS5$H[ind1])/sum(baseNS5$LE[ind1]) # Bowen ratio
    NS5$Ga[s]=median(GaNS5$Ga_h[ind], na.rm=T) # aerodynamic conductance
    NS5$LE[s]=median(baseNS5$LE[ind1], na.rm=T)
    NS5$VPD[s]=median(baseNS5$VPD[ind], na.rm=T)
    NS5$WS[s]=median(baseNS5$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    NS5$TA[s]=max(baseNS5$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    NS5$TS[s]=mean(baseNS5$TSURF[ind], na.rm=T)
    #NS5$LST[s]=mean(baseNS5$LST[ind], na.rm=T)
    # temperature gradient
    NS5$dT[s]=mean(baseNS5$TSURF[ind]-baseNS5$TA[ind], na.rm=T)
    NS5$dT10[s]=mean(baseNS5$TSURF[ind]-baseNS5$T10[ind], na.rm=T)
    
    # albedo
    NS5$ALB[s]=sum(baseNS5$SW_OUT[ind])/sum(baseNS5$SW_IN[ind])
    # afternoon evaporative fraction
    NS5$EF[s]=sum(baseNS5$LE[ind])/sum(baseNS5$LE[ind]+baseNS5$H[ind])
  }
}

NS5$BR[abs(NS5$BR)>6]=NA
NS5$ALB[abs(NS5$ALB)>1 | NS5$ALB<0.05]=NA
NS5$EF[NS5$EF<0]=NA
NS5$EF[NS5$EF>1]=NA
ind = tsoutliers(NS5$dT)
NS5$dT[ind$index]=NA
ind = tsoutliers(NS5$dT10)
NS5$dT10[ind$index]=NA
ind = tsoutliers(NS5$TS)
NS5$TS[ind$index]=NA
#####################################
### read CA-NS1 ###
#####################################
# read flux dataset
baseNS1 <- amf_read_base(file = ".../AMF_CA-NS1_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude of site
lat = 55.8792
lon = -98.4839
# extract time vector
baseNS1$Time <- format(as.POSIXct(baseNS1$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseNS1$Date <- format(as.Date(baseNS1$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseNS1$Time <- as.ITime(baseNS1$Time, tz = "Canada/Central")
baseNS1$Date <- as.IDate(baseNS1$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseNS1$EF = baseNS1$LE/(baseNS1$LE+baseNS1$H)
baseNS1$EF[baseNS1$EF<0 | baseNS1$EF>1]=NA 

# read  radiation fluxes
baseNS1$RN=baseNS1$NETRAD
baseNS1$PA=101.3

# air temperature and humidity
baseNS1$VPD=baseNS1$VPD_PI/10

# get aerodynamic conductance
dfNS1 <- data.frame(Tair=baseNS1$TA_1_1_1,pressure=baseNS1$PA,wind=baseNS1$WS,ustar=baseNS1$USTAR,H=baseNS1$H)
GaNS1 = aerodynamic.conductance(dfNS1,Rb_model="Thom_1972")

# get aerodynamic temperature
surfNS1 = surface.conditions(Tair=baseNS1$TA_1_1_1,pressure=baseNS1$PA,LE=baseNS1$LE,H=baseNS1$H,VPD=baseNS1$VPD,Ga=GaNS1$Ga_h)
baseNS1$TSURF=surfNS1$Tsurf
ind = tsoutliers(baseNS1$TSURF)
baseNS1$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS1']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS1']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-NS1']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseNS1$H )/(1.293*1004*(baseNS1$TA +273.15)*baseNS1$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseNS1$T10 = baseNS1$TSURF-(baseNS1$H/(0.4*1.293*1004*baseNS1$USTAR))*(log(10/z0h)-theta)

uYEARNS1 = unique(baseNS1$YEAR)
s=0
NS1 <- data.frame(matrix(ncol = 18, nrow = length(uYEARNS1)*365))
colnames(NS1) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARNS1)-1

for (y in 1:length(uYEARNS1))
{
  for (n in 1:365)
  {
    s=s+1
    NS1$DOY[s]=n
    NS1$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseNS1$DOY == n & baseNS1$YEAR == sYEAR+y & baseNS1$HOUR >= 12 & baseNS1$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseNS1$DOY == n & baseNS1$YEAR == sYEAR+y & baseNS1$HOUR >= 6 & baseNS1$HOUR <= 18);
    ind2 = which(baseNS1$DOY == n & baseNS1$YEAR == sYEAR+y & baseNS1$HOUR >= 6 & baseNS1$HOUR <= 12);
    all = which(baseNS1$DOY == n & baseNS1$YEAR == sYEAR+y)
    mnt = baseNS1$MONTH[ind[1]]
    year = baseNS1$YEAR[ind[1]]
    
    NS1$BR[s]=sum(baseNS1$H[ind1])/sum(baseNS1$LE[ind1]) # Bowen ratio
    NS1$Ga[s]=median(GaNS1$Ga_h[ind], na.rm=T) # aerodynamic conductance
    NS1$LE[s]=median(baseNS1$LE[ind1], na.rm=T)
    NS1$VPD[s]=median(baseNS1$VPD[ind], na.rm=T)
    NS1$WS[s]=median(baseNS1$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    NS1$TA[s]=max(baseNS1$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    NS1$TS[s]=mean(baseNS1$TSURF[ind], na.rm=T)
    #NS1$LST[s]=mean(baseNS1$LST[ind], na.rm=T)
    # temperature gradient
    NS1$dT[s]=mean(baseNS1$TSURF[ind]-baseNS1$TA[ind], na.rm=T)
    NS1$dT10[s]=mean(baseNS1$TSURF[ind]-baseNS1$T10[ind], na.rm=T)
    
    # albedo
    NS1$ALB[s]=sum(baseNS1$SW_OUT[ind])/sum(baseNS1$SW_IN[ind])
    # evaporative fraction
    NS1$EF[s]=sum(baseNS1$LE[ind])/sum(baseNS1$LE[ind]+baseNS1$H[ind])
  }
}

NS1$BR[abs(NS1$BR)>6]=NA
NS1$ALB[abs(NS1$ALB)>1 | NS1$ALB<0.05]=NA
NS1$EF[NS1$EF<0]=NA
NS1$EF[NS1$EF>1]=NA
ind = tsoutliers(NS1$dT)
NS1$dT[ind$index]=NA
ind = tsoutliers(NS1$dT10)
NS1$dT10[ind$index]=NA
ind = tsoutliers(NS1$TS)
NS1$TS[ind$index]=NA
#####################################
### read CA-NS2 ###
#####################################
# read flux dataset
baseNS2 <- amf_read_base(file = ".../AMF_CA-NS2_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude of site
lat = 55.9058
lon = -98.5247

# extract time vector
baseNS2$Time <- format(as.POSIXct(baseNS2$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseNS2$Date <- format(as.Date(baseNS2$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseNS2$Time <- as.ITime(baseNS2$Time, tz = "Canada/Central")
baseNS2$Date <- as.IDate(baseNS2$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseNS2$EF = baseNS2$LE/(baseNS2$LE+baseNS2$H)
baseNS2$EF[baseNS2$EF<0 | baseNS2$EF>1]=NA 

# read  radiation fluxes
baseNS2$RN=baseNS2$NETRAD
baseNS2$PA=101.3

# air temperature and humidity
baseNS2$TA=baseNS2$TA_1_1_1
baseNS2$RH[baseNS2$RH>100]=100
baseNS2$VPD=baseNS2$VPD_PI/10

# get aerodynamic conductance
dfNS2 <- data.frame(Tair=baseNS2$TA,pressure=baseNS2$PA,wind=baseNS2$WS,ustar=baseNS2$USTAR,H=baseNS2$H)
GaNS2 = aerodynamic.conductance(dfNS2,Rb_model="Thom_1972")

# get aerodynamic temperature
surfNS2 = surface.conditions(Tair=baseNS2$TA,pressure=baseNS2$PA,LE=baseNS2$LE,H=baseNS2$H,VPD=baseNS2$VPD,Ga=GaNS2$Ga_h)
baseNS2$TSURF=surfNS2$Tsurf
ind = tsoutliers(baseNS2$TSURF)
baseNS2$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS2']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS2']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-NS2']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseNS2$H )/(1.293*1004*(baseNS2$TA +273.15)*baseNS2$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseNS2$T10 = baseNS2$TSURF-(baseNS2$H/(0.4*1.293*1004*baseNS2$USTAR))*(log(10/z0h)-theta)

uYEARNS2 = unique(baseNS2$YEAR)
s=0
NS2 <- data.frame(matrix(ncol = 18, nrow = length(uYEARNS2)*365))
colnames(NS2) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARNS2)-1

for (y in 1:length(uYEARNS2))
{
  for (n in 1:365)
  {
    s=s+1
    NS2$DOY[s]=n
    NS2$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseNS2$DOY == n & baseNS2$YEAR == sYEAR+y & baseNS2$HOUR >= 12 & baseNS2$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseNS2$DOY == n & baseNS2$YEAR == sYEAR+y & baseNS2$HOUR >= 6 & baseNS2$HOUR <= 18);
    ind2 = which(baseNS2$DOY == n & baseNS2$YEAR == sYEAR+y & baseNS2$HOUR >= 6 & baseNS2$HOUR <= 12);
    all = which(baseNS2$DOY == n & baseNS2$YEAR == sYEAR+y)
    mnt = baseNS2$MONTH[ind[1]]
    year = baseNS2$YEAR[ind[1]]
    
    NS2$BR[s]=sum(baseNS2$H[ind1])/sum(baseNS2$LE[ind1]) # Bowen ratio
    NS2$Ga[s]=median(GaNS2$Ga_h[ind], na.rm=T) # aerodynamic conductance
    NS2$LE[s]=median(baseNS2$LE[ind1], na.rm=T)
    NS2$VPD[s]=median(baseNS2$VPD[ind], na.rm=T)
    NS2$WS[s]=median(baseNS2$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    NS2$TA[s]=max(baseNS2$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    NS2$TS[s]=mean(baseNS2$TSURF[ind], na.rm=T)
    #NS2$LST[s]=mean(baseNS2$LST[ind], na.rm=T)
    # temperature gradient
    NS2$dT[s]=mean(baseNS2$TSURF[ind]-baseNS2$TA[ind], na.rm=T)
    NS2$dT10[s]=mean(baseNS2$TSURF[ind]-baseNS2$T10[ind], na.rm=T)
    
    # albedo
    NS2$ALB[s]=sum(baseNS2$SW_OUT[ind])/sum(baseNS2$SW_IN[ind])
    # evaporative fraction
    NS2$EF[s]=sum(baseNS2$LE[ind])/sum(baseNS2$LE[ind]+baseNS2$H[ind])
  }
}

NS2$BR[abs(NS2$BR)>6]=NA
NS2$ALB[abs(NS2$ALB)>1 | NS2$ALB<0.05]=NA
NS2$EF[NS2$EF<0]=NA
NS2$EF[NS2$EF>1]=NA
ind = tsoutliers(NS2$dT)
NS2$dT[ind$index]=NA
ind = tsoutliers(NS2$dT10)
NS2$dT10[ind$index]=NA
ind = tsoutliers(NS2$TS)
NS2$TS[ind$index]=NA
#####################################
### read CA-NS3 ###
#####################################
# latitude and longitude of site
lat = 55.9117
lon = -98.3822

# read flux dataset
baseNS3 <- amf_read_base(file = ".../AMF_CA-NS3_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# extract time vector
baseNS3$Time <- format(as.POSIXct(baseNS3$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseNS3$Date <- format(as.Date(baseNS3$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseNS3$Time <- as.ITime(baseNS3$Time, tz = "Canada/Central")
baseNS3$Date <- as.IDate(baseNS3$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseNS3$EF = baseNS3$LE/(baseNS3$LE+baseNS3$H)
baseNS3$EF[baseNS3$EF<0 | baseNS3$EF>1]=NA 

# read  radiation fluxes
baseNS3$RN=baseNS3$NETRAD
baseNS3$PA=101.3

# air temperature and humidity
baseNS3$TA=baseNS3$TA_1_1_1
baseNS3$RH[baseNS3$RH>100]=100
baseNS3$VPD=baseNS3$VPD_PI/10

# get aerodynamic conductance
dfNS3 <- data.frame(Tair=baseNS3$TA,pressure=baseNS3$PA,wind=baseNS3$WS,ustar=baseNS3$USTAR,H=baseNS3$H)
GaNS3 = aerodynamic.conductance(dfNS3,Rb_model="Thom_1972")

# get aerodynamic temperature
surfNS3 = surface.conditions(Tair=baseNS3$TA,pressure=baseNS3$PA,LE=baseNS3$LE,H=baseNS3$H,VPD=baseNS3$VPD,Ga=GaNS3$Ga_h)
baseNS3$TSURF=surfNS3$Tsurf
ind = tsoutliers(baseNS3$TSURF)
baseNS3$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS3']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS3']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-NS3']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseNS3$H )/(1.293*1004*(baseNS3$TA +273.15)*baseNS3$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseNS3$T10 = baseNS3$TSURF-(baseNS3$H/(0.4*1.293*1004*baseNS3$USTAR))*(log(10/z0h)-theta)

uYEARNS3 = unique(baseNS3$YEAR)
s=0
NS3 <- data.frame(matrix(ncol = 18, nrow = length(uYEARNS3)*365))
colnames(NS3) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARNS3)-1

for (y in 1:length(uYEARNS3))
{
  for (n in 1:365)
  {
    s=s+1
    NS3$DOY[s]=n
    NS3$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseNS3$DOY == n & baseNS3$YEAR == sYEAR+y & baseNS3$HOUR >= 12 & baseNS3$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseNS3$DOY == n & baseNS3$YEAR == sYEAR+y & baseNS3$HOUR >= 6 & baseNS3$HOUR <= 18);
    ind2 = which(baseNS3$DOY == n & baseNS3$YEAR == sYEAR+y & baseNS3$HOUR >= 6 & baseNS3$HOUR <= 12);
    all = which(baseNS3$DOY == n & baseNS3$YEAR == sYEAR+y);
    mnt = baseNS3$MONTH[ind[1]]
    year = baseNS3$YEAR[ind[1]]
    
    NS3$BR[s]=sum(baseNS3$H[ind1])/sum(baseNS3$LE[ind1]) # Bowen ratio
    NS3$Ga[s]=median(GaNS3$Ga_h[ind], na.rm=T) # aerodynamic conductance
    NS3$LE[s]=median(baseNS3$LE[ind1], na.rm=T)
    NS3$VPD[s]=median(baseNS3$VPD[ind], na.rm=T)
    NS3$WS[s]=median(baseNS3$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    NS3$TA[s]=max(baseNS3$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    NS3$TS[s]=mean(baseNS3$TSURF[ind], na.rm=T)
    #NS3$LST[s]=mean(baseNS3$LST[ind],na.rm=T)
    # temperature gradient
    NS3$dT[s]=mean(baseNS3$TSURF[ind]-baseNS3$TA[ind], na.rm=T)
    NS3$dT10[s]=mean(baseNS3$TSURF[ind]-baseNS3$T10[ind], na.rm=T)
    
    # albedo
    NS3$ALB[s]=sum(baseNS3$SW_OUT[ind])/sum(baseNS3$SW_IN[ind])
    # evaporative fraction
    NS3$EF[s]=sum(baseNS3$LE[ind])/sum(baseNS3$LE[ind]+baseNS3$H[ind])
  }
}

NS3$BR[abs(NS3$BR)>6]=NA
NS3$ALB[abs(NS3$ALB)>1 | NS3$ALB<0.05]=NA
NS3$EF[NS3$EF<0]=NA
NS3$EF[NS3$EF>1]=NA
ind = tsoutliers(NS3$dT)
NS3$dT[ind$index]=NA
ind = tsoutliers(NS3$dT10)
NS3$dT10[ind$index]=NA
ind = tsoutliers(NS3$TS)
NS3$TS[ind$index]=NA
#####################################
### read CA-NS4 ###
#####################################
# read flux dataset
baseNS4 <- amf_read_base(file = ".../AMF_CA-NS4_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude of site
lat = 55.9144
lon = -98.3806

# extract time vector
baseNS4$Time <- format(as.POSIXct(baseNS4$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseNS4$Date <- format(as.Date(baseNS4$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseNS4$Time <- as.ITime(baseNS4$Time, tz = "Canada/Central")
baseNS4$Date <- as.IDate(baseNS4$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseNS4$EF = baseNS4$LE/(baseNS4$LE+baseNS4$H)
baseNS4$EF[baseNS4$EF<0 | baseNS4$EF>1]=NA 

# read  radiation fluxes
baseNS4$RN=baseNS4$NETRAD
baseNS4$PA=101.3

# air temperature and humidity
baseNS4$TA=baseNS4$TA_1_1_1
baseNS4$RH[baseNS4$RH>100]=100
baseNS4$VPD=baseNS4$VPD_PI/10

# get aerodynamic conductance
dfNS4 <- data.frame(Tair=baseNS4$TA,pressure=baseNS4$PA,wind=baseNS4$WS,ustar=baseNS4$USTAR,H=baseNS4$H)
GaNS4 = aerodynamic.conductance(dfNS4,Rb_model="Thom_1972")

# get aerodynamic temperature
surfNS4 = surface.conditions(Tair=baseNS4$TA,pressure=baseNS4$PA,LE=baseNS4$LE,H=baseNS4$H,VPD=baseNS4$VPD,Ga=GaNS4$Ga_h)
baseNS4$TSURF=surfNS4$Tsurf
ind = tsoutliers(baseNS4$TSURF)
baseNS4$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS4']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS4']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-NS4']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseNS4$H )/(1.293*1004*(baseNS4$TA +273.15)*baseNS4$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseNS4$T10 = baseNS4$TSURF-(baseNS4$H/(0.4*1.293*1004*baseNS4$USTAR))*(log(10/z0h)-theta)

uYEARNS4 = unique(baseNS4$YEAR)
s=0
NS4 <- data.frame(matrix(ncol = 18, nrow = length(uYEARNS4)*365))
colnames(NS4) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARNS4)-1

for (y in 1:length(uYEARNS4))
{
  for (n in 1:365)
  {
    s=s+1
    NS4$DOY[s]=n
    NS4$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseNS4$DOY == n & baseNS4$YEAR == sYEAR+y & baseNS4$HOUR >= 12 & baseNS4$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseNS4$DOY == n & baseNS4$YEAR == sYEAR+y & baseNS4$HOUR >= 6 & baseNS4$HOUR <= 18);
    ind2 = which(baseNS4$DOY == n & baseNS4$YEAR == sYEAR+y & baseNS4$HOUR >= 6 & baseNS4$HOUR <= 12);
    all = which(baseNS4$DOY == n & baseNS4$YEAR == sYEAR+y)
    mnt = baseNS4$MONTH[ind[1]]
    year = baseNS4$YEAR[ind[1]]
    
    NS4$BR[s]=sum(baseNS4$H[ind1])/sum(baseNS4$LE[ind1]) # Bowen ratio
    NS4$Ga[s]=median(GaNS4$Ga_h[ind], na.rm=T) # aerodynamic conductance
    NS4$LE[s]=median(baseNS4$LE[ind1], na.rm=T)
    NS4$VPD[s]=median(baseNS4$VPD[ind], na.rm=T)
    NS4$WS[s]=median(baseNS4$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    NS4$TA[s]=max(baseNS4$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    NS4$TS[s]=mean(baseNS4$TSURF[ind], na.rm=T)
    #NS4$LST[s]=mean(baseNS4$LST[ind], na.rm=T)
    # temperature gradient
    NS4$dT[s]=mean(baseNS4$TSURF[ind]-baseNS4$TA[ind], na.rm=T)
    NS4$dT10[s]=mean(baseNS4$TSURF[ind]-baseNS4$T10[ind], na.rm=T)
    
    # albedo
    NS4$ALB[s]=sum(baseNS4$SW_OUT[ind])/sum(baseNS4$SW_IN[ind])
    # evaporative fraction
    NS4$EF[s]=sum(baseNS4$LE[ind])/sum(baseNS4$LE[ind]+baseNS4$H[ind])
  }
}

NS4$BR[abs(NS4$BR)>6]=NA
NS4$ALB[abs(NS4$ALB)>1 | NS4$ALB<0.05]=NA
NS4$EF[NS4$EF<0]=NA
NS4$EF[NS4$EF>1]=NA
ind = tsoutliers(NS4$dT)
NS4$dT[ind$index]=NA
ind = tsoutliers(NS4$dT10)
NS4$dT10[ind$index]=NA
ind = tsoutliers(NS4$TS)
NS4$TS[ind$index]=NA
#####################################
### read CA-NS6 ###
#####################################
# read flux dataset
baseNS6 <- amf_read_base(file = ".../AMF_CA-NS6_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude
lat = 55.9167
lon = -98.9644

# extract time vector
baseNS6$Time <- format(as.POSIXct(baseNS6$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseNS6$Date <- format(as.Date(baseNS6$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseNS6$Time <- as.ITime(baseNS6$Time, tz = "Canada/Central")
baseNS6$Date <- as.IDate(baseNS6$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseNS6$EF = baseNS6$LE/(baseNS6$LE+baseNS6$H)
baseNS6$EF[baseNS6$EF<0 | baseNS6$EF>1]=NA 

# read  radiation fluxes
baseNS6$RN=baseNS6$NETRAD
baseNS6$PA=101.3

# air temperature and humidity
baseNS6$TA=baseNS6$TA_1_1_1
baseNS6$RH[baseNS6$RH>100]=100
baseNS6$VPD=baseNS6$VPD_PI/10

# get aerodynamic conductance
dfNS6 <- data.frame(Tair=baseNS6$TA,pressure=baseNS6$PA,wind=baseNS6$WS,ustar=baseNS6$USTAR,H=baseNS6$H)
GaNS6 = aerodynamic.conductance(dfNS6,Rb_model="Thom_1972")

# get aerodynamic temperature
surfNS6 = surface.conditions(Tair=baseNS6$TA,pressure=baseNS6$PA,LE=baseNS6$LE,H=baseNS6$H,VPD=baseNS6$VPD,Ga=GaNS6$Ga_h)
baseNS6$TSURF=surfNS6$Tsurf
ind = tsoutliers(baseNS6$TSURF)
baseNS6$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS6']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS6']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-NS6']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseNS6$H )/(1.293*1004*(baseNS6$TA +273.15)*baseNS6$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseNS6$T10 = baseNS6$TSURF-(baseNS6$H/(0.4*1.293*1004*baseNS6$USTAR))*(log(10/z0h)-theta)

uYEARNS6 = unique(baseNS6$YEAR)
s=0
NS6 <- data.frame(matrix(ncol = 18, nrow = length(uYEARNS6)*365))
colnames(NS6) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARNS6)-1

for (y in 1:length(uYEARNS6))
{
  for (n in 1:365)
  {
    s=s+1
    NS6$DOY[s]=n
    NS6$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseNS6$DOY == n & baseNS6$YEAR == sYEAR+y & baseNS6$HOUR >= 12 & baseNS6$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseNS6$DOY == n & baseNS6$YEAR == sYEAR+y & baseNS6$HOUR >= 6 & baseNS6$HOUR <= 18);
    ind2 = which(baseNS6$DOY == n & baseNS6$YEAR == sYEAR+y & baseNS6$HOUR >= 6 & baseNS6$HOUR <= 12);
    all = which(baseNS6$DOY == n & baseNS6$YEAR == sYEAR+y)
    mnt = baseNS6$MONTH[ind[1]]
    year = baseNS6$YEAR[ind[1]]
    
    NS6$BR[s]=sum(baseNS6$H[ind1])/sum(baseNS6$LE[ind1]) # Bowen ratio
    NS6$Ga[s]=median(GaNS6$Ga_h[ind], na.rm=T) # aerodynamic conductance
    NS6$LE[s]=median(baseNS6$LE[ind1], na.rm=T)
    NS6$VPD[s]=median(baseNS6$VPD[ind], na.rm=T)
    NS6$WS[s]=median(baseNS6$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    NS6$TA[s]=max(baseNS6$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    NS6$TS[s]=mean(baseNS6$TSURF[ind], na.rm=T)
    #NS6$LST[s]=mean(baseNS6$LST[ind], na.rm=T)
    # temperature gradient
    NS6$dT[s]=mean(baseNS6$TSURF[ind]-baseNS6$TA[ind], na.rm=T)
    NS6$dT10[s]=mean(baseNS6$TSURF[ind]-baseNS6$T10[ind], na.rm=T)

    # albedo
    NS6$ALB[s]=sum(baseNS6$SW_OUT[ind])/sum(baseNS6$SW_IN[ind])
    # evaporative fraction
    NS6$EF[s]=sum(baseNS6$LE[ind])/sum(baseNS6$LE[ind]+baseNS6$H[ind])
  }
}

NS6$BR[abs(NS6$BR)>6]=NA
NS6$ALB[abs(NS6$ALB)>1 | NS6$ALB<0.05]=NA
NS6$EF[NS6$EF<0]=NA
NS6$EF[NS6$EF>1]=NA
ind = tsoutliers(NS6$dT)
NS6$dT[ind$index]=NA
ind = tsoutliers(NS6$dT10)
NS6$dT10[ind$index]=NA
ind = tsoutliers(NS6$TS)
NS6$TS[ind$index]=NA
#####################################
### read CA-NS7 ###
#####################################
# read flux dataset
baseNS7 <- amf_read_base(file = ".../AMF_CA-NS7_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
# latitude and longitude
lat = 56.6358
lon = -99.9483

# extract time vector
baseNS7$Time <- format(as.POSIXct(baseNS7$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseNS7$Date <- format(as.Date(baseNS7$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseNS7$Time <- as.ITime(baseNS7$Time, tz = "Canada/Central")
baseNS7$Date <- as.IDate(baseNS7$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseNS7$EF = baseNS7$LE/(baseNS7$LE+baseNS7$H)
baseNS7$EF[baseNS7$EF<0 | baseNS7$EF>1]=NA 

# read  radiation fluxes
baseNS7$RN=baseNS7$NETRAD
baseNS7$PA=101.3

# air temperature and humidity
baseNS7$TA=baseNS7$TA_1_1_1
baseNS7$RH[baseNS7$RH>100]=100
baseNS7$VPD=baseNS7$VPD_PI/10

# get aerodynamic conductance
dfNS7 <- data.frame(Tair=baseNS7$TA,pressure=baseNS7$PA,wind=baseNS7$WS,ustar=baseNS7$USTAR,H=baseNS7$H)
GaNS7 = aerodynamic.conductance(dfNS7,Rb_model="Thom_1972")

# get aerodynamic temperature
surfNS7 = surface.conditions(Tair=baseNS7$TA,pressure=baseNS7$PA,LE=baseNS7$LE,H=baseNS7$H,VPD=baseNS7$VPD,Ga=GaNS7$Ga_h)
baseNS7$TSURF=surfNS7$Tsurf
ind = tsoutliers(baseNS7$TSURF)
baseNS7$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS7']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-NS7']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-NS7']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseNS7$H )/(1.293*1004*(baseNS7$TA +273.15)*baseNS7$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseNS7$T10 = baseNS7$TSURF-(baseNS7$H/(0.4*1.293*1004*baseNS7$USTAR))*(log(10/z0h)-theta)

uYEARNS7 = unique(baseNS7$YEAR)
s=0
NS7 <- data.frame(matrix(ncol = 18, nrow = length(uYEARNS7)*365))
colnames(NS7) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARNS7)-1

for (y in 1:length(uYEARNS7))
{
  for (n in 1:365)
  {
    s=s+1
    NS7$DOY[s]=n
    NS7$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseNS7$DOY == n & baseNS7$YEAR == sYEAR+y & baseNS7$HOUR >= 12 & baseNS7$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseNS7$DOY == n & baseNS7$YEAR == sYEAR+y & baseNS7$HOUR >= 6 & baseNS7$HOUR <= 18);
    ind2 = which(baseNS7$DOY == n & baseNS7$YEAR == sYEAR+y & baseNS7$HOUR >= 6 & baseNS7$HOUR <= 12);
    all = which(baseNS7$DOY == n & baseNS7$YEAR == sYEAR+y)
    mnt = baseNS7$MONTH[ind[1]]
    year = baseNS7$YEAR[ind[1]]
    
    NS7$BR[s]=sum(baseNS7$H[ind1])/sum(baseNS7$LE[ind1]) # Bowen ratio
    NS7$Ga[s]=median(GaNS7$Ga_h[ind], na.rm=T) # aerodynamic conductance
    NS7$LE[s]=median(baseNS7$LE[ind1], na.rm=T)
    NS7$VPD[s]=median(baseNS7$VPD[ind], na.rm=T)
    NS7$WS[s]=median(baseNS7$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    NS7$TA[s]=max(baseNS7$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    NS7$TS[s]=mean(baseNS7$TSURF[ind], na.rm=T)
    #NS7$LST[s]=mean(baseNS7$LST[ind], na.rm=T)
    # temperature gradient
    NS7$dT[s]=mean(baseNS7$TSURF[ind]-baseNS7$TA[ind], na.rm=T)
    NS7$dT10[s]=mean(baseNS7$TSURF[ind]-baseNS7$T10[ind], na.rm=T)
    
    # albedo
    NS7$ALB[s]=sum(baseNS7$SW_OUT[ind])/sum(baseNS7$SW_IN[ind])
    # evaporative fraction
    NS7$EF[s]=sum(baseNS7$LE[ind])/sum(baseNS7$LE[ind]+baseNS7$H[ind])
  }
}

NS7$BR[abs(NS7$BR)>6]=NA
NS7$ALB[abs(NS7$ALB)>1 | NS7$ALB<0.05]=NA
NS7$EF[NS7$EF<0]=NA
NS7$EF[NS7$EF>1]=NA
ind = tsoutliers(NS7$dT)
NS7$dT[ind$index]=NA
ind = tsoutliers(NS7$dT10)
NS7$dT10[ind$index]=NA
ind = tsoutliers(NS7$TS)
NS7$TS[ind$index]=NA
#####################################
### read CA-SF3 ###
#####################################
# read flux dataset
baseSF3 <- amf_read_base(file = ".../AMF_CA-SF3_BASE-BADM_2-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)

baseSF31 <- amf_read_base(".../AMF_CA-SF3_BASE-BADM_2-5/AMF_CA-SF3_BASE_HH_2-5.csv",
                          unzip = FALSE,
                          parse_timestamp = FALSE)

# calculate air temperature at eddy covariance measurement height from sonic temperature
baseSF31$H2O_2=baseSF31$H2O_2*18.02/1000/1000
baseSF31$H2Oalt = baseSF31$H2O.1*18.02/10^6/0.0224
baseSF31$H2O_2[is.na(baseSF31$H2O_2)]=baseSF31$H2Oalt[is.na(baseSF31$H2O_2)]
baseSF31$H2O.1=baseSF31$H2O_2
baseSF31$e = baseSF31$H2O.1*461.5*(baseSF31$TA_1_1_1+273.15)
baseSF31$Pd = baseSF31$PA*1000 - baseSF31$e
baseSF31$rohd = baseSF31$Pd/(287*(baseSF3$TA_1_1_1+273.15))
baseSF31$TA_1_1_2=baseSF31$TSONIC*(1+0.51*(baseSF31$H2O.1/(baseSF31$rohd+baseSF31$H2O.1)))^(-1) # use air temp from sonic
baseSF3$TA_1_1_2 = baseSF31$TA_1_1_2
rm(baseSF31)

# latitude and longitude
lat = 54.0916
lon = -106.0053

# extract time vector
baseSF3$Time <- format(as.POSIXct(baseSF3$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseSF3$Date <- format(as.Date(baseSF3$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseSF3$Time <- as.ITime(baseSF3$Time, tz = "Canada/Central")
baseSF3$Date <- as.IDate(baseSF3$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseSF3$EF = baseSF3$LE/(baseSF3$LE+baseSF3$H)
baseSF3$EF[baseSF3$EF<0 | baseSF3$EF>1]=NA 

# read  radiation fluxes
baseSF3$RN=baseSF3$NETRAD_1_1_1

# air temperature and humidity
baseSF3$TA=baseSF3$TA_1_1_2
baseSF3$RH=baseSF3$RH
baseSF3$RH[baseSF3$RH>100]=100
baseSF3$VPD=baseSF3$VPD_PI/10

# get aerodynamic conductance
dfSF3 <- data.frame(Tair=baseSF3$TA,pressure=baseSF3$PA,wind=baseSF3$WS,ustar=baseSF3$USTAR,H=baseSF3$H)
GaSF3 = aerodynamic.conductance(dfSF3,Rb_model="Thom_1972")

# get aerodynamic temperature
surfSF3 = surface.conditions(Tair=baseSF3$TA,pressure=baseSF3$PA,LE=baseSF3$LE,H=baseSF3$H,VPD=baseSF3$VPD,Ga=GaSF3$Ga_h)
baseSF3$TSURF=surfSF3$Tsurf
ind = tsoutliers(baseSF3$TSURF)
baseSF3$TSURF[ind$index]=NA

# zero height displacement and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-SF3']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-SF3']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-SF3']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseSF3$H )/(1.293*1004*(baseSF3$TA +273.15)*baseSF3$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseSF3$T10 = baseSF3$TSURF-(baseSF3$H/(0.4*1.293*1004*baseSF3$USTAR))*(log(10/z0h)-theta)

uYEARSF3 = unique(baseSF3$YEAR)
s=0
SF3 <- data.frame(matrix(ncol = 18, nrow = length(uYEARSF3)*365))

colnames(SF3) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARSF3)-1

# remove data before August 2002 when measurement height was only at 7.7 m
ind=which(baseSF3$YEAR==2002 & baseSF3$MONTH == 12 & baseSF3$DAY ==31  & baseSF3$HOUR == 23 & baseSF3$MINUTE==45)
baseSF3$TSURF[ind:length(baseSF3$dTS)]=NA
baseSF3$T10[ind:length(baseSF3$dTS)]=NA
GaSF3$Ga_h[ind:length(baseSF3$dTS)]=NA

for (y in 1:length(uYEARSF3))
{
  for (n in 1:365)
  {
    s=s+1
    SF3$DOY[s]=n
    SF3$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseSF3$DOY == n & baseSF3$YEAR == sYEAR+y & baseSF3$HOUR >= 12 & baseSF3$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseSF3$DOY == n & baseSF3$YEAR == sYEAR+y & baseSF3$HOUR >= 6 & baseSF3$HOUR <= 18);
    ind2 = which(baseSF3$DOY == n & baseSF3$YEAR == sYEAR+y & baseSF3$HOUR >= 6 & baseSF3$HOUR <= 12);
    all = which(baseSF3$DOY == n & baseSF3$YEAR == sYEAR+y)
    mnt = baseSF3$MONTH[ind[1]]
    year = baseSF3$YEAR[ind[1]]
    
    SF3$BR[s]=sum(baseSF3$H[ind1])/sum(baseSF3$LE[ind1]) # Bowen ratio
    SF3$Ga[s]=median(GaSF3$Ga_h[ind], na.rm=T) # aerodynamic conductance
    SF3$LE[s]=median(baseSF3$LE[ind1], na.rm=T)
    SF3$VPD[s]=median(baseSF3$VPD[ind], na.rm=T)
    SF3$WS[s]=median(baseSF3$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    SF3$TA[s]=max(baseSF3$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    SF3$TS[s]=mean(baseSF3$TSURF[ind], na.rm=T)
    SF3$LST[s]=mean(baseSF3$LST[ind], na.rm=T)
    # temperature gradient
    SF3$dT[s]=mean(baseSF3$TSURF[ind]-baseSF3$TA[ind], na.rm=T)
    SF3$dT10[s]=mean(baseSF3$TSURF[ind]-baseSF3$T10[ind], na.rm=T)
    
    # albedo
    SF3$ALB[s]=sum(baseSF3$SW_OUT[ind])/sum(baseSF3$SW_IN[ind])
    # evaporative fraction
    SF3$EF[s]=sum(baseSF3$LE[ind])/sum(baseSF3$LE[ind]+baseSF3$H[ind])
  }
}

SF3$BR[abs(SF3$BR)>6]=NA
SF3$ALB[abs(SF3$ALB)>1 | SF3$ALB<0.05]=NA
SF3$EF[SF3$EF<0]=NA
SF3$EF[SF3$EF>1]=NA
ind = tsoutliers(SF3$dT)
SF3$dT[ind$index]=NA
ind = tsoutliers(SF3$dT10)
SF3$dT10[ind$index]=NA
ind = tsoutliers(SF3$TS)
SF3$TS[ind$index]=NA
#####################################
### read CA-SF2 ###
#####################################
# read flux dataset
baseSF2 <- amf_read_base(file = ".../AMF_CA-SF2_BASE-BADM_3-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)

baseSF21 <- amf_read_base(".../AMF_CA-SF2_BASE-BADM_3-5/AMF_CA-SF2_BASE_HH_3-5.csv",
                          unzip = FALSE,
                          parse_timestamp = FALSE)

# calculate air temperature at eddy covariance measurement height from sonic temperature
baseSF21$H2O.1=baseSF21$H2O.1*18.01/1000/1000
baseSF21$e = baseSF21$H2O.1*461.5*(baseSF2$TA_1_1_1+273.15)
baseSF21$Pd = baseSF21$PA*1000 - baseSF21$e
baseSF21$rohd = baseSF21$Pd/(287*(baseSF2$TA_1_1_1+273.15))
baseSF21$TA_1_1_2=baseSF21$TSONIC*(1+0.51*(baseSF21$H2O.1/(baseSF21$rohd+baseSF21$H2O.1)))^(-1) # use air temp from sonic
baseSF2$TA_1_1_2 = baseSF21$TA_1_1_2
rm(baseSF21)

# latitude and longitude
lat = 54.2539
lon = -105.8775
# extract time vector
baseSF2$Time <- format(as.POSIXct(baseSF2$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseSF2$Date <- format(as.Date(baseSF2$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseSF2$Time <- as.ITime(baseSF2$Time, tz = "Canada/Central")
baseSF2$Date <- as.IDate(baseSF2$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseSF2$EF = baseSF2$LE/(baseSF2$LE+baseSF2$H)
baseSF2$EF[baseSF2$EF<0 | baseSF2$EF>1]=NA 

# read  radiation fluxes
baseSF2$RN=baseSF2$NETRAD_1_1_1
baseSF2$PA=101.3

# air temperature and humidity
baseSF2$TA=baseSF2$TA_1_1_2
baseSF2$RH=baseSF2$RH
baseSF2$RH[baseSF2$RH>100]=100
baseSF2$VPD=baseSF2$VPD_PI/10

# get aerodynamic conductance
dfSF2 <- data.frame(Tair=baseSF2$TA,pressure=baseSF2$PA,wind=baseSF2$WS,ustar=baseSF2$USTAR,H=baseSF2$H)
GaSF2 = aerodynamic.conductance(dfSF2,Rb_model="Thom_1972")

# get aerodynamic temperature
surfSF2 = surface.conditions(Tair=baseSF2$TA,pressure=baseSF2$PA,LE=baseSF2$LE,H=baseSF2$H,VPD=baseSF2$VPD,Ga=GaSF2$Ga_h)
baseSF2$TSURF=surfSF2$Tsurf
ind = tsoutliers(baseSF2$TSURF)
baseSF2$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-SF2']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-SF2']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-SF2']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseSF2$H )/(1.293*1004*(baseSF2$TA +273.15)*baseSF2$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseSF2$T10 = baseSF2$TSURF-(baseSF2$H/(0.4*1.293*1004*baseSF2$USTAR))*(log(10/z0h)-theta)

uYEARSF2 = unique(baseSF2$YEAR)
s=0

SF2 <- data.frame(matrix(ncol = 18, nrow = length(uYEARSF2)*365))
colnames(SF2) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARSF2)-1

for (y in 1:length(uYEARSF2))
{
  for (n in 1:365)
  {
    s=s+1
    SF2$DOY[s]=n
    SF2$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseSF2$DOY == n & baseSF2$YEAR == sYEAR+y & baseSF2$HOUR >= 12 & baseSF2$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseSF2$DOY == n & baseSF2$YEAR == sYEAR+y & baseSF2$HOUR >= 6 & baseSF2$HOUR <= 18);
    ind2 = which(baseSF2$DOY == n & baseSF2$YEAR == sYEAR+y & baseSF2$HOUR >= 6 & baseSF2$HOUR <= 12);
    all = which(baseSF2$DOY == n & baseSF2$YEAR == sYEAR+y)
    mnt = baseSF2$MONTH[ind[1]]
    year = baseSF2$YEAR[ind[1]]
    
    SF2$BR[s]=sum(baseSF2$H[ind1])/sum(baseSF2$LE[ind1]) # Bowen ratio
    SF2$Ga[s]=median(GaSF2$Ga_h[ind], na.rm=T) # aerodynamic conductance
    SF2$LE[s]=median(baseSF2$LE[ind1], na.rm=T)
    SF2$VPD[s]=median(baseSF2$VPD[ind], na.rm=T)
    SF2$WS[s]=median(baseSF2$WS[ind], na.rm=T)

    # maximum air temperature in the afternoon
    SF2$TA[s]=max(baseSF2$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    SF2$TS[s]=mean(baseSF2$TSURF[ind], na.rm=T)
    # temperature gradient
    SF2$dT[s]=mean(baseSF2$TSURF[ind]-baseSF2$TA[ind], na.rm=T)
    SF2$dT10[s]=mean(baseSF2$TSURF[ind]-baseSF2$T10[ind], na.rm=T)
    
    # albedo
    SF2$ALB[s]=sum(baseSF2$SW_OUT[ind])/sum(baseSF2$SW_IN[ind])
    # evaporative fraction
    SF2$EF[s]=sum(baseSF2$LE[ind])/sum(baseSF2$LE[ind]+baseSF2$H[ind])
  }
}

SF2$BR[abs(SF2$BR)>6]=NA
SF2$ALB[abs(SF2$ALB)>1 | SF2$ALB<0.05]=NA
SF2$EF[SF2$EF<0]=NA
SF2$EF[SF2$EF>1]=NA
ind = tsoutliers(SF2$dT)
SF2$dT[ind$index]=NA
ind = tsoutliers(SF2$dT10)
SF2$dT10[ind$index]=NA
ind = tsoutliers(SF2$TS)
SF2$TS[ind$index]=NA
#####################################
### read CA-SF1 ###
#####################################
# read flux dataset
baseSF1 <- amf_read_base(file = ".../AMF_CA-SF1_BASE-BADM_2-5.zip",
                         unzip = TRUE,
                         parse_timestamp = TRUE)
baseSF11 <- amf_read_base(".../AMF_CA-SF1_BASE-BADM_2-5/AMF_CA-SF1_BASE_HH_2-5.csv",
                         unzip = FALSE,
                         parse_timestamp = FALSE)

# calculate air temperature at eddy covariance measurement height from sonic temperature
baseSF11$H2O.1=baseSF11$H2O.1*18.01/1000/1000
baseSF11$e = baseSF11$H2O.1*461.5*(baseSF1$TA_1_1_1+273.15)
baseSF11$Pd = baseSF11$PA*1000 - baseSF11$e
baseSF11$rohd = baseSF11$Pd/(287*(baseSF1$TA_1_1_1+273.15))
baseSF11$TA_1_1_2=baseSF11$TSONIC*(1+0.51*(baseSF11$H2O.1/(baseSF11$rohd+baseSF11$H2O.1)))^(-1) # use air temp from sonic
baseSF1$TA_1_1_2 = baseSF11$TA_1_1_2
rm(baseSF11)

# latitude and longitude
lat = 54.4850 
lon = -105.8176
# extract time vector
baseSF1$Time <- format(as.POSIXct(baseSF1$TIMESTAMP, "%Y-%m-%d %H:%M:%S", tz = "Canada/Central"), format = "%H:%M")
baseSF1$Date <- format(as.Date(baseSF1$TIMESTAMP,"%Y-%m-%d",tz = "Canada/Central"), format = "%Y-%m-%d")
baseSF1$Time <- as.ITime(baseSF1$Time, tz = "Canada/Central")
baseSF1$Date <- as.IDate(baseSF1$Date, tz = "Canada/Central")

# get latent and sensible heat flux
baseSF1$EF = baseSF1$LE/(baseSF1$LE+baseSF1$H)
baseSF1$EF[baseSF1$EF<0 | baseSF1$EF>1]=NA 

# read  radiation fluxes
baseSF1$RN=baseSF1$NETRAD

# air temperature and humidity
baseSF1$TA=baseSF1$TA_1_1_2
baseSF1$RH=baseSF1$RH
baseSF1$RH[baseSF1$RH>100]=100
baseSF1$VPD=baseSF1$VPD_PI/10

# get aerodynamic conductance
dfSF1 <- data.frame(Tair=baseSF1$TA,pressure=baseSF1$PA,wind=baseSF1$WS,ustar=baseSF1$USTAR,H=baseSF1$H)
GaSF1 = aerodynamic.conductance(dfSF1,Rb_model="Thom_1972")

# get aerodynamic temperature
surfSF1 = surface.conditions(Tair=baseSF1$TA,pressure=baseSF1$PA,LE=baseSF1$LE,H=baseSF1$H,VPD=baseSF1$VPD,Ga=GaSF1$Ga_h)
baseSF1$TSURF=surfSF1$Tsurf
ind = tsoutliers(baseSF1$TSURF)
baseSF1$TSURF[ind$index]=NA

# zero plane displacement height and roughness length
d0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-SF1']*0.66
z0 = HEIGHTS$hc[HEIGHTS$SITE=='CA-SF1']*0.1
z0h = z0/exp(2)
z=HEIGHTS$z[HEIGHTS$SITE=='CA-SF1']

# ratio of convective to mechanical production of TKE (see Novick & Katul 2020: eq A1)
zeta = (-0.4*(z-d0)*9.81*baseSF1$H )/(1.293*1004*(baseSF1$TA +273.15)*baseSF1$USTAR^3)
# diabatic correction factors (see Novick & Katul 2020: eq A2)
theta=2*log((1+(1-16*zeta)^0.5)/2)
theta[zeta>0 & !is.na(zeta)]=-6*log(1+zeta[zeta>0 & !is.na(zeta)])
theta[zeta>1]=NA
# get 10 m air temperature (see Novick & Katul 2020: eq 4)
baseSF1$T10 = baseSF1$TSURF-(baseSF1$H/(0.4*1.293*1004*baseSF1$USTAR))*(log(10/z0h)-theta)

uYEARSF1 = unique(baseSF1$YEAR)
s=0
baseSF1$ClearSky = NA
SF1 <- data.frame(matrix(ncol = 18, nrow = length(uYEARSF1)*365))
colnames(SF1) <- c('DOY', 'YEAR', 'BR','dT','TA','ALB','EF','Ga','WS',"Gs","Gs_grd","LE","VPD","OMG","TS","TCAN","ClearSky","SWC")
sYEAR = min(uYEARSF1)-1

for (y in 1:length(uYEARSF1))
{
  for (n in 1:365)
  {
    s=s+1
    SF1$DOY[s]=n
    SF1$YEAR[s]=sYEAR+y
    # select afternoon hours
    ind = which(baseSF1$DOY == n & baseSF1$YEAR == sYEAR+y & baseSF1$HOUR >= 12 & baseSF1$HOUR <= 16)
    # select all daytime hours
    ind1 = which(baseSF1$DOY == n & baseSF1$YEAR == sYEAR+y & baseSF1$HOUR >= 6 & baseSF1$HOUR <= 18);
    ind2 = which(baseSF1$DOY == n & baseSF1$YEAR == sYEAR+y & baseSF1$HOUR >= 6 & baseSF1$HOUR <= 12);
    all = which(baseSF1$DOY == n & baseSF1$YEAR == sYEAR+y)
    mnt = baseSF1$MONTH[ind[1]]
    year = baseSF1$YEAR[ind[1]]
    
    SF1$BR[s]=sum(baseSF1$H[ind1])/sum(baseSF1$LE[ind1]) # Bowen ratio
    SF1$Ga[s]=median(GaSF1$Ga_h[ind], na.rm=T) # aerodynamic conductance
    SF1$LE[s]=median(baseSF1$LE[ind1], na.rm=T)
    SF1$VPD[s]=median(baseSF1$VPD[ind], na.rm=T)
    SF1$WS[s]=median(baseSF1$WS[ind], na.rm=T)
    
    # maximum air temperature in the afternoon
    SF1$TA[s]=max(baseSF1$TA[ind], na.rm=T)
    # mean aerodynamic temperature in the afternoon
    SF1$TS[s]=mean(baseSF1$TSURF[ind], na.rm=T)
    # temperature gradient
    SF1$dT[s]=mean(baseSF1$TSURF[ind]-baseSF1$TA[ind], na.rm=T)
    SF1$dT10[s]=mean(baseSF1$TSURF[ind]-baseSF1$T10[ind], na.rm=T)
    
    # albedo
    SF1$ALB[s]=sum(baseSF1$SW_OUT[ind])/sum(baseSF1$SW_IN[ind])
    # evaporative fraction
    SF1$EF[s]=sum(baseSF1$LE[ind])/sum(baseSF1$LE[ind]+baseSF1$H[ind])
  }
}

SF1$BR[abs(SF1$BR)>6]=NA
SF1$ALB[abs(SF1$ALB)>1 | SF1$ALB<0.05]=NA
SF1$EF[SF1$EF<0]=NA
SF1$EF[SF1$EF>1]=NA
ind = tsoutliers(SF1$dT)
SF1$dT[ind$index]=NA
ind = tsoutliers(SF1$dT10)
SF1$dT10[ind$index]=NA
ind = tsoutliers(SF1$TS)
SF1$TS[ind$index]=NA

################################################################################
################# plot chronosequence of ecosystem properties ##################
################################################################################
###################### air-surface temperature difference ######################
library(matrixStats)
library(patchwork)
library(gridExtra)
library(wesanderson)

MODIS <- read.csv(".../MODIS_EC_WINTER_dTEMP.csv")
# select seasons #
sFeAp = 32
eFeAp = 120
sJuSe = 182
eJuSe = 273

RPF$SEASON=NA
RPF$SEASON[RPF$DOY>=sFeAp & RPF$DOY<=eFeAp]=1 # Feb - Apr
RPF$SEASON[RPF$DOY>=sJuSe & RPF$DOY<=eJuSe]=2 # Jul - Sep

dTS_RPF = aggregate(dT10 ~ SEASON+YEAR, RPF[!is.na(RPF$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004

dTS_RPF$SITE[dTS_RPF$YEAR<=11]=as.factor(1)
dTS_RPF$SITE[dTS_RPF$YEAR>11]=2
dTS_RPF$SITE=as.factor(dTS_RPF$SITE)
dTS_RPF$ECO[dTS_RPF$YEAR<=11]="OSH"
dTS_RPF$ECO[dTS_RPF$YEAR>11]="DBF"
dTS_RPF=dTS_RPF[dTS_RPF$SEASON==1,]
for (y in 1:length(dTS_RPF$YEAR))
{
  dTS_RPF$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_RPF$YEAR[y] & MODIS$SITE=='RPF']
}

FCR$SEASON=NA
FCR$SEASON[FCR$DOY>=sFeAp & FCR$DOY<=eFeAp]=1
FCR$SEASON[FCR$DOY>=sJuSe & FCR$DOY<=eJuSe]=2

dTS_FCR = aggregate(dT10 ~ SEASON+YEAR, FCR[!is.na(FCR$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(3)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"
dTS_FCR=dTS_FCR[dTS_FCR$SEASON==1,]
for (y in 1:length(dTS_FCR$YEAR))
{
  dTS_FCR$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_FCR$YEAR[y] & MODIS$SITE=='FCR']
}

BN3$SEASON=NA
BN3$SEASON[BN3$DOY>=sFeAp & BN3$DOY<=eFeAp]=1
BN3$SEASON[BN3$DOY>=sJuSe & BN3$DOY<=eJuSe]=2

dTS_BN3 = aggregate(dT10 ~ SEASON+YEAR, BN3[!is.na(BN3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(4)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"
dTS_BN3=dTS_BN3[dTS_BN3$SEASON==1,]
for (y in 1:length(dTS_BN3$YEAR))
{
  dTS_BN3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN3$YEAR[y] & MODIS$SITE=='BN3']
}


BN2$SEASON=NA
BN2$SEASON[BN2$DOY>=sFeAp & BN2$DOY<=eFeAp]=1
BN2$SEASON[BN2$DOY>=sJuSe & BN2$DOY<=eJuSe]=2

dTS_BN2 = aggregate(dT10 ~ SEASON+YEAR, BN2[!is.na(BN2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN2$SEASON=as.factor(dTS_BN2$SEASON)
dTS_BN2$SITE=as.factor(5)
dTS_BN2$YEAR=dTS_BN2$YEAR-1987
dTS_BN2$ECO="DBF"
dTS_BN2=dTS_BN2[dTS_BN2$SEASON==1,]
for (y in 1:length(dTS_BN2$YEAR))
{
  dTS_BN2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN2$YEAR[y] & MODIS$SITE=='BN2']
}

NS2$SEASON=NA
NS2$SEASON[NS2$DOY>=sFeAp & NS2$DOY<=eFeAp]=1
NS2$SEASON[NS2$DOY>=sJuSe & NS2$DOY<=eJuSe]=2

dTS_NS2 = aggregate(dT10 ~ SEASON+YEAR, NS2[!is.na(NS2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(6)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"
dTS_NS2=dTS_NS2[dTS_NS2$SEASON==1,]
for (y in 1:length(dTS_NS2$YEAR))
{
  dTS_NS2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS2$YEAR[y] & MODIS$SITE=='NS2']
}

NS3$SEASON=NA
NS3$SEASON[NS3$DOY>=sFeAp & NS3$DOY<=eFeAp]=1
NS3$SEASON[NS3$DOY>=sJuSe & NS3$DOY<=eJuSe]=2

dTS_NS3 = aggregate(dT10 ~ SEASON+YEAR, NS3[!is.na(NS3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(7)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"
dTS_NS3=dTS_NS3[dTS_NS3$SEASON==1,]
for (y in 1:length(dTS_NS3$YEAR))
{
  dTS_NS3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS3$YEAR[y] & MODIS$SITE=='NS3']
}

NS4$SEASON=NA
NS4$SEASON[NS4$DOY>=sFeAp & NS4$DOY<=eFeAp]=1
NS4$SEASON[NS4$DOY>=sJuSe & NS4$DOY<=eJuSe]=2

dTS_NS4 = aggregate(dT10 ~ SEASON+YEAR, NS4[!is.na(NS4$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(8)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"
dTS_NS4=dTS_NS4[dTS_NS4$SEASON==1,]
for (y in 1:length(dTS_NS4$YEAR))
{
  dTS_NS4$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS4$YEAR[y] & MODIS$SITE=='NS4']
}

NS5$SEASON=NA
NS5$SEASON[NS5$DOY>=sFeAp & NS5$DOY<=eFeAp]=1
NS5$SEASON[NS5$DOY>=sJuSe & NS5$DOY<=eJuSe]=2

dTS_NS5 = aggregate(dT10 ~ SEASON+YEAR, NS5[!is.na(NS5$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(9)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"
dTS_NS5=dTS_NS5[dTS_NS5$SEASON==1,]
for (y in 1:length(dTS_NS5$YEAR))
{
  dTS_NS5$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS5$YEAR[y] & MODIS$SITE=='NS5']
}

NS6$SEASON=NA
NS6$SEASON[NS6$DOY>=sFeAp & NS6$DOY<=eFeAp]=1
NS6$SEASON[NS6$DOY>=sJuSe & NS6$DOY<=eJuSe]=2

dTS_NS6 = aggregate(dT10 ~ SEASON+YEAR, NS6[!is.na(NS6$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(10)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"
dTS_NS6=dTS_NS6[dTS_NS6$SEASON==1,]
for (y in 1:length(dTS_NS6$YEAR))
{
  dTS_NS6$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS6$YEAR[y] & MODIS$SITE=='NS6']
}

NS7$SEASON=NA
NS7$SEASON[NS7$DOY>=sFeAp & NS7$DOY<=eFeAp]=1
NS7$SEASON[NS7$DOY>=sJuSe & NS7$DOY<=eJuSe]=2

dTS_NS7 = aggregate(dT10 ~ SEASON+YEAR, NS7[!is.na(NS7$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(11)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"
dTS_NS7=dTS_NS7[dTS_NS7$SEASON==1,]
for (y in 1:length(dTS_NS7$YEAR))
{
  dTS_NS7$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS7$YEAR[y] & MODIS$SITE=='NS7']
}


SF3$SEASON=NA
SF3$SEASON[SF3$DOY>=sFeAp & SF3$DOY<=eFeAp]=1
SF3$SEASON[SF3$DOY>=sJuSe & SF3$DOY<=eJuSe]=2

dTS_SF3 = aggregate(dT10 ~ SEASON+YEAR, SF3[!is.na(SF3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(12)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"
dTS_SF3=dTS_SF3[dTS_SF3$SEASON==1,]
for (y in 1:length(dTS_SF3$YEAR))
{
  dTS_SF3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF3$YEAR[y] & MODIS$SITE=='SF3']
}

SF2$SEASON=NA
SF2$SEASON[SF2$DOY>=sFeAp & SF2$DOY<=eFeAp]=1
SF2$SEASON[SF2$DOY>=sJuSe & SF2$DOY<=eJuSe]=2

dTS_SF2 = aggregate(dT10 ~ SEASON+YEAR, SF2[!is.na(SF2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(13)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"
dTS_SF2=dTS_SF2[dTS_SF2$SEASON==1,]
for (y in 1:length(dTS_SF2$YEAR))
{
  dTS_SF2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF2$YEAR[y] & MODIS$SITE=='SF2']
}

SF1$SEASON=NA
SF1$SEASON[SF1$DOY>=sFeAp & SF1$DOY<=eFeAp]=1
SF1$SEASON[SF1$DOY>=sJuSe & SF1$DOY<=eJuSe]=2

dTS_SF1 = aggregate(dT10 ~ SEASON+YEAR, SF1[!is.na(SF1$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(14)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"
dTS_SF1=dTS_SF1[dTS_SF1$SEASON==1,]
for (y in 1:length(dTS_SF1$YEAR))
{
  dTS_SF1$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF1$YEAR[y] & MODIS$SITE=='SF1']
}

###### plot successional winter changes in air-surface temperature difference #####
# bind dataframes from individual sites (omit CA-SF3 and CA-NS7: sites with standing dead trees)
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) 

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$dT10[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# get ecotype by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$dT10[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC2 = group_mean
dT_SUCC2$ECO=group_mean2$x
names(dT_SUCC2)=c("SITE","dT","YEAR","SD","dY","ECO")

# bind dataframes from individual sites (all sites)
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF1,dTS_SF2,dTS_SF3) # ,dTS_BN2,dTS_SF3,dTS_SF2,dTS_SF1

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$dT10[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$dT10[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC1 = group_mean
dT_SUCC1$ECO=group_mean2$x
names(dT_SUCC1)=c("SITE","dT","YEAR","SD","dY","ECO")

# fit logarithmic function
fit <- nls(dT~ SSasymp(YEAR, yf, y0, log_alpha), data = dT_SUCC1)
r <- range(dT_SUCC1$YEAR)
x = predict(fit)
mat = matrix(ncol = 2, nrow = 200) 
fitted=data.frame(mat)
fitted$YEAR <- seq(1,95,length.out = 200)
fitted$fit <- predict(fit,newdata=fitted)
# get correlation between model fit and data
cor = cor.test(x, dT_SUCC1$dT, method=c("pearson", "kendall", "spearman"))

# plot successional change in air-surface temperature gradient
pal <- wes_palette("FantasticFox1", 14, type = "continuous")
fig1 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = YEAR, y = dT10[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = dT, ymin = dT - SD,
                                        ymax = dT + SD), width=0) + # plot site means and errors
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("Years since fire") +
  ylab("\u0394 T (surface - air) [\u00B0C]")+
  ggtitle("February to April")+
  annotate("text", x=77, y=2.9, label= "(a)")+
  theme_bw()+
  geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = dT, shape=ECO),size=3,colour='black', fill='grey')+
  geom_point(dT_SUCC2, mapping = aes(x = YEAR, y = dT, shape=ECO, fill=SITE),size=3,colour='black')+
  geom_line(fitted, mapping = aes(x = YEAR, y = fit))+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(fitted, mapping = aes(x = YEAR, y = fit))+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  scale_y_continuous(expand = expand_scale(),limits = c(-0.5,3))+
  scale_x_continuous(expand = expand_scale(),limits = c(0,81))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

######### summer data ##############
MODIS <- read.csv(".../MODIS_EC_SUMMER_dTEMP.csv")
dTS_RPF = aggregate(dT10 ~ SEASON+YEAR, RPF[!is.na(RPF$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004
dTS_RPF$SITE[dTS_RPF$YEAR<=11]=as.factor(1)
dTS_RPF$SITE[dTS_RPF$YEAR>11]=2
dTS_RPF$SITE=as.factor(dTS_RPF$SITE)
dTS_RPF$ECO[dTS_RPF$YEAR<=11]="OSH"
dTS_RPF$ECO[dTS_RPF$YEAR>11]="DBF"
dTS_RPF=dTS_RPF[dTS_RPF$SEASON==2,]
for (y in 1:length(dTS_RPF$YEAR))
{
  dTS_RPF$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_RPF$YEAR[y] & MODIS$SITE=='RPF']
}

FCR$SEASON=NA
FCR$SEASON[FCR$DOY>=sFeAp & FCR$DOY<=eFeAp]=1
FCR$SEASON[FCR$DOY>=sJuSe & FCR$DOY<=eJuSe]=2

dTS_FCR = aggregate(dT10 ~ SEASON+YEAR, FCR[!is.na(FCR$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(3)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"
dTS_FCR=dTS_FCR[dTS_FCR$SEASON==2,]
for (y in 1:length(dTS_FCR$YEAR))
{
  dTS_FCR$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_FCR$YEAR[y] & MODIS$SITE=='FCR']
}

BN3$SEASON=NA
BN3$SEASON[BN3$DOY>=sFeAp & BN3$DOY<=eFeAp]=1
BN3$SEASON[BN3$DOY>=sJuSe & BN3$DOY<=eJuSe]=2

dTS_BN3 = aggregate(dT10 ~ SEASON+YEAR, BN3[!is.na(BN3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(4)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"
dTS_BN3=dTS_BN3[dTS_BN3$SEASON==2,]
for (y in 1:length(dTS_BN3$YEAR))
{
  dTS_BN3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN3$YEAR[y] & MODIS$SITE=='BN3']
}


BN2$SEASON=NA
BN2$SEASON[BN2$DOY>=sFeAp & BN2$DOY<=eFeAp]=1
BN2$SEASON[BN2$DOY>=sJuSe & BN2$DOY<=eJuSe]=2

dTS_BN2 = aggregate(dT10 ~ SEASON+YEAR, BN2[!is.na(BN2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN2$SEASON=as.factor(dTS_BN2$SEASON)
dTS_BN2$SITE=as.factor(5)
dTS_BN2$YEAR=dTS_BN2$YEAR-1987
dTS_BN2$ECO="DBF"
dTS_BN2=dTS_BN2[dTS_BN2$SEASON==2,]
for (y in 1:length(dTS_BN2$YEAR))
{
  dTS_BN2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN2$YEAR[y] & MODIS$SITE=='BN2']
}

NS2$SEASON=NA
NS2$SEASON[NS2$DOY>=sFeAp & NS2$DOY<=eFeAp]=1
NS2$SEASON[NS2$DOY>=sJuSe & NS2$DOY<=eJuSe]=2

dTS_NS2 = aggregate(dT10 ~ SEASON+YEAR, NS2[!is.na(NS2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(6)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"
dTS_NS2=dTS_NS2[dTS_NS2$SEASON==2,]
for (y in 1:length(dTS_NS2$YEAR))
{
  dTS_NS2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS2$YEAR[y] & MODIS$SITE=='NS2']
}

NS3$SEASON=NA
NS3$SEASON[NS3$DOY>=sFeAp & NS3$DOY<=eFeAp]=1
NS3$SEASON[NS3$DOY>=sJuSe & NS3$DOY<=eJuSe]=2

dTS_NS3 = aggregate(dT10 ~ SEASON+YEAR, NS3[!is.na(NS3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(7)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"
dTS_NS3=dTS_NS3[dTS_NS3$SEASON==2,]
for (y in 1:length(dTS_NS3$YEAR))
{
  dTS_NS3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS3$YEAR[y] & MODIS$SITE=='NS3']
}

NS4$SEASON=NA
NS4$SEASON[NS4$DOY>=sFeAp & NS4$DOY<=eFeAp]=1
NS4$SEASON[NS4$DOY>=sJuSe & NS4$DOY<=eJuSe]=2

dTS_NS4 = aggregate(dT10 ~ SEASON+YEAR, NS4[!is.na(NS4$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(8)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"
dTS_NS4=dTS_NS4[dTS_NS4$SEASON==2,]
for (y in 1:length(dTS_NS4$YEAR))
{
  dTS_NS4$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS4$YEAR[y] & MODIS$SITE=='NS4']
}

NS5$SEASON=NA
NS5$SEASON[NS5$DOY>=sFeAp & NS5$DOY<=eFeAp]=1
NS5$SEASON[NS5$DOY>=sJuSe & NS5$DOY<=eJuSe]=2

dTS_NS5 = aggregate(dT10 ~ SEASON+YEAR, NS5[!is.na(NS5$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(9)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"
dTS_NS5=dTS_NS5[dTS_NS5$SEASON==2,]
for (y in 1:length(dTS_NS5$YEAR))
{
  dTS_NS5$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS5$YEAR[y] & MODIS$SITE=='NS5']
}

NS6$SEASON=NA
NS6$SEASON[NS6$DOY>=sFeAp & NS6$DOY<=eFeAp]=1
NS6$SEASON[NS6$DOY>=sJuSe & NS6$DOY<=eJuSe]=2

dTS_NS6 = aggregate(dT10 ~ SEASON+YEAR, NS6[!is.na(NS6$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(10)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"
dTS_NS6=dTS_NS6[dTS_NS6$SEASON==2,]
for (y in 1:length(dTS_NS6$YEAR))
{
  dTS_NS6$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS6$YEAR[y] & MODIS$SITE=='NS6']
}

NS7$SEASON=NA
NS7$SEASON[NS7$DOY>=sFeAp & NS7$DOY<=eFeAp]=1
NS7$SEASON[NS7$DOY>=sJuSe & NS7$DOY<=eJuSe]=2

dTS_NS7 = aggregate(dT10 ~ SEASON+YEAR, NS7[!is.na(NS7$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(11)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"
dTS_NS7=dTS_NS7[dTS_NS7$SEASON==2,]
for (y in 1:length(dTS_NS7$YEAR))
{
  dTS_NS7$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS7$YEAR[y] & MODIS$SITE=='NS7']
}


SF3$SEASON=NA
SF3$SEASON[SF3$DOY>=sFeAp & SF3$DOY<=eFeAp]=1
SF3$SEASON[SF3$DOY>=sJuSe & SF3$DOY<=eJuSe]=2

dTS_SF3 = aggregate(dT10 ~ SEASON+YEAR, SF3[!is.na(SF3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(12)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"
dTS_SF3=dTS_SF3[dTS_SF3$SEASON==2,]
for (y in 1:length(dTS_SF3$YEAR))
{
  dTS_SF3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF3$YEAR[y] & MODIS$SITE=='SF3']
}

SF2$SEASON=NA
SF2$SEASON[SF2$DOY>=sFeAp & SF2$DOY<=eFeAp]=1
SF2$SEASON[SF2$DOY>=sJuSe & SF2$DOY<=eJuSe]=2

dTS_SF2 = aggregate(dT10 ~ SEASON+YEAR, SF2[!is.na(SF2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(13)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"
dTS_SF2=dTS_SF2[dTS_SF2$SEASON==2,]
for (y in 1:length(dTS_SF2$YEAR))
{
  dTS_SF2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF2$YEAR[y] & MODIS$SITE=='SF2']
}

SF1$SEASON=NA
SF1$SEASON[SF1$DOY>=sFeAp & SF1$DOY<=eFeAp]=1
SF1$SEASON[SF1$DOY>=sJuSe & SF1$DOY<=eJuSe]=2

dTS_SF1 = aggregate(dT10 ~ SEASON+YEAR, SF1[!is.na(SF1$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(14)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"
dTS_SF1=dTS_SF1[dTS_SF1$SEASON==2,]
for (y in 1:length(dTS_SF1$YEAR))
{
  dTS_SF1$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF1$YEAR[y] & MODIS$SITE=='SF1']
}

###### plot successional winter changes in air-surface temperature difference #####
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) # omit sites with standing dead trees

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$dT10[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$dT10[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC2 = group_mean
dT_SUCC2$ECO=group_mean2$x
names(dT_SUCC2)=c("SITE","dT","YEAR","SD","dY","ECO")

# bind dataframes from individual sites
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF1,dTS_SF2,dTS_SF3) # all sites

# mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$dT10[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$dT10[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2

dT_SUCC1 = group_mean
dT_SUCC1$ECO=group_mean2$x
names(dT_SUCC1)=c("SITE","dT","YEAR","SD","dY","ECO")

# fit logarithmic function
fit <- nls(dT ~ SSasymp(YEAR, yf, y0, log_alpha), data = dT_SUCC1)
r <- range(dT_SUCC1$YEAR)
x = predict(fit)
mat = matrix(ncol = 2, nrow = 200) 
fitted=data.frame(mat)
fitted$YEAR <- seq(1,95,length.out = 200)
fitted$fit <- predict(fit,newdata=fitted)
# get correlation between model fit and data
cor = cor.test(x, dT_SUCC1$dT, method=c("pearson", "kendall", "spearman"))

# plot successional change
pal <- wes_palette("FantasticFox1", 14, type = "continuous")
fig2 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = YEAR, y = dT10[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = dT, ymin = dT - SD,
                                        ymax = dT + SD), width=0) + # plot site means and errors
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("Years since fire") +
  ylab(element_blank())+
  ggtitle("July to September")+
  annotate("text", x=77, y=4.2, label= "(b)")+
  theme_bw()+
  geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = dT, shape=ECO),size=3,colour='black', fill='grey')+
  geom_point(dT_SUCC2, mapping = aes(x = YEAR, y = dT, shape=ECO, fill=SITE),size=3,colour='black')+
  geom_line(fitted, mapping = aes(x = YEAR, y = fit))+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(fitted, mapping = aes(x = YEAR, y = fit))+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  scale_y_continuous(expand = expand_scale(),limits = c(0,4.3))+
  scale_x_continuous(expand = expand_scale(),limits = c(0,81))+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = "none",colour = "none")+
  labs(shape='Vegetation type') 

plot1 = grid.arrange(fig1, fig2, nrow = 1)
plot1
#ggsave("Fig5_Wildfire.svg", width =20, height = 10, units = "cm", plot1)

################################################################################
################# plot chronosequence of albedo changes ##################
################################################################################
ALBEDO <- read.csv(".../ALB_MODIS_FIRE_WINTER.csv")

library(matrixStats)
library(patchwork)

RPF$SEASON=NA
RPF$SEASON[RPF$DOY>=sFeAp & RPF$DOY<=eFeAp]=1 
RPF$SEASON[RPF$DOY>=sJuSe & RPF$DOY<=eJuSe]=2 

dTS_RPF = aggregate(ALB ~ SEASON+YEAR, RPF[!is.na(RPF$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004
dTS_RPF$SITE[dTS_RPF$YEAR<=11]=as.factor(1)
dTS_RPF$SITE[dTS_RPF$YEAR>11]=2
dTS_RPF$SITE=as.factor(dTS_RPF$SITE)
dTS_RPF$ECO[dTS_RPF$YEAR<=11]="OSH"
dTS_RPF$ECO[dTS_RPF$YEAR>11]="DBF"

FCR$SEASON=NA
FCR$SEASON[FCR$DOY>=sFeAp & FCR$DOY<=eFeAp]=1
FCR$SEASON[FCR$DOY>=sJuSe & FCR$DOY<=eJuSe]=2

dTS_FCR = aggregate(ALB ~ SEASON+YEAR, FCR[!is.na(FCR$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(3)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"

BN3$SEASON=NA
BN3$SEASON[BN3$DOY>=sFeAp & BN3$DOY<=eFeAp]=1
BN3$SEASON[BN3$DOY>=sJuSe & BN3$DOY<=eJuSe]=2

dTS_BN3 = aggregate(ALB ~ SEASON+YEAR, BN3[!is.na(BN3$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(4)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"

NS1$SEASON=NA
NS1$SEASON[NS1$DOY>=sFeAp & NS1$DOY<=eFeAp]=1
NS1$SEASON[NS1$DOY>=sJuSe & NS1$DOY<=eJuSe]=2

dTS_NS1 = aggregate(ALB ~ SEASON+YEAR, NS1[!is.na(NS1$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS1$SEASON=as.factor(dTS_NS1$SEASON)
dTS_NS1$SITE=as.factor(6)
dTS_NS1$YEAR=dTS_NS1$YEAR-1850
dTS_NS1$ECO="ENF"

NS2$SEASON=NA
NS2$SEASON[NS2$DOY>=sFeAp & NS2$DOY<=eFeAp]=1
NS2$SEASON[NS2$DOY>=sJuSe & NS2$DOY<=eJuSe]=2

dTS_NS2 = aggregate(ALB ~ SEASON+YEAR, NS2[!is.na(NS2$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(7)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"

NS3$SEASON=NA
NS3$SEASON[NS3$DOY>=sFeAp & NS3$DOY<=eFeAp]=1
NS3$SEASON[NS3$DOY>=sJuSe & NS3$DOY<=eJuSe]=2

dTS_NS3 = aggregate(ALB ~ SEASON+YEAR, NS3[!is.na(NS3$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(8)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"

NS4$SEASON=NA
NS4$SEASON[NS4$DOY>=sFeAp & NS4$DOY<=eFeAp]=1
NS4$SEASON[NS4$DOY>=sJuSe & NS4$DOY<=eJuSe]=2

dTS_NS4 = aggregate(ALB ~ SEASON+YEAR, NS4[!is.na(NS4$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(9)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"

NS5$SEASON=NA
NS5$SEASON[NS5$DOY>=sFeAp & NS5$DOY<=eFeAp]=1
NS5$SEASON[NS5$DOY>=sJuSe & NS5$DOY<=eJuSe]=2

dTS_NS5 = aggregate(ALB ~ SEASON+YEAR, NS5[!is.na(NS5$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(10)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"

NS6$SEASON=NA
NS6$SEASON[NS6$DOY>=sFeAp & NS6$DOY<=eFeAp]=1
NS6$SEASON[NS6$DOY>=sJuSe & NS6$DOY<=eJuSe]=2

dTS_NS6 = aggregate(ALB ~ SEASON+YEAR, NS6[!is.na(NS6$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(11)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"

NS7$SEASON=NA
NS7$SEASON[NS7$DOY>=sFeAp & NS7$DOY<=eFeAp]=1
NS7$SEASON[NS7$DOY>=sJuSe & NS7$DOY<=eJuSe]=2

dTS_NS7 = aggregate(ALB ~ SEASON+YEAR, NS7[!is.na(NS7$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(12)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"

SF3$SEASON=NA
SF3$SEASON[SF3$DOY>=sFeAp & SF3$DOY<=eFeAp]=1
SF3$SEASON[SF3$DOY>=sJuSe & SF3$DOY<=eJuSe]=2

dTS_SF3 = aggregate(ALB ~ SEASON+YEAR, SF3[!is.na(SF3$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(13)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"

SF2$SEASON=NA
SF2$SEASON[SF2$DOY>=sFeAp & SF2$DOY<=eFeAp]=1
SF2$SEASON[SF2$DOY>=sJuSe & SF2$DOY<=eJuSe]=2

dTS_SF2 = aggregate(ALB ~ SEASON+YEAR, SF2[!is.na(SF2$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(14)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"

SF1$SEASON=NA
SF1$SEASON[SF1$DOY>=sFeAp & SF1$DOY<=eFeAp]=1
SF1$SEASON[SF1$DOY>=sJuSe & SF1$DOY<=eJuSe]=2

dTS_SF1 = aggregate(ALB ~ SEASON+YEAR, SF1[!is.na(SF1$ALB),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(15)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"

# MODIS albedo for sites
MODIS_RPF=cbind(ALBEDO$FireYear[ALBEDO$SITE=="RPF"],ALBEDO$dALB[ALBEDO$SITE=="RPF"])
dTS_RPF$dALB=NA
for (x in 1:length(dTS_RPF$YEAR)) {
  dTS_RPF$dALB[x] = MODIS_RPF[MODIS_RPF[,1]==dTS_RPF$YEAR[x],2]
}

MODIS_FCR=cbind(ALBEDO$FireYear[ALBEDO$SITE=="FCR"],ALBEDO$dALB[ALBEDO$SITE=="FCR"])
dTS_FCR$dALB=NA
for (x in 1:length(dTS_FCR$YEAR)) {
  dTS_FCR$dALB[x] = MODIS_FCR[MODIS_FCR[,1]==dTS_FCR$YEAR[x],2]
}

MODIS_BN3=cbind(ALBEDO$FireYear[ALBEDO$SITE=="BN3"],ALBEDO$dALB[ALBEDO$SITE=="BN3"])
dTS_BN3$dALB=NA
for (x in 1:length(dTS_BN3$YEAR)) {
  dTS_BN3$dALB[x] = MODIS_BN3[MODIS_BN3[,1]==dTS_BN3$YEAR[x],2]
}

MODIS_BN2=cbind(ALBEDO$FireYear[ALBEDO$SITE=="BN2"],ALBEDO$dALB[ALBEDO$SITE=="BN2"])
dTS_BN2$dALB=NA
for (x in 1:length(dTS_BN2$YEAR)) {
  dTS_BN2$dALB[x] = MODIS_BN2[MODIS_BN2[,1]==dTS_BN2$YEAR[x],2]
}

MODIS_NS2=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS2"],ALBEDO$dALB[ALBEDO$SITE=="NS2"])
dTS_NS2$dALB=NA
for (x in 1:length(dTS_NS2$YEAR)) {
  dTS_NS2$dALB[x] = MODIS_NS2[MODIS_NS2[,1]==dTS_NS2$YEAR[x],2]
}

MODIS_NS3=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS3"],ALBEDO$dALB[ALBEDO$SITE=="NS3"])
dTS_NS3$dALB=NA
for (x in 1:length(dTS_NS3$YEAR)) {
  dTS_NS3$dALB[x] = MODIS_NS3[MODIS_NS3[,1]==dTS_NS3$YEAR[x],2]
}

MODIS_NS4=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS4"],ALBEDO$dALB[ALBEDO$SITE=="NS4"])
dTS_NS4$dALB=NA
for (x in 1:length(dTS_NS4$YEAR)) {
  dTS_NS4$dALB[x] = MODIS_NS4[MODIS_NS4[,1]==dTS_NS4$YEAR[x],2]
}

MODIS_NS5=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS5"],ALBEDO$dALB[ALBEDO$SITE=="NS5"])
dTS_NS5$dALB=NA
for (x in 1:length(dTS_NS5$YEAR)) {
  dTS_NS5$dALB[x] = MODIS_NS5[MODIS_NS5[,1]==dTS_NS5$YEAR[x],2]
}

MODIS_NS6=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS6"],ALBEDO$dALB[ALBEDO$SITE=="NS6"])
dTS_NS6$dALB=NA
for (x in 1:length(dTS_NS6$YEAR)) {
  dTS_NS6$dALB[x] = MODIS_NS6[MODIS_NS6[,1]==dTS_NS6$YEAR[x],2]
}

MODIS_NS7=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS7"],ALBEDO$dALB[ALBEDO$SITE=="NS7"])
dTS_NS7$dALB=NA
for (x in 1:length(dTS_NS7$YEAR)) {
  dTS_NS7$dALB[x] = MODIS_NS7[MODIS_NS7[,1]==dTS_NS7$YEAR[x],2]
}

MODIS_SF3=cbind(ALBEDO$FireYear[ALBEDO$SITE=="SF3"],ALBEDO$dALB[ALBEDO$SITE=="SF3"])
dTS_SF3$dALB=NA
for (x in 1:length(dTS_SF3$YEAR)) {
  dTS_SF3$dALB[x] = MODIS_SF3[MODIS_SF3[,1]==dTS_SF3$YEAR[x],2]
}

MODIS_SF2=cbind(ALBEDO$FireYear[ALBEDO$SITE=="SF2"],ALBEDO$dALB[ALBEDO$SITE=="SF2"])
dTS_SF2$dALB=NA
for (x in 1:length(dTS_SF2$YEAR)) {
  dTS_SF2$dALB[x] = MODIS_SF2[MODIS_SF2[,1]==dTS_SF2$YEAR[x],2]
}

MODIS_SF1=cbind(ALBEDO$FireYear[ALBEDO$SITE=="SF1"],ALBEDO$dALB[ALBEDO$SITE=="SF1"])
dTS_SF1$dALB=NA
for (x in 1:length(dTS_SF1$YEAR)) {
  dTS_SF1$dALB[x] = MODIS_SF1[MODIS_SF1[,1]==dTS_SF1$YEAR[x],2]
}

# bind dataframes
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) # omit sites with standing dead trees

# plot successional winter changes in albedo
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==1,]

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$ALB[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$ALB[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC2 = group_mean
dT_SUCC2$ECO=group_mean2$x
names(dT_SUCC2)=c("SITE","ALB","YEAR","SD","dY","ECO")

# bind dataframes
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF1,dTS_SF2,dTS_SF3) # all sites

dT_SUCC=dT_SUCC[dT_SUCC$SEASON==1,]

# mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$ALB[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$ALB[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC1 = group_mean
dT_SUCC1$ECO=group_mean2$x
names(dT_SUCC1)=c("SITE","ALB","YEAR","SD","dY","ECO")

# fit logarithmic function
fit <- nls(ALB ~ SSasymp(YEAR, yf, y0, log_alpha), data = dT_SUCC1)
r <- range(dT_SUCC1$YEAR)
x = predict(fit)
mat = matrix(ncol = 2, nrow = 200) 
fitted=data.frame(mat)
fitted$YEAR <- seq(1,95,length.out = 200)
fitted$fit <- predict(fit,newdata=fitted)
# get correlation between model fit and data
cor = cor.test(x, dT_SUCC1$ALB, method=c("pearson", "kendall", "spearman"))

# plot successional albedo change
pal <- wes_palette("FantasticFox1", 16, type = "continuous")
fig1 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = YEAR, y = ALB[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = ALB, ymin = ALB - SD,
                                        ymax = ALB + SD), width=0) + # plot site means and errors
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("Years since fire") +
  ylab("Albedo")+
  ggtitle("February to April")+
  annotate("text", x=6, y=0.75, label= "(a)")+
  theme_bw()+
  geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = ALB, shape=ECO),size=3,colour='black', fill='grey')+
  geom_point(dT_SUCC2, mapping = aes(x = YEAR, y = ALB, shape=ECO, fill=SITE),size=3,colour='black')+
  geom_line(fitted, mapping = aes(x = YEAR, y = fit))+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  scale_y_continuous(expand = expand_scale(),limits = c(0,0.8))+
  scale_x_continuous(expand = expand_scale(),limits = c(0,82))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

############ read summer albedo ############
ALBEDO <- read.csv(".../ALB_MODIS_FIRE_SUMMER.csv")

MODIS_RPF=cbind(ALBEDO$FireYear[ALBEDO$SITE=="RPF"],ALBEDO$dALB[ALBEDO$SITE=="RPF"])
dTS_RPF$dALB=NA
for (x in 1:length(dTS_RPF$YEAR)) {
  dTS_RPF$dALB[x] = MODIS_RPF[MODIS_RPF[,1]==dTS_RPF$YEAR[x],2]
}

MODIS_FCR=cbind(ALBEDO$FireYear[ALBEDO$SITE=="FCR"],ALBEDO$dALB[ALBEDO$SITE=="FCR"])
dTS_FCR$dALB=NA
for (x in 1:length(dTS_FCR$YEAR)) {
  dTS_FCR$dALB[x] = MODIS_FCR[MODIS_FCR[,1]==dTS_FCR$YEAR[x],2]
}

MODIS_BN3=cbind(ALBEDO$FireYear[ALBEDO$SITE=="BN3"],ALBEDO$dALB[ALBEDO$SITE=="BN3"])
dTS_BN3$dALB=NA
for (x in 1:length(dTS_BN3$YEAR)) {
  dTS_BN3$dALB[x] = MODIS_BN3[MODIS_BN3[,1]==dTS_BN3$YEAR[x],2]
}

MODIS_BN2=cbind(ALBEDO$FireYear[ALBEDO$SITE=="BN2"],ALBEDO$dALB[ALBEDO$SITE=="BN2"])
dTS_BN2$dALB=NA
for (x in 1:length(dTS_BN2$YEAR)) {
  dTS_BN2$dALB[x] = MODIS_BN2[MODIS_BN2[,1]==dTS_BN2$YEAR[x],2]
}

MODIS_NS2=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS2"],ALBEDO$dALB[ALBEDO$SITE=="NS2"])
dTS_NS2$dALB=NA
for (x in 1:length(dTS_NS2$YEAR)) {
  dTS_NS2$dALB[x] = MODIS_NS2[MODIS_NS2[,1]==dTS_NS2$YEAR[x],2]
}

MODIS_NS3=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS3"],ALBEDO$dALB[ALBEDO$SITE=="NS3"])
dTS_NS3$dALB=NA
for (x in 1:length(dTS_NS3$YEAR)) {
  dTS_NS3$dALB[x] = MODIS_NS3[MODIS_NS3[,1]==dTS_NS3$YEAR[x],2]
}

MODIS_NS4=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS4"],ALBEDO$dALB[ALBEDO$SITE=="NS4"])
dTS_NS4$dALB=NA
for (x in 1:length(dTS_NS4$YEAR)) {
  dTS_NS4$dALB[x] = MODIS_NS4[MODIS_NS4[,1]==dTS_NS4$YEAR[x],2]
}

MODIS_NS5=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS5"],ALBEDO$dALB[ALBEDO$SITE=="NS5"])
dTS_NS5$dALB=NA
for (x in 1:length(dTS_NS5$YEAR)) {
  dTS_NS5$dALB[x] = MODIS_NS5[MODIS_NS5[,1]==dTS_NS5$YEAR[x],2]
}

MODIS_NS6=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS6"],ALBEDO$dALB[ALBEDO$SITE=="NS6"])
dTS_NS6$dALB=NA
for (x in 1:length(dTS_NS6$YEAR)) {
  dTS_NS6$dALB[x] = MODIS_NS6[MODIS_NS6[,1]==dTS_NS6$YEAR[x],2]
}

MODIS_NS7=cbind(ALBEDO$FireYear[ALBEDO$SITE=="NS7"],ALBEDO$dALB[ALBEDO$SITE=="NS7"])
dTS_NS7$dALB=NA
for (x in 1:length(dTS_NS7$YEAR)) {
  dTS_NS7$dALB[x] = MODIS_NS7[MODIS_NS7[,1]==dTS_NS7$YEAR[x],2]
}

MODIS_SF3=cbind(ALBEDO$FireYear[ALBEDO$SITE=="SF3"],ALBEDO$dALB[ALBEDO$SITE=="SF3"])
dTS_SF3$dALB=NA
for (x in 1:length(dTS_SF3$YEAR)) {
  dTS_SF3$dALB[x] = MODIS_SF3[MODIS_SF3[,1]==dTS_SF3$YEAR[x],2]
}

MODIS_SF2=cbind(ALBEDO$FireYear[ALBEDO$SITE=="SF2"],ALBEDO$dALB[ALBEDO$SITE=="SF2"])
dTS_SF2$dALB=NA
for (x in 1:length(dTS_SF2$YEAR)) {
  dTS_SF2$dALB[x] = MODIS_SF2[MODIS_SF2[,1]==dTS_SF2$YEAR[x],2]
}

MODIS_SF1=cbind(ALBEDO$FireYear[ALBEDO$SITE=="SF1"],ALBEDO$dALB[ALBEDO$SITE=="SF1"])
dTS_SF1$dALB=NA
for (x in 1:length(dTS_SF1$YEAR)) {
  dTS_SF1$dALB[x] = MODIS_SF1[MODIS_SF1[,1]==dTS_SF1$YEAR[x],2]
}

# bind dataframes
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) # omit sites with standing trees

# plot successional summer changes in albedo
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$ALB[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$ALB[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC2 = group_mean
dT_SUCC2$ECO=group_mean2$x
names(dT_SUCC2)=c("SITE","ALB","YEAR","SD","dY","ECO")

# bind dataframes from individual sites
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF1,dTS_SF2,dTS_SF3) # all sites

# plot successional summer changes in albedo
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]

# mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$ALB[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$ALB[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC1 = group_mean
dT_SUCC1$ECO=group_mean2$x
names(dT_SUCC1)=c("SITE","ALB","YEAR","SD","dY","ECO")

# plot successional albedo change
pal <- wes_palette("FantasticFox1", 16, type = "continuous")
fig2 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = YEAR, y = ALB[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = ALB, ymin = ALB - SD,
                                        ymax = ALB + SD), width=0) + # plot site means and errors
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("Years since fire") +
  ggtitle("July to September")+
  annotate("text", x=6, y=0.75, label= "(b)")+
  theme_bw()+
  geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = ALB, shape=ECO),size=3,colour='black', fill='grey')+
  geom_point(dT_SUCC2, mapping = aes(x = YEAR, y = ALB, shape=ECO, fill=SITE),size=3,colour='black')+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  ylab(element_blank())+
  scale_y_continuous(expand = expand_scale(),limits = c(0,0.8))+
  scale_x_continuous(expand = expand_scale(),limits = c(0,82))+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = "none",colour = "none")+
  labs(shape='Vegetation type') 

plot1 = grid.arrange(fig1, fig2, nrow = 1)
plot1
#ggsave("Fig6_ALB_WILDFIRE.pdf", width =20, height = 10, units = "cm",plot1)

################################################################################
################# plot chronosequence of aerodynamic conductance ###############
################################################################################

library(matrixStats)
library(patchwork)
MODIS <- read.csv(".../LAI_MODIS_FIRE.csv")

RPF$SEASON=NA
RPF$SEASON[RPF$DOY>=sFeAp & RPF$DOY<=eFeAp]=1 
RPF$SEASON[RPF$DOY>=sJuSe & RPF$DOY<=eJuSe]=2 

dTS_RPF = aggregate(Ga ~ SEASON+YEAR, RPF[!is.na(RPF$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004
dTS_RPF$SITE[dTS_RPF$YEAR<=11]=as.factor(1)
dTS_RPF$SITE[dTS_RPF$YEAR>11]=2
dTS_RPF$SITE=as.factor(dTS_RPF$SITE)
dTS_RPF$ECO[dTS_RPF$YEAR<=11]="OSH"
dTS_RPF$ECO[dTS_RPF$YEAR>11]="DBF"
MODIS_RPF=cbind(MODIS$FireYear[MODIS$SITE=="US-RPF"],MODIS$LAI[MODIS$SITE=="US-RPF"])
dTS_RPF$LAI=NA
for (x in 1:length(dTS_RPF$YEAR)) {
  dTS_RPF$LAI[x] = MODIS_RPF[MODIS_RPF[,1]==dTS_RPF$YEAR[x],2]
}

FCR$SEASON=NA
FCR$SEASON[FCR$DOY>=sFeAp & FCR$DOY<=eFeAp]=1
FCR$SEASON[FCR$DOY>=sJuSe & FCR$DOY<=eJuSe]=2

dTS_FCR = aggregate(Ga ~ SEASON+YEAR, FCR[!is.na(FCR$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(3)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"
MODIS_FCR=cbind(MODIS$FireYear[MODIS$SITE=="US-FCR"],MODIS$LAI[MODIS$SITE=="US-FCR"])
dTS_FCR$LAI=NA
for (x in 1:length(dTS_FCR$YEAR)) {
  dTS_FCR$LAI[x] = MODIS_FCR[MODIS_FCR[,1]==dTS_FCR$YEAR[x],2]
}

BN3$SEASON=NA
BN3$SEASON[BN3$DOY>=sFeAp & BN3$DOY<=eFeAp]=1
BN3$SEASON[BN3$DOY>=sJuSe & BN3$DOY<=eJuSe]=2

dTS_BN3 = aggregate(Ga ~ SEASON+YEAR, BN3[!is.na(BN3$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(4)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"
MODIS_BN3=cbind(MODIS$FireYear[MODIS$SITE=="US-BN3"],MODIS$LAI[MODIS$SITE=="US-BN3"])
dTS_BN3$LAI=NA
for (x in 1:length(dTS_BN3$YEAR)) {
  dTS_BN3$LAI[x] = MODIS_BN3[MODIS_BN3[,1]==dTS_BN3$YEAR[x],2]
}

BN2$SEASON=NA
BN2$SEASON[BN2$DOY>=sFeAp & BN2$DOY<=eFeAp]=1
BN2$SEASON[BN2$DOY>=sJuSe & BN2$DOY<=eJuSe]=2

dTS_BN2 = aggregate(Ga ~ SEASON+YEAR, BN2[!is.na(BN2$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN2$SEASON=as.factor(dTS_BN2$SEASON)
dTS_BN2$SITE=as.factor(5)
dTS_BN2$YEAR=dTS_BN2$YEAR-1987
dTS_BN2$ECO="DBF"
MODIS_BN2=cbind(MODIS$FireYear[MODIS$SITE=="US-BN2"],MODIS$LAI[MODIS$SITE=="US-BN2"])
dTS_BN2$LAI=NA
for (x in 1:length(dTS_BN2$YEAR)) {
  dTS_BN2$LAI[x] = MODIS_BN2[MODIS_BN2[,1]==dTS_BN2$YEAR[x],2]
}

NS2$SEASON=NA
NS2$SEASON[NS2$DOY>=sFeAp & NS2$DOY<=eFeAp]=1
NS2$SEASON[NS2$DOY>=sJuSe & NS2$DOY<=eJuSe]=2

dTS_NS2 = aggregate(Ga ~ SEASON+YEAR, NS2[!is.na(NS2$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(6)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"
MODIS_NS2=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS2"],MODIS$LAI[MODIS$SITE=="CA-NS2"])
dTS_NS2$LAI=NA
for (x in 1:length(dTS_NS2$YEAR)) {
  dTS_NS2$LAI[x] = MODIS_NS2[MODIS_NS2[,1]==dTS_NS2$YEAR[x],2]
}

NS3$SEASON=NA
NS3$SEASON[NS3$DOY>=sFeAp & NS3$DOY<=eFeAp]=1
NS3$SEASON[NS3$DOY>=sJuSe & NS3$DOY<=eJuSe]=2

dTS_NS3 = aggregate(Ga ~ SEASON+YEAR, NS3[!is.na(NS3$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(7)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"
MODIS_NS3=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS3"],MODIS$LAI[MODIS$SITE=="CA-NS3"])
dTS_NS3$LAI=NA
for (x in 1:length(dTS_NS3$YEAR)) {
  dTS_NS3$LAI[x] = MODIS_NS3[MODIS_NS3[,1]==dTS_NS3$YEAR[x],2]
}

NS4$SEASON=NA
NS4$SEASON[NS4$DOY>=sFeAp & NS4$DOY<=eFeAp]=1
NS4$SEASON[NS4$DOY>=sJuSe & NS4$DOY<=eJuSe]=2

dTS_NS4 = aggregate(Ga ~ SEASON+YEAR, NS4[!is.na(NS4$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(8)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"
MODIS_NS4=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS4"],MODIS$LAI[MODIS$SITE=="CA-NS4"])
dTS_NS4$LAI=NA
for (x in 1:length(dTS_NS4$YEAR)) {
  dTS_NS4$LAI[x] = MODIS_NS4[MODIS_NS4[,1]==dTS_NS4$YEAR[x],2]
}

NS5$SEASON=NA
NS5$SEASON[NS5$DOY>=sFeAp & NS5$DOY<=eFeAp]=1
NS5$SEASON[NS5$DOY>=sJuSe & NS5$DOY<=eJuSe]=2

dTS_NS5 = aggregate(Ga ~ SEASON+YEAR, NS5[!is.na(NS5$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(9)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"
MODIS_NS5=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS5"],MODIS$LAI[MODIS$SITE=="CA-NS5"])
dTS_NS5$LAI=NA
for (x in 1:length(dTS_NS5$YEAR)) {
  dTS_NS5$LAI[x] = MODIS_NS5[MODIS_NS5[,1]==dTS_NS5$YEAR[x],2]
}

NS6$SEASON=NA
NS6$SEASON[NS6$DOY>=sFeAp & NS6$DOY<=eFeAp]=1
NS6$SEASON[NS6$DOY>=sJuSe & NS6$DOY<=eJuSe]=2

dTS_NS6 = aggregate(Ga ~ SEASON+YEAR, NS6[!is.na(NS6$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(10)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"
MODIS_NS6=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS6"],MODIS$LAI[MODIS$SITE=="CA-NS6"])
dTS_NS6$LAI=NA
for (x in 1:length(dTS_NS6$YEAR)) {
  dTS_NS6$LAI[x] = MODIS_NS6[MODIS_NS6[,1]==dTS_NS6$YEAR[x],2]
}

NS7$SEASON=NA
NS7$SEASON[NS7$DOY>=sFeAp & NS7$DOY<=eFeAp]=1
NS7$SEASON[NS7$DOY>=sJuSe & NS7$DOY<=eJuSe]=2

dTS_NS7 = aggregate(Ga ~ SEASON+YEAR, NS7[!is.na(NS7$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(11)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"
MODIS_NS7=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS7"],MODIS$LAI[MODIS$SITE=="CA-NS7"])
dTS_NS7$LAI=NA
for (x in 1:length(dTS_NS7$YEAR)) {
  dTS_NS7$LAI[x] = MODIS_NS7[MODIS_NS7[,1]==dTS_NS7$YEAR[x],2]
}

SF3$SEASON=NA
SF3$SEASON[SF3$DOY>=sFeAp & SF3$DOY<=eFeAp]=1
SF3$SEASON[SF3$DOY>=sJuSe & SF3$DOY<=eJuSe]=2

dTS_SF3 = aggregate(Ga ~ SEASON+YEAR, SF3[!is.na(SF3$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(12)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"
MODIS_SF3=cbind(MODIS$FireYear[MODIS$SITE=="CA-SF3"],MODIS$LAI[MODIS$SITE=="CA-SF3"])
dTS_SF3$LAI=NA
for (x in 1:length(dTS_SF3$YEAR)) {
  dTS_SF3$LAI[x] = MODIS_SF3[MODIS_SF3[,1]==dTS_SF3$YEAR[x],2]
}

SF2$SEASON=NA
SF2$SEASON[SF2$DOY>=sFeAp & SF2$DOY<=eFeAp]=1
SF2$SEASON[SF2$DOY>=sJuSe & SF2$DOY<=eJuSe]=2

dTS_SF2 = aggregate(Ga ~ SEASON+YEAR, SF2[!is.na(SF2$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(13)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"
MODIS_SF2=cbind(MODIS$FireYear[MODIS$SITE=="CA-SF2"],MODIS$LAI[MODIS$SITE=="CA-SF2"])
dTS_SF2$LAI=NA
for (x in 1:length(dTS_SF2$YEAR)) {
  dTS_SF2$LAI[x] = MODIS_SF2[MODIS_SF2[,1]==dTS_SF2$YEAR[x],2]
}

SF1$SEASON=NA
SF1$SEASON[SF1$DOY>=sFeAp & SF1$DOY<=eFeAp]=1
SF1$SEASON[SF1$DOY>=sJuSe & SF1$DOY<=eJuSe]=2

dTS_SF1 = aggregate(Ga ~ SEASON+YEAR, SF1[!is.na(SF1$Ga),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(14)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"
MODIS_SF1=cbind(MODIS$FireYear[MODIS$SITE=="CA-SF1"],MODIS$LAI[MODIS$SITE=="CA-SF1"])
dTS_SF1$LAI=NA
for (x in 1:length(dTS_SF1$YEAR)) {
  dTS_SF1$LAI[x] = MODIS_SF1[MODIS_SF1[,1]==dTS_SF1$YEAR[x],2]
}

dTS_RPF=dTS_RPF[!is.na(dTS_RPF$Ga[,2]),]
dTS_FCR=dTS_FCR[!is.na(dTS_FCR$Ga[,2]),]
dTS_BN3=dTS_BN3[!is.na(dTS_BN3$Ga[,2]),]
dTS_BN2=dTS_BN2[!is.na(dTS_BN2$Ga[,2]),]
dTS_NS2=dTS_NS2[!is.na(dTS_NS2$Ga[,2]),]
dTS_NS3=dTS_NS3[!is.na(dTS_NS3$Ga[,2]),]
dTS_NS4=dTS_NS4[!is.na(dTS_NS4$Ga[,2]),]
dTS_NS5=dTS_NS5[!is.na(dTS_NS5$Ga[,2]),]
dTS_NS6=dTS_NS6[!is.na(dTS_NS6$Ga[,2]),]
dTS_NS7=dTS_NS7[!is.na(dTS_NS7$Ga[,2]),]
dTS_SF3=dTS_SF3[!is.na(dTS_SF3$Ga[,2]),]
dTS_SF2=dTS_SF2[!is.na(dTS_SF2$Ga[,2]),]
dTS_SF1=dTS_SF1[!is.na(dTS_SF1$Ga[,2]),]

# bind dataframes from individual sites
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) # omit sites with standing dead trees

# plot successional changes in winter aerodynamic conductance 
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==1,]

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$Ga[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$Ga[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC2 = group_mean
dT_SUCC2$ECO=group_mean2$x
names(dT_SUCC2)=c("SITE","Ga","YEAR","SD","dY","ECO")

dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF1,dTS_SF2,dTS_SF3) # all sites

# plot successional changes in winter aerodynamic conductance
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==1,]

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$Ga[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$Ga[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC1 = group_mean
dT_SUCC1$ECO=group_mean2$x
names(dT_SUCC1)=c("SITE","Ga","YEAR","SD","dY","ECO")

# fit logarithmic function
fit <- nls(Ga ~ SSasymp(YEAR, yf, y0, log_alpha), data = dT_SUCC2)
r <- range(dT_SUCC2$YEAR)
x = predict(fit)
mat = matrix(ncol = 2, nrow = 200) 
fitted=data.frame(mat)
fitted$YEAR <- seq(1,95,length.out = 200)
fitted$fit <- predict(fit,newdata=fitted)
# get correlation between model fit and data
cor = cor.test(x, dT_SUCC2$Ga, method=c("pearson", "kendall", "spearman"))

# plot successional change
pal <- wes_palette("FantasticFox1", 16, type = "continuous")
fig1 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = YEAR, y = Ga[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = Ga, ymin = Ga - SD,
                                        ymax = Ga + SD), width=0) + # plot site means and errors
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("Years since fire") +
  ylab(expression(paste("Aerodynamic conductance [m ", s^-1, "]")))+
  ggtitle("February to April")+
  theme_bw()+
  geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = Ga, shape=ECO),size=3,colour='black', fill='grey')+
  geom_point(dT_SUCC2, mapping = aes(x = YEAR, y = Ga, shape=ECO, fill=SITE),size=3,colour='black')+
  #geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = Ga),size=3,shape=21,fill=NA,colour='black')+
  geom_line(fitted, mapping = aes(x = YEAR, y = fit))+
  scale_shape_manual(values=c(21, 22, 23))+
  annotate("text", x=3, y=0.073, label= "(a)")+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  scale_y_continuous(expand = expand_scale(),limits = c(0.01,0.08))+
  scale_x_continuous(expand = expand_scale(),limits = c(0,82))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

############# plot wintertime successional changes #############
# bind dataframes from individual sites
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) # omit sites with standing dead trees

# plot successional changes in summer aerodynamioc conductance
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$Ga[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$Ga[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC2 = group_mean
dT_SUCC2$ECO=group_mean2$x
names(dT_SUCC2)=c("SITE","Ga","YEAR","SD","dY","ECO")

# bind dataframes from individual sites
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF1,dTS_SF2,dTS_SF3) # all sites

# plot successional changes in summer aerodynamic conductance
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]

# calculate mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$Ga[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# calculate mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# calculate mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# calculate standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$Ga[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC1 = group_mean
dT_SUCC1$ECO=group_mean2$x
names(dT_SUCC1)=c("SITE","Ga","YEAR","SD","dY","ECO")

# fit logarithmic function
fit <- nls(Ga ~ SSasymp(YEAR, yf, y0, log_alpha), data = dT_SUCC2)
r <- range(dT_SUCC2$YEAR)
x = predict(fit)
mat = matrix(ncol = 2, nrow = 200) 
fitted=data.frame(mat)
fitted$YEAR <- seq(1,95,length.out = 200)
fitted$fit <- predict(fit,newdata=fitted)
# get correlation between model fit and data
cor = cor.test(x, dT_SUCC2$Ga, method=c("pearson", "kendall", "spearman"))

# plot successional change
pal <- wes_palette("FantasticFox1", 16, type = "continuous")
fig2 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = YEAR, y = Ga[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = Ga, ymin = Ga - SD,
                                        ymax = Ga + SD), width=0) + # plot site means and errors
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("Years since fire") +
  ylab(element_blank())+
  ggtitle("July to September")+
  theme_bw()+
  geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = Ga, shape=ECO),size=3,colour='black', fill='grey')+
  geom_point(dT_SUCC2, mapping = aes(x = YEAR, y = Ga, shape=ECO, fill=SITE),size=3,colour='black')+
  geom_line(fitted, mapping = aes(x = YEAR, y = fit))+
  scale_shape_manual(values=c(21, 22, 23))+
  annotate("text", x=3, y=0.073, label= "(b)")+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  scale_y_continuous(expand = expand_scale(),limits = c(0.01,0.08))+
  scale_x_continuous(expand = expand_scale(),limits = c(0,82))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = c(.95, .4),
        legend.justification = c("right", "top"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = "none",colour = "none")+
  labs(shape='Vegetation type') 

plot1 = grid.arrange(fig1, fig2, nrow = 1)
plot1
#ggsave("Fig7_Ga_Wildfire.pdf", width =20, height = 10, units = "cm",plot1)

################################################################################
################# plot chronosequence of ecosystem properties ##################
################################################################################
################################################################################
########################### evaporative fraction ###############################
################################################################################

library(matrixStats)
library(patchwork)
MODIS <- read.csv(".../LAI_MODIS_FIRE.csv")

RPF$SEASON=NA
RPF$SEASON[RPF$DOY>=sFeAp & RPF$DOY<=eFeAp]=1 
RPF$SEASON[RPF$DOY>=sJuSe & RPF$DOY<=eJuSe]=2 

dTS_RPF = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, RPF[!is.na(RPF$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004
dTS_RPF$SITE[dTS_RPF$YEAR<=11]=as.factor(1)
dTS_RPF$SITE[dTS_RPF$YEAR>11]=2
dTS_RPF$SITE=as.factor(dTS_RPF$SITE)
dTS_RPF$ECO[dTS_RPF$YEAR<=11]="OSH"
dTS_RPF$ECO[dTS_RPF$YEAR>11]="DBF"

FCR$SEASON=NA
FCR$SEASON[FCR$DOY>=sFeAp & FCR$DOY<=eFeAp]=1
FCR$SEASON[FCR$DOY>=sJuSe & FCR$DOY<=eJuSe]=2

dTS_FCR = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, FCR[!is.na(FCR$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(3)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"

BN3$SEASON=NA
BN3$SEASON[BN3$DOY>=sFeAp & BN3$DOY<=eFeAp]=1
BN3$SEASON[BN3$DOY>=sJuSe & BN3$DOY<=eJuSe]=2

dTS_BN3 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, BN3[!is.na(BN3$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(4)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"

BN2$SEASON=NA
BN2$SEASON[BN2$DOY>=sFeAp & BN2$DOY<=eFeAp]=1
BN2$SEASON[BN2$DOY>=sJuSe & BN2$DOY<=eJuSe]=2

dTS_BN2 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, BN2[!is.na(BN2$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN2$SEASON=as.factor(dTS_BN2$SEASON)
dTS_BN2$SITE=as.factor(5)
dTS_BN2$YEAR=dTS_BN2$YEAR-1987
dTS_BN2$ECO="DBF"

NS1$SEASON=NA
NS1$SEASON[NS1$DOY>=sFeAp & NS1$DOY<=eFeAp]=1
NS1$SEASON[NS1$DOY>=sJuSe & NS1$DOY<=eJuSe]=2

dTS_NS1 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, NS1[!is.na(NS1$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS1$SEASON=as.factor(dTS_NS1$SEASON)
dTS_NS1$SITE=as.factor(6)
dTS_NS1$YEAR=dTS_NS1$YEAR-1850
dTS_NS1$ECO="ENF"

NS2$SEASON=NA
NS2$SEASON[NS2$DOY>=sFeAp & NS2$DOY<=eFeAp]=1
NS2$SEASON[NS2$DOY>=sJuSe & NS2$DOY<=eJuSe]=2

dTS_NS2 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, NS2[!is.na(NS2$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(7)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"

NS3$SEASON=NA
NS3$SEASON[NS3$DOY>=sFeAp & NS3$DOY<=eFeAp]=1
NS3$SEASON[NS3$DOY>=sJuSe & NS3$DOY<=eJuSe]=2

dTS_NS3 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, NS3[!is.na(NS3$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(8)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"

NS4$SEASON=NA
NS4$SEASON[NS4$DOY>=sFeAp & NS4$DOY<=eFeAp]=1
NS4$SEASON[NS4$DOY>=sJuSe & NS4$DOY<=eJuSe]=2

dTS_NS4 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, NS4[!is.na(NS4$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(9)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"

NS5$SEASON=NA
NS5$SEASON[NS5$DOY>=sFeAp & NS5$DOY<=eFeAp]=1
NS5$SEASON[NS5$DOY>=sJuSe & NS5$DOY<=eJuSe]=2

dTS_NS5 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, NS5[!is.na(NS5$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(10)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"

NS6$SEASON=NA
NS6$SEASON[NS6$DOY>=sFeAp & NS6$DOY<=eFeAp]=1
NS6$SEASON[NS6$DOY>=sJuSe & NS6$DOY<=eJuSe]=2

dTS_NS6 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, NS6[!is.na(NS6$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(11)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"

NS7$SEASON=NA
NS7$SEASON[NS7$DOY>=sFeAp & NS7$DOY<=eFeAp]=1
NS7$SEASON[NS7$DOY>=sJuSe & NS7$DOY<=eJuSe]=2

dTS_NS7 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, NS7[!is.na(NS7$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(12)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"

SF3$SEASON=NA
SF3$SEASON[SF3$DOY>=sFeAp & SF3$DOY<=eFeAp]=1
SF3$SEASON[SF3$DOY>=sJuSe & SF3$DOY<=eJuSe]=2

dTS_SF3 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, SF3[!is.na(SF3$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(13)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"

SF2$SEASON=NA
SF2$SEASON[SF2$DOY>=sFeAp & SF2$DOY<=eFeAp]=1
SF2$SEASON[SF2$DOY>=sJuSe & SF2$DOY<=eJuSe]=2

dTS_SF2 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, SF2[!is.na(SF2$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(14)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"

SF1$SEASON=NA
SF1$SEASON[SF1$DOY>=sFeAp & SF1$DOY<=eFeAp]=1
SF1$SEASON[SF1$DOY>=sJuSe & SF1$DOY<=eJuSe]=2

dTS_SF1 = aggregate(cbind(EF,Ga) ~ SEASON+YEAR, SF1[!is.na(SF1$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(15)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"

MODIS_RPF=cbind(MODIS$FireYear[MODIS$SITE=="US-RPF"],MODIS$LAI[MODIS$SITE=="US-RPF"])
dTS_RPF$LAI=NA
for (x in 1:length(dTS_RPF$YEAR)) {
  dTS_RPF$LAI[x] = MODIS_RPF[MODIS_RPF[,1]==dTS_RPF$YEAR[x],2]
}

MODIS_FCR=cbind(MODIS$FireYear[MODIS$SITE=="US-FCR"],MODIS$LAI[MODIS$SITE=="US-FCR"])
dTS_FCR$LAI=NA
for (x in 1:length(dTS_FCR$YEAR)) {
  dTS_FCR$LAI[x] = MODIS_FCR[MODIS_FCR[,1]==dTS_FCR$YEAR[x],2]
}

MODIS_BN3=cbind(MODIS$FireYear[MODIS$SITE=="US-BN3"],MODIS$LAI[MODIS$SITE=="US-BN3"])
dTS_BN3$LAI=NA
for (x in 1:length(dTS_BN3$YEAR)) {
  dTS_BN3$LAI[x] = MODIS_BN3[MODIS_BN3[,1]==dTS_BN3$YEAR[x],2]
}

MODIS_BN2=cbind(MODIS$FireYear[MODIS$SITE=="US-BN2"],MODIS$LAI[MODIS$SITE=="US-BN2"])
dTS_BN2$LAI=NA
for (x in 1:length(dTS_BN2$YEAR)) {
  dTS_BN2$LAI[x] = MODIS_BN2[MODIS_BN2[,1]==dTS_BN2$YEAR[x],2]
}

MODIS_NS2=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS2"],MODIS$LAI[MODIS$SITE=="CA-NS2"])
dTS_NS2$LAI=NA
for (x in 1:length(dTS_NS2$YEAR)) {
  dTS_NS2$LAI[x] = MODIS_NS2[MODIS_NS2[,1]==dTS_NS2$YEAR[x],2]
}

MODIS_NS3=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS3"],MODIS$LAI[MODIS$SITE=="CA-NS3"])
dTS_NS3$LAI=NA
for (x in 1:length(dTS_NS3$YEAR)) {
  dTS_NS3$LAI[x] = MODIS_NS3[MODIS_NS3[,1]==dTS_NS3$YEAR[x],2]
}

MODIS_NS4=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS4"],MODIS$LAI[MODIS$SITE=="CA-NS4"])
dTS_NS4$LAI=NA
for (x in 1:length(dTS_NS4$YEAR)) {
  dTS_NS4$LAI[x] = MODIS_NS4[MODIS_NS4[,1]==dTS_NS4$YEAR[x],2]
}

MODIS_NS5=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS5"],MODIS$LAI[MODIS$SITE=="CA-NS5"])
dTS_NS5$LAI=NA
for (x in 1:length(dTS_NS5$YEAR)) {
  dTS_NS5$LAI[x] = MODIS_NS5[MODIS_NS5[,1]==dTS_NS5$YEAR[x],2]
}

MODIS_NS6=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS6"],MODIS$LAI[MODIS$SITE=="CA-NS6"])
dTS_NS6$LAI=NA
for (x in 1:length(dTS_NS6$YEAR)) {
  dTS_NS6$LAI[x] = MODIS_NS6[MODIS_NS6[,1]==dTS_NS6$YEAR[x],2]
}

MODIS_NS7=cbind(MODIS$FireYear[MODIS$SITE=="CA-NS7"],MODIS$LAI[MODIS$SITE=="CA-NS7"])
dTS_NS7$LAI=NA
for (x in 1:length(dTS_NS7$YEAR)) {
  dTS_NS7$LAI[x] = MODIS_NS7[MODIS_NS7[,1]==dTS_NS7$YEAR[x],2]
}

MODIS_SF3=cbind(MODIS$FireYear[MODIS$SITE=="CA-SF3"],MODIS$LAI[MODIS$SITE=="CA-SF3"])
dTS_SF3$LAI=NA
for (x in 1:length(dTS_SF3$YEAR)) {
  dTS_SF3$LAI[x] = MODIS_SF3[MODIS_SF3[,1]==dTS_SF3$YEAR[x],2]
}

MODIS_SF2=cbind(MODIS$FireYear[MODIS$SITE=="CA-SF2"],MODIS$LAI[MODIS$SITE=="CA-SF2"])
dTS_SF2$LAI=NA
for (x in 1:length(dTS_SF2$YEAR)) {
  dTS_SF2$LAI[x] = MODIS_SF2[MODIS_SF2[,1]==dTS_SF2$YEAR[x],2]
}

MODIS_SF1=cbind(MODIS$FireYear[MODIS$SITE=="CA-SF1"],MODIS$LAI[MODIS$SITE=="CA-SF1"])
dTS_SF1$LAI=NA
for (x in 1:length(dTS_SF1$YEAR)) {
  dTS_SF1$LAI[x] = MODIS_SF1[MODIS_SF1[,1]==dTS_SF1$YEAR[x],2]
}

dTS_RPF=dTS_RPF[!is.na(dTS_RPF$EF[,2]),]
dTS_FCR=dTS_FCR[!is.na(dTS_FCR$EF[,2]),]
dTS_BN3=dTS_BN3[!is.na(dTS_BN3$EF[,2]),]
dTS_BN2=dTS_BN2[!is.na(dTS_BN2$EF[,2]),]
dTS_NS2=dTS_NS2[!is.na(dTS_NS2$EF[,2]),]
dTS_NS3=dTS_NS3[!is.na(dTS_NS3$EF[,2]),]
dTS_NS4=dTS_NS4[!is.na(dTS_NS4$EF[,2]),]
dTS_NS5=dTS_NS5[!is.na(dTS_NS5$EF[,2]),]
dTS_NS6=dTS_NS6[!is.na(dTS_NS6$EF[,2]),]
dTS_NS7=dTS_NS7[!is.na(dTS_NS7$EF[,2]),]
dTS_SF3=dTS_SF3[!is.na(dTS_SF3$EF[,2]),]
dTS_SF2=dTS_SF2[!is.na(dTS_SF2$EF[,2]),]
dTS_SF1=dTS_SF1[!is.na(dTS_SF1$EF[,2]),]


# bind dataframes from individual sites
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) # omit sites with standing dead trees

# plot successional changes in summer evaporative fraction
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]
dT_SUCC=dT_SUCC[dT_SUCC$YEAR<80,]

# mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$EF[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$EF[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC2 = group_mean
dT_SUCC2$ECO=group_mean2$x
names(dT_SUCC2)=c("SITE","EF","YEAR","SD","dY","ECO")

dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF1,dTS_SF2,dTS_SF3) # all sites

dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]
dT_SUCC=dT_SUCC[dT_SUCC$YEAR<80,]

# mean temp difference by site
group_mean<- aggregate(x= dT_SUCC$EF[,1],
                       # Specify group indicator
                       by = list(dT_SUCC$SITE),      
                       # Specify function (i.e. mean)
                       FUN = mean)

# mean age by site
group_mean1<- aggregate(x= dT_SUCC$YEAR,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = mean)

# mean age by site
group_mean2<- aggregate(x= dT_SUCC$ECO,
                        # Specify group indicator
                        by = list(dT_SUCC$SITE),      
                        # Specify function (i.e. mean)
                        FUN = min)

# standard deviation in temp difference by site
group_sd<- aggregate(x= dT_SUCC$EF[,1],
                     # Specify group indicator
                     by = list(dT_SUCC$SITE),      
                     # Specify function (i.e. mean)
                     FUN = sd)

# range of years by site
group_y<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = min)

group_z<- aggregate(x= dT_SUCC$YEAR,
                    # Specify group indicator
                    by = list(dT_SUCC$SITE),      
                    # Specify function (i.e. mean)
                    FUN = max)

group_mean$YEAR = group_mean1$x
group_mean$SD = group_sd$x
group_mean$dY=(group_z$x-group_y$x)/2
dT_SUCC1 = group_mean
dT_SUCC1$ECO=group_mean2$x
names(dT_SUCC1)=c("SITE","EF","YEAR","SD","dY",'ECO')

pal <- wes_palette("FantasticFox1", 16, type = "continuous")
fig1 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = LAI, y = EF[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = EF, ymin = EF - SD,
                                        ymax = EF + SD), width=0)

# plot successional change
pal <- wes_palette("FantasticFox1", 16, type = "continuous")
fig1 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = YEAR, y = EF[,1], colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 3/10)+ # all points
  geom_errorbar(dT_SUCC1, mapping = aes(x = YEAR, y = EF, ymin = EF - SD,
                                        ymax = EF + SD), width=0) + # plot site means and errors
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("Years since fire") +
  ylab("Evaporative fraction")+
  ggtitle("July to September")+
  theme_bw()+
  geom_point(dT_SUCC1, mapping = aes(x = YEAR, y = EF, shape=ECO),size=3,colour='black', fill='grey')+
  geom_point(dT_SUCC2, mapping = aes(x = YEAR, y = EF, shape=ECO, fill=SITE),size=3,colour='black')+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  scale_x_continuous(expand = expand_scale(),limits = c(0,82))+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = "none",colour = "none")+
  labs(shape='Vegetation type') 

fig1
#ggsave("Fig8_EF_Wildfire.pdf", width =12.5, height = 10, units = "cm")

################################################################################
######## Partition the variation in air-surface temperature difference ########
################################################################################

library(vegan)
library(eulerr)
library(VennDiagram)
library(matrixStats)
library(patchwork)
library(matrixStats)
library(gridExtra)
library(ggpubr)
############

dTS_RPF = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, RPF[!is.na(FCR$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004
dTS_RPF$SITE=as.factor(1)
dTS_RPF$ECO="DBF"
dTS_RPF=dTS_RPF[dTS_RPF$SEASON==2,] # choose summer data
dTS_RPF = na.omit(dTS_RPF)

dTS_FCR = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, FCR[!is.na(FCR$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(2)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"
dTS_FCR=dTS_FCR[dTS_FCR$SEASON==2,]

dTS_BN3 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, BN3[!is.na(BN3$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(3)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"
dTS_BN3=dTS_BN3[dTS_BN3$SEASON==2,]

dTS_BN2 = aggregate(cbind(EF,Ga,dT10) ~ SEASON+YEAR, BN2[!is.na(BN2$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN2$SEASON=as.factor(dTS_BN2$SEASON)
dTS_BN2$SITE=as.factor(4)
dTS_BN2$YEAR=dTS_BN2$YEAR-1987
dTS_BN2$ECO="DBF"
dTS_BN2=dTS_BN2[dTS_BN2$SEASON==2,]

dTS_NS2 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, NS2[!is.na(NS2$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(7)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"
dTS_NS2=dTS_NS2[dTS_NS2$SEASON==2,]

dTS_NS3 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, NS3[!is.na(NS3$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(8)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"
dTS_NS3=dTS_NS3[dTS_NS3$SEASON==2,]

dTS_NS4 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, NS4[!is.na(NS4$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(9)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"
dTS_NS4=dTS_NS4[dTS_NS4$SEASON==2,]

dTS_NS5 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, NS5[!is.na(NS5$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(10)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"
dTS_NS5=dTS_NS5[dTS_NS5$SEASON==2,]

dTS_NS6 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, NS6[!is.na(NS6$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(11)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"
dTS_NS6=dTS_NS6[dTS_NS6$SEASON==2,]

dTS_NS7 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, NS7[!is.na(NS7$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(12)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"
dTS_NS7=dTS_NS7[dTS_NS7$SEASON==2,]

dTS_SF3 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, SF3[!is.na(SF3$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(13)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"
dTS_SF3=dTS_SF3[dTS_SF3$SEASON==2,]

dTS_SF2 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, SF2[!is.na(SF2$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(14)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"
dTS_SF2=dTS_SF2[dTS_SF2$SEASON==2,]

dTS_SF1 = aggregate(cbind(EF,Ga,ALB,dT10) ~ SEASON+YEAR, SF1[!is.na(SF1$EF),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(15)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"
dTS_SF1=dTS_SF1[dTS_SF1$SEASON==2,]

########### run variation partitioning ###########
# run it 1000 times to randomly sample individual years
library(wesanderson)
fractions = data.frame(matrix(ncol = 8, nrow = 1000))
nlc <- nls.control(maxiter = 1000)
for (y in 1:1000)
{
  # sample one year per site to avoid bias due to different sample size
  dTS_RPFs=dTS_RPF[sample(1:nrow(dTS_RPF), 3),]
  dTS_FCRs=dTS_FCR[sample(1:nrow(dTS_FCR), 3),]
  dTS_BN3s=dTS_BN3[sample(1:nrow(dTS_BN3), 3),]
  dTS_NS2s=dTS_NS2[sample(1:nrow(dTS_NS2), 3),]
  dTS_NS3s=dTS_NS3[sample(1:nrow(dTS_NS3), 3),]
  dTS_NS4s=dTS_NS4[sample(1:nrow(dTS_NS4), 3),]
  dTS_NS5s=dTS_NS5[sample(1:nrow(dTS_NS5), 3),]
  dTS_NS6s=dTS_NS6[sample(1:nrow(dTS_NS6), 3),]
  dTS_NS7s=dTS_NS7[sample(1:nrow(dTS_NS7), 3),]
  dTS_SF3s=dTS_SF3[sample(1:nrow(dTS_SF3), 3),]
  dTS_SF2s=dTS_SF2[sample(1:nrow(dTS_SF2), 3),]
  dTS_SF1s=dTS_SF1[sample(1:nrow(dTS_SF1), 3),]
  
  dT_SUCC1=rbind(dTS_RPFs,dTS_FCRs,dTS_BN3s,dTS_NS2s,dTS_NS3s,dTS_NS4s,dTS_NS5s,dTS_NS6s,dTS_NS7s,dTS_SF1s,dTS_SF2s,dTS_SF3s) # BN2 has no albedo measurements
  ts = dT_SUCC1$dT10[,1]
  alb = dT_SUCC1$ALB[,1]
  ga=dT_SUCC1$Ga[,1]
  ef=dT_SUCC1$EF[,1]
  ts.part.all <- varpart(ts, alb, ga, ef)
  fractions[y,]=ts.part.all$part$indfract$Adj.R.square
  rm(dT_SUCC1)
  rm(ts.part.all)
}

####### run variation partitioning without subsampling #######

dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF3,dTS_SF2,dTS_SF1)
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]
dT_SUCC <- na.omit(dT_SUCC)
ts = dT_SUCC$dT10[,1]
alb = dT_SUCC$ALB[,1]
ga = dT_SUCC$Ga[,1]
ef=dT_SUCC$EF[,1]
part = cbind(dT_SUCC$Ga[,1],dT_SUCC$EF[,1])
ts.part.all <- varpart(ts, alb, ga, ef)
ts.part.all$part

####### plot variation partitioning results #######
library(ggplot2)
library(ggforce)
alb_fr = median(fractions[,1])
ga_fr = median(fractions[,2])
ef_fr = median(fractions[,3])
alb_ga_fr = median(fractions[,4])
ga_ef_fr = median(fractions[,5])
ef_alb_fr = median(fractions[,6])
ga_ef_alb_fr = median(fractions[,7])
pal <- c("#B40F20","#E2D200", "#46ACC8")
df.venn <- data.frame(x = c(3, 1, 2),y = c(1, 1,2.8),labels = c('Albedo [<0]', 'Aerodyn. cond. [34]',"Evap. fraction [0]"))
p1 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +geom_circle(alpha = .5, size = 2, colour = 'grey',show.legend = FALSE ) +coord_fixed()+
  annotate("text", x = df.venn$x[2]-0.5 , y = df.venn$y[2]-0.2,label=df.venn$labels[2] ,size = 3.5)+annotate("text", x = df.venn$x[1]+0.5, y = df.venn$y[1]-0.2,label=df.venn$labels[1] ,size = 3.5)+
  annotate("text", x = df.venn$x[3], y = df.venn$y[3]+0.2,label=df.venn$labels[3] ,size = 3)+scale_fill_manual(values = pal)+ggtitle("(d)")
pf = p1+annotate("text", x = 1.35 , y =2,label="[5]" ,size = 3.5)+annotate("text", x = 2.02 , y =1.65,label="[4]" ,size = 3.5)+
  annotate("text", x = 2.7 , y =2,label="[0]" ,size = 3.5)+annotate("text", x = 2.02 , y =0.8,label="[20]" ,size = 3.5)+
  annotate("text", x = 2 , y =-0.75,label="Unexplained variance: 36" ,size = 3.5)+theme_void()
pf

####### plot relationship between aerodynamic conductance and dT #######
pd <- position_dodge(0.1) # move them .05 to the left and right
pal <- wes_palette("FantasticFox1", 16, type = "continuous")

p1=ggplot(dT_SUCC, aes(x=Ga[,1], y = dT10[,1], color=SITE, shape =ECO, fill=SITE))+
  geom_errorbar(aes(xmin=Ga[,1]-Ga[,2], xmax=Ga[,1]+Ga[,2]), width=.0, position=pd) +
  geom_errorbar(aes(ymin=dT10[,1]-dT10[,2], ymax=dT10[,1]+dT10[,2]), width=.0, position=pd) +
  geom_point(size=3,colour='black')+ #position=pd, 
  xlab(expression(Aerodynamic ~ conductance ~ "[m" ~ s^-1 ~ "]")) +
  annotate("text", x=0.02, y=5, label= "(c)")+
  ylab(expression(paste(Delta *"temperature " (surface-air) *" [" *~degree*C* "]")))+
  scale_shape_discrete(labels=c('deciduous broafleaf','evergreen needleleaf','open shrubland'))+
  guides(color="none")+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  guides(shape=guide_legend(title="Vegetation type"))+
  scale_color_manual(values = pal)+
  ggtitle("July - September")+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(fill = "none",colour = "none")+
  labs(shape='Vegetation type') 

p1 = p1 +  theme(
  # Hide panel borders and remove grid lines
  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  aspect.ratio=1,
  # Change axis line
  axis.line = element_line(colour = "black"),
  panel.background = element_blank(),
  text = element_text(size = 16)
)
#ggsave("Fig_dT_Ga.pdf", width =15, height = 15, units = "cm")
#ggsave("Fig_dT_Ga.png", width =15, height = 15, units = "cm")
#grid.draw(p);

############ variation partitioning for winter ############
dTS_RPF = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, RPF[!is.na(RPF$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004
dTS_RPF$SITE=as.factor(1)
dTS_RPF$ECO="DBF"

dTS_FCR = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, FCR[!is.na(FCR$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(2)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"

dTS_BN3 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, BN3[!is.na(BN3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(3)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"

dTS_NS2 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, NS2[!is.na(NS2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(7)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"

dTS_NS3 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, NS3[!is.na(NS3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(8)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"

dTS_NS4 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, NS4[!is.na(NS4$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(9)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"

dTS_NS5 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, NS5[!is.na(NS5$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(10)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"

dTS_NS6 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, NS6[!is.na(NS6$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(11)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"

dTS_NS7 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, NS7[!is.na(NS7$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(12)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"

dTS_SF3 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, SF3[!is.na(SF3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(13)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"

dTS_SF2 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, SF2[!is.na(SF2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(14)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"

dTS_SF1 = aggregate(cbind(Ga,ALB,dT10) ~ SEASON+YEAR, SF1[!is.na(SF1$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(15)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"

############ choose winter season ############
dTS_RPF=dTS_RPF[dTS_RPF$SEASON==1,,]
dTS_FCR=dTS_FCR[dTS_FCR$SEASON==1,,]
dTS_BN3=dTS_BN3[dTS_BN3$SEASON==1,,]
dTS_NS2=dTS_NS2[dTS_NS2$SEASON==1,,]
dTS_NS3=dTS_NS3[dTS_NS3$SEASON==1,,]
dTS_NS4=dTS_NS4[dTS_NS4$SEASON==1,,]
dTS_NS5=dTS_NS5[dTS_NS5$SEASON==1,,]
dTS_NS6=dTS_NS6[dTS_NS6$SEASON==1,,]
dTS_NS7=dTS_NS7[dTS_NS7$SEASON==1,,]
dTS_SF3=dTS_SF3[dTS_SF3$SEASON==1,,]
dTS_SF2=dTS_SF2[dTS_SF2$SEASON==1,,]
dTS_SF1=dTS_SF1[dTS_SF1$SEASON==1,,]
dTS_SF3=dTS_SF3[dTS_SF3$SEASON==1,,]

fractions = data.frame(matrix(ncol = 4, nrow = 1000))
nlc <- nls.control(maxiter = 1000)
for (y in 1:1000)
{
  # sample one year per site
  dTS_RPFs=dTS_RPF[sample(1:nrow(dTS_RPF), 2),]
  dTS_FCRs=dTS_FCR[sample(1:nrow(dTS_FCR), 2),]
  dTS_BN3s=dTS_BN3[sample(1:nrow(dTS_BN3), 2),]
  dTS_NS2s=dTS_NS2[sample(1:nrow(dTS_NS2), 2),]
  dTS_NS3s=dTS_NS3[sample(1:nrow(dTS_NS3), 2),]
  dTS_NS4s=dTS_NS4[sample(1:nrow(dTS_NS4), 2, replace = TRUE),]
  dTS_NS5s=dTS_NS5[sample(1:nrow(dTS_NS5), 2),]
  dTS_NS6s=dTS_NS6[sample(1:nrow(dTS_NS6), 2),]
  dTS_NS7s=dTS_NS7[sample(1:nrow(dTS_NS7), 2),]
  dTS_SF3s=dTS_SF3[sample(1:nrow(dTS_SF3), 2),]
  dTS_SF2s=dTS_SF2[sample(1:nrow(dTS_SF2), 2),]
  dTS_SF1s=dTS_SF1[sample(1:nrow(dTS_SF1), 2),]
  
  dT_SUCC1=rbind(dTS_RPFs,dTS_FCRs,dTS_BN3s,dTS_NS2s,dTS_NS3s,dTS_NS4s,dTS_NS5s,dTS_NS6s,dTS_NS7s,dTS_SF1s,dTS_SF2s,dTS_SF3s) 
  ts = dT_SUCC1$dT10[,1]
  alb = dT_SUCC1$ALB[,1]
  ga = dT_SUCC1$Ga[,1]
  part = cbind(dT_SUCC1$Ga[,1],dT_SUCC1$EF[,1])
  ts.part.all <- varpart(ts, alb, ga)
  fractions[y,]=ts.part.all$part$indfract$Adj.R.squared
  rm(dT_SUCC1)
  rm(ts.part.all)
}
alb_fr = median(fractions[,1])
alb_ga_fr = median(fractions[,2])
ga_fr = median(fractions[,3])
unexpl_fr = median(fractions[,4])

# run variation partitioning without subsampling
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_NS7,dTS_SF3,dTS_SF2,dTS_SF1)
dT_SUCC <- na.omit(dT_SUCC)

ts = dT_SUCC$dT10
alb = dT_SUCC$ALB[,1]
ga = dT_SUCC$Ga[,1]
part = cbind(dT_SUCC$Ga[,1],dT_SUCC$EF[,1])
ts.part.all <- varpart(ts, alb, ga)
ts.part.all$part

pal <- c("#B40F20","#E2D200")
df.venn <- data.frame(x = c(3, 1),y = c(1, 1),labels = c('Albedo [<0]', 'Aerodyn. cond. [7]'))
p3 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +geom_circle(alpha = .5, size = 2, colour = 'grey',show.legend = FALSE ) +coord_fixed()+
  annotate("text", x = df.venn$x[2]-0.5 , y = df.venn$y[2],label=df.venn$labels[2] ,size = 3.5)+
  annotate("text", x = df.venn$x[1]+0.5, y = df.venn$y[1],label=df.venn$labels[1] ,size = 3.5)+scale_fill_manual(values = pal)
pg = p3+annotate("text", x = 2 , y =1,label="[21]" ,size = 3.5)+
  annotate("text", x = -0.5 , y =2.3,label="(b)" ,size = 3.5)+
  annotate("text", x = 2 , y =-0.75,label="Unexplained variance: 71" ,size = 3.5)+theme_void()
#ggsave("Fig_VarPartWint.pdf", width =15, height = 10, units = "cm")
#ggsave("Fig_VarPartWint.png", width =15, height = 10, units = "cm")

pal <- wes_palette("FantasticFox1", 16, type = "continuous")
p2=ggplot(dT_SUCC, aes(x=ALB[,1], y = dT10[,1], color=SITE, shape =ECO, fill=SITE))+
  geom_errorbar(aes(xmin=ALB[,1]-ALB[,2], xmax=ALB[,1]+ALB[,2]), width=.0, position=pd) +
  geom_errorbar(aes(ymin=dT10[,1]-dT10[,2], ymax=dT10[,1]+dT10[,2]), width=.0, position=pd) +
  geom_point(size=3,colour='black')+
  xlab("Albedo") +
  annotate("text", x=0.02, y=3.7, label= "(a)")+
  ylab(expression(paste(Delta *"temperature " (surface-air) *" [" *~degree*C* "]")))+
  scale_shape_discrete(labels=c('deciduous broafleaf','evergreen needleleaf','open shrubland'))+
  guides(color="none")+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  guides(shape=guide_legend(title="Vegetation type"))+
  scale_color_manual(values = pal)+
  ggtitle("February - April")+
  theme(legend.position = c(.98, 0.99),
        legend.justification = c("right", "top"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.text=element_text(size=10),legend.title=element_text(size=10))+
  guides(fill = "none",colour = "none")+
  labs(shape='Vegetation type') 
p2 = p2 + theme(
  # Hide panel borders and remove grid lines
  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  aspect.ratio=1,
  # Change axis line
  axis.line = element_line(colour = "black"),
  panel.background = element_blank(),
  text = element_text(size = 16)
)

plot1 = grid.arrange(p2, pg, p1, pf, nrow = 2, ncol = 2)
plot1

#ggsave("Fig10_Partitioning.pdf", width =25, height = 25, units = "cm",plot1)

################################################################################
################ compare satellite versus EC temp differences ##################
################################################################################
library(matrixStats)
library(patchwork)
library(gridExtra)
library(ggpubr)
library(ggpmisc)

MODIS <- read.csv(".../MODIS_EC_WINTER_dTEMP.csv")
# select seasons #
sMaAp = 32
eMaAp = 120
sJuAu = 182
eJuAu = 273

RPF$SEASON=NA
RPF$SEASON[RPF$DOY>=sMaAp & RPF$DOY<=eMaAp]=1 # Jan - Mar
RPF$SEASON[RPF$DOY>=sJuAu & RPF$DOY<=eJuAu]=2 # Jul - Sep

# select clear-sky days (ClearSky == 1) or all days (ClearSky == 1 | 0)
dTS_RPF = aggregate(dT10 ~ SEASON+YEAR, RPF[!is.na(RPF$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004

dTS_RPF$SITE[dTS_RPF$YEAR<=11]=as.factor(1)
dTS_RPF$SITE[dTS_RPF$YEAR>11]=2
dTS_RPF$SITE=as.factor(dTS_RPF$SITE)
dTS_RPF$ECO[dTS_RPF$YEAR<=11]="OSH"
dTS_RPF$ECO[dTS_RPF$YEAR>11]="DBF"
dTS_RPF=dTS_RPF[dTS_RPF$SEASON==1,]
for (y in 1:length(dTS_RPF$YEAR))
{
  dTS_RPF$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_RPF$YEAR[y] & MODIS$SITE=='RPF']
}

FCR$SEASON=NA
FCR$SEASON[FCR$DOY>=sMaAp & FCR$DOY<=eMaAp]=1
FCR$SEASON[FCR$DOY>=sJuAu & FCR$DOY<=eJuAu]=2

dTS_FCR = aggregate(dT10 ~ SEASON+YEAR, FCR[!is.na(FCR$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(3)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"
dTS_FCR=dTS_FCR[dTS_FCR$SEASON==1,]
for (y in 1:length(dTS_FCR$YEAR))
{
  dTS_FCR$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_FCR$YEAR[y] & MODIS$SITE=='FCR']
}

BN3$SEASON=NA
BN3$SEASON[BN3$DOY>=sMaAp & BN3$DOY<=eMaAp]=1
BN3$SEASON[BN3$DOY>=sJuAu & BN3$DOY<=eJuAu]=2

dTS_BN3 = aggregate(dT10 ~ SEASON+YEAR, BN3[!is.na(BN3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(4)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"
dTS_BN3=dTS_BN3[dTS_BN3$SEASON==1,]
for (y in 1:length(dTS_BN3$YEAR))
{
  dTS_BN3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN3$YEAR[y] & MODIS$SITE=='BN3']
}


BN2$SEASON=NA
BN2$SEASON[BN2$DOY>=sMaAp & BN2$DOY<=eMaAp]=1
BN2$SEASON[BN2$DOY>=sJuAu & BN2$DOY<=eJuAu]=2

dTS_BN2 = aggregate(dT10 ~ SEASON+YEAR, BN2[!is.na(BN2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN2$SEASON=as.factor(dTS_BN2$SEASON)
dTS_BN2$SITE=as.factor(5)
dTS_BN2$YEAR=dTS_BN2$YEAR-1987
dTS_BN2$ECO="DBF"
dTS_BN2=dTS_BN2[dTS_BN2$SEASON==1,]
for (y in 1:length(dTS_BN2$YEAR))
{
  dTS_BN2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN2$YEAR[y] & MODIS$SITE=='BN2']
}

NS2$SEASON=NA
NS2$SEASON[NS2$DOY>=sMaAp & NS2$DOY<=eMaAp]=1
NS2$SEASON[NS2$DOY>=sJuAu & NS2$DOY<=eJuAu]=2

dTS_NS2 = aggregate(dT10 ~ SEASON+YEAR, NS2[!is.na(NS2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(6)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"
dTS_NS2=dTS_NS2[dTS_NS2$SEASON==1,]
for (y in 1:length(dTS_NS2$YEAR))
{
  dTS_NS2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS2$YEAR[y] & MODIS$SITE=='NS2']
}

NS3$SEASON=NA
NS3$SEASON[NS3$DOY>=sMaAp & NS3$DOY<=eMaAp]=1
NS3$SEASON[NS3$DOY>=sJuAu & NS3$DOY<=eJuAu]=2

dTS_NS3 = aggregate(dT10 ~ SEASON+YEAR, NS3[!is.na(NS3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(7)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"
dTS_NS3=dTS_NS3[dTS_NS3$SEASON==1,]
for (y in 1:length(dTS_NS3$YEAR))
{
  dTS_NS3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS3$YEAR[y] & MODIS$SITE=='NS3']
}

NS4$SEASON=NA
NS4$SEASON[NS4$DOY>=sMaAp & NS4$DOY<=eMaAp]=1
NS4$SEASON[NS4$DOY>=sJuAu & NS4$DOY<=eJuAu]=2

dTS_NS4 = aggregate(dT10 ~ SEASON+YEAR, NS4[!is.na(NS4$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(8)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"
dTS_NS4=dTS_NS4[dTS_NS4$SEASON==1,]
for (y in 1:length(dTS_NS4$YEAR))
{
  dTS_NS4$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS4$YEAR[y] & MODIS$SITE=='NS4']
}

NS5$SEASON=NA
NS5$SEASON[NS5$DOY>=sMaAp & NS5$DOY<=eMaAp]=1
NS5$SEASON[NS5$DOY>=sJuAu & NS5$DOY<=eJuAu]=2

dTS_NS5 = aggregate(dT10 ~ SEASON+YEAR, NS5[!is.na(NS5$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(9)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"
dTS_NS5=dTS_NS5[dTS_NS5$SEASON==1,]
for (y in 1:length(dTS_NS5$YEAR))
{
  dTS_NS5$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS5$YEAR[y] & MODIS$SITE=='NS5']
}

NS6$SEASON=NA
NS6$SEASON[NS6$DOY>=sMaAp & NS6$DOY<=eMaAp]=1
NS6$SEASON[NS6$DOY>=sJuAu & NS6$DOY<=eJuAu]=2

dTS_NS6 = aggregate(dT10 ~ SEASON+YEAR, NS6[!is.na(NS6$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(10)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"
dTS_NS6=dTS_NS6[dTS_NS6$SEASON==1,]
for (y in 1:length(dTS_NS6$YEAR))
{
  dTS_NS6$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS6$YEAR[y] & MODIS$SITE=='NS6']
}

NS7$SEASON=NA
NS7$SEASON[NS7$DOY>=sMaAp & NS7$DOY<=eMaAp]=1
NS7$SEASON[NS7$DOY>=sJuAu & NS7$DOY<=eJuAu]=2

dTS_NS7 = aggregate(dT10 ~ SEASON+YEAR, NS7[!is.na(NS7$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(11)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"
dTS_NS7=dTS_NS7[dTS_NS7$SEASON==1,]
for (y in 1:length(dTS_NS7$YEAR))
{
  dTS_NS7$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS7$YEAR[y] & MODIS$SITE=='NS7']
}


SF3$SEASON=NA
SF3$SEASON[SF3$DOY>=sMaAp & SF3$DOY<=eMaAp]=1
SF3$SEASON[SF3$DOY>=sJuAu & SF3$DOY<=eJuAu]=2

dTS_SF3 = aggregate(dT10 ~ SEASON+YEAR, SF3[!is.na(SF3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(12)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"
dTS_SF3=dTS_SF3[dTS_SF3$SEASON==1,]
for (y in 1:length(dTS_SF3$YEAR))
{
  dTS_SF3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF3$YEAR[y] & MODIS$SITE=='SF3']
}

SF2$SEASON=NA
SF2$SEASON[SF2$DOY>=sMaAp & SF2$DOY<=eMaAp]=1
SF2$SEASON[SF2$DOY>=sJuAu & SF2$DOY<=eJuAu]=2

dTS_SF2 = aggregate(dT10 ~ SEASON+YEAR, SF2[!is.na(SF2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(13)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"
dTS_SF2=dTS_SF2[dTS_SF2$SEASON==1,]
for (y in 1:length(dTS_SF2$YEAR))
{
  dTS_SF2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF2$YEAR[y] & MODIS$SITE=='SF2']
}

SF1$SEASON=NA
SF1$SEASON[SF1$DOY>=sMaAp & SF1$DOY<=eMaAp]=1
SF1$SEASON[SF1$DOY>=sJuAu & SF1$DOY<=eJuAu]=2

dTS_SF1 = aggregate(dT10 ~ SEASON+YEAR, SF1[!is.na(SF1$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(14)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"
dTS_SF1=dTS_SF1[dTS_SF1$SEASON==1,]
for (y in 1:length(dTS_SF1$YEAR))
{
  dTS_SF1$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF1$YEAR[y] & MODIS$SITE=='SF1']
}

# bind dataframes from individual sites
dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) 

# plot successional changes in air-surface temperature difference
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==1,]

library(wesanderson)
pal <- wes_palette("FantasticFox1", 14, type = "continuous")
fig3 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = dT10[,1], y = dTEMP, colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 8/10)+
  stat_poly_line(data=dT_SUCC,aes(dT10[,1], dTEMP), color = "chocolate3") +
  stat_poly_eq(data=dT_SUCC,aes(dT10[,1], dTEMP)) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  )+
  xlab("\u0394 T (surface - air) [\u00B0C]") +
  ylab("\u0394 LST (fire -control) [\u00B0C]")+
  ggtitle("February to April")+
  annotate("text", x=2.7, y=2, label= "(a)")+
  theme_bw()+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#########summer data##############
MODIS <- read.csv(".../MODIS_EC_SUMMER_dTEMP.csv")
dTS_RPF = aggregate(dT10 ~ SEASON+YEAR, RPF[!is.na(RPF$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_RPF$SEASON=as.factor(dTS_RPF$SEASON)
dTS_RPF$YEAR=dTS_RPF$YEAR-2004
dTS_RPF$SITE[dTS_RPF$YEAR<=11]=as.factor(1)
dTS_RPF$SITE[dTS_RPF$YEAR>11]=2
dTS_RPF$SITE=as.factor(dTS_RPF$SITE)
dTS_RPF$ECO[dTS_RPF$YEAR<=11]="OSH"
dTS_RPF$ECO[dTS_RPF$YEAR>11]="DBF"
dTS_RPF=dTS_RPF[dTS_RPF$SEASON==2,]
for (y in 1:length(dTS_RPF$YEAR))
{
  dTS_RPF$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_RPF$YEAR[y] & MODIS$SITE=='RPF']
}

FCR$SEASON=NA
FCR$SEASON[FCR$DOY>=sMaAp & FCR$DOY<=eMaAp]=1
FCR$SEASON[FCR$DOY>=sJuAu & FCR$DOY<=eJuAu]=2

dTS_FCR = aggregate(dT10 ~ SEASON+YEAR, FCR[!is.na(FCR$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_FCR$SEASON=as.factor(dTS_FCR$SEASON)
dTS_FCR$SITE=as.factor(3)
dTS_FCR$YEAR=dTS_FCR$YEAR-2010
dTS_FCR$ECO="OSH"
dTS_FCR=dTS_FCR[dTS_FCR$SEASON==2,]
for (y in 1:length(dTS_FCR$YEAR))
{
  dTS_FCR$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_FCR$YEAR[y] & MODIS$SITE=='FCR']
}

BN3$SEASON=NA
BN3$SEASON[BN3$DOY>=sMaAp & BN3$DOY<=eMaAp]=1
BN3$SEASON[BN3$DOY>=sJuAu & BN3$DOY<=eJuAu]=2

dTS_BN3 = aggregate(dT10 ~ SEASON+YEAR, BN3[!is.na(BN3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN3$SEASON=as.factor(dTS_BN3$SEASON)
dTS_BN3$SITE=as.factor(4)
dTS_BN3$YEAR=dTS_BN3$YEAR-1999
dTS_BN3$ECO="OSH"
dTS_BN3=dTS_BN3[dTS_BN3$SEASON==2,]
for (y in 1:length(dTS_BN3$YEAR))
{
  dTS_BN3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN3$YEAR[y] & MODIS$SITE=='BN3']
}


BN2$SEASON=NA
BN2$SEASON[BN2$DOY>=sMaAp & BN2$DOY<=eMaAp]=1
BN2$SEASON[BN2$DOY>=sJuAu & BN2$DOY<=eJuAu]=2

dTS_BN2 = aggregate(dT10 ~ SEASON+YEAR, BN2[!is.na(BN2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_BN2$SEASON=as.factor(dTS_BN2$SEASON)
dTS_BN2$SITE=as.factor(5)
dTS_BN2$YEAR=dTS_BN2$YEAR-1987
dTS_BN2$ECO="DBF"
dTS_BN2=dTS_BN2[dTS_BN2$SEASON==2,]
for (y in 1:length(dTS_BN2$YEAR))
{
  dTS_BN2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_BN2$YEAR[y] & MODIS$SITE=='BN2']
}

NS2$SEASON=NA
NS2$SEASON[NS2$DOY>=sMaAp & NS2$DOY<=eMaAp]=1
NS2$SEASON[NS2$DOY>=sJuAu & NS2$DOY<=eJuAu]=2

dTS_NS2 = aggregate(dT10 ~ SEASON+YEAR, NS2[!is.na(NS2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS2$SEASON=as.factor(dTS_NS2$SEASON)
dTS_NS2$SITE=as.factor(6)
dTS_NS2$YEAR=dTS_NS2$YEAR-1930
dTS_NS2$ECO="ENF"
dTS_NS2=dTS_NS2[dTS_NS2$SEASON==2,]
for (y in 1:length(dTS_NS2$YEAR))
{
  dTS_NS2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS2$YEAR[y] & MODIS$SITE=='NS2']
}

NS3$SEASON=NA
NS3$SEASON[NS3$DOY>=sMaAp & NS3$DOY<=eMaAp]=1
NS3$SEASON[NS3$DOY>=sJuAu & NS3$DOY<=eJuAu]=2

dTS_NS3 = aggregate(dT10 ~ SEASON+YEAR, NS3[!is.na(NS3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS3$SEASON=as.factor(dTS_NS3$SEASON)
dTS_NS3$SITE=as.factor(7)
dTS_NS3$YEAR=dTS_NS3$YEAR-1964
dTS_NS3$ECO="ENF"
dTS_NS3=dTS_NS3[dTS_NS3$SEASON==2,]
for (y in 1:length(dTS_NS3$YEAR))
{
  dTS_NS3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS3$YEAR[y] & MODIS$SITE=='NS3']
}

NS4$SEASON=NA
NS4$SEASON[NS4$DOY>=sMaAp & NS4$DOY<=eMaAp]=1
NS4$SEASON[NS4$DOY>=sJuAu & NS4$DOY<=eJuAu]=2

dTS_NS4 = aggregate(dT10 ~ SEASON+YEAR, NS4[!is.na(NS4$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS4$SEASON=as.factor(dTS_NS4$SEASON)
dTS_NS4$SITE=as.factor(8)
dTS_NS4$YEAR=dTS_NS4$YEAR-1964
dTS_NS4$ECO="ENF"
dTS_NS4=dTS_NS4[dTS_NS4$SEASON==2,]
for (y in 1:length(dTS_NS4$YEAR))
{
  dTS_NS4$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS4$YEAR[y] & MODIS$SITE=='NS4']
}

NS5$SEASON=NA
NS5$SEASON[NS5$DOY>=sMaAp & NS5$DOY<=eMaAp]=1
NS5$SEASON[NS5$DOY>=sJuAu & NS5$DOY<=eJuAu]=2

dTS_NS5 = aggregate(dT10 ~ SEASON+YEAR, NS5[!is.na(NS5$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS5$SEASON=as.factor(dTS_NS5$SEASON)
dTS_NS5$SITE=as.factor(9)
dTS_NS5$YEAR=dTS_NS5$YEAR-1981
dTS_NS5$ECO="ENF"
dTS_NS5=dTS_NS5[dTS_NS5$SEASON==2,]
for (y in 1:length(dTS_NS5$YEAR))
{
  dTS_NS5$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS5$YEAR[y] & MODIS$SITE=='NS5']
}

NS6$SEASON=NA
NS6$SEASON[NS6$DOY>=sMaAp & NS6$DOY<=eMaAp]=1
NS6$SEASON[NS6$DOY>=sJuAu & NS6$DOY<=eJuAu]=2

dTS_NS6 = aggregate(dT10 ~ SEASON+YEAR, NS6[!is.na(NS6$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS6$SEASON=as.factor(dTS_NS6$SEASON)
dTS_NS6$SITE=as.factor(10)
dTS_NS6$YEAR=dTS_NS6$YEAR-1989
dTS_NS6$ECO="OSH"
dTS_NS6=dTS_NS6[dTS_NS6$SEASON==2,]
for (y in 1:length(dTS_NS6$YEAR))
{
  dTS_NS6$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS6$YEAR[y] & MODIS$SITE=='NS6']
}

NS7$SEASON=NA
NS7$SEASON[NS7$DOY>=sMaAp & NS7$DOY<=eMaAp]=1
NS7$SEASON[NS7$DOY>=sJuAu & NS7$DOY<=eJuAu]=2

dTS_NS7 = aggregate(dT10 ~ SEASON+YEAR, NS7[!is.na(NS7$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_NS7$SEASON=as.factor(dTS_NS7$SEASON)
dTS_NS7$SITE=as.factor(11)
dTS_NS7$YEAR=dTS_NS7$YEAR-1998
dTS_NS7$ECO="OSH"
dTS_NS7=dTS_NS7[dTS_NS7$SEASON==2,]
for (y in 1:length(dTS_NS7$YEAR))
{
  dTS_NS7$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_NS7$YEAR[y] & MODIS$SITE=='NS7']
}


SF3$SEASON=NA
SF3$SEASON[SF3$DOY>=sMaAp & SF3$DOY<=eMaAp]=1
SF3$SEASON[SF3$DOY>=sJuAu & SF3$DOY<=eJuAu]=2

dTS_SF3 = aggregate(dT10 ~ SEASON+YEAR, SF3[!is.na(SF3$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF3$SEASON=as.factor(dTS_SF3$SEASON)
dTS_SF3$SITE=as.factor(12)
dTS_SF3$YEAR=dTS_SF3$YEAR-1998
dTS_SF3$ECO="OSH"
dTS_SF3=dTS_SF3[dTS_SF3$SEASON==2,]
for (y in 1:length(dTS_SF3$YEAR))
{
  dTS_SF3$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF3$YEAR[y] & MODIS$SITE=='SF3']
}

SF2$SEASON=NA
SF2$SEASON[SF2$DOY>=sMaAp & SF2$DOY<=eMaAp]=1
SF2$SEASON[SF2$DOY>=sJuAu & SF2$DOY<=eJuAu]=2

dTS_SF2 = aggregate(dT10 ~ SEASON+YEAR, SF2[!is.na(SF2$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF2$SEASON=as.factor(dTS_SF2$SEASON)
dTS_SF2$SITE=as.factor(13)
dTS_SF2$YEAR=dTS_SF2$YEAR-1989
dTS_SF2$ECO="ENF"
dTS_SF2=dTS_SF2[dTS_SF2$SEASON==2,]
for (y in 1:length(dTS_SF2$YEAR))
{
  dTS_SF2$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF2$YEAR[y] & MODIS$SITE=='SF2']
}

SF1$SEASON=NA
SF1$SEASON[SF1$DOY>=sMaAp & SF1$DOY<=eMaAp]=1
SF1$SEASON[SF1$DOY>=sJuAu & SF1$DOY<=eJuAu]=2

dTS_SF1 = aggregate(dT10 ~ SEASON+YEAR, SF1[!is.na(SF1$dT10),], function(x) c(mean = mean(x), sd =sd(x)/sqrt(length(x))))
dTS_SF1$SEASON=as.factor(dTS_SF1$SEASON)
dTS_SF1$SITE=as.factor(14)
dTS_SF1$YEAR=dTS_SF1$YEAR-1977
dTS_SF1$ECO="ENF"
dTS_SF1=dTS_SF1[dTS_SF1$SEASON==2,]
for (y in 1:length(dTS_SF1$YEAR))
{
  dTS_SF1$dTEMP[y]=MODIS$dTEMP[MODIS$YEAR==dTS_SF1$YEAR[y] & MODIS$SITE=='SF1']
}

dT_SUCC=rbind(dTS_RPF,dTS_FCR,dTS_BN3,dTS_BN2,dTS_NS2,dTS_NS3,dTS_NS4,dTS_NS5,dTS_NS6,dTS_SF1,dTS_SF2) # ,dTS_BN2,dTS_SF3,dTS_SF2,dTS_SF1

# plot successional changes in air-surface temperature difference
dT_SUCC=dT_SUCC[dT_SUCC$SEASON==2,]
dT_SUCC=dT_SUCC[dT_SUCC$YEAR<80,]
library(ggpmisc)

fig4 = ggplot() + 
  geom_point(dT_SUCC, mapping = aes(x = dT10[,1], y = dTEMP, colour=SITE, shape=ECO, fill=SITE),size=2,alpha = 8/10)+
  stat_poly_line(data=dT_SUCC,aes(dT10[,1], dTEMP), color = "chocolate3") +
  stat_poly_eq(data=dT_SUCC,aes(dT10[,1], dTEMP)) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  )+
  xlab("\u0394 T (surface - air) [\u00B0C]") +
  ylab("\u0394 LST (fire -control) [\u00B0C]")+
  ggtitle("July to September")+
  annotate("text", x=4.8, y=6, label= "(b)")+
  theme_bw()+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_shape_manual(values=c(21, 22, 23))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  #scale_y_continuous(expand = expand_scale(),limits = c(-3,4))+
  #scale_x_continuous(expand = expand_scale(),limits = c(0,81))+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot2 = grid.arrange(fig3, fig4, nrow = 1)
plot2
#ggsave("FigS3_dSAT.pdf", width =20, height = 10, units = "cm", plot2)
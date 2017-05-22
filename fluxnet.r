## 5/22/2017
## clean up code for calculating relationship between GPP/ER with precipitation and temperature

## use leave-one-out method and Mente Carlo to estimate the uncertainty of sensitivity


library(raster)
library(rgdal)
library(maptools)
library(sp)
library(rasterVis)

na.state = readOGR(dsn="F:\\zhihua\\dataset\\ecoregion", layer = "na.state2")
na.state = na.state[-which(na.state$NAME_1 == "Prince Edward Island"|
                           na.state$NAME_1 == "Hawaii"),]

Work.Dir <- "F:/zhihua/dataset/flux2015jul"
setwd(Work.Dir)

fluxsite.df = data.frame(read.csv("F:/zhihua/dataset/flux2015jul/selected_fluxtower_to_steve_3.csv"))
fluxsite.df = fluxsite.df[which(fluxsite.df$No_YEAR >= 4 & fluxsite.df$MOD_lc_percent > 0.5 & fluxsite.df$IGBP != "WET"),]

#REMOVE the following SITE because the daily gpp series is suspecious
#US-SRC, US-WCr, US-Wjs, US-Mpj
fluxsite.df = fluxsite.df[-which(fluxsite.df$SITE_ID == "US-SRC" | fluxsite.df$SITE_ID == "US-WCr" |
                                fluxsite.df$SITE_ID == "US-Wjs" | fluxsite.df$SITE_ID == "US-Mpj"),]

#further REMOVE the following SITE because the gpp-nee relationship is suspecious
#US-Blo [it is a plantation], US-Me2, US-Me6 [20-year old ponderosa pine trees]
fluxsite.df = fluxsite.df[-which(fluxsite.df$SITE_ID == "US-Blo" | fluxsite.df$SITE_ID == "US-Me6"),]
								
#further REMOVE the following SITE because they are irrigated cropland
#US-Ne1, US-Ne2, US-Twt
fluxsite.df = fluxsite.df[-which(fluxsite.df$SITE_ID == "US-Ne1" | fluxsite.df$SITE_ID == "US-Ne2" |
                                fluxsite.df$SITE_ID == "US-Twt"),]

flux.info2 = fluxsite.df
flux.info.sp2 = fluxsite.df						 
coordinates(flux.info.sp2) <- ~LOCATION_LONG + LOCATION_LAT
projection(flux.info.sp2) <- "+proj=longlat +datum=WGS84"

fluxnet.fn = list.files(path = "F:/zhihua/dataset/flux2015jul/download", pattern = "*.zip")
ameriflux.fn = list.files(path = "F:/zhihua/dataset/ameriflux/level2/unprocessed", pattern = "*.csv")
ameriflux.fn = ameriflux.fn[which(substr(ameriflux.fn, 8,14) == "monthly")]

plot(na.state[which(na.state$ISO == "USA" & na.state$NAME_1 != "Alaska"),])
plot(flux.info.sp2, add = T)
text(x = flux.info2$LOCATION_LONG, y = flux.info2$LOCATION_LAT, flux.info2$SITE_ID)
#############################################################################################################

#need to remove the following site and site/year
site.ex.df = data.frame(SITE_ID = c("US-AR2", rep("US-Cop",3),"US-Ha1","US-IB2","US-Me6","US-Ne1","US-Ne2","US-Ne3","US-NR1",
                                    "US-PFa", rep("US-Sta",2), rep("US-Syv",4),"US-Var","US-Whs", "US-Ced","US-MOz","US-Ton",rep("US-Wkg",2)), 
						year = c(2011,2001,2004:2005,1991,2004,2010,2013,2013,2013,1998,1995,2005:2006, 2008:2011,2000,2007,2005,2004, 2001, 2004:2005))

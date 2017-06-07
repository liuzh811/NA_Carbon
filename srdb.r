setwd("C:/zhihua/dataset/srdb_v3-2")

srdb <- read.csv( "srdb-data.csv" )
srdb_info <- read.table("srdb-data_fields.txt", sep="," )

library(rgdal)
library(raster)
library(ggplot2)

## read into USA map
usa.state = readOGR(dsn="C:\\zhihua\\dataset\\ecoregion", layer = "usa.state")
usa.state = usa.state[-which(usa.state$NAME_1 == "Hawaii" | usa.state$NAME_1 == "Alaska"),]
ext.usa = extent(usa.state)	

##examine the soil respitaion rate along the soc gradient
srdb1 <- subset( srdb, Rs_annual > 0 & Rs_annual < 4000 )

## only select years >= 2000 without manipulation
srdb4 <- subset( srdb1, Rs_annual>0 & Rs_annual < 2500 
						# & Ecosystem_state == "Natural"
						& MAP < 2000
						# & Stage == "Mature"						
						& Country == "USA" 
						& Biome == "Temperate" 
						& Ecosystem_type != "Desert" & Ecosystem_type != "Wetland" 
						& !is.na(Latitude) & !is.na(Longitude) 
						& Manipulation == "None" 
						& Study_midyear > 2000
						# & Annual_coverage
						# & Longitude > -121 & Longitude < -100 & Latitude > 30 & Latitude < 43 #select western 
						# & Longitude > -90 & Longitude < -70 & Latitude > 30 & Latitude < 43 #select eastern
						& Longitude > -130 & Latitude > 30 #select CONUS
						) #select western 
						
# plot to see distributions
world <- map_data("world")
srdb4$long <- srdb4$Longitude
srdb4$lat <- srdb4$Latitude
p1 <- ggplot( srdb4, aes( x=long, y=lat ) ) + 
        geom_point( aes( color=Ecosystem_type ) ) +
        geom_path( data=world, aes( group=group ) ) +
        scale_y_continuous( breaks=( -2:2 ) * 30 ) +
        scale_x_continuous( breaks=( -4:4 ) * 45 ) +
        # coord_fixed( xlim=c( -180, 180 ), ylim=c( -90, 90 ) ) 
		coord_fixed( xlim=c( ext.usa@xmin, ext.usa@xmax ), ylim=c( ext.usa@ymin, ext.usa@ymax ) )
print( p1 )

# plot Rs_annual and MAP correlationship
ggplot(srdb4, aes(x=MAP, y=Rs_annual)) +
    geom_point(shape=1, cex = 3) +    # Use hollow circles
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE)    # Don't add shaded confidence region

ggplot(srdb4, aes(x=MAP, y=Rs_annual, color=Ecosystem_type)) + geom_point(shape=1, cex = 3) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) # Extend regression lines

# select unique site based on lat and long
srdb5 = data.frame(aggregate(cbind(Rs_annual, MAP, MAT) ~ Latitude + Longitude, data = srdb4, "mean", na.rm = TRUE))
srdb5_sd = data.frame(aggregate(cbind(Rs_annual, MAP, MAT) ~ Latitude + Longitude, data = srdb4, "sd", na.rm = TRUE))

srdb5 = data.frame(srdb5, srdb5_sd$Rs_annual)

colnames(srdb5) <- c("Latitude","Longitude","Rs_annual","MAP", "MAT", "Rs_annual_sd")

srdb5.sp = srdb5
coordinates(srdb5.sp) <- ~Longitude + Latitude
projection(srdb5.sp) <- "+proj=longlat +datum=WGS84"

########### 5/31/2017 export srdb5.sp 
writeOGR(srdb5.sp, "F:/zhihua/dataset/srdb_v3-2", "srdb5.sp", driver="ESRI Shapefile")

############### 6/7/2017 include Ra/Rh in srdb5
srdb52 = data.frame(aggregate(cbind(Rs_annual, Ra_annual, Rh_annual, MAP, MAT) ~ Latitude + Longitude, data = srdb4, "mean", na.rm = TRUE))
srdb52_sd = data.frame(aggregate(cbind(Rs_annual, Ra_annual, Rh_annual, MAP, MAT) ~ Latitude + Longitude, data = srdb4, "sd", na.rm = TRUE))
srdb52 = data.frame(srdb52, srdb52_sd$Rs_annual,srdb52_sd$Ra_annual,srdb52_sd$Rh_annual)
colnames(srdb52) <- c("Latitude","Longitude","Rs_annual","Ra_annual", "Rh_annual", "MAP", "MAT","Rs_annual_sd","Ra_annual_sd", "Rh_annual_sd")

srdb5.sp2 = srdb52
coordinates(srdb5.sp2) <- ~Longitude + Latitude
projection(srdb5.sp2) <- "+proj=longlat +datum=WGS84"

writeOGR(srdb5.sp2, "F:/zhihua/dataset/srdb_v3-2", "srdb5.sp2", driver="ESRI Shapefile")

## assess the correlationship between Rs_annual and T/P
library(MASS)

dat.new1 = data.frame(MAT = seq(1,30, length.out = 30), MAP = mean(srdb5$MAP, na.rm = TRUE))
dat.new2 = data.frame(MAT = mean(srdb5$MAT, na.rm = TRUE), MAP = seq(50,1500, length.out = 30))

lm1.rh = lm(Rs_annual~poly(MAT,2)+poly(MAP,2), data = srdb5, na.action=na.exclude)
lm1.rh = stepAIC(lm1.rh)
#construct new data for prediction
pred1 = predict(lm1.rh, newdata = dat.new1)
pred2 = predict(lm1.rh, newdata = dat.new2)
reg.abs.mod.pred1.rh <- data.frame(MAT = dat.new1$MAT, predt = pred1, MAP = dat.new2$MAP, predp = pred2)

delt.rh = diff(reg.abs.mod.pred1.rh$predp)

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2000),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)
#plot(0, xlim = c(-5, 30), ylim = c(-2,2),bty='n',pch='',ylab='',xlab='')

lines(smooth.spline(reg.abs.mod.pred1.rh[,3], reg.abs.mod.pred1.rh[,4]), lty = 1, lwd = 6, col = "black")
points(srdb4_1$MAP, srdb4_1$Rs_annual, col = "black")

#plot spatial sensitity of gpp vs er in response to precipitation
library(oce)
plotInset(100, 800, 800, 2000,
          expr= {
		  
# plot(0, xlim = c(0, 1500), ylim = c(-20,150),bty='n',pch='',ylab='',xlab='',yaxt='n', ann=FALSE)
plot(0, xlim = c(0, 1500), ylim = c(-20,150), pch='', xlab = "Precipation (mm)", ylab = "Spatial Sensitivty (Per 50mm)", cex.lab = 1.5, cex.axis = 1.5)

#draw a box around the plot
box()
lines(smooth.spline(reg.abs.mod.pred1.rh[-1,3], delt.rh), lty = 1, lwd = 6, col = "blue")
abline(h = 0)

  },
		  mar=c(5.5,0,0,0))



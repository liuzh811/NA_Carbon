
###################################################################################
# 6/29/2017
# this code calculates the contribution of T/P to gross carbon flux (GPP/ER)


library(rgdal)
library(raster)
library(rasterVis)
require(ncdf4)

usa.state = readOGR(dsn="F:\\zhihua\\dataset\\ecoregion", layer = "usa.state")
usa.state = usa.state[-which(usa.state$NAME_1 == "Hawaii" | usa.state$NAME_1 == "Alaska"),]

######### section 2: read into MODIS GPP and aggregate into 1 degree resolution #########################################

ensm.nee1 = stack("F:/zhihua/dataset/co2_inversion/NEE.annual.mean.grd")

na.ext <- extent(-180,-48,15,85)
ensm.nee1 = crop(ensm.nee1, na.ext)
# read into modis 17 gpp
gpp.annual = stack("F:/zhihua/dataset/mod17a2/processed/c.gpp.annual.grd")
gpp.annual2 = list()

for(i in 1:15){
 
gpp.annual2[[i]] <- aggregate(gpp.annual[[i]], fact=20, fun=mean, expand=TRUE, na.rm=TRUE)
  
print(paste("Finish calculating year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.annual2 = stack(gpp.annual2)


# read into climate variables and aggregate into 1 degree resolution
# for annual temperature and precipitations
temp = list()
prep = list()

for (yr in 2000:2014){
temp1 = stack(paste0("F:/zhihua/dataset/cru_ts3.23/temp.", yr, ".grd"))
temp1 = crop(temp1, na.ext)
temp1 = calc(temp1, mean, na.rm = TRUE)
temp1 = aggregate(temp1, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)

prep1 = stack(paste0("F:/zhihua/dataset/cru_ts3.23/prep", yr, ".grd"))
prep1 = crop(prep1, na.ext)
prep1 = calc(prep1, sum, na.rm = TRUE)
prep1 = aggregate(prep1, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)

temp[[yr - 1999]] <- temp1
prep[[yr - 1999]] <- prep1

}

temp = stack(temp)
prep = stack(prep)

# calculate 3 by 3 regions
nee1 = list()
gpp1 = list()
temp1 = list()
prep1 = list()

for (i in 1:15){
nee1[[i]] <- focal(ensm.nee1[[i]], w=matrix(1/25, nc=5, nr=5))
gpp1[[i]] <- focal(gpp.annual2[[i]], w=matrix(1/25, nc=5, nr=5))
temp1[[i]] <- focal(temp[[i]], w=matrix(1/25, nc=5, nr=5))
prep1[[i]] <- focal(prep[[i]], w=matrix(1/25, nc=5, nr=5))

# nee1[[i]] <- aggregate(ensm.nee1[[i]], fact=3, fun=mean, expand=TRUE, na.rm=TRUE)
# gpp1[[i]] <- aggregate(gpp.annual2[[i]],fact=3, fun=mean, expand=TRUE, na.rm=TRUE)
# temp1[[i]] <- aggregate(temp[[i]], fact=3, fun=mean, expand=TRUE, na.rm=TRUE)
# prep1[[i]] <- aggregate(prep[[i]], fact=3, fun=mean, expand=TRUE, na.rm=TRUE)

}

temp1 = stack(temp1)
prep1 = stack(prep1)
nee1 = stack(nee1)
gpp1 = stack(gpp1)

# create mask to crop data
mask = crop(temp1[[1]], usa.state)
mask = rasterize(usa.state, mask)

mask = mask > 0
mask[mask == 0] = NA

temp1 = temp1*mask
prep1 = prep1*mask
nee1 = nee1*mask
gpp1 = gpp1*mask
er1 = gpp1 - nee1
# get each point
pts.sp = Ex.pts.all(nee1[[1]]) #get the point locations

nee.df = raster::extract(nee1, pts.sp)
gpp.df= raster::extract(gpp1, pts.sp)
temp.df = raster::extract(temp1, pts.sp)
prep.df= raster::extract(prep1, pts.sp)
er.df = gpp.df - nee.df

#calculate anomaly
Ano = function(x){x-mean(x)} 
nee.df.z = t(apply(nee.df, 1, Ano))
gpp.df.z = t(apply(gpp.df, 1, Ano))
temp.df.z = t(apply(temp.df, 1, Ano))
prep.df.z = t(apply(prep.df, 1, Ano))
er.df.z = t(apply(er.df, 1, Ano))

boot.rs1 = function(dat.gpp, dat.temp, dat.prep, ln = 12, n = 30){ 
# for each site, at least 12 year was used to get correlationship
Coef.df = c()
NofSim = 0	

if (length(which(is.na(dat.gpp))) > 4 | length(which(is.na(dat.temp))) > 4 | length(which(is.na(dat.prep))) > 4) {
Coef.df = cbind(rep(NA,n+1), rep(NA,n+1),rep(NA,n+1),rep(NA,n+1),rep(NA,n+1))
}

else {

repeat{

NofSim = NofSim+1
n1 = sample(ln:length(dat.gpp),1)
idx = sample(1:length(dat.gpp),n1)
# print(idx)
dat1 = dat.gpp[idx]
dat2 = dat.temp[idx]
dat3 = dat.prep[idx]

x = c(as.numeric(dat1),as.numeric(dat2),as.numeric(dat3))

lm1 = lm(as.numeric(x[1:length(idx)]) ~ as.numeric(x[(length(idx)+1):(length(idx)*2)])+as.numeric(x[(length(idx)*2+1):(length(idx)*3)]))
lm11 = summary(lm1)
Coef.df = rbind(Coef.df,  c(coef(lm1), lm11$r.squared, round(sqrt(mean(resid(lm11)^2)), 2)))
if(NofSim > n){break}
}

}

Coef.df = data.frame(Coef.df)
colnames(Coef.df) <- c("int","airtemp","prep","r2","rmse")
return(Coef.df)

}


dat1.list = list()
coef.gpp = list() # 
coef.er = list() # 
for (i in 1:nrow(gpp.df.z)){

gpp.df.z3 = gpp.df.z[i,]
er.df.z3 = er.df.z[i,]
temp.df.z3 = temp.df.z[i,]
prep.df.z3 = prep.df.z[i,]

coef.gpp1 = boot.rs1(dat.gpp = gpp.df.z3,dat.temp = temp.df.z3,dat.prep = prep.df.z3, n = 100)
coef.er1 = boot.rs1(dat.gpp = er.df.z3,dat.temp = temp.df.z3,dat.prep = prep.df.z3, n = 100)

coef.gpp[[i]] <- coef.gpp1
coef.er[[i]] <- coef.er1

print(paste("Finish calculating ", i, " of ", nrow(gpp.df.z), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
		  
}

#extract mean sensitivity

Rplc = function(x){
p1 = quantile(x, probs = c(0.1,0.9))
x[which(x < p1[1] | x > p1[2])] = NA
return(x)
}


er.temp = c()
gpp.temp = c()

er.prep = c()
gpp.prep = c()

gpp.r = c()
er.r = c()

for (i in 1:length(coef.er)){

er.prep = c(er.prep, median(coef.er[[i]][,3]))
gpp.prep = c(gpp.prep, median(coef.gpp[[i]][,3]))

er.temp = c(er.temp, median(coef.er[[i]][,2]))
gpp.temp = c(gpp.temp, median(coef.gpp[[i]][,2]))

gpp.r = c(gpp.r, median(coef.gpp[[i]][,4]))
er.r = c(er.r, median(coef.er[[i]][,4]))

}

#change to raster
er.tempr = Point2raster(er.temp, raster = nee1[[1]])
gpp.tempr = Point2raster(gpp.temp, raster = nee1[[1]])

er.prepr = Point2raster(er.prep, raster = nee1[[1]])
gpp.prepr = Point2raster(gpp.prep, raster = nee1[[1]])

gpp.rr = Point2raster(gpp.r, raster = nee1[[1]])
er.rr = Point2raster(er.r, raster = nee1[[1]])

# calculate contribution, follows Jung Methods
gpp1.ano = calc(gpp1, Ano)
er1.ano = calc(er1, Ano)

temp1.ano = calc(temp1, Ano)
prep1.ano = calc(prep1, Ano)

temp1.gpp = gpp.tempr*temp1.ano
prep1.gpp = gpp.prepr*prep1.ano
gpp1.ano.predict = temp1.gpp + prep1.gpp
temp1.gpp1 = calc(temp1.gpp, sd)/mean(as.matrix(calc(gpp1.ano.predict, sd)), na.rm = TRUE)
prep1.gpp1 = calc(prep1.gpp, sd)/mean(as.matrix(calc(gpp1.ano.predict, sd)), na.rm = TRUE)

temp1.er = er.tempr*temp1.ano
prep1.er = er.prepr*prep1.ano
er1.ano.predict = temp1.er + prep1.er
temp1.er1 = calc(temp1.er, sd)/mean(as.matrix(calc(er1.ano.predict, sd)), na.rm = TRUE)
prep1.er1 = calc(prep1.er, sd)/mean(as.matrix(calc(er1.ano.predict, sd)), na.rm = TRUE)

prep.mn = calc(prep1, mean)
q10 = quantile(prep.mn, probs = seq(0, 1, length.out= 20), na.rm = TRUE)
prep.mn2 = cut(prep.mn, breaks = q10)
		
se <- function(x) sqrt(var(x)/length(x))
	
dat1.df = data.frame(gpp.mean.p = zonal(prep1.gpp1, prep.mn2, mean)[,2],
						   gpp.sd.p = zonal(prep1.gpp1, prep.mn2, sd)[,2]/2,
						   gpp.mean.t = zonal(temp1.gpp1, prep.mn2, mean)[,2],
						   gpp.sd.t = zonal(temp1.gpp1, prep.mn2, sd)[,2]/2,
						   er.mean.p = zonal(prep1.er1, prep.mn2, mean)[,2],
						   er.sd.p = zonal(prep1.er1, prep.mn2, sd)[,2]/2,
						   er.mean.t = zonal(temp1.er1, prep.mn2, mean)[,2],
						   er.sd.t = zonal(temp1.er1, prep.mn2, sd)[,2]/2,
						   prep = zonal(prep.mn, prep.mn2, mean)[,2],
						   prep.sd = zonal(prep.mn, prep.mn2, sd)[,2]/2)

png("F:/zhihua/dataset/results2/contributionTP.rs.png",height = 2500, width = 2500, res = 300, units = "px")						   
						   
par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(5,5,0,0))					   
						   
# gpp plot
plot(gpp.mean.p~prep, data = dat1.df, type = "l", 
								  ylim = c(0,2), 
								  cex = 2, lwd = 4, 
								  col = "blue",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  ylab = "Mean grid-cell IAV (normalized)", 
								  xaxt='n', ann=FALSE,
								  cex.axis = 1.5, cex.lab = 1.6)
			
polygon(c(rev(dat1.df$prep), dat1.df$prep), 
		c(rev(dat1.df$gpp.mean.p-dat1.df$gpp.sd.p), dat1.df$gpp.mean.p + dat1.df$gpp.sd.p), 
        col=rgb(0, 0, 1,0.25),
		border = NA)

# overlay trendy r	
par(new=TRUE)
plot(gpp.mean.t~prep, data = dat1.df, type = "l", 
								  ylim = c(0,2), 
								  cex = 2, lwd = 4, 
								  col = "red",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  ylab = "Mean grid-cell IAV (normalized)", 
								  xaxt='n', ann=FALSE,
								  cex.axis = 1.5, cex.lab = 1.6)
			
polygon(c(rev(dat1.df$prep), dat1.df$prep), 
		c(rev(dat1.df$gpp.mean.t-dat1.df$gpp.sd.t), dat1.df$gpp.mean.t + dat1.df$gpp.sd.t), 
        col=rgb(1, 0, 0,0.25),
		border = NA)

text(x = 300, y = 1.9, "a):GPP",cex = 2)
					  
legend("topright", 
       # inset=0.05, 
	   # legend = c("ENF","DBF","MF","SHB","GRA","CRO"),
#	   legend = c("EC", "RS","TRENDY"),
       	   legend = c("Precipitation", "Temperature"),
	   horiz=F,
	   # pch = 0:2,
	   lwd = 4,
	   col = c(rgb(0, 0, 1,1), rgb(1,0,0,1)),
	   text.col = c(rgb(0, 0, 1,1), rgb(1,0,0,1)),
	   cex = 1.5,
	   box.col = "transparent",
       bg = "transparent")	

	   
# er plot
plot(er.mean.p~prep, data = dat1.df, type = "l", 
								  ylim = c(0,2), 
								  cex = 2, lwd = 4, 
								  col = "blue",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  ylab = "Mean grid-cell IAV (normalized)", cex.axis = 1.5, cex.lab = 1.6)
			
polygon(c(rev(dat1.df$prep), dat1.df$prep), 
		c(rev(dat1.df$er.mean.p-dat1.df$er.sd.p), dat1.df$er.mean.p + dat1.df$er.sd.p), 
        col=rgb(0, 0, 1,0.25),
		border = NA)

# overlay trendy r	
par(new=TRUE)
plot(er.mean.t~prep, data = dat1.df, type = "l", 
								  ylim = c(0,2), 
								  cex = 2, lwd = 4, 
								  col = "red",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  ylab = "Mean grid-cell IAV (normalized)", cex.axis = 1.5, cex.lab = 1.6)
			
polygon(c(rev(dat1.df$prep), dat1.df$prep), 
		c(rev(dat1.df$er.mean.t-dat1.df$er.sd.t), dat1.df$er.mean.t + dat1.df$er.sd.t), 
        col=rgb(1, 0, 0,0.25),
		border = NA)

text(x = 300, y = 1.9, "b):ER",cex = 2)
		
mtext(side = 1, line = 3, "Anuual Mean Precipitation (mm)", 
      outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))
				
mtext(side = 2, line = 3, "Mean grid-cell IAV (normalized)", 
outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))

dev.off()

## plot longitude/latitude gradient
temp1.gpp2 = as.matrix(temp1.gpp1)
prep1.gpp2 = as.matrix(prep1.gpp1)
lat = yFromRow(nee1[[1]], 1:nrow(nee1[[1]]))
long = xFromCol(nee1[[1]], 1:ncol(nee1[[1]]))

#latitude
temp1.gpp21 = data.frame(latmean = apply(temp1.gpp2, 1, mean, na.rm = TRUE),latsd = 0.5*apply(temp1.gpp2, 1, sd, na.rm = TRUE))
prep1.gpp21 = data.frame(latmean = apply(prep1.gpp2, 1, mean, na.rm = TRUE),latsd = 0.5*apply(prep1.gpp2, 1, sd, na.rm = TRUE))

temp1.gpp21 = data.frame(long = lat, temp1.gpp21)
prep1.gpp21 = data.frame(long = lat, prep1.gpp21)

plot(latmean~long,  type = "l", lwd = 2, data = temp1.gpp21,
                                  ylim = c(0,2), 
								  col="red",
                                  xlab = "Latitude", 
								  ylab = "Mean grid-cell IAV (normalized)", cex.axis = 1.5, cex.lab = 1.6)

polygon(c(rev(temp1.gpp21$long), temp1.gpp21$long), 
		c(rev(temp1.gpp21$latmean-temp1.gpp21$latsd), temp1.gpp21$latmean+temp1.gpp21$latsd), 
        col=rgb(1, 0, 0, 0.25),
		border = NA)
		
par(new=TRUE)
plot(latmean~long,  type = "l", lwd = 2, data = prep1.gpp21,
                                  ylim = c(0,2), 
								  col="blue",
                                  xlab = "Latitude", 
								  ylab = "Percent High Severity Fire", cex.axis = 1.5, cex.lab = 1.6)
			
polygon(c(rev(prep1.gpp21$long), prep1.gpp21$long), 
		c(rev(prep1.gpp21$latmean-prep1.gpp21$latsd), prep1.gpp21$latmean+prep1.gpp21$latsd), 
         col=rgb(0, 0, 1,0.25),
		border = NA)


# longitude
# prepared GPP data
temp1.gpp22 = data.frame(latmean = apply(temp1.gpp2, 2, mean, na.rm = TRUE),latsd = 0.5*apply(temp1.gpp2, 2, sd, na.rm = TRUE))
prep1.gpp22 = data.frame(latmean = apply(prep1.gpp2, 2, mean, na.rm = TRUE),latsd = 0.5*apply(prep1.gpp2, 2, sd, na.rm = TRUE))

temp1.gpp22[sapply(temp1.gpp22,is.na)] = NA 

temp1.gpp22 = data.frame(long = long, temp1.gpp22)
prep1.gpp22 = data.frame(long = long, prep1.gpp22)

temp1.gpp22 = temp1.gpp22[complete.cases(temp1.gpp22),]
prep1.gpp22 = prep1.gpp22[complete.cases(prep1.gpp22),]

# prepared ER data
temp1.er2 = as.matrix(temp1.er1)
prep1.er2 = as.matrix(prep1.er1)

temp1.er22 = data.frame(latmean = apply(temp1.er2, 2, mean, na.rm = TRUE),latsd = 0.5*apply(temp1.er2, 2, sd, na.rm = TRUE))
prep1.er22 = data.frame(latmean = apply(prep1.er2, 2, mean, na.rm = TRUE),latsd = 0.5*apply(prep1.er2, 2, sd, na.rm = TRUE))

temp1.er22[sapply(temp1.er22,is.na)] = NA 

temp1.er22 = data.frame(long = long, temp1.er22)
prep1.er22 = data.frame(long = long, prep1.er22)

temp1.er22 = temp1.er22[complete.cases(temp1.er22),]
prep1.er22 = prep1.gpp22[complete.cases(prep1.er22),]

png("F:/zhihua/dataset/results2/contributionTP.rs.long.png",height = 1500, width = 2000, res = 300, units = "px")						   

par(mfrow=c(2,1),mar=c(0,0,0,0),oma=c(5,5,0,0))					   

# gpp plot
plot(latmean~long,  type = "l", lwd = 2, data = temp1.gpp22,
                                  ylim = c(0,2), xlim = c(-125,-70),
				  col="red",
                                  xlab = "Longitude", 
				  ylab = "Mean grid-cell IAV (normalized)", 
                                  xaxt='n', ann=FALSE,
                                  cex.axis = 1.5, cex.lab = 1.6)

 polygon(c(rev(temp1.gpp22$long), temp1.gpp22$long), 
		 c(rev(temp1.gpp22$latmean-temp1.gpp22$latsd), temp1.gpp22$latmean+temp1.gpp22$latsd), 
         col=rgb(1, 0, 0, 0.25),
		 border = NA)
								  
par(new=TRUE)
plot(latmean~long,  type = "l", lwd = 2, data = prep1.gpp22,
                                  ylim = c(0,2), xlim = c(-125,-70),
								  col="blue",
                                  xlab = "Longitude", 
				  ylab = "Mean grid-cell IAV (normalized)", 
                                  xaxt='n', ann=FALSE,
                                  cex.axis = 1.5, cex.lab = 1.6)
			
 polygon(c(rev(prep1.gpp22$long), prep1.gpp22$long), 
		 c(rev(prep1.gpp22$latmean-prep1.gpp22$latsd), prep1.gpp22$latmean+prep1.gpp22$latsd), 
          col=rgb(0, 0, 1,0.25),
		 border = NA)
text(x = -120, y = 1.9, "a):GPP",cex = 2)

# er plot
plot(latmean~long,  type = "l", lwd = 2, data = temp1.er22,
                                  ylim = c(0,2), xlim = c(-125,-70),
								  col="red",
                                  xlab = "Longitude", 
								  ylab = "Mean grid-cell IAV (normalized)", cex.axis = 1.5, cex.lab = 1.6)

 polygon(c(rev(temp1.er22$long), temp1.er22$long), 
		 c(rev(temp1.er22$latmean-temp1.er22$latsd), temp1.er22$latmean+temp1.er22$latsd), 
         col=rgb(1, 0, 0, 0.25),
		 border = NA)
								  
par(new=TRUE)
plot(latmean~long,  type = "l", lwd = 2, data = prep1.er22,
                                  ylim = c(0,2), xlim = c(-125,-70),
								  col="blue",
                                  xlab = "Longitude", 
								  ylab = "Mean grid-cell IAV (normalized)", cex.axis = 1.5, cex.lab = 1.6)
			
 polygon(c(rev(prep1.er22$long), prep1.er22$long), 
		 c(rev(prep1.er22$latmean-prep1.er22$latsd), prep1.er22$latmean+prep1.er22$latsd), 
          col=rgb(0, 0, 1,0.25),
		 border = NA)

text(x = -120, y = 1.9, "b):ER",cex = 2)

mtext(side = 1, line = 3, "Longitude", 
      outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))
				
mtext(side = 2, line = 3, "Mean grid-cell IAV (normalized)", 
outer = TRUE, cex = 1.6, col = rgb(0, 0, 0,1))

dev.off()


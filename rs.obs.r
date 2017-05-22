## 5/22/2017
## using remote sensing and co2 inversion method to calculate the relationship between GPP/ER with precipitation and temperature
## using the median/mean NEE estimate from four co2 inversion methods

library(rgdal)
library(raster)
library(rasterVis)
require(ncdf4)

usa.state = readOGR(dsn="F:\\zhihua\\dataset\\ecoregion", layer = "usa.state")
usa.state = usa.state[-which(usa.state$NAME_1 == "Hawaii" | usa.state$NAME_1 == "Alaska"),]

######### section 1: read into CO2 inversion and calculate mean/median #########################################
#########             refer to code at jenoCarboScope.r                #########################################
# read into CT NEE: 2000 - 2014
ct.nee = stack("F:/zhihua/dataset/ct/ct2015flux/processed/c.flux.annual.grd")
# ct.nee = crop(ct.nee, na.state)
ct.nee = -1*ct.nee

#read and process CTE: 2001 - 2014
nc <- nc_open("F:/zhihua/dataset/co2_inversion/cte/flux1x1_all_years.nc")
print(paste("The file has",nc$nvars,"variables"))

#variable 5 is the bio_flux_opt
v2 <- nc$var[[5]]
data2 <- ncvar_get( nc, v2 ) #data2 is an array
print(paste("Var 2 has name",v2$name,"and is of shape",dim(data2),
	". Here are the values:"))

# float ocn_flux_opt[longitude,latitude,date]   (Chunking: [360,180,1])  
# comment: time-interval average, centered on times in the date axis
# long_name: Surface flux of carbon dioxide, open ocean , optimized 
# standard_name: surface_carbon_dioxide_mole_flux
# units: mol m-2 s-1	
days.mon = c(31,28,31,30,31,30,31,31,30,31,30,31)
sec.mon = days.mon*86400

data2 = data2*sum(sec.mon)*12.011 #convert to g C/m-2/yr

data2.t <- aperm(data2, c(2,1,3))
data2.r = brick(data2.t,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r,direction="y")

cte.nee = -1*data2.r

# add NA for 2000 estimate
cte.nee = stack(cte.nee[[1]], cte.nee)
cte.nee[[1]][] <- NA

# read into CAMS NEE: 2000-2014
cams.nee = stack("F:/zhihua/dataset/co2_inversion/cams/CAMS.NEE.2000-2014.grd")
# cams.nee = crop(cams.nee, na.state)
cams.nee = -1*cams.nee

# read into jena NEE:2000-2014
jena.nee = stack("F:/zhihua/dataset/co2_inversion/JenaCarboScope/jena.NEE.s99_v3.8.2000-2014.grd")
# jena.nee = crop(jena.nee, na.state)
jena.nee = -1*jena.nee

# compute the mean value from three co2 inversions
# change to 1 degree first
cams.nee2 = list()
jena.nee2 = list()

for (i in 1:15){
cams.nee2[[i]] <- resample(cams.nee[[i]], ct.nee[[1]], method = "ngb")
jena.nee2[[i]] <- resample(jena.nee[[i]], ct.nee[[1]], method = "ngb")
}

cams.nee2 = stack(cams.nee2)
jena.nee2 = stack(jena.nee2)


# calculate the mean/median value 
ensm.nee1 = list() # mean value
ensm.nee2 = list() # median value
for (i in 1:15){
temp.r = stack(cams.nee2[[i]], 
			   jena.nee2[[i]], 
			   ct.nee[[i]], 
			   cte.nee[[i]])
ensm.nee1[[i]] <- calc(temp.r, mean, na.rm = TRUE)
ensm.nee2[[i]] <- calc(temp.r, median, na.rm = TRUE)

}
ensm.nee1 = stack(ensm.nee1)
ensm.nee2 = stack(ensm.nee2)

######### section 2: read into MODIS GPP and aggregate into 1 degree resolution #########################################
na.ext <- extent(-180,-48,15,85)

# read into modis 17 gpp
gpp.annual = stack("F:/zhihua/dataset/mod17a2/processed/c.gpp.annual.grd")
gpp.annual2 = list()

for(i in 1:15){
t1 = aggregate(gpp.annual[[i]], fact=20, fun=mean, expand=TRUE, na.rm=TRUE)
gpp.annual2[[i]] <- t1
  
print(paste("Finish calculating year ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.annual2 = stack(gpp.annual2)


######### section 3: read into various grid and polygon #########################################
## read into cell-area maps
library(ncdf4)
nc <- nc_open("F:/zhihua/dataset/ecoregion/regions.nc")
print(nc)
v2 <- nc$var[[1]]
data2 <- ncvar_get( nc, v2 ) #data2 is an array

grid.area = raster(t(data2),xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

grid.area.gpp = resample(grid.area.nee, mask.gpp, method = "ngb"); grid.area.gpp = grid.area.gpp/400

## read into NEON domains
neon.sp = readOGR(dsn="F:/zhihua/dataset/ecoregion", layer="NEON_Domains")

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

# create mask to crop data
temp1 = crop(temp[[1]], usa.state)
mask = rasterize(usa.state, temp1)

mask = mask > 0
mask[mask == 0] = NA

# crop dataset into USA
nee.annual = crop(ensm.nee1, usa.state) # mean
nee.annual = nee.annual*mask

gpp.annual = gpp.annual2*mask
temp = temp*mask
prep = prep*mask

# change NEON into grid
noen.grd = rasterize(neon.sp,mask, field = "DomainID")

######### section 4: calculate correlationship #########################################
nee.df = zonal(nee.annual, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)] 
gpp.df = zonal(gpp.annual, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]
er.df = gpp.df - nee.df
temp.df = zonal(temp, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]
prep.df = zonal(prep, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]

temp.df.mn = zonal(calc(temp, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2]
prep.df.mn = zonal(calc(prep, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2]
nee.df.mn = zonal(calc(nee.annual, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2] 
gpp.df.mn = zonal(calc(gpp.annual, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2]
er.df.mn = gpp.df.mn - nee.df.mn

#calculate anomaly
Ano = function(x){x-mean(x)} 
nee.df.z = t(apply(nee.df, 1, Ano))
gpp.df.z = t(apply(gpp.df, 1, Ano))
temp.df.z = t(apply(temp.df, 1, Ano))
prep.df.z = t(apply(prep.df, 1, Ano))
er.df.z = t(apply(er.df, 1, Ano))

#calculate cofficient, like flux tower processing
# first value is coefficient, 2nd is intercept, 3rd is p value
fun.cor2 <- function(x) {
  if (length(which(is.na(x[1:15]))) > 8 | length(which(is.na(x[16:30])))>8 | length(which(is.na(x[31:45])))>8) 
  {return(c(NA,NA,NA,NA,NA))} 
  else 
  {
    lm1 = lm(as.numeric(x[1:15]) ~ as.numeric(x[16:30])+as.numeric(x[31:45]))
	lm11 = summary(lm1)
    return(c(lm11$coefficients[2:3,1],lm11$coefficients[2:3,4], lm11$r.squared))
  }
}

#plot coefficient, use as temproal sensitivity
# need to remove region 4
gpp.df.z = cbind(gpp.df.z, temp.df.z, prep.df.z)
er.df.z = cbind(er.df.z, temp.df.z, prep.df.z)

gpp.df.z.r = apply(gpp.df.z, 1, fun.cor2) 
er.df.z.r = apply(er.df.z, 1, fun.cor2) 

mean(gpp.df.z.r[5,]) # 0.40+-0.21
mean(er.df.z.r[5,]) # 0.37+-0.19


dat1 = data.frame(prep = prep.df.mn[-4],gpp.p = gpp.df.z.r[2,-4],er.p = er.df.z.r[2,-4], 
				  temp = temp.df.mn[-4],gpp.t = gpp.df.z.r[1,-4],er.t = er.df.z.r[1,-4])
dat1$p.dif = dat1$gpp.p - dat1$er.p
dat1$lc = 1
dat1$lc[dat1$prep > 750] = 2

#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.mn[-4], nee = nee.df.mn[-4],er = er.df.mn[-4],airtemp = temp.df.mn[-4],prep = prep.df.mn[-4])

lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.er = stepAIC(lm1.er)

#construct new data for prediction
dat.new1 = data.frame(airtemp = seq(1,23, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

pred1 = predict(lm1.er, newdata = dat.new1)
pred2 = predict(lm1.er, newdata = dat.new2)
reg.abs.mod.pred1.er <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.gpp = stepAIC(lm1.gpp)

pred1 = predict(lm1.gpp, newdata = dat.new1)
pred2 = predict(lm1.gpp, newdata = dat.new2)
reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

lm1.nee = lm(nee~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
#lm1.nee = stepAIC(lm1.nee)

pred1 = predict(lm1.nee, newdata = dat.new1)
pred2 = predict(lm1.nee, newdata = dat.new2)
reg.abs.mod.pred1.nee <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

delt.er = diff(reg.abs.mod.pred1.er$predp)
delt.gpp = diff(reg.abs.mod.pred1.gpp$predp)


##### plot

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(-0.2,2000),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)

lines(smooth.spline(reg.abs.mod.pred1.gpp[,3], reg.abs.mod.pred1.gpp[,4]), lty = 1, lwd = 6, col = "green")
points(dat.df.annual.mean$prep, dat.df.annual.mean$gpp, col = "green")

lines(smooth.spline(reg.abs.mod.pred1.er[,3], reg.abs.mod.pred1.er[,4]), lty = 1, lwd = 6, col = "red")
points(dat.df.annual.mean$prep, dat.df.annual.mean$er, col = "red")

library(oce)
plotInset(750,0, 1450, 800,
          expr= {
		  
# plot(0, xlim = c(0, 1500), ylim = c(0,100), pch='', xlab = "Precipation (mm)", ylab = "g C m-2 yr-1", cex.lab = 1.5, cex.axis = 1.5)
plot(0, xlim = c(0, 1500), ylim = c(25,75)*2, pch='', xlab = "Precipation (mm)", 
     # ylab = "Spatial Sensitivty (Per 50mm)", 
	 ylab = expression(paste(beta ["spatial"])),
	 cex.lab = 1.5, cex.axis = 1.5)

#draw a box around the plot
box()
lines(smooth.spline(reg.abs.mod.pred1.gpp[-1,3], delt.gpp*2), lty = 1, lwd = 6, col = "green")
lines(smooth.spline(reg.abs.mod.pred1.er[-1,3], delt.er*2), lty = 1, lwd = 6, col = "red")

abline(h = 0)
# abline(v = reg.abs.mod.pred1.gpp[-1,3][which.max(which(delt.gpp > delt.er))])
abline(v = 750)

	  },	  
		  mar=c(5,0,0,0))
		  
library(oce)
plotInset(200, 700, 1000, 2000,
          expr= {
		 		  
#plot temproal sensitivity
# plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',xlab='Temproal Sensitivity (Per 50mm)',yaxt='n', ann=FALSE)

plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',
   # xlab='Temproal Sensitivity (Per 50mm)',
   xlab = expression(paste(beta ["temporal"])),
   yaxt='n', 
   # ann=FALSE
   cex.lab = 1.5, cex.axis = 1.5)
   

#draw a box around the plot
box()
abline(v = 0)

# high prep : MAP > 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 2, x1 = x1[4], y1 = 2, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 2, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 1.5, "MAP > 750 mm")
# paste("MAP > 750 mm ( n = ", length(which(dat1$lc == 2)), " )", sep = ""),  cex = 1)
  

# low prep : : MAP < 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, "MAP < 750 mm")
# paste("MAP < 750 mm ( n = ", length(which(dat1$lc == 1)), " )", sep = ""),  cex = 1)

################## add ER         #################################
# high prep : MAP > 750 mm
x11 = dat1$er.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 1, x1 = x1[4], y1 = 1, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 1, type = "p", pch = 17, cex = 2, col = "red")

# low prep : : MAP < 750 mm
x11 = dat1$er.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 3, x1 = x1[4], y1 = 3, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 3, type = "p", pch = 17, cex = 2, col = "red")

	  },
	  
		  mar=c(5,0,0,0))


######### section 4: calculate correlationship: leave one out method #########################################
############################## leave-one-out methods ###################################

# first value is coefficient, 2nd is intercept, 3rd is p value
# the order of the 5 value, coef for temp and prep, p value for temp and prep, r squared

fun.cor3 <- function(x) {
  if (length(which(is.na(x[1:14]))) > 7 | length(which(is.na(x[15:28])))> 7 | length(which(is.na(x[29:42])))>8) 
  {return(c(NA,NA,NA,NA,NA))} 
  else 
  {
    lm1 = lm(as.numeric(x[1:14]) ~ as.numeric(x[15:28])+as.numeric(x[29:42]))
	lm11 = summary(lm1)
    return(c(lm11$coefficients[2:3,1],lm11$coefficients[2:3,4], lm11$r.squared))
  }
}

#calculate anomaly
Ano = function(x){x-mean(x)} 
nee.df.z = t(apply(nee.df, 1, Ano))
gpp.df.z = t(apply(gpp.df, 1, Ano))
temp.df.z = t(apply(temp.df, 1, Ano))
prep.df.z = t(apply(prep.df, 1, Ano))
er.df.z = t(apply(er.df, 1, Ano))

# temproal sensitivity need to remove region 4
gpp.df.z2 = gpp.df.z[-4,]
er.df.z2 = er.df.z[-4,]
temp.df.z2 = temp.df.z[-4,]
prep.df.z2 = prep.df.z[-4,]

dat1.list = list()

for (i in 1:ncol(gpp.df.z2)){

gpp.df.z1 = cbind(gpp.df.z2[,-i], temp.df.z2[,-i], prep.df.z2[,-i])
er.df.z1 = cbind(er.df.z2[,-i], temp.df.z2[,-i], prep.df.z2[,-i])

gpp.df.z.r = apply(gpp.df.z1, 1, fun.cor3) 
er.df.z.r = apply(er.df.z1, 1, fun.cor3) 

dat1 = data.frame(prep = prep.df.mn[-4],gpp.p = gpp.df.z.r[2,],er.p = er.df.z.r[2,], 
				  temp = temp.df.mn[-4],gpp.t = gpp.df.z.r[1,],er.t = er.df.z.r[1,])

dat1.list[[i]] <- dat1		  
}

#get into one table, and calculate the mean
gpp.temporal2 = c()
er.temporal2 = c()

for(i in 1:length(dat1.list)){

gpp.temporal2 = cbind(gpp.temporal2, dat1.list[[i]][,2])
er.temporal2 = cbind(er.temporal2, dat1.list[[i]][,3])

}

apply(gpp.temporal2, 1, mean)
dat1 = data.frame(prep = prep.df.mn[-4],gpp.p = apply(gpp.temporal2, 1, mean),er.p = apply(er.temporal2, 1, mean))
dat1$lc = 1
dat1$lc[dat1$prep > 750] = 2


#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.mn[-4], nee = nee.df.mn[-4],er = er.df.mn[-4],airtemp = temp.df.mn[-4],prep = prep.df.mn[-4])

reg.abs.mod.pred1.er.list = list()
reg.abs.mod.pred1.gpp.list = list()

delt.er = c()
delt.gpp = c()

for (i in 1:nrow(dat.df.annual.mean)){

dat.df.annual.mean1 = dat.df.annual.mean[-i,]

lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)

dat.new1 = data.frame(airtemp = seq(1,23, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

#construct new data for prediction
pred1 = predict(lm1.er, newdata = dat.new1)
pred2 = predict(lm1.er, newdata = dat.new2)
reg.abs.mod.pred1.er <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)
pred1 = predict(lm1.gpp, newdata = dat.new1)
pred2 = predict(lm1.gpp, newdata = dat.new2)
reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

delt.er1 = diff(reg.abs.mod.pred1.er$predp)
delt.gpp1 = diff(reg.abs.mod.pred1.gpp$predp)

reg.abs.mod.pred1.er.list[[i]] <- reg.abs.mod.pred1.er
reg.abs.mod.pred1.gpp.list[[i]] <- reg.abs.mod.pred1.gpp

delt.er = cbind(delt.er, delt.er1)
delt.gpp = cbind(delt.gpp, delt.gpp1)

}

pred.er = c()
pred.gpp = c()

for (i in 1:length(reg.abs.mod.pred1.er.list)){
pred.er = cbind(pred.er, reg.abs.mod.pred1.er.list[[i]][,4])
pred.gpp = cbind(pred.gpp, reg.abs.mod.pred1.gpp.list[[i]][,4])
}

## plot using the confidence intervals
# functions
cifun <- function(data, ALPHA = 0.05){
  c(mean(data) - qnorm(1-ALPHA/2) * sd(data)/sqrt(length(data)),
    mean(data),
    mean(data) + qnorm(1-ALPHA/2) * sd(data)/sqrt(length(data)))
}

cifun(pred.er[1,], 0.1) 

cifun2 <- function(data){
  c(mean(data) - sd(data),
    mean(data),
    mean(data) + sd(data))
}

pred.er2 = apply(pred.er, 1, cifun)
pred.gpp2 = apply(pred.gpp, 1, cifun)

delt.gpp2 = apply(delt.gpp, 1, cifun)
delt.er2 = apply(delt.er, 1, cifun)

prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(-0.2,2000),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)

### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.gpp2[3,], lty = 'dashed', col = 'green')
lines(prep.grd,pred.gpp2[1,], lty = 'dashed', col = 'green')
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")
points(dat.df.annual.mean$prep, dat.df.annual.mean$gpp, col = "green")


### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.er2[3,], lty = 'dashed', col = 'red')
lines(prep.grd,pred.er2[1,], lty = 'dashed', col = 'red')
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
points(dat.df.annual.mean$prep, dat.df.annual.mean$er, col = "red")

# plot spatial sensitivity
library(oce)
plotInset(750,0, 1450, 800,
          expr= {
		  
# plot(0, xlim = c(0, 1500), ylim = c(0,100), pch='', xlab = "Precipation (mm)", ylab = "g C m-2 yr-1", cex.lab = 1.5, cex.axis = 1.5)
plot(0, xlim = c(0, 1500), ylim = c(25,75)*2, pch='', xlab = "Precipation (mm)", 
     # ylab = "Spatial Sensitivty (Per 50mm)", 
	 ylab = expression(paste(beta ["spatial"])),
	 cex.lab = 1.5, cex.axis = 1.5)

#draw a box around the plot
box()
#######  GPP
lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# add fill
# polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
        # c(rev(delt.gpp2[3,]*2), delt.gpp2[1,]*2), col = 'grey80', border = NA)
# lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# intervals
lines(prep.grd[-1],delt.gpp2[3,]*2, lty = 'dashed', col = 'green')
lines(prep.grd[-1],delt.gpp2[1,]*2, lty = 'dashed', col = 'green')

#######  ER
lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# add fill
# polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
        # c(rev(delt.er2[3,]*2), delt.er2[1,]*2), col = 'grey80', border = NA)
# lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# intervals
lines(prep.grd[-1],delt.er2[3,]*2, lty = 'dashed', col = 'red')
lines(prep.grd[-1],delt.er2[1,]*2, lty = 'dashed', col = 'red')

abline(h = 0)
# abline(v = reg.abs.mod.pred1.gpp[-1,3][which.max(which(delt.gpp > delt.er))])
abline(v = 750)

	  },	  
		  mar=c(5,0,0,0))

# temproal sensitivity
	  
library(oce)
plotInset(200, 700, 1000, 2000,
          expr= {
		 		  
#plot temproal sensitivity
# plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',xlab='Temproal Sensitivity (Per 50mm)',yaxt='n', ann=FALSE)

plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',
   # xlab='Temproal Sensitivity (Per 50mm)',
   xlab = expression(paste(beta ["temporal"])),
   yaxt='n', 
   # ann=FALSE
   cex.lab = 1.5, cex.axis = 1.5)
   

#draw a box around the plot
box()
abline(v = 0)

# high prep : MAP > 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 2, x1 = x1[4], y1 = 2, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 2, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 1.5, "MAP > 750 mm")
# paste("MAP > 750 mm ( n = ", length(which(dat1$lc == 2)), " )", sep = ""),  cex = 1)
  

# low prep : : MAP < 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, "MAP < 750 mm")
# paste("MAP < 750 mm ( n = ", length(which(dat1$lc == 1)), " )", sep = ""),  cex = 1)

################## add ER         #################################
# high prep : MAP > 750 mm
x11 = dat1$er.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 1, x1 = x1[4], y1 = 1, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 1, type = "p", pch = 17, cex = 2, col = "red")

# low prep : : MAP < 750 mm
x11 = dat1$er.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 3, x1 = x1[4], y1 = 3, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 3, type = "p", pch = 17, cex = 2, col = "red")

	  },
	  
		  mar=c(5,0,0,0))

############################## END of leave-one-out methods ###################################

  
##############################       bootstrap methods            ##############################
## use bootstrap methodsm: produce the same results as leave one out method above
## leave 1, 2, or 3 out, and calculate the coefficient
## repeat 999 times
## create 95% CI

gpp.df.z1 = cbind(gpp.df.z, temp.df.z, prep.df.z)
er.df.z1 = cbind(er.df.z, temp.df.z, prep.df.z)

# for temporal sensitivity
gpp.beta.temporal = list() # each element contains one region
er.beta.temporal = list()

for (i in 1:nrow(gpp.df.z)){
## calcualte GPP

# select regions
dat11 = data.frame(gpp = gpp.df.z[i,], temp = temp.df.z[i,], prep = prep.df.z[i,])

coef1 = c()
N = 0
repeat{
N = N+1
# select number of element to removed random
n1 = sample(c(1:3),1)
#select index of the element to be used for regression
l1 = nrow(dat11) - n1
dat12 = dat11[sample(c(1:nrow(dat11)), l1), ]
lm1 = lm(gpp ~ temp + prep, data = dat12)
coef1 = rbind(coef1, coef(lm1))
if(N >= 100){break}

	}
	
gpp.beta.temporal[[i]] = apply(coef1, 2, cifun)
	
## calcualte ER

# select regions
dat11 = data.frame(er = er.df.z[i,], temp = temp.df.z[i,], prep = prep.df.z[i,])

coef1 = c()
N = 0
repeat{
N = N+1
# select number of element to removed random
n1 = sample(c(1:3),1)
#select index of the element to be used for regression
l1 = nrow(dat11) - n1
dat12 = dat11[sample(c(1:nrow(dat11)), l1), ]
lm1 = lm(er ~ temp + prep, data = dat12)
coef1 = rbind(coef1, coef(lm1))
if(N >= 100){break}

	}
	
er.beta.temporal[[i]] = apply(coef1, 2, cifun)

print(paste("Finish calculating for region ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

gpp.beta.temporal2 = c()
er.beta.temporal2 = c()

for(i in 1:length(gpp.beta.temporal)){

gpp.beta.temporal2 = rbind(gpp.beta.temporal2, gpp.beta.temporal[[i]][2,])
er.beta.temporal2 = rbind(er.beta.temporal2, er.beta.temporal[[i]][2,])

}

dat1 = data.frame(prep = prep.df.mn[-4],gpp.p = gpp.beta.temporal2[-4,3],er.p = er.beta.temporal2[-4,3], 
				  temp = temp.df.mn[-4],gpp.t = gpp.beta.temporal2[-4,2],er.t = er.beta.temporal2[-4,2])
dat1$lc = 1
dat1$lc[dat1$prep > 750] = 2


#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.mn[-4], nee = nee.df.mn[-4],er = er.df.mn[-4],airtemp = temp.df.mn[-4],prep = prep.df.mn[-4])

reg.abs.mod.pred1.er.list2 = list()
reg.abs.mod.pred1.gpp.list2 = list()

delt.er = c()
delt.gpp = c()

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

## calcualte GPP
N = 0
repeat{
N = N+1
# select number of element to removed random
n1 = sample(c(1:3),1)
#select index of the element to be used for regression
l1 = nrow(dat.df.annual.mean) - n1
dat.df.annual.mean1 = dat.df.annual.mean[sample(c(1:nrow(dat.df.annual.mean)), l1), ]

lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)

#construct new data for prediction
pred1 = predict(lm1.er, newdata = dat.new1)
pred2 = predict(lm1.er, newdata = dat.new2)
reg.abs.mod.pred1.er <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)
pred1 = predict(lm1.gpp, newdata = dat.new1)
pred2 = predict(lm1.gpp, newdata = dat.new2)
reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

delt.er1 = diff(reg.abs.mod.pred1.er$predp)
delt.gpp1 = diff(reg.abs.mod.pred1.gpp$predp)

reg.abs.mod.pred1.er.list2[[N]] <- reg.abs.mod.pred1.er
reg.abs.mod.pred1.gpp.list2[[N]] <- reg.abs.mod.pred1.gpp

delt.er = cbind(delt.er, delt.er1)
delt.gpp = cbind(delt.gpp, delt.gpp1)

if(N >= 100){break}

}


pred.er = c()
pred.gpp = c()

for (i in 1:length(reg.abs.mod.pred1.er.list2)){
pred.er = cbind(pred.er, reg.abs.mod.pred1.er.list2[[i]][,4])
pred.gpp = cbind(pred.gpp, reg.abs.mod.pred1.gpp.list2[[i]][,4])
}


pred.er2 = apply(pred.er, 1, cifun2)
pred.gpp2 = apply(pred.gpp, 1, cifun2)

delt.gpp2 = apply(delt.gpp, 1, cifun2)
delt.er2 = apply(delt.er, 1, cifun2)

prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(-0.2,2000),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)

### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 1, col = "green")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.gpp2[3,], lty = 'dashed', col = 'green')
lines(prep.grd,pred.gpp2[1,], lty = 'dashed', col = 'green')
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 1, col = "green")
points(dat.df.annual.mean$prep, dat.df.annual.mean$gpp, col = "green")


### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.er2[3,], lty = 'dashed', col = 'red')
lines(prep.grd,pred.er2[1,], lty = 'dashed', col = 'red')
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
points(dat.df.annual.mean$prep, dat.df.annual.mean$er, col = "red")

# plot spatial sensitivity
library(oce)
plotInset(750,0, 1450, 800,
          expr= {
		  
# plot(0, xlim = c(0, 1500), ylim = c(0,100), pch='', xlab = "Precipation (mm)", ylab = "g C m-2 yr-1", cex.lab = 1.5, cex.axis = 1.5)
plot(0, xlim = c(0, 1500), ylim = c(25,75)*2, pch='', xlab = "Precipation (mm)", 
     # ylab = "Spatial Sensitivty (Per 50mm)", 
	 ylab = expression(paste(beta ["spatial"])),
	 cex.lab = 1.5, cex.axis = 1.5)

#draw a box around the plot
box()
#######  GPP
lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# add fill
# polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
        # c(rev(delt.gpp2[3,]*2), delt.gpp2[1,]*2), col = 'grey80', border = NA)
# lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# intervals
lines(prep.grd[-1],delt.gpp2[3,]*2, lty = 'dashed', col = 'green')
lines(prep.grd[-1],delt.gpp2[1,]*2, lty = 'dashed', col = 'green')

#######  ER
lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# add fill
# polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
        # c(rev(delt.er2[3,]*2), delt.er2[1,]*2), col = 'grey80', border = NA)
# lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# intervals
lines(prep.grd[-1],delt.er2[3,]*2, lty = 'dashed', col = 'red')
lines(prep.grd[-1],delt.er2[1,]*2, lty = 'dashed', col = 'red')

abline(h = 0)
# abline(v = reg.abs.mod.pred1.gpp[-1,3][which.max(which(delt.gpp > delt.er))])
abline(v = 750)

	  },	  
		  mar=c(5,0,0,0))
		  
#plot temporal sensitivity
		  
library(oce)
plotInset(200, 700, 1000, 2000,
          expr= {
		 		  
#plot temproal sensitivity
# plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',xlab='Temproal Sensitivity (Per 50mm)',yaxt='n', ann=FALSE)

plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',
   # xlab='Temproal Sensitivity (Per 50mm)',
   xlab = expression(paste(beta ["temporal"])),
   yaxt='n', 
   # ann=FALSE
   cex.lab = 1.5, cex.axis = 1.5)
   

#draw a box around the plot
box()
abline(v = 0)

# high prep : MAP > 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 2, x1 = x1[4], y1 = 2, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 2, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 1.5, "MAP > 750 mm")
# paste("MAP > 750 mm ( n = ", length(which(dat1$lc == 2)), " )", sep = ""),  cex = 1)
  

# low prep : : MAP < 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, "MAP < 750 mm")
# paste("MAP < 750 mm ( n = ", length(which(dat1$lc == 1)), " )", sep = ""),  cex = 1)

################## add ER         #################################
# high prep : MAP > 750 mm
x11 = dat1$er.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 1, x1 = x1[4], y1 = 1, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 1, type = "p", pch = 17, cex = 2, col = "red")

# low prep : : MAP < 750 mm
x11 = dat1$er.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 3, x1 = x1[4], y1 = 3, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 3, type = "p", pch = 17, cex = 2, col = "red")

	  },
	  
		  mar=c(5,0,0,0))

####################### END of bootstrap methods       ##############################


##############################       Mento Carlo Method      ##############################

gpp.df.z1 = cbind(gpp.df.z, temp.df.z, prep.df.z)
er.df.z1 = cbind(er.df.z, temp.df.z, prep.df.z)

# for temporal sensitivity
gpp.beta.temporal = list() # each element contains one region
er.beta.temporal = list()

for (i in 1:nrow(gpp.df.z)){
## calcualte GPP

# select regions
dat11 = data.frame(gpp = gpp.df.z[i,], temp = temp.df.z[i,], prep = prep.df.z[i,])

coef1 = c()
N = 0
repeat{
N = N+1
# select number of element to removed random
n1 = sample(c(1:3),1)
#select index of the element to be used for regression
l1 = nrow(dat11) - n1
dat12 = dat11[sample(c(1:nrow(dat11)), l1), ]
lm1 = lm(gpp ~ temp + prep, data = dat12)
coef1 = rbind(coef1, coef(lm1))
if(N >= 100){break}

	}
	
gpp.beta.temporal[[i]] = apply(coef1, 2, cifun)
	
## calcualte ER

# select regions
dat11 = data.frame(er = er.df.z[i,], temp = temp.df.z[i,], prep = prep.df.z[i,])

coef1 = c()
N = 0
repeat{
N = N+1
# select number of element to removed random
n1 = sample(c(1:3),1)
#select index of the element to be used for regression
l1 = nrow(dat11) - n1
dat12 = dat11[sample(c(1:nrow(dat11)), l1), ]
lm1 = lm(er ~ temp + prep, data = dat12)
coef1 = rbind(coef1, coef(lm1))
if(N >= 100){break}

	}
	
er.beta.temporal[[i]] = apply(coef1, 2, cifun)

print(paste("Finish calculating for region ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

gpp.beta.temporal2 = c()
er.beta.temporal2 = c()

for(i in 1:length(gpp.beta.temporal)){

gpp.beta.temporal2 = rbind(gpp.beta.temporal2, gpp.beta.temporal[[i]][2,])
er.beta.temporal2 = rbind(er.beta.temporal2, er.beta.temporal[[i]][2,])

}

dat1 = data.frame(prep = prep.df.mn[-4],gpp.p = gpp.beta.temporal2[-4,3],er.p = er.beta.temporal2[-4,3], 
				  temp = temp.df.mn[-4],gpp.t = gpp.beta.temporal2[-4,2],er.t = er.beta.temporal2[-4,2])
dat1$lc = 1
dat1$lc[dat1$prep > 750] = 2


#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.mn[-4], nee = nee.df.mn[-4],er = er.df.mn[-4],airtemp = temp.df.mn[-4],prep = prep.df.mn[-4])


# step 1: fit model
lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)

# step 2: get residual
resd1.gpp = lm1.gpp$residuals
resd1.mean.gpp = mean(resd1.gpp)
resd1.sd.gpp = sd(resd1.gpp)

resd1.er = lm1.er$residuals
resd1.mean.er = mean(resd1.er)
resd1.sd.er = sd(resd1.er)

reg.abs.mod.pred1.er.list2 = list()
reg.abs.mod.pred1.gpp.list2 = list()

delt.er = c()
delt.gpp = c()

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

## calcualte GPP
N = 0
repeat{
N = N+1

resd1.value = rnorm(nrow(dat.df.annual.mean), mean = resd1.mean.gpp, sd = resd1.sd.gpp)
dat.df.annual.mean$gpp1 = dat.df.annual.mean$gpp + resd1.value
lm1.gpp1 = lm(gpp1~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)

resd1.value = rnorm(nrow(dat.df.annual.mean), mean = resd1.mean.er, sd = resd1.sd.er)
dat.df.annual.mean$er1 = dat.df.annual.mean$er + resd1.value
lm1.er1 = lm(er1~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)

#construct new data for prediction
pred1 = predict(lm1.er1, newdata = dat.new1)
pred2 = predict(lm1.er1, newdata = dat.new2)
reg.abs.mod.pred1.er <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)
pred1 = predict(lm1.gpp1, newdata = dat.new1)
pred2 = predict(lm1.gpp1, newdata = dat.new2)
reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

delt.er1 = diff(reg.abs.mod.pred1.er$predp)
delt.gpp1 = diff(reg.abs.mod.pred1.gpp$predp)

reg.abs.mod.pred1.er.list2[[N]] <- reg.abs.mod.pred1.er
reg.abs.mod.pred1.gpp.list2[[N]] <- reg.abs.mod.pred1.gpp

delt.er = cbind(delt.er, delt.er1)
delt.gpp = cbind(delt.gpp, delt.gpp1)

if(N >= 100){break}

}

pred.er = c()
pred.gpp = c()

for (i in 1:length(reg.abs.mod.pred1.er.list2)){
pred.er = cbind(pred.er, reg.abs.mod.pred1.er.list2[[i]][,4])
pred.gpp = cbind(pred.gpp, reg.abs.mod.pred1.gpp.list2[[i]][,4])
}


pred.er2 = apply(pred.er, 1, cifun, ALPHA = 0.01)
pred.gpp2 = apply(pred.gpp, 1, cifun, ALPHA = 0.01)

delt.gpp2 = apply(delt.gpp, 1, cifun, ALPHA = 0.01)
delt.er2 = apply(delt.er, 1, cifun, ALPHA = 0.01)

prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(-0.2,2000),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)

### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 1, col = "green")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.gpp2[3,], lty = 'dashed', col = 'green')
lines(prep.grd,pred.gpp2[1,], lty = 'dashed', col = 'green')
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 1, col = "green")
points(dat.df.annual.mean$prep, dat.df.annual.mean$gpp, col = "green")


### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.er2[3,], lty = 'dashed', col = 'red')
lines(prep.grd,pred.er2[1,], lty = 'dashed', col = 'red')
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
points(dat.df.annual.mean$prep, dat.df.annual.mean$er, col = "red")

# plot spatial sensitivity
library(oce)
plotInset(750,0, 1450, 800,
          expr= {
		  
# plot(0, xlim = c(0, 1500), ylim = c(0,100), pch='', xlab = "Precipation (mm)", ylab = "g C m-2 yr-1", cex.lab = 1.5, cex.axis = 1.5)
plot(0, xlim = c(0, 1500), ylim = c(25,75)*2, pch='', xlab = "Precipation (mm)", 
     # ylab = "Spatial Sensitivty (Per 50mm)", 
	 ylab = expression(paste(beta ["spatial"])),
	 cex.lab = 1.5, cex.axis = 1.5)

#draw a box around the plot
box()
#######  GPP
lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# add fill
# polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
        # c(rev(delt.gpp2[3,]*2), delt.gpp2[1,]*2), col = 'grey80', border = NA)
# lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# intervals
lines(prep.grd[-1],delt.gpp2[3,]*2, lty = 'dashed', col = 'green')
lines(prep.grd[-1],delt.gpp2[1,]*2, lty = 'dashed', col = 'green')

#######  ER
lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# add fill
# polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
        # c(rev(delt.er2[3,]*2), delt.er2[1,]*2), col = 'grey80', border = NA)
# lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# intervals
lines(prep.grd[-1],delt.er2[3,]*2, lty = 'dashed', col = 'red')
lines(prep.grd[-1],delt.er2[1,]*2, lty = 'dashed', col = 'red')

abline(h = 0)
# abline(v = reg.abs.mod.pred1.gpp[-1,3][which.max(which(delt.gpp > delt.er))])
abline(v = 750)

	  },	  
		  mar=c(5,0,0,0))
	
# writeRaster(nee.annual, filename="carbontrackerNEE.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
# writeRaster(gpp.annual, filename="mod17gpp.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

#plot temporal sensitivity
		  
library(oce)
plotInset(200, 700, 1000, 2000,
          expr= {
		 		  
#plot temproal sensitivity
# plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',xlab='Temproal Sensitivity (Per 50mm)',yaxt='n', ann=FALSE)

plot(0, xlim = c(-100, 100), ylim = c(0,5),bty='n',pch='',ylab='',
   # xlab='Temproal Sensitivity (Per 50mm)',
   xlab = expression(paste(beta ["temporal"])),
   yaxt='n', 
   # ann=FALSE
   cex.lab = 1.5, cex.axis = 1.5)
   

#draw a box around the plot
box()
abline(v = 0)

# high prep : MAP > 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 2, x1 = x1[4], y1 = 2, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 2, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 1.5, "MAP > 750 mm")
# paste("MAP > 750 mm ( n = ", length(which(dat1$lc == 2)), " )", sep = ""),  cex = 1)
  

# low prep : : MAP < 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, "MAP < 750 mm")
# paste("MAP < 750 mm ( n = ", length(which(dat1$lc == 1)), " )", sep = ""),  cex = 1)

################## add ER         #################################
# high prep : MAP > 750 mm
x11 = dat1$er.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 1, x1 = x1[4], y1 = 1, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 1, type = "p", pch = 17, cex = 2, col = "red")

# low prep : : MAP < 750 mm
x11 = dat1$er.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 3, x1 = x1[4], y1 = 3, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 3, type = "p", pch = 17, cex = 2, col = "red")

	  },
	  
		  mar=c(5,0,0,0))

##############################  END of Mento Carlo Method      ##############################








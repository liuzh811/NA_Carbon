## using LPJ model to calculate the relationship between GPP/ER with precipitation and temperature
## at 0.5 degree resolution
## using the median/mean NEE estimate from four co2 inversion methods

library(rgdal)
library(raster)
library(rasterVis)
require(ncdf4)

usa.state = readOGR(dsn="F:\\zhihua\\dataset\\ecoregion", layer = "usa.state")
usa.state = usa.state[-which(usa.state$NAME_1 == "Hawaii" | usa.state$NAME_1 == "Alaska"),]

######### section 1: read into LPJ simulations #########################################
na.ext <- extent(-180,-48,15,85)

##look at the lpj models
## nee.annual in the unit of Kg C m-1 yr-1
nee.annual.lpj = stack("F:/zhihua/dataset/lpj/processed/LPJ_nbp.2000.2014.annual.grd")
nee.annual.lpj = crop(nee.annual.lpj, usa.state)
nee.annual.lpj = nee.annual.lpj*1000 # change to g C m-1 yr-1 to compare with Carbon Tracker

# read into ER
ra.annual = stack("F:/zhihua/dataset/lpj/processed/LPJ_ra..2000.2014.annual.grd")
ra.annual = crop(ra.annual, na.state)
ra.annual = ra.annual*1000 # change to g C m-2 yr-1 

rh.annual = stack("F:/zhihua/dataset/lpj/processed/LPJ_rh..2000.2014.annual.grd")
rh.annual = crop(rh.annual, na.state)
rh.annual = rh.annual*1000 # change to g C m-2 yr-1 


# read into GPP
gpp.annual.lpj = stack("F:/zhihua/dataset/lpj/processed/LPJ_gpp.2000.2014.annual.grd")
gpp.annual.lpj = crop(gpp.annual.lpj, na.state)
gpp.annual.lpj = gpp.annual.lpj*1000 # change to g C m-2 yr-1 

# gpp.annual.lpj = gpp.annual.lpj*grid.area.gpp
# er.annual = er.annual*grid.area.gpp
# nee.annual.lpj = nee.annual.lpj*grid.area.gpp
# rh.annual = rh.annual*grid.area.gpp
# ra.annual = ra.annual*grid.area.gpp

######### section 2: read into various grid and polygon #########################################
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

prep1 = stack(paste0("F:/zhihua/dataset/cru_ts3.23/prep", yr, ".grd"))
prep1 = crop(prep1, na.ext)
prep1 = calc(prep1, sum, na.rm = TRUE)

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
gpp.annual.lpj = gpp.annual.lpj*mask
nee.annual.lpj = nee.annual.lpj*mask

ra.annual = ra.annual*mask
rh.annual = rh.annual*mask
er.annual = rh.annual+ra.annual

temp = temp*mask
prep = prep*mask

# change NEON into grid
noen.grd.gpp = rasterize(neon.sp,mask, field = "DomainID")				

###########  section 3, calculating the correlationship between GPP/ER and prep/temp

nee.df.lpj = zonal(nee.annual.lpj, noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]
gpp.df.lpj = zonal(gpp.annual.lpj, noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]
er.df.lpj = zonal(er.annual, noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]

rh.df.lpj = zonal(rh.annual, noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]
ra.df.lpj = zonal(ra.annual, noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]

nee.df.lpj.mn = zonal(calc(nee.annual.lpj, mean),noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,2]
gpp.df.lpj.mn = zonal(calc(gpp.annual.lpj, mean),noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,2]
er.df.lpj.mn = zonal(calc(er.annual, mean),noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,2] 
rh.df.lpj.mn = zonal(calc(rh.annual, mean),noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,2] 
ra.df.lpj.mn = zonal(calc(ra.annual, mean),noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,2] 

temp.df = zonal(temp, noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]
prep.df = zonal(prep, noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,c(2:16)]
temp.df.mn = zonal(calc(temp, mean),noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,2]
prep.df.mn = zonal(calc(prep, mean),noen.grd.gpp, fun='mean', digits=0, na.rm=TRUE)[,2]


#find the dominant land cover within each neon domain
lc_us2_rc3 = raster("C:/zhihua/dataset/ecoregion/LCTypeUS_rc3_0.5_2.tif")
newclas_names <- c("EvergreenForest", "BroadForest","Shrubland","Grass", "Crop")
lc_us2_rc3[lc_us2_rc3 == 0] = NA
extent(lc_us2_rc3) <- extent(mask.gpp)
lc_us2_rc3 = lc_us2_rc3*mask.gpp
lc_us2_rc3.df = data.frame(zonal(lc_us2_rc3, noen.grd.gpp1, fun='modal', digits=0, na.rm=TRUE))
lc_us2_rc3.df[2,2] = 2
lc_us2_rc3.df[4,2] = 2
lc_us2_rc3.df[6,2] = 2
lc_us2_rc3.df[7,2] = 2
lc_us2_rc3.df[11,2] = 1
lc_us2_rc3.df$modal = factor(lc_us2_rc3.df$modal)
levels(lc_us2_rc3.df$modal) <- c("ENF","DBF","Nonforest","Nonforest","Nonforest")

## temproal correlationship
## clculate Z-score
Zscore = function(x){(x-mean(x))/sd(x)} 
nee.df.lpj.z = t(apply(nee.df.lpj, 1, Zscore))
gpp.df.lpj.z = t(apply(gpp.df.lpj, 1, Zscore))
er.df.lpj.z = t(apply(er.df.lpj, 1, Zscore))
rh.df.lpj.z = t(apply(rh.df.lpj, 1, Zscore))
ra.df.lpj.z = t(apply(ra.df.lpj, 1, Zscore))

##use anomaly
nee.df.lpj.z = t(apply(nee.df.lpj, 1, Ano))
gpp.df.lpj.z = t(apply(gpp.df.lpj, 1, Ano))
er.df.lpj.z = t(apply(er.df.lpj, 1, Ano))
rh.df.lpj.z = t(apply(rh.df.lpj, 1, Ano))
ra.df.lpj.z = t(apply(ra.df.lpj, 1, Ano))

prep.df.z = t(apply(prep.df, 1, Ano))
temp.df.z = t(apply(temp.df, 1, Ano))

#plot coefficient, use as temproal sensitivity
# need to remove region 4

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

gpp.df.lpj.z = cbind(gpp.df.lpj.z, temp.df.z, prep.df.z)
er.df.lpj.z = cbind(er.df.lpj.z, temp.df.z, prep.df.z)
rh.df.lpj.z = cbind(rh.df.lpj.z, temp.df.z, prep.df.z)
ra.df.lpj.z = cbind(ra.df.lpj.z, temp.df.z, prep.df.z)

gpp.df.lpj.z.r = apply(gpp.df.lpj.z, 1, fun.cor2) 
er.df.lpj.z.r = apply(er.df.lpj.z, 1, fun.cor2) 
rh.df.lpj.z.r = apply(rh.df.lpj.z, 1, fun.cor2) 
ra.df.lpj.z.r = apply(ra.df.lpj.z, 1, fun.cor2) 

mean(gpp.df.lpj.z.r[5,]) # 0.43 +-0.19
mean(er.df.lpj.z.r[5,]) # 0.50 +- 0.27
mean(rh.df.lpj.z.r[5,]) # 0.56 +- 0.21
mean(ra.df.lpj.z.r[5,]) # 0.50 +- 0.33


dat1 = data.frame(prep = prep.df.mn[-4],gpp.p = gpp.df.lpj.z.r[2,-4],er.p = er.df.lpj.z.r[2,-4], rh.p = rh.df.lpj.z.r[2,-4],ra.p = ra.df.lpj.z.r[2,-4],
				  temp = temp.df.mn[-4],gpp.t = gpp.df.lpj.z.r[1,-4],er.t = er.df.lpj.z.r[1,-4], ra.t = ra.df.lpj.z.r[2,-4],ra.t = ra.df.lpj.z.r[2,-4])
dat1$p.dif = dat1$gpp.p - dat1$er.p
dat1$lc = 1
dat1$lc[dat1$prep > 750] = 2

#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.lpj.mn[-4], nee = nee.df.lpj.mn[-4],er = er.df.lpj.mn[-4],
								airtemp = temp.df.mn[-4],prep = prep.df.mn[-4],rh = rh.df.lpj.mn[-4],ra = ra.df.lpj.mn[-4])

#construct new data for prediction
dat.new1 = data.frame(airtemp = seq(1,23, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.er = stepAIC(lm1.er)

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

lm1.rh = lm(rh~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.rh = stepAIC(lm1.rh)
pred1 = predict(lm1.rh, newdata = dat.new1)
pred2 = predict(lm1.rh, newdata = dat.new2)
reg.abs.mod.pred1.rh <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

lm1.ra = lm(ra~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.ra = stepAIC(lm1.ra)
pred1 = predict(lm1.ra, newdata = dat.new1)
pred2 = predict(lm1.ra, newdata = dat.new2)
reg.abs.mod.pred1.ra <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

# png("C:/zhihua/dataset/results/er.prep01032017.lpj.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2500),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)

lines(smooth.spline(reg.abs.mod.pred1.gpp[,3], reg.abs.mod.pred1.gpp[,4]), lty = 1, lwd = 6, col = "green")
points(dat.df.annual.mean$prep,dat.df.annual.mean$gpp, cex = 1, col = "green")

lines(smooth.spline(reg.abs.mod.pred1.er[,3], reg.abs.mod.pred1.er[,4]), lty = 1, lwd = 6, col = "red")
points(dat.df.annual.mean$prep,dat.df.annual.mean$er, cex = 1, col = "red")

lines(smooth.spline(reg.abs.mod.pred1.ra[,3], reg.abs.mod.pred1.ra[,4]), lty = 1, lwd = 6, col = "blue")
lines(smooth.spline(reg.abs.mod.pred1.rh[,3], reg.abs.mod.pred1.rh[,4]), lty = 2, lwd = 6, col = "blue")


lines(smooth.spline(reg.abs.mod.pred1.nee[,3], reg.abs.mod.pred1.nee[,4]), lty = 1, lwd = 6, col = "red")
lines(smooth.spline(reg.abs.mod.pred1.nee[,3],reg.abs.mod.pred1.gpp[,4]/reg.abs.mod.pred1.er[,4]), lty = 1, lwd = 6, col = "red")

delt.er = diff(reg.abs.mod.pred1.er$predp)
delt.gpp = diff(reg.abs.mod.pred1.gpp$predp)
delt.ra = diff(reg.abs.mod.pred1.ra$predp)
delt.rh = diff(reg.abs.mod.pred1.rh$predp)

library(oce)
plotInset(900, 0, 1500, 900,
          expr= {
		  
# plot(0, xlim = c(0, 1500), ylim = c(-0.2,0.2), pch='', xlab = "Precipation (mm)", ylab = "Spatial Sensitivty (Per 50mm)", cex.lab = 1.5, cex.axis = 1.5)
plot(0, xlim = c(0, 1500), ylim = c(-5,160)*2, pch='', xlab = "Precipation (mm)", 
		# ylab = "g C m-2 yr-1", 
		ylab = expression(paste(beta ["spatial"])),
		cex.lab = 1.5, cex.axis = 1.5)

#draw a box around the plot
box()
lines(smooth.spline(reg.abs.mod.pred1.gpp[-1,3], delt.gpp*2), lty = 1, lwd = 6, col = "green")
lines(smooth.spline(reg.abs.mod.pred1.er[-1,3], delt.er*2), lty = 1, lwd = 6, col = "red")
# lines(smooth.spline(reg.abs.mod.pred1.ra[-1,3], delt.ra*2), lty = 1, lwd = 6, col = "blue")
# lines(smooth.spline(reg.abs.mod.pred1.rh[-1,3], delt.rh*2), lty = 2, lwd = 6, col = "blue")

abline(h = 0)
abline(v = reg.abs.mod.pred1.gpp[-1,3][which.max(which(delt.gpp > delt.er))])

 },
		  mar=c(2.5,0,0,0))

##plot temproal sensitivity
	  
library(oce)
plotInset(50, 1200, 650, 2400,
          expr= {
		 		  
#plot temproal sensitivity
# plot(0, xlim = c(-100, 200), ylim = c(0,5),bty='n',pch='',ylab='',xlab='Temproal Sensitivity (Per 50mm)',yaxt='n', ann=FALSE)
plot(0, xlim = c(-100, 200), ylim = c(0,5),bty='n',pch='',ylab='',
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
# paste("MAP > 750 mm ( n = ", length(which(dat1$lc == 2)), " )", sep = ""), cex = 1)
  

# low prep : : MAP < 750 mm
x11 = dat1$gpp.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, "MAP < 750 mm")
# paste("MAP < 750 mm ( n = ", length(which(dat1$lc == 1)), " )", sep = ""), cex = 1)

################## add ER         #################################
# high prep : MAP > 750 mm
x11 = dat1$er.p[which(dat1$lc == 2)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 1, x1 = x1[4], y1 = 1, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 1, type = "p", pch = 17, cex = 2, col = "red")

#ra
# x11 = dat1$ra.p[which(dat1$lc == 2)]*100
# x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
# segments(x0 = x1[2], y0 = 0.75, x1 = x1[4], y1 = 0.75, lty = 1, lwd = 2, col = "blue")
# points(x = mean(x11), y = 0.75, type = "p", pch = 17, cex = 2, col = "blue")
# rh
# x11 = dat1$rh.p[which(dat1$lc == 2)]*100
# x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
# segments(x0 = x1[2], y0 = 0.5, x1 = x1[4], y1 = 0.5, lty = 1, lwd = 2, col = "blue")
# points(x = mean(x11), y = 0.5, type = "p", pch = 15, cex = 2, col = "blue")

# low prep : : MAP < 750 mm
x11 = dat1$er.p[which(dat1$lc == 1)]*100
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 3, x1 = x1[4], y1 = 3, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 3, type = "p", pch = 17, cex = 2, col = "red")

# ra
# x11 = dat1$ra.p[which(dat1$lc == 1)]*100
# x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
# segments(x0 = x1[2], y0 = 2.75, x1 = x1[4], y1 = 2.75, lty = 1, lwd = 2, col = "blue")
# points(x = mean(x11), y = 2.75, type = "p", pch = 17, cex = 2, col = "blue")
# rh
# x11 = dat1$rh.p[which(dat1$lc == 1)]*100
# x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
# segments(x0 = x1[2], y0 = 2.5, x1 = x1[4], y1 = 2.5, lty = 1, lwd = 2, col = "blue")
# points(x = mean(x11), y = 2.5, type = "p", pch = 15, cex = 2, col = "blue")

	  },
	  
		  mar=c(3,0,0,0))

# dev.off()
	
############################## leave-one-out methods ###################################

#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.lpj.mn[-4], nee = nee.df.lpj.mn[-4],er = er.df.lpj.mn[-4],
								airtemp = temp.df.mn[-4],prep = prep.df.mn[-4],rh = rh.df.lpj.mn[-4],ra = ra.df.lpj.mn[-4])

reg.abs.mod.pred1.er.list = list()
reg.abs.mod.pred1.gpp.list = list()

delt.er = c()
delt.gpp = c()

reg.abs.mod.pred1.gpp.model.coef = c() #store model coefficients
reg.abs.mod.pred1.er.model.coef = c() #store model coefficients

rmse.gpp = c() #store model efficients, including rmse, mae, aic, r-squired
rmse.er = c() #store model efficients, including rmse, mae, aic, r-squired


for (i in 1:nrow(dat.df.annual.mean)){

dat.df.annual.mean1 = dat.df.annual.mean[-i,]

lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)
lm1.er = stepAIC(lm1.er)

rmse <- round(sqrt(mean(resid(lm1.er)^2)), 2)
mae <- round(mean(abs(resid(lm1.er))), 2)
aic = AIC(lm1.er)
r2 = round(summary(lm1.er)$r.squared, 2)
rmse.er = rbind(rmse.er, c(rmse, mae,aic,r2))

reg.abs.mod.pred1.er.model.coef = rbind(reg.abs.mod.pred1.gpp.model.coef, coef(lm1.gpp))


dat.new1 = data.frame(airtemp = seq(1,23, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

#construct new data for prediction
pred1 = predict(lm1.er, newdata = dat.new1)
pred2 = predict(lm1.er, newdata = dat.new2)
reg.abs.mod.pred1.er <- data.frame(airtemp = seq(1,23, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)
lm1.gpp = stepAIC(lm1.gpp)
rmse <- round(sqrt(mean(resid(lm1.gpp)^2)), 2)
mae <- round(mean(abs(resid(lm1.gpp))), 2)
aic = AIC(lm1.gpp)
r2 = round(summary(lm1.gpp)$r.squared, 2)
rmse.gpp = rbind(rmse.gpp, c(rmse, mae,aic,r2))

reg.abs.mod.pred1.gpp.model.coef = rbind(reg.abs.mod.pred1.gpp.model.coef, coef(lm1.gpp))


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

png("F:/zhihua/dataset/results2/lpj.sensitivity_ci-2.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,6,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), 
		ylim = c(-0.2,2000),bty='n',pch='',
		 xlab='Precipation (mm)',
			ylab=expression("GPP/TER" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")), 
		cex.axis = 2, cex.lab = 2)

### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.gpp2[3,], lty = 'dashed', col = 'green')
lines(prep.grd,pred.gpp2[1,], lty = 'dashed', col = 'green')
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")
# points(dat.df.annual.mean$prep, dat.df.annual.mean$gpp, col = "green")

### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.er2[3,], lty = 'dashed', col = 'red')
lines(prep.grd,pred.er2[1,], lty = 'dashed', col = 'red')
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# points(dat.df.annual.mean$prep, dat.df.annual.mean$er, col = "red")

# plot spatial sensitivity
library(oce)
plotInset(750,0, 1450, 800,
          expr= {
		  
plot(0, xlim = c(0, 1500), ylim = c(0,200)*2, pch='', xlab = "Precipation (mm)", 
     # ylab = "Spatial Sensitivty (Per 50mm)", 
	 ylab = expression(paste(beta ["spatial"])),
	 cex.lab = 1.5, cex.axis = 1.5)

#draw a box around the plot
box()
#######  GPP
lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# add fill
 polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
         c(rev(delt.gpp2[3,]*2), delt.gpp2[1,]*2), col = 'grey90', border = NA)
 lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "green")
# intervals
lines(prep.grd[-1],delt.gpp2[3,]*2, lty = 'dashed', col = 'green')
lines(prep.grd[-1],delt.gpp2[1,]*2, lty = 'dashed', col = 'green')

#######  ER
lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# add fill
 polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
         c(rev(delt.er2[3,]*2), delt.er2[1,]*2), col = 'grey90', border = NA)
 lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# intervals
lines(prep.grd[-1],delt.er2[3,]*2, lty = 'dashed', col = 'red')
lines(prep.grd[-1],delt.er2[1,]*2, lty = 'dashed', col = 'red')

# abline(v = c(750-100, 750, 750 + 100), lty = c(2,1,2))

	  },	  
		  mar=c(5,0,0,0))

# temproal sensitivity
x11 = dat1$gpp.p[which(dat1$lc == 2)]*100
x12 = dat1$gpp.p[which(dat1$lc == 1)]*100
x13 = dat1$er.p[which(dat1$lc == 2)]*100
x14 = dat1$er.p[which(dat1$lc == 1)]*100
 

se <- function(x) sqrt(var(x)/length(x))

D1 = data.frame(mean = c(mean(x11),mean(x12),mean(x13),mean(x14)),
				se = c(se(x11),se(x12),se(x13),se(x14)),
				flux = c("GPP","GPP","TER","TER"),
				region = c("> 750 mm", "< 750 mm","> 750 mm","< 750 mm"),
				type = rep("temproal",4))
	

ylab = 	expression(paste(beta ["temporal"]) ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "/100mm"))
color1 = c("green", "red")

p1 <- ggplot(data=D1, aes(x=region, y=mean, fill=flux)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 	
	# ylab(ylab) + 	 
    ylab(expression(paste(beta ["temporal"]))) + 
	xlab("Site by MAP") + # Set axis labels
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
	theme(legend.position=c(.5, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 12)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=12),axis.text.x  = element_text(colour="black",size=10))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=12),axis.text.y  = element_text(colour="black",size=12))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) +
    scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(D1$flux),
                    labels=levels(D1$flux)) +
	# geom_text(aes(label = "2"), vjust = "inward", hjust = "inward") + 
	geom_text(x = 0.5, y = 70, label = "2)",size=6)
	  
	  
library(ggplot2)
require(grid)
print(p1, vp=viewport(.35, .8, .3, .3)) #viewport, first two is the x,y coor; last 2 is the inset size

dev.off()

############################## END of leave-one-out methods ###################################

w.idx = which(prep.grd < 750)
e.idx = which(prep.grd > 750)

delt = data.frame(mean = c(mean(delt.gpp[e.idx-1,]),mean(delt.gpp[w.idx,]),mean(delt.er[e.idx-1,]),mean(delt.er[w.idx,])),
					se = c(sd(delt.gpp[e.idx-1,]),sd(delt.gpp[w.idx,]),sd(delt.er[e.idx-1,]),sd(delt.er[w.idx,])),
				    flux = c("GPP","GPP","TER","TER"),
				    region = c("> 750 mm", "< 750 mm","> 750 mm","< 750 mm"),
					type = rep("spatial",4))

D3 = rbind(delt, D1)			
write.csv(D3, "F:/zhihua/dataset/results2/lpj.sensitivity.csv")				

p1 <- ggplot(data=D3, aes(x=region, y=mean, fill=flux)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 
    facet_grid(.~type)+			 
	# ylab(ylab) + 	 
    ylab(expression(paste(beta ["temporal"]))) + 
	xlab("Site by MAP") + # Set axis labels
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
	theme(legend.position=c(.5, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 12)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=12),axis.text.x  = element_text(colour="black",size=10))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=12),axis.text.y  = element_text(colour="black",size=12))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) +
    scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(D1$flux),
                    labels=levels(D1$flux)) +
	# geom_text(aes(label = "2"), vjust = "inward", hjust = "inward") + 
geom_text(x = 0.5, y = 70, label = "2)",size=6)					

colnames(rmse.er) = c("rmse","mae","aic","r2") 
colnames(rmse.gpp) = c("rmse","mae","aic","r2") 
rmse1 = rbind(data.frame(rmse.gpp, flux = "gpp"),data.frame(rmse.er, flux = "er"))
write.csv(rmse1, "F:/zhihua/dataset/results2/lpj.rmse1.csv")

colnames(reg.abs.mod.pred1.er.model.coef) = c("int","p2","p") 
colnames(reg.abs.mod.pred1.gpp.model.coef) = c("int","p2","p") 
coef1 = rbind(data.frame(reg.abs.mod.pred1.gpp.model.coef, flux = "gpp"),data.frame(reg.abs.mod.pred1.er.model.coef, flux = "er"))
write.csv(coef1, "F:/zhihua/dataset/results2/lpj.coef1.csv")


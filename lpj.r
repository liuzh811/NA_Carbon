## 5/22/2017
## using LPJ model to calculate the relationship between GPP/ER with precipitation and temperature



##################################################
##look at the lpj models
## nee.annual in the unit of Kg C m-1 yr-1
nee.annual.lpj = stack("D:/zhihua/dataset/lpj/processed/LPJ_nbp.2000.2014.annual.grd")
nee.annual.lpj = crop(nee.annual.lpj, usa.state)
nee.annual.lpj = nee.annual.lpj*1000 # change to g C m-1 yr-1 to compare with Carbon Tracker
nee.annual.lpj = nee.annual.lpj*mask.gpp

# read into ER
ra.annual = stack("D:/zhihua/dataset/lpj/processed/LPJ_ra..2000.2014.annual.grd")
ra.annual = crop(ra.annual, na.state)
ra.annual = ra.annual*1000 # change to g C m-2 yr-1 
ra.annual = ra.annual*mask.gpp

rh.annual = stack("D:/zhihua/dataset/lpj/processed/LPJ_rh..2000.2014.annual.grd")
rh.annual = crop(rh.annual, na.state)
rh.annual = rh.annual*1000 # change to g C m-2 yr-1 
rh.annual = rh.annual*mask.gpp
er.annual = rh.annual+ra.annual

# read into GPP
gpp.annual.lpj = stack("D:/zhihua/dataset/lpj/processed/LPJ_gpp.2000.2014.annual.grd")
gpp.annual.lpj = crop(gpp.annual.lpj, na.state)
gpp.annual.lpj = gpp.annual.lpj*1000 # change to g C m-2 yr-1 
gpp.annual.lpj = gpp.annual.lpj*mask.gpp

# gpp.annual.lpj = gpp.annual.lpj*grid.area.gpp
# er.annual = er.annual*grid.area.gpp
# nee.annual.lpj = nee.annual.lpj*grid.area.gpp
# rh.annual = rh.annual*grid.area.gpp
# ra.annual = ra.annual*grid.area.gpp


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
lculate Z-score
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

png("C:/zhihua/dataset/results/er.prep01032017.lpj.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2500),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)

lines(smooth.spline(reg.abs.mod.pred1.gpp[,3], reg.abs.mod.pred1.gpp[,4]), lty = 1, lwd = 6, col = "green")
points(dat.df.annual.mean$prep,dat.df.annual.mean$gpp, cex = 1, col = "green")

lines(smooth.spline(reg.abs.mod.pred1.er[,3], reg.abs.mod.pred1.er[,4]), lty = 1, lwd = 6, col = "red")
points(dat.df.annual.mean$prep,dat.df.annual.mean$er, cex = 1, col = "red")

# lines(smooth.spline(reg.abs.mod.pred1.ra[,3], reg.abs.mod.pred1.ra[,4]), lty = 1, lwd = 6, col = "blue")
# lines(smooth.spline(reg.abs.mod.pred1.rh[,3], reg.abs.mod.pred1.rh[,4]), lty = 2, lwd = 6, col = "blue")


# lines(smooth.spline(reg.abs.mod.pred1.nee[,3], reg.abs.mod.pred1.nee[,4]), lty = 1, lwd = 6, col = "red")
# lines(smooth.spline(reg.abs.mod.pred1.nee[,3],reg.abs.mod.pred1.gpp[,4]/reg.abs.mod.pred1.er[,4]), lty = 1, lwd = 6, col = "red")

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

dev.off()
		  
## plot gpp/er sensitive to temperature 
delt.er2 = diff(reg.abs.mod.pred1.er$predt)
delt.gpp2 = diff(reg.abs.mod.pred1.gpp$predt)
plot(0, xlim = c(1, 25), ylim = c(0,2500), pch='', xlab = "Temperature (degree)", ylab = "Spatial Sensitivty (Per 50mm)", cex.lab = 1.5, cex.axis = 1.5)

lines(smooth.spline(reg.abs.mod.pred1.gpp[,1], reg.abs.mod.pred1.gpp[,2]), lty = 1, lwd = 6, col = "green")
points(dat.df.annual.mean$temp,dat.df.annual.mean$gpp, cex = 1, col = "green")

lines(smooth.spline(reg.abs.mod.pred1.er[,1], reg.abs.mod.pred1.er[,2]), lty = 1, lwd = 6, col = "red")
points(dat.df.annual.mean$prep,dat.df.annual.mean$er, cex = 1, col = "red")


#draw a box around the plot
box()
lines(smooth.spline(reg.abs.mod.pred1.gpp[-1,1], delt.gpp2), lty = 1, lwd = 6, col = "green")
lines(smooth.spline(reg.abs.mod.pred1.er[-1,1], delt.er2), lty = 1, lwd = 6, col = "red")

abline(h = 0)
abline(v = reg.abs.mod.pred1.gpp[-1,1][which.max(which(delt.gpp2 > delt.er2))+1])

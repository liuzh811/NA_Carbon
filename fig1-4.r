## plot fig1-4 and fig1s

Work.Dir <- "F:/zhihua/dataset/fpar"
setwd(Work.Dir)


#load library
library(rgdal)
library(raster)
library(rasterVis)

na.state = readOGR(dsn="F:\\zhihua\\dataset\\ecoregion", layer = "na.state2")
na.state = na.state[-which(na.state$NAME_1 == "Prince Edward Island"|
                           na.state$NAME_1 == "Hawaii"),]

na.ext <- extent(-180,-48,15,85)

usa.state = na.state[which(na.state$ISO == "USA" & na.state$NAME_1 != "Alaska"), ]

##################################################################		
######plot ensemble mean for 8 maps :::: RS
R.stack1 = stack("F:/zhihua/dataset/ct/cor.stack.ct.grd")
P.stack1 = stack("F:/zhihua/dataset/ct/P.stack.ct.grd")

# R.stack2 = stack("C:/zhihua/dataset/cte/cor.stack.cte.grd")
# P.stack2 = stack("C:/zhihua/dataset/cte/P.stack.cte.grd")

R.stack2 = stack("F:/zhihua/dataset/ct/cor.stack.ecmod.grd")
P.stack2 = stack("F:/zhihua/dataset/ct/p.stack.ecmod.grd")
	
R = calc(stack(R.stack1, R.stack2), mean, na.rm = TRUE)	
P = calc(stack(P.stack1, P.stack2), function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 4] = 1
P1[P > 4] = 0.05
pts.sp.sig2 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points

#plot individually

R.stack = stack(R.stack1, R.stack2)
P.stack = stack(P.stack1, P.stack2)

##################################################################		
######plot ensemble mean for 10 maps :::: TRENDY model

R.trendy = stack("F:/zhihua/dataset/trendy/R.trendy.gpp.nee.stack.grd")
P.trendy = stack("F:/zhihua/dataset/trendy/P.trendy.gpp.nee.stack.grd")
	
R.trendy1 = calc(R.trendy, mean, na.rm = TRUE)	
P.trendy1 = calc(P.trendy, function(x){length(which(x <= 0.1))})	
P2 = P.trendy1
P2[P.trendy1 <= 5] = 1
P2[P.trendy1 > 5] = 0.05
pts.sp.sig2 = Ex.pts(P2, sig.level = 0.1) #extract significant relation points


#############################################################################################
### 	read into MOD17 GPP/ CO2 NEE

R.grd = raster()   
R.grd[] <- 1						   
R.grd2 = crop(R.grd, usa.state)
R.grd = rasterize(usa.state, R.grd2, field = 1)
 
gpp.annual = stack("F:/zhihua/dataset/mod17a2/processed/c.gpp.annual.grd")

#aggreate into 1 degree
gpp.annual2 = list()
for (i in 1:nlayers(gpp.annual)){
gpp.annual2[[i]] <- aggregate(gpp.annual[[i]], fact=20, fun=mean, expand=TRUE, na.rm=TRUE)
gpp.annual2[[i]] <- focal(gpp.annual2[[i]], w=matrix(1/9,nrow=3,ncol=3), na.rm = TRUE)

}

gpp.annual2 = stack(gpp.annual2)
gpp.annual2 = gpp.annual2*R.grd

# for annual temperature and precipitations
prep = list()

for (yr in 2000:2014){

prep1 = stack(paste0("F:/zhihua/dataset/cru_ts3.23/prep", yr, ".grd"))
prep1 = crop(prep1, R.grd)
prep1 = calc(prep1, sum, na.rm = TRUE)
prep1 = aggregate(prep1, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)

prep[[yr - 1999]] <- prep1

}

prep = stack(prep)
prep.mn = calc(prep, mean, na.rm = TRUE)

#plot correlationship along the GPP gradient
Q1 = function(x){quantile(x, probs = c(0.1,0.9))}

se <- function(x) sqrt(var(x)/length(x))

gpp.mn = crop(calc(gpp.annual2, mean, na.rm = TRUE), usa.state)
gpp.mn = gpp.mn*R.grd
#remove some very low values
gpp.mn[gpp.mn < 200] = NA
q10 = quantile(gpp.mn, probs = seq(0, 1, length.out= 15), na.rm = TRUE)
gpp.mn2 = cut(gpp.mn, breaks = q10)

prep.mn[prep.mn < 10] = NA
q10 = quantile(prep.mn, probs = seq(0, 1, length.out= 25), na.rm = TRUE)
prep.mn2 = cut(prep.mn, breaks = q10)



dat1.corr2.df = data.frame(r = zonal(R, prep.mn2, mean)[,2],
						   sd = zonal(R, prep.mn2, sd)[,2],
						   r.trendy = zonal(R.trendy1, prep.mn2, mean)[,2],
						   sd.trendy = zonal(R.trendy1, prep.mn2, sd)[,2],						   
						   gpp = zonal(gpp.mn, prep.mn2, mean)[,2],
						   gpp.sd = zonal(gpp.mn, prep.mn2, sd)[,2],
						   prep = zonal(prep.mn, prep.mn2, mean)[,2],
						   prep.sd = zonal(prep.mn, prep.mn2, sd)[,2])

						   

###calculate correlation between region > 750 and region < 750
require(ncdf4)

nc <- nc_open("F:/zhihua/dataset/ecoregion/regions.nc")
print(nc)
v2 <- nc$var[[1]]
data2 <- ncvar_get( nc, v2 ) #data2 is an array

grid.area = raster(t(data2),xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

grid.area.nee = grid.area*R.grd		

## read into nee
nee.annual = stack("F:/zhihua/dataset/co2_inversion/NEE.annual.mean.grd")
nee.annual = nee.annual*R.grd

# area unit: units: m2
usa.grd.nee = prep.mn
usa.grd.nee[prep.mn > 750] = 2
usa.grd.nee[prep.mn < 750] = 1

usa.grd.nee[prep.mn > 900] = 2
usa.grd.nee[prep.mn < 900] = 1

usa.grd.nee = usa.grd.nee*R.grd

# usa.state = readOGR(dsn="C:\\zhihua\\dataset\\ecoregion", layer = "usa.state")
# usa.state = usa.state[-which(usa.state$NAME_1 == "Hawaii" | usa.state$NAME_1 == "Alaska"),]
# usa.grd.nee = rasterize(usa.state, R.grd, field = "region")
# usa.grd.nee = usa.grd.nee*R.grd

nee.annual2 = nee.annual*grid.area.nee
gpp.annual3 = gpp.annual2*grid.area.nee
 
nee.mn = calc(nee.annual, mean, na.rm = TRUE) 
gpp.mn = calc(gpp.annual2, mean, na.rm = TRUE) 
 
nee.df2 = zonal(nee.annual2, usa.grd.nee, fun='sum', digits=0, na.rm=TRUE)[,c(2:16)]/1000000000000000 
gpp.df2 = zonal(gpp.annual3, usa.grd.nee, fun='sum', digits=0, na.rm=TRUE)[,c(2:16)]/1000000000000000 
er.df2 = gpp.df2 - nee.df2
 
nee.df3 = sweep(nee.df2, MARGIN = 1, STATS = apply(nee.df2,1,mean), FUN = "-")
gpp.df3 = sweep(gpp.df2, MARGIN = 1, STATS = apply(gpp.df2,1,mean), FUN = "-")
er.df3 = sweep(er.df2, MARGIN = 1, STATS = apply(er.df2,1,mean), FUN = "-")

########read into trendy
trendy.nee = stack("F:/zhihua/dataset/trendy/Trendy.NBP.mean.grd")
trendy.gpp = stack("F:/zhihua/dataset/trendy/Trendy.GPP.mean.grd")
trendy.ra = stack("F:/zhihua/dataset/trendy/Trendy.Ra.mean.grd")
trendy.rh = stack("F:/zhihua/dataset/trendy/Trendy.Rh.mean.grd")

trendy.nee = trendy.nee*R.grd
trendy.gpp = trendy.gpp*R.grd
trendy.ra = trendy.ra*R.grd
trendy.rh = trendy.rh*R.grd
trendy.er = trendy.ra + trendy.rh

nee.mn.trendy = calc(trendy.nee, mean, na.rm=TRUE)
gpp.mn.trendy = calc(trendy.gpp, mean, na.rm=TRUE)

trendy.nee = trendy.nee*grid.area.nee
trendy.gpp = trendy.gpp*grid.area.nee
trendy.er = trendy.er*grid.area.nee

nee.df2.t = zonal(trendy.nee, usa.grd.nee, fun='sum', digits=0, na.rm=TRUE)[,c(2:12)]/1000000000000000 
gpp.df2.t = zonal(trendy.gpp, usa.grd.nee, fun='sum', digits=0, na.rm=TRUE)[,c(2:12)]/1000000000000000 
er.df2.t = zonal(trendy.er, usa.grd.nee, fun='sum', digits=0, na.rm=TRUE)[,c(2:12)]/1000000000000000 
 
nee.df3.t = sweep(nee.df2.t, MARGIN = 1, STATS = apply(nee.df2.t,1,mean), FUN = "-")
gpp.df3.t = sweep(gpp.df2.t, MARGIN = 1, STATS = apply(gpp.df2.t,1,mean), FUN = "-")
er.df3.t = sweep(er.df2.t, MARGIN = 1, STATS = apply(er.df2.t,1,mean), FUN = "-")


Year = 2000:2010
plot(Year, gpp.df3.t[1,], ylim = c(-0.5, 0.5), type = "b", pch = 1, col = 1, cex = 1.5, lwd = 3, 
     cex.axis = 1.5, xlab = "Year", ylab = "Anomaly (Pg C)",cex.lab = 1.5)

abline(h = 0, lty = 2)	 
legend("topleft",legend=c("GPP", "NEE","ER"),lwd = 3,
				  col = 1:3) 
	 	 
par(new=T)
plot(Year,nee.df3.t[1,], ylim = c(-0.5, 0.5), type = "b", pch = 1, col = 2, bty='n',cex = 1.5, lwd = 3, 
     cex.axis = 1.5, xlab = "Year", ylab = "Anomaly (Pg C)",cex.lab = 1.5)

par(new=T)
plot(Year,er.df3.t[1,], ylim = c(-0.5, 0.5), type = "b", pch = 1, col = 3, bty='n',cex = 1.5, lwd = 3, 
     cex.axis = 1.5, xlab = "Year", ylab = "Anomaly (Pg C)",cex.lab = 1.5)
	 
	 
r1 = cor.test(gpp.df3.t[1,], nee.df3.t[1,])
r2 = cor.test(er.df3.t[1,], nee.df3.t[1,])

text(x = 2006, y = 0.43, paste("GPP:NEE: ", 
                                "r = ", round(r1$estimate, 2),"; ",
								"p = ", round(r1$p.value, 2)), cex = 1.5)  
text(x = 2006, y = 0.35, paste("ER:NEE:", 
                                "r = ", round(r2$estimate, 2),"; ",
								"p = ", round(r2$p.value, 2)), cex = 1.5)   


Year = 2000:2010
plot(Year, gpp.df3.t[2,], ylim = c(-0.5, 0.5), type = "b", pch = 1, col = 1, cex = 1.5, lwd = 3, 
     cex.axis = 1.5, xlab = "Year", ylab = "Anomaly (Pg C)",cex.lab = 1.5)

abline(h = 0, lty = 2)	 
	 
par(new=T)
plot(Year,nee.df3.t[2,], ylim = c(-0.5, 0.5), type = "b", pch = 1, col = 2, bty='n',cex = 1.5, lwd = 3, 
     cex.axis = 1.5, xlab = "Year", ylab = "Anomaly (Pg C)",cex.lab = 1.5)

par(new=T)
plot(Year,er.df3.t[2,], ylim = c(-0.5, 0.5), type = "b", pch = 1, col = 3, bty='n',cex = 1.5, lwd = 3, 
     cex.axis = 1.5, xlab = "Year", ylab = "Anomaly (Pg C)",cex.lab = 1.5)
	 

r1 = cor.test(gpp.df3.t[2,], nee.df3.t[2,])
r2 = cor.test(er.df3.t[2,], nee.df3.t[2,])
	 
text(x = 2006, y = 0.43, paste("GPP:NEE: ", 
                                "r = ", round(r1$estimate, 2),"; ",
								"p = ", round(r1$p.value, 2)), cex = 1.5)  
text(x = 2006, y = 0.35, paste("ER:NEE:", 
                                "r = ", round(r2$estimate, 2),"; ",
								"p = ", round(r2$p.value, 2)), cex = 1.5)  


# random sampling to calculate r and p
boot.r = function(x,y,ln = 5, n = 100){  #ln: minimum to calculate r, n: number of sampling
# sampling at least 5  of the data
D = c()
NofSim = 0	
repeat{
NofSim = NofSim+1
n1 = sample(ln:length(x),1)
idx = sample(1:length(x),n1)
r1 = cor.test(x[idx], y[idx])
D = rbind(D, c(r1$estimate,r1$p.value))
if(NofSim > n){break}
}
D = data.frame(D)
colnames(D) <- c("r","p")
return(D)

}								
								
r1 = boot.r(x = gpp.df3[1,], y = nee.df3[1,],ln = 7, n = 100)
r2 = boot.r(x = er.df3[1,], y = nee.df3[1,],ln = 7, n = 100)
r3 = boot.r(x = gpp.df3[2,], y = nee.df3[2,],ln = 7, n = 100)
r4 = boot.r(x = er.df3[2,], y = nee.df3[2,],ln = 7, n = 100)

r5 = boot.r(x = gpp.df3.t[1,], y = nee.df3.t[1,],ln = 7, n = 100)
r6 = boot.r(x = er.df3.t[1,], y = nee.df3.t[1,],ln = 7, n = 100)
r7 = boot.r(x = gpp.df3.t[2,], y = nee.df3.t[2,],ln = 7, n = 100)
r8 = boot.r(x = er.df3.t[2,], y = nee.df3.t[2,],ln = 7, n = 100)

se = function(x) sqrt(var(x)/length(x))
							
library(ggplot2)
	

D2 = data.frame(mean.r = c(mean(r1$r), mean(r2$r),mean(r3$r),mean(r4$r),mean(r5$r),mean(r6$r),mean(r7$r),mean(r8$r)),
     sd.r = c(sd(r1$r), sd(r2$r),sd(r3$r),sd(r4$r),sd(r5$r),sd(r6$r),sd(r7$r),sd(r8$r)), 								
     mean.p = c(mean(r1$p), mean(r2$p),mean(r3$p),mean(r4$p),mean(r5$p),mean(r6$p),mean(r7$p),mean(r8$p)),
     sd.p = c(sd(r1$p), sd(r2$p),sd(r3$p),sd(r4$p),sd(r5$p),sd(r6$p),sd(r7$p),sd(r8$p)), 
	 method = c(rep("RS",4), rep("TRENDY",4)),
	 region = c("< 750","< 750","> 750","> 750","< 750","< 750","> 750","> 750"),
	 R = c("GPP:NEE","ER:NEE","GPP:NEE","ER:NEE","GPP:NEE","ER:NEE","GPP:NEE","ER:NEE"),
	 R2 = c("GPP:NEE1","ER:NEE2","GPP:NEE1","ER:NEE2","GPP:NEE3","ER:NEE4","GPP:NEE3","ER:NEE4")
	 )


D2$sym = c("**", "","","","***","","*","")
# write.csv(D1, "F:/zhihua/dataset/results2/corr.csv")
# D1 = read.csv("F:/zhihua/dataset/results2/corr.csv")
levels(D2$method) <- c("Constrained Global Obs", "TRENDY")

library(ggplot2)
color1 = c(rgb(0,0,1,0.95),rgb(178/255,178/255,0,0.95),rgb(0,0,1,0.25),rgb(178/255,178/255,0,0.25))
ylab=expression("" ~ g ~ C ~ m^{-2} ~ yr ^{-1}~ " Per 100 mm")

p1 <- ggplot(data=D2, aes(x=region, y=mean.r, fill=R2)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 
    facet_grid(.~method)+			 
	# ylab(ylab) + 	 
    # ylab(expression(paste(beta ["temporal"]))) + 
	# xlab("Region by MAP") + 
	ylab("Pearson's r") + 
	xlab("MAP (mm)") + 
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	geom_errorbar(aes(ymin=mean.r-sd.r, ymax=mean.r+sd.r),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
	# theme(legend.position=c(.35, .8)) + 	
	# theme(legend.title=element_blank()) +
	# theme(legend.text = element_text(size = 4)) +
	theme(legend.position="none") +
	theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
    theme(strip.text.x = element_text(size=18))+
    theme(strip.text.y = element_text(size=18)) +
    scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(D2$R2),
                    labels=levels(D2$R2)) + 
    geom_text(aes(label=sym), vjust=0.5,hjust=0,size = 18)		
		
ggsave("F:/zhihua/dataset/results2/fig1s.png", width = 9, height = 6, units = "in")
	
### fig 1, add GPP/NEE

png("F:/zhihua/dataset/results2/fig1-4.png",height = 2700, width = 3300, res = 300, units = "px")

par(mfrow=c(3,2),mar=c(0,0,0,0)+.0)

# plot GPP

####### Fpar
plot(gpp.mn, 
                    # zlim=c(-1,1),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
plot(usa.state, add = TRUE, lty = 2)
# tx = expression("Mean Annual FPAR by MODIS, 2007 - 2014")
tx = "a)"
# text(x = 0, y = -85, "Mean Annual GPP by MODIS (g C/m-2*yr, 2000-2014)", cex = 2)

#add legend
plot(gpp.mn, 
         col = rainbow(n = 100), 
         legend.only=TRUE,
		 legend.width=1, legend.shrink=0.75,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
        smallplot=c(0.8,0.85, 0.2,0.45))
	 par(mar = par("mar"))
text(x = -123, y = 50, "a)", cex = 2)

# add unit
tx = expression("" ~ gC ~ m^{-2} ~ yr ^{-1}~ "")
text(x = -70, y = 38, tx, cex = 1.5)
# text(x = -69, y = 38, "GPP\n gC m-2*yr-1", cex = 1.5)

# plot nee

####### Fpar
plot(nee.mn, 
                    # zlim=c(-1,1),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
plot(usa.state, add = TRUE, lty = 2)
# tx = expression("Mean Annual FPAR by MODIS, 2007 - 2014")
tx = "a)"
# text(x = 0, y = -85, "Mean Annual GPP by MODIS (g C/m-2*yr, 2000-2014)", cex = 2)

#add legend
plot(nee.mn, 
         col = rainbow(n = 100), 
         legend.only=TRUE,
		 legend.width=1, legend.shrink=0.75,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
        smallplot=c(0.8,0.85, 0.2,0.45))
	 par(mar = par("mar"))
text(x = -123, y = 50, "b)", cex = 2)

# add unit
tx = expression("" ~ gC ~ m^{-2} ~ yr ^{-1}~ "")
text(x = -70, y = 38, tx, cex = 1.5)
# text(x = -69, y = 38, "NEE\n gC m-2*yr-1", cex = 1.5)

# plot mean r from rs model
my.colors = colorRampPalette(c("blue", "white", "red"))

plot(R, zlim=c(-1,1),col = my.colors(100), 
					# main = "Obeserved Mean r Between GPP and NEE",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE,
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)

P = calc(stack(P.stack1, P.stack2), function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 4] = 1
P1[P > 4] = 0.05
pts.sp.sig2 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points

plot(pts.sp.sig2, add = TRUE, cex = 0.1)

# plot.new()
plot(R, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=1, legend.shrink=0.75,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
        smallplot=c(0.8,0.85, 0.2,0.45))
	 par(mar = par("mar"))
text(x = -123, y = 50, "c)", cex = 2)

# text(x = -70, y = 38, "r", cex = 1.6)

# plot mean r from trendy model

plot(R.trendy1, zlim=c(-1,1),col = my.colors(100), 
					# main = "Obeserved Mean r Between GPP and NEE",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE,
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)

P.trendy1 = calc(P.trendy, function(x){length(which(x <= 0.1))})	
P2 = P.trendy1
P2[P.trendy1 <= 5] = 1
P2[P.trendy1 > 5] = 0.05
pts.sp.sig2 = Ex.pts(P2, sig.level = 0.1) #extract significant relation points

plot(pts.sp.sig2, add = TRUE, cex = 0.1)

text(x = -123, y = 50, "d)", cex = 2)

# box()
#######plot r response to prep
par(mar=c(3,3,0,3)+2)

plot(r~prep, data = dat1.corr2.df, type = "l", 
								  ylim = c(-0.2,1), 
								  cex = 1.5, lwd = 4, 
								  col = "blue",
								  xlab = "Mean Anuual Precipitation (MAP: mm)", 
								  ylab = "Pearson's r", cex.axis = 1.5, cex.lab = 1.3)
			
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$r-dat1.corr2.df$sd), dat1.corr2.df$r+dat1.corr2.df$sd), 
        col=rgb(0, 0, 1,0.25),
		border = NA)

# overlay trendy r	
par(new=TRUE)
plot(dat1.corr2.df$prep, dat1.corr2.df$r.trendy, type="l", 
								  ylim = c(-0.2,1), 
								  col=rgb(178/255, 178/255, 0,1),lwd = 4,
								  bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)
								  
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$r.trendy-dat1.corr2.df$sd.trendy), dat1.corr2.df$r.trendy+dat1.corr2.df$sd.trendy), 
        col=rgb(178/255, 178/255, 0,0.25),
		border = NA)
								  
								  
# overlay EC r	
# points(cor.df$gpp,cor.df$r,ylim = c(-0.2,1), col="red",pch = 1, cex = 3)						  
# text(cor.df$gpp,cor.df$r, cor.df$SITE_ID)	
		
# overlay GPP
par(new=TRUE)
plot(dat1.corr2.df$prep, dat1.corr2.df$gpp, type="l", 
								  ylim = c(0, 1800), 
								  col=rgb(128/255, 255/255, 204/255,1),
								  lwd = 4, bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$gpp-dat1.corr2.df$gpp.sd), dat1.corr2.df$gpp+dat1.corr2.df$gpp.sd), 
        col=rgb(128/255, 255/255, 204/255,0.25),
		border = NA)
										  

legend("topleft", 
	   # inset=0.05, 
	   # legend = c("ENF","DBF","MF","SHB","GRA","CRO"),
	   legend = c("Constrained Global Obs","TRENDY", "GPP"),
	   horiz=F,
	   lwd = 4,
	   col = c(rgb(0, 0, 1,1), rgb(178/255, 178/255, 0,1),rgb(128/255, 255/255, 204/255,1)),
	   text.col = c(rgb(0, 0, 1,1), rgb(178/255, 178/255, 0,1),rgb(128/255, 255/255, 204/255,1)),
	   cex = 1,
	   box.col = "transparent",
	   bg = "transparent")

# text(x = 0, y = 1600, "c)", cex = 2)
# box()	
								  
axis(side = 4,cex.axis = 1.3, col.axis = rgb(128/255, 255/255, 204/255,1))
mtext(side = 4, line = 3, expression("GPP" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")),cex = 1.3, col = rgb(128/255, 255/255, 204/255,1))


dev.off()

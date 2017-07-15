
## plot figure 1
Work.Dir <- "C:/zhihua/dataset/fpar"
setwd(Work.Dir)


#load library
library(rgdal)
library(raster)
library(rasterVis)

na.state = readOGR(dsn="C:\\zhihua\\dataset\\ecoregion", layer = "na.state2")
na.state = na.state[-which(na.state$NAME_1 == "Prince Edward Island"|
                           na.state$NAME_1 == "Hawaii"),]

na.ext <- extent(-180,-48,15,85)

usa.state = na.state[which(na.state$ISO == "USA" & na.state$NAME_1 != "Alaska"), ]

##################################################################		
######plot ensemble mean for 8 maps :::: RS
R.stack1 = stack("C:/zhihua/dataset/ct/cor.stack.ct.grd")
P.stack1 = stack("C:/zhihua/dataset/ct/P.stack.ct.grd")

# R.stack2 = stack("C:/zhihua/dataset/cte/cor.stack.cte.grd")
# P.stack2 = stack("C:/zhihua/dataset/cte/P.stack.cte.grd")

R.stack2 = stack("C:/zhihua/dataset/ct/cor.stack.ecmod.grd")
P.stack2 = stack("C:/zhihua/dataset/ct/p.stack.ecmod.grd")
	
R = calc(stack(R.stack1, R.stack2), mean, na.rm = TRUE)	
P = calc(stack(P.stack1, P.stack2), function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 4] = 1
P1[P > 4] = 0.05
pts.sp.sig2 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points

#plot individually

R.stack = stack(R.stack1, R.stack2)
P.stack = stack(P.stack1, P.stack2)

#Col.name = c("GPP:NEE","NDVI:NEE","FPAR:NEE","SIF:NEE","GPP:NEE","NDVI:NEE","FPAR:NEE","SIF:NEE")
Col.name = c("GPP:CT-NEE","NDVI:CT-NEE","FPAR:CT-NEE","SIF:CT-NEE","GPP:EC-NEE","NDVI:EC-NEE","FPAR:EC-NEE","SIF:EC-NEE")
let.name = c("a)","b)","c)","d)","e)","f)","g)","h)")

png("F:/zhihua/dataset/results2/productivity.nee.cor.rs-2.png",height = 2500, width = 2500, res = 300, units = "px")

par(mfrow=c(4,2),mar=c(0,0,0,0)+.1)
my.colors = colorRampPalette(c("blue", "white", "red"))

for (i in 1:8){
pts.sp.sig1 = Ex.pts(P.stack[[i]], sig.level = 0.1) #extract significant relation points

plot(R.stack[[i]], zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)

# text(x = -114.5, y = 27.5, Col.name[i], cex = 2)
text(x = -116, y = 28, let.name[i], cex = 3)
plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)
}


dev.off()

R = calc(R.stack, mean, na.rm = TRUE)		
#R <- focal(R, w=matrix(1/9,nrow=3,ncol=3), na.rm = TRUE)
		
png("C:/zhihua/dataset/results/productivity.nee.cor.annual.usa.mean.png",height = 1500, width = 2500, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
levelplot(R, zlim = c(-1,1),main = "Obeserved Mean r Between GPP and NEE",
maxpixels = nrow(R)*ncol(R),
margin = FALSE,
# FUN.margin=median,
col.regions= colorRampPalette(c("blue", "white", "red"))(255),
colorkey= list(labels= list(cex = 1.5)),
scales=list(x=list(cex=1),y=list(cex=1)),
xlab=list(label = "Longtitude", cex=1),ylab=list(label = "Latitude", cex=1),
par.strip.text=p.strip) +
latticeExtra::layer(sp.polygons(na.state, col = "black", lwd = 1.5)) + 
latticeExtra::layer(sp.points(pts.sp.sig2,pch=20, cex=0.25, col="black"))

dev.off()		

#read into flux tower and overlay on the mean r map
flux.info.sp3 <- readOGR("C:/zhihua/dataset/flux2015jul", "flux.info.sp3")

png("C:/zhihua/dataset/results/productivity.nee.cor.annual.usa.mean-3.png",height = 1500, width = 2500, res = 300, units = "px")
	
my.colors = colorRampPalette(c("#053061", "#f7f7f7", "#67001f"))
my.colors = colorRampPalette(c("blue", "white", "red"))

plot(R, zlim=c(-1,1),col = my.colors(100), 
					main = "Obeserved Mean r Between GPP and NEE",
					legend=FALSE,
                    # axes=FALSE,
                    # box=FALSE
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig2, add = TRUE, cex = 0.1)
					
plot(R, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=0.5, legend.shrink=0.5,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
		 smallplot=c(0.73,0.76, 0.25,0.45))
	 par(mar = par("mar"))

	 
plot(flux.info.sp3[which(flux.info.sp3$IGBP == "DBF" | flux.info.sp3$IGBP == "MF"), ], pch = 21, cex = 2, add = T)
plot(flux.info.sp3[which(flux.info.sp3$IGBP == "ENF"), ], pch = 22, cex = 2, add = T)
plot(flux.info.sp3[which(flux.info.sp3$IGBP == "GRA" | flux.info.sp3$IGBP == "CSH" | flux.info.sp3$IGBP == "OSH" 
                    | flux.info.sp3$IGBP == "SAV" | flux.info.sp3$IGBP == "WSA" | flux.info.sp3$IGBP == "CRO"), ], pch = 23, cex = 2, add = T)

legend("bottomleft", 
	   # inset=0.05, 
	   legend = c("DBF/MF (7)","ENF (3)","GRA/Shurb/CRO(6)" ),
	   horiz=F,
	   pch = c(21:23),
	   #col = 1:6,
	   #text.col = 1:6,
	   cex = 1,
	   box.col = "transparent",
	   bg = "transparent")	 

dev.off()

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

#plot individually
Mod.name = c("CLM4C", "CLM4CN", "HYLAND", "LPJ", "LPJ-GUESS", "OCN", "ORCHIDEE","TRIFFID","VEGAS","SDGVM")

png("C:/zhihua/dataset/results2/trendy.cor.png",height = 2500, width = 2500, res = 300, units = "px")

par(mfrow=c(4,3),mar=c(0,0,0,0)+.1)
my.colors = colorRampPalette(c("blue", "white", "red"))

for (i in 1:10){
pts.sp.sig1 = Ex.pts(P.trendy[[i]], sig.level = 0.1) #extract significant relation points

plot(R.trendy[[i]], zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, Mod.name[i], cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)
}


plot(R.trendy1, zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, "Mean r", cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)

P.trendy1 = calc(P.trendy, function(x){length(which(x <= 0.1))})	
P1 = P.trendy1
P1[P <= 5] = 1
P1[P > 5] = 0.05
pts.sp.sig2 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points

plot(pts.sp.sig1, add = TRUE, cex = 0.1)

plot.new()
plot(R.trendy1, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=1, legend.shrink=0.75,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="r", side=4, font=2, line=2.5, cex=1.2),
    smallplot=c(0.3,0.5, 0.2,0.7))
	 par(mar = par("mar"))

legend(x = -100, y = 43, legend = "Correlation between NBP and GPP (2000-2010)",  cex = 1.5,  box.lwd = 0,box.col = "white",bg = "white")

dev.off()

# use ggplot2
# https://nrelscience.org/2013/05/30/this-is-how-i-did-it-mapping-in-r-with-ggplot2/
#load libraries
library(raster)
library(ggplot2)

#convert the raster to points for plotting
R.trendy1.p <- rasterToPoints(R.trendy1)

#Make the points a dataframe for ggplot
R.trendy1.df <- data.frame(R.trendy1.p)
#Make appropriate column headings
colnames(R.trendy1.df) <- c("Longitude", "Latitude", "r")

#Get the significant points
pts.sp.sig2.df <- data.frame(coordinates(pts.sp.sig2))
colnames(pts.sp.sig2.df) <- c("x","y")

#Now make the map
ggplot(data=R.trendy1.df, aes(y=Latitude, x=Longitude)) +
geom_raster(aes(fill=r)) +
geom_point(data=pts.sp.sig2.df, aes(x=x, y=y), color="red", size=1, shape=2) +
theme_bw() +
coord_equal() +
scale_fill_gradient2("r", 
                     breaks = seq(-1,1, by = 0.2),
					 limits=c(-1,1),
					 low = "red", 
					 mid = "yellow",
					 high = "blue", 
					 midpoint = 0) +
theme(axis.title.x = element_text(size=16),
axis.title.y = element_text(size=16, angle=90),
axis.text.x = element_text(size=14),
axis.text.y = element_text(size=14),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
legend.position = "right",
legend.key = element_blank()
)


##################################################################		
#
#############################################################################################
### 	

R.grd = raster()   
R.grd[] <- 1						   
R.grd2 = crop(R.grd, usa.state)
R.grd = rasterize(usa.state, R.grd2, field = 1)
 
gpp.annual = stack("C:/zhihua/dataset/mod17a2/processed/c.gpp.annual.grd")

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

						   
##### calculate EC correlationships
# use code until 241: NA_Carbon/fluxnet2.r 
dat1 = data.frame(dat.df.ano.rel2) 
SITE_ID = unique(dat1$site_id)

cor.df = c()
for (i in SITE_ID){
dat11 = subset(dat1, site_id == i)
r = cor.test(x=dat11$gpp, y=dat11$nee, method = 'spearman')
r2 = cor.test(x=dat11$er, y=dat11$nee, method = 'spearman')

cor.df = rbind(cor.df,c(r$estimate, r$p.value,r2$estimate, r2$p.value))

}

cor.df = data.frame(SITE_ID = SITE_ID, cor.df)
colnames(cor.df) <- c("SITE_ID","r","p","r2","p2")
cor.df = merge(cor.df, flux.info.sp3@data[,c("SITE_ID","gpp","prep")], by.x = "SITE_ID",by.y = "SITE_ID")

cor.df$r3 <- abs(cor.df$r)-abs(cor.df$r2)

						   
plot(usa.state)
plot(flux.info.sp3, add = T)
text(x = flux.info3$LOCATION_LONG, y = flux.info3$LOCATION_LAT, flux.info3$SITE_ID)	

	
png("F:/zhihua/dataset/results2/productivity.nee.cor.annual.usa2.mean.png",height = 1500, width = 3000, res = 300, units = "px")

# http://htmlcolorcodes.com/
par(mar=c(3,3,0,3)+2)

plot(r~prep, data = dat1.corr2.df, type = "l", 
								  ylim = c(-0.2,1), 
								  cex = 2, lwd = 4, 
								  col = "green",
								  xlab = "Anuual Mean Precipitation (MAP: mm)", 
								  ylab = "Pearson's r", cex.axis = 1.5, cex.lab = 1.6)
			
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$r-dat1.corr2.df$sd), dat1.corr2.df$r+dat1.corr2.df$sd), 
        col=rgb(0, 0.5, 0,0.25),
		border = NA)

# overlay trendy r	
par(new=TRUE)
plot(dat1.corr2.df$prep, dat1.corr2.df$r.trendy, type="l", 
								  ylim = c(-0.2,1), 
								  col="blue",lwd = 4,
								  bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)
								  
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$r.trendy-dat1.corr2.df$sd.trendy), dat1.corr2.df$r.trendy+dat1.corr2.df$sd.trendy), 
        col=rgb(0, 0, 0.5,0.25),
		border = NA)
								  
								  
# overlay EC r	
# points(cor.df$gpp,cor.df$r,ylim = c(-0.2,1), col="red",pch = 1, cex = 3)						  
# text(cor.df$gpp,cor.df$r, cor.df$SITE_ID)	
		
# overlay GPP
par(new=TRUE)
plot(dat1.corr2.df$prep, dat1.corr2.df$gpp, type="l", 
								  ylim = c(0, 1800), 
								  col=rgb(0.5, 0.5, 0.5,1),
								  lwd = 4, bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$gpp-dat1.corr2.df$gpp.sd), dat1.corr2.df$gpp+dat1.corr2.df$gpp.sd), 
        col=rgb(0.5, 0.5, 0.5,0.25),
		border = NA)
										  

legend("topleft", 
	   # inset=0.05, 
	   # legend = c("ENF","DBF","MF","SHB","GRA","CRO"),
	   legend = c("RS","TRENDY"),
	   horiz=F,
	   lwd = 4,
	   col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1)),
	   text.col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1)),
	   cex = 1.5,
	   box.col = "transparent",
	   bg = "transparent")

										  
axis(side = 4)
#mtext(side = 4, line = 3, 'precipitation (mm/yr)',cex.axis = 1.5, cex.lab = 1.6, col = "blue")
mtext(side = 4, line = 3, expression("GPP" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")),cex = 1.6, col = rgb(0.5, 0.5, 0.5,1))

dev.off()
		
###calculate correlation between region > 750 and region < 750
require(ncdf4)

nc <- nc_open("C:/zhihua/dataset/ecoregion/regions.nc")
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
							
D1 = data.frame(mean.r = c(mean(r1$r), mean(r2$r),mean(r3$r),mean(r4$r),mean(r5$r),mean(r6$r),mean(r7$r),mean(r8$r)),
     sd.r = c(sd(r1$r), sd(r2$r),sd(r3$r),sd(r4$r),sd(r5$r),sd(r6$r),sd(r7$r),sd(r8$r)), 								
     mean.p = c(mean(r1$p), mean(r2$p),mean(r3$p),mean(r4$p),mean(r5$p),mean(r6$p),mean(r7$p),mean(r8$p)),
     sd.p = c(sd(r1$p), sd(r2$p),sd(r3$p),sd(r4$p),sd(r5$p),sd(r6$p),sd(r7$p),sd(r8$p)), 
	 method = c(rep("RS",4), rep("TRENDY",4)),
	 region = c("< 750","< 750","> 750","> 750","< 750","< 750","> 750","> 750"),
	 R = c("GPP:NEE","ER:NEE","GPP:NEE","ER:NEE","GPP:NEE","ER:NEE","GPP:NEE","ER:NEE"))

D1$sym = c("**", "","","","***","","*","")
# write.csv(D1, "F:/zhihua/dataset/results2/corr.csv")
# D1 = read.csv("F:/zhihua/dataset/results2/corr.csv")
levels(D1$method) <- c("Constrained Global Obs", "TRENDY")


color1 = c("grey20", "grey90")
ylab=expression("" ~ g ~ C ~ m^{-2} ~ yr ^{-1}~ " Per 100 mm")

p1 <- ggplot(data=D1, aes(x=region, y=mean.r, fill=R)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 
    facet_grid(.~method)+			 
	# ylab(ylab) + 	 
    # ylab(expression(paste(beta ["temporal"]))) + 
	# xlab("Region by MAP") + 
	ylab("Pearson's r") + 
	xlab("Region by MAP (mm)") + 
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	geom_errorbar(aes(ymin=mean.r-sd.r, ymax=mean.r+sd.r),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) + 
	theme(legend.position=c(.35, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 6)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=12),axis.text.x  = element_text(colour="black",size=12))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=12),axis.text.y  = element_text(colour="black",size=12))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) +
    scale_fill_manual(values=color1, 
                    name="",
                    breaks=levels(D1$R),
                    labels=levels(D1$R)) + 
    geom_text(aes(label=sym), vjust=0.5,hjust=0,size = 8)		
	

##plot figure1


########

png("F:/zhihua/dataset/results2/fig1.png",height = 2000, width = 3300, res = 300, units = "px")

par(mfrow=c(2,2),mar=c(0,0,0,0)+.0)


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
text(x = -123, y = 51, "a)", cex = 2)

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

text(x = -123, y = 51, "b)", cex = 2)


#######plot r response to prep
par(mar=c(3,3,0,3)+2)

plot(r~prep, data = dat1.corr2.df, type = "l", 
								  ylim = c(-0.2,1), 
								  cex = 1.5, lwd = 4, 
								  col = "green",
								  xlab = "Mean Anuual Precipitation (MAP: mm)", 
								  ylab = "Pearson's r", cex.axis = 1.5, cex.lab = 1.3)
			
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$r-dat1.corr2.df$sd), dat1.corr2.df$r+dat1.corr2.df$sd), 
        col=rgb(0, 0.5, 0,0.25),
		border = NA)

# overlay trendy r	
par(new=TRUE)
plot(dat1.corr2.df$prep, dat1.corr2.df$r.trendy, type="l", 
								  ylim = c(-0.2,1), 
								  col="blue",lwd = 4,
								  bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)
								  
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$r.trendy-dat1.corr2.df$sd.trendy), dat1.corr2.df$r.trendy+dat1.corr2.df$sd.trendy), 
        col=rgb(0, 0, 0.5,0.25),
		border = NA)
								  
								  
# overlay EC r	
# points(cor.df$gpp,cor.df$r,ylim = c(-0.2,1), col="red",pch = 1, cex = 3)						  
# text(cor.df$gpp,cor.df$r, cor.df$SITE_ID)	
		
# overlay GPP
par(new=TRUE)
plot(dat1.corr2.df$prep, dat1.corr2.df$gpp, type="l", 
								  ylim = c(0, 1800), 
								  col=rgb(0.5, 0.5, 0.5,1),
								  lwd = 4, bty='n',pch='',ylab='',xlab='',yaxt='n',xaxt='n', ann=FALSE)
polygon(c(rev(dat1.corr2.df$prep), dat1.corr2.df$prep), 
		c(rev(dat1.corr2.df$gpp-dat1.corr2.df$gpp.sd), dat1.corr2.df$gpp+dat1.corr2.df$gpp.sd), 
        col=rgb(0.5, 0.5, 0.5,0.25),
		border = NA)
										  

legend("topleft", 
	   # inset=0.05, 
	   # legend = c("ENF","DBF","MF","SHB","GRA","CRO"),
	   legend = c("Constrained Global Obs","TRENDY"),
	   horiz=F,
	   lwd = 4,
	   col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1)),
	   text.col = c(rgb(0, 1, 0,1), rgb(0, 0, 1,1)),
	   cex = 1.1,
	   box.col = "transparent",
	   bg = "transparent")

# text(x = 100, y = 100, "c)", cex = 2)
									  
axis(side = 4,cex = 1.2)
mtext(side = 4, line = 3, expression("GPP" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")),cex = 1.3, col = rgb(0.5, 0.5, 0.5,1))

## plot 4
require(grid)
par(new=TRUE)
print(p1, vp=viewport(.73, .25, .45, .45)) #viewport, first two is the x,y coor; last 2 is the inset size

text(x = 0.5, y = 0.5, "c)", cex = 2)

dev.off()


###plot 17 EC flux tower on the land cover map
# classific land cover into forest and non-forest[shrubland, grassland, and crop]
lc = raster("C:/zhihua/dataset/ecoregion/LCTypeNA.tif")
#only for US
lc_us = crop(lc, usa.state)
# aggregate to 0.05 degree
lc_us2 = aggregate(lc_us, fact=12, fun=modal, expand=TRUE, na.rm=TRUE)
m1 = rasterize(usa.state, lc_us2, field = 1)
lc_us2 = lc_us2*m1
oldclas <- unique(lc_us2)
 # oldclas
 # [1]  0  1  2  4  5  6  7  8  9 10 11 12 13 14 15 16
# http://www.landcover.org/data/lc/ see IGBP classification
# oldclas <-     c(0         1      2    4     5    6     7      8      9      10    11   12      13      14   15     16
oldclas_names <- c("Water","ENF","EBF","DBF","MF","CSH","OSH", "WSA", "SAV", "GRA","WET","CRO", "URBAN", "MOSAIC","Snow and ice", "Rock/Sand/Clay")
      newclas <- c(1,        2,    3,    3,    4,   5,    5,     5,      5,     6,    7,    8,      9,      8,      10, 11)
newclas_names <- c("Water","ENF", "DF", "MF", "SHB", "Grass", "Wetland", "Crop", "Urban", "Ice/Snow", "Barren")
#     newclas <- c(   "1",   "2",   "3",  "4",  "5",    "6",     "7",       "8",     "9",  "10",       "11")
newclas_color <- c(rgb(0,0,255/255, 1), rgb(0,102/255,0,1),rgb(0,178/255,0,1), rgb(0,178/255,178/255,1),
		   rgb(178/255,178/255,0,1),rgb(229/255,204/255,153/255,1),rgb(128/255,255/255,204/255,1),
		   rgb(255/255,179/255,204/255,1), rgb(255/255,0,0,1), rgb(1,1,1,1), rgb(229/255,229/255,204/255,1))
rclastab.df <- data.frame(oldclas, newclas)
lc_us2_rc2 = subs(lc_us2, rclastab.df)

library(rgdal)
flux.info.sp3 <- readOGR("C:/zhihua/dataset/flux2015jul", "flux.info.sp3")

png("F:/zhihua/dataset/results2/land cover with ec site.png",height = 2000, width = 3300, res = 300, units = "px")

plot(lc_us2_rc2, col = newclas_color,
      legend=FALSE, 
	 axes=FALSE,
	 box=FALSE)
	 
plot(usa.state, add = T)

legend(x = -75, y = 40, legend = newclas_names, fill = newclas_color, cex = 1, box.lwd = 0,box.col = "transparent",bg = "transparent")

plot(flux.info.sp3[which(flux.info.sp3$IGBP == "DBF" | flux.info.sp3$IGBP == "MF"), ], pch = 21, cex = 2, add = T)
plot(flux.info.sp3[which(flux.info.sp3$IGBP == "ENF"), ], pch = 22, cex = 2, add = T)
plot(flux.info.sp3[which(flux.info.sp3$IGBP == "GRA" | flux.info.sp3$IGBP == "CSH" | flux.info.sp3$IGBP == "OSH" 
                    | flux.info.sp3$IGBP == "SAV" | flux.info.sp3$IGBP == "WSA" | flux.info.sp3$IGBP == "CRO"), ], pch = 23, cex = 2, add = T)

text(x = coordinates(flux.info.sp3)[,1]+1, y = coordinates(flux.info.sp3)[,2]+1, flux.info.sp3$SITE_ID)

legend("bottomleft", 
	   # inset=0.05, 
	   legend = c("DBF/MF (n = 7)","ENF (n = 3)","GRA/Shurb/CRO (n = 7)" ),
	   horiz=F,
	   pch = c(21:23),
	   #col = 1:6,
	   #text.col = 1:6,
	   cex = 1,
	   box.col = "transparent",
	   bg = "transparent")	 

text()

dev.off()









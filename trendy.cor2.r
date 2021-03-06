rm(list = ls())

## 5/22/2017
## using remote sensing and co2 inversion method to calculate the relationship between GPP/ER with precipitation and temperature
## using the median/mean NEE estimate from four co2 inversion methods

library(rgdal)
library(raster)
library(rasterVis)
require(ncdf4)

usa.state = readOGR(dsn="F:\\zhihua\\dataset\\ecoregion", layer = "usa.state")
usa.state = usa.state[-which(usa.state$NAME_1 == "Hawaii" | usa.state$NAME_1 == "Alaska"),]

na.ext <- extent(-180,-48,15,85)

#### read into trendy, 
trendy.ra.mean = stack("F:/zhihua/dataset/trendy/Trendy.Ra.mean.grd")
trendy.rh.mean = stack("F:/zhihua/dataset/trendy/Trendy.Rh.mean.grd")

trendy.er.mean = trendy.ra.mean + trendy.rh.mean
trendy.gpp.mean = stack("F:/zhihua/dataset/trendy/Trendy.gpp.mean.grd")

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
er.annual = crop(trendy.er.mean, usa.state) # mean
er.annual = er.annual*mask

gpp.annual = trendy.gpp.mean*mask

temp = temp*mask
prep = prep*mask

temp = temp[[c(1:11)]]
prep = prep[[c(1:11)]]

# change NEON into grid
noen.grd = rasterize(neon.sp,mask, field = "DomainID")

######### section 4: calculate correlationship #########################################
er.df = zonal(er.annual, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:12)] 
gpp.df = zonal(gpp.annual, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:12)]
temp.df = zonal(temp, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:12)]
prep.df = zonal(prep, noen.grd, fun='mean', digits=0, na.rm=TRUE)[,c(2:12)]

temp.df.mn = zonal(calc(temp, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2]
prep.df.mn = zonal(calc(prep, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2]
er.df.mn = zonal(calc(er.annual, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2] 
gpp.df.mn = zonal(calc(gpp.annual, mean),noen.grd, fun='mean', digits=0, na.rm=TRUE)[,2]

#calculate anomaly
Ano = function(x){x-mean(x)} 
gpp.df.z = t(apply(gpp.df, 1, Ano))
temp.df.z = t(apply(temp.df, 1, Ano))
prep.df.z = t(apply(prep.df, 1, Ano))
er.df.z = t(apply(er.df, 1, Ano))

############################## bootstrapping  ###################################

# first value is coefficient, 2nd is intercept, 3rd is p value
# the order of the 5 value, coef for temp and prep, p value for temp and prep, r squared

boot.rs1 = function(dat.gpp, dat.temp, dat.prep, ln = 8, n = 100){ 
# for each site, at least 12 year was used to get correlationship
Coef.df = c()
NofSim = 0	

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
Coef.df = rbind(Coef.df, c(coef(lm1), lm11$r.squared, round(sqrt(mean(resid(lm11)^2)), 2)))

if(NofSim > n){break}
}

Coef.df = data.frame(Coef.df)
colnames(Coef.df) <- c("int","airtemp","prep","r2","rmse")
return(Coef.df)

}

boot.rs1(dat.gpp = gpp.df.z2[1,],dat.temp = temp.df.z2[1,],dat.prep = prep.df.z2[1,], n = 100)


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
coef.gpp = list() # 
coef.er = list() # 

for (i in 1:nrow(gpp.df.z2)){

gpp.df.z3 = gpp.df.z2[i,]
er.df.z3 = er.df.z2[i,]
temp.df.z3 = temp.df.z2[i,]
prep.df.z3 = prep.df.z2[i,]

coef.gpp1 = boot.rs1(dat.gpp = gpp.df.z3,dat.temp = temp.df.z3,dat.prep = prep.df.z3, n = 100)
coef.er1 = boot.rs1(dat.gpp = er.df.z3,dat.temp = temp.df.z3,dat.prep = prep.df.z3, n = 100)

coef.gpp[[i]] <- coef.gpp1
coef.er[[i]] <- coef.er1
		  
}

# before calculating, replace value > 90 percentile and < 10 percentile with NA
Rplc = function(x){
p1 = quantile(x, probs = c(0.1,0.9))
x[which(x < p1[1] | x > p1[2])] = NA
return(x)
}

pred.er.t = c()
pred.gpp.t = c()
dif.t = c()

for (i in 1:length(coef.er)){

x1 = coef.er[[i]][,3]
x1 = Rplc(x1)
pred.er.t = c(pred.er.t, x1)

x2 = coef.gpp[[i]][,3]
x2 = Rplc(x2)
pred.gpp.t = c(pred.gpp.t, x2)

dif.t = rbind(dif.t, c(x2 - x1))
}
	
dif.t = data.frame(prep = prep.df.mn[-4], mean = 100*apply(dif.t,1,mean, na.rm = TRUE),
			                              sd = 100*apply(dif.t,1,sd, na.rm = TRUE))	

plot(mean~prep, data = dif.t, cex = 3)
abline(lm(mean~prep, data = dif.t))
										  
dif.t = dif.t[complete.cases(dif.t),]
				
ggplot(dif.t, aes(x=prep, y=coef1)) + 
	geom_point(shape=1, cex = 5) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
	# geom_smooth(method="auto",   # Add linear regression lines	
               se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) 
	
pred.er.t2 = data.frame(coef1=pred.er.t, prep = rep(prep.df.mn[-4], each = 101), flux = "TER")					
pred.gpp.t2 = data.frame(coef1=pred.gpp.t, prep = rep(prep.df.mn[-4], each = 101), flux = "GPP")					

d1 = rbind(pred.er.t2, pred.gpp.t2)	
d1 = d1[complete.cases(d1),]
d1$coef1 = d1$coef1*100
	
library(ggplot2)		
		
p1 = ggplot(d1, aes(x=prep, y=coef1, color=flux)) + 
	geom_point(shape=1, cex = 1) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) +
    coord_cartesian(xlim=c(100, 1500), ylim=c(-50, 150))	+ 	 
  #  ylab(expression(paste(beta ["temporal"]))) + 
 ylab(expression("" ~ delta ^{t} ~ "")) + 
	xlab("MAP (mm)") + # Set axis labels
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	theme(legend.position=c(.5, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 12)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) 
   
print(p1)
ggsave("F:/zhihua/dataset/results2/fig2.trendy1-2.png", width = 4, height = 3, units = "in")
	
#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.mn[-4], er = er.df.mn[-4],airtemp = temp.df.mn[-4],prep = prep.df.mn[-4])

list.gpp = boot.ec2.gpp(dat = dat.df.annual.mean, ln = 13, n = 100)
reg.abs.mod.pred1.gpp.list = list.gpp[[1]]
reg.abs.mod.pred1.gpp.model.coef = list.gpp[[2]]
delt.gpp = list.gpp[[3]]
rmse.gpp = list.gpp[[4]]
delt.gpp.temp = list.gpp[[5]]

list.er = boot.ec2.er(dat = dat.df.annual.mean, ln = 13, n = 100)
reg.abs.mod.pred1.er.list = list.er[[1]]
reg.abs.mod.pred1.er.model.coef = list.er[[2]]
delt.er = list.er[[3]]
rmse.er = list.er[[4]]
delt.er.temp = list.er[[5]]

# calculate the plot spatial coefficient
prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)


delt.gpp_1 = apply(delt.gpp, 1, Rplc)
delt.gpp_2 = as.vector(delt.gpp_1)

delt.er_1 = apply(delt.er, 1, Rplc)
delt.er_2 = as.vector(delt.er_1)


delt.spatial = delt.gpp_1 - delt.er_1
delt.spatial = data.frame(prep = prep.grd[-1], mean = 2*apply(delt.spatial, 2, mean, na.rm = TRUE), 
			  sd = 2*apply(delt.spatial, 2, sd, na.rm = TRUE))

plot(mean~prep, data = delt.spatial, cex = 3)


delt.gpp_3 = data.frame(coef1=delt.gpp_2, prep = rep(prep.grd[-1], each = 101), flux = "GPP")		
delt.er_3 = data.frame(coef1=delt.er_2, prep = rep(prep.grd[-1], each = 101), flux = "TER")		
d2 = rbind(delt.er_3,delt.gpp_3)

d2 = d2[complete.cases(d2),]
d2$coef1 = d2$coef1*2		
library(ggplot2)		
		
p2 = ggplot(d2, aes(x=prep, y=coef1, color=flux)) + 
	geom_point(shape=1, cex = 1) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) + 
	coord_cartesian(xlim=c(100, 1500), ylim=c(0, 250))	+ 	 
  #  ylab(expression(paste(beta ["Spatial"]))) + 
 ylab(expression("" ~ delta ^{s} ~ "")) + 
	xlab("MAP (mm)") + # Set axis labels
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	theme(legend.position=c(.5, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 12)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) 

print(p2)
ggsave("F:/zhihua/dataset/results2/fig2.trendy2-2.png", width = 4, height = 3, units = "in")

		   
d3 = rbind(data.frame(d1, doman = "Temporal"),data.frame(d2, doman = "Spatial"))			   
ggplot(d3, aes(x=prep, y=coef1, color=flux)) + 
	geom_point(shape=1, cex = 3) +
	facet_grid(.~doman)
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) # Extend regression lines

write.csv(d3, file = "F:/zhihua/dataset/results2/TRENDY.sensitivity2.csv")
write.csv(dif.t, file = "F:/zhihua/dataset/results2/delt.temporal.trendy.csv")
write.csv(delt.spatial, file = "F:/zhihua/dataset/results2/delt.spatial.trendy.csv")

##plot sensitivity to temperature
####### USE ggplot2 to map 
prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)

# calculate the plot spatial coefficient
delt.gpp_1t = apply(delt.gpp.temp, 1, Rplc)
delt.gpp_2t = as.vector(delt.gpp_1t)

delt.er_1t = apply(delt.er.temp, 1, Rplc)
delt.er_2t = as.vector(delt.er_1t)

delt.gpp_3t = data.frame(coef1=delt.gpp_2t, prep = rep(airtemp.grd[-1], each = 101), flux = "GPP")		
delt.er_3t = data.frame(coef1=delt.er_2t, prep = rep(airtemp.grd[-1], each = 101), flux = "ER")		
d2t = rbind(delt.er_3t, delt.gpp_3t)

d2t = d2t[complete.cases(d2t),]
library(ggplot2)		
		
p2t = ggplot(d2t, aes(x=prep, y=coef1, color=flux)) + 
	geom_point(shape=1, cex = 1) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) + 
	coord_cartesian(xlim=c(1, 25), ylim=c(-250, 250))	+ 	 
    ylab(expression(paste(beta ["Spatial"]))) + 
	xlab("Temperature (degree)") + # Set axis labels
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	theme(legend.position=c(.5, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 12)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=12),axis.text.x  = element_text(colour="black",size=10))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=12),axis.text.y  = element_text(colour="black",size=12))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) 

# plot the difference between spatial sensitivity between spatial and temporal
delt.spatial.temp = delt.gpp_1t - delt.er_1t
delt.spatial.temp = data.frame(prep = airtemp.grd[-1], mean = apply(delt.spatial.temp, 2, mean, na.rm = TRUE), 
			  sd = apply(delt.spatial.temp, 2, sd, na.rm = TRUE))

plot(mean~prep, data = delt.spatial.temp, cex = 3)

# temporal

pred.er.t.temp = c()
pred.gpp.t.temp = c()
delt.temporal.temp = c()
for (i in 1:length(coef.er)){

x1 = coef.er[[i]][,2]
x1 = Rplc(x1)
pred.er.t.temp = c(pred.er.t.temp, x1)

x2 = coef.gpp[[i]][,2]
x2 = Rplc(x2)
pred.gpp.t.temp = c(pred.gpp.t.temp, x2)
	
delt.temporal.temp = rbind(delt.temporal.temp, c(x2 - x1))	
}
					
pred.er.t2.temp = data.frame(coef1=pred.er.t.temp, prep = rep(dat.df.annual.mean[,c("airtemp")], each = 101), flux = "ER")					
pred.gpp.t2.temp = data.frame(coef1=pred.gpp.t.temp, prep = rep(dat.df.annual.mean[,c("airtemp")], each = 101), flux = "GPP")					

d1.temp = rbind(pred.er.t2.temp, pred.gpp.t2.temp)	
d1.temp = d1.temp[complete.cases(d1.temp),]
	
library(ggplot2)		
		
p1.temp = ggplot(d1.temp, aes(x=prep, y=coef1, color=flux)) + 
	geom_point(shape=1, cex = 1) +
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method="lm",   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) +
    coord_cartesian(xlim=c(1, 25), ylim=c(-50, 150))	+ 	 
    ylab(expression(paste(beta ["temporal"]))) + 
	xlab("Temperature (degree)") + # Set axis labels
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	theme(legend.position=c(.5, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 12)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=12),axis.text.x  = element_text(colour="black",size=10))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=12),axis.text.y  = element_text(colour="black",size=12))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) 

	
delt.temporal.temp = data.frame(prep = dat.df.annual.mean[,c("airtemp")], mean = apply(delt.temporal.temp, 1, mean, na.rm = TRUE), 
			   sd = apply(delt.temporal.temp, 1, sd, na.rm = TRUE)) 


plot(mean~prep, data = delt.temporal.temp)	

d3t1 = rbind(data.frame(d1.temp, doman = "Temporal"),data.frame(d2t, doman = "Spatial"))			   

write.csv(d3t1, file = "F:/zhihua/dataset/results2/TRENDY.sensitivity2.temp.csv")
write.csv(delt.temporal.temp, file = "F:/zhihua/dataset/results2/delt.temporal.trendy.temp.csv")
write.csv(delt.spatial.temp, file = "F:/zhihua/dataset/results2/delt.spatial.trendy.temp.csv")

print(p1.temp)
ggsave("F:/zhihua/dataset/results2/fig2.trendy1.temp.png", width = 4, height = 3, units = "in")
print(p2t)
ggsave("F:/zhihua/dataset/results2/fig2.trendy2.temp.png", width = 4, height = 3, units = "in")


# statistics, temporal
# x1 = mean(delt.gpp_2t, na.rm = TRUE);x1sd = sd(delt.gpp_2t, na.rm = TRUE)
# x2 = mean(delt.er_2t, na.rm = TRUE);x2sd = sd(delt.er_2t, na.rm = TRUE)

r1 = c()
r2 = c()
r3 = c()
r4 = c()

r5 = c()
r6 = c()
r7 = c()
r8 = c()


for (i in 1:length(coef.er)){
r1 = c(r1, coef.gpp[[i]][,4]) #r-squired
r2 = c(r2, coef.er[[i]][,4])
r3 = c(r3, coef.gpp[[i]][,5]) #rmse
r4 = c(r4, coef.er[[i]][,5])
r5 = c(r5, 100*coef.gpp[[i]][,3]) #gpp prep sensitivity
r6 = c(r6, 100*coef.er[[i]][,3]) #er prep sensitivity
r7 = c(r7, coef.gpp[[i]][,2]) #gpp temperature sensitivity
r8 = c(r8, coef.er[[i]][,2]) #er temperature sensitivity

}

x31 = mean(r5, na.rm = TRUE);x31sd = sd(r5, na.rm = TRUE)
x41 = mean(r6, na.rm = TRUE);x41sd = sd(r6, na.rm = TRUE)
x51 = mean(r7, na.rm = TRUE);x51sd = sd(r7, na.rm = TRUE)
x61 = mean(r8, na.rm = TRUE);x61sd = sd(r8, na.rm = TRUE)

x3 = mean(r1, na.rm = TRUE);x3sd = sd(r1, na.rm = TRUE)
x4 = mean(r2, na.rm = TRUE);x4sd = sd(r2, na.rm = TRUE)
x5 = mean(r3, na.rm = TRUE);x5sd = sd(r3, na.rm = TRUE)
x6 = mean(r4, na.rm = TRUE);x6sd = sd(r4, na.rm = TRUE)


# statistics, spatial
x7 = mean(delt.gpp_2, na.rm = TRUE);x7sd = sd(delt.gpp_2, na.rm = TRUE)
x8 = mean(delt.er_2, na.rm = TRUE);x8sd = sd(delt.er_2, na.rm = TRUE)
  # r2
x9 = mean(rmse.gpp[,4]); x9sd = sd(rmse.gpp[,4])
x10 = mean(rmse.er[,4]); x10sd = sd(rmse.er[,4])
  #rmse
x11 = mean(rmse.gpp[,1]); x11sd = sd(rmse.gpp[,1])
x12 = mean(rmse.er[,1]); x12sd = sd(rmse.er[,1])

mean(delt.gpp)*2; sd(delt.gpp)*2
mean(delt.er)*2; sd(delt.er)*2

mean(delt.gpp.temp); sd(delt.gpp.temp)
mean(delt.er.temp); sd(delt.er.temp)

# r-squired
mean(rmse.gpp[,4]);sd(rmse.gpp[,4])
mean(rmse.er[,4]);sd(rmse.er[,4])

# rmse
mean(rmse.gpp[,1]);sd(rmse.gpp[,1])
mean(rmse.er[,1]);sd(rmse.er[,1])

data.frame(name = c("T.p", "T.t", "T.rmse", "T.r2","S.p", "S.t", "S.rmse", "S.r2"),
           gppmean = c(mean(r5, na.rm = TRUE), mean(r7, na.rm = TRUE), mean(r3, na.rm = TRUE), mean(r1, na.rm = TRUE),
                       mean(delt.gpp)*2,mean(delt.gpp.temp), mean(rmse.gpp[,1]), mean(rmse.gpp[,4])),
		   gppsd = c(sd(r5, na.rm = TRUE), sd(r7, na.rm = TRUE), sd(r3, na.rm = TRUE), sd(r1, na.rm = TRUE),
                       sd(delt.gpp)*2,sd(delt.gpp.temp), sd(rmse.gpp[,1]), sd(rmse.gpp[,4])),
		   ermean = c(mean(r6, na.rm = TRUE), mean(r8, na.rm = TRUE), mean(r4, na.rm = TRUE), mean(r2, na.rm = TRUE),
                       mean(delt.er)*2, mean(delt.er.temp), mean(rmse.er[,1]), mean(rmse.er[,4])),
		   ersd = c(sd(r6, na.rm = TRUE), sd(r8, na.rm = TRUE), sd(r4, na.rm = TRUE), sd(r2, na.rm = TRUE),
                       sd(delt.er)*2, sd(delt.er.temp), sd(rmse.er[,1]), sd(rmse.er[,4])))

    name     gppmean       gppsd      ermean        ersd
1    T.p  45.2701440 31.00044691  27.4067538 18.49194341
2    T.t   9.3286961 52.30448346  38.8916796 36.38927090
3 T.rmse  37.1397339 18.09029305  27.3492884 14.52325176
4   T.r2   0.5216908  0.20405700   0.5209353  0.22354785
5    S.p 111.6168958 24.13136325 113.5866939 27.28597975
6    S.t  16.6337353 25.02555875  22.7403048 33.65501057
7 S.rmse 141.7775248 13.78281832 153.1353465 13.15385932
8   S.r2   0.9033663  0.02055613   0.8971287  0.02046151
>    
					   
save.image("F:/zhihua/dataset/results2/trendy.cor2.RData")

load("F:/zhihua/dataset/results2/trendy.cor2.RData")

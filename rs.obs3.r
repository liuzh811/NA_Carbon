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

writeRaster(ensm.nee1,"F:/zhihua/dataset/co2_inversion/NEE.annual.mean.grd",overwrite=TRUE) 
writeRaster(ensm.nee2,"F:/zhihua/dataset/co2_inversion/NEE.annual.median.grd",overwrite=TRUE) 

######### section 2: read into MODIS GPP and aggregate into 1 degree resolution #########################################

ensm.nee1 = stack("F:/zhihua/dataset/co2_inversion/NEE.annual.mean.grd")

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

############################## bootstrapping  ###################################

# first value is coefficient, 2nd is intercept, 3rd is p value
# the order of the 5 value, coef for temp and prep, p value for temp and prep, r squared

boot.rs1 = function(dat.gpp, dat.temp, dat.prep, ln = 14, n = 100){ 
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
Coef.df = rbind(Coef.df,  c(coef(lm1), lm11$r.squared, round(sqrt(mean(resid(lm11)^2)), 2)))
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
delt.temporal = c()

for (i in 1:length(coef.er)){

x1 = coef.er[[i]][,3]
x1 = Rplc(x1)
pred.er.t = c(pred.er.t, x1)

x2 = coef.gpp[[i]][,3]
x2 = Rplc(x2)
pred.gpp.t = c(pred.gpp.t, x2)
	
delt.temporal = rbind(delt.temporal, c(x2 - x1))
}
	

pred.er.t2 = data.frame(coef1=pred.er.t, prep = rep(prep.df.mn[-4], each = 101), flux = "ER")					
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
    ylab(expression(paste(beta ["temporal"]))) + 
	xlab("MAP (mm)") + # Set axis labels
    # ggtitle("Average bill for 2 people") +     # Set title
    theme_bw() + 
	theme(legend.position=c(.5, .8)) + 	
	theme(legend.title=element_blank()) +
	theme(legend.text = element_text(size = 12)) +
	theme(axis.title.x = element_text(face="bold", colour="black", size=12),axis.text.x  = element_text(colour="black",size=10))+
    theme(axis.title.y = element_text(face="bold", colour="black", size=12),axis.text.y  = element_text(colour="black",size=12))+
    theme(strip.text.x = element_text(size=12))+
    theme(strip.text.y = element_text(size=12)) 
  

# plot spatial sensitivity difference
delt.temporal = data.frame(prep = dat.df.annual.mean[,c("prep")], mean = 100*apply(delt.temporal, 1, mean, na.rm = TRUE), 
sd = 100*apply(delt.temporal, 1, sd, na.rm = TRUE)) 

plot(mean~prep, data = delt.temporal, cex = 3)
abline(lm(mean~prep, data = delt.temporal))
  	
#spatial relationship
dat.df.annual.mean = data.frame(gpp = gpp.df.mn[-4], nee = nee.df.mn[-4],er = er.df.mn[-4],airtemp = temp.df.mn[-4],prep = prep.df.mn[-4])

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
delt.gpp_1 = apply(delt.gpp, 1, Rplc)
delt.gpp_2 = as.vector(delt.gpp_1)

delt.er_1 = apply(delt.er, 1, Rplc)
delt.er_2 = as.vector(delt.er_1)

prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)

delt.gpp_3 = data.frame(coef1=delt.gpp_2, prep = rep(prep.grd[-1], each = 101), flux = "GPP")		
delt.er_3 = data.frame(coef1=delt.er_2, prep = rep(prep.grd[-1], each = 101), flux = "ER")		
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
    ylab(expression(paste(beta ["Spatial"]))) + 
	xlab("MAP (mm)") + # Set axis labels
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
delt.spatial = delt.gpp_1 - delt.er_1
delt.spatial = data.frame(prep = prep.grd[-1], mean = 2*apply(delt.spatial, 2, mean, na.rm = TRUE), 
			  sd = 2*apply(delt.spatial, 2, sd, na.rm = TRUE))

plot(mean~prep, data = delt.spatial, cex = 3)


print(p1)
ggsave("F:/zhihua/dataset/results2/fig2.rs1.png", width = 4, height = 3, units = "in")
print(p2)
ggsave("F:/zhihua/dataset/results2/fig2.rs2.png", width = 4, height = 3, units = "in")


d3 = rbind(data.frame(d1, doman = "Temporal"),data.frame(d2, doman = "Spatial"))			   
ggplot(d3, aes(x=prep, y=coef1, color=flux)) + 
	geom_point(shape=1, cex = 3) +
	facet_grid(.~doman)
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) # Extend regression lines



write.csv(d3, file = "F:/zhihua/dataset/results2/RS.sensitivity2.csv")
write.csv(delt.temporal, file = "F:/zhihua/dataset/results2/delt.temporal.rs.csv")
write.csv(delt.spatial, file = "F:/zhihua/dataset/results2/delt.spatial.rs.csv")

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

save.image("F:/zhihua/dataset/results2/rs.obs3.RData")

















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

png("F:/zhihua/dataset/results2/rs.sensitivity_ci-2.png",height = 2500, width = 3000, res = 300, units = "px")

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
points(dat.df.annual.mean$prep, dat.df.annual.mean$er, col = "red")

# plot spatial sensitivity
library(oce)
plotInset(750,0, 1450, 800,
          expr= {
		  
plot(0, xlim = c(0, 1500), ylim = c(25,75)*2, pch='', xlab = "Precipation (mm)", 
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

abline(v = c(900-100, 900, 900 + 100), lty = c(2,1,2))

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
write.csv(D3, "F:/zhihua/dataset/results2/RS.sensitivity.csv")				

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
write.csv(rmse1, "F:/zhihua/dataset/results2/RS.rmse1.csv")

colnames(reg.abs.mod.pred1.er.model.coef) = c("int","p2","p") 
colnames(reg.abs.mod.pred1.gpp.model.coef) = c("int","p2","p") 
coef1 = rbind(data.frame(reg.abs.mod.pred1.gpp.model.coef, flux = "gpp"),data.frame(reg.abs.mod.pred1.er.model.coef, flux = "er"))
write.csv(coef1, "F:/zhihua/dataset/results2/RS.coef1.csv")

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
## only select number of year > 4, homogenous land cover type and non-wetland 
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

##change to spatial point
flux.info2 = fluxsite.df
flux.info.sp2 = fluxsite.df						 
coordinates(flux.info.sp2) <- ~LOCATION_LONG + LOCATION_LAT
projection(flux.info.sp2) <- "+proj=longlat +datum=WGS84"

#get the name of fluxnet from FLUXNET2015 dataset
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
					
#######################################################################################
#select the following site
SITE_ID_sel = data.frame(SITE_ID = c("US-Cop", "US-ARM", "US-MMS","US-Ne3","US-NR1","US-Oho","US-PFa","US-Ha1","US-Me2","US-SRM","US-Syv","US-Var","US-Whs","US-Wkg"))
SITE_ID_sel_yrexl = data.frame(SITE_ID = c("US-ARM", "US-Ha1",rep("US-Me2",2),rep("US-Var",4)),
                                year = c(2012, 2010, c(2013,2014), c(2001,2005, 2006,2011)))

#select whose within the database
SITE_ID_sel2 = data.frame(SITE_ID = c("US-MOz","US-Slt","US-Vcm","US-Vcp","US-Wjs"))
SITE_ID_sel_yrex2 = data.frame(SITE_ID = c("US-Slt"),
                                year = c(2007))


SITE_select = rbind(SITE_ID_sel, SITE_ID_sel2)
SITE_exclude = rbind(SITE_ID_sel_yrexl, SITE_ID_sel_yrex2)			
################################################################################################################

##################### read into annual data, and compile it into one data sheet #################################
dat.annual.df = c()

for (i in 1:nrow(flux.info2)){

if (flux.info2$DataSource[i] == "FLUXNET2015"){
file.name <- fluxnet.fn[which(substr(fluxnet.fn, 5,10) == flux.info2$SITE_ID[i])]
WCr.df = data.frame(read.csv(paste("F:/zhihua/dataset/flux2015jul/download/", 
                                             paste(substr(file.name,1,30),"YY",substr(file.name,30,44),"csv", sep = ""), 
											 sep = "")))
WCr.df[WCr.df == -9999] = NA	
## gpp and er was calculated as mean value of two methods (night time and day time method)	
WCr.df = data.frame(site_id = substr(file.name,5,10),
					year = WCr.df$TIMESTAMP,
					#month =  WCr.df$TIMESTAMP - floor(WCr.df$TIMESTAMP/100)*100,
					nee = WCr.df$NEE_VUT_50,
					# gpp = WCr.df$GPP_NT_VUT_50,
					# er = WCr.df$RECO_NT_VUT_50,
					# gpp = WCr.df$GPP_DT_VUT_50,
					# er = WCr.df$RECO_DT_VUT_50,
					gpp = apply(data.frame(WCr.df$GPP_NT_VUT_50,WCr.df$GPP_DT_VUT_50), 1, mean, na.rm = TRUE),
					er = apply(data.frame(WCr.df$RECO_NT_VUT_50,WCr.df$RECO_DT_VUT_50), 1, mean, na.rm = TRUE),		
					airtemp = WCr.df$TA_F,
					prep = WCr.df$P_F,					
					vpd = WCr.df$VPD_F,
					par = 0.45*WCr.df$SW_IN_F, #PAR = 0.45*incoming shortwave length
					#soiltemp = WCr.df$TS_F_MDS_1,
					#swc = WCr.df$SWC_F_MDS_1,
					igbp = flux.info2$IGBP[which(flux.info2$SITE_ID == substr(file.name,5,10))],
					country = substr(file.name,5,6)
					)


} else {


WCr.df = data.frame(read.csv(paste("F:/zhihua/dataset/ameriflux/level2/unprocessed/", flux.info2$SITE_ID[i],".annual.csv", sep = "")))

if("NEE_st_fMDS" %in% colnames(WCr.df)){ #ckeck if NEE_st_fMDS if one of the collumn name in the database
								
WCr.df = data.frame(site_id = as.character(flux.info2$SITE_ID[i]),
					year = WCr.df$Year,
					nee = WCr.df$NEE_st_fMDS,
					gpp = WCr.df$GPP_st_MDS,
					er = WCr.df$Reco_st,
					airtemp = WCr.df$Ta_f,
					prep = WCr.df$Precip,					
					vpd = WCr.df$VPD_f,
					par = 0.45*WCr.df$Rg_f, #PAR = 0.45*incoming shortwave length
					igbp = as.character(flux.info2$IGBP[i]),
					country = substr(as.character(flux.info2$SITE_ID[i]),1,2)
					)
} else {

WCr.df = data.frame(site_id = as.character(flux.info2$SITE_ID[i]),
					year = WCr.df$Year,
					nee = WCr.df$NEE_WithUstar_f,
					gpp = WCr.df$GPP_WithUstar_f,
					er = WCr.df$Reco_WithUstar,
					airtemp = WCr.df$Tair_f,
					prep = WCr.df$Prec,					
					vpd = WCr.df$VPD,
					par = 0.45*WCr.df$PotRad_WithUstar, #PAR = 0.45*incoming shortwave length
					igbp = as.character(flux.info2$IGBP[i]),
					country = substr(as.character(flux.info2$SITE_ID[i]),1,2)
					)
					
}


}

dat.annual.df = rbind(dat.annual.df, WCr.df)

print(paste("Finish for ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

#remove bad years, based on site.ex.df
ind.exd = c()
for(i in 1:nrow(site.ex.df)){
ind.exd = c(ind.exd, which(as.character(dat.annual.df$site_id) == as.character(site.ex.df$SITE_ID[i]) & dat.annual.df$year == site.ex.df$year[i]))
}

dat.annual.df = dat.annual.df[-ind.exd, ]

#select good site and years
dat.annual.df = merge(dat.annual.df, SITE_select, by.x = "site_id", by.y = "SITE_ID")
#remove bad site/years

ind.exd = c()
for(i in 1:nrow(SITE_exclude)){
ind.exd = c(ind.exd, which(as.character(dat.annual.df$site_id) == as.character(SITE_exclude$SITE_ID[i]) & dat.annual.df$year == SITE_exclude$year[i]))
}

dat.annual.df = dat.annual.df[-ind.exd, ]
 

dat.annual.df = dat.annual.df[with(dat.annual.df, order(site_id, year)), ]

dat.df = dat.annual.df
dat.df[dat.df == -9999] = NA

dat.df$nee = -1*dat.df$nee

dat.df$er2 = dat.df$er/dat.df$gpp
dat.df$er2[which(dat.df$er2 < 0)] = NA
dat.df$er2[which(dat.df$er2 > 2)] = NA

dat.df$b = dat.df$nee/dat.df$gpp
dat.df$b[which(dat.df$b < -1)] = NA

dat.df$r = dat.df$nee/dat.df$er

###
dat.df.annual.mean = data.frame(aggregate(cbind(nee, gpp, er, airtemp, prep, vpd, par, er2, b, r) ~ site_id, data = dat.df, mean, na.rm = TRUE))
dat.df.annual.mean = merge(dat.df.annual.mean, data.frame(site_id = flux.info2$SITE_ID, igbp = flux.info2$IGBP,country = substr(flux.info2$SITE_ID,1,2)),
                             by = "site_id", all.x = TRUE)

dat.df.annual.sd = data.frame(aggregate(cbind(nee, gpp, er, airtemp, prep, vpd, par, er2, b,r) ~ site_id, data = dat.df, sd, na.rm = TRUE))
###################################################################################################################

###########################  calculate absolute and relative annual anomaly   ######################################
#step 2 and 3: 
dat.df.ano.abs2 = c()
dat.df.ano.rel2 = c()

for (i in 1:length(dat.df.annual.mean$site_id)){
#select the raw data

tmp = dat.df[which(toupper(dat.df$site_id) == toupper(dat.df.annual.mean$site_id[i])),]
tmp[which(tmp$nee == -9999),] = NA
tmp$b[which(is.infinite(tmp$b))] = NA

tmp11 = tmp
tmp12 = tmp

#select annual mean
tmp2 = dat.df.annual.mean[which(toupper(dat.df.annual.mean$site_id) == toupper(dat.df.annual.mean$site_id[i])),]

#calculate absolute anomaly
tmp11[, c("nee", "gpp", "er", "airtemp", "prep", "vpd", "par","er2","b","r")] = sweep(
															 tmp[, c("nee", "gpp", "er", "airtemp", "prep", "vpd", "par","er2","b","r")], 
															 MARGIN=2, 
                                                             STATS = as.numeric(tmp2[, c("nee", "gpp", "er", "airtemp", "prep", "vpd", "par","er2","b","r")]), 
															 FUN="-")
#calculate relative anomaly													 
tmp12[, c("nee", "gpp", "er", "airtemp", "prep", "vpd", "par","er2","b","r")] = sweep(
															 tmp11[, c("nee", "gpp", "er", "airtemp", "prep", "vpd", "par","er2","b","r")], #here is the absoulate anomaly
															 MARGIN=2, 
                                                             STATS = as.numeric(tmp2[, c("nee", "gpp", "er", "airtemp", "prep", "vpd", "par","er2","b","r")]), 
															 FUN="/")
	
dat.df.ano.abs2 = rbind(dat.df.ano.abs2, tmp11)
dat.df.ano.rel2 = rbind(dat.df.ano.rel2, tmp12)

}

######################################################################################################################


############ step 4: temporal sensitivity
############ calculate the GPP and Rainfall/temperature relationship for EACH individual site

### use Mento Carlo methods to estimate the temporal sensitivity, BUT NOT USED FOR PLOTTING
reg.abs.gpp = c()
reg.rel.gpp = c()

reg.abs.mod.gpp = list()
reg.rel.mod.gpp = list()

reg.abs.mod.pred.gpp = list()
reg.rel.mod.pred.gpp = list()

coef.gpp = list() # each element stores mento carlo coeficient for each site

N = 100 #number of Mento Carlo run

for(i in 1:length(dat.df.annual.mean$site_id)) {

tmp = dat.df.ano.abs2[which(toupper(dat.df.ano.abs2$site_id) == toupper(dat.df.annual.mean$site_id[i])),]
lm1 = lm(gpp~airtemp+prep, data = tmp,na.action=na.exclude)
#lm1 = lm(nee~gpp, data = tmp,na.action=na.exclude)

#construct new data for prediction
dat.new1 = data.frame(airtemp = tmp$airtemp, prep = mean(tmp$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(tmp$airtemp, na.rm = TRUE), prep = tmp$prep)

pred1 = predict(lm1, newdata = dat.new1)
pred2 = predict(lm1, newdata = dat.new2)

reg.abs.mod.pred.gpp[[i]] <- data.frame(airtemp = tmp$airtemp, predt = pred1, prep = tmp$prep, predp = pred2)

## mento carlo from here ##
# get residual
resd1 = lm1$residuals
resd1.mean = mean(resd1)
resd1.sd = sd(resd1)


NofSim = 0	
coef.gpp.df = c()
repeat {
NofSim = NofSim+1

resd1.value = rnorm(nrow(tmp), mean = resd1.mean, sd = resd1.sd )
tmp$gpp1 = tmp$gpp + resd1.value

lm11 = lm(gpp1~airtemp+prep, data = tmp,na.action=na.exclude)
coef.gpp.df = rbind(coef.gpp.df, coef(lm11))

if(NofSim > N){break}
}

coef.gpp[[i]] <- coef.gpp.df
## mento carlo stop here ##


tmp = dat.df.ano.rel2[which(toupper(dat.df.ano.rel2$site_id) == toupper(dat.df.annual.mean$site_id[i])),]
lm2 = lm(gpp~airtemp+prep, data = tmp)

#construct new data for prediction
pred1 = predict(lm2, newdata = dat.new1)
pred2 = predict(lm2, newdata = dat.new2)

reg.rel.mod.pred.gpp[[i]] <- data.frame(airtemp = tmp$airtemp, predt = pred1, prep = tmp$prep, predp = pred2)


reg.abs.gpp <- rbind(reg.abs.gpp, c(summary(lm1)$coefficients[,1], 
						  summary(lm1)$coefficients[,4],
			              summary(lm1)$r.squared))
reg.rel.gpp <- rbind(reg.rel.gpp, c(summary(lm2)$coefficients[,1], 
						  summary(lm2)$coefficients[,4],
			              summary(lm2)$r.squared))

reg.abs.mod.gpp[[i]] <- lm1						  
reg.rel.mod.gpp[[i]] <- lm2							  
}


##################  plot to see results ################################################################################
#plot partial plot [hold airtemp or prep as constant, plot nee response]
flux.info3 = dat.df.annual.mean[,c("site_id", "igbp")]
flux.info3$IGBP2 = factor(flux.info3$igbp)
levels(flux.info3$IGBP2) <- c("CRO", "DBF","ENF","GRA","MF", "SHB","SHB")

# plot biome level coefficient							
reg.abs.gpp = data.frame(reg.abs.gpp)
colnames(reg.abs.gpp) <- c("Int.coef","airtemp.coef", "prep.coef","Int.p","airtemp.p", "prep.p","r2")				
reg.abs.gpp = data.frame(reg.abs.gpp, dat.df.annual.mean)													
reg.abs.gpp$IGBP2 = factor(reg.abs.gpp$igbp)
levels(reg.abs.gpp$IGBP2) <- c("CRO", "DBF","ENF","GRA","MF", "SHB","SHB")
							
boxplot(airtemp.coef~IGBP2, data = reg.abs.gpp)
boxplot(prep.coef~IGBP2, data = reg.abs.gpp)

reg.abs.gpp$dif = reg.abs.gpp$airtemp.coef + reg.abs.gpp$prep.coef
boxplot(dif~IGBP2, data = reg.abs.gpp, subset = (IGBP2 != "Mixed Forest" | IGBP2 != "Wetland"))

#################################################################################################

############ step 5: spatial sensitivity
## two methods were used here, leave one out methods AND Mento Carlo Method
library(MASS)

########################## 5.1 LEAVE ONE OUT METHOD   ########################
## leave one out methods start from here
reg.abs.mod.pred1.gpp.list = list()
delt.gpp = c()

for (i in 1:nrow(dat.df.annual.mean)){

dat.df.annual.mean1 = dat.df.annual.mean[-i,]

lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)
lm1.gpp = stepAIC(lm1.gpp)

#construct new data for prediction

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

pred1 = predict(lm1.gpp, newdata = dat.new1)
pred2 = predict(lm1.gpp, newdata = dat.new2)

reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)
reg.abs.mod.pred1.gpp.list[[i]] <- reg.abs.mod.pred1.gpp

delt.gpp1 = diff(reg.abs.mod.pred1.gpp$predp)
delt.gpp = cbind(delt.gpp, delt.gpp1)

}
## leave one out methods ends here

##########################  5.2 MENTO CARLO METHOD   ########################
## mento carlo start from here ##
reg.abs.mod.pred1.gpp.list = list()
delt.gpp = c()

# fit model
lm1.gpp = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.gpp = stepAIC(lm1.gpp)

# get residual
resd1 = lm1.gpp$residuals
resd1.mean = mean(resd1)
resd1.sd = sd(resd1)

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

NofSim = 0	
N = 100
repeat {
NofSim = NofSim+1

resd1.value = rnorm(nrow(dat.df.annual.mean), mean = resd1.mean, sd = resd1.sd )
dat.df.annual.mean$gpp1 = dat.df.annual.mean$gpp + resd1.value

lm1.gpp1 = lm(gpp1~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.gpp1 = stepAIC(lm1.gpp1)

pred1 = predict(lm1.gpp1, newdata = dat.new1)
pred2 = predict(lm1.gpp1, newdata = dat.new2)

reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)
reg.abs.mod.pred1.gpp.list[[NofSim]] <- reg.abs.mod.pred1.gpp

delt.gpp1 = diff(reg.abs.mod.pred1.gpp$predp)
delt.gpp = cbind(delt.gpp, delt.gpp1)


if(NofSim > N){break}
}

## mento carlo end here ##

###############################################################################################

#################################################################################
## repeat step 4 and 5
## step 4: calculate the NEE and Rainfall/temperature relationship for each individual site
reg.abs.er = c()
reg.rel.er = c()

reg.abs.mod.er = list()
reg.rel.mod.er = list()

reg.abs.mod.pred.er = list()
reg.rel.mod.pred.er = list()

coef.er = list() # each element stores mento carlo coeficient for each site

for(i in 1:length(dat.df.annual.mean$site_id)) {

tmp = dat.df.ano.abs2[which(toupper(dat.df.ano.abs2$site_id) == toupper(dat.df.annual.mean$site_id[i])),]
lm1 = lm(er~airtemp+prep, data = tmp,na.action=na.exclude)

#construct new data for prediction
dat.new1 = data.frame(airtemp = tmp$airtemp, prep = mean(tmp$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(tmp$airtemp, na.rm = TRUE), prep = tmp$prep)

pred1 = predict(lm1, newdata = dat.new1)
pred2 = predict(lm1, newdata = dat.new2)

reg.abs.mod.pred.er[[i]] <- data.frame(airtemp = tmp$airtemp, predt = pred1, prep = tmp$prep, predp = pred2)

## mento carlo from here ##
# get residual
resd1 = lm1$residuals
resd1.mean = mean(resd1)
resd1.sd = sd(resd1)

NofSim = 0	
coef.er.df = c()
repeat {
NofSim = NofSim+1

resd1.value = rnorm(nrow(tmp), mean = resd1.mean, sd = resd1.sd )
tmp$er1 = tmp$er + resd1.value

lm11 = lm(er1~airtemp+prep, data = tmp,na.action=na.exclude)
coef.er.df = rbind(coef.er.df, coef(lm11))

if(NofSim > N){break}
}

coef.er[[i]] <- coef.er.df
## mento carlo stop here ##


tmp = dat.df.ano.rel2[which(toupper(dat.df.ano.rel2$site_id) == toupper(dat.df.annual.mean$site_id[i])),]
lm2 = lm(er~airtemp+prep, data = tmp)

pred1 = predict(lm2, newdata = dat.new1)
pred2 = predict(lm2, newdata = dat.new2)

reg.rel.mod.pred.er[[i]] <- data.frame(airtemp = tmp$airtemp, predt = pred1, prep = tmp$prep, predp = pred2)


reg.abs.er <- rbind(reg.abs.er, c(summary(lm1)$coefficients[,1], 
						  summary(lm1)$coefficients[,4],
			              summary(lm1)$r.squared))
reg.rel.er <- rbind(reg.rel.er, c(summary(lm2)$coefficients[,1], 
						  summary(lm2)$coefficients[,4],
			              summary(lm2)$r.squared))

reg.abs.mod.er[[i]] <- lm1						  
reg.rel.mod.er[[i]] <- lm2							  
}

##############################################################################################################
#plot partial plot [hold airtemp or prep as constant, plot nee response]
flux.info3 = dat.df.annual.mean[,c("site_id", "igbp")]
flux.info3$IGBP2 = factor(flux.info3$igbp)
levels(flux.info3$IGBP2) <- c("CRO", "DBF","ENF","GRA","MF", "SHB","SHB")

####################################################################################################
# plot biome level coefficient							
reg.abs.er = data.frame(reg.abs.er)

#relative influence
#reg.abs = data.frame(reg.rel)

colnames(reg.abs.er) <- c("Int.coef","airtemp.coef", "prep.coef","Int.p","airtemp.p", "prep.p","r2")	
				
reg.abs.er = data.frame(reg.abs.er, dat.df.annual.mean)						
							
reg.abs.er$IGBP2 = factor(reg.abs.er$igbp)

levels(reg.abs.er$IGBP2) <- c("CRO", "DBF","ENF","GRA","MF", "SHB","SHB")
							
boxplot(airtemp.coef~IGBP2, data = reg.abs.er)
boxplot(prep.coef~IGBP2, data = reg.abs.er)

plot(0, xlim = c(0, 1500), ylim = c(-1,1),bty='n',pch='',xlab='Precipation (mm)',ylab=' Coefficient', cex.axis = 2, cex.lab = 2)
box()
points(reg.abs.gpp$prep, reg.abs.gpp$prep.coef, col = "green", cex = 2)
points(reg.abs.er$prep, reg.abs.er$prep.coef, col = "red", cex = 2)

reg.abs.er$dif = reg.abs.er$airtemp.coef + reg.abs.er$prep.coef
boxplot(dif~IGBP2, data = reg.abs.er, subset = (IGBP2 != "Mixed Forest" | IGBP2 != "Wetland"))

###########################################################################################

# step 5: calculate the spatial correlationship

##########################      LEAVE ONE OUT METHOD   ########################
reg.abs.mod.pred1.er.list = list()
delt.er = c()
for (i in 1:nrow(dat.df.annual.mean)){

dat.df.annual.mean1 = dat.df.annual.mean[-i,]

lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean1, na.action=na.exclude)
lm1.er = stepAIC(lm1.er)

#construct new data for prediction

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

pred1 = predict(lm1.er, newdata = dat.new1)
pred2 = predict(lm1.er, newdata = dat.new2)

reg.abs.mod.pred1.er <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)
reg.abs.mod.pred1.er.list[[i]] <- reg.abs.mod.pred1.er

delt.er1 = diff(reg.abs.mod.pred1.er$predp)
delt.er = cbind(delt.er, delt.er1)

}

# LEAVE-ONE-OUT METHODS end here

##########################      MENTO CARLO METHOD   ########################
## mento carlo start from here ##
reg.abs.mod.pred1.er.list = list()
delt.er = c()

# fit model
lm1.er = lm(er~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.er = stepAIC(lm1.er)

# get residual
resd1 = lm1.er$residuals
resd1.mean = mean(resd1)
resd1.sd = sd(resd1)

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

NofSim = 0	
N = 100
repeat {
NofSim = NofSim+1

resd1.value = rnorm(nrow(dat.df.annual.mean), mean = resd1.mean, sd = resd1.sd )
dat.df.annual.mean$er1 = dat.df.annual.mean$er + resd1.value

lm1.er1 = lm(er1~poly(airtemp,2)+poly(prep,2), data = dat.df.annual.mean, na.action=na.exclude)
lm1.er1 = stepAIC(lm1.er1)

pred1 = predict(lm1.er1, newdata = dat.new1)
pred2 = predict(lm1.er1, newdata = dat.new2)

reg.abs.mod.pred1.er <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)
reg.abs.mod.pred1.er.list[[NofSim]] <- reg.abs.mod.pred1.er

delt.er1 = diff(reg.abs.mod.pred1.er$predp)
delt.er = cbind(delt.er, delt.er1)


if(NofSim > N){break}
}

## mento carlo end here ##

######################################################################################

###########################  plot the results in one figure  ####################
# define functions to calculate confidence intervals
cifun <- function(data, ALPHA = 0.05){
  c(mean(data) - qnorm(1-ALPHA/2) * sd(data)/sqrt(length(data)),
    mean(data),
    mean(data) + qnorm(1-ALPHA/2) * sd(data)/sqrt(length(data)))
}

cifun(pred.er[1,], 0.1) 

# define functions to calculate mean +- 1 S.E.
cifun2 <- function(data){
  c(mean(data) - sd(data),
    mean(data),
    mean(data) + sd(data))
}

## plot gpp/er sensitive to precipitaion 
pred.er = c()
pred.gpp = c()

for (i in 1:length(reg.abs.mod.pred1.er.list)){
pred.er = cbind(pred.er, reg.abs.mod.pred1.er.list[[i]][,4])
pred.gpp = cbind(pred.gpp, reg.abs.mod.pred1.gpp.list[[i]][,4])
}

## plot using the confidence intervals

pred.er2 = apply(pred.er, 1, cifun)
pred.gpp2 = apply(pred.gpp, 1, cifun)

delt.gpp2 = apply(delt.gpp, 1, cifun)
delt.er2 = apply(delt.er, 1, cifun)

prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)

# dat.df.annual.mean$igbp = CRO GRA DBF ENF DBF DBF CRO ENF DBF MF  DBF WSA MF  GRA ENF OSH GRA
# DBF/MF:1, ENF:2, others:3
pch =   c(3,  3,  1, 2,   1, 1,  3,  2,  1,  1,   1,  3, 1,  3,  2,  3,3)

png("D:/zhihua/dataset/results/er.prep01032017_2_ci.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,5,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2000),bty='n',pch='',xlab='Precipation (mm)',ylab=' GPP or ER (g C m-2 yr-1)', cex.axis = 2, cex.lab = 2)

### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.gpp2[3,], lty = 'dashed', col = 'green')
lines(prep.grd,pred.gpp2[1,], lty = 'dashed', col = 'green')
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")

arrows(dat.df.annual.mean$prep,dat.df.annual.mean$gpp-dat.df.annual.sd$gpp,
        dat.df.annual.mean$prep,dat.df.annual.mean$gpp+dat.df.annual.sd$gpp, code=3, length=0.05, angle = 90, lwd = 0.5, col = "green")
arrows(dat.df.annual.mean$prep-dat.df.annual.sd$prep,dat.df.annual.mean$gpp,
        dat.df.annual.mean$prep+dat.df.annual.sd$prep,dat.df.annual.mean$gpp, code=3, length=0.05, angle = 90, lwd = 0.5, col = "green")

### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.er2[3,], lty = 'dashed', col = 'red')
lines(prep.grd,pred.er2[1,], lty = 'dashed', col = 'red')
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")

arrows(dat.df.annual.mean$prep,dat.df.annual.mean$er-dat.df.annual.sd$er,
        dat.df.annual.mean$prep,dat.df.annual.mean$er+dat.df.annual.sd$er, code=3, length=0.05, angle = 90, lwd = 0.5, col = "red")
arrows(dat.df.annual.mean$prep-dat.df.annual.sd$prep,dat.df.annual.mean$er,
        dat.df.annual.mean$prep+dat.df.annual.sd$prep,dat.df.annual.mean$er, code=3, length=0.05, angle = 90, lwd = 0.5, col = "red")


#plot spatial sensitity of gpp vs er in response to precipitation
library(oce)
plotInset(1000, 50, 1500, 1000,
          expr= {
		  
# plot(0, xlim = c(0, 1500), ylim = c(-20,150),bty='n',pch='',ylab='',xlab='',yaxt='n', ann=FALSE)
# change to per 100 mm
plot(0, xlim = c(0, 1500), ylim = c(-20,150)*2, pch='', xlab = "Precipation (mm)", 
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

#read into soil resprition sensitivity to precipitaion
# delt.rs.srdb = read.csv("C:/zhihua/dataset/srdb_v3-2/delt.rs.srdb.csv")[,-1]
# lines(smooth.spline(reg.abs.mod.pred1.er[-1,3], delt.rs.srdb*2), lty = 2, lwd = 6, col = "blue")

abline(h = 0)
abline(v = reg.abs.mod.pred1.gpp[-1,3][which.max(which(delt.gpp > delt.er))])

  },
		  mar=c(5.5,0,0,0))

#plot temporal sensitity of gpp vs er in response to precipitation
################ plot GPP 	  ####################################

library(oce)
plotInset(50, 1250, 500, 2000,
          expr= {
		  
mar=c(5,5,5,5)+.1
		  
plot(0, xlim = c(-2, 2)*50*2, ylim = c(0,7),bty='n',pch='',ylab='',
   # xlab='Temproal Sensitivity (Per 50mm)',
   xlab = expression(paste(beta ["temporal"])),
   yaxt='n', 
   # ann=FALSE
   cex.lab = 1.5, cex.axis = 1.5)
#draw a box around the plot
box()
abline(v = 0)

# Deciduous Forest
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "DBF" | reg.abs.gpp$IGBP2 == "MF")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 2, x1 = x1[4], y1 = 2, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 2, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 1.5, paste("DBF/MF ( n = ", length(which(reg.abs.gpp$IGBP2 == "DBF"| reg.abs.gpp$IGBP2 == "MF")), " )", sep = ""),
  cex = 1)
  

# Evergreen Forest
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "ENF")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, paste("ENF ( n = ", length(which(reg.abs.gpp$IGBP2 == "ENF")), " )", sep = ""),
  cex = 1)

# SHB/GRA/CRO
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 6, x1 = x1[4], y1 = 6, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 6, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 5.5, paste("SHB/GRA/CRO ( n = ", length(which(reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")), " )", sep = ""),
  cex = 1)

################## add ER         #################################
# Deciduous Forest
x11 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "DBF" | reg.abs.er$IGBP2 == "MF")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 1, x1 = x1[4], y1 = 1, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 1, type = "p", pch = 17, cex = 2, col = "red")

# Evergreen Forest
x11 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "ENF")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 3, x1 = x1[4], y1 = 3, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 3, type = "p", pch = 17, cex = 2, col = "red")

# SHB/GRA/CRO
x11 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "SHB" | reg.abs.er$IGBP2 == "GRA" | reg.abs.er$IGBP2 == "CRO")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 5, x1 = x1[4], y1 = 5, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 5, type = "p", pch = 17, cex = 2, col = "red")

  },
		  mar=c(5,1,0,0))

############## devide the site based on annual rainfall		  
	
library(oce)
plotInset(500, 1250, 1000, 2000,
          expr= {
		  
		  
plot(0, xlim = c(-2, 2)*50*2, ylim = c(0,7),bty='n',pch='',ylab='',
   # xlab='Temproal Sensitivity (Per 50mm)',
   xlab = expression(paste(beta ["temporal"])),
   yaxt='n', 
   # ann=FALSE
   cex.lab = 1.5, cex.axis = 1.5)
   
#draw a box around the plot
box()
abline(v = 0)

# DBF/MF
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "DBF" | reg.abs.gpp$IGBP2 == "MF")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 2, x1 = x1[4], y1 = 2, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 2, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 1.5, paste("DBF/MF ( n = ", length(which(reg.abs.gpp$IGBP2 == "DBF"| reg.abs.gpp$IGBP2 == "MF")), " )", sep = ""),
  cex = 1)
  
# ENF/SHB/GRA/CRO
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "ENF" | reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, paste("ENF/SHB/GRA/CRO ( n = ", length(which(reg.abs.gpp$IGBP2 == "ENF" |reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")), " )", sep = ""),
  cex = 1)

################## add ER         #################################
# DBF/MF
x11 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "DBF" | reg.abs.gpp$IGBP2 == "MF")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 1, x1 = x1[4], y1 = 1, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 1, type = "p", pch = 17, cex = 2, col = "red")

# ENF/SHB/GRA/CRO
x11 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "ENF" | reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 3, x1 = x1[4], y1 = 3, lty = 1, lwd = 2, col = "red")
points(x = mean(x11), y = 3, type = "p", pch = 17, cex = 2, col = "red")

  },
		  mar=c(5,1,0,0))
	  
		  
dev.off()

### DOES NOT WORK FROM HERE
## plot temporal sensitivity from Mento Carlo results
length(coef.gpp) 
length(coef.er) 
# each element represent one flux site, and contains 100 replicates
prep1.id = which(dat.df.annual.mean$prep < 750) 
coef.er1 = c()
coef.gpp1 = c()
for (i in prep1.id){
coef.er1 = c(coef.er1, coef.er[[i]][,3])
coef.gpp1 = c(coef.gpp1, coef.gpp[[i]][,3])
}

mean(coef.er1)
sd(coef.er1)

mean(coef.gpp1)
sd(coef.gpp1)

coef.prep1 = data.frame(coef1 = 100*c(coef.gpp1, coef.er1), var1 = c(rep("GPP", length(coef.gpp1)),rep("ER", length(coef.er1))))

library(ggplot2)

ggplot(coef.prep1, aes(x=coef1, fill=var1)) +
    geom_histogram(binwidth=.5, alpha=.5, position="identity")

ggplot(coef.prep1, aes(x=coef1, fill=var1)) + geom_density(alpha=.3)

ggplot(coef.prep1, aes(x=coef1, fill=var1)) + geom_density()


mean(reg.abs.gpp$prep.coef[prep1.id])
mean(reg.abs.er$prep.coef[prep1.id])

### DOES NOT WORK STOP HERE






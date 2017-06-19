## final script to produce figures in foler results2

## use bootstrapping methods

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

# plot final points
flux.info3 = merge(fluxsite.df, dat.df.annual.mean, by.x = "SITE_ID", by.y = "site_id")
flux.info.sp3 = flux.info3			 
coordinates(flux.info.sp3) <- ~LOCATION_LONG + LOCATION_LAT
projection(flux.info.sp3) <- "+proj=longlat +datum=WGS84"

plot(na.state[which(na.state$ISO == "USA" & na.state$NAME_1 != "Alaska"),])
plot(flux.info.sp3, add = T)
text(x = flux.info3$LOCATION_LONG, y = flux.info3$LOCATION_LAT, flux.info3$SITE_ID)
#############################################################################################################

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

# random sampling to calculate temporal sensitivity
# dat is the data frame, ln: minimum to calculate, n: number of sampling

boot.ec1.gpp = function(dat, n = 100){  
# sampling at least 5  of the data

if (nrow(dat) <= 6) {
ln = nrow(dat) - 1
} else if (nrow(dat) >= 7 & nrow(dat) <= 10) {
ln = nrow(dat) - 2
} 
else {
ln = nrow(dat) - 3
}

Coef.df = c()
NofSim = 0	
repeat{
NofSim = NofSim+1
n1 = sample(ln:nrow(dat),1)
idx = sample(1:nrow(dat),n1)
# print(idx)
dat1 = dat[idx,]

lm11 = lm(gpp~airtemp+prep, data = dat1,na.action=na.exclude)

Coef.df = rbind(Coef.df, c(coef(lm11),summary(lm11)$r.squared))
if(NofSim > n){break}
}

Coef.df = data.frame(Coef.df)
colnames(Coef.df) <- c("int","airtemp","prep","r2")
return(Coef.df)

}								

boot.ec1.gpp(dat = tmp, n = 100)

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

## bootstraping from here ##
coef.gpp.df = boot.ec1.gpp(dat = tmp, n = 100)
coef.gpp[[i]] <- coef.gpp.df

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

# use bootstrapping methods here

boot.ec2.gpp = function(dat = dat.df.annual.mean, ln = 13, n = 100){  
require(MASS)
# sampling at least 5  of the data

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

reg.abs.mod.pred1.gpp.list = list()
reg.abs.mod.pred1.gpp.model.coef = c() #store model coefficients
delt.gpp = c()
rmse.gpp = c() #store model efficients, including rmse, mae, aic, r-squired

NofSim = 0	

repeat{

NofSim = NofSim+1
n1 = sample(ln:nrow(dat),1)
idx = sample(1:nrow(dat),n1)
# print(idx)
dat1 = dat.df.annual.mean[idx,]

lm1 = lm(gpp~poly(airtemp,2)+poly(prep,2), data = dat1, na.action=na.exclude)
lm1 = stepAIC(lm1)

rmse <- round(sqrt(mean(resid(lm1)^2)), 2)
mae <- round(mean(abs(resid(lm1))), 2)
aic = AIC(lm1)
r2 = round(summary(lm1)$r.squared, 2)

pred1 = predict(lm1, newdata = dat.new1)
pred2 = predict(lm1, newdata = dat.new2)

reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

delt.gpp1 = diff(reg.abs.mod.pred1.gpp$predp)

rmse.gpp = rbind(rmse.gpp, c(rmse, mae,aic,r2))
reg.abs.mod.pred1.gpp.list[[NofSim]] <- reg.abs.mod.pred1.gpp
delt.gpp = cbind(delt.gpp, delt.gpp1)
reg.abs.mod.pred1.gpp.model.coef = rbind(reg.abs.mod.pred1.gpp.model.coef, coef(lm1))

if(NofSim > n){break}
}


return(list(reg.abs.mod.pred1.gpp.list, reg.abs.mod.pred1.gpp.model.coef, delt.gpp, rmse.gpp))

}								

list.gpp = boot.ec2.gpp(dat = dat.df.annual.mean, ln = 13, n = 100)
reg.abs.mod.pred1.gpp.list = list.gpp[[1]]
reg.abs.mod.pred1.gpp.model.coef = list.gpp[[2]]
delt.gpp = list.gpp[[3]]
rmse.gpp = list.gpp[[4]]

###############################################################################################

#################################################################################
## repeat step 4 and 5

boot.ec1.er = function(dat, n = 100){  
# sampling at least 5  of the data

if (nrow(dat) <= 6) {
ln = nrow(dat) - 1
} else if (nrow(dat) >= 7 & nrow(dat) <= 10) {
ln = nrow(dat) - 2
} 
else {
ln = nrow(dat) - 3
}

Coef.df = c()
NofSim = 0	
repeat{
NofSim = NofSim+1
n1 = sample(ln:nrow(dat),1)
idx = sample(1:nrow(dat),n1)
# print(idx)
dat1 = dat[idx,]

lm11 = lm(er~airtemp+prep, data = dat1,na.action=na.exclude)

Coef.df = rbind(Coef.df, c(coef(lm11),summary(lm11)$r.squared))
if(NofSim > n){break}
}

Coef.df = data.frame(Coef.df)
colnames(Coef.df) <- c("int","airtemp","prep","r2")
return(Coef.df)

}								

boot.ec1.er(dat = tmp, n = 100)


## step 4: calculate the ER and Rainfall/temperature relationship for each individual site
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

## bootstraping from here ##
coef.er.df = boot.ec1.er(dat = tmp, n = 100)
coef.er[[i]] <- coef.er.df

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

boot.ec2.er = function(dat = dat.df.annual.mean, ln = 13, n = 100){  
require(MASS)
# sampling at least 5  of the data

dat.new1 = data.frame(airtemp = seq(1,25, length.out = 25), prep = mean(dat.df.annual.mean$prep, na.rm = TRUE))
dat.new2 = data.frame(airtemp = mean(dat.df.annual.mean$airtemp, na.rm = TRUE), prep = seq(100,1300, length.out = 25))

reg.abs.mod.pred1.gpp.list = list()
reg.abs.mod.pred1.gpp.model.coef = c() #store model coefficients
delt.gpp = c()
rmse.gpp = c() #store model efficients, including rmse, mae, aic, r-squired

NofSim = 0	

repeat{

NofSim = NofSim+1
n1 = sample(ln:nrow(dat),1)
idx = sample(1:nrow(dat),n1)
# print(idx)
dat1 = dat.df.annual.mean[idx,]

lm1 = lm(er~poly(airtemp,2)+poly(prep,2), data = dat1, na.action=na.exclude)
lm1 = stepAIC(lm1)

rmse <- round(sqrt(mean(resid(lm1)^2)), 2)
mae <- round(mean(abs(resid(lm1))), 2)
aic = AIC(lm1)
r2 = round(summary(lm1)$r.squared, 2)

pred1 = predict(lm1, newdata = dat.new1)
pred2 = predict(lm1, newdata = dat.new2)

reg.abs.mod.pred1.gpp <- data.frame(airtemp = seq(1,25, length.out = 25), predt = pred1, prep = seq(100,1300, length.out = 25), predp = pred2)

delt.gpp1 = diff(reg.abs.mod.pred1.gpp$predp)

rmse.gpp = rbind(rmse.gpp, c(rmse, mae,aic,r2))
reg.abs.mod.pred1.gpp.list[[NofSim]] <- reg.abs.mod.pred1.gpp
delt.gpp = cbind(delt.gpp, delt.gpp1)
reg.abs.mod.pred1.gpp.model.coef = rbind(reg.abs.mod.pred1.gpp.model.coef, coef(lm1))

if(NofSim > n){break}
}


return(list(reg.abs.mod.pred1.gpp.list, reg.abs.mod.pred1.gpp.model.coef, delt.gpp, rmse.gpp))

}								


list.er = boot.ec2.er(dat = dat.df.annual.mean, ln = 13, n = 100)
reg.abs.mod.pred1.er.list = list.er[[1]]
reg.abs.mod.pred1.er.model.coef = list.er[[2]]
delt.er = list.er[[3]]
rmse.er = list.er[[4]]

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


# before calculating, replace value > 90 percentile and < 10 percentile with NA
Rplc = function(x){
p1 = quantile(x, probs = c(0.1,0.9))
x[which(x < p1[1] | x > p1[2])] = NA
return(x)
}

####### USE ggplot 2 to map 
prep.grd = prep = seq(100,1300, length.out = 25)
airtemp.grd = seq(1,25, length.out = 25)

# calculate the plot spatial coefficient
delt.gpp_1 = apply(delt.gpp, 1, Rplc)
delt.gpp_2 = as.vector(delt.gpp_1)

delt.er_1 = apply(delt.er, 1, Rplc)
delt.er_2 = as.vector(delt.er_1)

delt.gpp_3 = data.frame(coef1=delt.gpp_2, prep = rep(prep.grd[-1], each = 101), flux = "GPP")		
delt.er_3 = data.frame(coef1=delt.er_2, prep = rep(prep.grd[-1], each = 101), flux = "ER")		
d2 = rbind(delt.gpp_3,delt.er_3)

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

r2s = data.frame(r = c(rmse.er[,4],rmse.gpp[,4]), flux = c(rep("ER", 101),rep("GPP", 101)))

library(plyr)
cdat <- ddply(r2s, "flux", summarise, rating.mean=mean(r))
ggplot(r2s, aes(x=r, fill=flux)) +
    geom_histogram(binwidth=.02, alpha=.5, position="identity") +
    geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=flux),
               linetype="dashed", size=1)

	
# calculate and plot temporal coeffifiecnt
# for annual temperature and precipitations

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
					
pred.er.t2 = data.frame(coef1=pred.er.t, prep = rep(dat.df.annual.mean[,c("prep")], each = 101), flux = "ER")					
pred.gpp.t2 = data.frame(coef1=pred.gpp.t, prep = rep(dat.df.annual.mean[,c("prep")], each = 101), flux = "GPP")					

d1 = rbind(pred.er.t2, pred.gpp.t2)	
d1 = d1[complete.cases(d1),]
d1$coef1 = d1$coef1*100
	
# reclassifiy prep <200, 200-400, 400-600,600-800,800-1000,1000-1200,1200-1400,1400-1600
# d1$prep[which(d1$prep < 200)] = mean(d1$prep[which(d1$prep < 200)])
# d1$prep[which( d1$prep > 200 & d1$prep < 400)] = mean(d1$prep[which( d1$prep > 200 & d1$prep < 400)])
# d1$prep[which( d1$prep > 400 & d1$prep < 600)] = mean(d1$prep[which( d1$prep > 400 & d1$prep < 600)])
# d1$prep[which( d1$prep > 600 & d1$prep < 800)] = mean(d1$prep[which( d1$prep > 600 & d1$prep < 800)])
# d1$prep[which( d1$prep > 800 & d1$prep < 1000)] = mean(d1$prep[which( d1$prep > 800 & d1$prep < 1000)])
# d1$prep[which( d1$prep > 1000 & d1$prep < 1200)] = mean(d1$prep[which( d1$prep > 1000 & d1$prep < 1200)])
	
# d1$coef1[which((d1$prep < 200) & (d1$flux == "GPP"))] = mean(d1$coef1[which((d1$prep < 200) & (d1$flux == "GPP"))])
# d1$coef1[which( (d1$prep > 200 & d1$prep < 400)& (d1$flux == "GPP"))] = mean(d1$coef1[which( (d1$prep > 200 & d1$prep < 400)& (d1$flux == "GPP"))])
# d1$coef1[which( (d1$prep > 400 & d1$prep < 600)& (d1$flux == "GPP"))] = mean(d1$coef1[which( (d1$prep > 400 & d1$prep < 600)& (d1$flux == "GPP"))])
# d1$coef1[which( (d1$prep > 600 & d1$prep < 800)& (d1$flux == "GPP"))] = mean(d1$coef1[which( (d1$prep > 600 & d1$prep < 800)& (d1$flux == "GPP"))])
# d1$coef1[which( (d1$prep > 800 & d1$prep < 1000)& (d1$flux == "GPP"))] = mean(d1$coef1[which( (d1$prep > 800 & d1$prep < 1000)& (d1$flux == "GPP"))])
# d1$coef1[which( (d1$prep > 1000 & d1$prep < 1200)& (d1$flux == "GPP"))] = mean(d1$coef1[which( (d1$prep > 1000 & d1$prep < 1200)& (d1$flux == "GPP"))])

# d1$coef1[which((d1$prep < 200) & (d1$flux == "ER"))] = mean(d1$coef1[which((d1$prep < 200) & (d1$flux == "ER"))])
# d1$coef1[which( (d1$prep > 200 & d1$prep < 400)& (d1$flux == "ER"))] = mean(d1$coef1[which( (d1$prep > 200 & d1$prep < 400)& (d1$flux == "ER"))])
# d1$coef1[which( (d1$prep > 400 & d1$prep < 600)& (d1$flux == "ER"))] = mean(d1$coef1[which( (d1$prep > 400 & d1$prep < 600)& (d1$flux == "ER"))])
# d1$coef1[which( (d1$prep > 600 & d1$prep < 800)& (d1$flux == "ER"))] = mean(d1$coef1[which( (d1$prep > 600 & d1$prep < 800)& (d1$flux == "ER"))])
# d1$coef1[which( (d1$prep > 800 & d1$prep < 1000)& (d1$flux == "ER"))] = mean(d1$coef1[which( (d1$prep > 800 & d1$prep < 1000)& (d1$flux == "ER"))])
# d1$coef1[which( (d1$prep > 1000 & d1$prep < 1200)& (d1$flux == "ER"))] = mean(d1$coef1[which( (d1$prep > 1000 & d1$prep < 1200)& (d1$flux == "ER"))])
		
	
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

# multiplot(p1, p2, cols=1)					
print(p1)
ggsave("F:/zhihua/dataset/results2/fig2.ec1.png", width = 4, height = 3, units = "in")
print(p2)
ggsave("F:/zhihua/dataset/results2/fig2.ec2.png", width = 4, height = 3, units = "in")
					
# plot spatial sensitivity difference
delt.temporal = data.frame(prep = dat.df.annual.mean[,c("prep")], mean = 100*apply(delt.temporal, 1, mean, na.rm = TRUE), 
			   sd = 100*apply(delt.temporal, 1, sd, na.rm = TRUE)) 



r1 = c()
r2 = c()

for (i in 1:length(coef.er)){
r1 = c(r1, coef.er[[i]][,4])
r2 = c(r2, coef.gpp[[i]][,4])
}

r2t = data.frame(r = c(r1,r2), flux = c(rep("ER", length(r1)),rep("GPP", length(r2))))

library(plyr)
cdat <- ddply(r2t, "flux", summarise, rating.mean=mean(r))

ggplot(r2t, aes(x=r, colour=flux)) +
    geom_density() +
    geom_vline(data=cdat, aes(xintercept=rating.mean,  colour=flux),
               linetype="dashed", size=1)

			   
d3 = rbind(data.frame(d1, doman = "Temporal"),data.frame(d2, doman = "Spatial"))			   
ggplot(d3, aes(x=prep, y=coef1, color=flux)) + 
	geom_point(shape=1, cex = 3) +
	facet_grid(.~doman)
    scale_colour_hue(l=50) + # Use a slightly darker palette than normal
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,    # Don't add shaded confidence region
                fullrange=FALSE) # Extend regression lines



write.csv(d3, file = "F:/zhihua/dataset/results2/EC.sensitivity2.csv")
write.csv(delt.temporal, file = "F:/zhihua/dataset/results2/delt.temporal.ec.csv")
write.csv(delt.spatial, file = "F:/zhihua/dataset/results2/delt.spatial.ec.csv")

####### END of USE ggplot 2 to map 					

	
##############################################################################################################
## old mapping 

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

# for two regions
w.idx = which(prep.grd < 750)
e.idx = which(prep.grd > 750)

delt = data.frame(mean = c(mean(delt.gpp[e.idx-1,]),mean(delt.gpp[w.idx,]),mean(delt.er[e.idx-1,]),mean(delt.er[w.idx,])),
					se = c(sd(delt.gpp[e.idx-1,]),sd(delt.gpp[w.idx,]),sd(delt.er[e.idx-1,]),sd(delt.er[w.idx,])),
				    flux = c("GPP","GPP","TER","TER"),
				    region = c("> 750 mm", "< 750 mm","> 750 mm","< 750 mm"),
					type = rep("spatial",4))
	
# dat.df.annual.mean$igbp = CRO GRA DBF ENF DBF DBF CRO ENF DBF MF  DBF WSA MF  GRA ENF OSH GRA
# DBF/MF:1, ENF:2, others:3
pch =   c(3,  3,  1, 2,   1, 1,  3,  2,  1,  1,   1,  3, 1,  3,  2,  3,3)

png("F:/zhihua/dataset/results2/ec.sensitivity_ci-2.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,6,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2000),bty='n',pch='',
        xlab='Precipation (mm)',
		ylab=expression("GPP/TER" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")),
		cex.axis = 2, cex.lab = 2)

text(x = 1450, y = 1900, "a)",   cex = 3)
		
### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "blue")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), 
   col = rgb(0,0,1,0.3), 
   border = NA)
# intervals
# lines(prep.grd,pred.gpp2[3,], lty = 'dashed', col = 'blue')
# lines(prep.grd,pred.gpp2[1,], lty = 'dashed', col = 'blue')
# lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "blue")

# arrows(dat.df.annual.mean$prep,dat.df.annual.mean$gpp-dat.df.annual.sd$gpp,
        # dat.df.annual.mean$prep,dat.df.annual.mean$gpp+dat.df.annual.sd$gpp, code=3, length=0.05, angle = 90, lwd = 0.5, col = "green")
# arrows(dat.df.annual.mean$prep-dat.df.annual.sd$prep,dat.df.annual.mean$gpp,
        # dat.df.annual.mean$prep+dat.df.annual.sd$prep,dat.df.annual.mean$gpp, code=3, length=0.05, angle = 90, lwd = 0.5, col = "green")

### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), 
  col = rgb(1,0,0,0.3), 
  border = NA)
# intervals
# lines(prep.grd,pred.er2[3,], lty = 'dashed', col = 'red')
# lines(prep.grd,pred.er2[1,], lty = 'dashed', col = 'red')
# lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")

# arrows(dat.df.annual.mean$prep,dat.df.annual.mean$er-dat.df.annual.sd$er,
        # dat.df.annual.mean$prep,dat.df.annual.mean$er+dat.df.annual.sd$er, code=3, length=0.05, angle = 90, lwd = 0.5, col = "red")
# arrows(dat.df.annual.mean$prep-dat.df.annual.sd$prep,dat.df.annual.mean$er,
        # dat.df.annual.mean$prep+dat.df.annual.sd$prep,dat.df.annual.mean$er, code=3, length=0.05, angle = 90, lwd = 0.5, col = "red")


#plot spatial sensitity of gpp vs er in response to precipitation
library(oce)
plotInset(1000, 50, 1500, 1000,
          expr= {
		  
# change to per 100 mm
plot(0, xlim = c(0, 1500), ylim = c(0,150)*2, pch='', 
     xlab = "Precipation (mm)", 
	 ylab = expression(paste(italic(beta) ["spatial"])),
	 cex.lab = 1.5, cex.axis = 1.5)

text(x = 20, y = 290, "1)",   cex = 1.5) 
	 
#draw a box around the plot
box()
#######  GPP
lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "blue")
# add fill
polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
        c(rev(delt.gpp2[3,]*2), delt.gpp2[1,]*2), 
		# col = 'grey90', 
		col = rgb(0,0,1,0.3), 
		border = NA)
lines(prep.grd[-1], delt.gpp2[2,]*2, lty = 1, lwd = 3, col = "blue")
# intervals
# lines(prep.grd[-1],delt.gpp2[3,]*2, lty = 'dashed', col = 'blue')
# lines(prep.grd[-1],delt.gpp2[1,]*2, lty = 'dashed', col = 'blue')

#######  ER
lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# add fill
 polygon(c(rev(prep.grd[-1]), prep.grd[-1]), 
         c(rev(delt.er2[3,]*2), delt.er2[1,]*2), 
		 # col = 'grey90', 
		 col = rgb(1,0,0,0.3), 
		border = NA)
 lines(prep.grd[-1], delt.er2[2,]*2, lty = 1, lwd = 3, col = "red")
# intervals
# lines(prep.grd[-1],delt.er2[3,]*2, lty = 'dashed', col = 'red')
# lines(prep.grd[-1],delt.er2[1,]*2, lty = 'dashed', col = 'red')

#read into soil resprition sensitivity to precipitaion
# delt.rs.srdb = read.csv("C:/zhihua/dataset/srdb_v3-2/delt.rs.srdb.csv")[,-1]
# lines(smooth.spline(reg.abs.mod.pred1.er[-1,3], delt.rs.srdb*2), lty = 2, lwd = 6, col = "blue")

# abline(h = 0)
# add threshold regions

x1 = c(650,800)
y1 = c(500,500)
y2 = c(-100,-100)
# abline(v = c(650, 750, 700 + 100), lty = c(2,1,2))
abline(v = c(650, 700 + 100), lty = c(2,2))


polygon(c(x1, rev(x1)), c(y1, rev(y2)),
      col = rgb(0.5,0.5,0.5,0.3), 
	  border = NA)

  },
		  mar=c(5.5,0,0,0))

#plot temporal sensitity of gpp vs er in response to precipitation

############## devide the site based on annual rainfall		  
	
# use bar plot
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "DBF" | reg.abs.gpp$IGBP2 == "MF")]*50*2
x12 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "ENF" | reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")]*50*2
x13 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "ENF" | reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")]*50*2
x14 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "DBF" | reg.abs.gpp$IGBP2 == "MF")]*50*2

se <- function(x) sqrt(var(x)/length(x))

D1 = data.frame(mean = c(mean(x11),mean(x12),mean(x13),mean(x14)),
				se = c(se(x11),se(x12),se(x13),se(x14)),
				flux = c("GPP","GPP","TER","TER"),
				region = c("> 750 mm", "< 750 mm","> 750 mm","< 750 mm"),
				type = rep("temproal",4))
	

# write.csv(D1, "F:/zhihua/dataset/results2/D1.csv")
	
ylab = 	expression(paste(beta ["temporal"]) ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "/100mm"))
color1 = c("blue", "red")

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
	
	
library(oce)
library(ggplot2)
require(grid)
print(p1, vp=viewport(.3, .8, .3, .3)) #viewport, first two is the x,y coor; last 2 is the inset size
   
dev.off()


png("F:/zhihua/dataset/results2/ec.sensitivity_ci.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,6,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2000),bty='n',pch='',
        xlab='Precipation (mm)',
		ylab=expression("GPP/TER" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")),
		cex.axis = 2, cex.lab = 2)

text(x = 1450, y = 1900, "a)",   cex = 3)
		
### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.gpp2[3,], lty = 'dashed', col = 'green')
lines(prep.grd,pred.gpp2[1,], lty = 'dashed', col = 'green')
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "green")

# arrows(dat.df.annual.mean$prep,dat.df.annual.mean$gpp-dat.df.annual.sd$gpp,
        # dat.df.annual.mean$prep,dat.df.annual.mean$gpp+dat.df.annual.sd$gpp, code=3, length=0.05, angle = 90, lwd = 0.5, col = "green")
# arrows(dat.df.annual.mean$prep-dat.df.annual.sd$prep,dat.df.annual.mean$gpp,
        # dat.df.annual.mean$prep+dat.df.annual.sd$prep,dat.df.annual.mean$gpp, code=3, length=0.05, angle = 90, lwd = 0.5, col = "green")

### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), col = 'grey80', border = NA)
# intervals
lines(prep.grd,pred.er2[3,], lty = 'dashed', col = 'red')
lines(prep.grd,pred.er2[1,], lty = 'dashed', col = 'red')
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")

# arrows(dat.df.annual.mean$prep,dat.df.annual.mean$er-dat.df.annual.sd$er,
        # dat.df.annual.mean$prep,dat.df.annual.mean$er+dat.df.annual.sd$er, code=3, length=0.05, angle = 90, lwd = 0.5, col = "red")
# arrows(dat.df.annual.mean$prep-dat.df.annual.sd$prep,dat.df.annual.mean$er,
        # dat.df.annual.mean$prep+dat.df.annual.sd$prep,dat.df.annual.mean$er, code=3, length=0.05, angle = 90, lwd = 0.5, col = "red")


#plot spatial sensitity of gpp vs er in response to precipitation
library(oce)
plotInset(1000, 50, 1500, 1000,
          expr= {
		  
# change to per 100 mm
plot(0, xlim = c(0, 1500), ylim = c(0,150)*2, pch='', 
     xlab = "Precipation (mm)", 
	 ylab = expression(paste(italic(beta) ["spatial"])),
	 cex.lab = 1.5, cex.axis = 1.5)

text(x = 20, y = 290, "1)",   cex = 1.5) 
	 
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

#read into soil resprition sensitivity to precipitaion
# delt.rs.srdb = read.csv("C:/zhihua/dataset/srdb_v3-2/delt.rs.srdb.csv")[,-1]
# lines(smooth.spline(reg.abs.mod.pred1.er[-1,3], delt.rs.srdb*2), lty = 2, lwd = 6, col = "blue")

# abline(h = 0)
# add threshold regions

abline(v = c(730-100, 730, 730 + 100), lty = c(2,1,2))

  },
		  mar=c(5.5,0,0,0))

#plot temporal sensitity of gpp vs er in response to precipitation
	
library(oce)

plotInset(50, 1250, 500, 2000,
          expr= {
		  
mar=c(5,5,5,5)+.1
		  
plot(0, xlim = c(-2, 2)*50*2, ylim = c(0,5),bty='n',pch='',ylab='',
   xlab = expression(paste(beta ["temporal"])),
   yaxt='n', 
   # ann=FALSE
   cex.lab = 1.5, cex.axis = 1.5)

text(x = -190, y = 4.7, "2)",   cex = 1.5)   
 
#draw a box around the plot
box()
abline(v = 0)

# DBF/MF
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "DBF" | reg.abs.gpp$IGBP2 == "MF")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 2, x1 = x1[4], y1 = 2, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 2, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 1.5, paste("> 750 mm ( n = ", length(which(reg.abs.gpp$IGBP2 == "DBF"| reg.abs.gpp$IGBP2 == "MF")), " )", sep = ""),
  cex = 1)
# text(x = x1[3], y = 1.5, paste("MAP > 750 mm"),   cex = 1)

# ENF/SHB/GRA/CRO
x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "ENF" | reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")]*50*2
x1 = quantile(x11,probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
segments(x0 = x1[2], y0 = 4, x1 = x1[4], y1 = 4, lty = 1, lwd = 2, col = "green")
points(x = mean(x11), y = 4, type = "p", pch = 15, cex = 2, col = "green")
text(x = x1[3], y = 3.5, paste("< 750 mm ( n = ", length(which(reg.abs.gpp$IGBP2 == "ENF" |reg.abs.gpp$IGBP2 == "SHB" | reg.abs.gpp$IGBP2 == "GRA" | reg.abs.gpp$IGBP2 == "CRO")), " )", sep = ""),
  cex = 1)
# text(x = x1[3], y = 3.5, paste("MAP < 750 mm"),   cex = 1)
  
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
		  mar=c(5.5,0,0,0))


dev.off()
	  
## plot individual site

png("F:/zhihua/dataset/results2/gpp.prep.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,6,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2300),bty='n',pch='',
     xlab='Precipation (mm)',
	 ylab=expression("GPP" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")),
     cex.axis = 2, cex.lab = 2)

#draw a box around the plot
box()
#plot Evergreen Forest
for(i in 1:length(flux.info3$IGBP2)){

if (flux.info3$IGBP2[i] == "ENF" ){
lines(reg.abs.mod.pred.gpp[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.gpp[[i]][,4]+dat.df.annual.mean$gpp[i],lwd = 1.5, lty = 1)

} else if (flux.info3$IGBP2[i] == "DBF") {
lines(reg.abs.mod.pred.gpp[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.gpp[[i]][,4]+dat.df.annual.mean$gpp[i],lwd = 1.5,col = 2)

} else if (flux.info3$IGBP2[i] == "MF") {
lines(reg.abs.mod.pred.gpp[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.gpp[[i]][,4]+dat.df.annual.mean$gpp[i],lwd = 1.5,col = 2)

} else if (flux.info3$IGBP2[i] == "SHB") {
lines(reg.abs.mod.pred.gpp[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.gpp[[i]][,4]+dat.df.annual.mean$gpp[i],lwd = 1.5, col=4)

} else if (flux.info3$IGBP2[i] == "GRA") {
lines(reg.abs.mod.pred.gpp[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.gpp[[i]][,4]+dat.df.annual.mean$gpp[i],lwd = 1.5, col=4)

} else if (flux.info3$IGBP2[i] == "CRO") {
lines(reg.abs.mod.pred.gpp[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.gpp[[i]][,4]+dat.df.annual.mean$gpp[i],lwd = 1.5, col=4)

}

# text(x = dat.df.annual.mean$prep[i],
     # y = dat.df.annual.mean$gpp[i],
	 # paste(factor(flux.info3[i,1]), factor(flux.info3[i,2])))


}

legend("topleft", 
	   # inset=0.05, 
	   # legend = c("ENF","DBF","MF","SHB","GRA","CRO"),
	   legend = c("ENF","DBF/MF","Shrub/GRA/CRO"),
	   horiz=F,
	   lwd = 2,
	   col = c(1,2,4),
	   text.col = c(1,2,4),
	   cex = 1.5,
	   box.col = "transparent",
	   bg = "transparent")


### plot gpp
lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3, col = "blue")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.gpp2[3,]), pred.gpp2[1,]), 
     col = rgb(0,0,1,0.3), 
	 border = NA)
# intervals
# lines(prep.grd,pred.gpp2[3,], lty = 'dashed')
# lines(prep.grd,pred.gpp2[1,], lty = 'dashed')
# lines(prep.grd, pred.gpp2[2,], lty = 1, lwd = 3)


#plot inset for precipitation sensitity for each type

x11 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "DBF" | reg.abs.gpp$IGBP2 == "MF")]*50*2
x12 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "ENF")]*50*2
x13 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "SHB")]*50*2
x14 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "GRA")]*50*2
x15 = reg.abs.gpp$prep.coef[which(reg.abs.gpp$IGBP2 == "CRO")]*50*2


se <- function(x) sqrt(var(x)/length(x))

D2 = data.frame(mean = c(mean(x11),mean(x12),mean(x13),mean(x14),mean(x15)),
				se = c(se(x11),se(x12),se(x13),se(x14),se(x15)),
				type = c("DBF/MF","ENF","Shrub","Grass","Crop"))
	

ylab = 	expression(paste(beta ["temporal"]) ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "/100mm"))
color1 = c("green", "red")

p1 <- ggplot(data=D2, aes(x=type, y=mean)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 	
	# ylab(ylab) + 	 
    ylab(expression(paste(beta ["temporal"]))) + 
	xlab("") + # Set axis labels
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
    theme(strip.text.y = element_text(size=12)) 
	
library(ggplot2)
require(grid)
print(p1, vp=viewport(.75, .3, .4, .3)) #viewport, first two is the x,y coor; last 2 is the inset size
  
dev.off()

# plot prepcipation sensitivity
png("F:/zhihua/dataset/results2/er.prep.png",height = 2500, width = 3000, res = 300, units = "px")

par(mar=c(5,6,1,1))
#plot the relationship
plot(0, xlim = c(0, 1500), ylim = c(0,2300),bty='n',pch='',
     xlab='Precipation (mm)',
	 ylab=expression("TER" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "")),
     cex.axis = 2, cex.lab = 2)

#draw a box around the plot
box()
#plot Evergreen Forest
for(i in 1:length(flux.info3$IGBP2)){

if (flux.info3$IGBP2[i] == "ENF" ){
lines(reg.abs.mod.pred.er[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.er[[i]][,4]+dat.df.annual.mean$er[i],lwd = 1.5, lty = 2)

} else if (flux.info3$IGBP2[i] == "DBF") {
lines(reg.abs.mod.pred.er[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.er[[i]][,4]+dat.df.annual.mean$er[i],lwd = 1.5,col = 2)

} else if (flux.info3$IGBP2[i] == "MF") {
lines(reg.abs.mod.pred.er[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.er[[i]][,4]+dat.df.annual.mean$er[i],lwd = 1.5,col = 2)

} else if (flux.info3$IGBP2[i] == "SHB") {
lines(reg.abs.mod.pred.er[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.er[[i]][,4]+dat.df.annual.mean$er[i],lwd = 1.5, col=4)

} else if (flux.info3$IGBP2[i] == "GRA") {
lines(reg.abs.mod.pred.er[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.er[[i]][,4]+dat.df.annual.mean$er[i],lwd = 1.5, col=4)

} else if (flux.info3$IGBP2[i] == "CRO") {
lines(reg.abs.mod.pred.er[[i]][,3]+dat.df.annual.mean$prep[i], reg.abs.mod.pred.er[[i]][,4]+dat.df.annual.mean$er[i],lwd = 1.5, col=4)

}

# text(x = reg.abs.mod.pred.er[[i]][1,3]+dat.df.annual.mean$prep[i],
     # y = reg.abs.mod.pred.er[[i]][1,4]+dat.df.annual.mean$er[i],
	 # paste(factor(flux.info3[i,1]), factor(flux.info3[i,2])))


}


legend("bottomright", 
	   # inset=0.05, 
	   # legend = c("ENF","DBF","MF","SHB","GRA","CRO"),
	   legend = c("ENF","DBF/MF","Shrub/GRA/CRO"),
	   horiz=F,
	   lwd = 2,
	   col = c(1,2,4),
	   text.col = c(1,2,4),
	   cex = 1.5,
	   box.col = "transparent",
	   bg = "transparent")

### plot er
lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3, col = "red")
# add fill
polygon(c(rev(prep.grd), prep.grd), c(rev(pred.er2[3,]), pred.er2[1,]), 
  col = rgb(1,0,0,0.3), 
  border = NA)
# intervals
# lines(prep.grd,pred.er2[3,], lty = 'dashed')
# lines(prep.grd,pred.er2[1,], lty = 'dashed')
# lines(prep.grd, pred.er2[2,], lty = 1, lwd = 3)


#plot inset for precipitation sensitity for each type
x11 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "DBF" | reg.abs.er$IGBP2 == "MF")]*50*2
x12 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "ENF")]*50*2
x13 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "SHB")]*50*2
x14 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "GRA")]*50*2
x15 = reg.abs.er$prep.coef[which(reg.abs.er$IGBP2 == "CRO")]*50*2


se <- function(x) sqrt(var(x)/length(x))

D2 = data.frame(mean = c(mean(x11),mean(x12),mean(x13),mean(x14),mean(x15)),
				se = c(se(x11),se(x12),se(x13),se(x14),se(x15)),
				type = c("DBF/MF","ENF","Shrub","Grass","Crop"))
	

ylab = 	expression(paste(beta ["temporal"]) ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ "/100mm"))
color1 = c("green", "red")

p1 <- ggplot(data=D2, aes(x=type, y=mean)) + 
    geom_bar(colour="black", stat="identity",
             position=position_dodge(),
             size=.3) + 	
	# ylab(ylab) + 	 
    ylab(expression(paste(beta ["temporal"]))) + 
	xlab("") + # Set axis labels
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
    theme(strip.text.y = element_text(size=12)) 
	
library(ggplot2)
require(grid)
print(p1, vp=viewport(.4, 0.8, .4, .3)) #viewport, first two is the x,y coor; last 2 is the inset size

dev.off()



## DO NOT USE THIS ONE

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




### DOES NOT WORK FROM HERE
D3 = rbind(delt, D1)
write.csv(D3, "F:/zhihua/dataset/results2/EC.sensitivity.csv")

reg.abs.mod.pred1.er.model.coef = c() #store model coefficients
colnames(rmse.er) = c("rmse","mae","aic","r2") 
colnames(rmse.gpp) = c("rmse","mae","aic","r2") 
rmse1 = rbind(data.frame(rmse.gpp, flux = "gpp"),data.frame(rmse.er, flux = "er"))
write.csv(rmse1, "F:/zhihua/dataset/results2/EC.rmse1.csv")

colnames(reg.abs.mod.pred1.er.model.coef) = c("int","t2","t","p2","p") 
colnames(reg.abs.mod.pred1.gpp.model.coef) = c("int","t2","t","p2","p") 
coef1 = rbind(data.frame(reg.abs.mod.pred1.gpp.model.coef, flux = "gpp"),data.frame(reg.abs.mod.pred1.er.model.coef, flux = "er"))
write.csv(coef1, "F:/zhihua/dataset/results2/EC.coef1.csv")


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
geom_text(x = 0.5, y = 70, label = "2)",size=6

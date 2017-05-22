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





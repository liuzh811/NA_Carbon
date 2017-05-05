
# 5/3/2017 plotting CO2 NEE/MOD GPP/GOME2 SIF at global scale

library(rgdal)
library(raster)
library(rasterVis)
require(ncdf4)

usa.state = readOGR(dsn="D:\\zhihua\\dataset\\ecoregion", layer = "usa.state")
usa.state = usa.state[-which(usa.state$NAME_1 == "Hawaii" | usa.state$NAME_1 == "Alaska"),]

# read into CO2 inversion
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

# add 2000 estimate
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


# calculate the median value ##WITHOUT JENA
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

# read into MOD 17 GPP

Work.Dir <- "F:/zhihua/dataset/mod17a2"
setwd(Work.Dir)

rasterOptions(tmpdir="F:/zhihua/dataset/mod17a2/TempDir")

gpp.annual = list()

gpp.files = list.files(path = "./GeoTIFF_0.05degree", pattern = "*.tif$")

for (i in 2000:2014){

gpp.file = which(as.numeric(substr(gpp.files, 13,16)) == i)

gpp.r <- stack(paste("./GeoTIFF_0.05degree/", gpp.files[gpp.file], sep = ""))
# gpp.r = crop(gpp.r, na.ext)
gpp.r[gpp.r>30000] = NA #remove filled value

gpp.r = gpp.r*0.1
names(gpp.r) <- paste("M", 1:12, sep = "")
writeRaster(gpp.r,paste("./processed/global.c.gpp.monthly.", i, ".grd", sep = ""),overwrite=TRUE) 

gpp.r = calc(gpp.r, sum)

gpp.annual[[i - 1999]] <- gpp.r

print(paste("Finish calculating annual GPP for ", i, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

gpp.annual= stack(gpp.annual)
names(gpp.annual) <- paste("Y", 2000:2014, sep = "")
writeRaster(gpp.annual,paste("./processed/global.c.gpp.annual", ".grd", sep = ""),overwrite=TRUE) 

sort(sapply(ls(), function(x){object.size(get(x))}))

# read into GOME SIF product.

sif = stack("F:/zhihua/dataset/gome2_sif/sif.annual.2007-2015.grd") 

# plot
gpp.annual.mean = calc(gpp.annual, mean)

mask.land = gpp.annual.mean > 0
mask.land1 = aggregate(mask.land, fact=20, fun=mean)

ensm.nee1.mean = calc(ensm.nee1, mean)*mask.land1
ensm.nee2.mean = calc(ensm.nee2, mean)*mask.land1

ct.nee.mean = calc(ct.nee, mean)*mask.land1
cte.nee.mean = calc(cte.nee, mean)*mask.land1

mask.land2 = aggregate(mask.land, fact=10, fun=mean)
sif.mean = calc(sif, mean)*mask.land2

# get global adminstrate boundary
install.packages("rworldmap")
library(rworldmap)

newmap <- getMap(resolution = "coarse")

png("F:/zhihua/dataset/results/plot20170505.png",height = 4500, width = 3000, res = 300, units = "px")

# plot start from here
par(mfrow=c(3,1),mar=c(0,0,0,0)+.1)

####### GPP plot 
plot(gpp.annual.mean, 
                    # zlim=c(-1,1),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
plot(newmap, add = TRUE, lty = 2)
tx = expression("Mean Annual GPP by MODIS" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ ", 2000 - 2014"))
text(x = 0, y = -85, tx, cex = 2)
# text(x = 0, y = -85, "Mean Annual GPP by MODIS (g C/m-2*yr, 2000-2014)", cex = 2)

#add legend
plot(gpp.annual.mean, 
                    # zlim=c(-1,1),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					legend.only=TRUE,
					legend.width=0.5, legend.shrink=0.5,
					axis.args=list(cex.axis=1.5),
					legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
					smallplot=c(0.13,0.16, 0.35,0.55))
				par(mar = par("mar"))

####### SIF plot 
my.colors = colorRampPalette(c("red", "white", "blue"))
plot( sif.mean, 
                    # zlim=c(-100,100),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
plot(newmap, add = TRUE, lty = 2)
tx = expression("Mean Annual SIF by GOME2" ~ (mu~W ~ m^{-2} ~ mu~m^{-1}~ ", 2007 - 2014"))
text(x = 0, y = -85, tx, cex = 2)

#add legend
plot(sif.mean, 
                    # zlim=c(-1,1),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					legend.only=TRUE,
					legend.width=0.5, legend.shrink=0.5,
					axis.args=list(cex.axis=1.5),
					legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
					smallplot=c(0.13,0.16, 0.35,0.55))
				par(mar = par("mar"))
				
		

####### NEE plot 
my.colors = colorRampPalette(c("red", "white", "blue"))

ensm.nee1.mean1 = ensm.nee1.mean
ensm.nee1.mean1[ensm.nee1.mean1 > 200] = 200
ensm.nee1.mean1[ensm.nee1.mean1 < -100] = -100

plot(               
					# ct.nee.mean,
					ensm.nee1.mean1, 
                    # zlim=c(-100,100),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
plot(newmap, add = TRUE, lty = 2)

tx = expression("Mean Annual NEE by CO2 Inversions" ~ (g ~ C ~ m^{-2} ~ yr ^{-1}~ ", 2000 - 2014"))
text(x = 0, y = -85, tx, cex = 2)
# text(x = 0, y = -85, "Mean Annual NEE by CO2 Inversions (g C/m-2*yr, 2000-2014)", cex = 2)

#add legend
plot(ensm.nee1.mean1, 
                    # zlim=c(-1,1),
					# col = my.colors(100), 
					col = rainbow(n = 100), 
					legend.only=TRUE,
					legend.width=0.5, legend.shrink=0.5,
					axis.args=list(cex.axis=1.5),
					legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
					smallplot=c(0.13,0.16, 0.35,0.55))
				par(mar = par("mar"))
				
dev.off()
		

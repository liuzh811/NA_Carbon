
##### TRENDY results analysis
##### calcuate mean GPP/NEE/ER, and calculate their relationship with T/P

##### In 2017/3/20, change the unit to g C/m2/s, all the model can be compared with each other
library(rgdal)
library(raster)
library(rasterVis)
require(ncdf4)

Work.Dir <- "F:/zhihua/dataset/trendy"
setwd(Work.Dir)

na.state = readOGR(dsn="F:\\zhihua\\dataset\\ecoregion", layer = "na.state2")
na.state = na.state[-which(na.state$NAME_1 == "Prince Edward Island"|
                           na.state$NAME_1 == "Hawaii"),]
# create a raster mask for North America
ext = extent(na.state)
R.grd = raster()   
R.grd[] <- 1						   
R.grd2 = crop(R.grd, na.state)
R.grd2 = rasterize(na.state, R.grd2, field = 1)
#

usa.state = na.state[which(na.state$ISO == "USA" & na.state$NAME_1 != "Alaska"), ]

# correlation between GPP:NEE
R1.stack <- list()
P1.stack <- list()
# correlation between GPP:Ra
R2.stack <- list()
P2.stack <- list()
# correlation between GPP:Rh
R3.stack <- list()
P3.stack <- list()

GPP.stack2 <- list()
NBP.stack2 <- list()
ER.stack2 <- list()
Ra.stack2 <- list()
Rh.stack2 <- list()

days.mon = c(31,28,31,30,31,30,31,31,30,31,30,31)
sec.mon = days.mon*86400

##########################1: CLM4C #################################
#days since 1901-1, monthly, to 2010-12
#time  Size:1368
#so 190101 == 1; 190102 == 2
# calculate 1985-2014 (30 years), 198501 = (1985-1901)*12+1
# year 1985 is (1985-1901)*12 + 1:12
# year x is (x-1901)*12 + 1:12

gpp.list = list()
nc1.file = nc_open("./download/clm4c/gpp.nc")

v2 <- nc1.file$var[[1]]

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

gpp.list[[yr - 1999]] <- r1
# gpp.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[1]] <- gpp.list

nee.list = list()
nc1.file = nc_open("./download/clm4c/nbp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
NBP.stack2[[1]] <- nee.list


nee.list = list()
nc1.file = nc_open("./download/clm4c/ra.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Ra.stack2[[1]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/clm4c/rh.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Rh.stack2[[1]] <- nee.list

ER.stack2[[1]] <- Rh.stack2[[1]] + Ra.stack2[[1]]

#calcuate correlationship

###### calculate relationship between GPP/NDVI/EVI and NEE
pts.sp = Ex.pts.all(nee.list[[1]]) #get the point locations

#extract NEE
nee.df = raster::extract(NBP.stack2[[1]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[1]], pts.sp)
ra.df = raster::extract(Ra.stack2[[1]], pts.sp)
rh.df= raster::extract(Rh.stack2[[1]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[1]] <- dat1.corr
P1.stack[[1]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[1]] <- dat1.corr
P2.stack[[1]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[1]] <- dat1.corr
P3.stack[[1]] <- dat1.p

##########################2: CLM4CN #################################
#days since 1901-1, monthly, to 2010-12
#time  Size:1368
#so 190101 == 1; 190102 == 2
# calculate 1985-2014 (30 years), 198501 = (1985-1901)*12+1
# year 1985 is (1985-1901)*12 + 1:12
# year x is (x-1901)*12 + 1:12

gpp.list = list()
nc1.file = nc_open("./download/clm4cn/gpp.nc")
print(nc1.file)
v2 <- nc1.file$var[[1]]

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

gpp.list[[yr - 1999]] <- r1
# gpp.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[2]] <- gpp.list

nee.list = list()
nc1.file = nc_open("./download/clm4cn/nbp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
NBP.stack2[[2]] <- nee.list


nee.list = list()
nc1.file = nc_open("./download/clm4cn/ra.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Ra.stack2[[2]] <- nee.list


nee.list = list()
nc1.file = nc_open("./download/clm4cn/rh.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
# change unit, try umolCO2 m-2 s-1 first # NOT THIS UNIT
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Rh.stack2[[2]] <- nee.list

#calcuate correlationship
###### calculate relationship between GPP/NDVI/EVI and NEE
pts.sp = Ex.pts.all(nee.list[[1]]) #get the point locations

nee.df = raster::extract(NBP.stack2[[2]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[2]], pts.sp)
ra.df = raster::extract(Ra.stack2[[2]], pts.sp)
rh.df= raster::extract(Rh.stack2[[2]], pts.sp)

dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])
R1.stack[[2]] <- dat1.corr
P1.stack[[2]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])
R2.stack[[2]] <- dat1.corr
P2.stack[[2]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[2]] <- dat1.corr
P3.stack[[2]] <- dat1.p

##########################3: HYLAND #################################
## unknow unit, do not use

gpp.list = list()
nc1.file = nc_open("./download/hyland/gpp.nc")
print(nc1.file)

v2 <- nc1.file$var[[1]]

# Lat = v2$dim[[2]]$vals
# Lon = v2$dim[[1]]$vals
# res.x = Lon[2]-Lon[1]
# res.y = Lat[2]-Lat[1]

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- aperm(data2[,,c(100:110)], c(2,1,3))
data2.r = brick(data2,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r = flip(data2.r, direction = "y")
# Change the xlim from 0~360 to -180 ~ 180
data2.r1 = list()
for (i in 1:nlayers(data2.r)){
r1 = data2.r[[i]]
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

data2.r1[[i]] <- r1
}
data2.r1 = stack(data2.r1)
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r1 = data2.r1*60*60*24*365*1000
gpp.list = stack(data2.r1)
GPP.stack2[[3]] <- gpp.list

nee.list = list()
nc1.file = nc_open("./download/hyland/nbp.nc")
v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- aperm(data2[,,c(100:110)], c(2,1,3))
data2.r = brick(data2,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r = flip(data2.r, direction = "y")
# Change the xlim from 0~360 to -180 ~ 180
data2.r1 = list()
for (i in 1:nlayers(data2.r)){
r1 = data2.r[[i]]
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

data2.r1[[i]] <- r1
}
data2.r1 = stack(data2.r1)
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r1 = data2.r1*60*60*24*365*1000
nee.list = -1*stack(data2.r1)
NBP.stack2[[3]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/hyland/rh.nc")
v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- aperm(data2[,,c(100:110)], c(2,1,3))
data2.r = brick(data2,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r = flip(data2.r, direction = "y")
# Change the xlim from 0~360 to -180 ~ 180
data2.r1 = list()
for (i in 1:nlayers(data2.r)){
r1 = data2.r[[i]]
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

data2.r1[[i]] <- r1
}
data2.r1 = stack(data2.r1)
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r1 = data2.r1*60*60*24*365*1000
Rh.stack2[[3]] <- data2.r1

Ra.stack2[[3]] = GPP.stack2[[3]] - Rh.stack2[[3]] - NBP.stack2[[3]]

#calcuate correlationship

###### calculate relationship between GPP/NDVI/EVI and NEE
nee.df = raster::extract(NBP.stack2[[3]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[3]], pts.sp)
ra.df = raster::extract(Ra.stack2[[3]], pts.sp)
rh.df= raster::extract(Rh.stack2[[3]], pts.sp)

dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
dat1.corr = Point2raster(dat1.df1[1,], raster = gpp.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = gpp.list[[1]])
R1.stack[[3]] <- dat1.corr
P1.stack[[3]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = gpp.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = gpp.list[[1]])

R2.stack[[3]] <- dat1.corr
P2.stack[[3]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
dat1.corr = Point2raster(dat1.df1[1,], raster = gpp.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = gpp.list[[1]])

R3.stack[[3]] <- dat1.corr
P3.stack[[3]] <- dat1.p
##########################4: LPJ #################################
# units: kg C m-2 month-1

gpp.list = list()
nc1.file = nc_open("./download/lpj/gpp.nc")
print(nc1.file)

v2 <- nc1.file$var[[1]]

Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r[data2.r == -99999] = NA
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000

# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2
	
gpp.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[4]] <- gpp.list

nee.list = list()
nc1.file = nc_open("./download/lpj/nbp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r[data2.r == -99999] = NA
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2
		
nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
NBP.stack2[[4]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/lpj/LPJ_S2_ra.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r[data2.r == -99999] = NA
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2
		
nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Ra.stack2[[4]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/lpj/LPJ_S2_rh.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r[data2.r == -99999] = NA
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2
		
nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Rh.stack2[[4]] <- nee.list

#calcuate correlationship
pts.sp = Ex.pts.all(nee.list[[1]]) #get the point locations

#extract NEE
nee.df = raster::extract(NBP.stack2[[4]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[4]], pts.sp)
ra.df = raster::extract(Ra.stack2[[4]], pts.sp)
rh.df= raster::extract(Rh.stack2[[4]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[4]] <- dat1.corr
P1.stack[[4]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[4]] <- dat1.corr
P2.stack[[4]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[4]] <- dat1.corr
P3.stack[[4]] <- dat1.p
##########################5: LPJ-GUESS #################################
## unknown units, do not use

gpp.list = list()
nc1.file = nc_open("./download/lpj_guess/gpp.nc")
print(nc1.file)

v2 <- nc1.file$var[[1]]

Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
		
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2

gpp.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[5]] <- gpp.list

nee.list = list()
nc1.file = nc_open("./download/lpj_guess/nbp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
	
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2
		
nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
NBP.stack2[[5]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/lpj_guess/LPJ_GUESS_s2_ra.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
	
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2
		
nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Ra.stack2[[5]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/lpj_guess/LPJ_GUESS_s2_rh.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=-180, xmx=180, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
	
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# aggregate to 1 degree
data2.r	= aggregate(data2.r, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
data2.r = data2.r*R.grd2
		
nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Rh.stack2[[5]] <- nee.list
#calcuate correlationship
pts.sp = Ex.pts.all(nee.list[[1]]) #get the point locations

#extract NEE
nee.df = raster::extract(NBP.stack2[[5]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[5]], pts.sp)
ra.df = raster::extract(Ra.stack2[[5]], pts.sp)
rh.df= raster::extract(Rh.stack2[[5]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[5]] <- dat1.corr
P1.stack[[5]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[5]] <- dat1.corr
P2.stack[[5]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[5]] <- dat1.corr
P3.stack[[5]] <- dat1.p

##########################6: OCN #################################
## units: kg C m^-2 s^-1

gpp.list = list()
nc1.file = nc_open("./download/ocn/gpp.nc")
print(nc1.file)

v2 <- nc1.file$var[[1]]

Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

gpp.list[[yr - 1999]] <- r1
# gpp.list[[yr - 1999]] <- data2.r
	
print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[6]] <- nee.list


nee.list = list()
nc1.file = nc_open("./download/ocn/nbp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
NBP.stack2[[6]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/ocn/OCN_S2_ra.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Ra.stack2[[6]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/ocn/OCN_S2_rh.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Rh.stack2[[6]] <- nee.list

#calcuate correlationship

#extract NEE
nee.df = raster::extract(NBP.stack2[[6]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[6]], pts.sp)
ra.df = raster::extract(Ra.stack2[[6]], pts.sp)
rh.df= raster::extract(Rh.stack2[[6]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[6]] <- dat1.corr
P1.stack[[6]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[6]] <- dat1.corr
P2.stack[[6]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[6]] <- dat1.corr
P3.stack[[6]] <- dat1.p

##########################7: ORCHIDEE #################################
## unknown units, do not use

gpp.list = list()
nc1.file = nc_open("./download/orchidee/gpp.nc")
print(nc1.file)

v2 <- nc1.file$var[[1]]
Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- aperm(data2[,,c(100:110)], c(2,1,3))
data2.r = brick(data2,xmn=min(Lon), xmx=max(Lon)+res.x, ymn=min(Lat), ymx=max(Lat) + res.y, 
        crs = "+proj=longlat +datum=WGS84")
		
data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
data2.r = data2.r*60*60*24*365*1000

#aggregate resulution from 0.5 to 1, and crop to North America
data2.r1 = list()
for(i in 1:nlayers(data2.r)){
r1 = aggregate(data2.r[[i]], fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
r1 = r1*R.grd2
data2.r1[[i]] <- r1	
}
gpp.list = stack(data2.r1)
GPP.stack2[[7]] <- gpp.list

nee.list = list()
nc1.file = nc_open("./download/orchidee/nbp.nc")
v2 <- nc1.file$var[[1]]
Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- aperm(data2[,,c(100:110)], c(2,1,3))
data2.r = brick(data2,xmn=min(Lon), xmx=max(Lon)+res.x, ymn=min(Lat), ymx=max(Lat) + res.y, 
        crs = "+proj=longlat +datum=WGS84")
		
data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
data2.r = data2.r*60*60*24*365*1000

#aggregate resulution from 0.5 to 1, and crop to North America
data2.r1 = list()
for(i in 1:nlayers(data2.r)){
r1 = aggregate(data2.r[[i]], fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
r1 = r1*R.grd2
data2.r1[[i]] <- r1	
}
nee.list = stack(data2.r1)
NBP.stack2[[7]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/orchidee/ORCHIDEE_S2_rh_y.nc")
v2 <- nc1.file$var[[1]]
Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- aperm(data2[,,c(100:110)], c(2,1,3))
data2.r = brick(data2,xmn=min(Lon), xmx=max(Lon)+res.x, ymn=min(Lat), ymx=max(Lat) + res.y, 
        crs = "+proj=longlat +datum=WGS84")
		
data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA
data2.r = data2.r*60*60*24*365*1000

#aggregate resulution from 0.5 to 1, and crop to North America
data2.r1 = list()
for(i in 1:nlayers(data2.r)){
r1 = aggregate(data2.r[[i]], fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
r1 = r1*R.grd2
data2.r1[[i]] <- r1	
}
nee.list = stack(data2.r1)
Rh.stack2[[7]] <- nee.list

Ra.stack2[[7]] = GPP.stack2[[7]] - Rh.stack2[[7]] - NBP.stack2[[7]]  
#extract NEE
nee.df = raster::extract(NBP.stack2[[7]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[7]], pts.sp)
ra.df= raster::extract(Ra.stack2[[7]], pts.sp)
rh.df= raster::extract(Rh.stack2[[7]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[7]] <- dat1.corr
P1.stack[[7]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[7]] <- dat1.corr
P2.stack[[7]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[7]] <- dat1.corr
P3.stack[[7]] <- dat1.p
##########################8: triffid #################################
## units: kgC/m2/s
gpp.list = list()
nc1.file = nc_open("./download/triffid/gpp.nc")

v2 <- nc1.file$var[[1]]
Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

gpp.list[[yr - 1999]] <- r1
# gpp.list[[yr - 1999]] <- data2.r
	
print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[8]] <- gpp.list


nee.list = list()
nc1.file = nc_open("./download/triffid/nbp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
NBP.stack2[[8]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/triffid/ra.nc")
v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Ra.stack2[[8]] <- nee.list


nee.list = list()
nc1.file = nc_open("./download/triffid/rh.nc")
v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# resample to 1 degree
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Rh.stack2[[8]] <- nee.list

#calcuate correlationship
#extract NEE
nee.df = raster::extract(NBP.stack2[[8]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[8]], pts.sp)
ra.df = raster::extract(Ra.stack2[[8]], pts.sp)
rh.df= raster::extract(Rh.stack2[[8]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[8]] <- dat1.corr
P1.stack[[8]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[8]] <- dat1.corr
P2.stack[[8]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[8]] <- dat1.corr
P3.stack[[8]] <- dat1.p
##########################9: vegas #################################
## unit: kg/m2/s

gpp.list = list()
nc1.file = nc_open("./download/vegas/VEGAS_S2_Monthly_gpp.nc")

v2 <- nc1.file$var[[1]]

Lat = v2$dim[[2]]$vals
Lon = v2$dim[[1]]$vals
res.x = abs(Lon[2]-Lon[1])
res.y = abs(Lat[2]-Lat[1])

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# aggregate to 1 degree
r1	= aggregate(r1, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
r1= r1*R.grd2
	
gpp.list[[yr - 1999]] <- r1
# gpp.list[[yr - 1999]] <- data2.r


print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[9]] <- gpp.list


nee.list = list()
nc1.file = nc_open("./download/vegas/VEGAS_S2_Monthly_totresp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# aggregate to 1 degree
r1	= aggregate(r1, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
r1= r1*R.grd2
	
nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
nee.list = gpp.list-nee.list
NBP.stack2[[9]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/vegas/VEGAS_S2_Monthly_autoresp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# aggregate to 1 degree
r1	= aggregate(r1, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
r1= r1*R.grd2
	
nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Ra.stack2[[9]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/vegas/VEGAS_S2_Monthly_hetresp.nc")

v2 <- nc1.file$var[[1]]
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
for (yr in 2000:2010){

data2.tmp = data2[,,c((yr-1901)*12 + 1:12)]
data2.tmp <- aperm(data2.tmp, c(2,1,3))
data2.r = brick(data2.tmp,xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")

data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == -99999] = NA	
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# aggregate to 1 degree
r1	= aggregate(r1, fact=2, fun=mean, expand=TRUE, na.rm=TRUE)
r1= r1*R.grd2
	
nee.list[[yr - 1999]] <- r1
# nee.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
nee.list = stack(nee.list)
Rh.stack2[[9]] <- nee.list

#calcuate correlationship
#extract NEE
nee.df = raster::extract(NBP.stack2[[9]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[9]], pts.sp)
ra.df = raster::extract(Ra.stack2[[9]], pts.sp)
rh.df= raster::extract(Rh.stack2[[9]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[9]] <- dat1.corr
P1.stack[[9]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[9]] <- dat1.corr
P2.stack[[9]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[9]] <- dat1.corr
P3.stack[[9]] <- dat1.p
##########################10: sdgvm #################################
## unknown units, do not use

gpp.list = list()
nc1.file = nc_open("./download/sdgvm/gpp.nc")

v2 <- nc1.file$var[[1]]

data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array

for (yr in 2000:2010){

data2.tmp = data2[c((yr-1901)*12 + 1:12),,]

tmp.list = list()
for (i in 1:12){
data2.r = raster(data2.tmp[i,,],xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r = flip(data2.r, direction = "y")
data2.r[data2.r == 99999] = NA
tmp.list[[i]] <- data2.r
}
data2.r = stack(tmp.list)
# Seem to be Kg C/m2/s: multiple *60*60*24*365*1000, change to g C/m2/yr
data2.r = data2.r*sec.mon
data2.r = calc(data2.r, sum, na.rm = TRUE)*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees
# crop to North America
r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

gpp.list[[yr - 1999]] <- r1
# gpp.list[[yr - 1999]] <- data2.r

print(paste("Finish calculating for year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}
gpp.list = stack(gpp.list)
GPP.stack2[[10]] <- gpp.list

nee.list = list()
nc1.file = nc_open("./download/sdgvm/nbp.nc")
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- data2[c(100:110),,]

for (i in 1:11){
data2.r = raster(data2[i,,],xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r[data2.r == -99999] = NA	
data2.r = flip(data2.r, direction = "y")
data2.r = data2.r*60*60*24*365*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[i]] <- r1
# nee.list[[yr - 1999]] <- data2.r

}

nee.list = stack(nee.list)
NBP.stack2[[10]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/sdgvm/SDGVM_S2_ra.nc")
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- data2[c(100:110),,]

for (i in 1:11){
data2.r = raster(data2[i,,],xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r[data2.r == -99999] = NA	
data2.r = flip(data2.r, direction = "y")
data2.r = data2.r*60*60*24*365*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[i]] <- r1
# nee.list[[yr - 1999]] <- data2.r

}

nee.list = stack(nee.list)
Ra.stack2[[10]] <- nee.list

nee.list = list()
nc1.file = nc_open("./download/sdgvm/SDGVM_S2_rh.nc")
data2 <- ncvar_get( nc1.file, v2 ) #data2 is an 3-d array
data2 <- data2[c(100:110),,]

for (i in 1:11){
data2.r = raster(data2[i,,],xmn=0, xmx=360, ymn=-90, ymx=90, 
        crs = "+proj=longlat +datum=WGS84")
data2.r[data2.r == -99999] = NA	
data2.r = flip(data2.r, direction = "y")
data2.r = data2.r*60*60*24*365*1000
# Change the xlim from 0~360 to -180 ~ 180
r1 = data2.r
r2 = r1
nr = nrow(r1); nc = ncol(r1)
r1[,c(1:(nc/2))] = r2[,c((nc/2+1):nc)]
r1[,c((nc/2+1):nc)] = r2[,c(1:(nc/2))]
r1 <- shift(r1, x= -180, y=0)#shift by 180 degrees

r1 = resample(r1, R.grd, method = "ngb")
r1 = r1*R.grd2

nee.list[[i]] <- r1
# nee.list[[yr - 1999]] <- data2.r

}

nee.list = stack(nee.list)
Rh.stack2[[10]] <- nee.list

#calcuate correlationship
#extract NEE
nee.df = raster::extract(NBP.stack2[[10]], pts.sp)
gpp.df= raster::extract(GPP.stack2[[10]], pts.sp)
ra.df = raster::extract(Ra.stack2[[10]], pts.sp)
rh.df= raster::extract(Rh.stack2[[10]], pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R1.stack[[10]] <- dat1.corr
P1.stack[[10]] <- dat1.p

#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R2.stack[[10]] <- dat1.corr
P2.stack[[10]] <- dat1.p
#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

R3.stack[[10]] <- dat1.corr
P3.stack[[10]] <- dat1.p

################### plot  ####################
# for correlationship

usa.state = na.state[which(na.state$ISO == "USA" & na.state$NAME_1 != "Alaska"), ]
na.ext <- extent(-180,-48,15,85)

mask1 = R.grd2
pts.sp = Ex.pts.all(mask1) #get the point locations

mask2 = rasterize(usa.state, crop(mask1,usa.state))
mask2 = mask2 > 0
mask2[mask2 == 0] = NA

## plot GPP/NEE r
R1.stack2 = crop(stack(R1.stack),usa.state)
P1.stack2 = crop(stack(P1.stack),usa.state)

R1.stack2 = R1.stack2*mask2
P1.stack2 = P1.stack2*mask2

# writeRaster(R1.stack2,"D:/zhihua/dataset/trendy/R.stack3.grd",overwrite=TRUE) 
# writeRaster(R1.stack2,"D:/zhihua/dataset/trendy/P.stack3.grd",overwrite=TRUE) 

#plot individual correlationship map
Mod.name = c("CLM4C", "CLM4CN", "HYLAND", "LPJ", "LPJ-GUESS", "OCN", "ORCHIDEE","TRIFFID","VEGAS","SDGVM")

png("D:/zhihua/dataset/results/productivity.nee.cor.trendy.png",height = 2500, width = 3500, res = 300, units = "px")

par(mfrow=c(4,3),mar=c(0,0,0,0)+.1)
my.colors = colorRampPalette(c("blue", "white", "red"))

for (i in 1:10){
pts.sp.sig1 = Ex.pts(P1.stack2[[i]], sig.level = 0.1) #extract significant relation points

plot(R1.stack2[[i]], zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, Mod.name[i], cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)
}


plot(calc(R1.stack2, mean, na.rm = TRUE), zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, "Mean r", cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)

P = calc(P1.stack2, function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 7] = 1
P1[P > 7] = 0.05
pts.sp.sig1 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points
plot(pts.sp.sig1, add = TRUE, cex = 0.1)


plot.new()
plot(R1.stack2[[i]], zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=1, legend.shrink=0.75,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="r", side=4, font=2, line=2.5, cex=1.2),
    smallplot=c(0.3,0.5, 0.2,0.7))
	 par(mar = par("mar"))

legend(x = -100, y = 43, legend = "Correlation between NBP and GPP (2000-2010)",  cex = 1.5,  box.lwd = 0,box.col = "white",bg = "white")


dev.off()

## plot GPP/Rh r
R3.stack2 = crop(stack(R3.stack),usa.state)
P3.stack2 = crop(stack(P3.stack),usa.state)

R3.stack2 = R3.stack2*mask2
P3.stack2 = P3.stack2*mask2

# writeRaster(R1.stack2,"D:/zhihua/dataset/trendy/R.stack3.grd",overwrite=TRUE) 
# writeRaster(R1.stack2,"D:/zhihua/dataset/trendy/P.stack3.grd",overwrite=TRUE) 

#plot individual correlationship map
Mod.name = c("CLM4C", "CLM4CN", "HYLAND", "LPJ", "LPJ-GUESS", "OCN", "ORCHIDEE","TRIFFID","VEGAS","SDGVM")

png("F:/zhihua/dataset/results/productivity.rh.cor.trendy.png",height = 2500, width = 3500, res = 300, units = "px")

par(mfrow=c(4,3),mar=c(0,0,0,0)+.1)
my.colors = colorRampPalette(c("blue", "white", "red"))

for (i in 1:10){
pts.sp.sig1 = Ex.pts(P3.stack2[[i]], sig.level = 0.1) #extract significant relation points

plot(R3.stack2[[i]], zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, Mod.name[i], cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)
}


plot(calc(R3.stack2, mean, na.rm = TRUE), zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, "Mean r", cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)

P = calc(P3.stack2, function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 7] = 1
P1[P > 7] = 0.05
pts.sp.sig1 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points
plot(pts.sp.sig1, add = TRUE, cex = 0.1)


plot.new()
plot(R3.stack2[[i]], zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=1, legend.shrink=0.75,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="r", side=4, font=2, line=2.5, cex=1.2),
    smallplot=c(0.3,0.5, 0.2,0.7))
	 par(mar = par("mar"))

legend(x = -100, y = 43, legend = "Correlation between NBP and GPP (2000-2010)",  cex = 1.5,  box.lwd = 0,box.col = "white",bg = "white")

dev.off()

## plot GPP/Ra r
R2.stack2 = crop(stack(R2.stack),usa.state)
P2.stack2 = crop(stack(P2.stack),usa.state)

R2.stack2 = R2.stack2*mask2
P2.stack2 = P2.stack2*mask2

# writeRaster(R1.stack2,"D:/zhihua/dataset/trendy/R.stack3.grd",overwrite=TRUE) 
# writeRaster(R1.stack2,"D:/zhihua/dataset/trendy/P.stack3.grd",overwrite=TRUE) 

#plot individual correlationship map
Mod.name = c("CLM4C", "CLM4CN", "HYLAND", "LPJ", "LPJ-GUESS", "OCN", "ORCHIDEE","TRIFFID","VEGAS","SDGVM")

png("F:/zhihua/dataset/results/productivity.ra.cor.trendy.png",height = 2500, width = 3500, res = 300, units = "px")

par(mfrow=c(4,3),mar=c(0,0,0,0)+.1)
my.colors = colorRampPalette(c("blue", "white", "red"))

for (i in 1:10){
pts.sp.sig1 = Ex.pts(P2.stack2[[i]], sig.level = 0.1) #extract significant relation points

plot(R2.stack2[[i]], zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, Mod.name[i], cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)

}


plot(calc(R2.stack2, mean, na.rm = TRUE), zlim=c(-1,1),col = my.colors(100), 
					main = "",
					legend=FALSE,
                    axes=FALSE,
                    box=FALSE)
text(x = -115, y = 28, "Mean r", cex = 2)
plot(usa.state, lwd = 1.5, add = TRUE)

P = calc(P2.stack2, function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 7] = 1
P1[P > 7] = 0.05
pts.sp.sig1 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points
plot(pts.sp.sig1, add = TRUE, cex = 0.1)


plot.new()
plot(R2.stack2[[i]], zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=1, legend.shrink=0.75,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="r", side=4, font=2, line=2.5, cex=1.2),
    smallplot=c(0.3,0.5, 0.2,0.7))
	 par(mar = par("mar"))

legend(x = -100, y = 43, legend = "Correlation between NBP and GPP (2000-2010)",  cex = 1.5,  box.lwd = 0,box.col = "white",bg = "white")

dev.off()


### plot mean r, GPP/RA
R2.stack3 = calc(R2.stack2, mean, na.rm = TRUE)
P = calc(P2.stack2, function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 7] = 1
P1[P > 7] = 0.05
pts.sp.sig2 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points
	
		
png("F:/zhihua/dataset/results/productivity.ra.cor.annual.usa.mean.trendy2.png",height = 1500, width = 2500, res = 300, units = "px")
	
my.colors = colorRampPalette(c("blue", "white", "red"))
plot(R2.stack3, zlim=c(-1,1),col = my.colors(100), 
					main = "TRENDY Mean r Between GPP and Ra",
					legend=FALSE,
                    # axes=FALSE,
                    # box=FALSE
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig2, add = TRUE, cex = 0.1)
					
plot(R2.stack3, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=0.5, legend.shrink=0.5,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
		 smallplot=c(0.73,0.76, 0.25,0.45))
	 par(mar = par("mar"))


dev.off()

### plot mean r, GPP/NEE
R1.stack3 = calc(R1.stack2, mean, na.rm = TRUE)
P = calc(P1.stack2, function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 7] = 1
P1[P > 7] = 0.05
pts.sp.sig2 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points
	
		
png("F:/zhihua/dataset/results/productivity.nee.cor.annual.usa.mean.trendy2.png",height = 1500, width = 2500, res = 300, units = "px")
	
my.colors = colorRampPalette(c("blue", "white", "red"))
plot(R1.stack3, zlim=c(-1,1),col = my.colors(100), 
					main = "TRENDY Mean r Between GPP and NEE",
					legend=FALSE,
                    # axes=FALSE,
                    # box=FALSE
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig2, add = TRUE, cex = 0.1)
					
plot(R1.stack3, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=0.5, legend.shrink=0.5,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
		 smallplot=c(0.73,0.76, 0.25,0.45))
	 par(mar = par("mar"))


dev.off()

### plot mean r, GPP/RH
R3.stack3 = calc(R3.stack2, mean, na.rm = TRUE)
P = calc(P3.stack2, function(x){length(which(x <= 0.1))})	
P1 = P
P1[P <= 7] = 1
P1[P > 7] = 0.05
pts.sp.sig2 = Ex.pts(P1, sig.level = 0.1) #extract significant relation points
	
		
png("F:/zhihua/dataset/results/productivity.rh.cor.annual.usa.mean.trendy2.png",height = 1500, width = 2500, res = 300, units = "px")
	
my.colors = colorRampPalette(c("blue", "white", "red"))
plot(R3.stack3, zlim=c(-1,1),col = my.colors(100), 
					main = "TRENDY Mean r Between GPP and Rh",
					legend=FALSE,
                    # axes=FALSE,
                    # box=FALSE
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig2, add = TRUE, cex = 0.1)
					
plot(R3.stack3, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=0.5, legend.shrink=0.5,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
		 smallplot=c(0.73,0.76, 0.25,0.45))
	 par(mar = par("mar"))


dev.off()


##### use mean GPP/NBP/Ra/Rh to calculate correlationship
GPP.mean = list()
NBP.mean = list()
Ra.mean = list()
Rh.mean = list()

for (i in 1:11){

GPP1 = stack(GPP.stack2[[1]][[i]],GPP.stack2[[2]][[i]],GPP.stack2[[3]][[i]],GPP.stack2[[4]][[i]],GPP.stack2[[5]][[i]],
            GPP.stack2[[6]][[i]],GPP.stack2[[7]][[i]],GPP.stack2[[8]][[i]],GPP.stack2[[9]][[i]],GPP.stack2[[10]][[i]])
NBP1 = stack(NBP.stack2[[1]][[i]],NBP.stack2[[2]][[i]],NBP.stack2[[3]][[i]],NBP.stack2[[4]][[i]],NBP.stack2[[5]][[i]],
            NBP.stack2[[6]][[i]],NBP.stack2[[7]][[i]],NBP.stack2[[8]][[i]],NBP.stack2[[9]][[i]],NBP.stack2[[10]][[i]])
Ra1 = stack(Ra.stack2[[1]][[i]],Ra.stack2[[2]][[i]],Ra.stack2[[3]][[i]],Ra.stack2[[4]][[i]],Ra.stack2[[5]][[i]],
            Ra.stack2[[6]][[i]],Ra.stack2[[7]][[i]],Ra.stack2[[8]][[i]],Ra.stack2[[9]][[i]],Ra.stack2[[10]][[i]])
Rh1 = stack(Rh.stack2[[1]][[i]],Rh.stack2[[2]][[i]],Rh.stack2[[3]][[i]],Rh.stack2[[4]][[i]],Rh.stack2[[5]][[i]],
            Rh.stack2[[6]][[i]],Rh.stack2[[7]][[i]],Rh.stack2[[8]][[i]],Rh.stack2[[9]][[i]],Rh.stack2[[10]][[i]])

GPP.mean[[i]] <- calc(GPP1, mean, na.rm = TRUE)	
NBP.mean[[i]] <- calc(NBP1, mean, na.rm = TRUE)			
Ra.mean[[i]] <- calc(Ra1, mean, na.rm = TRUE)			
Rh.mean[[i]] <- calc(Rh1, mean, na.rm = TRUE)			
		
}

GPP.mean = stack(GPP.mean)
NBP.mean = stack(NBP.mean)
Ra.mean = stack(Ra.mean)
Rh.mean = stack(Rh.mean)

#extract NEE
nee.df = raster::extract(NBP.mean, pts.sp)
gpp.df= raster::extract(GPP.mean, pts.sp)
ra.df = raster::extract(Ra.mean, pts.sp)
rh.df= raster::extract(Rh.mean, pts.sp)

#combine gpp and nee
dat1.df = data.frame(gpp.df, nee.df)
dat2.df = data.frame(gpp.df, ra.df)
dat3.df = data.frame(gpp.df, rh.df)

#calculate the relationship
dat1.df1 = apply(dat1.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

dat1.corr = dat1.corr*mask2
dat1.p = dat1.p*mask2
pts.sp.sig1 = Ex.pts(dat1.p, sig.level = 0.1) #extract significant relation points

png("F:/zhihua/dataset/results/productivity.nee.cor.annual.usa.mean.trendy.png",height = 1500, width = 2500, res = 300, units = "px")

my.colors = colorRampPalette(c("blue", "white", "red"))
plot(dat1.corr, zlim=c(-1,1),col = my.colors(100), 
					main = "r Between Mean GPP and NEE",
					legend=FALSE,
                    # axes=FALSE,
                    # box=FALSE
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)
					
plot(dat1.corr, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=0.5, legend.shrink=0.5,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
		 smallplot=c(0.73,0.76, 0.25,0.45))
	 par(mar = par("mar"))

dev.off()


#calculate the relationship
dat1.df1 = apply(dat2.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

dat1.corr = dat1.corr*mask2
dat1.p = dat1.p*mask2
pts.sp.sig1 = Ex.pts(dat1.p, sig.level = 0.1) #extract significant relation points

png("F:/zhihua/dataset/results/productivity.ra.cor.annual.usa.mean.trendy.png",height = 1500, width = 2500, res = 300, units = "px")

my.colors = colorRampPalette(c("blue", "white", "red"))
plot(dat1.corr, zlim=c(-1,1),col = my.colors(100), 
					main = "r Between Mean GPP and Ra",
					legend=FALSE,
                    # axes=FALSE,
                    # box=FALSE
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)
					
plot(dat1.corr, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=0.5, legend.shrink=0.5,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
		 smallplot=c(0.73,0.76, 0.25,0.45))
	 par(mar = par("mar"))

dev.off()



#calculate the relationship
dat1.df1 = apply(dat3.df, 1, cor.test1)                         
# change to raster: GPP:NEE
dat1.corr = Point2raster(dat1.df1[1,], raster = nee.list[[1]])
dat1.p = Point2raster(dat1.df1[2,], raster = nee.list[[1]])

dat1.corr = dat1.corr*mask2
dat1.p = dat1.p*mask2
pts.sp.sig1 = Ex.pts(dat1.p, sig.level = 0.1) #extract significant relation points

png("F:/zhihua/dataset/results/productivity.rh.cor.annual.usa.mean.trendy.png",height = 1500, width = 2500, res = 300, units = "px")

my.colors = colorRampPalette(c("blue", "white", "red"))
plot(dat1.corr, zlim=c(-1,1),col = my.colors(100), 
					main = "r Between Mean GPP and Rh",
					legend=FALSE,
                    # axes=FALSE,
                    # box=FALSE
					# 
					xlab = "Longtitude", ylab = "Latitude", cex.lab = 1.5,cex = 1.5
					)

plot(usa.state, lwd = 1.5, add = TRUE)
plot(pts.sp.sig1, add = TRUE, cex = 0.1)
					
plot(dat1.corr, zlim=c(-1,1),col = my.colors(100),
         legend.only=TRUE,
		 legend.width=0.5, legend.shrink=0.5,
		 axis.args=list(cex.axis=1.5),
		 legend.args=list(text="", side=4, font=2, line=2.5, cex=1.2),
		 smallplot=c(0.73,0.76, 0.25,0.45))
	 par(mar = par("mar"))

dev.off()


	

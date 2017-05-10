# functions used in the analysis

################ define functions ##################
#define a function to get all points for a raster
Ex.pts.all = function(x){
  proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
  pts = xyFromCell(x, seq(1, (nrow(x)*ncol(x))))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}

#define a function to get significant change points
Ex.pts = function(x, sig.level = 0.1){
  proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
  dry.evi.sig2 = x <= sig.level
  dry.evi.sig2[dry.evi.sig2 ==0] = NA
  dry.evi.sig2.xy = rasterToPoints(dry.evi.sig2) #get raster coordinate, from left to right, top to bottom
  dry.evi.sig2.xy = data.frame(dry.evi.sig2.xy)
  dry.evi.sig2.sp <- SpatialPoints(coords = cbind(dry.evi.sig2.xy$x,dry.evi.sig2.xy$y),proj4string = CRS(proj.geo))
  return(dry.evi.sig2.sp)
}


# define a function to find the correlationship between NEE and GPP

cor.test1 <- function(x) {
  x = as.numeric(x)
  n = length(x)
  if (length(which(!is.na(x[1:(n/2)]))) < 4 | length(which(!is.na(x[((n/2)+1):n]))) < 4) 
  {
    c(NA,NA)
  } else 
  {
    test = cor.test(as.numeric(x[c(1:(n/2))]), as.numeric(x[c(((n/2)+1):n)]),na.action = na.omit,method = "pearson")
    c(test$estimate,test$p.value)
  }
}

# define a function to change points to raster
Point2raster = function(points, raster){ #points is vector, raster the template
  nr = nrow(raster)
  nc = ncol(raster)
  ext = extent(raster)
  proj.geo = projection(raster)
  r = raster(matrix(points, nrow = nr, ncol = nc, byrow = TRUE),
             xmn=ext@xmin, xmx=ext@xmax,
             ymn=ext@ymin, ymx=ext@ymax, 
             crs=CRS(proj.geo))
  return(r)
}

#reclassification
recls2 = function(x, threshold = threshold){
#reclassify
x.list = list()
for(layer in 1:nlayers(x))
{
x2 = x[[layer]]
x2.rc = reclassify(x2, threshold)
x.list[[layer]] <- x2.rc
}
x.list = stack(x.list)
return(x.list)
}

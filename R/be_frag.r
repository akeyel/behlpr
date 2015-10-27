#**# Note to self - I think Bastian might be doing something similar - maybe talk to him & share code.

# Code out steps needed for CORINE Data

# Load required packages
library(raster)
library(rgdal)

# Load CORINE Data
# Reference site: http://neondataskills.org/R/Raster-Data-In-R/
corine = raster("C:/docs/beplants/datasets/CORINE_2006/g100_06.tif")

# Load BE Exploratories plot locations
plots.dir = "C:/docs/beplants/datasets/GIS/BExIS_data" #**# NOTE: Trailing slash crashes readOGR. I'd totally write a patch for this.
plots = readOGR(dsn= plots.dir, layer="Grassland_EPs")

# Project plots to UTM
plots.utm = sp::spTransform(plots, CRS("+proj=utm +zone=32 ellps=WGS84"))


#**# Trying to solve spatialpixels problem
#?raster
x = seq(479629,867282)
y = seq(5301768, 5940455)
vals = length(x) * length(y)
z = runif(vals,1,10)

raster(nrows = length(x), ncols = length(y), xmn = min(x), xmx = max(x), ymn = min(y), ymx = max(y), vals = z, crs = CRS("+proj=utm +zone=32 ellps=WGS84"))

test = matrix(c(x,y,z), nrow = 3)
raster(test)

test1 = seq(1,100)
test2 =

test1 = rep(seq(1,100),2)
test2 = sort(rep(seq(1,100),2))
test0 = seq(1,200)
test3 = matrix(c(test1,test2,test0), ncol = 3)

test = raster(test3)
test.sp = as(test, "SpatialPixels")

# Get buffer distances around plots
require(adehabitatMA)

adehabitatMA::buffer(plots.utm, plots.utm, 50)
#Error in adehabitatMA::buffer(plots, plots, 50) : 
#  x should inherit the class SpatialPixels



# Extract Corine data by buffer
#**# This step will overlap with what Bastian is doing.

#http://stackoverflow.com/questions/13982773/crop-for-spatialpolygonsdataframe

# Reclassify Corine data to grassland, forest, & Other


# Estimate % cover for grassland & forest


# Estimate contagion for grassland habitat
#**# Rodolphe had suggested doing this by edge density. I'm not exactly sure how to implement this in the raster context
#**# but if I do, I can add it as an output to spatialdemography at the same time (highly desirable)

# output data to be read in by R


## Helper functions (or is Rodolphe's edge density measure already coded in R?)
#**# Start with internet search - try to avoid recreating the wheel. & talk w/Bastian


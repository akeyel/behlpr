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

# Get buffer distances around plots


# Extract Corine data by buffer


# Reclassify Corine data to grassland, forest, & Other


# Estimate % cover for grassland & forest


# Estimate contagion for grassland habitat
#**# Rodolphe had suggested doing this by edge density. I'm not exactly sure how to implement this in the raster context
#**# but if I do, I can add it as an output to spatialdemography at the same time (highly desirable)

# output data to be read in by R


## Helper functions (or is Rodolphe's edge density measure already coded in R?)
#**# Start with internet search - try to avoid recreating the wheel. & talk w/Bastian


#**# Note to self - I think Bastian might be doing something similar - maybe talk to him & share code.

# Code out steps needed for CORINE Data

# Load required packages
library(raster)

# Load CORINE Data
corine.file = "C:/docs/beplants/datasets/CORINE_2006/g100_06.tif"
#**# Read in as raster here

# LOad BE Exploratories plot locations
#**# Looks like they are in kml format, will need to figure out how to read it.
plots.file = "C:/docs/beplants/datasets/BExIS/google_earth_files/AllPlots.kml"
#**# Read in as points here

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


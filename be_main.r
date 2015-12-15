## Supporting analysis for Keyel et al. in prep: Predicting grassland community
#    changes in response to changes in climate and land management

## Code Overview
# See Analysis_Methods.docx for details on what the analysis is, and why each step is taken.
# The below code is organized as follows:
##  load required packages
##  Set up path information & indicator variables
##  Set up species data
##  Do optional analyses (Sandra, you can skip these)
##  Set up predictor variables
##  Do some basic descriptive analyses
### Run Species Distribution Model
### Run SpatialDemography model

## Set up path information
# Main working directory
spath = "C:/docs/beplants/Scripts/"
setwd(spath)

# Install custom packages (only need to do this 1x, or after updating packages)
install.packages("behlpr", type = "source", repos = NULL)
#install.packages("myr", type = "source", repos = NULL)

## load required packages
library(behlpr)
library(spatialdemography)
library(myr)

# Set up paths
in.path = "C:/docs/beplants/drafts/BE_data_paper/R_inputs/"  # Path for inputs for spatialdemography model
out.path = "C:/docs/beplants/drafts/BE_data_paper/R_outputs/" # Base path for outputs for spatialdemography model
opath = sprintf("%s/analysis/", out.path)  # Path for main analysis results
desc.path = sprintf("%sdescriptive_information/", out.path) # Path to descriptive analysis results
DispPath = sprintf("%sdispersal_tables/", spath) # Path to dispersal tables for spatialdemography model
dpath = "C:/docs/beplants/datasets/" # Base path for data sets
template.dir = sprintf("%sinput_templates/", in.path) # Path to templates to provide the basis for the analyses in SpatialDemography
sp.base = sprintf("%s/Species/spfile_base.csv", template.dir) # Path containing species data to use as a template for input into the SpatialDemography model
spfile = sprintf("%s/Species/spfile_spdem.csv", template.dir) # This file will be generated, the actual species file to be used by the spatialdemography model

# Set up indicators
give.messages = 0 # 0 = no messages, 1 = messages; # Indicator for whether to give messages about possible concerns about the code
optional.analyses = 0 # Do optional analyses
sdm.type = 1 # 1 = simple, 2 = with biomod2, 3 = with previous.occupancy data; # Set up SDM (species distribution model) settings

# Set up spatialdemography settings
file.ending = "BE"  
start.year = 2008 # Start with first year of consistent plant data
end.year = 2013   # End with last year of plant data
pool.type = "plot" # region & all are other anticipated options
landscape.extent = 5 #**# can try different landscape sizes
n.cells = landscape.extent ^ 2 # get number of cells in landscape
focal.cell.position = floor(median(seq(1,n.cells))) # get position of focal cell # median of the sequence gets the middle number. Floor ensures that it is a whole number.
n.sim = 1 # number of simulations
n.plots = 150 # number of exploratory plots

## Species Data
# Reformat species data as needed, then read in species data file
sp.file = sprintf("%s/plants/PlantDataBE.csv", dpath)
sp.lookup.file = sprintf("%s/plants/PlantDataBE_names.csv", dpath)
sp.file.reformatted = sprintf("%s/plants/PlantDataBE_rf.csv", dpath)
# only reformat data file if a reformatted file does not already exist.
if (!file.exists(sp.file.reformatted)){
  sp.data = read.sp.data(sp.file, sp.lookup.file)  
}else{
  sp.data = read.csv(sp.file.reformatted)
}

# Get names of species
sp.names = names(sp.data[2:length(sp.data)]) # Drop plotyear column, but retain all other species names.

# Number of species
n.sp = length(sp.names) 

# Set timestep to use for model output verification
model.timestep = 2013

# Convert the observed data into a format that can be used for model output verification
if (give.messages == 1){ message("Non-optimal approach - convert sp.data to a SpeciesData.csv format to allow me to use an existing function.")}
actual.presences.lst = convert.sp.data(sp.data, model.timestep) # Create the model object based on sp.data


## Do optional analyses
if (optional.analyses == 1){
  # Do metacom analysis
  metacom.analysis(sp.data, 2008) #**# NOT WORKING
  if (give.messages == 1){ metacom.analysis.notes() }

  # Read in trait data
  trait.data = read.traits()  #**# This fails now, I moved the files somewhere else.
}

#This point saved as SavePoint1.RData
if (give.messages == 1){message("SavePoint1.RData is out of date")}

## Set up predictor variables
if (give.messages == 1){message("Note: to use distance to most recently occupied patch, I think we have to restrict the dataset to colonizations")}

# Give a descriptor of the analysis.
analysis.type = "abiotic"

# Abiotic variables
lu.path = sprintf("%sBExIS/landuse/LUI_tool_cleaned.csv", dpath)
frag.path = sprintf("%sfragmentation/frag_data_for_analysis.csv", dpath)
predictors.vec = c("landuse", "frag") #**# Sandra, you will want to change this to a vector with a single element "climate"
predictors.paths = c(lu.path, frag.path) #**# Sandra, you will want to change this to be a single element vector called climate.path

# Actually set up the predictor variables
out.info = setup.predictors(sp.data, predictors.vec, predictors.paths, start.year, end.year)
my.data = out.info[[1]]   # This will be the dataset merging species and your predictor variables
my.index = out.info[[2]]  # This will contain an index identifying properties of the variables in my.data
col.index = out.info[[3]] # This is an index used by the code for selectively removing columns from the analysis.

# Save Point 1 Updated 2015-12-12

## Do basic descriptive analysis of data
if(give.messages == 1){message("Descriptive analysis still needs cleaning & proper outputting of results")}
do.descriptives(my.data, my.index, desc.path)

## Drop NA records from dataset
my.data = na.omit(my.data)
#str(my.data) # 11, 91, 325, 326, 439, 798 - rows dropped from my.data due to missing values

### Run Species Distribution model
sdm.out = do.sdm(my.data, my.index, col.index, sdm.type)
validation.diagnostics.lst = sdm.out[[1]]
thresholds = sdm.out[[2]]
if (give.messages == 1){message("Need a way to make the output of thresholds more general if we start comparing linear models too.")}

# Create table output of results for further analysis
model.type = "SDM-threshold"
model.name = analysis.type
results.table = sprintf("%sMain_Results.csv", opath)
dir.create(opath, recursive = T)
my.sdm.results = output.results.to.table(model.type, model.name, results.table, sp.names, validation.diagnostics.lst) 

#**# NEED TO FIX my.sdm.results for plausibility name to plausible

# This location saved as SavePoint2.RData
#if(give.messages == 1){message("SavePoint2.RData is currently out of date")} #Updated 2015-12-12

# Visualize SDM analysis (show landscape with SDM results)
#visualize.sdm.main(sdm.out)

# Visually examine quality of sdm results for all species (via plots)
check.sdm.quality.v2(my.sdm.results) #**# Add pdf here?



### Run SpatialDemography  
# Set up species file
prepare.spfile(sp.base, spfile, thresholds)

# Restrict data set to starting year
spdem.out = spdem.prep(my.data, start.year, col.index)
sp.dat = spdem.out[[1]]
land.dat = spdem.out[[2]]

#**# LEFT OFF HERE IN UPDATING CODE FOR SANDRA
#**# Needs to be flexible to take data sets other than LUI.
#**# NOTE: Sandra can start with the SDM part of the code, and come to this later,
#so you don't actually need to have this by Sunday.

# Create a list to hold results from ALL simulations
sim.lst = list()
# loop through simulations (note: for pool.type = plot, it should be fairly deterministic?)
for (k in 1:n.sim){
  k = 1 
  predicted.presences.lst = rep(list(NA), n.sp)
  # Loop through plots in the exploratories
  for (p in 1:n.plots){
    p = 1
    # Restrict data set to a specific plot
    if (nrow(sp.dat) != n.plots) { stop("Species dataset does not match plot dataset. Something is wrong.")}
    if (nrow(land.dat) != n.plots) {stop("Landscape data dataset does not match plot dataset. Something is wrong.")}
    plot = sp.dat$plots[p]
    sp.plot.dat = sp.dat[sp.dat$plots == plot, ]
    land.plot.dat = land.dat[land.dat$plots == plot, ]
    # Set up landscape configuration (hacky way for now, will come back to this)
    # Currently using fixed values which may or may not correpsond to the focal cell. This stupid, because we need the value for the focal cell. In addition to species pool, we need an environmental pool!). Maybe I should just use LUI and not worry about mowing, grazing, and fertilization to start with??
    # customize to focal plot 
    fert.plot = land.plot.dat[[2]]
    mowing.plot = land.plot.dat[[3]]
    grazing.plot = land.plot.dat[[4]]
    lui.plot = land.plot.dat[[5]]
    
    lui = rep(lui.plot, n.cells) 
    fert = rep(fert.plot, n.cells) 
    mowing = rep(mowing.plot, n.cells) 
    grazing = rep(grazing.plot, n.cells)
    landscape = list(lui, fert, mowing, grazing) #**# make this correspond to an input for spatialdemography
    
    # Create landscape layer files (can you just input landscape objects to spatialdemography? I think so???)
    lpath = landscape # Hidden functionality that needs better documentation: this can also be the landscape itself, if properly formatted!
    
    #Run name & run path
    run.path = sprintf("%splot_%s_sim_%s/", opath, p, k)
    run.name = sprintf("plot_%s_sim_%s", p, k)
    
    # Set up locations file
    sp.locs = sprintf("%s/locations_%s_%s.csv", template.dir, pool.type, Digit(n.sim, 4))
    prepare.sp.locations(sp.locs, sp.plot.dat, pool.type, n.cells) # Create sp.locs file
    
    scale.vec = "landscape" # not really interested in this - the actual model results are irrelevant for my purposes.
    
    # Run main model
    spdem.out = do.spatial.demography(run.name, run.path, spfile, sp.locs,template.dir, file.ending, lpath, DisPath, opath, scale.vec)    
    
    # Extract result from focal cell & append to list
    output.file = sprintf("%sLanduse/SpeciesData.csv", run.path)
    validation.timestep = 2 #**# Don't run it for 6 timesteps then!
    predicted.presences.lst = extract.cell.values(predicted.presences.lst, output.file, focal.cell.position, validation.timestep)
  }
  
  # Compare model-produced results to observed data
  validation.diagnostics.lst = evaluate.spatialdemography.v2(predicted.presences.lst, actual.presences.lst, n.sp)
  
  # Add this runs accuracy stats to the running total.
  sim.lst = append(sim.lst, list(validation.diagnostics.lst))
  
}

# Convert validation.diagnostics from multiple simulations into an average result & output the average result
validation.diagnostics.lst = summarize.validation(sim.lst, n.sp)

# Write average accuracy to file
model.type = "metacommunity-agg-mean"
model.name = analysis
my.spdem.results = output.results.to.table(model.type, model.name, results.table, sp.names, validation.diagnostics.lst)

# Examine results table
combined.results = read.csv(results.table)
combined.results

#} # End bracket if code is being run in a for loop for multiple analysis sets.


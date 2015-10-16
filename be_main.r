## Supporting analysis for Keyel et al. in prep: "TITLE GOES HERE"

## Index to data analyses
# One analysis with just abiotic factors
# One analysis with abiotic factors & proxies for dispersal and biotic filters.

## load required packages
library(behlpr)
library(spatialdemography)
library(myr)

## 

# Set up path information
spath = "C:/docs/beplants/Scripts/"
in.path = "C:/docs/beplants/drafts/BE_data_paper/R_inputs/"
out.path = "C:/docs/beplants/drafts/BE_data_paper/R_outputs/"
DispPath = sprintf("%sdispersal_tables/", spath) # Path to dispersal tables
dpath = "C:/docs/beplants/datasets/"
opath = sprintf("%s/proof_of_concept/", out.path)
desc.path = sprintf("%sdescriptive_information/", out.path)

#Set up paths for SpatialDemography
template.dir = sprintf("%sinput_templates/", in.path) # Path to templates to provide the basis for the analyses
sp.base = sprintf("%s/Species/spfile_base.csv", template.dir) 
spfile = sprintf("%s/Species/spfile_spdem.csv", template.dir) # This file will be generated

setwd(spath)

# Set up inputs
# Species Data
sp.file = sprintf("%s/combined/SpeciesDataBE.csv", dpath)
sp.lookup.file = sprintf("%s/combined/species_lookup.csv", dpath)
sp.data = read.sp.data(sp.file, sp.lookup.file)
sp.names = names(sp.data[2:length(sp.data)]) # Drop plotyear column, but retain all other species names.

# Do metacom analysis
message("Metacom analysis code in progress!")
metacom.analysis(sp.data)

# Read in trait data
trait.data = read.traits()

#**# Consider dropping rare species (I dealt with this in the descriptives, too, so need to figure out if, where, and when to code this.)

# Subset species data based on chosen criteria for trait data
#link = etwas
#trait.lst = etwas
#trait.vals.lst = etwas
#sp.data = subset.by.traits(sp.data, trait.data, link, trait.lst, trait.vals.lst)

# Set up SDM settings
sdm.type = 1 # 1 = simple, 2 = with biomod2, 3 = with previous data

# Set up spatialdemography settings
file.ending = "BE"  
start.year = 2009 #**# Arbitrarily chosen starting year
pool.type = "plot" # region & all are other anticipated options
landscape.extent = 5 #**# can try different landscape sizes
n.cells = landscape.extent ^ 2 # get number of cells in landscape
focal.cell.position = floor(median(seq(1,n.cells))) # get position of focal cell # median of the sequence gets the middle number. Floor ensures that it is a whole number.
n.sim = 1
n.plots = 150
n.sp = length(sp.names)

# Set up validation data for spatialdemography model
#**# Non-optimal approach - convert my.data to a SpeciesData.csv format to allow me to use an existing function.
model.timestep = 2010
actual.presences.lst = convert.sp.data(sp.data, model.timestep) # Create the model object based on sp.data

# Set up predictor variables
analysis.vec = c("abiotic", "all")

#This point saved as SavePoint1.RData

# Run once with only abiotic factors, then a second time with biotic factors
#for (analysis.type in analysis.vec){
  analysis.type = analysis.vec[1]
  # Potential predictor variables: c("climate", "landuse", "biotic", "isolation")
  
  # Abiotic variables
  lu.path = sprintf("%sBExIS/landuse/16029.csv", dpath)
  config.path = stop("This variable is not yet available")
  #abiotic.vec = c("landuse", "configuration")
  abiotic.vec = c("landuse")
  #abiotic.paths = c(lu.path, config.path)
  abiotic.paths= c(lu.path)
  
  # Biotic variables
  biotic.vec = biotic.paths = c() # Create empty variables in case no biotic variables
  if (analysis.type == "all"){
    iso.path = stop("This variable has not yet been configured")
    prev.path = stop("This variable has not yet been configured")
    biotic.path = stop("This variable has not yet been configured")
    biotic.vec = c("isolation", "prior occupancy", "biotic filter")
    biotic.paths = c(iso.path, prev.path, biotic.path)
  }
  
  # Set up predictors
  predictors.vec = c(abiotic.vec, biotic.vec)
  predictors.paths = c(abiotic.paths, biotic.paths)
  
  out.info = setup.predictors(sp.data, predictors.vec, predictors.paths)
  my.data = out.info[[1]]
  my.index = out.info[[2]]
  col.index = out.info[[3]]
      
  ### Do basic descriptive analysis of data
  message("Descriptive analysis still needs cleaning & proper outputting of results")
  do.descriptives(my.data, my.index, desc.path)
  
  ### Run Species Distribution model
  sdm.out = do.sdm(my.data, my.index, col.index, sdm.type)
  validation.diagnostics.lst = sdm.out[[1]]
  thresholds = sdm.out[[2]] #**# Need a way to make this more general if we start comparing linear models too.

  # Visually examine quality of sdm results for all species
  #check.sdm.quality(validation.diagnostics.lst, "Land Use")
  
  # Visualize SDM analysis
  #visualize.sdm.main(sdm.out)
  
  # Create table output of results for further analysis
  model.type = "SDM-threshold"
  model.name = analysis.type
  results.table = sprintf("%sMain_Results.csv", opath)
  my.sdm.results = output.results.to.table(model.type, model.name, results.table, sp.names, validation.diagnostics.lst) 
  
  # This location saved as SavePoint2.RData

  ### Run SpatialDemography  
  # Set up species file
  prepare.spfile(sp.base, spfile, thresholds)
    
  # Restrict data set to starting year #**# Note: species pools should be set up prior to this step
  spdem.out = spdem.prep(my.data, start.year, col.index) # Initialize a new data object (should I just use my.data?)
  sp.dat = spdem.out[[1]]
  land.dat = spdem.out[[2]]
      
  # Create a list to hold results from ALL simulations
  sim.lst = list()
  # loop through simulations (note: for pool.type = plot, it should be fairly deterministic?)
  #for (k in 1:n.sim){
   k = 1 
    predicted.presences.lst = rep(list(NA), n.sp)
    # Loop through plots in the exploratories
    #for (p in 1:n.plots){
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
  
}


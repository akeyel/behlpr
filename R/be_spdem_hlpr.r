# This script contains functions related to the implementation of
# the spatialdemography analyses

#' Main function for generating the spatial demography outputs
#' 
#' 
#' 
#' @param inputs.ending If the input files need a different ending than the output files, this allows that to happen.
#' @export do.spatial.demography
do.spatial.demography = function(run.name, run.path, spfile, sp.locs,template.dir,
                                 file.ending, lpath, DisPath, opath, scale.vec, inputs.ending = "none"){
  
  # Check if run.path exists, if not, create it (should you add this to spatialdemography instead??? Maybe with an option to turn it on?)
  
  # Set up path to inputs (if not specified)
  if (inputs.ending == "none"){
    inputs.ending = file.ending
  }
  
  # Get cells for which to calculate diversity metrics
  scale.cells.lst = get.scale.lst(scale.vec)
  
  ## Add in data for environmental layers & environmental layer change (separate concepts, single file for now)
  
  env.file = sprintf("%sEnvironmental_layers_%s.csv", template.dir, inputs.ending)
  # landcover, mowing, temperature. No change
  
  #Set up initial conditions and settings. Using basic settings for simplicity
  initial.conditions.file = sprintf("%sInitial_conditions_%s.csv", template.dir, inputs.ending)
  settings.file = sprintf("%sSettings_%s.csv", template.dir, inputs.ending)
  
  ## The species file here is generated external to spatialdemography, so these options are not used
  sp.instr.file = sp.resp.instr.file = NA
    
  ## Create a path to store results
  ResultsFile = sprintf("%sResults%s.csv",opath, run.name) # Create a file to contain results information & write the header for the file
  if (file.exists(ResultsFile))  file.remove(ResultsFile)    # Delete any existing results file
  
  scn = 1
  s.lbl = sprintf("scn%s", scn)
    
  #**# Should be irrelevant, the entire code ought to be deterministic at this point
  #my.seed = 987654 + sc #Arbitrary seed, but get a different set of random numbers for each scenario and make it so that running scenarios out of order still gives the same seed value
  #set.seed(my.seed)
    
  #Run main model
  results = SpatialDemography(scn,s.lbl,file.ending,DispPath,run.path,opath,ResultsFile,initial.conditions.file,settings.file,env.file,spfile, sp.instr.file,sp.resp.instr.file, landscape.dir = lpath, locations.file = sp.locs, scale.vec = scale.vec, scale.cells.lst = scale.cells.lst)
  
  # Return outputs
  return(results)
}

#' Set up scale.cells.lst for the biodiveristy exploratories
#' 
get.scale.lst = function(scale.vec){

  scale.cells.lst = list()
  
  if (scale.vec[1] == "landscape" & length(scale.vec) == 1){ scale.cells.lst = "default" }
  
  for (i in 1:length(scale.vec)){
    expl = scale.vec[i]
    #scale.cells.lst = list(c(T,F,F,T)) #T indicates that the cell values should be included in the calculations
    #Set up subsets for all and for each exploratory (the all subset is needed because the model requires a square landscape, so some irrelevant cells have been added for this analysis.
    if (expl == "all"){ scale.cells.lst = append(scale.cells.lst, list(rep(T,150),rep(F,19)))           }
    if (expl == "alb"){ scale.cells.lst = append(scale.cells.lst, list(rep(T,50),rep(F,119)))           }
    if (expl == "hai"){ scale.cells.lst = append(scale.cells.lst, list(rep(F,50),rep(T,50),rep(F,69)))  }
    if (expl == "sch"){ scale.cells.lst = append(scale.cells.lst, list(rep(F,100),rep(T,50),rep(F,19))) }
  }
  
  return(scale.cells.lst)
}

#' Extract data for a species distribution model from spatialdemography model
#' 
#' I think I wrote a function doing exactly this within spatialdemography, to produce species locations files?
#' Look to see if this is redundant code!
#' 
#' @param data.file This input should be in the format of a species locations file from spatialdemography, or in the format of the SpeciesData.csv file created as a model output
#' @param timestep Will allow a user to restrict the data extraction to a particular timestep.
#' 
spdem.extract = function(data.file, timestep = "all", adults.only = TRUE){
  
  require(myr)
  
  if (timestep != "all"){ stop("Functionality for the extraction of a specific timestep hasn't been coded yet. Sorry for the inconvenience") }
  
  if (adults.only != TRUE){ stop("Currently the code can only extract adult presences. Merge this function with existing code in spatialdemography for added options.")}
  
  raw.data = read.csv(data.file)
  
  if (adults.only == TRUE){
    raw.data = raw.data[raw.data$LifeStage == "Adults", ]
  }
  
  # Reformat to have plotyear
  plots = names(raw.data)[4:length(raw.data)] #Cols 1 - 3 are not plots
  timesteps = unique(raw.data$TimeStep)
  timesteps = sapply(timesteps, Digit, 4)
  placeholders = rep(NA, (length(plots) * length(timesteps))) # Calculate the size of the dataframe to be produced
  plots.expanded = rep(plots, length(timesteps))
  timesteps.expanded = sort(rep(timesteps, length(plots)))
  plotyears = sprintf("%s_%s", plots.expanded, timesteps.expanded)
  my.data = data.frame(plotyear = plotyears)
  
  # Add each species to the dataframe
  species = unique(raw.data$Species)
  sp.len = ceiling(length(species) / 10 ) # Get number of digits to encompass all species
  for (sp in 1:length(species)){
    this.sp = sprintf("sp_%s", Digit(species[sp], sp.len))
    my.data[[this.sp]] = placeholders
  }
  
  for (i in 1:nrow(raw.data)){
    for (j in 1:length(plots)){
      plot = plots[j]
      this.rec = raw.data[i , plot]
      this.sp.num = raw.data[i, "Species"]
      this.sp = sprintf("sp_%s", Digit(this.sp.num, sp.len))
      this.time = raw.data[i, "TimeStep"]
      plotyear = sprintf("%s_%s", plot, Digit(this.time, 4))
      # Update dataframe with correct value
      my.data[my.data$plotyear == plotyear, this.sp] = this.rec
    }
  }
  return(my.data)
}

#' Prepare species file
#' 
#' Add additional response traits to an existing species file containing base
#' vital rates, dispersal parameters, and (any) known response traits.
#' 
#' @param spbase The base file to append response traits to
#' @param spfile The new species file to be created
#' @param thresholds Additional thresholds to be appended to the species file
#' 
#' @details thresholds should be a set of nested lists.
#' The first list should have one entry for each species
#' Each species list entry should be a list of 3 vectors
#' vector 1 gives variable names, vector 2 lower thresholds, vector 3 = upper thresholds.
#' @export prepare.spfile
prepare.spfile = function(spbase, spfile, thresholds){
  
  # Update the file based on the observed thresholds
  
  # Loop through species (Organizing by variable would make way more sense.)
  n.sp = length(thresholds)
  for (i in 1:n.sp){
    # 
    variables = thresholds[[i]]
    
    var.names = variables[[1]]    
    low.thresh = variables[[2]]
    up.thresh = variables[[3]]
    
    # Create check to watch for duplicate first letters
    check.var.names(var.names)
    
    # Create dataframe of response traits
    if (i == 1){
      # Create one row for each species
      response.traits = data.frame(index = seq(1:n.sp))
    }
    
    for (j in 1:length(var.names)){
      # Extract information for this variable
      this.var = var.names[j]
      this.low = low.thresh[j]
      this.up = up.thresh[j]

      var.abbr = toupper(substr(this.var,1,1))
      #**# For now, assuming that p23 is the response trait of interest.
      # Add to list of dubious assumptions that I can fix in the future!
      response.trait = sprintf("%s.p23", var.abbr) 
      
      # For first species, create a new object for the variable
      if (i == 1){
        # Create a variable abbreviation (first letter, uppercase)
        # Add response trait to response.traits dataframe
        response.traits[[response.trait]] = rep(NA, n.sp)
      }
      # update this particular response trait value
      this.response = sprintf("NA;106;%s;%s",this.low, this.up)
      response.traits[[response.trait]][i] = this.response
    }  
  }
  
  # Read in the base file
  sp.data = read.csv(spbase)
  
  # Check that number of rows in sp.data matches number of species for which there are thresholds. If not, drop excess species.
  if (nrow(sp.data) != n.sp){
    warning("Number of species in base file exceeds number of species for which thresholds are given. Excess species will be automatically dropped. This step should be checked carefully for errors!!!! Or the base file should be revised to avoid this issue.")
    sp.data = sp.data[1:n.sp, ]
  }
  
  # Add new species traits to existing species data
  sp.data = cbind(sp.data, response.traits)
  
  # Write finished species file (should I also the species file as an R object?)
  if (file.exists(spfile)){
    file.remove(spfile)
  }
  write.table(sp.data, file = spfile, row.names = F, col.names = T, sep = ",")
}

#' Check variable names for spatialdemography compliance
#' 
#' Check that each variable begins with a different first letter.
#' @param var.names the list of variable names to be checked
check.var.names = function(var.names){
    
  var.firsts = c()
  for (k in 1:length(var.names)){
    var.first = substr(var.names[k],1,1)
    var.first = toupper(var.first)
    var.firsts = c(var.firsts, var.first)
  }
  
  if (length(var.names) != length(unique(var.firsts))){
    warning("Only the response trait p23 is currently being set up")
    stop("Variable names must start with unique letters, otherwise SpatialDemography will have issues")
  }
  
}

#' Evaluate results of spatialdemography model runs
#'
#' Goal is to compare output of spatialdemography model to a validation dataset
#' 
#' @param validation.file A file containing the Validation data set
#' @param model.file A file containing the data set from model output
#' @param n.sp The number of species (probably could be extracted from the datasets)
#' 
#' @export evaluate.spatialdemography
#' 
evaluate.spatialdemography = function(validation.file, model.file, validation.timestep, model.timestep, n.sp){
  
  # Read in validation data & subset to presence/absence from specified timestep
  v.dat = prep.data(validation.file, validation.timestep)
  
  # Read in model data & subset to presence/absence from specified timestep
  model.dat = prep.data(model.file, model.timestep)
  
  require(myr)
  
  validation.diagnostics.lst = list()
  
  for (sp in 1:n.sp){
    #Subset data to particular species in question
    v.dat.sp = v.dat[v.dat$Species == sp, ]
    model.dat.sp = model.dat[model.dat$Species == sp, ]
    
    # Drop species column
    v.dat.sp$Species = NULL
    model.dat.sp$Species = NULL
    
    # Recode to presence/absence
    v.dat.presences = sapply(v.dat.sp, pa.recode)
    model.dat.presences = sapply(model.dat.sp, pa.recode)
   
    #validate.sdm comes from be_sdm_hlpr.r and could be renamed & moved to be_hlpr.r
    validation.diagnostics = validate.sdm(v.dat.presences, model.dat.presences)
    #sensitivity = validation.diagnostics[[1]]
    #specificity = validation.diagnostics[[2]]
    #accuracy = validation.diagnostics[[3]]
    #true.positive = validation.diagnostics[[4]]
    #false.positive = validation.diagnostics[[5]] 
    #true.negative = validation.diagnostics[[6]]
    #false.negative = validation.diagnostics[[7]]
    #positive.predictive.value = validation.diagnostics[[8]]
    #negative.predictive.value = validation.diagnostics[[9]]
    
    # Append to list
    validation.diagnostics.lst = append(validation.diagnostics.lst, list(validation.diagnostics))
  }
  
  return(validation.diagnostics.lst)  
}

#' Evaluate results of spatialdemography model runs using list objects as inputs.
#'
#' Goal is to compare output of spatialdemography model to a validation dataset
#' 
#' @param estimated.presences.lst A list with one entry for every species, giving predicted occupancy for each exploratory plot
#' @param actual.presences.lst A list with one entry for every species, giving observed occupancy for each exploratory plot
#' @param n.sp The number of species (equal to length of the two above lists!)
#' 
#' @export evaluate.spatialdemography.v2
#' 
evaluate.spatialdemography.v2 = function(estimated.presences.lst, actual.presences.lst, n.sp){
  
  # Check that all lists are correct
  if (n.sp != length(estimated.presences.lst) & n.sp != length(actual.presences.lst)){
    stop("Number of species must equal length of estimated presences list and length of actual presences list")
  }

  validation.diagnostics.lst = list()
  for (sp in 1:n.sp){
    
    estimated.presences = estimated.presences.lst[[sp]]
    actual.presences = actual.presences.lst[[sp]]
    
    #validate.sdm comes from be_sdm_hlpr.r and could be renamed & moved to be_hlpr.r
    validation.diagnostics = validate.sdm(actual.presences, estimated.presences)
    #sensitivity = validation.diagnostics[[1]]
    #specificity = validation.diagnostics[[2]]
    #accuracy = validation.diagnostics[[3]]
    #true.positive = validation.diagnostics[[4]]
    #false.positive = validation.diagnostics[[5]] 
    #true.negative = validation.diagnostics[[6]]
    #false.negative = validation.diagnostics[[7]]
    #positive.predictive.value = validation.diagnostics[[8]]
    #negative.predictive.value = validation.diagnostics[[9]]
    
    # Append to list
    validation.diagnostics.lst = append(validation.diagnostics.lst, list(validation.diagnostics))
  }
  
  return(validation.diagnostics.lst)
}

#' Create a summary from multiple simulations
#' 
#' Think about whether there is a better way to code this!
#' Add options for median & quantiles later
#' @export summarize.validation
#' 
summarize.validation = function(sim.lst, n.sp){ # summary.type,
    
  summary.type = "mean" #**# currently hardwired, can change this
  n.sim = length(sim.lst)
  
  if (n.sim > 1){
    message("This function needs testing/checking for more than one simulation")  
  }
  
  # Set up for calculation of mean
  sp.sen = sp.spec = sp.acc = sp.tp = sp.fp = sp.tn = sp.fn = 0
  vd.objects = c(sp.sen, sp.spec, sp.acc, sp.tp, sp.fp, sp.tn, sp.fn)
  # Create a list to contain results
  out.lst = rep(list(vd.objects), n.sp)
  
  #**# Follow up on this later if you actually want to use it!
  #sp.sen.min = sp.spec.min = sp.acc.min = 0
  #sp.sen.max = sp.spec.max = sp.spec.max = 0
  
  # Loop through simulations
  for (i in 1:length(sim.lst)){
    vd.lst = sim.lst[[i]]
    if (length(vd.lst) != n.sp){ stop("number of species must match the length of the validation diagnostics lists") }
    
    # Loop through species
    for (sp in 1:n.sp){
      # Pull out the current values from validation.diagnostics.lst
      validation.diagnostics = vd.lst[[sp]]
      sensitivity = validation.diagnostics[[1]]
      specificity = validation.diagnostics[[2]]
      accuracy = validation.diagnostics[[3]]
      true.positive = validation.diagnostics[[4]]
      false.positive = validation.diagnostics[[5]] 
      true.negative = validation.diagnostics[[6]]
      false.negative = validation.diagnostics[[7]]

      # Create vector of current values
      cur.vals = c(sensitivity, specificity, accuracy, true.positive, false.positive, true.negative, false.negative)
      
      # Extract previous values
      old.vals = out.lst[[sp]]
      
      # Add value to running total & update main list
      new.vals = old.vals + cur.vals #**# How to handle NA's?? 
      out.lst[[sp]] = new.vals
    }
  }
  
  # Loop through species one last time & divide by sample size
  for (sp in 1:n.sp){
    total.vals = out.lst[[sp]]
    mean.vals = total.vals / n.sim
    mv = mean.vals # Make it shorter
    
    
    if (summary.type == "mean"){
      # Reformat to be in list format needed by downstream code
      out.stuff = list(sensitivity = mv[1], specificity = mv[2], accuracy = mv[3],true.positive = mv[4], false.positive = mv[5], true.negative = mv[6], false.negative = mv[7])
      out.lst[[sp]] = out.stuff  
    }     
  }
  
  return(out.lst)
}

#' prepare data for validation
#' 
#' 
prep.data = function(in.file, timestep){
  my.dat = read.csv(in.file)
  # Subset to adults only
  my.dat = my.dat[my.dat$LifeStage == "Adults", ]
  my.dat = my.dat[my.dat$TimeStep == timestep, ]
  my.dat$LifeStage = NULL
  my.dat$TimeStep = NULL
  
  return(my.dat)
}

#' prepare data set for spatialdemography model
#' 
#' @export spdem.prep
spdem.prep = function(my.data, start.year, col.index){

  # Add year column
  plotyears = my.data$plotyear
  years = sapply(plotyears, substr, 7, 10)
  years = as.num(years)
  my.data$year = years
  
  # Add plots column
  plots = sapply(plotyears, substr, 1, 5)
  my.data$plots = plots
  
  # Subset my.data to just the chosen starting year
  my.data = my.data[my.data$year == start.year, ]
  
  # update col.index to include additional year and plot columns (labeled -1)
  col.index = c(col.index, -1, -1) #**# These columns will be dropped in the next step
  # Get plot as an object to join back after subsetting
  plots = my.data$plots # create new plots variable after subsetting

  # Get data on landscape parameters
  land.dat = matrix.subset(my.data, col.index, 1) #**# May need to be >0 if more than 1 is used as an index.
  land.dat = cbind(plots, land.dat)
  land.dat$plotyear = NULL
  
  # use col.index to drop non-species columns
  sp.dat = matrix.subset(my.data, col.index, 0)
  sp.dat = cbind(plots, sp.dat) #Would the plots be more useful as row names?
  sp.dat$plotyear = NULL
  
  return(list(sp.dat, land.dat))
  
}

#' Create species locations file
#'
#' Given a landscape, create initial species locations. 
#' 
#' Note to self: think about whether you want species initialization to depend on landscape cells
#' I.e. to what extent are species drawn randomly, and to what extent are they drawn based on the
#' environmental conditions present? (we could test some interesting hypotheses with this model
#' setup.)
#' 
#' @export prepare.sp.locations
#' 
prepare.sp.locations = function(sp.locs, sp.plot.dat, pool.type, n.cells){

  # set up locations in landscape
  # If pool.type == plot, just make every cell in the landscape exactly like the plot
  if (pool.type == "plot"){
    cell.vals = sp.plot.dat[2:length(sp.plot.dat)] #Drop plot (unless you drop it sooner?)
    cell.vals = matrix(cell.vals, ncol = 1)
    loc.mat = cell.vals
    
    # Go through cells in landscape, and add the columns to the species locations matrix
    for (i in 2:n.cells){
      loc.mat = cbind(loc.mat, cell.vals)
    }
    
    #Add header to matrix
    colnames(loc.mat) = sprintf("Cell%s", seq(1,n.cells))
    
    #Add header columns to matrix
    species = names(sp.plot.dat)[2:length(sp.plot.dat)]
    n.sp = length(species)
    Species = seq(1,n.sp) #Convert species to number, without a lookup table. Need to deal with this!
    LifeStage = rep("Adults", n.sp)
    TimeStep = rep(1, n.sp)
    loc.mat = cbind(LifeStage, Species, TimeStep, loc.mat)
  }
    
  # Replace focal cell info with focal cell data #**# for other pool types
  
  # Write data to file (are you sure you don't want to just pass it as an R object?)
    #**# I think the flatfile approach you are using is leading to a lot of confusion and frustration
    #**# You might want to try to switch to a more "pure R" approach. (like the reviewers suggested, I might add)
  write.table(loc.mat, file = sp.locs, col.names = T, row.names = F, append = F, sep = ",")
  
}

#' Set up species locations
#' 
#' 
setup.locs = function(this.plotyear, my.data, pool.type){
  
  # Check that predictor variables have already been dropped (otherwise do this before merging them in??)
  
  # Select data just from this plot
  fc.sp = my.data[plotyear == this.plotyear, ]
  #**# LEFT OFF HERE - NEED TO EXTRACT A CELL's WORTH OF INFOMRATION IN SPECIES.LOCS format
  # WIll also need to add leading columns at some point!
  
  
  return(cell.vals)
}

#' Create species locations file
#' 
#' This approach was designed to convert biodiversity exploratories data into
#' the spatialdemography model. Problem is, plots are not adjacent to one another
#' so that really isn't an appropriate analysis technique for this dataset.
#' 
prepare.sp.locations.old = function(sp.locs, my.data, start.year, col.index){
  
  # use col.index to drop non-species columns
  my.data = matrix.subset(my.data, col.index, 0)
  
  # Add year column
  plotyears = my.data$plotyear
  years = sapply(plotyears, substr, 7, 10)
  years = as.num(years)
  my.data$year = years

  # Add plots column
  plots = sapply(plotyears, substr, 1, 5)
  my.data$plots = plots
  
  
  # Subset my.data to just the chosen starting year
  my.data = my.data[my.data$year == start.year, ]
  
  #Drop unnecessary columns
  my.data$plotyear = NULL
  my.data$year = NULL  
  
  #Convert plots to cells #**# I don't like this step
  
  cells = sprintf("Cell%s", seq(1, length(my.data$plots)))
  my.data$cell = cells
  my.data$plots = NULL #**# Drop plots column without making a lookup file for re-matching later. (bad idea)
  
  # Convert species names to numeric.
  sp.vec = names(my.data[1:(length(my.data) - 1)])
  n.sp = length(sp.vec)
  
  # Create dataframe with required columns (LifeStage, Species, TimeStep)  
  loc.data = data.frame(LifeStage = rep("Adults", n.sp), Species = seq(1,n.sp), TimeStep = 1)
  
  # add an entry for each plot
  for (i in 1:nrow(my.data)){
    cell = my.data$cell[i]
    loc.data[[cell]] = rep(NA, n.sp)
    #plot = my.data$plots[i]
    #loc.data[[plot]] = rep(NA, n.sp)
  }
  
  # Go through rows & columns of my.data & put into loc.data data.frame
  for (i in 1:nrow(my.data)){
    for (j in 1:(ncol(my.data) - 1)){
      #plot = my.data$plots[i]
      cell = my.data$cell[i]
      
      # Get the value corresponding to a species and plot
      this.val = my.data[i, j] 
      
      # Add it to the dataframe
      loc.data[[cell]][j] = this.val
    }
  }
  
  # Write locations file
  write.table(loc.data, file = sp.locs, col.names = T, row.names = F, append = F, sep = ",")
}

#' Junk - alternative approach that I may come back to
#' 
my.junk = function(NOTHING){
  # Create a list object to hold results for each species
  sp.results.lst = rep(list(NA), n.sp)
  
  # Loop through species & calculate for each species
  for (sp in 1:n.sp){
    #**# Kind of a quick and dirty approach here - something more sophisticated is desirable.
    # Calculate percent of simulation runs where each species is present
    percents.vec = calculate.percent.present(STUFF)
    
    # Get the number of presences from reality to compare to estimates
    actual.vec = extract.actual(my.data, 2010) #**# 6 years is too long - this is one year later! #**# need to adjust initial conditions file.
    
    # Compare percent predictions to actual observations
    # One number approach that I made up (that someone may have made up before me??)
    
    # Calculate expected number of presences based on probabilities
    est.prevalence = sum(percents.vec) # Treat probabilities as partial presences, and sum predicted probabilities. This completely ignores the discrete nature of binary data, and is probably a mathematical "no-no", but it has an intuitive appeal to me. Problem is there is no associated confidence interval to this approach. So maybe some sort of sampling of the probabilities to get a range of prevelance estimates would be more accurate??? Except with a deterministic model, that won't help much! Alternative would be to draw distributions based on the probabilities, and evaluate those (in aggregate). That actually makes a bit more sense.
    
    # Calculate expected distribution based on simulation approach based on probabilities
    probs.out = probability.draws(percents.vec)
    #min, 25%, 50% 75% max? Add 5% in there?
    
    
    # Calculate actual number of presences based on data
    total = sum(actual.vec) # Get total number of presences
    perc.prevalence = total / length(actual.vec) # Scale by total number of plots
    
  }
}

#' Extract actual results
#' 
#' @export extract.actual
extract.actual = function(val.data, year, col.index){
  
  # Drop non-species columns
  val.data = matrix.subset(val.data, col.index, 0)
  
  # Convert val.data to have year column
  plotyears = val.data$plotyear
  years = sapply(plotyears, substr,7,10)
  val.data$year = years
  
  # Restrict val.data to selected timestep
  ext.data = val.data[val.data$year == year, ]
  
  # Drop unnecessary columns
  plotyear = NULL
  year = NULL
  
  pres.vec = sapply(val.vec, pa.recode)
  return(pres.vec)
}

#' Extract the value of a single cell from the SpeciesData.csv
#' 
#' @export extract.cell.values
extract.cell.values = function(predicted.presences.lst, output.file, focal.cell.position, validation.timestep){
  
  #Identify the column to extract
  target.col = focal.cell.position + 3 # File has LifeStage, SPecies, and TimeStep columns
  
  # Extract result from focal cell & append to list
  all.dat = read.csv(output.file)
  
  # Subset all.dat to just adults
  all.dat = all.dat[all.dat$LifeStage == "Adults", ]
  
  # Subset all.dat to just the timestep of interest
  all.dat = all.dat[all.dat$TimeStep == validation.timestep, ]
  
  # Loop through species and extract values
  for (i in 1:length(predicted.presences.lst)){
    # Get value of focal cell
    rec.val = all.dat[all.dat$Species == i, target.col]
    
    # Convert to presences/absence
    rec.val = pa.recode(rec.val)
    
    old.val = predicted.presences.lst[[i]]
    if (is.na(old.val[1])){
      new.val = rec.val
    }else{
      new.val = c(old.val, rec.val)
    }
    # Update predicted.presences.lst
    predicted.presences.lst[[i]] = new.val
  }

  return(predicted.presences.lst)
}


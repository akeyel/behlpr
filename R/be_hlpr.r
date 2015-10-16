### R Code in support of:

## Predicting grassland community changes in response to changes in climate and
# land management
# Keyel et al. in prep

## Supporting functions for be_overview.r

#Author: A.C. Keyel <skeyel@gmail.com>

# Comments:
# This package is not CRAN-worthy.
# - unstated dependencies
# - call functions without identifying the external package they are called from
# - arguments present without documentation
# - some oddities in functions that never became legitimate functions

#' Check plausibility of model
#' 
check.plausibility = function(accuracy, sensitivity){
  plaus = 0
  
  # Currently 0.85 is an arbitrary cutoff roughly corresponding to the threshold used
  # in AIC (except that's a p threshold, and this isn't a probability, well sort of it.
  if (!is.na(accuracy) & !is.na(sensitivity)){
    if (accuracy > 0.85 & sensitivity > 0.85){
      plaus = 1
    }
  }
  return(plaus)
}

#' Output Results to Table
#' 
#' Purpose is to output results in table format, appending to any existing table, unless overwrite is set to TRUE.
#' @export output.results.to.table
output.results.to.table = function(model.type, model.name, results.table,sp.names, validation.diagnostics.lst, do.append = TRUE){
  
  # Format validation diagnostics for output
  # "False Positive Rate", "False Negative Rate", "Positive predictive value", "Negative predictive value",
  out.hdr = c("Model Type", "Model", "Species", "Plausible Model", "Accuracy", "Sensitivity", "Specificity", "True Presences", "False Presences", "True Absences", "False Absences")
  
  # Extract results for first species
  vd = validation.diagnostics.lst[[1]]
  plaus = check.plausibility(vd$accuracy, vd$sensitivity)
  #vd$fpr, vd$fnr, vd$ppv, vd$npv,
  species = sp.names[1]
  out.dat = c(model.type,model.name, species, plaus, vd$accuracy,  vd$sensitivity, vd$specificity, vd$true.positive, vd$false.positive, vd$true.negative, vd$false.negative)
  out.dat = matrix(out.dat, nrow = 1)
  colnames(out.dat) = out.hdr
  
  for (i in 2:length(sp.names)){
    vd = validation.diagnostics.lst[[i]]
    plaus = check.plausibility(vd$accuracy, vd$sensitivity)
    species = sp.names[i]
    add.dat = c(model.type, model.name, species, plaus, vd$accuracy,vd$fpr, vd$fnr, vd$ppv, vd$npv, vd$sensitivity, vd$specificity, vd$true.positive, vd$false.positive, vd$true.negative, vd$false.negative)
    out.dat = rbind(out.dat, add.dat)
  }
  
  col.names = TRUE
  # Don't re-write column names if appending to an existing file
  if (file.exists(results.table) & do.append == TRUE){
    col.names = FALSE
  }
  
  #Write the table
  write.table(out.dat, file = results.table, sep = ",", row.names = F, col.names = col.names, append = do.append)
  
  # Return results to workspace (not planning on using this, but it's an option)
  return(out.dat)
}

#' Main function for reading in x variables
#' 
#' @export setup.predictors
setup.predictors = function(sp.data, predictors.vec, predictors.paths){
  
  #**# Can move this into a sub-function to save space
  # Convert predictors.vec to numeric vector (should I just start with a numeric vector?)
  run.vec = c(0,0,0,0,0,0,0)
  path.vec = rep(NA, length(run.vec))
  
  for (i in 1:length(predictors.vec)){
    item = predictors.vec[i]
    item.path = predictors.paths[i]
    # For inputs from raw data
    if (item == "landuse"  ) { run.vec[1] = 1 ; path.vec[1] = item.path }
    if (item == "climate"  ) { run.vec[2] = 1 ; path.vec[2] = item.path }
    if (item == "biotic"   ) { run.vec[3] = 1 ; path.vec[3] = item.path }
    if (item == "frag"     ) { run.vec[4] = 1 ; path.vec[4] = item.path }
    if (item == "isolation") { run.vec[5] = 1 ; path.vec[5] = item.path }
    if (item == "survival" ) { run.vec[6] = 1 ; path.vec[6] = item.path }
    
    # For inputs taken from the metacommunity model files
    if (item == "landuse-sim") { run.vec[7] = 1 ; path.vec[7] = item.path }
  }
  
  # Create list starting with species data
  data.lst = list(sp.data)
  
  # Set up land use data
  if (run.vec[1] == 1){
    # Read in land use data (grazing, mowing, fertilziation, LUI)
    ludata = read.ludata(path.vec[1])  
    
    # Add dataset to list of datasets
    data.lst = append(data.lst, list(ludata))
  }
  
  # Set up climate data
  if (run.vec[2] == 1){
    # Read in climate data from Exploratories (will this actually be useful???)
    # Redownload the climate data using the new interface?
    # or just use the data that you have already?
    cddata = read.climate("climate/climate_data_year_qc3/plots_year_qc3.csv")
    
    # Add dataset to list of datasets
    data.lst = append(data.lst, cddata)
  }
  
  # Set up biotic data (note: this is going to be species data - so will probably be in the species file!)
  if (run.vec[3] == 1){
    stop("Need to calculate this based on lag terms in the data")
    
    data.lst = append(data.lst, biodata)
  }
  
  if (run.vec[4] == 1){
    stop("This data set has not yet been set up")
    # Read in fragmentation data #**# CORINE or something from Thomas Nauss (or both?)
    # Corine is available, but Nauss's data that I don't have is theoretically better.
    #frag = ETWAS
    
    data.lst = append(data.lst, frag)
  }
  
  if (run.vec[5] == 1){
    stop("Need to calculate this based on the biotic data")
    data.lst = append(data.lst, isolation)
  }
  
  if (run.vec[6] == 1){
    # Read in survival data
    survival.file = "C:/docs/beplants/datasets/EJB/Rechenmatrix_AHS.csv"
    sudata = get.survival(survival.file, ludata)
    
    data.lst = append(data.lst, sudata)
  }

  # Read in metacommunity landscape layer
  if (run.vec[7] == 1){
    raw.ludata = read.csv(path.vec[7])  
    
    if (nrow(raw.ludata) != 1){ stop("Only a single timestep is currently supported. Sorry for the inconvenience") }
    
    #Reformat to be useful
    plots = names(raw.ludata)[2:length(raw.ludata)]
    plotyears = sprintf("%s_%s",plots, Digit(1, 4)) #Only set up for one year of data!
    
    # Poor extraction method because more sophisticated approaches aren't working for me today
    values = c()
    for (val in 1:(ncol(raw.ludata) - 1)){
      this.value = raw.ludata[1,(val + 1)]
      values = c(values, this.value)
    }
    
    # This is a hacky approach that depends on there only being a single row in the input data.
    ludata = data.frame(plotyear = plotyears, landcover = values)
    
    # Add dataset to list of datasets
    data.lst = append(data.lst, list(ludata))
    
  }
  
  # Merge data sets into a single (giant) data object
  id.field = "plotyear"
  my.data = merge.data(data.lst, id.field)
  # (From a data management standpointis this a bad idea?)
  # If it works, it's fine, if you have memory problems, nothing you are doing at
  # this stage is especially multivariate, so you can read things in in chunks
  # and work with the data that way. Just need the plot-year ID for merging purposes.
  
  # Get index of which data columns correspond to which dataset for analysis purposes
  my.index = get.index(data.lst)
  
  # Create column index, this will indicate which data columns to include as predictor variables
  #**# Is this outdated - I should just be able to use my.index, right??
  col.index = rep(0, length(my.index)) # Initialize the variable with no variables as predictor variables
  col.index[ my.index > 1 ] = 1 # Set to use any columns that have a value > 1 in my.index May want to edit this if including other species as covariate.
  
  return(list(my.data, my.index, col.index))
}


#' Read in be_data file
#'
#' Reads the species data file
#' Reformat the file so that species are columns and plot-years are rows.
#' NOTE: species are coded by number #**# NEED TO LOOK UP FILE THAT HAS species ID's
#' NOTE: Not every species has an entry for every year!
#'
#' @param inpath the input path for the data
#' @param dropcol Whether or not to drop the added empty columns to make the landscape square for spatialdemography
#' @export read.sp.data
read.sp.data = function(inpath, sp.lookup.file, dropcol = TRUE, plotyear.rows = TRUE){
  bedata = read.csv(inpath)
  if (dropcol == TRUE){
    bedata = bedata[ ,1:153] #Drop extra blank columns for spatialdemography
  }
  
  sp.lookup = read.csv(sp.lookup.file)
  sp.lookup$Species = as.char(sp.lookup$Species)
  
  out.dat = bedata
  
  #Reformat to have plot years as rows, instead of columns
  if (plotyear.rows == TRUE){
    
    #**# DO I WANT PLOT AND YEAR COLUMNS??? Or just plotyear?
    
    #Get minimum and maximum year values
    min.year = min(bedata$TimeStep)
    max.year = max(bedata$TimeStep)
    
    #Set up plotyear.id
    #**# Think about making this its own function
    plotid.vec = names(bedata[4:153])
    years.vec = seq(min.year,max.year)
    plotid.vec2 = rep(plotid.vec, length(years.vec))
    years.vec2 = rep(years.vec, length(plotid.vec))
    years.vec2 = sort(years.vec2) # Put in order for joining purposes
    plotyear.id = sprintf("%s_%s", plotid.vec2, years.vec2)
    
    out.dat = data.frame(plotyear = plotyear.id)
    
    #Loop through species
    for (sp in 1:max(bedata$Species)){
      
      # Convert species from numbers to names
      sp.lbl = sp.lookup[sp, 1]
      
      #Create a label, since R doesn't like columns that start with numbers
      #sp.lbl = sprintf("sp_%s", sp)
      
      #Create a vector to hold output values
      val.vec = c()
      
      #Loop through timesteps
      for (ts in min.year:max.year){
        #Subset to data for this species and timestep AND drop Lifestage, Species and Timestep columns
        subdat = bedata[bedata$TimeStep == ts & bedata$Species == sp, 4:153]
        
        # Convert subdat to vector format
        #subdat = matrix(subdat, nrow = 1) #**# This does not seem to work - I'm still getting a list out
        #loop through stupid list that should be a vector and extract values
        sub.vec = c()
        for (val in subdat){
          #Check that the vector has valid values. If not, add 0's.
          if (length(val) == 0){
            sub.vec = c(sub.vec, 0)
          }else{
            sub.vec = c(sub.vec, val)  
          }
          
        }
        
        #Add to vector of values (assumes data are in the expected order (plots grouped by year, in order of ALB, HAI, SCH))
        val.vec = c(val.vec, sub.vec)
      }
      
      #Add val.vec to data frame
      out.dat[[sp.lbl]] = val.vec
    }
  }
  
  return(out.dat)
}

#' Read in land use data
#'
#' Reads in and reformats land use data to make for easy merging
#'
#' @param inpath input path for the land use file
#'
read.ludata = function(inpath){
  if (!require(myr)){
    stop("please intstall myr")
  }
  ludata = read.csv(inpath)
  
  # Add leading 0 to EPID
  plotpart = substr(ludata$EP,1,3)
  numpart = substr(ludata$EP,4,5)
  numpart = sapply(numpart, myr::Digit, 2)
  ludata$EP = sprintf("%s%s", plotpart, numpart)
  
  # Add plotyear field  
  ludata$plotyear = sprintf("%s_%s",ludata$EP, ludata$Year)
  
  #**# RECALCULATE LUI
  
  # Drop extra columns
  ludata$Year = NULL
  ludata$EP = NULL
  ludata$Exploratory = NULL
  ludata$Area = NULL
  ludata$N_org = NULL #Drop organic fertilizer
  ludata$N_min = NULL #Drop mineral fertilizer
  
  #Drop standardized values
  ludata$Fstd = NULL
  ludata$Mstd = NULL
  ludata$ Gstd = NULL
  
  return(ludata)  
}


#' Read in (yearly) climate data
#'
#' Read in climate data that has been quality controlled.
#' Note that 2008 data should probably be dropped:
#' it looks like some further quality control is needed.
read.climate = function(inpath){
  # Read in data
  cd = read.csv(inpath)
  
  # Create year field
  cd$year = substring(cd$datetime,1,4)
  
  # Drop 2008 data because I don't trust it
  cddata = cd[cd$year != 2008, ]
  
  # Create plotyear variable
  cddata$plotyear = sprintf("%s_%s", cddata$plotID, cddata$year)
  
  # Drop extra columns - for temperature, additional sensors do not appear to be contributing
  # novel information (~98% shared variation)
  cddata$datetime = NULL
  cddata$Ta_10 = NULL
  cddata$Ts_10 = NULL
  cddata$Ts_20 = NULL
  cddata$Ts_5 = NULL
  cddata$Ts_50 = NULL
  
  return(cddata)
}

#' Read in survival data
#' 
#' Also relevant as response traits to the metacommunity model
#' Data from Helge, Ute, & Eva
#' 
#' Hmm, I've been using plotyear as a joining method, but these data aren't associated with a particular year.
#' Need to add a second join to the main code to add the datasets that are not year specific
#' 
#' @param survival.data The survival data from Helge, Ute, & Eva
#' @param ludata Land use data from the Biodiversity Exploratories
get.survival = function(survival.file, ludata){
  # Read in survival data set
  
  sudata = read.csv(survival.data)
  head(sudata)
  
  # Note: not sure how many replicates per species there are - I think the sample design works for their
  # question, but may not be suitable for mine (at least not with a logistic exposure analysis)
  
  # Adjust file as needed (e.g., add plotyear, compute survivals, etc.)
  # Keep only relevant columns
  # For now, keeping it simple and only taking those columns I know how to use.
  # Note: there is a lot more information in this file!
  keeps = rep(0,ncol(sudata)) #Create a vector with the default to drop all columns
  keeps[3] = 1 # Keep species name for qc purposes
  keeps[4] = 1 # Keep Gsl_Nr for species joining purposes
  keeps[8] = 1 # Keep plot number
  keeps[ncol(sudata) - 11] = 1 # Keep survival at last time point
  su.sub = sudata[ ,keeps == 1]
  
  # Merge with landuse file (this requires remembering what year the study was carried out in - note some of the land use may have changed between years)
  
  # Planted 2012 & monitored in 2013 - one of the landuse files I have is missing these data (but the other has it? Just not LUI)
  
  # Estimate survival probabilities based on landuse type (logistic regression?)
  
  # output survival probabilities
  
  
}

#' Create fragmentation data set for use
#' 
#' Right now, using CORINE, because it is easy and available, and I'm running out of time
#' Alternatives will be to follow up with Thomas Nauss about options or to contact Catrin Westfall
#' Catrin has digitized about a 2km buffer around each plot and classified the land use, so would
#' have much more precise fragmentation statistics.
#' (if that is really the case, why am I re-creating the wheel?
#' Because rumor has it she's not particularly interested in sharing the data.
#' NOTE: anything you set up for CORINE will also work with her data set, so set it up for CORINE
#' And then follow up with her - the comparison could make an interesting note.
#' I wonder if there is already an R tool for this?
#' NOTE: Fragmentation and heterogeneity are two words for effectively the same thing.
#' Or you could just try FRAGSTATS - that might be faster and easier.
#' Quick look doesn't look like there is an easy way to do what I want - just do what you know (even if you haven't done it in R.)
#' 
#' @param points an input of points (even though the exporatories are plots, I'm using the center point)
#' @param dataset A data set to calculate fragmentation for
#' @param buffer.vec a vector of spatial scales at which to calculate fragmentation
#' 
compute.frag = function(points, dataset, buffer.vec){
  
  require(sp)
  require(raster)
  require(rgdal)
  
  # Points layer was converted from kml to shapefile in QGIS (something seems wrong about that - the kml format is so much better!)
  # But... there are only 122 EP - Grassland plots when I searched by EP - Grassland ! Uh oh. Maybe there is something more useful on BExIS? Forgot the VIP's.
  # Load points layer
  points = "GIS/BExIS_data/Grassalnd_EPs.shp" #Note, it's currently in WGS 1984, so we'll have to think about how to deal with that.
  in.points = rgdal::readShapePoints(points)
  
  #**# Left off here - time flies when you're doing things the stupid way (2 hours should have been more than enough for this!)
  
  # Load dataset layer
  if (dataset.type == "raster"){
    #in.data = scriptme!
  }
  
  # Add polygon function when you actually want to read in polygons.
  if (dataset.type == "polygon"){
    stop("Polygon support is not yet scripted")
  }
  
  # Create buffer polygons
  
  # Check if dataset is raster or polygon
  
  # loop through points 
  # & use buffer polygons to extract the information of interest
  
  # if dataset is raster, convert to polygon (only for a region surrounding the points - otherwise it will take forever)
  # or use crop raster tool
  
  # Convert from clipped raster to summary statistics of interest.
  
}


#' Merge data
#' 
#' Function to handle merging of data sets for species distribution modeling
#' 
#' @param data.lst a list of data files to be merged
#' @param id.field the field to use for merging the data
#' 
merge.data = function(data.lst, id.field){
  
  #Create a base data set
  my.data = data.lst[[1]]
  
  #Loop through and merge data sets
  if (length(data.lst) > 1){
    for (i in 2:length(data.lst)){
      set = data.lst[[i]]
      my.data = merge(my.data, set, by = id.field, all = FALSE) #all = FALSE drops incomplete rows! 
    }
  }

  return(my.data)
}

#' Get Index
#' 
#' Provides an index to which columns come from which data source for later data manipulation
#' Note that the first column in the combined data set needs to be plotyear, the ID field
#' Also note that the species data must be first in data.lst
#' 
#' @param A list of data sets requiring the index.
#' Each data set must have one, and only one identical ID field.
#' 
#' 
get.index = function(data.lst){

  # First column needs to be plotyear, the ID field
  my.index = c(0)
  
  # Species data needs to be the first data set in data.lst
  n.sp = length(data.lst[[1]]) - 1 # -1 is to drop plotyear field
  my.index = c(my.index, rep(1, n.sp)) # Add a 1 indicator for each species
  
  for (i in 2:length(data.lst)){
    this.set = data.lst[[i]]
    n.col = length(this.set) - 1 # Get number of columns in this data set (minus plotyear)
    this.val = i # Get a value to indicate this data set #**# I may want something more sophisticated here.
    my.index = c(my.index, rep(this.val, n.col))    
  }
  
  return(my.index)
}

#' Do descriptives
#'
#' Basic descriptive analysis on the assembled dataset. Highlight key insights
#' 
#' @param my.data A data set with data for descriptive analysis
#' @param my.index an index coding data for different types of descriptive analysis
#' @param desc.path a folder to put the descriptive analysis results
#' 
#' @export do.descriptives
do.descriptives = function(my.data, my.index, desc.path){
  message("This function needs testing")
  
  ## Add columns for plot, year, and exploratory
  plotyears = my.data$plotyear
  plots = sapply(plotyears, substr, 1,5)
  years = sapply(plotyears, substr, 7,10)
  expl = sapply(plotyears, substr, 1,3)
  
  my.data$plots = plots
  my.data$years = years
  my.data$expl = expl
  
  # update my.index (new columns are indicated with 0)
  my.index = c(my.index, rep(0, 3))
  
  my.subset = "all"
  # Look at data all together #**# Add subsets for by year, and by exploratory
  if (my.subset == "all") { in.data = my.data }
  
  # Do descriptive analysis for species data
  species.descriptive.analysis(in.data, my.index, desc.path)
  
  # Do descriptive analysis for predictor variables
  predictor.descriptive.analysis(in.data, my.index, desc.path)
      
}


#' descriptive analysis for species
#' 
#' 
species.descriptive.analysis = function(in.data, my.index, desc.path){

  # Drop non-species, non-identifier columns
  #sub.dat = in.data[  , my.index <= 1 ]
  # Drop non-species columns
  rownames(in.data) = in.data$plotyear
  sub.dat = in.data[  , my.index == 1 ]
  sub.dat = apply(sub.dat, c(1,2), as.num)
  
  
  # Species summary
  
  #**# Note: Rabinowitz used above/below median for rarity, at least according to Yu & Dobson 2000 who used the approach for mammals (http://www.jstor.org/stable/2655991?seq=1#page_scan_tab_contents).
  
  # Plot species occupancy distributions
  # Recode to presence/absence
  occ.dat = apply(sub.dat, c(1,2), pa.recode)
  # Get species richness per plot
  species.richness = get.sr(occ.dat)
  max.possible.richness = ncol(occ.dat)

  # Identify plots with a species richness of 0
  #**# Need to think about what to do with this - definitely should exclude from analysis (and need to figure out how to do that)
  plots.with.no.sp = species.richness[ species.richness == 0 ]
  
  
  # Get number of occupied plots per species
  plots.occupied = get.occ.plots(occ.dat)
  max.possible.plots = nrow(occ.dat)
  sp.not.present = plots.occupied[ plots.occupied == 0 ]
  length(sp.not.present) #39 species never present in any of the 3 time steps examined.
  
  # Set a threshold (arbitrary) for rarity by occupancy standards (10% = 15 plots for all exploratories, or 5 for single exploratories)
  occupancy.rarity.threshold = 0.1 # Use present on <= 10% of plots as "rare"
  occ.thresh.val = max.possible.plots * occupancy.rarity.threshold
  # Drop species present less than threshold
  message("Would you like to indicate your proposed threshold on the occupancy distributions above?")
  occ.plots.no.rare = plots.occupied[plots.occupied > occ.thresh.val]
  #**# Drops species number to 72 most common species.
  message("Is this really what you want to do? Also, what are you doing with this number after dropping the rare species? Outputting it somewhere, presumably")
  
  ## Plot species abundance distributions
  #**# How do you want to do this? Abundance per plot for each plot? Total abundance in landscape? All of the above?
  # to start - mirror the occupancy analysis
  plot.coverage = get.sr(sub.dat)
  message("This might not work. Revise function so that your axis legends at least make sense.")
  #**# Abundance is still right-skewed, but not nearly so much as species richness. -
  message("wait, why is there a cumulative score of at least 50 for occupied plots?? something's not right here")
  #**# they are coming in as factors - not as abundances! - need to fix this!
  
  species.abundance = get.occ.plots(sub.dat)
  
  # Set a threshold (arbitrary) for rarity by abundance standards (is at least 5% cover reasonable?)
  abundance.rarity.threshold = 65 #**# for cumulative abundance across all plots. Chosen because it looked like a natural break to me. Arguably 20 or 55 could also be chosen.
  
  # number rare species, number common species, number intermediate species
  
  # rare in terms of occupancy  
  # rare in terms of abundance
  # Number of species that occur: only ALB, only HAI, only SCH, ALB & HAI, ALB & SCH, HAI & SCH, ALB, HAI, SCH.
  
  # Check for irregularities:
  # does a widespread species disappear across a wide area just in one year then reappear the next year (saw this already, formally want to flag it)
  # Ask Helge about other potential irregularities?? He might criticize you for using this data, period.
  # - should look at alternative results if you just restrict to Helge's data, but I don't think they have enough timesteps to be useful.
  # - what other irregularities should there be? Is there a way to test for mis-identification of species? or other observer biases?  
  
  # Species plots (parasitize code from the visualize sdm for this) (or maybe just do this with the SDM??)
  #**# Should have kept those results in an easier to locate place! - may need to re-create them.
  # but this is the pre-subset data, so this might be the place to do it instead??
  
}

#' Get species richness
#' 
#' Helper function to streamline the above code
get.sr = function(occ.dat){
  species.richness = apply(occ.dat, c(1), sum)
  max.richness = max(species.richness)

  message("Figure out where to output the histogram")
  interval = 1
  hist(species.richness, breaks = seq(0,(max.richness + interval), interval))
  
  # Check data for plausibility
  if (min(species.richness) == 0){
    warning("Some plots have no species. This strikes me as unlikely. ")
    message("Figure out how to output a list of problem plots")
  }
  
  message("Think about how to deal with two upper-bound outliers: looks like they've got more than 10 more species than the next closest plots!")
  
  # Check shape of distribution
  message("Please visually inspect the histogram for occupancy distribution, as I do not know how to do this in R.")
  message("Google to see if there is a simple check for distribution shape!")
  #**# Overall dataset truncates unexpectedly on the left, and extends farther than expected on right.
  
  # Drop 0's and then log-transform, examine that distribution
  species.richness.drop.zeros = species.richness[species.richness > 0]
  ln.species.richness = log(species.richness)
  hist(ln.species.richness)
  #**# Note: more normally distributed, but still a short tail on the left with a broad shoulder on the right. Data are not log-normal either (but are more normal in the log-distribution)
  #**# So, the interquantile distance above the median is greater than the interquantile distance below the median.
  # I wonder if anyone uses the distance between the median and 75% quantile, or 95% quantile as a way of describing a distribution.
  message("Think about talking to Thomas Knieb about distributions")
  
  return(species.richness)
}

#' Get number of plots occupied by each species
get.occ.plots = function(occ.data){
  plots.occupied = apply(occ.dat, c(2), sum)
  max.plots = max(plots.occupied)
  
  interval = 1
  hist(plots.occupied, breaks = seq(0,(max.plots + interval), interval))
  #hist(plots.occupied, breaks = seq(0,(max.plots + interval), interval), xlim = c(0,100))
  
  message("What am I looking for here? What would be a problem?")
  #**# Distribution is zero-inflated - most species have no presences, some have one, and so forth, with some odd intermediate peaks (that may or may not be biological)
  message("How do I identify this distribution? Seems like this is something Preston or Hubbell worked with")
  
  message("Did you want to cut rare species? And what are you doing with the data set?")
  
  return(plots.occupied)
}

#' Predictor descriptive analysis
#' 
#' 
predictor.descriptive.analysis = function(in.data, my.index, desc.path){

  message("Needs testing" )
  pred.dat = in.data[, my.index > 1 ]
  
  a = summary(pred.dat)
  a
  # Max fertilization = 298, max mowing = 3, max grazing = 1401, max LUI =3.7166 (min = 0.3505)
  # medians: Fertilization = 0, mowing = 1, grzing = 59.41, LUI = 1.5224
  
  cor(pred.dat)
  require(psych)
  psych::pairs.panels(pred.dat, smooth = F, ellipses = F, scale = F)  
  
  ## For predictor data
  # -  Are all values plausible? #Maybe input a list of plausible ranges?
  # examine range & distribution of data & report fits
  # Also give plots of each dataset.
  # Examine multicollinearity
  
  ## Run princple component analysis, to get a simple overview of which variables are
  # strongly correlated with one another (does this really matter? Two variables can
  # be very strongly correlated, but still have very different mechanistic value!)
  #**# Maybe link to R script that already did this?
  #**# Or incorporate it as a helper script.
  # Not sure a PCA is really appropriate for the landuse data?? But I'm not really sure what will happen if you do one! Because the LUI is kind of a sort of composite variable as it is - so I'd say don't do a PCA for landuse.
  return(a)
  
}


#' Drop rare species
#'
#'
drop.rare.sp = function(my.data, my.index){
  #Default is keep all species
  # get number of species columns
  sp.cols = ncol(my.data) - length(id.keeps) - length(lu.keeps)
  sp.to.keep = rep(1, sp.cols)
  
  #Drop species with fewer than 10 occurrences overall
  if (drop.rare == 1){
    # Find out which species occur fewer than 10 times in the dataset
    message("Check to make sure no extremely abundant species are being excluded")
    
    ## Drop extremely rare species for which a species distribution model would be non-sensical (with fewer than 10 presences?)
    message("10 presences arbitrarily used as minimum threshold for species distribution modeling")
    message("Think about whether that's a good idea.")
    my.threshold = 10
    
    sp.recode = apply(my.data[2:(ncol(my.data) - 4)], c(1,2), myr::threshold.recode, 0.1)
    sp.count = apply(sp.recode, c(2), sum)
    plot(sp.count)
    sp.to.keep = sapply(sp.count, myr::threshold.recode, my.threshold)
    sum(sp.to.keep) #163 species with more than 10 occurrences. So, about 1/2.
    
    warning("How do you predict the behavior of really rare species? This is a key concept that is being left out here!")
    message("One thought would be to model rare species as a group and see if anything interesting shakes out.")
  }
  
  
  keeps = c(id.keeps, sp.to.keep, lu.keeps)
  sdm.data = my.data[ , keeps == 1] #Just keep columns where keeps == 1. Note: leaving out the == 1 part leads to very weird behavior!
  
}

#' Convert my.data to a format more easily used by the code
#' 
#' @export convert.sp.data
convert.sp.data = function(sp.data, year){
  
  # Subset data to chosen year (#**# I should have made this a function!)
  plotyears = sp.data$plotyear
  my.years = sapply(plotyears, substr, 7, 10)
  sp.data$year = my.years
  c.data = sp.data[sp.data$year == year, ]
  
  # Drop columns without species data
  c.data$plotyear = NULL
  c.data$year = NULL
  
  # Reformat data
  #LifeStage, TimeStep, Species, Cells1 - 150
  n.sp = ncol(c.data)
  
  actual.presences.lst = rep(list(NA), n.sp)
  # Loop through columns
  for (sp in 1:n.sp){
    plot.vals = c()
    # Loop through plots
    for (i in 1:nrow(c.data)){
      # Extract values
      this.val = c.data[i, sp]
      plot.vals = c(plot.vals, this.val) 
    }
    
    # Convert to presence absence
    plot.vals = sapply(plot.vals, pa.recode)
    actual.presences.lst[[sp]] = plot.vals
  }
  
  return(actual.presences.lst)   
}


#' Convert my.data to a format more easily used by the code
#' 
#' This version converted my.data to a SpeciesData.csv. Some of this code may be useful, but is not actually what I wanted in the end.
#' 
convert.sp.data.old = function(sp.data, model.file, year){
  
  # Subset data to chosen year (#**# I should have made this a function!)
  plotyears = sp.data$plotyear
  my.years = sapply(plotyears, substr, 7, 10)
  sp.data$year = my.years
  c.data = sp.data[sp.data$year == year, ]
  
  # Drop columns without species data
  c.data$plotyear = NULL
  c.data$year = NULL
  
  # Reformat data
  #LifeStage, TimeStep, Species, Cells1 - 150
  n.sp = ncol(c.data)
  lifestages = rep("Adults", n.sp)
  timesteps = rep(year, n.sp)
  species = seq(1,n.sp) # Create numeric sequential species ID's. Should think about adding look up tables for this! #**# this is something to watch out for in the future!
  loc.dat = data.frame(LifeStage = lifestages, TimeStep = timesteps, Species = species)
  
  # Add cells to dataframe
  for (i in 1:nrow(c.data)){
    this.cell = sprintf("Cell%s", i)
    #these.vals = c.data[i , ] # Take the ith row of species data #**# THis wasn't working
    #these.vals = matrix(these.vals, ncol = 1)
    
    these.vals = c()
    for (j in 1:ncol(c.data)){
      this.val = c.data[i, j]
      these.vals = c(these.vals, this.val)
    }
    
    loc.dat[[this.cell]] = these.vals
  }  
  
  # Save to file
  write.table(loc.dat, file = model.file, sep = ",", row.names = F, col.names = T, append = F)
}

#' Subset by traits
#' 
#' Plan is to give an option to subset the dataset based on preselected traits
#' 
#' @name subset.by.traits
subset.by.traits = function(sp.data, trait.data, link, trait.lst, trait.vals.lst){
  stop("Scripting in progress")
  
  # Drop extraneous columns from trait data
  
  # Join trait data to species data
  
  # Subset species data to just species with the required trait values
  
  # Return modified species data
  return(sp.data)
}

#' Read traits
#' 
#' Read in trait dataset from Helge. Note: in the future, you could make the trait source an input, and then use this for any chosen trait set.
#' 
#' @name read.traits
read.traits = function(){
  library(foreign)
  
  # Set hardwired paths
  bpath = "C:/docs/beplants/datasets/"
  hpath = sprintf("%sEJB/BeRich/TRAITS/", bpath)
  
  ## Read in raw data from EVA and INA
  #Eva #These trait data correspond to measurments taken on plants in the plots
  sp.tr.eva.file = sprintf("%sTRAITS_Eva_Plots_54_2011_eval2\\SPECIES_TRAITS.txt", hpath)
  sp.tr.eva.meta.file =  sprintf("%sTRAITS_Eva_Plots_54_2011_eval2\\SPECIES_TRAITS_DESCR.txt", hpath)        
  species.traits.eva<-read.table(sp.tr.eva.file, sep = ";", dec = ".", header=T, stringsAsFactors=F)
  species.traits.eva_<-read.table(sp.tr.eva.meta.file, sep = ";", dec = ".", header=T, stringsAsFactors=F)    
  
  #**# Okay, looks like the INA trait data are in the EVA columns. This makes things a bit confusing.
  #Ina #These trait data correspond to measurements taken on the seedlings to be planted(?)
  sp.tr.ina.file = sprintf("%sTRAITS_Ina_Species_3_2011_eval2\\SPECIES_TRAITS.txt", hpath)
  sp.tr.ina.meta.file = sprintf("%sTRAITS_Ina_Species_3_2011_eval2\\SPECIES_TRAITS_DESCR.txt", hpath)
  species.traits.ina<-read.table(sp.tr.ina.file, sep = ";", dec = ".", header=T, stringsAsFactors=F)
  species.traits.ina_<-read.table(sp.tr.ina.meta.file, sep = ";", dec = ".", header=T, stringsAsFactors=F)
  
  #Merge Ina and Eva 
  # Note: where these two datasets overlap, the values are identical, hence the merge
  spec.traits.all = merge(species.traits.eva, species.traits.ina, all = TRUE)
  
  return(spec.traits.all)
  
}

#**# NEEDS ADJuStMENT TO ADAPT TO CODE CHANGES
#-should be minor - just to the EVA/INA section?
#- do I really want to do this? Well, at least now I have an option
#' Split aggregated trait file into separate files
#' 
#' Trait data provided by Helge Bruelheide, and is subject to restrictions
hb.split = function(spec.traits.all){
  
  #Create a lookup variable for getting the species name that corresponds to the gsl_number
  gsl.lookup = spec.traits.all[ ,1:2]
  gsl.lookup = gsl.lookup[!duplicated(gsl.lookup), ] #Just get unique values
  names(gsl.lookup) = c("gsl_number","Species") #Rename ABBREVIAT to be species.
  
  #These allow adding the id and origin (eva vs. ina) for each dataset
  id = spec.traits.all[ ,1:2]
  origin = spec.traits.all[ ,1212]
  
  biolflor = spec.traits.all[ , 1:622]
  
  biopop = spec.traits.all[ , 623:1006] 
  biopop = cbind(id, biopop) #Add gsl & Abbrev 
  
  ellenberg = spec.traits.all[ ,1007:1089]
  ellenberg = cbind(id, ellenberg)
  
  evaina = spec.traits.all[ ,1090:1149]
  evaina = cbind(id, evaina, origin) #Add gsl, abbrev, and origin
  
  #Split the eva and ina datasets into separate datasets
  eva = evaina[evaina$origin == "EVA", ] #196 rows
  ina = evaina[evaina$origin == "INA", ] #150 rows
  #Drop the origin column now that all data have common origin
  eva$origin = NULL 
  ina$origin = NULL
  
  leda = spec.traits.all[ ,1150:1201]
  leda = cbind(id, leda)
  
  rothmaler = spec.traits.all[ ,1202:1211]
  rothmaler = cbind(id, rothmaler)
  
  #Create list of datasets
  data.sets = list(eva,ina,biolflor,biopop,leda,ellenberg,rothmaler)
  set.id = c("eva","ina","biolflor","biopop","leda","ellenberg","rothmaler")
  
  to.trans = c()
  #For each dataset, drop the dataset name from the column header
  for (i in 1:length(data.sets)){
    set = data.sets[[i]]
    set = set[!duplicated(set), ] #Drop duplicates
    set = col.strip(set) #Remove database name from trait title
    set$data.source = set.id[i]
    data.sets[[i]] = set
    
    #Get list of traits present
    to.trans = c(to.trans, names(set))
    
    # Write out each dataset
    out.set = sprintf("combined/%s_ebj.csv", set.id[i])
    write.table(set, file = out.set, sep = ",", row.names = F)
  }
  
  #Write list of names to translate
  #write.table(to.trans, file = "combined/totrans.csv", sep = ",", row.names = F)
  
  #Return data sets as a list
  return(list(data.sets,set.id,gsl.lookup))    
}

#' Year subset
#' 
#' Restrict to one year of the exploratories data, convert plots to row.names and drop plotyear column.
#' May want to back update other parts of the code that do essentially this!
#'
#' @param in.data Data set (e.g. sp.data) with a plotyear column containing plots and years and the rest as data columns
#' @param year the year of data to subset to. (what happens if you do it across years instead of across plots?)
#' @name year.subset
year.subset = function(in.data, year){
  # Subset data to chosen year (#**# I should have made this a function!)
  plotyears = in.data$plotyear
  my.years = sapply(plotyears, substr, 7, 10)
  sub.data = in.data[my.years == year, ]

  my.plots = sapply(sub.data$plotyear, substr, 1, 5)
  
  # Add plots as row.names
  rownames(sub.data) = my.plots
  
  # Drop columns without species data
  sub.data$plotyear = NULL

  return(sub.data)
}


#' Use metacom package to assess the biotic filter for Biodiversity Exploratories data
#' 
#' #**# Which data does it take as an input? Need to remember how I was doing this.
#' NOTE: requires an updated version of vegan to work.
#' @param sp.data Species data with a plotyear column containing plots and years
#' @param year the year of data to subset to. (what happens if you do it across years instead of across plots?)
#' 
#' @name metacom.analysis
#' @export metacom.analysis
metacom.analysis = function(sp.data, year){
  require(metacom)
  require(myr)

  # subset data to year & drop plotyear column
  my.data = year.subset(sp.data, year)
  
  #Convert data to presence/absence
  my.data = apply(my.data, c(1,2), myr::pa.recode)
  
  # Drop sites with no species
  # metacom does not allow empty sites. Interesting.Wait, what? how am I getting empty sites? I can see having species that are absent, but not sites?
  # Get list of rows with row sum of 0
  test = apply(my.data, c(1), sum) # if test == 0, the site has no species
  #Drop sites with no species
  plot.data = my.data[test != 0, ]
  
  # if sites were dropped, issue warning with site ID's
  if (nrow(plot.data != nrow(my.data))){
    dropped.sites = grep('\\b0\\b', test, value = TRUE)
    dropped.site.names = myr::listtotext(names(dropped.sites), separator = " ")
    drop.message = sprintf("The following sites had no species (!) and were dropped: %s; year == %s", dropped.site.names, year)
    message(drop.message)
  }
  
  # Drop missing species from plot
  col.test = apply(plot.data, c(2), sum)
  plot.data2 = plot.data[ , col.test != 0 ]
    
  # if species were dropped, issue warning with site ID's
  if (ncol(plot.data2) != ncol(plot.data)){
    dropped.sp = grep('\\b0\\b', col.test, value = TRUE)
    dropped.sp.names = myr::listtotext(names(dropped.sp), separator = ", ")
    drop.message = sprintf("The following species had no presences and were dropped: %s; year == %s", dropped.sp.names, year)
    message(drop.message)
  }
  
  # Plot with metacom  
  metacom::Imagine(plot.data2, fill = F)
  
  # SavePoint metacom
  
  #**# MOVE THIS STUFF
  # try principle component analysis on species data?
  a = princomp(plot.data2)
  #Error in princomp.default(plot.data2) : 
  #  'princomp' can only be used with more units than variables
  # Fail - but if I transpose it, have species be the rows and sites as variables, then what does it tell me? which species load onto which sites.
  # but I want to know which plots load on which species. Why do you need more units than variables for a principal component analysis? Isn't this just looking at collinearity??
  
  #library(psych) #**# Need this package to get this to work!
  #a = pairs.panel(plot.data2)
  a = cor(plot.data2)
  diag(a) = 0 # Diagonal == 1 by definition & is uninformative. Remove so I can use max.
  max(a) # still 1. Which species pairs???
  
  b = apply(a, c(1), min) #-0.61 is the strongest negative correlation
  d = apply(a, c(1), max) # 1 is the strongest positive correlation, but I don't really see it in the data as visualized with imagine.
  # Now I see it - a couple of sites have a couple of joint presences of just one species. Wait, why are there sites with just one species? Do we believe that?
  #**# This goes really slow in R studio & effectivly crashes it.
  #for (a.col in 1:ncol(a)){
  #  this.col = a[ ,a.col]
  #  val = max(this.col)
  #  if (val > 0.5){ message(sprintf("max r = %s, col = %s", val, a.col))}
  #  min.val = min(this.col)
  # if (min.val < -0.5){message(sprintf("min r = %s, col = %s", min.val, a.col))}
  #}
  
  a = pairs(plot.data2) # FAIL - figure margins too large (i.e., it wants to make a BIG graph)
  # Do metacom's metacommunity analysis
    # This takes a long time to run!
  outputs = metacom::Metacommunity(plot.data2)

  #**# Maybe try the metacom stuff with Helge's data?? Maybe it'll give better results? I'm a little concerned about data quailty, given some of the perfect correlations.
  
  # Looks like more likely to get species co-occuring than anti-occurring. But... biotic filter? Or facilitation? or just convergent habitat requirements?
  # and how do you separate those possibilities??
  
  # what happens if you exclude species that are only in one plot? Is that a legitimate exclusion?
  
  #**# MOVE THE ABOVE STUFF
  
  # Trying for a better display option, but this looks much worse!
  #Reimagine(plot.data2, fill = F)
  
}

#' Function to fix problems with Imagine function in metacom package
#' 
#' Copied the Imagine function from metacom, and then adapted it to my needs
#' First real use of the Open Source part of open source!
#' Well, I've got it plotting my labels, but they're so tiny, that they're illegible!
#' I think I need staggered labels = have some above and some below
#' 
Reimagine = function (comm, col = c(0, 1), ordinate = TRUE, scores = 1, fill = TRUE, 
                      xlab = "Species", ylab = "Site", yline = 2, xline = 2, sitenames = rownames(comm), 
                      speciesnames = colnames(comm), binary = TRUE, cex.x = 1.5, cex.y = 1.5,
                      staggered.x.labels = FALSE, staggered.y.labels = FALSE) 
{
  require(metacom)
  if (ordinate == TRUE) {
    comm = OrderMatrix(comm, binary = binary, scores = scores)
  }
  if (fill == TRUE) {
    for (i in 1:dim(comm)[2]) {
      temp = comm[, i]
      if (sum(temp) < 2) {
        comm[, i] = temp
      }
      else {
        first = min(which(temp > 0))
        last = max(which(temp > 0))
        comm[first:last, i] <- max(temp)
      }
    }
  }
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]
  par(mar = c(2, 6, 6, 1))
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm), col = col, 
        xlab = "", ylab = "", axes = FALSE)
  box()
  if (length(sitenames) > 1) {
    offset = rep(c("-","---------------"), length(sitenames))
    offset = offset[1:length(sitenames)]
    offset.sitenames = sprintf("%s%s",sitenames, offset)
    axis(2, at = 1:dim(comm)[1], labels = offset.sitenames, las = 1, 
         cex.axis = cex.y, lwd.ticks = 0)
  }
  if (length(speciesnames) > 1) {
    
    if (staggered.x.labels == FALSE){
      text(1:dim(comm)[2], par("usr")[4] + 1, labels = speciesnames, cex = cex.x, srt = 45, adj = 0, xpd = T)  
    }else{
      #Make this a function and apply it to the y labels too (although there srt probably won't be needed. And the side will be different)
      row1 = seq(1,dim(comm)[2],2)
      row2 = seq(2,dim(comm)[2],2)
      index = rep(c(1,0),dim(comm)[2]) # create a repeating vector of 0,1 longer than the rows (not just using dim(comm[2])/ 2 because it might be odd and be fractional)
      index = index[1:dim(comm)[2]] # Truncate to length of dim(comm[2])
      sp.names1 = speciesnames[index == 1]
      sp1.offset = rep(c("-","-------"),length(sp.names1))
      sp1.offset = sp1.offset[1:length(sp.names1)] #Truncate to length of sp.names1
      sp.names1 = sprintf("%s%s", sp1.offset, sp.names1)
      sp.names2 = speciesnames[index == 0]
      sp2.offset = rep(c("-","-------"), length(sp.names2))
      sp2.offset = sp2.offset[1:length(sp.names2)]
      sp.names2 = sprintf("%s%s", sp.names2, sp2.offset)
      #**# what if you added ---- before the label? Then they'd be staggered, but there would still be the line for continuity
      #**# looks like you need double-staggering of the labels - having them on top and bottom, and also staggered on top and bottom!
      # But it looks like I'm moving towards something usable? Maybe?
      # can set this based on the maximum label length - I think keeping it short with coding would help!
      # but you can also recode yourself, prior to putting the labels into the function. Then it should work.
      
      #**# Might consider a fourth layer of offsets? But that will make it hard to identify species on the bottom left.
      # I think it's always going to be a little bit painful!
      
      # Also, replace with number codes, and then return the look up for the number codes to the actual species.
      #**# Note: you should be getting actual species names, but you appear still not to be!
      text(row1, par("usr")[4] + 1, labels = sp.names1, cex = cex.x, srt = 45, adj = 0, xpd = T)  
      text(row2, par("usr")[3] - 1, labels = sp.names2, cex = cex.x, srt = 45, adj = 1, xpd = T)  
    }
    
    #axis(3, at = 1:dim(comm)[2], labels = speciesnames, cex.axis = cex.x, 
    #     lwd.ticks = 0)
  }
  mtext(xlab, 3, cex = 1.5, line = xline)
  mtext(ylab, 2, cex = 1.5, line = yline)
}


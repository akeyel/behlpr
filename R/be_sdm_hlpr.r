#### HELPER FUNCTIONS FOR BE_DATA PROCESSING ####

#' A main function to run the Species Distribution Model
#' 
#' 
#'  
#' @param my.data
#' @param my.index an index coding different column types.
#' type 0 corresponds to the ID field
#' type 1 corresponds to the species columns
#' type 2+ correspond to different suites of independent variables
#' @param col.index Contains indicators for whether or not to include the column as a predictor variable.
#' @param sdm.type 1 = basic sdm (nothing fancy), 2 = uses biomod
#' @param colonization.only An indicator for whether to restrict the data set to only
#' newly colonized plots #**# is this really a good idea? Seemed like one on the surface
#' but not it seems stupid - you lose stable populations and parts of the core distribution.
#' But you gain a weak attempt to look at dispersal limitation.
#' 
#' @export do.sdm
do.sdm = function(sdm.data, my.index, col.index, sdm.type, do.sdm.visualization, out.pdf, start.year, end.year, colonization.only = 0){
  
  index.table = table(my.index)
  n.sp = index.table[["1"]] # gives the count of columns that pertain to species
  
  #Change column to character to make later steps work
  sdm.data$plotyear = as.char(sdm.data$plotyear)

  # Create lists to contain outputs
  sp.thresholds.lst = validation.diagnostics.lst = rep(list(NA),n.sp)
  
  if (sdm.type == 2 | sdm.type == "biomod"){ env.lyrs = setup.env.lyrs(sdm.data) }

  # Set up plotting options
  if (out.pdf != "none" & do.sdm.visualization == 1){
    pdf(file = out.pdf) 
  }  
  
  ## Loop through species, and create species distribution model for each species
  for (i in 1:n.sp){  #note: species numbers are not sequential
    #i = 1
    # Get species column & name
    sp.col = i + 1 # + 1 adjusts for plotyear
    sp.name = names(sdm.data[sp.col]) 
    
    # Get data for focal species
    focal.species.all = sdm.data[[sp.col]] 
    
    #Recode to pres/abs for the focal species
    fs.presences.all = sapply(focal.species.all, threshold.recode, 0.1)
    
    # Basic sdm
    if (sdm.type == 1 | sdm.type == "basic"){
      basic.out = do.basic.sdm(sp.col, sdm.data, col.index, fs.presences.all, do.sdm.visualization, start.year, end.year)
      validation.diagnostics.lst[[i]] = basic.out[[1]]
      sp.thresholds.lst[[i]] = basic.out[[2]]
      other.info = sp.thresholds.lst
    }
    
    # Biomod2 sdm
    if (sdm.type == 2 | sdm.type == "biomod"){
      try.biomod(sp.name, fs.presences.all, env.lyrs)
    }
    
    if (sdm.type == 3 | sdm.type == "prefit"){
      do.sdm.prefit(STUFF)
    }
    
  }
  
  if (out.pdf != "none" & do.sdm.visualization == 1){
    dev.off()
  }
  
  return(list(validation.diagnostics.lst, other.info))
}

#' Basic SDM approach
#' 
#' Mainly for getting started - will want to do something more sophisticated for the final analysis
#' 
do.basic.sdm = function(sp.col, sdm.data, col.index, fs.presences.all, do.sdm.visualization, start.year, end.year, training.proportion = 2/3){
  
  # Make sure the species' own data won't be used as a predictor
  col.index[sp.col] = 0 # This does not affect the global col.index, just the local one in this function
  
  # Set number of training observations (I bet I can do the pre-fit just by setting training proportion to 0!)
  n.train = length(sdm.data) * training.proportion
  
  #Select subset of data for training and for validation
  training.sample = sample(sdm.data$plotyear, n.train)
  training.index = (sdm.data$plotyear %in% training.sample)
  #validation.index = !(sdm.data$plotyear %in% training.sample) #Can skip this and just use training.index != TRUE as in the command to create train.dat.
  
  train.dat = sdm.data[training.index == TRUE, ] #Training data set
  v.dat = sdm.data[training.index == FALSE, ]    #Validation data set
    
  focal.species.training = train.dat[[sp.col]] # For fitting the sdm parameters
  focal.species.validation = v.dat[[sp.col]] # For validating the sdm parameters
  
  fs.presences.training = sapply(focal.species.training, pa.recode)
  fs.presences.validation = sapply(focal.species.validation, pa.recode)  
  
  # Get thresholds for each species
  sp.thresholds = get.thresholds(train.dat, sp.col, col.index, sp.name, fs.presences.training)
  # var.names = sp.thresholds[1]
  # thresh.mins = sp.thresholds[2]
  # thresh.maxs = sp.thresholds[3]
  
  # Use thresholds to classify remaining plots as suitable or unsuitable
  suitability.vec = classify.plots(v.dat, sp.thresholds)  
  
  # Evaluate classification using validation data set
  validation.diagnostics = validate.sdm(fs.presences.validation, suitability.vec)  
  #sensitivity = validation.diagnostics[[1]]
  #specificity = validation.diagnostics[[2]]
  #accuracy = validation.diagnostics[[3]]
  #true.positive = validation.diagnostics[[4]]
  #false.positive = validation.diagnostics[[5]] 
  #true.negative = validation.diagnostics[[6]]
  #false.negative = validation.diagnostics[[7]]
  #positive.predictive.value = validation.diagnostics[[8]]
  #negative.predictive.value = validation.diagnostics[[9]]
  
  if (do.sdm.visualization == 1){
    # Create visualization of SDM for this species
    suitability.vec.all = classify.plots(sdm.data, sp.thresholds)
    
    # Calculate range of values used in training set, and mark values as extraploation or not
    extrapolation.index = extrapolation.test(sdm.data, train.dat)
    #extrapolation.index = rep("not calculated",length(sdm.data$plotyear)) #Temporary fix to make ploting work. This part of the code needs to be scripted.
    
    sp.name = names(sdm.data)[sp.col]
    
    # Create visual "map" of SDM for Exploratories
    visualize.sdm(sp.name, validation.diagnostics, sdm.data$plotyear, suitability.vec.all,
                  fs.presences.all, training.index, extrapolation.index, start.year, end.year)
    #**# the below does not appear to match the current function!
    #visualize.sdm(sdm.data$plotyear, suitability.vec.all, fs.presences.all, training.index, extrapolation.index)
    
  }
  
  return(list(validation.diagnostics, sp.thresholds)) # Combine sp.thresholds & visualization.info.lst, will be returned as an other.info object
  
}

#' Check SDM quality
#' 
#' 
check.sdm.quality.old = function(validation.diagnostics.lst, plot.title, out.pdf = "none"){
  
  #**# Add ability to export to pdf.
  
  # Extract accuracy, sensitivity, and specificity
  n.sp = length(validation.diagnostics.lst)
  acc.vec = sen.vec = spec.vec = rep(NA, n.sp)
  for (i in 1:n.sp){
    this.info = validation.diagnostics.lst[[i]]
    acc = this.info[[3]]
    sen = this.info[[1]]
    spec = this.info[[2]]
    
    acc.vec[i] = acc
    sen.vec[i] = sen
    spec.vec[i] = spec
    
  }
  
  # Plot accuracy, sensitivity, and specificity
  par(mfrow = c(1,2))  
  
  boxplot(acc.vec, ylim = c(0,1), ylab = "Accuracy") #**# check if this should be 1!
  plot(sen.vec, spec.vec, xlab = "Sensitivity", ylab = "Specificity", xlim = c(0,1), ylim = c(0,1))
  
}


#' Quick recode (move to myr!)
#' 
#' Recode a value if it matches a particular value
#' 
quick.recode = function(x, target, new.value){
  
  y = x
  if (x == target){
    y = new.value
  }
  
  return(y)
}

#' Check SDM quality
#' 
#' Plot results of the SDM, and color code SDM results that are considered "plausible"
#' 
check.sdm.quality.v2 = function(my.sdm.results, out.pdf = "none"){
  
  # Export to pdf
  if (out.pdf != "none"){
    pdf("sdm.quality.test.pdf") #**# This can be improved!
  }

  my.sdm.results = as.data.frame(my.sdm.results)
  
  # Plot sensitivity vs. specificity, with colors giving accuracy. Use symbol to denote plausibility
  sen = as.num(my.sdm.results$Sensitivity)
  spec = as.num(my.sdm.results$Specificity)
  acc = as.num(my.sdm.results$Accuracy)
  plaus = as.num(my.sdm.results$Plausible)
  
  # Create color ramp based on accuracy
  acc.col = sapply(acc, threshold.recode, 0.85) #**# FIX THIS LATER
  acc.col = sapply(acc.col,quick.recode, 1, "green")
  acc.col = sapply(acc.col, quick.recode, 0, 1)
  
  # Recode plaus to give the symbols I want
  plaus = sapply(plaus, quick.recode, 0, 5)
  plaus = sapply(plaus, quick.recode, 1, 16)
  
  #X axis is typically 1 - specificity, with y axis as sensitivity. Oops, I got taht switched.
  plot(spec, sen, xlab = "Specificity (True Negative Rate)", ylab = "Sensitivity (True Positive Rate)", col = acc.col, pch = plaus)
  
  
}


#' SDM Evaluate the fit of an externally defined relationship
#' 
#' Here, all data is used as validation, as the suitability relationship
#' was defined based on an external data set
#' 
#' 
do.sdm.prefit = function(STUFF){
  suitability.vec = my.data[years == "2009" ,sp.name] # Data from the previous year
  suitability.vec = sapply(suitability.vec, threshold.recode, 0.1)
  
  # Validate predictions based on previous year
  validation.diagnostics = validate.sdm(fs.presences.all, suitability.vec)  
  
  # Set up extrapolation.index
  extrapolation.index = rep(0,length(suitability.vec)) # For use when not extrapolating
  #extrapolation.index = rep("not calculated",length(sdm.data$plotyear)) # for use when extrapolation is uncertain.
  
  # Define all data to be validation data (none of it was used to parameterize the model)
  training.index = rep(0,length(suitability.vec))
  
  # Create visualization maps of the SDM data
  visualize.sdm(sp.name, sp.image, validation.diagnostics, sdm.data$plotyear, suitability.vec, fs.presences.all, training.index, extrapolation.index)
  
}
  

#' Helper function to set data up for biomod
#'
setup.env.lyrs = function(sdm.data){
  
  stop("This function needs to be recoded in the context of my.index and with critical thought")
  
  # subset & set up environmental layers for biomod (later in species loop)
  # Start with 4? Not sure how many I should do? Maybe just the non-biotic ones?
  env1 = sdm.data$LUI
  env2 = sdm.data$Fertilization_intensity
  env3 = sdm.data$Mowing_frequency
  env4 = sdm.data$Grazing_intensity
  env5 = my.data[years == "2009" ,"sp_1"] # Data from the previous year for species 1.
  env6 = sapply(env5, threshold.recode, 0.1)
  
  env.lyrs = list(env1,env2,env3,env4,env5,env6)
}

#' Get thresholds
#' 
#' Calculate a simple bioclimatic envelope based on extremes of known occurrence.\
#' Changed from previous version in that only extremes from presences are used as thresholds.
#' This is a dubious approach at best.
#' 
get.thresholds = function(train.dat, sp.col, col.index, sp.name, fs.presences){
  
  require(myr)
  target.value = 1
  train.dat = myr::matrix.subset(train.dat,col.index, target.value)
  
  var.names = colnames(train.dat)
  var.mins = rep(NA, length(var.names))
  var.maxs = rep(NA, length(var.names))
  
  # Calculate thresholds for each remaining variable in train.dat
  for (j in 1:ncol(train.dat)){
    
    # Calculate range of values of plots with and without species
    this.var = train.dat[ , j]
    var.name = var.names[j]
    with.sp = this.var[fs.presences == 1]
    
    # Only calculate mins & maxes if species is not entirely absent
    if (length(with.sp) != 0){
      w.min = min(with.sp)
      w.max = max(with.sp)    
      data.min = min(this.var)
      data.max = max(this.var)
      
      # Compare min/max values from with to data range
      
      # Create lower threshold (if applicable)
      if (w.min > data.min){  var.mins[j] = w.min  }
      
      # Create upper threshold (if applicable)
      if (w.max < data.max){  var.maxs[j] = w.max  }   
    }
  }
  
  # May want to reorganize into having all three elements together. But I think this is ok.
  return(list(var.names, var.mins, var.maxs))
}



#' Get thresholds
#' 
#' Calculate a simple bioclimatic envelope based on extremes of known occurrence
#' 
get.thresholds.old = function(train.dat, sp.col, col.index, sp.name, fs.presences){
  
  require(myr)
  target.value = 1
  train.dat = matrix.subset(train.dat,col.index, target.value)
  
  var.names = colnames(train.dat)
  var.mins = rep(NA, length(var.names))
  var.maxs = rep(NA, length(var.names))
  
  # Calculate thresholds for each remaining variable in train.dat
  for (j in 1:ncol(train.dat)){
    
    # Calculate range of values of plots with and without species
    this.var = train.dat[ , j]
    var.name = var.names[j]
    with.sp = this.var[fs.presences == 1]
    without.sp = this.var[fs.presences == 0]  
    
    # Only calculate mins & maxes if species is not present everywhere or nowhere (causes issues in the code)
    if (length(with.sp) != 0 & length(without.sp) != 0){
      w.min = min(with.sp)
      w.max = max(with.sp)    
      wo.min = min(without.sp)
      wo.max = max(without.sp)
      
      # Compare min/max values from with/without
      
      # Create lower threshold (if applicable)
      if (w.min < wo.min){  var.mins[j] = w.min  }
      
      # Create upper threshold (if applicable)
      if (wo.max > w.max){  var.maxs[j] = w.max  }   
    }
  }
    
  # May want to reorganize into having all three elements together. But I think this is ok.
  return(list(var.names, var.mins, var.maxs))
}

#' Classify Plots
#' 
#' Classify plots as either suitable or unsuitable based on the
#' thresholds previously identified
#' 
classify.plots = function(v.dat, sp.thresholds){
  
  # unpack lsts.out
  var.names = sp.thresholds[[1]]
  thresh.mins = sp.thresholds[[2]]
  thresh.maxs = sp.thresholds[[3]]
  
  suitability.vec = c()
  
  # For each plot in v.dat, classify it as suitable or not
  # based on the thresholds given in lsts.out
  for (j in 1:nrow(v.dat)){
    this.row = v.dat[j, ]
    
    # Make the default that the cell is suitable
    is.suitable = 1
    
    # Check upper thresholds
    is.suitable = check.thresholds(this.row, thresh.maxs, var.names, "upper", is.suitable)
    
    # Check lower thresholds (will bypass checks if input cell is already unsuitable)
    is.suitable = check.thresholds(this.row, thresh.mins, var.names, "lower", is.suitable)
    
    # Add this suitability record to the suitability.vec
    suitability.vec = c(suitability.vec, is.suitable)    
  }
  
  return(suitability.vec)
}

#' Check Thresholds
#' 
#' This function can be used to check EITHER an upper or a lower threshold.
#'  
check.thresholds = function(this.row, vals.lst, lbls.lst, threshold.type, is.suitable = 1){

  # If no thresholds to test, maintain input suitability.
  if (length(lbls.lst) != 0){
  
    # Loop through fields that have a threshold
    for (i in 1:length(lbls.lst)){
      
      # End for loop if cell is already unsuitable (is.suitable == 0) #**# Should I put a return within the function? I kind of like having the functions always exit at the end.
      if (is.suitable == 0){
        break
      }
      
      #Extract this record
      thresh.lbl = lbls.lst[i]
      thresh.val = vals.lst[i]
      
      #If there is no threshold, the threshold will be given as NA. No threshold will result in a classifciation of suitable.
      if (!is.na(thresh.val)){
        test.val = this.row[[thresh.lbl]]
        
        # Check that value is appropriate for the threshold
        
        #Check if value is greater than an upper threshold
        if (threshold.type == "upper"){
          if (test.val > thresh.val){
            is.suitable = 0
          }
        }
        
        #Check if value is lower than a lower threshold
        if (threshold.type == "lower"){
          if (test.val < thresh.val){
            is.suitable = 0
          }
        }  
      }
    }
  }
  return(is.suitable)
}

#' Validate SDM
#' 
#' Compare the predicted suitability from the classification with the observed
#' distribution in the validation data set
#' 
validate.sdm = function(fs.presences.validation, suitability.vec){
  
  if (length(fs.presences.validation) != length(suitability.vec)){
    stop("List of presences/absences must match list of predicted suitability!")
  }
  
  # Compare presence/absence list to predicted suitability list.
  true.positive = 0 # Predicted present and actually present
  false.positive = 0 # Predicted present but actually absent # I'm okay with a higher error here - may be dispersal limitation or competitive exclusion (i.e. habitat is suitable, but this species lost the lottery)
  true.negative = 0 # Predicted absent and actually absent
  false.negative = 0 # Predicted absent but actually present # Trying to minimize this
  
  for (i in 1:length(fs.presences.validation)){
    predicted = suitability.vec[i]
    observed = fs.presences.validation[i]
    
    # increment the appropriate category
    if (predicted == 1 & observed == 1){  true.positive = true.positive + 1    }
    if (predicted == 1 & observed == 0){  false.positive = false.positive + 1  }
    if (predicted == 0 & observed == 0){  true.negative = true.negative + 1    }
    if (predicted == 0 & observed == 1){  false.negative = false.negative + 1  }
  }
  
  # Calculate summary statistics
  N = length(fs.presences.validation)
  total.positive = sum(fs.presences.validation)
  total.negative = N - total.positive
  
  #positive.predictive.value = true.positive / (true.positive + false.positive)
  #negative.predictive.value = true.negative / (true.negative + false.negative)

  #false.positive.rate = false.positive / (false.positive + true.negative)
  #false.negative.rate = false.negative / (false.negative + true.positive)
  
  sensitivity = true.positive / total.positive
  specificity = true.negative / total.negative
  accuracy = (true.positive + true.negative) / N

  #, ppv = positive.predictive.value, npv = negative.predictive.value, fpr = false.positive.rate, fnr = false.negative.rate#**# Why do you care about the AUC score? Why not just care about the best choice of threshold's ability to classify?
  return(list(sensitivity = sensitivity, specificity = specificity, accuracy = accuracy, true.positive = true.positive, false.positive = false.positive, true.negative = true.negative, false.negative = false.negative))
}

#' Extrapolation test
#' 
#' Test whether examined values fall outside the range of the training data set
#' NOTE: AT PRESENT THIS DOES NOT TEST WHETHER *COMBINIATIONS* fall outside the range of the training data set
#' But that probably could be done with a convex hull or with the UTC approach
#' 
#' @param eval.dat the data to be evaluated for whether they fall outside the predicted range
#' @param train.dat the data set used in training
#' 
extrapolation.test = function(eval.dat, train.dat){
  
  #Check that number of columns are the same in the two data sets
  if (ncol(eval.dat) != ncol(train.dat)){
    stop("Evaluation data set must have the same columns as the training dataset!")
  }
  
  #start at 2 to avoid plotyear (#**# Should tie this to keep.ids so it is more general and robust!)
  col.start = 2
  
  min.vec = c(NA) #Start with NA to account for starting offset
  max.vec = c(NA)
  # Loop through variables in datasets
    
  for (a.col in col.start:ncol(train.dat)){
    this.col = train.dat[ ,a.col]
    train.min = min(this.col)
    train.max = max(this.col)
    
    min.vec = c(min.vec, train.min)
    max.vec = c(max.vec, train.max)
  }
  
  extrapolation.index = c()
  # Loop through records in evaluation data set & check whether the values are within the min/max
  for (a.row in 1:nrow(eval.dat)){
    
    is.extrapolation = 0
    #Loop through columns to see if they fit within the min/max
    for (a.col in col.start:ncol(eval.dat)){
      # If an extrapolation is found, break loop and mark the cell as an extrapolation
      #**# may want to change this to identify all the variables contributing to the extrapolation.
      if (is.extrapolation == 1){
        break
      }
      
      this.val = eval.dat[a.row, a.col] #Take the kth row, and mth column
      test.min = min.vec[a.col] 
      test.max = max.vec[a.col]
      
      if (this.val > test.max | this.val < test.min){
        is.extrapolation = 1
      }
      
      
    }
    extrapolation.index = c(extrapolation.index, is.extrapolation)
  }
  
  return(extrapolation.index)
}


#' Create a visual of the SDM predictions
#' 
#' For visually assessing the quality of the results. Plan is to plot suitability as the color of the cell,
#' plot a symbol to indicate presence or absence, color the cell border to indicate training or validation data
#' and finally to have a second cell border that indicates whether or not the cell contains extrapolated values beyond the validation set
#' The approach to changing the background color (later adapted to producing additional rectangles) was from Joran and stackoverflow.
#' Originally from an R help thread.
#'
#' @param plotlabels A list of the names of each plot, for use in identifying particular plots
#' @param suitability.vec.all
#' @param training.index Index to indicate which cells were in the training (T) and which were in the validation sets (F).
#' 
visualize.sdm = function(sp.name, validation.diagnostics, plotlabels, suitability.vec.all, fs.presences.all, training.index, extrapolation.index, start.year, end.year){
  
  #message("validation.diagnostics part has been temporarily turned off.")
  
  # Check that all input vectors are the correct length
  t1 = length(plotlabels)
  t2 = length(suitability.vec.all)
  t3 = length(fs.presences.all)
  t4 = length(training.index)
  t5 = length(extrapolation.index)
  
  if (t1 != t2 | t1 != t3 | t1 != t4 | t1 != t5){
    message(sprintf("Input vector lengths: %s, %s, %s, %s, %s", t1,t2,t3,t4,t5))
    stop("Input vectors are not all of the same length")
  }
  
  plots = substr(plotlabels,1,5)
  years = substr(plotlabels,7,10)
  
  #Loop through years
  for (year in start.year:end.year){

    #Create index for subsetting vectors
    year.index = grep(year, years)
    
    # Subset vectors to a single year (plots are for 1 year only, otherwise, I think it would be too much)
    pltlbls = plots[year.index]
    t1 = length(pltlbls) # Reassign t1
    suit.vec.all = suitability.vec.all[year.index]
    fs.pres.all = fs.presences.all[year.index]
    train.i = training.index[year.index]
    e.i = extrapolation.index[year.index]
    
    # Set up multipanel plot (NOT plotting this as a raster, more flexibility this way, I think)
    rows = 15
    cols = 10
    if (t1 != rows * cols){
      #**# temporarily disabling this
      #warning("rows * cols does not match data length. Incorrect plotting behavior may occur.")
    }
    
    # Plot one blank row to use for title purposes & add plant image
    # (I couldn't figure out how to add an image to the upper margin!)
    rows = rows + 1
    
    par(mfrow = c(rows,cols))
    par(oma = c(0,0,0,0))
    par(mar = c(0,0,0,0))
    
    # Plot blank row as header
    for (blank in 1:cols){
      plot(1,1, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", bty = "n", col = 0)
      
    }
    
    #my.image = matrix(seq(0.1,0.9,.1), nrow = 3)
    # Add image of plant in last cell
    my.image = sprintf("../datasets/plants/species_images/%s.jpg", sp.name)
    if (file.exists(my.image)){
      require(jpeg)
      my.image = readJPEG(my.image)  
      rasterImage(my.image, par("usr")[1], par("usr")[3],par("usr")[2],par("usr")[4])    
      #xleft, ybottom, xright, ytop
    }
    
    ##**# These are size sensitive - so a better solution or some adjustment may be necessary
    ## Add species name & year
    sp.name.year = sprintf("%s %s", sp.name, year)
    mtext(sp.name.year, side = 3, outer = TRUE, at = 0.5, line = -3)

    #, validation information, and species photo to plot
    #validation.info = sprintf("Accuracy = %s\nSensitivity = %s\nSpecificity = %s",
    #                          validation.diagnostics$accuracy,
    #                          validation.diagnostics$sensitivity,
    #                          validation.diagnostics$specificity)
    #
    #mtext(validation.info, side = 3, outer = TRUE, at = 0.1, line = -4.5)
    
    #Loop through cells in landscape
    for (i in 1:t1){
      
      #Unpack cell values
      suitability = suit.vec.all[i]
      presence = fs.pres.all[i]
      training = train.i[i]
      extrap = e.i[i]
      label = pltlbls[i]
      
      # Set suitability color
      if (suitability == 1){
        suit.col = "darkgreen"
      }else{
        suit.col = "darkorange3" #"firebrick3" #Not using red because of red/green colorblindness
      }
      
      # Set color for extrapolation border
      if (extrap == 1){
        extr.col = "gold"
      }else{
        extr.col = "gray"
      }
      
      # Set color for whether data is training or validation
      if (training == 1){
        train.col = "dodgerblue4"
      }else{
        train.col = "mediumpurple2"
      }
      
      # Plot this cell
      #1,1 makes a single point. When absent, this point is set to the background color.
      plot(1,1, xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
      
      #par(xpd = NA)
      
      #par("usr") is a vector (x1,x2,y1,y2)
      
      # Set up & plot training rectangle
      xleft = par("usr")[1]
      ybottom = par("usr")[3]
      xright = par("usr")[2]
      ytop = par("usr")[4]
      rect(xleft,ybottom,xright,ytop,col = train.col)
      
      #**# Watch for problems with plotting with the increments (not sure how that will change, depending on the size of the plot. but seems to work, so if it ain't broke, don't fix it!)
      
      # Set up & plot extrapolation rectangle
      xleft = par("usr")[1] + 0.05
      ybottom = par("usr")[3] + 0.05
      xright = par("usr")[2] - 0.05
      ytop = par("usr")[4] - 0.05
      rect(xleft,ybottom,xright,ytop,col = extr.col)
      
      # Set up & plot suitability rectangle
      xleft = par("usr")[1] + 0.1
      ybottom = par("usr")[3] + 0.1
      xright = par("usr")[2] - 0.1
      ytop = par("usr")[4] - 0.1
      rect(xleft,ybottom,xright,ytop,col = suit.col)
      
      # Replot the points over the rectangles if the species is present
      if (presence == 1){
        points(1,1, cex = 3, pch = 16, col = "deepskyblue")
      }
      
      # Add plot label
      x.center = (par("usr")[1] + par("usr")[2]) / 2
      y.height = par("usr")[3] + 0.15
      text(x.center, y.height, labels = sprintf("%s",label), col = "white")
      # have a second cell border indicate training or validation data
      # have cell border indicate whether values represent extrapolation
      # Color code cell as suitable or not
      # plot an symbol to indicate plant presences/absences
      
    }
  }
  
}



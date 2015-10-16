# Code for using biomod2 package with this data set

# I don't know that any of this is functional.

#' Test biomod2 with Biodiversity Exploratories data
#' 
#' NOTE: some things are being plotted at odd times and in
#' odd ways. This needs improvement
#' 
#' @param sp.name The name of the focal species
#' @param fs.presences.all List of presences and absences of the species
#' @param env.lyrs a list of environmental layers to be examined.
#' 
try.biomod = function(sp.name, fs.presences.all, env.lyrs){
  # load data
  ###################################################
  # load the library
  
  library(biomod2)
  
  # Name the focal species
  myRespName = sp.name
  
  # load our species data
  myResp = as.num(fs.presences.all)
  
  # Construct artificial X,Y coordinates
  # x is between 1 & 10
  # y is between 1 & 15 (displays the exploratories as blocks of 5 rows each
  # These have to be offset by 1/2 to get the points to be IN the cells!
  x.base = seq(0.5,9.5)
  y.base = seq(0.5,14.5)
  x = rep(x.base, length(y.base))
  y = rep(y.base, length(x.base))
  y = sort(y)
  myRespXY = data.frame(X = x, Y = y)
  
  #convert environmental layers to raster & add to raster stack
  first.lyr = env.lyrs[[1]]
  first.lyr = matrix(first.lyr, ncol = length(x.base), byrow = T)
  first.raster = raster(first.lyr, xmn = 0, xmx = 10, ymn = 0, ymx = 15)
  myExpl = stack(first.raster)
  
  # Loop through other layers & add them to the stack
  if (length(env.lyrs) > 1){
    for (lyr in 2:length(env.lyrs)){
      this.lyr = env.lyrs[[lyr]]
      this.lyr = matrix(this.lyr, ncol = length(x.base), byrow = T)
      this.lyr =  raster(this.lyr, xmn = 0, xmx = 10, ymn = 0, ymx = 15)
      myExpl = stack(myExpl, this.lyr)
    }
  }
  
  message("Figure out how to split the data set (I think this means chaning presence/absence to NA, and then having those records in the evaluation data set")
  
  # Format Data
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName)
  # Check data                         
  myBiomodData
  
  # Plot data
  plot(myBiomodData)
  
  # Set model options
  message("Blindly using default biomod2 options. (BAD IDEA)")
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  # Do modeling
  message("Blindly using biomod2 model setup. (OK IDEA?)")
  myBiomodModelOut <- BIOMOD_Modeling( 
    myBiomodData, 
    models = c('SRE','CTA','RF','MARS','FDA'), 
    models.options = myBiomodOption, 
    NbRunEval=3, 
    DataSplit=80, 
    Prevalence=0.5, 
    VarImport=3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
  
  # Examine model summary
  myBiomodModelOut 
  
  # get all models evaluation                                     
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # print the dimnames of this object
  dimnames(myBiomodModelEval)
  
  # let's print the TSS scores of Random Forest
  myBiomodModelEval["TSS","Testing.data","RF",,]
  
  # let's print the ROC scores of all selected models
  myBiomodModelEval["ROC","Testing.data",,,]
  
  # Examine variable importance
  get_variables_importance(myBiomodModelOut)
  
  # Construct an ensemble model
  message("Again, blindly using biomod2 defaults for the ensemble model")
  myBiomodEM <- BIOMOD_EnsembleModeling( 
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = c(0.7),
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  #**# And this fails because none of the models meet the ensemble threshold!
  
  # print summary                     
  myBiomodEM
  
  # get evaluation scores
  get_evaluations(myBiomodEM)
  
  # Project current predictions 
  #**# do I want to try to do this with my version?
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  
  # summary of created oject
  myBiomodProj
  
  # Look at files on harddrive (I think this needs updating, but I don't know where it will be!
  list.files("sp.1/proj_current/")
  
  # make some plots sub-selected by str.grep argument
  plot(myBiomodProj, str.grep = 'MARS')
  
  # if you want to make custom plots, you can also get the projected map
  myCurrentProj <- get_predictions(myBiomodProj)
  myCurrentProj
  
  # Plot ensemble forecast
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
    EM.output = myBiomodEM,
    projection.output = myBiomodProj)
  
  myBiomodEF
  plot(myBiomodEF)
  
  #**# NOTE I HAVE NOT DONE FUTURE PROJECTIONS YET 
  # BUT THAT IS SOMETHING TO CONSIDER ADDING!
  
}


#' Test biomod
#' 
#' Runs the Gulo gulo vignette from biomod2
#'
test.biomod2 = function(){
  
  # load data
  ###################################################
  # load the library
  library(biomod2)
  
  # load our species data
  DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                      package="biomod2"))
  
  head(DataSpecies)
  
  # the name of studied species
  myRespName <- 'GuloGulo'
  
  # the presence/absences data for our species 
  myResp <- as.numeric(DataSpecies[,myRespName])
  
  # the XY coordinates of species data
  myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
  
  
  # load the environmental raster layers (could be .img, ArcGIS 
  # rasters or any supported format by the raster package)
  
  # Environmental variables extracted from Worldclim (bio_3, bio_4, 
  # bio_7, bio_11 & bio_12)
  myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
                               package="biomod2"),
                  system.file( "external/bioclim/current/bio4.grd", 
                               package="biomod2"), 
                  system.file( "external/bioclim/current/bio7.grd", 
                               package="biomod2"),  
                  system.file( "external/bioclim/current/bio11.grd", 
                               package="biomod2"), 
                  system.file( "external/bioclim/current/bio12.grd", 
                               package="biomod2"))
  
  
  ###################################################
  ### code chunk number 3: formating_data
  ###################################################
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName)
  
  
  ###################################################
  ### code chunk number 4: print_formating_data
  ###################################################
  myBiomodData
  
  
  ###################################################
  ### code chunk number 5: plot_formating_data
  ###################################################
  plot(myBiomodData)
  
  
  ###################################################
  ### code chunk number 6: modeling_options
  ###################################################
  # 2. Defining Models Options using default options.
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  
  ###################################################
  ### code chunk number 7: modeling
  ###################################################
  # 3. Computing the models 
  
  myBiomodModelOut <- BIOMOD_Modeling( 
    myBiomodData, 
    models = c('SRE','CTA','RF','MARS','FDA'), 
    models.options = myBiomodOption, 
    NbRunEval=3, 
    DataSplit=80, 
    Prevalence=0.5, 
    VarImport=3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
  
  
  
  ###################################################
  ### code chunk number 8: modeling_summary
  ###################################################
  myBiomodModelOut 
  
  
  ###################################################
  ### code chunk number 9: modeling_model_evaluation
  ###################################################
  # get all models evaluation                                     
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # print the dimnames of this object
  dimnames(myBiomodModelEval)
  
  # let's print the TSS scores of Random Forest
  myBiomodModelEval["TSS","Testing.data","RF",,]
  
  # let's print the ROC scores of all selected models
  myBiomodModelEval["ROC","Testing.data",,,]
  
  
  
  ###################################################
  ### code chunk number 10: modeling_variable_importance
  ###################################################
  # print variable importances                                    
  get_variables_importance(myBiomodModelOut)
  
  
  ###################################################
  ### code chunk number 11: ensemble_modeling
  ###################################################
  myBiomodEM <- BIOMOD_EnsembleModeling( 
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = c(0.7),
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  
  ###################################################
  ### code chunk number 12: ensemble_modeling_outputs
  ###################################################
  # print summary                     
  myBiomodEM
  
  # get evaluation scores
  get_evaluations(myBiomodEM)
  
  
  ###################################################
  ### code chunk number 13: projection_curent
  ###################################################
  # projection over the globe under current conditions
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = 'current',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = F,
    output.format = '.grd')
  
  # summary of crated oject
  myBiomodProj
  
  # files created on hard drive
  list.files("GuloGulo/proj_current/")
  
  
  
  ###################################################
  ### code chunk number 14: projection_curent_plot
  ###################################################
  # make some plots sub-selected by str.grep argument
  plot(myBiomodProj, str.grep = 'MARS')
  
  
  ###################################################
  ### code chunk number 15: projection_curent_getProj
  ###################################################
  # if you want to make custom plots, you can also get the projected map
  myCurrentProj <- get_predictions(myBiomodProj)
  myCurrentProj
  
  
  ###################################################
  ### code chunk number 16: projection_future
  ###################################################
  # load environmental variables for the future. 
  myExplFuture = stack( system.file( "external/bioclim/future/bio3.grd",
                                     package="biomod2"),
                        system.file( "external/bioclim/future/bio4.grd",
                                     package="biomod2"),
                        system.file( "external/bioclim/future/bio7.grd",
                                     package="biomod2"),
                        system.file( "external/bioclim/future/bio11.grd",
                                     package="biomod2"),
                        system.file( "external/bioclim/future/bio12.grd",
                                     package="biomod2"))
  
  myBiomodProjFuture <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExplFuture,
    proj.name = 'future',
    selected.models = 'all',
    binary.meth = 'TSS',
    compress = 'xz',
    clamping.mask = T,
    output.format = '.grd')
  
  
  
  
  ###################################################
  ### code chunk number 17: projection_current_plot
  ###################################################
  # make some plots, sub-selected by str.grep argument
  plot(myBiomodProjFuture, str.grep = 'MARS')
  
  
  ###################################################
  ### code chunk number 18: EnsembleForecasting_current
  ###################################################
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
    EM.output = myBiomodEM,
    projection.output = myBiomodProj)
  
  
  ###################################################
  ### code chunk number 19: EnsembleForecasting_loading_res
  ###################################################
  myBiomodEF
  
  
  ###################################################
  ### code chunk number 20: EnsembleForecasting_plotting_res
  ###################################################
  # reduce layer names for plotting convegences
  plot(myBiomodEF)
  
  
}

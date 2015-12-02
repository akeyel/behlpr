# Goal is to run metacom function on the cluster, as it seems to take forever.
# Give it a week to finish, if no luck, re-evaluate!

cluster.stuff = function(){ 
  # Set up cluster stuff. Change number of cores if running in parallel.
  library(parallel)
  cores        <- 1#6                # The number of cores I would like to use in paralell
  workers      <- makeCluster(getOption("cl.cores", cores))     # Start 16 worker process
  #clusterEvalQ(workers, library(spatstat))   # load package spatstat on each worker(FD in your case)
  
  
  #Main function to run on the cluster
  #Note that outfile is NOT defined here - it is defined in the specific workspaces above.
  do.metacom = function(x){
    # Load workspace
    library(metacom) # Note: this will give warnings, because for some reason this package depends on devtools!
    load("metacom/SavePoint_preMetacommunity_v2.RData")
    outputs = metacom::Metacommunity(plot.data2)
    save.image("metacom/Output_2009.RData")
  }
  
  #Save the function so you can load it into each of the workers
  save.image("metacom/do_metacom_function.RData")
  x = c(1)
  
  clusterEvalQ(workers, load("metacom/do_metacom_function.RData"))          # Adds my Analysis function to the workspaces of the 16 workers
  clusterApplyLB(cl=workers, x, do.metacom)
  # This is the main call for the analysis
  stopCluster(cl=workers)    # Once all analyses are done the 16 workers are closed
}




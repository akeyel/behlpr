# First task is to 

library(parallel)
cores        <- 1#6                # The number of cores I would like to use in paralell
workers      <- makeCluster(getOption("cl.cores", cores))     # Start 16 worker process
#clusterEvalQ(workers, library(spatstat))   # load package spatstat on each worker(FD in your case)

#As I understand Daniel's approach to parallelization, you
  #1. create workspaces with unique settings
  #2. call a function to execute in each workspace (oh, name is the name of the workspace - it's specifically loaded in the function!

#Create unique settings for each workspace
for (i in 1:16){
    outfile = sprintf("outputs/test%s.csv",i)
    save.image(sprintf("testfiles/input%s.RData",i))
    }

#Get a listing of the unique workspace settings
my.files = dir("testfiles")

#Main function to run on the cluster
  #Note that outfile is NOT defined here - it is defined in the specific workspaces above.
ParTest = function(nam){
    load(sprintf("testfiles/%s",nam))
    a = runif(100000000,0,1)
    cat(sum(a),file = outfile)
    }
    
#Save the function so you can load it into each of the workers
save.image("partest.rData")


clusterEvalQ(workers, load("partest.rData"))          # Adds my Analysis function to the workspaces of the 16 workers
clusterApplyLB(cl=workers, x=my.files, fun=function(x) ParTest(nam=x))
                                                              # This is the main call for the analysis
stopCluster(cl=workers)    # Once all analyses are done the 16 workers are closed





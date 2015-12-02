#Convert CORINE data set to buffers that can be used in analyses
#Author: A.C. Keyel
#Last Updated: 2015 Nov 11
#Requirements: ArcGIS 10.1

#Load required modules
from numpy import *
import arcpy

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

#Turn Overwrite on - keep from getting 1,000,000 files!
arcpy.env.overwriteOutput = True

arcpy.env.pyramid = "NONE"

#' Download CORINE 2006 (COR06) Data
COR06 = "C:/docs/beplants/datasets/GIS/CORINE_2006/g100_06.tif"

#' Converted CORINE 2006 to ESRI GRID format (was having alignment issues with the tif)
# Set snapping to the created COR06 GRID layer.

exploratories = ["ALB","HAI","SCH"]
ending = "_CORQ"
#' Cut COR06 data to match each exploratory
# *_CORQ
COR06_cut = [expl + ending for expl in exploratories] # Umm, that was way harder than what I do in R!

#' Reproject COR06 to match exploratories (UTM 32 for HAI & ALB, 33 for SCH). Chose cubic resampling - oh, I bet I know where I got the all forest cells from last time... /tab Better to have everything in same coord. system
#**# For now, skip this and see whether I get hung up.

#' Start with plot shapefiles (UTM 32 for HAI & ALB, 33 for SCH)  /tab I made these with ArcGIS based on the .kml files /cr
basepath = "C:/docs/beplants/datasets/GIS/"
plotpath = basepath + "Plots/"
bufferpath = basepath + "buffers/"
corinepath = basepath + "CORINE_2006/"
plotfiles = ["ALB_UTM32.shp","HAI_UTM32.shp","SCH_UTM33.shp"]
buffer_sizes = [500, 1000, 2500, 5000]

make_buffers = 0
subset_plots = 0
reproject_raster = 0
reclass_cor = 0
extract_cor = 1

# Create buffers
for plot in plotfiles:
    for my_buffer in buffer_sizes:
        #plot = plotfiles[0]
        #my_buffer = buffer_sizes[0]

        # Set up UTM zone & projection
        utmzone = "32"
        if plot[0:3] == "SCH":
            utmzone = "33"

        prj_file = basepath + "prj/" + "WGS 1984 UTM Zone %sN.prj" % utmzone
            

        plot_file = plotpath + plot
        buffer_file = bufferpath + plot[0:4] + "_%s" % my_buffer + ".shp" #**# change to plot[0:3] to drop extra whitespace
        buffer_subset = bufferpath + plot[0:3] + "_%s" % my_buffer + "EP" + ".shp"
        expl_raster = corinepath + "%s_CORQ.tif" % plot[0:3]
        expl_raster_proj = corinepath + "%s_utm%s" % (plot[0:3].lower(), utmzone)
        rc_raster_file = expl_raster_proj + "rc"
        
        if make_buffers == 1:
            # Create buffer file
            print "Warning: This code has been modified since it was actually run to produce data. Check for bugs"
            #NOTE: arcpy.Buffer_arc DOES NOT WORK!!!
            arcpy.Buffer_analysis(plot_file, buffer_file, my_buffer)
        
        if subset_plots == 1:        
            # Create a subset of plots to only include the EP's (why I didn't do this before making the buffers is beyond me)
            print "Warning: This code has been modified since it was actually run to produce data"
            whereClause = "\"PopupInfo\" = 'EP - Grassland' OR \"PopupInfo\" = 'VIP - Grassland'"
            my_lyr = "bufferlyr"
            arcpy.MakeFeatureLayer_management(buffer_file, my_lyr)
            arcpy.SelectLayerByAttribute_management(my_lyr, "NEW_SELECTION", whereClause)
            arcpy.CopyFeatures_management(my_lyr, buffer_subset) #If the input is a layer which has a selection, only the selected features will be copied.

        if reproject_raster == 1:
            ## slow step!
            raise ValueError("This wasn't working right and needs troubleshooting")
            #**# this step hung up for about 3 hours and I aborted it. It took under 3 seconds in ArcMap, so I just did the reprojections manually in ArcMap.
            # Reproject CORINE data to be in the appropriate coordinate system for calculations
            arcpy.ProjectRaster_management(expl_raster, expl_raster_proj, prj_file, "NEAREST")
            
        if reclass_cor == 1:
            # Reclassify Corine Subset to be Grassland/Arable, Forest, Urban and Other
            reclass_field = "Value" #**# Check that this is true!
            remap = arcpy.sa.RemapRange([[1,11,3],[12,22,1],[23,34,2],[35,255,4]])
            rc_raster = arcpy.sa.Reclassify(expl_raster_proj, reclass_field, remap, "NODATA")
            rc_raster.save(rc_raster_file)

        # Extract CORINE data by buffer
        if extract_cor == 1:
            
            # Loop through each polygon and extract the raster uniquely
            # probably a better approach, but this should work

            # Loop through plots (50 - known number)
            for i in range(0,1):#50:
                # Add ID to buffer subset file
                # Just using FID - check that this works.

                # Select each plot in the buffer & export to own shapefile (will be manually deleted at some point)
                this_lyr = "single_plot"
                this_plot = basepath + "temp/" + "%s_%s.shp" % (plot[0:3], i)
                whereClause = "\"FID\" = %s" % i
                arcpy.MakeFeatureLayer_management(buffer_subset, this_lyr)
                arcpy.SelectLayerByAttribute_management(this_lyr, "NEW_SELECTION", whereClause)
                arcpy.CopyFeatures_management(this_lyr, this_plot) #If the input is a layer which has a selection, only the selected features will be copied.

                # Extract from raster based on new shapefile
                plot_raster = basepath + "temp/" + "%s_%sr" % (plot[0:3], i)
                new_raster = arcpy.sa.ExtractByMask(rc_raster_file, this_plot)
                new_raster.save(plot_raster)

                # Compute summary info based on raster
                # umm, how do I do this? Should just be value, count - should be in the raster itself.
                #**# search cursor?? #**# Look up previous code for this after dinner. Need to get this running before bed.

                # Append summary info to main summary document


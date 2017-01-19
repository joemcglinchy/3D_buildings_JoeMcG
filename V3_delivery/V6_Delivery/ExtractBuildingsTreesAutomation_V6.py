# ---------------------------------------------------------------------------
# Name: ExtractBuildingsTreesAutomation.py
# Purpose:
# Usage:
# Description:
# Author:       Yiqun Xie, Joseph McGlinchy
# Organization: Esri Inc.
#
# Created:     06/25/2015 Yiqun Xie, Joseph McGlinchy
# Modified:    12/22/2015 Roslyn Dunn Version 1
#              Expand the algorithm to run over large areas in production
# Modified:    01/05/2016 Roslyn Dunn Version 2
#              Remove describe of file gdb, since gdb's are corrupted during
#                file deletion.
#              Output .tif file from Zonal Stats to avoid having SRS modified
#                on output.
#              Require output folder to exist (changed to input parameter instead
#                of output.
#              Fix issues with deletion of scratch workspaces (bug in core)
# Modified:    02/12/2016 Roslyn Dunn Version 3
#              Added functionality to construct surfaces directly from LAS
# Modified:    04/05/2016 Joe McGlinchy Version 4.0
#              Ensured correct execution from Desktop 10.4 and Pro 1.2 environments by including
#              BuildRasterAttributeTable around line 340
# Modified:    04/27/2016 Joe McGlinchy Version 5.0
#              Modified LAS Dataset creation to add LAS files with full paths and not use arcpy.env.workspace variable
#              Added try-except around arcpy.AddGeometryAttributes to add fields as GEODESIC. This may result in
#              other errors if neither condition is appropriate.
#              removed print() statement in Except block and changed to arcpy.AddError
# Modified:    08/25/2016 Roslyn Dunn Version 6.0
#              Checked for availability of Spatial Analyst and 3D Analyst licenses
#              Check these licenses in if returning mid-stream (for example, if no buildings found in current area)
#              Check to see if DSM and DTM have any data. Return if either surface is all NoData.
#              Modify boundary simplification to NOT check for errors, since this sometimes hangs. Changed Handle
#                topological errors from RESOLVE_ERRORS to NO_CHECK.
#              Take SR units into account when setting cell size in LAS_Raster_type.art.xml
#              Send most raster output to .tif files instead of Esri GRID to avoid a problem with spontaneous
#                 reprojection cause because the Esri GRID files don't support certain projections. The only
#                 exception is the CleanTreeR<n> raster, since the cleanfast routine gets different results otherwise.
#              Work around a bug in Create LAS Dataset Layer (in Pro) which requires that we only specify existing
#                 classification codes.
# ---------------------------------------------------------------------------

import sys
import os
import inspect
import traceback
import math
import multiprocessing
import time
import shutil
# import resourceLogger
import arcpy
import csv
from arcpy.sa import *


def extract_buildings_trees(pp_params):
    # In order to invoke this method as a separate process using python multiprocessing,
    #   import arcpy and other site packages
    import arcpy
    import sys
    import os
    import traceback
    import time
    try:
        start = time.time()
        # check out Spatial Analyst license
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            arcpy.AddMessage("\n*** No Spatial Analyst License Available."
                             "Exiting extract_buildings_trees for OID: {0} ***".format(pp_params[1]))
            sys.exit()
        arcpy.env.overwriteOutput = True
        # Unpack pp_params
        out_folder_path = pp_params[0]
        # arcpy.AddMessage("\n out_folder_path: {0}".format(out_folder_path))
        oid = pp_params[1]
        arcpy.AddMessage("\n** Processing oid: {0} **\n".format(oid))
        # set a flag to indicate if input data type is LAS or DSM/DTM rasters
        las_input = False
        if pp_params[3] == "ALL RETURNS" or pp_params[3] == "LAST RETURNS":
            las_input = True
            lasd = pp_params[2]
            # arcpy.AddMessage("\n lasd: {0}".format(lasd))
            las_return_type = pp_params[3]
            # arcpy.AddMessage("\n las_return_type: {0}".format(las_return_type))
            if arcpy.CheckExtension("3D") == "Available":
                arcpy.CheckOutExtension("3D")
            else:
                arcpy.AddMessage("\n*** No 3D Analyst License Available."
                                 "Exiting extract_buildings_trees for OID: {0} ***".format(pp_params[1]))
                sys.exit()
        else:
            dsm_md = pp_params[2]
            # arcpy.AddMessage("\n dsm_md: {0}".format(dsm_md))
            dtm_md = pp_params[3]
            # arcpy.AddMessage("\n dtm_md: {0}".format(dtm_md))
        xmin_str = pp_params[4]
        # arcpy.AddMessage("\n xmin_str: {0}".format(xmin_str))
        xmax_str = pp_params[5]
        # arcpy.AddMessage("\n xmax_str: {0}".format(xmax_str))
        ymin_str = pp_params[6]
        # arcpy.AddMessage("\n ymin_str: {0}".format(ymin_str))
        ymax_str = pp_params[7]
        # arcpy.AddMessage("\n ymax_str: {0}".format(ymax_str))
        featureextraction = pp_params[8]
        # arcpy.AddMessage("\n featureextraction: {0}".format(featureextraction))
        elevation_meter_scalefactor_str = pp_params[9]
        # arcpy.AddMessage("\n elevation_meter_scalefactor_str: {0}".format(elevation_meter_scalefactor_str))
        use_pos_terrain_method = pp_params[10]
        # arcpy.AddMessage("\n use_pos_terrain_method: {0}".format(use_pos_terrain_method))
        delete_intermediate_files = pp_params[11]
        # arcpy.AddMessage("\n delete_intermediate_files: {0}".format(delete_intermediate_files))
        height_path = pp_params[12]
        # arcpy.AddMessage("\n height_path: {0}".format(height_path))
        las_point_spacing = pp_params[13]   # only utilized when las_input = True
        # arcpy.AddMessage("\n las_point_spacing: {0}".format(las_point_spacing))
        dsm_class_codes = pp_params[14]     # only utilized when las_input = True
        dtm_class_codes = pp_params[15]     # only utilized when las_input = True

        # Convert some passed parameters to float (all parameters were passed as strings)
        elevation_meter_scalefactor = float(elevation_meter_scalefactor_str)
        xmin = float(xmin_str)
        xmax = float(xmax_str)
        ymin = float(ymin_str)
        ymax = float(ymax_str)

        # Each fishnet section gets a separate directory (named after the OID) to contain results
        out_subfolder_path = os.path.join(out_folder_path, oid)

        # The results file gdb will store trees and buildings
        results_gdb_name = r"Results.gdb"
        # Entire path of the file gdb
        results_file_gdb_path = os.path.join(out_subfolder_path, results_gdb_name)

        trees_output = ""
        if "TREES" in featureextraction.upper():
            trees_output = os.path.join(results_file_gdb_path, r"trees_to_merge" + oid)
        buildings_output = os.path.join(results_file_gdb_path, r"buildings_to_merge" + oid)

        # arcpy.env.workspace = out_subfolder_path
        # arcpy.env.workspace = "in_memory"
        if not os.path.exists(out_subfolder_path):
            arcpy.AddMessage("Creating Results sub-Folder:    " + out_subfolder_path)
            os.makedirs(out_subfolder_path)
        elif arcpy.Exists(buildings_output):
            # if desired outputs already exist (buildings_to_merge* and possibly trees_to_merge*),
            #  then don't continue, since this run has be re-started
            if "TREES" in featureextraction.upper():
                # If trees_to_merge* feature class exists then this oid has already been processed
                if arcpy.Exists(trees_output):
                    arcpy.AddMessage("OID {0} already processed...skipping".format(oid))
                    arcpy.CheckInExtension('3D')
                    arcpy.CheckInExtension('Spatial')
                    return
            else:
                arcpy.AddMessage("OID {0} already processed...skipping".format(oid))
                arcpy.CheckInExtension('3D')
                arcpy.CheckInExtension('Spatial')
                return

        # If the results file gdb doesn't exist, then create it
        if not os.path.exists(results_file_gdb_path):
            # arcpy.AddMessage("Creating Results File GDB:  {0}".format(results_file_gdb_path))
            arcpy.CreateFileGDB_management(out_subfolder_path, results_gdb_name, out_version="CURRENT")

        # Create a scratch file gdb for intermediate file output
        scratch_gdb_name = r"TempWorkArea.gdb"
        # Entire path of the file gdb
        scratch_file_gdb_path = os.path.join(out_subfolder_path, scratch_gdb_name)
        # If the file gdb doesn't exist, then create it
        if not os.path.exists(scratch_file_gdb_path):
            # arcpy.AddMessage("Creating Scratch File GDB:  {0}\n".format(scratch_file_gdb_path))
            arcpy.CreateFileGDB_management(out_subfolder_path, scratch_gdb_name, out_version="CURRENT")

        # send intermediate files to the scratch file gdb
        arcpy.env.workspace = scratch_file_gdb_path

        # clip dsm raster from either LAS Dataset or DSM Mosaic Dataset
        clip_dsm_raster = os.path.join(out_subfolder_path, "dsm{}.tif".format(oid))
        arcpy.AddMessage("\nDSM clip_raster: {0}".format(clip_dsm_raster))
        # Extend the clip footprint by 1 unit in each direction so resulting features overlap
        fishnet_rectangle = "{} {} {} {}".format(xmin - 1, ymin - 1, xmax + 1, ymax + 1)
        if not arcpy.Exists(clip_dsm_raster):
            if las_input:
                if las_return_type == "ALL RETURNS":
                    return_val = ""
                else:
                    # Las Return is better for defining buildings
                    return_val = "'Last Return'"
                # use LAS Class codes 0-6 to create the DSM
                # 0 - Never Classified, 1 - Unassigned, 2 - Ground, 3 - Low Vegetation,
                # 4 - Medium Vegetation, 5 - High Vegetation, 6 - Building
                arcpy.MakeLasDatasetLayer_management(lasd, out_layer="DSM_LASD_Layer",
                                                     class_code=dsm_class_codes,
                                                     return_values=return_val, no_flag="INCLUDE_UNFLAGGED",
                                                     synthetic="INCLUDE_SYNTHETIC", keypoint="INCLUDE_KEYPOINT",
                                                     withheld="EXCLUDE_WITHHELD", surface_constraints="")

                # messages = arcpy.GetMessages()
                # arcpy.AddMessage("\nResults output from MakeLasDatasetLayer of DSM is: \n{0}".format(messages))

                arcpy.env.extent = fishnet_rectangle
                arcpy.LasDatasetToRaster_conversion("DSM_LASD_Layer", clip_dsm_raster, value_field="ELEVATION",
                                                    interpolation_type="BINNING MAXIMUM NATURAL_NEIGHBOR",
                                                    data_type="FLOAT", sampling_type="CELLSIZE",
                                                    sampling_value=las_point_spacing, z_factor="1")
                # messages = arcpy.GetMessages()
                # arcpy.AddMessage("\nResults output from LasDatasetToRaster of DSM is: \n{0}".format(messages))
            else:
                # Input is DSM and DTM Mosaic datasets
                try:
                    arcpy.Clip_management(dsm_md, fishnet_rectangle, clip_dsm_raster, "#", "-3.40282346639e+038")
                    # messages = arcpy.GetMessages()
                    # arcpy.AddMessage("\nResults output from Clip of DSM is: \n{0}".format(messages))
                except:
                    arcpy.AddMessage("No data in DSM clip of area with oid: {0}".format(oid))
                    arcpy.CheckInExtension('3D')
                    arcpy.CheckInExtension('Spatial')
                    return

        # clip dtm raster from either LAS Dataset or DTM Mosaic Dataset
        clip_dtm_raster = os.path.join(out_subfolder_path, "dtm{}.tif".format(oid))
        arcpy.AddMessage("\nDTM clip_raster: {0}".format(clip_dtm_raster))
        if not arcpy.Exists(clip_dtm_raster):
            if las_input:
                # use LAS Class codes 2 and 8 to create the DTM
                # 2 - Ground, 8 - Model Key
                arcpy.MakeLasDatasetLayer_management(lasd, out_layer="DTM_LASD_Layer", class_code=dtm_class_codes,
                                                     return_values="", no_flag="INCLUDE_UNFLAGGED",
                                                     synthetic="INCLUDE_SYNTHETIC", keypoint="INCLUDE_KEYPOINT",
                                                     withheld="EXCLUDE_WITHHELD", surface_constraints="")
                # messages = arcpy.GetMessages()
                # arcpy.AddMessage("\nResults output from MakeLasDatasetLayer of DTM is: \n{0}".format(messages))
                arcpy.env.extent = fishnet_rectangle
                arcpy.LasDatasetToRaster_conversion("DTM_LASD_Layer", clip_dtm_raster, value_field="ELEVATION",
                                                    interpolation_type="BINNING AVERAGE NATURAL_NEIGHBOR",
                                                    data_type="FLOAT", sampling_type="CELLSIZE",
                                                    sampling_value=las_point_spacing, z_factor="1")
                # messages = arcpy.GetMessages()
                # arcpy.AddMessage("\nResults output from LasDatasetToRaster of DTM is: \n{0}".format(messages))
            else:
                # Input is DSM and DTM Mosaic datasets
                try:
                    arcpy.Clip_management(dtm_md, fishnet_rectangle, clip_dtm_raster, "#", "-3.40282346639e+038")
                    # messages = arcpy.GetMessages()
                    # arcpy.AddMessage("\nResults output from Clip of DTM is: \n{0}".format(messages))
                except:
                    arcpy.AddMessage("No data in DTM clip of area with oid: {0}".format(oid))
                    if delete_intermediate_files == "true":
                        remove_rasters(out_subfolder_path)
                        remove_shapefiles(out_subfolder_path)
                        arcpy.env.workspace = ""   # so we can successfully delete scratch_file_gdb_path
                        remove_filegdb(scratch_file_gdb_path)
                    arcpy.CheckInExtension('3D')
                    arcpy.CheckInExtension('Spatial')
                    return

        # Check the DTM and DSM for data. If either surface is entirely NoData then exit
        # allnodata_DTM_result = arcpy.GetRasterProperties_management(clip_dtm_raster, property_type="ALLNODATA",
        #                                                         band_index="Band_1")
        # messages = arcpy.GetMessages()
        # arcpy.AddMessage("\nMessages ALLNODATA output from GetRasterProperties are: \n{0}\n".format(messages))
        # all_nodata_DTM = int(allnodata_DTM_result.getOutput(0))
        # arcpy.AddMessage("\nResults ALLNODATA output from GetRasterProperties are: \n{0}\n".format(all_nodatavalue))
        # allnodata_DSM_result = arcpy.GetRasterProperties_management(clip_dsm_raster, property_type="ALLNODATA",
        #                                                             band_index="Band_1")
        # all_nodata_DSM = int(allnodata_DSM_result.getOutput(0))

        try:
            clip_dsm_ras = arcpy.Raster(clip_dsm_raster)
            dist_dsm = arcpy.RasterToNumPyArray(clip_dsm_ras, "", "", "", nodata_to_value=0)
            dsm_locs = dist_dsm.nonzero()
            arcpy.AddMessage("len(dsm_locs): {0}".format(len(dsm_locs[0])))

            clip_dtm_ras = arcpy.Raster(clip_dtm_raster)
            dist_dtm = arcpy.RasterToNumPyArray(clip_dtm_ras, "", "", "", nodata_to_value=0)
            dtm_locs = dist_dtm.nonzero()
            arcpy.AddMessage("len(dtm_locs): {0}".format(len(dtm_locs[0])))
        except Exception as e:
            arcpy.AddMessage(e.message)
            arcpy.CheckInExtension('3D')
            arcpy.CheckInExtension('Spatial')
            # Return any Python specific errors and any error returned by the geoprocessor
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            arcpy.AddError("PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
            messages = arcpy.GetMessages()
            arcpy.AddError("GP Messages:\n {0}".format(messages))
            raise

        # if all_nodata_DTM == 1 or all_nodata_DSM == 1:
        if len(dsm_locs[0]) == 0 or len(dtm_locs[0]) == 0:
            arcpy.AddMessage("Either the DTM or DSM is entirely NoData in region {0}".format(oid))
            if delete_intermediate_files == "true":
                remove_rasters(out_subfolder_path)
                remove_shapefiles(out_subfolder_path)
                arcpy.env.workspace = ""   # so we can successfully delete scratch_file_gdb_path
                remove_filegdb(scratch_file_gdb_path)
            arcpy.CheckInExtension('3D')
            arcpy.CheckInExtension('Spatial')
            return

        # Set output coordinate system (need this to keep it from changing)
        desc = arcpy.Describe(clip_dtm_raster)
        sr = desc.spatialReference
        arcpy.env.outputCoordinateSystem = sr
        linear_unit = sr.linearUnitName.upper()
        # dtm_SR = sr.exportToString()
        # arcpy.AddMessage("Spatial ref of all outputs should be equivalent to the DTM SR: \n{0}\n".format(dtm_SR))

        trees_output = ""
        if "TREES" in featureextraction.upper():
            trees_output = os.path.join(results_file_gdb_path, r"trees_to_merge" + oid)
        buildings_output = os.path.join(results_file_gdb_path, r"buildings_to_merge" + oid)

        # Minus - compute difference to get heights relative to the ground)
        # Note the Height raster is put into a separate folder (height_path)
        # (Used to create a Mosaic Dataset later in zonal stats step)
        diff_ori = os.path.join(height_path, "Minus_img" + oid + r".tif")
        arcpy.gp.Minus_sa(clip_dsm_raster, clip_dtm_raster, diff_ori)

        # Divide by a scale factor to convert all heights to meters
        # divide_minus1 = os.path.join(out_subfolder_path, "Divide_Minus1" + oid)
        flt_meter_elev_unit = float(elevation_meter_scalefactor)

        divide_minus1 = Raster(diff_ori) / flt_meter_elev_unit
        arcpy.Delete_management(flt_meter_elev_unit)

        # arcpy.AddMessage("divide_minus1.maximum: {0}".format(divide_minus1.maximum))
        # If there are no elevation values > 2 meters then return
        if divide_minus1.maximum <= 2.0:
            arcpy.AddMessage("No Buildings or Trees were identified in region {0}".format(oid))
            # print "No Buildings or Trees were identified in region: " + oid
            if delete_intermediate_files == "true":
                remove_rasters(out_subfolder_path)
                remove_shapefiles(out_subfolder_path)
                arcpy.env.workspace = ""  # so we can successfully delete scratch_file_gdb_path
                remove_filegdb(scratch_file_gdb_path)
            del desc, sr
            arcpy.CheckInExtension('3D')
            arcpy.CheckInExtension('Spatial')
            return

        # create a raster with all floating point values of 1.0
        # Note: Done because SetNull & Con won't take a constant 1.0 value for the false raster parameter
        flt_raster_one = divide_minus1 / divide_minus1

        setnull_divi1 = SetNull(divide_minus1, flt_raster_one, "VALUE<2")
        arcpy.Delete_management(flt_raster_one)

        arcpy.env.outputCoordinateSystem = sr
        # arcpy.AddMessage("Current SR (before IsNull): {}".format(sr.exportToString()))
        # Is Null - create a mask to indicate the areas where heights are < 2 (1: heights < 2, 0: heights >= 2)
        isnull_newma1 = os.path.join(out_subfolder_path, "IsNullnewm" + oid + r".tif")
        arcpy.gp.IsNull_sa(setnull_divi1, isnull_newma1)

        # Make Raster Layer from the previously created mask
        mask_null = "mask_null" + oid
        arcpy.MakeRasterLayer_management(isnull_newma1, mask_null, "", "", "")
        # arcpy.AddMessage("mask_null Layer SR {}".format(arcpy.Describe(mask_null).spatialReference.exportToString()))
        # Select Layer By Attribute - Select those pixels in the layer where heights are < 2
        arcpy.SelectLayerByAttribute_management(mask_null, "NEW_SELECTION", "VALUE=1")

        # Euclidean Distance - disBD represents the distance of each 'tall' pixel to the closest 'short' pixel
        #  If a tree exists, the center (trunk) of the tree will have a larger value (it's toward the center)
        #  If a building, the center of the building will have a larger value
        disbd = os.path.join(out_subfolder_path, "disBD" + oid + r".tif")
        arcpy.gp.EucDistance_sa(mask_null, disbd, "", mask_null, "")
        arcpy.Delete_management(mask_null)

        # Set Null - set all ZERO (and very small distance) values to NoData
        disbdnull = os.path.join(out_subfolder_path, "disBDnull" + oid + r".tif")
        arcpy.gp.SetNull_sa(disbd, disbd, disbdnull, "VALUE<0.0001")

        # Negate the Distances to create wells (as opposed to peaks)
        # Now the peaks we have in the middle of the objects become basins and the boundary pixels become ridges
        rdsm = os.path.join(out_subfolder_path, "Negatediff" + oid + r".tif")
        arcpy.gp.Negate_sa(disbdnull, rdsm)

        # Flow Direction
        flowdir_filt1 = os.path.join(out_subfolder_path, "FlowDirFlt" + oid + r".tif")
        arcpy.gp.FlowDirection_sa(rdsm, flowdir_filt1, "NORMAL", "")

        # Basin
        basin_flowdi5 = os.path.join(out_subfolder_path, "BasinFlowD" + oid + r".tif")
        arcpy.gp.Basin_sa(flowdir_filt1, basin_flowdi5)

        # Times
        diff = os.path.join(out_subfolder_path, "diff" + oid + r".tif")
        arcpy.gp.Times_sa(divide_minus1, setnull_divi1, diff)
        arcpy.Delete_management(divide_minus1)

        # Focal Statistics
        focalst_diff1 = os.path.join(out_subfolder_path, "FocalStdif" + oid + r".tif")
        arcpy.gp.FocalStatistics_sa(diff, focalst_diff1, "Rectangle 3 3 CELL", "MEAN", "DATA")

        # Times
        mean_diff = os.path.join(out_subfolder_path, "mean_diff" + oid + r".tif")
        arcpy.gp.Times_sa(focalst_diff1, setnull_divi1, mean_diff)
        arcpy.Delete_management(setnull_divi1)

        # Minus
        diff_minus_avg = os.path.join(out_subfolder_path, "diffMinusA" + oid + r".tif")
        arcpy.gp.Minus_sa(diff, mean_diff, diff_minus_avg)

        # Send the output to a different file name so results can be compared for the 2 methods
        if not use_pos_terrain_method == "true":
            # use the 'slope' method by default
            positive = os.path.join(out_subfolder_path, "Slope" + oid + r".tif")
            arcpy.gp.Slope_sa(diff_minus_avg, positive, "DEGREE", "1")
        else:
            # the 'Positive Terrains' method
            input_raster_or_constant_value_2 = "0.3"
            positive = os.path.join(out_subfolder_path, "GreaterMin" + oid + r".tif")
            arcpy.gp.GreaterThan_sa(diff_minus_avg, input_raster_or_constant_value_2, positive)

        # Set output coordinate system (need this to keep it from changing in Zonal Stats step)
        # desc = arcpy.Describe(clip_dtm_raster)
        # sr = desc.spatialReference
        # arcpy.env.outputCoordinateSystem = sr
        # linear_unit = sr.linearUnitName.upper()
        # arcpy.AddMessage("Spatial ref of DSM clipped raster is: \n{0}\n".format(sr.exportToString()))

        # Zonal Statistics
        # first run build raster attribute table on basin_flowdi5
        arcpy.BuildRasterAttributeTable_management(basin_flowdi5, "Overwrite")

        # Output a .tif file to avoid modification of SRS (bug which applies only to GRID)
        zonalst_basi4 = os.path.join(out_subfolder_path, "ZonalStBas" + oid + r".tif")
        arcpy.gp.ZonalStatistics_sa(basin_flowdi5, "VALUE", positive, zonalst_basi4, "MEAN", "DATA")

        # Iso Cluster Unsupervised Classification
        # arcpy.AddMessage("Classifying...")
        isocluster2 = os.path.join(out_subfolder_path, "isocluster" + oid + r".tif")
        # Write out a permanent signature file to avoid conflicts with other concurrent processes
        iso_sig_file = os.path.join(out_subfolder_path, r"iso_sig" + oid + r".gsg")
        # wrap this in a try block in case zero (0) classes are found (this happens in desert areas)
        try:
            arcpy.gp.IsoClusterUnsupervisedClassification_sa(zonalst_basi4, "2", isocluster2, "20", "10", iso_sig_file)
            # messages = arcpy.GetMessages()
            # arcpy.AddMessage("\nResults from IsoClusterUnsupervisedClassification_sa of oid {0} are:"
            #                  " \n{1}\n".format(oid, messages))

        except:
            arcpy.AddMessage("No Buildings or Trees were identified in region {0}".format(oid))
            if delete_intermediate_files == "true":
                remove_rasters(out_subfolder_path)
                remove_shapefiles(out_subfolder_path)
                arcpy.env.workspace = ""   # so we can successfully delete scratch_file_gdb_path
                remove_filegdb(scratch_file_gdb_path)
            del desc, sr
            arcpy.CheckInExtension('3D')
            arcpy.CheckInExtension('Spatial')
            return
        # check actual number of classes obtained
        raster_iso = arcpy.Raster(isocluster2)
        buildingsfound = True
        if int(raster_iso.maximum) == int(raster_iso.minimum):
            # Only one class found = assumed to be trees for now, but this might change
            arcpy.AddMessage("No buildings were identified in region {0}".format(oid))
            # arcpy.AddMessage("No records will be output to: {0}".format(buildings_output))
            buildingsfound = False

        # Always look for buildings regardless if asked for, since we look for trees where no buildings exist
        if buildingsfound:
            # Set Null
            setnullbd = os.path.join(out_subfolder_path, "SetNullbd" + oid + r".tif")
            arcpy.gp.SetNull_sa(isocluster2, isocluster2, setnullbd, "VALUE>1")

            # Raster to Polygon
            bdiso = os.path.join(scratch_file_gdb_path, "bdiso" + oid)
            arcpy.RasterToPolygon_conversion(setnullbd, bdiso, "NO_SIMPLIFY", "VALUE")

            # Add Geometry Attributes
            try:
                arcpy.AddGeometryAttributes_management(bdiso, "AREA;PERIMETER_LENGTH", "METERS", "SQUARE_METERS", "")
            
            except Exception as e:
                arcpy.AddMessage(e.message)
                arcpy.AddWarning("*** Dataset is not planar, accounting for geodesic area and perimeter ***")
                arcpy.AddGeometryAttributes_management(bdiso, "AREA_GEODESIC;PERIMETER_LENGTH_GEODESIC", "METERS",
                                                       "SQUARE_METERS", "")
                
                # create and copy area and perimeter fields
                arcpy.AddField_management(bdiso, "POLY_AREA", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
                arcpy.AddField_management(bdiso, "PERIMETER", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
                arcpy.CalculateField_management(bdiso, "POLY_AREA", "!AREA_GEO!", "PYTHON_9.3", "")
                arcpy.CalculateField_management(bdiso, "PERIMETER", "!PERIM_GEO!", "PYTHON_9.3", "")

            # Select
            buildings_sel = "in_memory\\building_area50" + oid
            if "FOOT" in linear_unit:
                # xy linear units are foot or foot_us
                arcpy.Select_analysis(bdiso, buildings_sel, "\"POLY_AREA\" >= 538.19")
            else:
                # xy linear units are meter
                arcpy.Select_analysis(bdiso, buildings_sel, "\"POLY_AREA\" >= 50")

            # Add Field
            if len(arcpy.ListFields(buildings_sel, "ratio")) < 1:
                arcpy.AddField_management(buildings_sel, "ratio", "DOUBLE", "", "", "", "",
                                          "NULLABLE", "NON_REQUIRED", "")
            # Calculate Field
            arcpy.CalculateField_management(buildings_sel, "ratio", "!POLY_AREA! / !PERIMETER!",
                                            "PYTHON_9.3", "")
            # Select
            # bd_075 = "in_memory\\bd_075"
            bd_075 = os.path.join(scratch_file_gdb_path, "bd_075" + oid)
            arcpy.Select_analysis(buildings_sel, bd_075, "\"ratio\" >=0.75")
            arcpy.Delete_management(buildings_sel)

            # Aggregate Polygons
            # bdagg_tbl = "in_memory\\tbl"
            bdagg_tbl = os.path.join(scratch_file_gdb_path, "bdagg_tbl" + oid)
            arcpy.AggregatePolygons_cartography(bd_075, buildings_output, "1.5 Meters", "50 SquareMeters",
                                                "250 SquareMeters", "ORTHOGONAL", "", bdagg_tbl)

            # Repair the building geometries in case they have self intersecting geometries
            arcpy.RepairGeometry_management(buildings_output, "DELETE_NULL")
            # Note: can do zonal stats here but buildings along borders are dissolved later,
            #       which invalidates the gathered stats for those particular buildings
            # zonalstats(buildings_output, diff_ori, results_file_gdb_path)

        # If the user asks for trees, then the process is different if buildings are also found
        if "TREES" in featureextraction.upper():
            isnull_setnu1 = os.path.join(out_subfolder_path, "IsNullSetN" + oid + r".tif")
            if buildingsfound:
                # Set the nodata mapping method environment to promote the value.
                arcpy.env.nodata = "PROMOTION"
                # Is Null
                isnull_isocl1 = os.path.join(out_subfolder_path, "IsNullisoc" + oid + r".tif")
                arcpy.gp.IsNull_sa(isocluster2, isnull_isocl1)
                # Con
                ras_isnull_isocl1 = Raster(isnull_isocl1)
                int_raster_one = ras_isnull_isocl1 / ras_isnull_isocl1
                intRasOne = os.path.join(out_subfolder_path, "intRasOne" + oid + r".tif")
                int_raster_one.save(intRasOne)
                con_isnull_i1 = os.path.join(out_subfolder_path, "ConIsNulli" + oid + r".tif")
                # arcpy.gp.Con_sa(isnull_isocl1, int_raster_one, con_isnull_i1, isocluster2, "VALUE=1")
                arcpy.gp.Con_sa(isnull_isocl1, intRasOne, con_isnull_i1, isocluster2, "VALUE=1")
                arcpy.Delete_management(ras_isnull_isocl1)
                arcpy.Delete_management(int_raster_one)
                # Focal Statistics
                focalst_isoc1 = os.path.join(out_subfolder_path, "FocalStiso" + oid + r".tif")
                # focalst_isoc1 = os.path.join(out_subfolder_path, "FocalStiso" + oid)
                intCon = os.path.join(out_subfolder_path, "intCon" + oid + r".tif")
                intConIsNull = Int(con_isnull_i1)
                intConIsNull.save(intCon)
                # arcpy.gp.FocalStatistics_sa(Int(con_isnull_i1), focalst_isoc1, "Rectangle 3 3 CELL","MAJORITY","DATA")
                arcpy.gp.FocalStatistics_sa(intCon, focalst_isoc1, "Rectangle 3 3 CELL", "MAJORITY", "DATA")
                # Set Null
                setnull_isot = os.path.join(out_subfolder_path, "SetNulliso" + oid + r".tif")
                # setnull_isot = os.path.join(out_subfolder_path, "SetNulliso" + oid)
                arcpy.gp.SetNull_sa(focalst_isoc1, 0, setnull_isot, "VALUE = 1")

                # Is Null
                # arcpy.gp.IsNull_sa(setnull_isot, isnull_setnu1)
                # call the IsNull function instead of the gp tool in order to output to .tif and avoid data projection
                isnull = IsNull(setnull_isot)
                isnull.save(isnull_setnu1)
                # arcpy.gp.Con_sa(setnull_isot, intRasOne, isnull_setnu1, setnull_isot, "VALUE<>0")
                # arcpy.Delete_management(int_raster_one)
            else:
                # arcpy.gp.IsNull_sa(isocluster2, isnull_setnu1)
                # call the IsNull function instead of the gp tool in order to output to .tif and avoid data projection
                isnull = IsNull(isocluster2)
                isnull.save(isnull_setnu1)
                # arcpy.gp.Con_sa(isocluster2, intRasOne, isnull_setnu1, isocluster2, "VALUE<>0")
                # arcpy.Delete_management(int_raster_one)
            # Make Raster Layer
            treenull = "treeNULL" + oid
            arcpy.MakeRasterLayer_management(isnull_setnu1, treenull, "", "", "")

            # Select Layer By Attribute
            arcpy.SelectLayerByAttribute_management(treenull, "NEW_SELECTION", "VALUE=1")

            # Euclidean Distance
            eucdist_make1 = os.path.join(out_subfolder_path, "EucDistMak" + oid + r".tif")
            # arcpy.gp.EucDistance_sa(treenull, eucdist_make1, "", treenull, "")
            eucDist = EucDistance(treenull, "", treenull, "")
            eucDist.save(eucdist_make1)
            arcpy.Delete_management(treenull)

            # Set Null
            setnull_eucd1 = os.path.join(out_subfolder_path, "SetNulEucD" + oid + r".tif")
            if "FOOT" in linear_unit:
                arcpy.gp.SetNull_sa(eucdist_make1, eucdist_make1, setnull_eucd1, "VALUE<2.95")
            else:
                arcpy.gp.SetNull_sa(eucdist_make1, eucdist_make1, setnull_eucd1, "VALUE<0.9")

            setnulldistest = arcpy.Raster(setnull_eucd1)
            if setnulldistest.maximum <= -1:
                arcpy.AddMessage("No Trees were identified in region {0}".format(oid))
                # print "No Trees were identified in region: " + oid
                if delete_intermediate_files == "true":
                    remove_rasters(out_subfolder_path)
                    remove_shapefiles(out_subfolder_path)
                    arcpy.env.workspace = ""   # so we can successfully delete scratch_file_gdb_path
                    remove_filegdb(scratch_file_gdb_path)
                if arcpy.Exists(raster_iso):
                    arcpy.Delete_management(raster_iso)
                del desc, sr
                arcpy.CheckInExtension('3D')
                arcpy.CheckInExtension('Spatial')
                return

            # Focal Statistics
            focalst_setn1 = os.path.join(out_subfolder_path, "FocalStStN" + oid + r".tif")
            arcpy.gp.FocalStatistics_sa(setnull_eucd1, focalst_setn1, "Circle 3 CELL", "MAXIMUM", "DATA")

            # Minus
            minus_focals1 = os.path.join(out_subfolder_path, "MinusFoclS" + oid + r".tif")
            arcpy.gp.Minus_sa(focalst_setn1, setnull_eucd1, minus_focals1)
            arcpy.Delete_management(setnulldistest)  # deleted after done w/ setnull_eucd1 & setnulldistest

            # Equal To
            equalto_minu1 = os.path.join(out_subfolder_path, "EqualToMin" + oid + r".tif")
            arcpy.gp.EqualTo_sa(minus_focals1, "0", equalto_minu1)

            # Set Null
            setnull_equa1 = os.path.join(out_subfolder_path, "SetNulEqua" + oid + r".tif")
            arcpy.gp.SetNull_sa(equalto_minu1, eucdist_make1, setnull_equa1, "VALUE=0")

            setnulleqtst = arcpy.Raster(setnull_equa1)
            if setnulleqtst.maximum <= -1:
                arcpy.AddMessage("No Trees were identified in region {0}".format(oid))
                # print "No Trees were identified in region: " + oid
                if delete_intermediate_files == "true":
                    remove_rasters(out_subfolder_path)
                    remove_shapefiles(out_subfolder_path)
                    arcpy.env.workspace = ""   # so we can successfully delete scratch_file_gdb_path
                    remove_filegdb(scratch_file_gdb_path)
                if arcpy.Exists(raster_iso):
                    arcpy.Delete_management(raster_iso)
                del desc, sr
                arcpy.CheckInExtension('3D')
                arcpy.CheckInExtension('Spatial')
                return

            # Plus
            # plus_int_set1 = os.path.join(out_subfolder_path, "Plus_Int_Set1" + oid)
            flt_meter_elev_unit = float(elevation_meter_scalefactor)
            plus_int_set1 = Raster(setnull_equa1) + flt_meter_elev_unit
            arcpy.Delete_management(flt_meter_elev_unit)

            # Int
            int_setnull_1 = os.path.join(out_subfolder_path, "IntSetNull" + oid + r".tif")
            arcpy.gp.Int_sa(plus_int_set1, int_setnull_1)
            arcpy.Delete_management(plus_int_set1)
            arcpy.Delete_management(setnulleqtst)  # deleted after done with setnull_equal & setnulleqtst

            #  Clean Tree Raster
            # NOTE: Don't output the CleanTreeR raster as .tif or sometimes Gridcode is set to 255, resulting
            #       in very large trees.
            cleantreeras1 = os.path.join(out_subfolder_path, "CleanTreeR" + oid)
            cleanfast(int_setnull_1, cleantreeras1)

            # Raster to Polygon
            rastert_setnull1 = os.path.join(scratch_file_gdb_path, "RastrTSetN" + oid)
            arcpy.RasterToPolygon_conversion(cleantreeras1, rastert_setnull1, "NO_SIMPLIFY", "Value")
            arcpy.Delete_management(cleantreeras1)

            # Feature To Point
            rastert_setnull1_featuretopo = os.path.join(scratch_file_gdb_path,
                                                        "RasterT_SetNull1_FeatureToPo" + oid)
            arcpy.FeatureToPoint_management(rastert_setnull1, rastert_setnull1_featuretopo, "INSIDE")

            # Make Feature Layer of gridcode > 0 (to eliminate erroneous features)
            rastert_setnull1_featuretopolyr = "RasterT_SetNull1_FeatureToPoLyr" + oid
            arcpy.MakeFeatureLayer_management(rastert_setnull1_featuretopo, rastert_setnull1_featuretopolyr,
                                              "GRIDCODE > 0", "", "")

            if buildingsfound:
                # Make Feature Layer moved to above (out of If statement, since it applies in either case)
                # rastert_setnull1_featuretopolyr = "RasterT_SetNull1_FeatureToPoLyr" + oid
                # arcpy.MakeFeatureLayer_management(rastert_setnull1_featuretopo, rastert_setnull1_featuretopolyr,
                #                                   "", "", "")

                # Copy Features
                bdagg_copyfeatures = os.path.join(scratch_file_gdb_path, "bdAgg_CopyFeatures" + oid)
                arcpy.CopyFeatures_management(buildings_output, bdagg_copyfeatures, "", "0", "0", "0")

                # Select Layer By Location
                arcpy.SelectLayerByLocation_management(rastert_setnull1_featuretopolyr, "WITHIN_A_DISTANCE",
                                                       bdagg_copyfeatures, "3.3 Feet", "NEW_SELECTION",
                                                       "NOT_INVERT")

                # Select Layer By Attribute
                arcpy.SelectLayerByAttribute_management(rastert_setnull1_featuretopolyr, "SWITCH_SELECTION", "")

                # Buffer
                arcpy.Buffer_analysis(rastert_setnull1_featuretopolyr, trees_output,
                                      "GRIDCODE", "FULL", "ROUND", "NONE", "", "PLANAR")
                # arcpy.AddMessage("Trees exported to:     {0}".format(trees_output))

            else:
                # Buffer
                # arcpy.Buffer_analysis(rastert_setnull1_featuretopo, trees_output, "GRIDCODE", "FULL", "ROUND",
                #                       "NONE", "", "PLANAR")
                arcpy.Buffer_analysis(rastert_setnull1_featuretopolyr, trees_output, "GRIDCODE", "FULL", "ROUND",
                                      "NONE", "", "PLANAR")
                # arcpy.AddMessage("Trees exported to:     {0}".format(trees_output))

            arcpy.Delete_management(rastert_setnull1_featuretopolyr)
        if arcpy.Exists(raster_iso):
            arcpy.Delete_management(raster_iso)
        arcpy.Delete_management("in_memory")
        # Clean up the rasters in out_subfolder_path
        if delete_intermediate_files == "true":
            remove_rasters(out_subfolder_path)
            remove_shapefiles(out_subfolder_path)
            arcpy.env.workspace = ""   # so we can successfully delete scratch_file_gdb_path
            remove_filegdb(scratch_file_gdb_path)
            # arcpy.env.workspace = out_subfolder_path
        del desc, sr
        end = time.time()
        delta = end - start
        arcpy.AddMessage("Elapsed time for OID {0} is {1} seconds".format(oid, delta))
        arcpy.CheckInExtension('3D')
        arcpy.CheckInExtension('Spatial')
        return

    except:
        # Get the traceback object
        #
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))
        raise


def remove_shapefiles(directory):
    # Delete all of the shapefiles in directory
    arcpy.env.workspace = directory
    featureclasses = arcpy.ListFeatureClasses()
    for featureclass in featureclasses:
        arcpy.Delete_management(featureclass)
    arcpy.ClearEnvironment("workspace")
    return


def remove_rasters(directory):
    # Delete all of the rasters in directory
    arcpy.env.workspace = directory
    rasters = arcpy.ListRasters("*", "ALL")
    for raster in rasters:
        arcpy.Delete_management(raster)
    arcpy.ClearEnvironment("workspace")
    return


def remove_filegdb(filegdb):
    # Delete all of the feature classes in filegdb, then delete filegdb
    arcpy.env.workspace = filegdb
    featureclasses = arcpy.ListFeatureClasses()
    for featureclass in featureclasses:
        arcpy.Delete_management(featureclass)
    arcpy.env.workspace = ""
    try:
        arcpy.Delete_management(filegdb)
    except:
        arcpy.AddMessage("Unable to delete file GDB: {0}".format(filegdb))
    # This extra code is needed because of a bug which results in the deleted file gdb being
    # changed into a folder (i.e. loses it's designation as a file gdb, but the folder still exists)
    if arcpy.Exists(filegdb):
        try:
            arcpy.Delete_management(filegdb)
        except:
            arcpy.AddMessage("Unable to delete file GDB the second time: {0}".format(filegdb))
    arcpy.ClearEnvironment("workspace")
    return



def remove_tables(filegdb):
    # Remove all of the tables in filegdb
    arcpy.env.workspace = filegdb
    tables = arcpy.ListTables()
    for table in tables:
        arcpy.Delete_management(table)
    arcpy.ClearEnvironment("workspace")
    return


def delete(array, i, j, cellsize):
    try:
        # This cleans overlapping trees
        ext = array[i][j]
        for x in range(i - ext, i + ext):
            for y in range(j - ext, j + ext):
                if x > 0 and x < (array.shape[0] - 1) and y > 0 and y < (array.shape[1] - 1) and (x != i or y != j):
                    r = array[x][y]
                    if r > 0 and ext >= r:
                        if (float(r) / float(ext)) > 0.5:
                            distance = math.sqrt(math.pow(x - i, 2) + math.pow(y - j, 2)) * cellsize
                            if distance < (1.0 * float(ext)):  # change threshold here
                                array[x][y] = 0
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))
        raise


def cleanfast(inputraster, outputraster):
    try:
        # get raster information
        desc = arcpy.Describe(inputraster)
        sr = desc.spatialReference
        point = arcpy.Point(desc.Extent.XMin, desc.Extent.YMin)

        # iterate through raster array
        dist = arcpy.RasterToNumPyArray(in_raster=inputraster, nodata_to_value=0)
        locs = dist.nonzero()
        # part = int(len(locs[0]) / 10)

        for x in range(0, len(locs[0])):
            # if (int(x) % part) == 0:
            #     arcpy.AddMessage(str(float(x) / float(part) * 10) + "% completed")
            locx = locs[0][x]
            locy = locs[1][x]
            delete(dist, locx, locy, desc.meanCellWidth)

        # output
        distraster = arcpy.NumPyArrayToRaster(in_array=dist, lower_left_corner=point,
                                              x_cell_size=desc.meanCellWidth,
                                              y_cell_size=desc.meanCellWidth,
                                              value_to_nodata=0)
        arcpy.DefineProjection_management(distraster, sr)
        distraster.save(outputraster)
        del desc, sr, point
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))
        raise


def zonalstats(vectorfc, rasterheights, resultstablegdb):
    # Gather zonal statistics for the features in vectorfc
    try:
        start = time.time()
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            arcpy.AddMessage("\n*** No Spatial Analyst License Available.")
            return
        table_prefix = os.path.splitext(os.path.basename(vectorfc))[0]
        zonaltable = os.path.join(resultstablegdb, table_prefix + r"_zonalTbl")
        oid_fieldname = arcpy.Describe(vectorfc).OIDFieldName
        arcpy.gp.ZonalStatisticsAsTable_sa(vectorfc, oid_fieldname, rasterheights, zonaltable, "DATA", "ALL")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from ZonalStatisticsAsTable_sa are: \n{0}\n".format(messages))
        arcpy.MakeFeatureLayer_management(vectorfc, "vectorfc_layer")
        arcpy.JoinField_management("vectorfc_layer", oid_fieldname, zonaltable, join_field="OBJECTID",
                                   fields="COUNT;AREA;MIN;MAX;RANGE;MEAN;STD;SUM")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from JoinField are: \n{0}\n".format(messages))
        arcpy.Delete_management("vectorfc_layer")
        end = time.time()
        delta = end - start
        arcpy.AddMessage("Elapsed time for ZonalStats and JoinField is {0} seconds".format(delta))
        arcpy.CheckInExtension('Spatial')
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))
        raise


def create_md_from_raster(raster_folder, file_gdb, md_name, md_boundary, build_footprints, export_boundary):
    arcpy.env.overwriteOutput = True
    # Create and populate a mosaic dataset with elevation rasters (32-bit float)
    try:
        arcpy.env.workspace = raster_folder
        if not arcpy.Exists(file_gdb):
            arcpy.AddMessage("\n*** Exiting create_md_from_raster...File GDB Does not exist: {0} ***".format(file_gdb))
            return

        if not os.path.exists(raster_folder):
            arcpy.AddMessage("\n*** Exiting create_md_from_raster..."
                             "Raster Folder Does not exist: {0} ***".format(raster_folder))
            return

        full_md_path = os.path.join(file_gdb, md_name)
        arcpy.AddMessage("\nMD to be created:  {0}".format(full_md_path))

        # Don't re-create the Mosaic Dataset if it already exists
        if not arcpy.Exists(full_md_path):
            # Get the spatial reference string of the first raster (to use in creation of MD)
            rasters = arcpy.ListRasters("*", "All")
            # Make sure there's at least one raster in raster_folder
            # If not, then exit the script
            # If so, get the raster's Spatial Reference
            if len(rasters) > 0:
                desc_firstraster = arcpy.Describe(rasters[0])
                spatref_firstraster = desc_firstraster.SpatialReference.exportToString()
                arcpy.AddMessage("Spatial ref of 1st raster in {0} is: \n{1}\n".format(raster_folder,
                                                                                       spatref_firstraster))
                arcpy.AddMessage("Number of rasters in {0}:  {1}".format(raster_folder, len(rasters)))
            else:
                arcpy.AddMessage("\n*** Exiting create_md_from_raster..."
                                 "No rasters found in {0} ***".format(raster_folder))
                return
            # Create a Mosaic Dataset
            arcpy.CreateMosaicDataset_management(file_gdb, md_name,
                                                 coordinate_system=spatref_firstraster,
                                                 num_bands="1", pixel_type="32_BIT_FLOAT", product_definition="NONE",
                                                 product_band_definitions="#")
            del desc_firstraster, spatref_firstraster, rasters
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from CreateMosaicDataset are: \n{0}\n".format(messages))

            # set the data_type to ELEVATION
            arcpy.SetRasterProperties_management(full_md_path, data_type="ELEVATION", statistics="",
                                                 stats_file="#", nodata="")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from SetRasterProperties are: \n{0}\n".format(messages))

            # Add rasters from Raster folder to MD
            arcpy.AddRastersToMosaicDataset_management(full_md_path, raster_type="Raster Dataset",
                                                       input_path=raster_folder,
                                                       update_cellsize_ranges="UPDATE_CELL_SIZES",
                                                       update_boundary="UPDATE_BOUNDARY",
                                                       update_overviews="NO_OVERVIEWS", maximum_pyramid_levels="",
                                                       maximum_cell_size="0",
                                                       minimum_dimension="1500", spatial_reference="", filter="",
                                                       sub_folder="SUBFOLDERS",
                                                       duplicate_items_action="ALLOW_DUPLICATES",
                                                       build_pyramids="NO_PYRAMIDS",
                                                       calculate_statistics="NO_STATISTICS",
                                                       build_thumbnails="NO_THUMBNAILS",
                                                       operation_description="#",
                                                       force_spatial_reference="NO_FORCE_SPATIAL_REFERENCE")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from AddRastersToMosaicDataset are: \n{0}\n".format(messages))

            # re-calculate cell size ranges so export will work at various scales
            arcpy.CalculateCellSizeRanges_management(full_md_path, where_clause="", do_compute_min="MIN_CELL_SIZES",
                                                     do_compute_max="MAX_CELL_SIZES", max_range_factor="100",
                                                     cell_size_tolerance_factor="0.8", update_missing_only="UPDATE_ALL")

            if build_footprints == "true":
                arcpy.BuildFootprints_management(full_md_path, "", "RADIOMETRY", "-10", "4294967295", "80", "0",
                                                 "NO_MAINTAIN_EDGES", "SKIP_DERIVED_IMAGES", "UPDATE_BOUNDARY",
                                                 "2000", "100", "NONE", "", "20", "0.05")
                messages = arcpy.GetMessages()
                arcpy.AddMessage("\nResults output from BuildFootprints are: \n{0}\n".format(messages))
        else:
            arcpy.AddMessage("\n*** Mosaic Dataset already exists: {0} ***".format(full_md_path))

        # default_compression_type="LERC"
        # clip_to_footprints="NOT_CLIP"
        # data_source_type="ELEVATION"
        # rows_maximum_imagesize="15000"
        arcpy.SetMosaicDatasetProperties_management(full_md_path, "15000", "15000", "None;JPEG;LZ77;LERC", "LERC",
                                                    "75", "0.01", "BILINEAR", "NOT_CLIP",
                                                    "FOOTPRINTS_MAY_CONTAIN_NODATA",
                                                    "NOT_CLIP", "NOT_APPLY", "#", "NONE",
                                                    "NorthWest;Center;LockRaster;ByAttribute;Nadir;Viewpoint;"
                                                    "Seamline;None",
                                                    "NorthWest", "", "", "ASCENDING", "FIRST", "10", "600", "300", "20",
                                                    "0.8", "", "BASIC",
                                                    "Name;MinPS;MaxPS;LowPS;HighPS;Tag;GroupName;ProductName;"
                                                    "CenterX;CenterY;ZOrder;Shape_Length;Shape_Area", "DISABLED", "",
                                                    "", "", "", "20", "1000", "ELEVATION", "1", "None", "None")

        # Get a record count just to be sure we found raster products to ingest
        result = arcpy.GetCount_management(full_md_path)
        count_rasters = int(result.getOutput(0))

        if count_rasters == 0:
            arcpy.AddMessage("\n*** Exiting: {0} Mosaic Dataset has no raster products ***".format(full_md_path))
            sys.exit()
        else:
            arcpy.AddMessage("{0} has {1} raster product(s).".format(full_md_path, count_rasters))

        # boundary = os.path.join(file_gdb, md_boundary)
        if export_boundary == "true":
            if not arcpy.Exists(md_boundary):
                # Export Boundary to the file GDB which holds the final results
                arcpy.ExportMosaicDatasetGeometry_management(full_md_path, md_boundary, "", "BOUNDARY")
                messages = arcpy.GetMessages()
                arcpy.AddMessage("Results output from ExportMosaicDatasetGeometry are: \n{0}\n".format(messages))
            else:
                arcpy.AddMessage("Exported boundary already exists: {}".format(md_boundary))
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))
    finally:
        arcpy.ClearEnvironment("workspace")


def create_md_from_las(lasd, raster_type_file, file_gdb, md_name, md_boundary, build_footprints,
                       export_boundary):

    arcpy.env.overwriteOutput = True
    # Create and populate a mosaic dataset with a LAS dataset (32-bit float)
    try:
        # arcpy.env.workspace = raster_folder
        if not arcpy.Exists(file_gdb):
            arcpy.AddMessage("\n*** Exiting create_md_from_las...File GDB Does not exist: {0} ***".format(file_gdb))
            return
        if not os.path.exists(lasd):
            arcpy.AddMessage("\n*** Exiting create_md_from_las...LAS dataset Does not exist: {0} ***".format(lasd))
            return

        full_md_path = os.path.join(file_gdb, md_name)
        arcpy.AddMessage("\nMD to be created:  {0}".format(full_md_path))
        # md_boundary = full_md_path + boundary_append

        # Don't re-create the Mosaic Dataset if it already exists
        if not arcpy.Exists(full_md_path):
            # Get the spatial reference string of the LAS Dataset (to use in creation of MD)
            desc_lasd = arcpy.Describe(lasd)
            spat_ref_lasd = desc_lasd.SpatialReference
            # Create a Mosaic Dataset
            arcpy.CreateMosaicDataset_management(file_gdb, md_name,
                                                 coordinate_system=spat_ref_lasd,
                                                 num_bands="1", pixel_type="32_BIT_FLOAT", product_definition="NONE",
                                                 product_band_definitions="#")
            del desc_lasd, spat_ref_lasd
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from CreateMosaicDataset are: \n{0}\n".format(messages))

            # set the NoData value to -3.40282346639e+038
            arcpy.SetRasterProperties_management(full_md_path, data_type="ELEVATION", statistics="",
                                                 stats_file="#", nodata="")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from SetRasterProperties are: \n{0}\n".format(messages))

            # Add rasters from Raster folder to MD
            arcpy.AddRastersToMosaicDataset_management(full_md_path, raster_type=raster_type_file,
                                                       input_path=lasd,
                                                       update_cellsize_ranges="UPDATE_CELL_SIZES",
                                                       update_boundary="UPDATE_BOUNDARY",
                                                       update_overviews="NO_OVERVIEWS", maximum_pyramid_levels="",
                                                       maximum_cell_size="0",
                                                       minimum_dimension="1500", spatial_reference="", filter="#",
                                                       sub_folder="SUBFOLDERS",
                                                       duplicate_items_action="ALLOW_DUPLICATES",
                                                       build_pyramids="NO_PYRAMIDS",
                                                       calculate_statistics="NO_STATISTICS",
                                                       build_thumbnails="NO_THUMBNAILS",
                                                       operation_description="#",
                                                       force_spatial_reference="NO_FORCE_SPATIAL_REFERENCE")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from AddRastersToMosaicDataset are: \n{0}\n".format(messages))

            if build_footprints == "true":
                arcpy.BuildFootprints_management(full_md_path, "", "RADIOMETRY", "-100", "4294967295", "300", "0",
                                                 "MAINTAIN_EDGES", "SKIP_DERIVED_IMAGES", "UPDATE_BOUNDARY",
                                                 "2000", "20", "NONE", "", "20", "0.05")
                messages = arcpy.GetMessages()
                arcpy.AddMessage("\nResults output from BuildFootprints are: \n{0}\n".format(messages))
        else:
            arcpy.AddMessage("\n*** Mosaic Dataset already exists: {0} ***".format(full_md_path))

        # Set the desired MD properties (non-default parameters are listed below):
        #   default_compression_type="LERC"
        #   clip_to_boundary="CLIP"
        #   data_source_type="ELEVATION"
        #   rows_maximum_imagesize="25000"
        #   columns_maximum_imagesize="25000"
        arcpy.SetMosaicDatasetProperties_management(full_md_path, rows_maximum_imagesize="25000",
                                                    columns_maximum_imagesize="25000",
                                                    allowed_compressions="None;JPEG;LZ77;LERC",
                                                    default_compression_type="LERC", JPEG_quality="75",
                                                    LERC_Tolerance="0.01", resampling_type="BILINEAR",
                                                    clip_to_footprints="NOT_CLIP",
                                                    footprints_may_contain_nodata="FOOTPRINTS_MAY_CONTAIN_NODATA",
                                                    clip_to_boundary="CLIP",
                                                    color_correction="NOT_APPLY",
                                                    allowed_mensuration_capabilities="Basic",
                                                    default_mensuration_capabilities="Basic",
                                                    allowed_mosaic_methods="NorthWest;Center;LockRaster;ByAttribute;"
                                                                           "Nadir;Viewpoint;Seamline;None",
                                                    default_mosaic_method="NorthWest", order_field="", order_base="#",
                                                    sorting_order="ASCENDING", mosaic_operator="FIRST", blend_width="0",
                                                    view_point_x="600", view_point_y="300", max_num_per_mosaic="20",
                                                    cell_size_tolerance="0.8", cell_size="#", metadata_level="BASIC",
                                                    transmission_fields="",
                                                    use_time="DISABLED", start_time_field="", end_time_field="#",
                                                    time_format="#", geographic_transform="#",
                                                    max_num_of_download_items="20", max_num_of_records_returned="1000",
                                                    data_source_type="ELEVATION", minimum_pixel_contribution="1")

        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from SetMosaicDatasetProperties are: \n{0}\n".format(messages))

        # Get a record count just to be sure we found raster products to ingest
        result = arcpy.GetCount_management(full_md_path)
        count_lasd = int(result.getOutput(0))

        if count_lasd == 0:
            arcpy.AddMessage("\n*** Exiting: {0} Mosaic Dataset has no LASD contents ***".format(full_md_path))
            sys.exit()
        else:
            arcpy.AddMessage("{0} has {1} LASD(s).".format(full_md_path, count_lasd))

        # boundary = os.path.join(file_gdb, md_boundary)
        if export_boundary == "true":
            if not arcpy.Exists(md_boundary):
                # Export Boundary to the file GDB which holds the final results
                arcpy.ExportMosaicDatasetGeometry_management(full_md_path, md_boundary, "", "BOUNDARY")
                messages = arcpy.GetMessages()
                arcpy.AddMessage("Results output from ExportMosaicDatasetGeometry are: \n{0}\n".format(messages))
            else:
                arcpy.AddMessage("Exported boundary already exists: {}".format(md_boundary))
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))


def regularize_some_buildings(building, bd_to_reg_shp, rastert_mask_mb1_symdiff_shp, expression, cellsize, scratch_ws):
    # This routine will regularize smaller buildings that are oriented North-South
    # The remainder of the buildings will be placed in bd_to_reg_shp for regularization by the 'Regularize Building
    # Footprint' gp tool in ArcGIS Pro
    try:
        loc = building.rfind(".gdb")
        bdfilename = building[loc + 5:]

        buildings_lyr = "buildings"
        # ratio_vr = Buildings_output__2_
        # ratio_occ = ratio_vr
        # ori_vex = ratio_occ
        bdmbrv_shp = scratch_ws + "\\" + "bdMBRv.shp"
        # MBR_vex = bdMBRv_shp
        # all = ori_vex
        bdmbrr_shp = scratch_ws + "\\" + "bdMBRr.shp"
        # MBR = bdMBRr_shp
        # building1 = all
        # building2 = building1
        bd_mbr_output_shp = scratch_ws + "\\" + "bd_MBR_output.shp"
        bd_mbr_final_shp = scratch_ws + "\\" + "bd_MBR_final.shp"
        bd_mbr_final_polygontoraster_tif = scratch_ws + "\\" + "bd_MBR_final_PolygonToRaster.tif"
        # Input_raster_or_constant_value_2__2_ = "-1"
        lesstha_buil1 = scratch_ws + "\\" + "LessTha_buil1"
        # Input_raster_or_constant_value_2__3_ = "1"
        mask_mbr = scratch_ws + "\\" + "mask_mbr"
        rastert_mask_mb1_shp = scratch_ws + "\\" + "RasterT_mask_mb1.shp"
        buildings_polygontoraster_tif = scratch_ws + "\\" + "buildings_PolygonToRaster.tif"
        # Input_raster_or_constant_value_2 = "-1"
        greater_bd_m1 = scratch_ws + "\\" + "Greater_bd_M1"
        isnull_great1 = scratch_ws + "\\" + "IsNull_Great1"
        plus_isnull_1 = scratch_ws + "\\" + "Plus_IsNull_1"
        # Input_true_raster_or_constant_value = "1"
        # Input_false_raster_or_constant_value__2_ = "0"
        con_plus_isn1 = scratch_ws + "\\" + "Con_Plus_IsN1"
        isnull_con_p1 = scratch_ws + "\\" + "IsNull_Con_P1"
        # Input_false_raster_or_constant_value__3_ = "0"
        con_plus_isn2 = scratch_ws + "\\" + "Con_Plus_IsN2"
        shrink_plus_1 = scratch_ws + "\\" + "Shrink_Plus_1"
        times_shrink1 = scratch_ws + "\\" + "Times_Shrink1"
        # Input_false_raster_or_constant_value = "0"
        setnull_time1 = scratch_ws + "\\" + "SetNull_Time1"
        rastert_setnull1_shp = scratch_ws + "\\" + "RasterT_SetNull1.shp"
        rastert_setnull1_minimumboun_shp = scratch_ws + "\\" + "RasterT_SetNull1_MinimumBoun.shp"
        mbr_sel_shp = scratch_ws + "\\" + "mbr_sel.shp"
        # RasterT_SetNull1_MinimumBoun = "mbr_sel_Layer"
        mbr_sel_layer = "mbr_sel_Layer"
        # RasterT_SetNull1_MinimumBoun__2_ = RasterT_SetNull1_MinimumBoun
        # RasterT_SetNull1_MinimumBoun__3_ = RasterT_SetNull1_MinimumBoun__2_

        arcpy.MakeFeatureLayer_management(building, buildings_lyr)
        # Process: Add Geometry Attributes
        arcpy.AddGeometryAttributes_management(building, "AREA", "", "SQUARE_METERS", "")

        # Process: Add Field
        if len(arcpy.ListFields(buildings_lyr, "ratio_vr")) < 1:
            arcpy.AddField_management(buildings_lyr, "ratio_vr", "DOUBLE", "", "", "", "",
                                      "NULLABLE", "NON_REQUIRED", "")

        # Process: Add Field
        if len(arcpy.ListFields(buildings_lyr, "ratio_occ")) < 1:
            arcpy.AddField_management(buildings_lyr, "ratio_occ", "DOUBLE", "", "", "", "",
                                      "NULLABLE", "NON_REQUIRED", "")

        # Process: Minimum Bounding Geometry
        arcpy.MinimumBoundingGeometry_management(building, bdmbrv_shp, "CONVEX_HULL", "NONE", "", "NO_MBG_FIELDS")

        # Process: Add Geometry Attributes
        arcpy.AddGeometryAttributes_management(bdmbrv_shp, "AREA", "", "SQUARE_METERS", "")

        # Process: Add Join
        arcpy.AddJoin_management(buildings_lyr, "OBJECTID", bdmbrv_shp, "ORIG_FID", "KEEP_ALL")

        # Process: Minimum Bounding Geometry
        arcpy.MinimumBoundingGeometry_management(building, bdmbrr_shp, "ENVELOPE", "NONE", "", "NO_MBG_FIELDS")

        # Process: Add Geometry Attributes
        arcpy.AddGeometryAttributes_management(bdmbrr_shp, "AREA", "", "SQUARE_METERS", "")

        # Process: Add Join
        arcpy.AddJoin_management(buildings_lyr, "bdMBRv.ORIG_FID", bdmbrr_shp, "ORIG_FID", "KEEP_ALL")

        # Process: Calculate Field
        arcpy.CalculateField_management(buildings_lyr, "ratio_vr", "!bdMBRv.POLY_AREA! / !bdMBRr.POLY_AREA!",
                                        "PYTHON_9.3", "")

        # Process: Calculate Field
        arcpy.CalculateField_management(buildings_lyr, "ratio_occ",
                                        "!" + bdfilename + ".POLY_AREA! / !bdMBRr.POLY_AREA!", "PYTHON_9.3", "")

        # Process: Select
        arcpy.Select_analysis(buildings_lyr, bd_to_reg_shp,
                              "\"" + bdfilename + ".POLY_AREA\" >= 500 OR \"" + bdfilename + ".ratio_vr\"<0.70")

        # Process: Select
        arcpy.Select_analysis(buildings_lyr, bd_mbr_output_shp,
                              "\"" + bdfilename + ".POLY_AREA\" <= 500 AND \"" + bdfilename + ".ratio_vr\">=0.70")

        # Process: Minimum Bounding Geometry
        arcpy.MinimumBoundingGeometry_management(bd_mbr_output_shp, bd_mbr_final_shp, "ENVELOPE", "NONE", "",
                                                 "NO_MBG_FIELDS")

        # Process: Polygon to Raster
        arcpy.PolygonToRaster_conversion(bd_mbr_final_shp, "bdMBRv_rat", bd_mbr_final_polygontoraster_tif,
                                         "CELL_CENTER", "NONE", cellsize)

        # Process: Less Than
        arcpy.gp.LessThan_sa(bd_mbr_final_polygontoraster_tif, "-1", lesstha_buil1)

        # Process: Plus
        # arcpy.gp.Plus_sa(LessTha_buil1, "1", mask_mbr)
        lesstha_buil1_ras = Raster(lesstha_buil1)
        int_raster_one = CreateConstantRaster(1, "INTEGER", lesstha_buil1_ras.meanCellWidth, lesstha_buil1_ras.extent)
        arcpy.gp.Plus_sa(lesstha_buil1, int_raster_one, mask_mbr)

        # Process: Raster to Polygon
        arcpy.RasterToPolygon_conversion(mask_mbr, rastert_mask_mb1_shp, "NO_SIMPLIFY", "VALUE")

        # Process: Extra Add Field (since building is now a feature class instead of a shapefile)
        if len(arcpy.ListFields(building, "Id")) < 1:
            arcpy.AddField_management(building, "Id", "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED", "")
        arcpy.CalculateField_management(building, "Id", "0", "PYTHON_9.3", "")

        # Process: Polygon to Raster
        arcpy.PolygonToRaster_conversion(building, "Id", buildings_polygontoraster_tif, "CELL_CENTER", "NONE",
                                         cellsize)

        # Process: Greater Than
        arcpy.gp.GreaterThan_sa(buildings_polygontoraster_tif, "-1", greater_bd_m1)

        # Process: Is Null
        arcpy.gp.IsNull_sa(greater_bd_m1, isnull_great1)

        # Process: Plus
        arcpy.gp.Plus_sa(isnull_great1, lesstha_buil1, plus_isnull_1)
        arcpy.Delete_management(lesstha_buil1_ras)  # Delete after we're done with lesstha_buil1_ras and lesstha_buil1

        # Process: Con
        arcpy.gp.Con_sa(plus_isnull_1, int_raster_one, con_plus_isn1, "0", "VALUE=0")
        arcpy.Delete_management(int_raster_one)

        # Process: Is Null
        arcpy.gp.IsNull_sa(con_plus_isn1, isnull_con_p1)

        # Process: Con
        arcpy.gp.Con_sa(isnull_con_p1, con_plus_isn1, con_plus_isn2, "0", "VALUE=0")

        # Process: Shrink
        arcpy.gp.Shrink_sa(con_plus_isn2, shrink_plus_1, "3", "0")

        # Process: Times
        arcpy.gp.Times_sa(shrink_plus_1, mask_mbr, times_shrink1)

        # Process: Set Null
        arcpy.gp.SetNull_sa(times_shrink1, "0", setnull_time1, "VALUE=1")

        # Process: Raster to Polygon
        arcpy.RasterToPolygon_conversion(setnull_time1, rastert_setnull1_shp, "NO_SIMPLIFY", "VALUE")

        # Process: Minimum Bounding Geometry
        arcpy.MinimumBoundingGeometry_management(rastert_setnull1_shp, rastert_setnull1_minimumboun_shp, "ENVELOPE",
                                                 "NONE", "", "MBG_FIELDS")

        # Process: Select
        arcpy.Select_analysis(rastert_setnull1_minimumboun_shp, mbr_sel_shp, expression)

        # Process: Make Feature Layer
        arcpy.MakeFeatureLayer_management(mbr_sel_shp, mbr_sel_layer, "", "", "")

        # Process: Select Layer By Location
        arcpy.SelectLayerByLocation_management(mbr_sel_layer, "COMPLETELY_WITHIN", rastert_mask_mb1_shp, "",
                                               "NEW_SELECTION", "NOT_INVERT")

        # Process: Select Layer By Attribute
        arcpy.SelectLayerByAttribute_management(mbr_sel_layer, "SWITCH_SELECTION", "")

        # Process: Symmetrical Difference
        arcpy.SymDiff_analysis(rastert_mask_mb1_shp, mbr_sel_layer, rastert_mask_mb1_symdiff_shp, "ALL", "")

        arcpy.Delete_management(buildings_lyr)
        arcpy.Delete_management(mbr_sel_layer)
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))


def set_pixel_size(raster_type_file, las_point_spacing, sr_units):
    # The template raster type file has PUT_PIXEL_SIZE_HERE in place of
    # the pixel size. This is because each set of LAS has a unique point
    # spacing. This module will set pixel size in the art.xml so the
    # resulting Mosaic Dataset will have an appropriate pixel size.
    # Note: the Mosaic Dataset pixel size can't be too small or footprints aren't
    # generated properly.
    try:
        # Get a good number for pixel size of the MD. It can't be too small or
        # footprints either won't be generated at all or they will be incorrect.
        if "METER" in sr_units:
            md_pixel_size = max(3.0, round(2 * las_point_spacing + 0.5))
        else:
            # assumed that xy units are either FOOT or US_FOOT
            md_pixel_size = max(10.0, round(6.5 * las_point_spacing))

        arcpy.AddMessage("Mosaic Dataset pixel size will be: {0}".format(md_pixel_size))
        search_text = r"PUT_PIXEL_SIZE_HERE"
        # Read in the file
        filedata = None
        with open(raster_type_file, 'r') as file:
            filedata = file.read()
        # Replace the target string
        filedata = filedata.replace(search_text, str(md_pixel_size))
        # Write the file out again
        with open(raster_type_file, 'w') as file:
            file.write(filedata)
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))


def create_lasd(las_path, lasd):
    try:
        arcpy.env.workspace = las_path
        lasfiles = arcpy.ListFiles("*.las")
        arcpy.AddMessage("Entire LAS list in {0} is: \n{1}\n".format(las_path, lasfiles))
        if len(lasfiles) > 0:
            arcpy.AddMessage("Creating LAS Dataset: {0}".format(lasd))
            # Create a LAS Dataset and add the LAS files in las_path
            # Compute stats (lasx) if they don't already exist
            lasfiles_list = [os.path.join(las_path, fi) for fi in lasfiles]
            arcpy.CreateLasDataset_management(lasfiles_list, lasd, folder_recursion="NO_RECURSION",
                                              in_surface_constraints="#", spatial_reference="#",
                                              compute_stats="COMPUTE_STATS", relative_paths="ABSOLUTE_PATHS")
                                              # compute_stats="COMPUTE_STATS", relative_paths="RELATIVE_PATHS")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from CreateLasDataset are: \n{0}\n".format(messages))
        del lasfiles
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))
    finally:
        arcpy.ClearEnvironment("workspace")


def get_las_point_spacing(las_path, lasd, output_path, output_gdb_name, lasdatasetstatstext):
    # This module determines the point spacing of the LAS files, and can
    # be used to determine a reasonable raster product pixel size.
    try:
        if arcpy.Exists(lasd):
            arcpy.AddMessage("Calculating point spacing of LASD: {0}".format(lasd))
            # lasdatasetstatstext = os.path.join(output_path, "lasDatasetStatsText.txt")
            if not arcpy.Exists(lasdatasetstatstext):
                arcpy.Delete_management(lasdatasetstatstext)
                arcpy.LasDatasetStatistics_management(lasd, "true", lasdatasetstatstext, "LAS_FILES", "COMMA",
                                                      "DECIMAL_POINT")
            else:
                arcpy.AddMessage("lasDatasetStatsText already exists: {0}".format(lasdatasetstatstext))

            ptFileInfoFC = os.path.join(output_gdb_name, 'ptFileInfoFC')
            if not arcpy.Exists(ptFileInfoFC):
                # Note: This step is optional, so if it takes too long it's safe to remove it
                # get lasd sr
                descLASD = arcpy.Describe(lasd)
                SpatRefLASD = descLASD.SpatialReference
                # SpatRefStringLASD = SpatRefLASD.SpatialReference.exportToString()
                arcpy.CheckOutExtension("3D")
                arcpy.PointFileInformation_3d(las_path, ptFileInfoFC, "LAS", "las", "", "NO_RECURSION", "NO_EXTRUSION",
                                              "DECIMAL_POINT", "NO_SUMMARIZE", "LAS_SPACING")
                arcpy.CheckInExtension("3D")
                messages = arcpy.GetMessages()
                arcpy.AddMessage("\nResults output from PointFileInformation_3d are: \n{0}\n".format(messages))
                del descLASD, SpatRefLASD
            else:
                arcpy.AddMessage("ptFileInfoFC already exists:  {0}".format(ptFileInfoFC))

            rows = arcpy.SearchCursor(ptFileInfoFC,
                                      fields="FileName; Pt_Spacing; Z_Min; Z_Max",
                                      sort_fields="FileName; Pt_Spacing; Z_Min; Z_Max")
            # Iterate through the rows in the cursor and store the
            # "FileName; Pt_Spacing; Z_Min; Z_Max"
            ptFileInfoList = []
            PtSpacing = []
            # Z Min & Z Max added for auto-detecting LiDAR tiles with potential artifacts in LiDAR in the future.
            for row in rows:
                formattedfields = ("{0}, {1}, {2}, {3}".format(
                    row.getValue("FileName"),
                    row.getValue("Pt_Spacing"),
                    row.getValue("Z_Min"),
                    row.getValue("Z_Max")))
                ptFileInfoList.append(formattedfields)
                ptspacinglist = float("{0}".format(row.getValue("Pt_Spacing")))
                PtSpacing.append(ptspacinglist)
                del row
            arcpy.AddMessage("ptFileInfoList:    {0}".format(str(ptFileInfoList)))
            arcpy.AddMessage("ptSpacing:    {0}".format(str(PtSpacing)))
            avgPtSpacing = sum(PtSpacing)/float(len(PtSpacing))
            arcpy.AddMessage("returning avgPtSpacing of:    {0}".format(str(avgPtSpacing)))
        else:
            arcpy.AddMessage("\nExiting get_las_point_spacing, since no LASD found: \n{0}\n".format(lasd))
            return ""
        del rows
        return avgPtSpacing

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))


def get_las_boundary(las_path, lasd, md_name, las_raster_type, file_gdb, surface_constraint_fc):
    try:
        # Ensure the LAS Raster type file exists
        if not os.path.exists(las_raster_type):
            arcpy.AddError("\nExiting: LAS Raster type file doesn't exist {0}".format(las_raster_type))
            return
        # Ensure the file gdb exists
        if not os.path.exists(file_gdb):
            arcpy.AddError("\nExiting: Geodatabase (in which to place boundary) doesn't exist {0}".format(file_gdb))
            return
        # Ensure the lasd exists
        if not arcpy.Exists(lasd):
            arcpy.AddError("\nExiting: LAS Dataset (from which to get boundary) doesn't exist {0}".format(lasd))
            return
        desc_lasd = arcpy.Describe(lasd)
        spat_ref_lasd = desc_lasd.SpatialReference
        spat_ref_lasd_str = desc_lasd.SpatialReference.exportToString()
        if spat_ref_lasd.PCSCode == 0:
            arcpy.AddWarning("\n*** NOTE: One or more LAS files has a PCSCode of 0.                  ***")
            arcpy.AddWarning("\n*** PCSCode = 0 indicates a non-standard datum or unit of measure.   ***")

        arcpy.AddMessage("\nSpatial reference of LASD is: \n{0}\n".format(spat_ref_lasd_str))
        #arcpy.AddMessage("Length of SR string is {0}:".format(len(SpatRefStringFirstLAS)))
        arcpy.AddMessage("Spatial Reference name of LAS Dataset:  {0}".format(spat_ref_lasd.name))
        arcpy.AddMessage("Spatial Reference XY Units of LAS Dataset: {0}".format(spat_ref_lasd.linearUnitName))

        loc = md_name.rfind(".gdb")
        # #arcpy.AddMessage("loc = {0}".format(loc))
        MD_ShortName = md_name[loc+5:]
        arcpy.AddMessage("Temp MD Short Name:  {0}".format(MD_ShortName))

        # Create a MD in same SR as LAS Dataset
        arcpy.CreateMosaicDataset_management(file_gdb, MD_ShortName,
                                             coordinate_system=spat_ref_lasd, num_bands="1", pixel_type="32_BIT_FLOAT",
                                             product_definition="NONE", product_band_definitions="#")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from CreateMosaicDataset are: \n{0}\n".format(messages))

        # Add the LAS files to the Mosaic Dataset and don't update the boundary yet.
        # The cell size of the Mosaic Dataset is determined by the art.xml file chosen by the user.
        arcpy.AddRastersToMosaicDataset_management(md_name, las_raster_type, las_path,
                                                   update_cellsize_ranges="UPDATE_CELL_SIZES",
                                                   update_boundary="NO_BOUNDARY", update_overviews="NO_OVERVIEWS",
                                                   maximum_pyramid_levels="#", maximum_cell_size="0",
                                                   minimum_dimension="1500", spatial_reference=spat_ref_lasd_str,
                                                   filter="*.las", sub_folder="NO_SUBFOLDERS",
                                                   duplicate_items_action="ALLOW_DUPLICATES",
                                                   build_pyramids="NO_PYRAMIDS", calculate_statistics="NO_STATISTICS",
                                                   build_thumbnails="NO_THUMBNAILS", operation_description="#",
                                                   force_spatial_reference="NO_FORCE_SPATIAL_REFERENCE")

        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from AddRastersToMosaicDataset are: \n{0}\n".format(messages))

        # Get a count of the number of LAS ingested
        result = arcpy.GetCount_management(md_name)
        countRowsWithLAS = int(result.getOutput(0))
        if countRowsWithLAS == 0:
            arcpy.AddMessage("\nNo LAS rows were ingested into {0}".format(md_name))
            return
        else:
            arcpy.AddMessage("{0} has {1} LAS row(s).".format(md_name, countRowsWithLAS))

        # Build Footprints with these non-standard parameters:
        #    min_region_size="20"
        #    approx_num_vertices="2000"
        #    Update the Boundary
        #    min_data_value="-1000"
        arcpy.BuildFootprints_management(md_name,where_clause="#", reset_footprint="RADIOMETRY", min_data_value="-1000",
                                         max_data_value="4294967295", approx_num_vertices="2000", shrink_distance="0",
                                         maintain_edges="MAINTAIN_EDGES", skip_derived_images="SKIP_DERIVED_IMAGES",
                                         update_boundary="UPDATE_BOUNDARY", request_size="2000", min_region_size="20",
                                         simplification_method="NONE", edge_tolerance="#", max_sliver_size="20",
                                         min_thinness_ratio="0.05")

        messages = arcpy.GetMessages()
        arcpy.AddMessage("Results output from BuildFootprints are: \n{0}\n".format(messages))

        # The boundary will potentially have lots of vertices, so simplify the feature after exporting.
        boundary_detailed = surface_constraint_fc + r"_detail"
        arcpy.ExportMosaicDatasetGeometry_management(md_name, boundary_detailed, where_clause="#",
                                                     geometry_type="BOUNDARY")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("Results output from ExportMosaicDatasetGeometry are: \n{0}\n".format(messages))

        arcpy.SimplifyPolygon_cartography(boundary_detailed, surface_constraint_fc,
                                          algorithm="POINT_REMOVE", tolerance="5 Meters", minimum_area="0 SquareMeters",
                                          error_option = "NO_CHECK", collapsed_point_option = "NO_KEEP")
                                          # error_option="RESOLVE_ERRORS", collapsed_point_option="KEEP_COLLAPSED_POINTS")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("Results output from SimplifyPolygon are: \n{0}\n".format(messages))
        del desc_lasd, spat_ref_lasd
        return

    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))


def main(first_parameter, second_parameter, z_unit, featureextraction, out_folder_path, processing_unit_length,
         use_pos_terrain_method, delete_intermediate_files, regularize_buildings):
    try:
        start = time.time()
        # resourceLogger.log()
        executed_from = sys.executable.upper()
        # Check out Spatial Analyst license
        arcpy.CheckOutExtension("Spatial")

        # This python code can be invoked by two different gp script tools, depending upon
        # whether the input is LAS or raster. Therefore the first two parameters can either be
        # LAS folder and DSM Creation Method or DSM (raster) Folder and DTM (raster) Folder.
        #
        # If the second parameter is DSM Creation Method ("ALL Returns" or "Last Returns"), then
        # the first parameter is assumed to be the LAS folder. Otherwise the first and second
        # parameters are assumed to be DSM folder and DTM folder.

        if second_parameter == "ALL RETURNS" or second_parameter == "LAST RETURNS":
            # Input is LAS
            las_path = first_parameter
            if not os.path.exists(las_path):
                arcpy.AddMessage("*** Exiting...the LAS Path does not exist: {0} ***".format(las_path))
                sys.exit()
        else:
            # Input is DSM & DTM Rasters
            dsm_path = first_parameter
            dtm_path = second_parameter
            if not os.path.exists(dsm_path):
                arcpy.AddMessage("*** Exiting...the DSM Path does not exist: {0} ***".format(dsm_path))
                sys.exit()
            if not os.path.exists(dtm_path):
                arcpy.AddMessage("*** Exiting...the DTM Path does not exist: {0} ***".format(dtm_path))
                sys.exit()

        if z_unit == "METER":
            elevation_meter_scalefactor = "1.0"
        elif z_unit == "FOOT":
            elevation_meter_scalefactor = "3.2808"  # Feet per Meter
        else:
            elevation_meter_scalefactor = "1.0"
        elevation_meter_scalefactor_str = str(elevation_meter_scalefactor)

        if os.path.exists(out_folder_path):
            arcpy.AddMessage("Results Folder exists:               {0}".format(out_folder_path))
        else:
            arcpy.AddMessage("Creating Results Folder:             {0}".format(out_folder_path))
            os.makedirs(out_folder_path)

        height_path = os.path.join(out_folder_path, r"HeightRasters")
        # If the path doesn't exist, then create the folder
        if not os.path.isdir(height_path):
            arcpy.AddMessage("Creating folder for Height Rasters:  {0}".format(height_path))
            os.makedirs(height_path)

        # Create a results file gdb to store trees and buildings results
        results_gdb_name = r"Results.gdb"
        # Entire path of the file gdb
        results_file_gdb_path = os.path.join(out_folder_path, results_gdb_name)
        # If the file gdb doesn't exist, then create it
        if not os.path.exists(results_file_gdb_path):
            arcpy.AddMessage("Creating Results File GDB:           {0}".format(results_file_gdb_path))
            arcpy.CreateFileGDB_management(out_folder_path, results_gdb_name, out_version="CURRENT")
        else:
            arcpy.AddMessage("\nResults File GDB already exists:  {0}".format(results_file_gdb_path))

        # feature classes to be created
        fishnet = os.path.join(results_file_gdb_path, r"aFishnet")
        tmpfishnet = os.path.join(results_file_gdb_path, r"tmpFishnet")
        all_trees_final = os.path.join(results_file_gdb_path, r"all_trees_final")
        all_buildings_final = os.path.join(results_file_gdb_path, r"all_buildings_final")

        # Exit process if all_buildings_final (and all_trees_final, if requested) already exist
        if arcpy.Exists(all_buildings_final):
            if "TREES" in featureextraction.upper():
                if arcpy.Exists(all_trees_final):
                    arcpy.AddMessage("\nExiting process...Buildings and Trees output products already "
                                     "exist: \n {0} \n {1}".format(all_buildings_final, all_trees_final))
                    sys.exit()
            else:
                arcpy.AddMessage("\nExiting process...Buildings output product already "
                                 "exists: \n {0}".format(all_buildings_final))
                sys.exit()

        # the following two fc's are only created if the user checks on
        #   "Regularize north-south oriented buildings"
        buildings_to_reg = os.path.join(results_file_gdb_path, r"partial_buildings_to_regularize")
        buildings_reg = os.path.join(results_file_gdb_path, r"partial_buildings_regularized")

        # the scratch gdb name
        # scratch_gdb_name = r"Scratch.gdb"
        scratch_gdb_name = r"TempWorkArea.gdb"
        # Create a gdb to store the intermediate products
        mp_gdb_name = r"MiscProducts.gdb"
        # Entire path of the Miscellaneous Intermediate Products file gdb
        mp_file_gdb_path = os.path.join(out_folder_path, mp_gdb_name)
        # If the Miscellanous Intermediate Products file gdb doesn't exist, then create it
        if not os.path.exists(mp_file_gdb_path):
            arcpy.AddMessage("Creating Miscellanous Intermediate Products File GDB:    {0}".format(mp_file_gdb_path))
            arcpy.CreateFileGDB_management(out_folder_path, mp_gdb_name, out_version="CURRENT")
        else:
            arcpy.AddMessage("\nMiscellanous Intermediate Products File GDB "
                             "already exists:  {0}".format(mp_file_gdb_path))

        # more Mosaic Datasets and Feature classes to be created...
        dsm_md_name = r"DSM"
        dsm_md = os.path.join(mp_file_gdb_path, dsm_md_name)
        dsm_boundary = dsm_md + r"_boundary"
        dtm_md_name = r"DTM"
        dtm_md = os.path.join(mp_file_gdb_path, dtm_md_name)
        dtm_boundary = dtm_md + r"_boundary"
        las_md_name = r"LAS"
        las_md = os.path.join(mp_file_gdb_path, las_md_name)
        # las_boundary = las_md + r"_boundary"
        height_md_name = r"HeightAboveGround"
        height_md = os.path.join(mp_file_gdb_path, height_md_name)
        buildings_merged = os.path.join(results_file_gdb_path, r"buildings_merged")
        las_point_spacing = 0.0
        lasd_name = r"LasDataset.lasd"
        lasd = os.path.join(out_folder_path, lasd_name)

        processing_area = os.path.join(mp_file_gdb_path, r"ProcessingArea")
        arcpy.AddMessage("General processing area FC:          {0}".format(processing_area))
        dsmClassCodes = ""
        dtmClassCodes = ""
        if second_parameter == "ALL RETURNS" or second_parameter == "LAST RETURNS":
            # Input is LAS
            if not arcpy.Exists(lasd):
                arcpy.AddMessage("Creating LAS Dataset: {}".format(lasd))
                create_lasd(las_path, lasd)
                if not arcpy.Exists(lasd):
                    arcpy.AddError("\nExiting...LAS Dataset not created: {0}".format(lasd))
                    sys.exit()
            lasdatasetstatstext = os.path.join(out_folder_path, "lasDatasetStatsText.txt")
            las_point_spacing = get_las_point_spacing(las_path, lasd, out_folder_path, mp_file_gdb_path,
                                                      lasdatasetstatstext)
            arcpy.AddMessage("\nLAS point spacing: {0}".format(str(las_point_spacing)))
            # # Set point spacing to 0.5 if get_las_point_spacing is unable to determine point spacing
            # if las_point_spacing_str == "":
            #     las_point_spacing_str = "0.5"

            # Get a list of both dsm and dtm class codes found in the lasdatasetstatstext file
            # Note: this is necessary to circumvent a bug in Make LAS Dataset Layer when run from ArcGIS Pro
            #       If invoked from Pro, Make LAS Dataset layer will bomb if passed a class code that does not
            #       exist in the input LAS Dataset.
            classCodesList = []
            if arcpy.Exists(lasdatasetstatstext):
                with open(lasdatasetstatstext) as csvfile:
                    reader = csv.DictReader(csvfile, delimiter=',', fieldnames=['File_Name','Item','Category'])
                    for row in reader:
                        # arcpy.AddMessage("row[Category] : {0}".format(row['Category']))
                        if row['Category'] == "ClassCodes":
                            classCodesList.append(int(row["Item"].split("_")[0]))
            # arcpy.AddMessage("Class Code List length: {0}".format(len(classCodesList)))
            # arcpy.AddMessage("Class Code List: {0}".format(classCodesList))
            uniqueClassCodesList = []
            for i in classCodesList:
                if i not in uniqueClassCodesList:
                    uniqueClassCodesList.append(i)
            arcpy.AddMessage("Unique Class Code List: {0}".format(uniqueClassCodesList))
            for i in [0, 1, 2, 3, 4, 5, 6]:
                if i in uniqueClassCodesList:
                    dsmClassCodes += str(i) + ";"
            if len(dsmClassCodes) > 0:
                dsmClassCodes = dsmClassCodes[:-1]
                arcpy.AddMessage("dsmClassCodes: {0}".format(dsmClassCodes))
            else:
                arcpy.AddMessage("Exiting...No DSM Class codes found.")
                sys.exit()
            for i in [2, 8]:
                if i in uniqueClassCodesList:
                    dtmClassCodes += str(i) + ";"
            if len(dtmClassCodes) > 0:
                dtmClassCodes = dtmClassCodes[:-1]
                arcpy.AddMessage("dtmClassCodes: {0}".format(dtmClassCodes))
            else:
                arcpy.AddMessage("Exiting...No DTM Class codes found.")
                sys.exit()
            # Use a template raster type to create a las_raster_type_file with a custom pixel size
            las_raster_type_template = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())),
                                                    r"LAS_Template_Raster_Type.art.xml")
            las_raster_type_file = os.path.join(out_folder_path, r"LAS_Raster_Type.art.xml")

            # Create a copy of the Template LAS Raster type for this execution, since the art.xml file
            # needs to be edited (in set_pixel_size) to modify the desired pixel size of the mosaic dataset.
            if not arcpy.Exists(las_raster_type_file):
                shutil.copyfile(las_raster_type_template, las_raster_type_file)
                desc_lasd = arcpy.Describe(lasd)
                spat_ref_lasd = desc_lasd.SpatialReference
                set_pixel_size(las_raster_type_file, las_point_spacing, spat_ref_lasd.linearUnitName)
                del desc_lasd, spat_ref_lasd
            else:
                arcpy.AddMessage("LAS Raster type file already exists: {0}".format(las_raster_type_file))
            # Load the LAS into a Mosaic dataset and return the boundary in processing_area
            if not arcpy.Exists(processing_area):
                get_las_boundary(las_path, lasd, las_md, las_raster_type_file, mp_file_gdb_path, processing_area)
            else:
                arcpy.AddMessage("LAS Boundary file already exists: {0}".format(processing_area))
            if not arcpy.Exists(processing_area):
                arcpy.AddError("\nExiting...Surface constraint FC not created: {0}".format(processing_area))
                sys.exit()

            # assign the boundary as a hard clip constraint
            surface_constraints = "'" + processing_area + "'" + r" <None> Hard_Clip"
            arcpy.AddFilesToLasDataset_management(lasd, "", "NO_RECURSION", surface_constraints)
        else:
            # Input is DSM & DTM Rasters
            if not arcpy.Exists(dsm_md):
                create_md_from_raster(dsm_path, mp_file_gdb_path, dsm_md_name, dsm_boundary, "true", "true")
                arcpy.AddMessage("Created MD:  {0}\n".format(dsm_md))
            else:
                arcpy.AddMessage("\nMD already exists:  {0}".format(dsm_md))
            if not arcpy.Exists(dtm_md):
                create_md_from_raster(dtm_path, mp_file_gdb_path, dtm_md_name, dtm_boundary, "true", "true")
                arcpy.AddMessage("Created MD:  {0}\n".format(dtm_md))
            else:
                arcpy.AddMessage("\nMD already exists:  {0}".format(dtm_md))

            # Find the intersection of the DSM and DTM datasets to determine the general processing area
            # Subsequent processing will further eliminate areas that don't need to be processed
            intersect_input_list = [dsm_boundary, dtm_boundary]
            arcpy.AddMessage("\nIntersection of:     {0}".format(intersect_input_list))
            if not arcpy.Exists(processing_area):
                arcpy.Intersect_analysis(intersect_input_list, processing_area,
                                         join_attributes="ONLY_FID", cluster_tolerance="", output_type="INPUT")
                messages = arcpy.GetMessages()
                arcpy.AddMessage("\nResults output from Intersect_analysis are: \n{0}\n".format(messages))
            else:
                arcpy.AddMessage("\nDSM & DTM intersection fc already exists:  {0}".format(processing_area))
            # desc_bound = arcpy.Describe(processing_area)

        # Get the bounds of the processing area to create a fishnet
        # If input is DSM and DTM folders, then processing area is the intersection of the two MD's (DSM & DTM)
        # If input is a LAS folder, then processing area is the boundary of the MD that contains the LAS
        desc_bound = arcpy.Describe(processing_area)
        xmin = desc_bound.Extent.XMin
        ymin = desc_bound.Extent.YMin
        # xmax = desc_bound.Extent.XMax
        ymax = desc_bound.Extent.YMax
        origin_coord = str(xmin) + " " + str(ymin)
        y_axis_coord = str(xmin) + " " + str(ymax)

        if not arcpy.Exists(fishnet):
            arcpy.CreateFishnet_management(tmpfishnet, origin_coord, y_axis_coord, processing_unit_length,
                                           processing_unit_length, number_rows="", number_columns="", corner_coord="",
                                           labels="NO_LABELS",
                                           template=processing_area, geometry_type="POLYGON")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from CreateFishnet are: \n{0}\n".format(messages))

            # Now keep only fishnet features that intersect the processing area
            fishnet_lyr = r"fishnetLyr"
            arcpy.MakeFeatureLayer_management(tmpfishnet, fishnet_lyr, "", "", "")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from MakeFeatureLayer are: \n{0}\n".format(messages))
            arcpy.SelectLayerByLocation_management(fishnet_lyr, "INTERSECT", processing_area,
                                                   "", "NEW_SELECTION", "NOT_INVERT")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from SelectLayerByLocation are: \n{0}\n".format(messages))
            arcpy.CopyFeatures_management(fishnet_lyr, fishnet, config_keyword="", spatial_grid_1="0",
                                          spatial_grid_2="0", spatial_grid_3="0")
            messages = arcpy.GetMessages()
            arcpy.AddMessage("\nResults output from CopyFeatures are: \n{0}\n".format(messages))
            arcpy.Delete_management(fishnet_lyr)
            arcpy.Delete_management(tmpfishnet)
        else:
            arcpy.AddMessage("\nFishnet fc already exists:  {0}".format(fishnet))

        # Intentionally set overwriteOutput here (as opposed to earlier in the script) because this might be a re-start
        #  of the script, in which case we don't want to re-do the pre-processing steps.  Subsequent logic
        #  manually checks for the existence of each folder (one per fishnet oid) before overwriting folder contents,
        #  so if the tool is being re-started the newest folder(s) may need to be deleted before restarting, especially
        # if the tool was forced to stop.
        arcpy.env.overwriteOutput = True

        # Initialize
        buildings_merged_list = []
        trees_merged_list = []
        # inititalize Lists used for passing arguments to extract_buildings_trees
        out_folder_path_list = []
        oid_list = []
        dsm_md_list = []    # only used if DSM & DTM Raster input
        dtm_md_list = []    # only used if DSM & DTM Raster input
        lasd_list = []      # only used if LAS input (i.e. if second_parameter == "ALL RETURNS" or "LAST RETURNS")
        dsm_type_list = []  # only used if LAS input (i.e. if second_parameter == "ALL RETURNS" or "LAST RETURNS")
        xmin_list = []
        xmax_list = []
        ymin_list = []
        ymax_list = []
        featureextraction_list = []
        elevation_meter_scalefactor_str_list = []
        use_pos_terrain_method_list = []
        delete_intermediate_files_list = []
        height_path_list = []
        point_spacing_list = []   # only used if LAS input (second_parameter == "ALL RETURNS" or "LAST RETURNS")
        dsm_ClassCodes_list = []
        dtm_ClassCodes_list = []

        fields = ["OID@", "SHAPE@"]
        with arcpy.da.SearchCursor(fishnet, fields) as sc:
            # iterate through the fishnet features to populate the arguments for extract_buildings_trees
            for row in sc:
                oid = str(row[0])
                geom = row[1]
                xmin = str(geom.extent.XMin)
                ymin = str(geom.extent.YMin)
                xmax = str(geom.extent.XMax)
                ymax = str(geom.extent.YMax)
                # populate lists for each parameter
                out_folder_path_list.append(out_folder_path)
                oid_list.append(oid)
                if second_parameter == "ALL RETURNS" or second_parameter == "LAST RETURNS":
                    lasd_list.append(lasd)
                    dsm_type_list.append(second_parameter)
                else:
                    dsm_md_list.append(dsm_md)
                    dtm_md_list.append(dtm_md)
                xmin_list.append(xmin)
                xmax_list.append(xmax)
                ymin_list.append(ymin)
                ymax_list.append(ymax)
                featureextraction_list.append(featureextraction)
                elevation_meter_scalefactor_str_list.append(elevation_meter_scalefactor_str)
                use_pos_terrain_method_list.append(use_pos_terrain_method)
                delete_intermediate_files_list.append(delete_intermediate_files)
                height_path_list.append(height_path)
                point_spacing_list.append(str(las_point_spacing))  # only applicable if LAS input
                dsm_ClassCodes_list.append(dsmClassCodes)          # only applicable if LAS input
                dtm_ClassCodes_list.append(dtmClassCodes)          # only applicable if LAS input
                del row, geom
            del sc

        num_iterations_str = str(len(out_folder_path_list))
        arcpy.AddMessage("\n** Number of iterations (fishnet features) is: {0} **\n".format(num_iterations_str))

        if second_parameter == "ALL RETURNS" or second_parameter == "LAST RETURNS":
            pp_params = [[out_folder_path_list[i], oid_list[i], lasd_list[i], dsm_type_list[i], xmin_list[i],
                          xmax_list[i], ymin_list[i], ymax_list[i], featureextraction_list[i],
                          elevation_meter_scalefactor_str_list[i], use_pos_terrain_method_list[i],
                          delete_intermediate_files_list[i], height_path_list[i], point_spacing_list[i],
                          dsm_ClassCodes_list[i], dtm_ClassCodes_list[i]]
                    for i in range(len(out_folder_path_list))]
        else:
            pp_params = [[out_folder_path_list[i], oid_list[i], dsm_md_list[i], dtm_md_list[i], xmin_list[i],
                          xmax_list[i], ymin_list[i], ymax_list[i], featureextraction_list[i],
                          elevation_meter_scalefactor_str_list[i], use_pos_terrain_method_list[i],
                          delete_intermediate_files_list[i], height_path_list[i], point_spacing_list[i],
                          dsm_ClassCodes_list[i], dtm_ClassCodes_list[i]]
                    for i in range(len(out_folder_path_list))]
        # arcpy.AddMessage("\n pp_params: {0}".format(pp_params))

        # If executing from the gp User Interface, then extract_buildings_trees will be run serially.
        # If executing from the command line, then extract_buildings_trees will be run in parallel.
        arcpy.AddMessage(executed_from)
        if "ARCMAP" in executed_from or "ARCCATALOG" in executed_from or "RUNTIME" in executed_from or \
                        "ARCGISPRO" in executed_from:
            list(map(extract_buildings_trees, pp_params))
        elif "PYTHON" in executed_from:
            # Number of cores to use (max will be 3 for now, otherwise we're I/O bound)
            cpu_num = min(multiprocessing.cpu_count(), 3)
            # Create the pool object
            pool = multiprocessing.Pool(processes=cpu_num, maxtasksperchild=1)
            arcpy.AddMessage("\nCPUs utilized: {0}".format(cpu_num))
            # Start Multiprocessing
            arcpy.AddMessage("Start Multiprocessing")
            pool.map(extract_buildings_trees, pp_params, chunksize=1)

            # Close the pool
            pool.close()
            pool.join()

        # clear extent for the remainder of processing - important step (or mosaic dataset functionality doesn't work)
        arcpy.env.extent = None

        # Create a mosaic dataset of all of the height rasters and use it later to get
        #   zonal stats on buildings and/or trees
        if not arcpy.Exists(height_md):
            # Don't need to build footprints, since this takes a while and isn't necessary
            create_md_from_raster(height_path, mp_file_gdb_path, height_md, "", "false", "false")
            arcpy.AddMessage("Created MD:  {0}".format(height_md))
        else:
            arcpy.AddMessage("\nMD already exists:  {0}".format(height_md))

        # If the user wanted trees output, then merge all of the tree feature classes
        if "TREES" in featureextraction.upper() and not arcpy.Exists(all_trees_final):
            for i in range(1, len(out_folder_path_list) + 1):
                sub_path = os.path.join(out_folder_path, str(i))
                # arcpy.AddMessage("sub Path: {0}".format(sub_path))
                if os.path.exists(sub_path):
                    sub_results_gdb = os.path.join(sub_path, results_gdb_name)
                    # arcpy.AddMessage("sub Results gdb: {0}".format(sub_results_gdb))
                    if arcpy.Exists(sub_results_gdb):
                        sub_fc = os.path.join(sub_results_gdb, r"trees_to_merge" + str(i))
                        # arcpy.AddMessage("sub fc: {0}".format(sub_fc))
                        if arcpy.Exists(sub_fc):
                            # Construct a semicolon delimited list of tree feature classes, for subsequent merging
                            # trees_merged_list = trees_merged_list + sub_fc + ";"
                            trees_merged_list.append(sub_fc)
            # arcpy.AddMessage("trees_merged_list: {0}".format(trees_merged_list))

            if len(trees_merged_list) > 0:
                # Merge all of the tree feature classes into one.
                arcpy.Merge_management(inputs=trees_merged_list, output=all_trees_final)
                messages = arcpy.GetMessages()
                arcpy.AddMessage("\nResults output from Merge are: \n{0}\n".format(messages))
                # Gather elevation statistics for trees
                zonalstats(all_trees_final, height_md, results_file_gdb_path)

        # merge all of the building feature classes
        if not arcpy.Exists(all_buildings_final):
            # Merge all of the building feature classes into one.
            for i in range(1, len(out_folder_path_list) + 1):
                sub_path = os.path.join(out_folder_path, str(i))
                # arcpy.AddMessage("sub Path: {0}".format(sub_path))
                if os.path.exists(sub_path):
                    sub_results_gdb = os.path.join(sub_path, results_gdb_name)
                    # arcpy.AddMessage("sub Results gdb: {0}".format(sub_results_gdb))
                    if arcpy.Exists(sub_results_gdb):
                        sub_fc = os.path.join(sub_results_gdb, r"buildings_to_merge" + str(i))
                        # arcpy.AddMessage("sub fc: {0}".format(sub_fc))
                        if arcpy.Exists(sub_fc):
                            # Construct a semicolon delimited list of building feature classes, for subsequent merging
                            # buildings_merged_list = buildings_merged_list + sub_fc + ";"
                            buildings_merged_list.append(sub_fc)
            # arcpy.AddMessage("buildings_merged_list: {0}".format(buildings_merged_list))

            if len(buildings_merged_list) > 0:
                arcpy.Merge_management(inputs=buildings_merged_list, output=buildings_merged)
                messages = arcpy.GetMessages()
                # arcpy.AddMessage("\nResults output from Merge are: \n{0}\n".format(messages))

                # Dissolve buildings into one feature class so that buildings at the borders of
                # each building feature class can be dissolved into one building feature.
                arcpy.Dissolve_management(buildings_merged, all_buildings_final, dissolve_field="",
                                          statistics_fields="",
                                          multi_part="SINGLE_PART", unsplit_lines="DISSOLVE_LINES")
                messages = arcpy.GetMessages()
                arcpy.AddMessage("\nResults output from Dissolve are: \n{0}\n".format(messages))
                # Gather elevation statistics for area under buildings
                zonalstats(all_buildings_final, height_md, results_file_gdb_path)

        # If user wants to regularize, then regularize those buildings that can best be regularized
        # (these are the buildings that are oriented North/South)
        # Instead of using this option, try using the Building Regularization gp tool in ArcGIS Pro or ArcMap 10.4
        if regularize_buildings == "true":
            # Don't bother if both feature classes already exist
            if not arcpy.Exists(buildings_to_reg) or not arcpy.Exists(buildings_reg):
                expression = "\"MBG_Width\"* \"MBG_Length\">1"
                # cellsize = Raster(diff_ori).meanCellHeight
                # Determine the cell size of the DSM Mosaic Dataset
                cellsize_result = arcpy.GetRasterProperties_management(dsm_md, property_type="CELLSIZEX", band_index="")
                cellsize = float(cellsize_result.getOutput(0))
                # arcpy.AddMessage("Cell size of MD:  {0}".format(cellsize))
                regularize_some_buildings(all_buildings_final, buildings_to_reg, buildings_reg, expression,
                                          cellsize, out_folder_path)
            else:
                arcpy.AddMessage("\npartial_buildings_regularized and partial_buildings_to_regularize already exist")

        # Delete all intermediate files if user checked on "Delete all intermediate files"
        if delete_intermediate_files == "true":
            if arcpy.Exists(buildings_merged):
                arcpy.Delete_management(buildings_merged)
            # Clean up the rasters in height_path and delete the height_path directory
            try:
                remove_rasters(height_path)
                arcpy.Delete_management(height_path)
            except:
                arcpy.AddMessage("Unable to clean up directory: {0}".format(height_path))

            # Delete files created during building regularization
            if regularize_buildings == "true":
                remove_rasters(out_folder_path)
                remove_shapefiles(out_folder_path)

            # Delete tables created during Zonal Statistics creation (in def zonalstats)
            remove_tables(results_file_gdb_path)

            # In each oid sub-folder, delete the individual Results.gdb & Scratch.gdb and all of their feature classes
            for i in range(1, len(out_folder_path_list) + 1):
                sub_path = os.path.join(out_folder_path, str(i))
                # arcpy.AddMessage("sub Path: {0}".format(sub_path))
                if os.path.exists(sub_path):
                    sub_results_gdb = os.path.join(sub_path, results_gdb_name)
                    if arcpy.Exists(sub_results_gdb):
                        try:
                            remove_filegdb(sub_results_gdb)
                        except:
                            arcpy.AddMessage("Unable to delete file GDB: {0}".format(sub_results_gdb))
                    sub_scratch_gdb = os.path.join(sub_path, scratch_gdb_name)
                    if arcpy.Exists(sub_scratch_gdb):
                        try:
                            remove_filegdb(sub_scratch_gdb)
                        except:
                            arcpy.AddMessage("Unable to delete file GDB: {0}".format(sub_scratch_gdb))

            # Delete Mosaic Datasets in MosaicDatasets.gdb, then delete remaining fc's and then MosaicDatasets.gdb
            if arcpy.Exists(dsm_md):
                arcpy.Delete_management(dsm_md)
            if arcpy.Exists(dtm_md):
                arcpy.Delete_management(dtm_md)
            if arcpy.Exists(height_md):
                arcpy.Delete_management(height_md)
            # if arcpy.Exists(las_md):
            #     arcpy.Delete_management(las_md)

            # remove_filegdb(mp_file_gdb_path)

        end = time.time()
        delta = end - start
        # This is useful if the tool is run at the command line
        arcpy.AddMessage("***** Total elapsed time is {0} hours *****".format(delta/3600))

    # except arcpy.ExecuteError:
    #     print(arcpy.GetMessages())
    except Exception as e:
        arcpy.AddMessage(e.message)
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        arcpy.AddError(
            "PYTHON ERRORS:\nTraceback Info:\n {0}\nError Info:\n {1}: {2}".format(tbinfo, sys.exc_type, sys.exc_value))
        messages = arcpy.GetMessages()
        arcpy.AddError("GP Messages:\n {0}".format(messages))
    finally:
        arcpy.CheckInExtension('3D')
        arcpy.CheckInExtension('Spatial')

if __name__ == '__main__':
    arcpy.AddMessage(inspect.getfile(inspect.currentframe()))
    arcpy.AddMessage(os.path.dirname(inspect.getfile(inspect.currentframe())))
    arcpy.AddMessage(sys.version)
    arcpy.AddMessage(sys.executable)
    executed_from = sys.executable.upper()

    PYTHON_EXE = os.path.join(sys.exec_prefix, 'pythonw.exe')
    # use pythonw for multiprocessing
    multiprocessing.set_executable(PYTHON_EXE)

    # This python code can be invoked by two different script tools, depending upon
    # whether the input is LAS or raster. Therefore the first two parameters can either be
    # LAS folder and DSM Creation Method or DSM Folder and DTM Folder.
    #
    # If the second parameter is either "ALL Returns" or "Last Returns" then
    # the first parameter is assumed to be the LAS folder.

    # if "ARCGISPRO" in executed_from:
    #     arcpy.AddMessage("Exiting...this tool does not yet run from ArcGIS Pro")
    #     # sys.exit(0)

    if not ("ARCMAP" in executed_from or "ARCCATALOG" in executed_from or
            "RUNTIME" in executed_from or "ARCGISPRO" in executed_from):
        arcpy.AddMessage("Getting parameters from command line...")

        second_parameter = sys.argv[2]
        second_parameter = second_parameter.strip()
        if second_parameter.upper() == "ALL RETURNS" or second_parameter.upper() == "LAST RETURNS":
            # Input is LAS
            first_parameter = sys.argv[1]
            first_parameter = first_parameter.strip()
            arcpy.AddMessage("LAS Path:                            {0}".format(first_parameter))
            second_parameter = second_parameter.upper()
            arcpy.AddMessage("DSM Creation Method:                 {0}".format(second_parameter))
        else:
            # Input is DSM & DTM Rasters
            first_parameter = sys.argv[1]
            first_parameter = first_parameter.strip()
            arcpy.AddMessage("DSM Path:                            {0}".format(first_parameter))
            arcpy.AddMessage("DTM Path:                            {0}".format(second_parameter))

        # dsm_path = sys.argv[1]
        # dsm_path = dsm_path.strip()
        # arcpy.AddMessage("DSM Path:                            {0}".format(dsm_path))
        # dtm_path = sys.argv[2]
        # dtm_path = dtm_path.strip()
        # arcpy.AddMessage("DTM Path:                            {0}".format(dtm_path))

        z_unit = sys.argv[3]
        z_unit = z_unit.upper()
        arcpy.AddMessage("Z Unit:                              {0}".format(z_unit))

        featureextraction = sys.argv[4]
        arcpy.AddMessage("Desired Features extracted:          {0}".format(featureextraction))

        out_folder_path = sys.argv[5]
        out_folder_path = out_folder_path.strip()
        arcpy.AddMessage("Output Folder Path:                  {0}".format(out_folder_path))

        processing_unit_length = sys.argv[6]
        arcpy.AddMessage("Processing Unit Distance:            {0}".format(processing_unit_length))

        use_pos_terrain_method = sys.argv[7]
        arcpy.AddMessage("Use alternative positive terrain method for buildings: {0}".format(use_pos_terrain_method))

        delete_intermediate_files = sys.argv[8]
        arcpy.AddMessage("Delete all intermediate files:       {0}".format(delete_intermediate_files))

        regularize_buildings = sys.argv[9]
        arcpy.AddMessage("Regularize some buildings:           {0}".format(regularize_buildings))

    else:
        arcpy.AddMessage("Getting parameters from GetParameterAsText...")

        second_parameter = arcpy.GetParameterAsText(1)
        second_parameter = second_parameter.strip()
        if second_parameter.upper() == "ALL RETURNS" or second_parameter.upper() == "LAST RETURNS":
            # Input is LAS
            first_parameter = arcpy.GetParameterAsText(0)
            first_parameter = first_parameter.strip()
            arcpy.AddMessage("LAS Path:                            {0}".format(first_parameter))
            second_parameter = second_parameter.upper()
            arcpy.AddMessage("DSM Creation Method:                 {0}".format(second_parameter))
        else:
            # Input is DSM & DTM Rasters
            first_parameter = arcpy.GetParameterAsText(0)
            first_parameter = first_parameter.strip()
            arcpy.AddMessage("DSM Path:                            {0}".format(first_parameter))
            arcpy.AddMessage("DTM Path:                            {0}".format(second_parameter))

        z_unit = arcpy.GetParameterAsText(2)
        z_unit = z_unit.upper()
        arcpy.AddMessage("Z Unit:                              {0}".format(z_unit))

        featureextraction = arcpy.GetParameterAsText(3)
        arcpy.AddMessage("Desired Features extracted:          {0}".format(featureextraction))

        out_folder_path = arcpy.GetParameterAsText(4)
        out_folder_path = out_folder_path.strip()
        arcpy.AddMessage("Output Folder Path:                  {0}".format(out_folder_path))

        processing_unit_length = arcpy.GetParameterAsText(5)
        arcpy.AddMessage("Processing Unit Distance:            {0}".format(processing_unit_length))

        use_pos_terrain_method = arcpy.GetParameterAsText(6)
        arcpy.AddMessage("Use alternative positive terrain method for buildings: {0}".format(use_pos_terrain_method))

        delete_intermediate_files = arcpy.GetParameterAsText(7)
        arcpy.AddMessage("Delete all intermediate files:       {0}".format(delete_intermediate_files))

        regularize_buildings = arcpy.GetParameterAsText(8)
        arcpy.AddMessage("Regularize some buildings:           {0}".format(regularize_buildings))

    main(first_parameter, second_parameter, z_unit, featureextraction, out_folder_path, processing_unit_length,
         use_pos_terrain_method, delete_intermediate_files, regularize_buildings)

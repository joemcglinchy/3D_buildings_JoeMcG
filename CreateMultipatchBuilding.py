#-------------------------------------------------------------------------------
# Name:        CreateMultipatchBuilding.py
# Purpose:     process for extracting LOD2 multipatch building from Classified &
#              Unclassifierd LiDAR.
#
#              Building Footprints and LiDAR in .las Folder, .las file or
#              .las dataset are required as inputs.
#              a spatial reference output folder must be designated as well as
#              the LiDAR building class codes.
#
# Author:      Geoff Taylor, Joe McGlinchy & Dennis Jarrard
#
# Created:     28/09/2015
# Copyright:   (c) Esri 2015
# Licence:     Internal Use Only!
#-------------------------------------------------------------------------------


import arcpy
import os
import os.path
import tempfile
import glob
from datetime import datetime

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('spatial')
arcpy.CheckOutExtension('3d')

inLASD = r''
inLAS = r'C:\Users\geof7015\PycharmProjects\testData\Boulder\LAS\329.las'
DTMRaster = r''  # C:\workspace\data\testdata\bh12TVK1800084000.img
DSMRaster = r''  # 'C:\workspace\data\testdata\hh12TVK1800084000.img'
buildingFootprints = r'C:\Users\geof7015\PycharmProjects\testData\Boulder\Data.gdb\Building329_small'
sr = "PROJCS['NAD_1983_HARN_StatePlane_Colorado_North_FIPS_0501_Feet',GEOGCS['GCS_North_American_1983_HARN',DATUM['D_North_American_1983_HARN',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',3000000.000316083],PARAMETER['False_Northing',999999.999996],PARAMETER['Central_Meridian',-105.5],PARAMETER['Standard_Parallel_1',39.71666666666667],PARAMETER['Standard_Parallel_2',40.78333333333333],PARAMETER['Latitude_Of_Origin',39.33333333333334],UNIT['Foot_US',0.3048006096012192]]"
outputWS = r'C:\Users\geof7015\PycharmProjects\testData\Boulder\Workspace.gdb'

#inLASD = r''
#inLAS = r'E:\3D_City_Data\United States\North Carolina\Charlotte\AIC\LiDAR'
#DTMRaster = r''  # C:\workspace\data\testdata\bh12TVK1800084000.img
#DSMRaster = r''  # 'C:\workspace\data\testdata\hh12TVK1800084000.img'
#buildingFootprints = r'E:\3D_City_Data\United States\North Carolina\Charlotte\AIC\Data.gdb\BuildingFootprints'
#sr = "PROJCS['NAD_1983_StatePlane_North_Carolina_FIPS_3200_Feet',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',2000000.002616666],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-79.0],PARAMETER['Standard_Parallel_1',34.33333333333334],PARAMETER['Standard_Parallel_2',36.16666666666666],PARAMETER['Latitude_Of_Origin',33.75],UNIT['Foot_US',0.3048006096012192]]"
#outputWS = r'E:\3D_City_Data\United States\North Carolina\Charlotte\AIC\Workspace1.gdb'

scratchGDB = arcpy.env.scratchGDB
tempFolder = tempfile.mkdtemp()

buildingClassCode = 6
groundClassCode = 2

groundReturn = ""
buildingReturn = "Last Return"

beginOnFeatureNumber = 0
pointSpacingCorrectionFactor = 0.75

# For Point-Cloud to raster process
interpolateBetweenPoints = True
rasterExtractionApproach = True

# For Raster Only Extraction Approach

#Currently Disconnected...
reduceTesselations = True

## TODO Check and ensure Con raster in_memory bug is resolved in arcpy python v 3.4 before enabling!
optimizeRaster = True
optimizeRasterFactor = 0.5

###############
# Definitions #
###############


def createlasdataset(inLAS, sr):
    global inLASD
    inLASD = os.path.join(tempFolder, "LASDataSet.lasd")
    if arcpy.Exists(inLASD):
        arcpy.Delete_management(inLASD)
    arcpy.CreateLasDataset_management(inLAS, inLASD, False, "", sr, "COMPUTE_STATS")
    if arcpy.Exists(inLASD):
        arcpy.AddMessage("LASD File Created @ Location: " + inLASD)
        return inLASD
        # for multiples:  return inLASD,output2,output3,etc...
    else:
        arcpy.AddMessage("Could Not Create LASD DataSet. Check LAS inputs for errors")


def findLasDatasetStatisticsfilePaths(file):
    file_object = open(file, 'r')
    lines = file_object.readlines()
    file_object.close()
    cleanLines = []
    for line in lines:
        if len(line) > 1:
            path = line.split(",")[0]
            if os.path.isabs(path) is True and path not in cleanLines:
                cleanLines.append(path)
    return cleanLines


# Create Lists with LiDAR Statistical Information.  Pt Spacing etc... Process only used in other modules.
def obtainLiDARInfo(inLASD,lasList):
    if arcpy.Exists(inLASD):
        arcpy.AddMessage("Calculating Necessary Statistics for Feature Extraction Process")
        lasDatasetStatsText = os.path.join(tempFolder, "lasDatasetStatsText.txt")
        if arcpy.Exists(lasDatasetStatsText):
            arcpy.Delete_management(lasDatasetStatsText)
        arcpy.LasDatasetStatistics_management(inLASD, "true", lasDatasetStatsText, "LAS_FILES", "COMMA",
                                              "DECIMAL_POINT")

        # TODO DJARRARD obtain a LiDAR file from list and parse the point_spacing to building footprints.
        # TODO DJARRARD if multiple LiDAR tiles overlap building footprints then point_spacing = pt_spacing_average
        #if recursivelyCreateAndClipRastersFromLasd:
            #pass

        # run arcpy.PointFileInfo_3d on the single tile (no recursion)
        ptFileInfoFile = os.path.join(outputWS, 'ptFileInfoFile')
        if arcpy.Exists(ptFileInfoFile):
            arcpy.Delete_management(ptFileInfoFile)
        arcpy.PointFileInformation_3d(lasList, ptFileInfoFile, "LAS", None, sr, "false", "false", "DECIMAL_POINT",
                                      "false", "false")

        rows = arcpy.SearchCursor(ptFileInfoFile,
                                  fields="FileName; Pt_Spacing; Z_Min; Z_Max",
                                  sort_fields="FileName; Pt_Spacing; Z_Min; Z_Max")
        # Iterate through the rows in the cursor and store the
        # "FileName; Pt_Spacing; Z_Min; Z_Max"
        ptFileInfoList = []
        PtSpacing = []
        for row in rows:
            formattedfields = ("{0}, {1}, {2}, {3}".format(
                row.getValue("FileName"),
                row.getValue("Pt_Spacing"),
                row.getValue("Z_Min"),
                row.getValue("Z_Max")))
            ptFileInfoList.append(formattedfields)
            ptspacinglist = float("{0}".format(row.getValue("Pt_Spacing")))
            PtSpacing.append(ptspacinglist)
        print(ptFileInfoList)
        print(PtSpacing)
        avgPtSpacing = sum(PtSpacing)/float(len(PtSpacing))
        print(avgPtSpacing)
        return ptFileInfoFile, ptFileInfoList, PtSpacing, avgPtSpacing


# Omitted in current process (for Surface Rasters Only)
# Raster process will be added for Non-LiDAR DSM version of tool when Con(Raster) GP tool in_memory bug is resolved.
def interpolateBetweenLasPts(LrDSM):
    # Run raster interpolation algorithm on LiDAR derived rasters if interpolateAdditionalPoints is True & <> Recursive
    TimesRaster = os.path.join(tempFolder, "TimesRaster.tif")
    if arcpy.Exists(TimesRaster):
        arcpy.Delete_management(TimesRaster)
    arcpy.Times_3d(LrDSM, 100, TimesRaster)
    arcpy.AddMessage("Times Raster Complete")

    IntegerRaster = os.path.join(tempFolder, "IntRaster.tif")
    if arcpy.Exists(IntegerRaster):
        arcpy.Delete_management(IntegerRaster)
    arcpy.Int_3d(TimesRaster, IntegerRaster)
    arcpy.AddMessage("Integer Raster Complete")

    BoundaryCleanRaster = os.path.join(tempFolder, "BoundaryClean.tif")
    if arcpy.Exists(BoundaryCleanRaster):
        arcpy.Delete_management(BoundaryCleanRaster)
    BC = arcpy.sa.BoundaryClean(IntegerRaster, "NO_SORT", "true")
    BC.save(BoundaryCleanRaster)
    arcpy.AddMessage("BoundaryClean Raster Complete")

    FloatRaster = os.path.join(tempFolder, "FloatRaster.tif")
    if arcpy.Exists(FloatRaster):
        arcpy.Delete_management(FloatRaster)
    arcpy.Float_3d(BoundaryCleanRaster, FloatRaster)
    arcpy.AddMessage("Float Raster Complete")

    if arcpy.Exists(LrDSM):
        arcpy.Delete_management(LrDSM)
    arcpy.Divide_3d(FloatRaster, 100, LrDSM)
    arcpy.AddMessage("Divide Raster Complete")
    return LrDSM


# Omitted in current process (for Surface Rasters Only)
# Raster process will be added for Non-LiDAR DSM version of tool when Con(Raster) GP tool in_memory bug is resolved.
def slopedAreaRasters(SlopeRaster, slopedAreasNullRaster):
    # TODO Fix Memory Leak 1
    slopedAreasRaster = os.path.join(tempFolder, "slopedAreasRaster.tif")
    if arcpy.Exists(slopedAreasRaster):
        arcpy.Delete_management(slopedAreasRaster)
    slopedAreasRasterProcess = arcpy.sa.Con(SlopeRaster, 1, 0, "VALUE >= 20")
    slopedAreasRasterProcess.save(slopedAreasRaster)
    # TODO Fix Memory Leak 2
    if arcpy.Exists(slopedAreasNullRaster):
        arcpy.Delete_management(slopedAreasNullRaster)
    slopedAreasNullRasterProcess = arcpy.sa.SetNull(slopedAreasRaster, 1, "Value = 0")
    slopedAreasNullRasterProcess.save(slopedAreasNullRaster)

    arcpy.Delete_management(slopedAreasRaster)

    return slopedAreasNullRaster


# Omitted in current process (for Surface Rasters Only)
# Raster process will be added for Non-LiDAR DSM version of tool when Con(Raster) GP tool in_memory bug is resolved.
def reduceTesselationProcess(LrDSM, SlopedAreasPolygonBuffered):
    SlopeRaster = os.path.join(tempFolder, "SlopeRaster.tif")
    if arcpy.Exists(SlopeRaster):
        arcpy.Delete_management(SlopeRaster)
    arcpy.Slope_3d(LrDSM, SlopeRaster, "DEGREE", 1)

    slopedAreasNullRaster = os.path.join(tempFolder, "slopedAreasNullRaster.tif")
    slopedAreaRasters(SlopeRaster=SlopeRaster, slopedAreasNullRaster=slopedAreasNullRaster)

    SlopedAreasPolygon = os.path.join(tempFolder, "SlopedAreasPolygon.shp")
    if arcpy.Exists(SlopedAreasPolygon):
        arcpy.Delete_management(SlopedAreasPolygon)
    arcpy.RasterToPolygon_conversion(slopedAreasNullRaster, SlopedAreasPolygon, "false", "Value")

    if arcpy.Exists(SlopedAreasPolygonBuffered):
        arcpy.Delete_management(SlopedAreasPolygonBuffered)
    arcpy.Buffer_analysis(SlopedAreasPolygon, SlopedAreasPolygonBuffered, "2 Feet", "FULL", "ROUND", "ALL", None, "PLANAR")

    arcpy.Delete_management(slopedAreasNullRaster)

    return SlopedAreasPolygonBuffered


# Automatically removes Artifacts from point-clouds.
# Resolves issues where building sides were triangulated and other geometric flaws
def cleanupArtifacts(single_bldg_pts, single_bldg_pts_cleaned):

    arcpy.Near3D_3d(single_bldg_pts, single_bldg_pts, str(1.4 * pointSpace), "LOCATION", "ANGLE", "DELTA")

    bldgpts = os.path.join("in_memory", "bldgPoints")
    arcpy.MakeFeatureLayer_management(single_bldg_pts, bldgpts)

    arcpy.SelectLayerByAttribute_management(bldgpts, "NEW_SELECTION", "NEAR_DELTZ >= -1 Or NEAR_ANG_V <> 0")

    arcpy.CopyFeatures_management(bldgpts, single_bldg_pts_cleaned)
    print("Artifacts Removed")
    if arcpy.Exists(bldgpts):
        arcpy.Delete_management(bldgpts)
    return single_bldg_pts_cleaned


# assigns nearest point value using planar distance to building footprint points from an input point feature class.
def interpolatePointsToBoundary(input_bldg_points, input_bldg_fp, output_bldg_points_with_border):

    # explode input multipoint FC to single part
    # it is understood the input point FC will be multipoint. Need to convert to single part features
    #single_bldg_pts = os.path.join(outputWS, "singlepts")
    single_bldg_pts = os.path.join("in_memory", "singlepts")
    arcpy.MultipartToSinglepart_management(input_bldg_points, single_bldg_pts)
    print("created Single Pts")

    # Cleanup Artifacts
    # single_bldg_pts_cleaned = os.path.join(outputWS, "single_bldg_pts_cleaned")
    single_bldg_pts_cleaned = os.path.join("in_memory", "single_bldg_pts_cleaned")
    cleanupArtifacts(single_bldg_pts=single_bldg_pts, single_bldg_pts_cleaned=single_bldg_pts_cleaned)
    print("artifacts Cleaned")

    # add geometry attributes
    arcpy.AddGeometryAttributes_management(Input_Features=single_bldg_pts_cleaned, Geometry_Properties="POINT_X_Y_Z_M",
                                           Length_Unit="", Area_Unit="", Coordinate_System="")
    print("added attr to geometry")

    # process the building footprint
    footprintBuffer = os.path.join("in_memory", "footprintBuffer")
    arcpy.Buffer_analysis(input_bldg_fp, footprintBuffer, "0.5 Feet", "FULL", "FLAT", "NONE", None, "GEODESIC")

    # convert to line
    # bldg_line = os.path.join(outputWS, "bldgline")
    bldg_line = os.path.join("in_memory", "bldgline")
    arcpy.FeatureToLine_management(in_features=footprintBuffer, out_feature_class=bldg_line, cluster_tolerance=None,
                                   attributes="NO_ATTRIBUTES")
    if arcpy.Exists(footprintBuffer):
        arcpy.Delete_management(footprintBuffer)

    # Densify
    arcpy.Densify_edit(in_features=bldg_line, densification_method="DISTANCE", distance="1 Feet",
                       max_deviation="0.33 Feet", max_angle="10")

    # convert to points
    # bldg_ln_pts = os.path.join(outputWS, "bldglinepts")
    bldg_ln_pts = os.path.join("in_memory", "bldglinepts")
    arcpy.FeatureVerticesToPoints_management(in_features=bldg_line, out_feature_class=bldg_ln_pts, point_location="ALL")

    # use Near tool to identify point FID from building points to the boundary points
    arcpy.Near_analysis(in_features=bldg_ln_pts, near_features=single_bldg_pts_cleaned, search_radius="5 Feet",
                        location="NO_LOCATION", angle="NO_ANGLE", method="PLANAR")

    # now, grab the NEARI_FID field and assign that feature's z-value to the building footprint point z value
    arcpy.AddField_management(bldg_ln_pts, "z_val", "DOUBLE")
    tbl_fp = arcpy.da.FeatureClassToNumPyArray(bldg_ln_pts, ["NEAR_FID"])
    tbl_pts = arcpy.da.FeatureClassToNumPyArray(single_bldg_pts_cleaned, ["POINT_Z"])

    # update the z_val attribute
    with arcpy.da.UpdateCursor(bldg_ln_pts, ["z_val"]) as Pointsc:
        for i, row in enumerate(Pointsc):
            fid = tbl_fp[i][0]
            row[0] = tbl_pts[fid-1][0]
            # print(row[0])
            Pointsc.updateRow(row)

    # convert to 3D and copy
    #bldg_ln_pts_z = os.path.join(outputWS, "bldg_ln_pts_Z")
    bldg_ln_pts_z = os.path.join("in_memory", "bldg_ln_pts_Z")
    arcpy.FeatureTo3DByAttribute_3d(bldg_ln_pts, bldg_ln_pts_z, "z_val")

    # pointsMerged = os.path.join("in_memory", "pointsMerged")
    arcpy.Merge_management([bldg_ln_pts_z, single_bldg_pts_cleaned], output_bldg_points_with_border)

    # Remove Intermediate Data
    if arcpy.Exists(single_bldg_pts):
        arcpy.Delete_management(single_bldg_pts)
    if arcpy.Exists(bldg_line):
        arcpy.Delete_management(bldg_line)
    if arcpy.Exists(bldg_ln_pts):
       arcpy.Delete_management(bldg_ln_pts)
    if arcpy.Exists(bldg_ln_pts_z):
        arcpy.Delete_management(bldg_ln_pts_z)

    return output_bldg_points_with_border


def extractMultipatchFromPts(fullextent, row):
    try:
        arcpy.env.extent = fullextent
        # get raster extent
        geom = row[0]
        print("geom = ", geom)

        # copy the feature temporarily
        #tp = os.path.join("in_memory", "tp{0}".format(i))
        tp = os.path.join(outputWS, "tp{0}".format(i))
        #tp = os.path.join("in_memory", "tp{0}".format(i))
        arcpy.CopyFeatures_management(geom, tp)

        extentgeom = arcpy.Describe(tp)
        arcpy.env.mask = tp
        print("extentgeom = ", extentgeom)
        extent = "{0} {1} {2} {3}".format(extentgeom.extent.XMin, extentgeom.extent.YMin, extentgeom.extent.XMax, extentgeom.extent.YMax)
        print("extent = ", extent)
        arcpy.env.extent = extent

        print("Begin Raster Creation Process")

        DTM = os.path.join("in_memory", "DTM")
        # Delete terrain rasters if existing.
        if arcpy.Exists(DTM):
            arcpy.Delete_management(DTM)

        arcpy.LasDatasetToRaster_conversion("DTMLASD", DTM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                            valueField, "CELLSIZE", pointSpacing, heightValue)
        print("Created DTM Raster at location: " + DTM)

        arcpy.env.mask = tp

        arcpy.env.extent = extent

         # clip the DTM
        print('clipping DTM')
        dtmClipRast = os.path.join("in_memory", 'tempDEMclip{0}'.format(i + beginOnFeatureNumber))
        arcpy.Clip_management(DTM, extent, dtmClipRast, tp, "true", "false")
        # convert DEM to Int
        #dtmClipRastInt = Int(dtmClipRast)

        # add Min Height to Building Footprints
        print('determining Minimum Building Elevation')
        arcpy.AddField_management(tp, "ID", "SHORT", None, None, None, "ID", "true", "true", None)
        arcpy.CalculateField_management(tp, "ID", 1, "PYTHON_9.3", None)
        minMaxElevTable = os.path.join("in_memory", "minMaxElevTable")
        arcpy.sa.ZonalStatisticsAsTable(tp, "ID", DTM, minMaxElevTable, "true", "MIN_MAX_MEAN")
        arcpy.JoinField_management(tp, "ID", minMaxElevTable, "ID", "MIN;MAX")

        if arcpy.Exists(DTM):
            arcpy.Delete_management(DTM)

        # then, move building footprints to MIN Z Height
        #out_poly3d = os.path.join("in_memory", "out_poly3d")
        out_poly3d = os.path.join(outputWS, "out_poly3d_{0}".format(i))
        arcpy.FeatureTo3DByAttribute_3d(tp, out_poly3d, "MIN", "")

        # multipoint = os.path.join("in_memory", "multipoint")
        #multipoint = os.path.join("in_memory", "multipoint")
        multipoint = os.path.join(outputWS, "multipoint{0}".format(i))
        if arcpy.Exists(multipoint):
            arcpy.Delete_management(multipoint)
        arcpy.LASToMultipoint_3d(lasList, multipoint, pointSpacing, buildingClassCode, buildingReturn, None, sr, "las",
                                 1, "false")
        print("Las to Multipoint complete")

        roofPoints = os.path.join(outputWS, "roofPoints{0}".format(i))
        if arcpy.Exists(roofPoints):
            arcpy.Delete_management(roofPoints)
        arcpy.Clip_analysis(multipoint, tp, roofPoints, None)

        # Delete Mulipoint shp
        #if arcpy.Exists(multipoint):
        #    arcpy.Delete_management(multipoint)

        # Interpolate Points to Boundary
        buildingInsideAndBorderPoints = os.path.join(outputWS, "buildingBorderPoints{0}".format(i))
        interpolatePointsToBoundary(input_bldg_points=roofPoints, input_bldg_fp=tp,
                                    output_bldg_points_with_border=buildingInsideAndBorderPoints)

        if arcpy.Exists(roofPoints):
            arcpy.Delete_management(roofPoints)
        roofPoints = os.path.join("in_memory", "roofPoints")

        arcpy.Dissolve_management(buildingInsideAndBorderPoints, roofPoints, None, None, "true", "false")

        #if arcpy.Exists(buildingInsideAndBorderPoints):
        #    arcpy.Delete_management(buildingInsideAndBorderPoints)

        # TODO: Resolve issue where roof-tin will not process
        # Check to ensure that paths can have spaces.  may be the problem.
        roofTin = os.path.join(tempFolder, "roofTin")
        arcpy.CreateTin_3d(roofTin, sr, "{0} Shape.Z Mass_Points <None>".format(roofPoints), "DELAUNAY")
        print("roof Tin Created")

        #if arcpy.Exists(roofPoints):
        #    arcpy.Delete_management(roofPoints)

        # make ground TIN
        gnd_feats_tin = "{} Shape.Z Hard_Clip <None>;".format(out_poly3d)
        out_gnd_tin = os.path.join(tempFolder, "gndTin")
        arcpy.CreateTin_3d(out_gnd_tin, sr, gnd_feats_tin, "DELAUNAY")

        # extrude polygon between TINs
        print('creating Multipatch')
        this_MP = os.path.join(outputWS, "bldgMP_{0}".format(i))
        arcpy.ExtrudeBetween_3d(roofTin, out_gnd_tin, out_poly3d, this_MP)

        # add feature name to list
        mp_list.append(this_MP)

        # Delete Unnecessary files
        #arcpy.Delete_management(tp)
        arcpy.Delete_management(minMaxElevTable)
        #arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(dtmClipRast)
        arcpy.Delete_management(out_gnd_tin)
        arcpy.Delete_management(roofTin)
        arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(buildingInsideAndBorderPoints)
        del row, tp

        print("Multipatch {0} Process complete @ ".format(i + beginOnFeatureNumber), str(datetime.now()))


        # TODO Geoff7015 Incorporate Cleanup Building CGA from Geof7015 rule into tool here:
        ''' every multipatch building must have LiDAR point spacing as a attribute and "Units: feet/meters
            will need to update CGA cleanup rules settings with a conditional calculator operation where
            it leverages these attributes and changes the cleanupGeometry operations optimally based on input features
            final output will be two file geodatabases. one with original buildings and other with cleaned.'''
        ''' Other Cleanup Utility tools/processes may be required to optimize building faces and roof geometries'''
    except:
        print("Unable to process feature {0}".format(i + beginOnFeatureNumber))
        print("Multipatch {0} Process failed @ ".format(i + beginOnFeatureNumber), str(datetime.now()))


def extractMultipatchRasterToPts(fullextent, row):
    try:
        arcpy.env.extent = fullextent
        # get raster extent
        geom = row[0]
        print("geom = ", geom)

        # copy the feature temporarily
        #tp = os.path.join("in_memory", "tp{0}".format(i))
        tp = os.path.join(outputWS, "tp{0}".format(i))
        #tp = os.path.join("in_memory", "tp{0}".format(i))
        arcpy.CopyFeatures_management(geom, tp)

        extentgeom = arcpy.Describe(tp)
        arcpy.env.mask = tp
        print("extentgeom = ", extentgeom)
        extent = "{0} {1} {2} {3}".format(extentgeom.extent.XMin, extentgeom.extent.YMin, extentgeom.extent.XMax, extentgeom.extent.YMax)
        print("extent = ", extent)
        arcpy.env.extent = extent

        print("Begin Raster Creation Process")

        DTM = os.path.join("in_memory", "DTM")
        # Delete terrain rasters if existing.
        if arcpy.Exists(DTM):
            arcpy.Delete_management(DTM)

        arcpy.LasDatasetToRaster_conversion("DTMLASD", DTM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                            valueField, "CELLSIZE", pointSpacing, heightValue)
        print("Created DTM Raster at location: " + DTM)

        # Set mask to building footprint geometry outline
        arcpy.env.mask = tp

         # clip the DTM
        print('clipping DTM')
        dtmClipRast = os.path.join("in_memory", 'tempDEMclip{0}'.format(i + beginOnFeatureNumber))
        arcpy.Clip_management(DTM, extent, dtmClipRast, tp, "true", "false")
        # convert DEM to Int
        #dtmClipRastInt = Int(dtmClipRast)

        # add Min Height to Building Footprints
        print('determining Minimum Building Elevation')
        arcpy.AddField_management(tp, "ID", "SHORT", None, None, None, "ID", "true", "true", None)
        arcpy.CalculateField_management(tp, "ID", 1, "PYTHON_9.3", None)
        minMaxElevTable = os.path.join("in_memory", "minMaxElevTable")
        arcpy.sa.ZonalStatisticsAsTable(tp, "ID", DTM, minMaxElevTable, "true", "MIN_MAX_MEAN")
        arcpy.JoinField_management(tp, "ID", minMaxElevTable, "ID", "MIN;MAX")

        # Delete the DTM Raster
        if arcpy.Exists(DTM):
            arcpy.Delete_management(DTM)

        # then, move building footprints to MIN Z Height
        #out_poly3d = os.path.join("in_memory", "out_poly3d")
        out_poly3d = os.path.join("in_memory", "out_poly3d_{0}".format(i))
        arcpy.FeatureTo3DByAttribute_3d(tp, out_poly3d, "MIN", "")

        # Create DSM Raster
        LrDSM = os.path.join("in_memory", "LrDSM{0}".format(i))
        # Delete terrain rasters if existing.
        if arcpy.Exists(LrDSM):
            arcpy.Delete_management(LrDSM)

        if optimizeRaster:
            arcpy.LasDatasetToRaster_conversion("LRDSMLASD", LrDSM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                        valueField, "CELLSIZE", pointSpacing * optimizeRasterFactor, heightValue)
            print("Created LrDSM Raster at location: " + LrDSM)

            interpolateBetweenLasPts(LrDSM)

        if not optimizeRaster:
            arcpy.LasDatasetToRaster_conversion("LRDSMLASD", LrDSM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                        valueField, "CELLSIZE", pointSpacing, heightValue)

        #nbr = arcpy.sa.NbrRectangle(3, 3, "CELL")

        # Filter LrDSM PointCloud to remove Artifacts
        LrDSMFilter = os.path.join("in_memory", "LrDSMFilter{0}".format(i))
        if arcpy.Exists(LrDSMFilter):
            arcpy.Delete_management(LrDSMFilter)

        filterOut = arcpy.sa.Filter(LrDSM, "LOW", "DATA")
        filterOut.save(LrDSMFilter)

        if arcpy.Exists(LrDSM):
            arcpy.Delete_management(LrDSM)

        # Remove Artifact Pits
        fillRaster = os.path.join("in_memory", "fillRaster{0}".format(i))
        if arcpy.Exists(fillRaster):
            arcpy.Delete_management(fillRaster)
        outputFill = arcpy.sa.Fill(LrDSMFilter, None)
        outputFill.save(fillRaster)

        if arcpy.Exists(LrDSMFilter):
            arcpy.Delete_management(LrDSMFilter)

        # slopedAreaRasters(SlopeRaster, slopedAreasNullRaster):

        # convert raster to points
        print('converting raster to points')
        extractedPoints = os.path.join("in_memory", "extractedPoints{0}".format(i))
        if arcpy.Exists(extractedPoints):
            arcpy.Delete_management(extractedPoints)
        arcpy.RasterToPoint_conversion(fillRaster, extractedPoints, "Value")

        if arcpy.Exists(fillRaster):
            arcpy.Delete_management(fillRaster)

        # convert points to 3D
        roofPoints = os.path.join("in_memory", "roofPoints{0}".format(i))
        arcpy.FeatureTo3DByAttribute_3d(in_features=extractedPoints, out_feature_class=roofPoints,
                                        height_field="grid_code", to_height_field="")

        if arcpy.Exists(extractedPoints):
            arcpy.Delete_management(extractedPoints)

        # Interpolate Points to Boundary
        buildingInsideAndBorderPoints = os.path.join("in_memory", "buildingBorderPoints{0}".format(i))
        interpolatePointsToBoundary(input_bldg_points=roofPoints, input_bldg_fp=tp,
                                    output_bldg_points_with_border=buildingInsideAndBorderPoints)

        if arcpy.Exists(roofPoints):
            arcpy.Delete_management(roofPoints)
        roofPoints = os.path.join("in_memory", "roofPoints{0}".format(i))

        arcpy.Dissolve_management(buildingInsideAndBorderPoints, roofPoints, None, None, "true", "false")
        print("Dissolved Points")

        #if arcpy.Exists(buildingInsideAndBorderPoints):
        #    arcpy.Delete_management(buildingInsideAndBorderPoints)

        # TODO: Resolve issue where roof-tin will not process
        # Check to ensure that paths can have spaces.  may be the problem.
        roofTin = os.path.join(tempFolder, "roofTin{0}".format(i))
        arcpy.CreateTin_3d(roofTin, sr, "{0} Shape.Z Mass_Points <None>".format(roofPoints), "DELAUNAY")
        print("roof Tin Created")

        #if arcpy.Exists(roofPoints):
        #    arcpy.Delete_management(roofPoints)

        # make ground TIN
        gnd_feats_tin = "{} Shape.Z Hard_Clip <None>;".format(out_poly3d)
        out_gnd_tin = os.path.join(tempFolder, "gndTin")
        arcpy.CreateTin_3d(out_gnd_tin, sr, gnd_feats_tin, "DELAUNAY")

        # extrude polygon between TINs
        print('creating Multipatch')
        this_MP = os.path.join(outputWS, "bldgMP_{0}".format(i))
        arcpy.ExtrudeBetween_3d(roofTin, out_gnd_tin, out_poly3d, this_MP)

        # add feature name to list
        mp_list.append(this_MP)

        # Delete Unnecessary files
        arcpy.Delete_management(tp)
        arcpy.Delete_management(minMaxElevTable)
        arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(dtmClipRast)
        arcpy.Delete_management(out_gnd_tin)
        arcpy.Delete_management(roofTin)
        arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(buildingInsideAndBorderPoints)

        print("Multipatch {0} Process complete @ ".format(i + beginOnFeatureNumber), str(datetime.now()))


        # TODO Geoff7015 Incorporate Cleanup Building CGA from Geof7015 rule into tool here:
        ''' every multipatch building must have LiDAR point spacing as a attribute and "Units: feet/meters
            will need to update CGA cleanup rules settings with a conditional calculator operation where
            it leverages these attributes and changes the cleanupGeometry operations optimally based on input features
            final output will be two file geodatabases. one with original buildings and other with cleaned.'''
        ''' Other Cleanup Utility tools/processes may be required to optimize building faces and roof geometries'''
    except:
        print("Unable to process feature {0}".format(i + beginOnFeatureNumber))
        print("Multipatch {0} Process failed @ ".format(i + beginOnFeatureNumber), str(datetime.now()))

def extractMultipatchFromRaster(fullextent, row):
    try:
        arcpy.env.extent = fullextent
        # get raster extent
        geom = row[0]
        print("geom = ", geom)

        # copy the feature temporarily
        #tp = os.path.join("in_memory", "tp{0}".format(i))
        tp = os.path.join(outputWS, "tp{0}".format(i))
        #tp = os.path.join("in_memory", "tp{0}".format(i))
        arcpy.CopyFeatures_management(geom, tp)

        extentgeom = arcpy.Describe(tp)
        arcpy.env.mask = tp
        print("extentgeom = ", extentgeom)
        extent = "{0} {1} {2} {3}".format(extentgeom.extent.XMin, extentgeom.extent.YMin, extentgeom.extent.XMax, extentgeom.extent.YMax)
        print("extent = ", extent)
        arcpy.env.extent = extent

        print("Begin Raster Creation Process")

        DTM = DTMRaster

        # Set mask to building footprint geometry outline
        arcpy.env.mask = tp

        # clip the DTM
        print('clipping DTM')
        dtmClipRast = os.path.join("in_memory", 'tempDEMclip{0}'.format(i + beginOnFeatureNumber))
        arcpy.Clip_management(DTM, extent, dtmClipRast, tp, "true", "false")

        # add Min Height to Building Footprints
        print('determining Minimum Building Elevation')
        arcpy.AddField_management(tp, "ID", "SHORT", None, None, None, "ID", "true", "true", None)
        arcpy.CalculateField_management(tp, "ID", 1, "PYTHON_9.3", None)
        minMaxElevTable = os.path.join("in_memory", "minMaxElevTable")
        arcpy.sa.ZonalStatisticsAsTable(tp, "ID", DTM, minMaxElevTable, "true", "MIN_MAX_MEAN")
        arcpy.JoinField_management(tp, "ID", minMaxElevTable, "ID", "MIN;MAX")

        # Delete the DTM Raster
        if arcpy.Exists(dtmClipRast):
            arcpy.Delete_management(dtmClipRast)

        # then, move building footprints to MIN Z Height
        #out_poly3d = os.path.join("in_memory", "out_poly3d")
        out_poly3d = os.path.join("in_memory", "out_poly3d_{0}".format(i))
        arcpy.FeatureTo3DByAttribute_3d(tp, out_poly3d, "MIN", "")

        # Create DSM Raster
        LrDSM = DSMRaster

        # clip the DTM
        print('clipping DTM')
        LrDSMClipRast = os.path.join("in_memory", 'tempLrDSMclip{0}'.format(i + beginOnFeatureNumber))
        arcpy.Clip_management(LrDSM, extent, dtmClipRast, tp, "true", "false")

        # Filter LrDSM PointCloud to remove Artifacts
        LrDSMFilter = os.path.join("in_memory", "LrDSMFilter{0}".format(i))
        if arcpy.Exists(LrDSMFilter):
            arcpy.Delete_management(LrDSMFilter)

        filterOut = arcpy.sa.Filter(LrDSMClipRast, "LOW", "DATA")
        filterOut.save(LrDSMFilter)

        if arcpy.Exists(LrDSMClipRast):
            arcpy.Delete_management(LrDSMClipRast)

        # Remove Artifact Pits
        fillRaster = os.path.join("in_memory", "fillRaster{0}".format(i))
        if arcpy.Exists(fillRaster):
            arcpy.Delete_management(fillRaster)
        outputFill = arcpy.sa.Fill(LrDSMFilter, None)
        outputFill.save(fillRaster)

        if arcpy.Exists(LrDSMFilter):
            arcpy.Delete_management(LrDSMFilter)

        # slopedAreaRasters(SlopeRaster, slopedAreasNullRaster):

        # convert raster to points
        print('converting raster to points')
        extractedPoints = os.path.join("in_memory", "extractedPoints{0}".format(i))
        if arcpy.Exists(extractedPoints):
            arcpy.Delete_management(extractedPoints)
        arcpy.RasterToPoint_conversion(fillRaster, extractedPoints, "Value")

        if arcpy.Exists(fillRaster):
            arcpy.Delete_management(fillRaster)

        # convert points to 3D
        roofPoints = os.path.join("in_memory", "roofPoints{0}".format(i))
        arcpy.FeatureTo3DByAttribute_3d(in_features=extractedPoints, out_feature_class=roofPoints,
                                        height_field="grid_code", to_height_field="")

        if arcpy.Exists(extractedPoints):
            arcpy.Delete_management(extractedPoints)

        # Interpolate Points to Boundary
        buildingInsideAndBorderPoints = os.path.join("in_memory", "buildingBorderPoints{0}".format(i))
        interpolatePointsToBoundary(input_bldg_points=roofPoints, input_bldg_fp=tp,
                                    output_bldg_points_with_border=buildingInsideAndBorderPoints)

        if arcpy.Exists(roofPoints):
            arcpy.Delete_management(roofPoints)
        roofPoints = os.path.join("in_memory", "roofPoints{0}".format(i))

        arcpy.Dissolve_management(buildingInsideAndBorderPoints, roofPoints, None, None, "true", "false")
        print("Dissolved Points")

        #if arcpy.Exists(buildingInsideAndBorderPoints):
        #    arcpy.Delete_management(buildingInsideAndBorderPoints)

        # TODO: Resolve issue where roof-tin will not process
        # Check to ensure that paths can have spaces.  may be the problem.
        roofTin = os.path.join(tempFolder, "roofTin{0}".format(i))
        arcpy.CreateTin_3d(roofTin, sr, "{0} Shape.Z Mass_Points <None>".format(roofPoints), "DELAUNAY")
        print("roof Tin Created")

        #if arcpy.Exists(roofPoints):
        #    arcpy.Delete_management(roofPoints)

        # make ground TIN
        gnd_feats_tin = "{} Shape.Z Hard_Clip <None>;".format(out_poly3d)
        out_gnd_tin = os.path.join(tempFolder, "gndTin")
        arcpy.CreateTin_3d(out_gnd_tin, sr, gnd_feats_tin, "DELAUNAY")

        # extrude polygon between TINs
        print('creating Multipatch')
        this_MP = os.path.join(outputWS, "bldgMP_{0}".format(i))
        arcpy.ExtrudeBetween_3d(roofTin, out_gnd_tin, out_poly3d, this_MP)

        # add feature name to list
        mp_list.append(this_MP)

        # Delete Unnecessary files
        arcpy.Delete_management(tp)
        arcpy.Delete_management(minMaxElevTable)
        arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(dtmClipRast)
        arcpy.Delete_management(out_gnd_tin)
        arcpy.Delete_management(roofTin)
        arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(buildingInsideAndBorderPoints)

        print("Multipatch {0} Process complete @ ".format(i + beginOnFeatureNumber), str(datetime.now()))


        # TODO Geoff7015 Incorporate Cleanup Building CGA from Geof7015 rule into tool here:
        ''' every multipatch building must have LiDAR point spacing as a attribute and "Units: feet/meters
            will need to update CGA cleanup rules settings with a conditional calculator operation where
            it leverages these attributes and changes the cleanupGeometry operations optimally based on input features
            final output will be two file geodatabases. one with original buildings and other with cleaned.'''
        ''' Other Cleanup Utility tools/processes may be required to optimize building faces and roof geometries'''
    except:
        print("Unable to process feature {0}".format(i + beginOnFeatureNumber))
        print("Multipatch {0} Process failed @ ".format(i + beginOnFeatureNumber), str(datetime.now()))
##############
# Begin Code #
##############
print("Starting Process at ", str(datetime.now()))
# If LiDAR Input is a LASD DataSet then count number of LAS files in LAS Dataset and List LAS files as input 4 GP tools
if arcpy.Exists(inLASD):
    arcpy.AddMessage("detected LASD Dataset as input: " + inLASD)
    lasDatasetStatsText = os.path.join(tempFolder, "lasDatasetStatsText.txt")
    arcpy.LasDatasetStatistics_management(inLASD, "true", lasDatasetStatsText, "LAS_FILES", "COMMA", "DECIMAL_POINT")
    filenames = findLasDatasetStatisticsfilePaths(lasDatasetStatsText)

    if len(filenames) == 0:
        arcpy.AddMessage("1 LAS file detected in LASD DATASET")
    else:
        arcpy.AddMessage("{0} LAS files detected in LASD DATASET".format(len(filenames)))

    # Process lasList into Esri GP tool friendly input format
    newstr = str(filenames)[1:-1].replace("', ", ";")
    lasList = '"' + newstr.replace("'", "") + '"'

# If the LiDAR Input is a a single LAS then return 1 of LAS files and format file to string for GP tools input.
if inLAS.lower().endswith('.las') and ";" not in inLAS:
    arcpy.AddMessage("1 LAS file detected")
    lasList = '"' + inLAS + '"'

# If the LiDAR Input is a string of LAS files then count number of LAS files and create List LAS files input 4 GP tools.
if inLAS.lower().endswith('.las') and ";" in inLAS:
    numberLASFiles = (inLAS.count(';')+1)
    arcpy.AddMessage(str(numberLASFiles) + " LAS file detected")
    lasList = '"' + inLAS + '"'

# If the LiDAR Input is a LAS Directory then count number of LAS files and create List of LAS files as input 4 GP tools.
if os.path.isdir(inLAS):
    lasSearchPathDirectory = inLAS + "/*.las"
    for name in glob.glob(lasSearchPathDirectory):
        filename = name
        file_extension = ".las"
        filename, file_extension = os.path.splitext(name)
        if file_extension == ".las":
            # Find all LAS files in input folder. Optionally search recursively
            recursive = True
            lasList = []
            if recursive:
                for root, dirs, files in os.walk(inLAS):
                    for file in files:
                        if file.endswith(".las") and file not in lasList:
                            lasList.append((os.path.join(root, file)))
            else:
                for file in os.listdir(inLAS):
                    if file.endswith(".las") and file not in lasList:
                        lasList.append((os.path.join(inLAS, file)))

    # Print Number of Las files
    if len(lasList) == 0:
        arcpy.AddMessage("1 LAS file detected in Directory")
    else:
        arcpy.AddMessage("{0} LAS files detected in Directory".format(len(lasList)))

    # Process lasList into Esri GP tool friendly input format
    newstr = str(lasList)[1:-1].replace("', ", ";")
    lasList = '"' + newstr.replace("'", "") + '"'

# Convert Las file List as String and format for GP tool input
# Create LASDataset from LAS files.
if inLAS.lower().endswith('.las') or os.path.isdir(inLAS):
    createlasdataset(inLAS=inLAS, sr=sr)


DTMLASD = "DTMLASD"
LRDSMLASD = "LRDSMLASD"

if arcpy.Exists(DTMLASD):
    arcpy.Delete_management(DTMLASD)
if arcpy.Exists(LRDSMLASD):
    arcpy.Delete_management(LRDSMLASD)
arcpy.MakeLasDatasetLayer_management(inLASD, DTMLASD, str(groundClassCode), groundReturn, "", "", "", "", "", "")
arcpy.MakeLasDatasetLayer_management(inLASD, LRDSMLASD, str(buildingClassCode), buildingReturn, "", "", "", "", "", "")
if arcpy.Exists(DTMLASD):
    if arcpy.Exists(LRDSMLASD):
        arcpy.AddMessage("LASD Layers Created")
else:
    arcpy.AddMessage("Could Not Create LASD Layers")
# selector determining whether or not to interpolateAdditional points for tin creation. Helps with Terribe LiDAR.
''' if interpolateAdditionalPoints is enabled then input correct heightValue & valueField for raster processing '''
''' Determine the correct point spacing settings based on raster processing algorithm requirements '''

# TODO make point spacing support recursions. obtainLiDARInfo(inLASD, lasList)[3] is currently placeholder
''' calculate average pt spacing of LiDAR tiles building footprints intersect'''
pointSpace = obtainLiDARInfo(inLASD, lasList)[3]
pointSpacing = pointSpace * pointSpacingCorrectionFactor
heightValue = 1
valueField = "FLOAT"

result = arcpy.GetCount_management(buildingFootprints)
FootprintCount = int(result.getOutput(0))
print("number of building Footprints to process = " + str(FootprintCount))

fullextent = arcpy.Describe(buildingFootprints).extent

# create list for multiPatch features
mp_list = []

# make search cursor for footprint polygons
fields = ["SHAPE@"]
with arcpy.da.SearchCursor(buildingFootprints, fields) as sc:
    for i, row in enumerate(sc):
        if (i + beginOnFeatureNumber) < FootprintCount:
            # if i is a multiple of 50 compact the gdb
            print("on BuildingFootprint {0}".format(i + beginOnFeatureNumber) + " of " + str(FootprintCount))
            if not i % 50:
                print("Began Compacting GDB @ ", str(datetime.now()))
                arcpy.Compact_management(outputWS)
                arcpy.Compact_management(scratchGDB)
                print("Complete Compacting GDB @ ", str(datetime.now()))
            if arcpy.Exists(DTMRaster) and arcpy.Exists(DSMRaster):
                extractMultipatchRasterToPts(fullextent=fullextent, row=row)
            if rasterExtractionApproach:
                extractMultipatchFromRaster(fullextent=fullextent, row=row)
            if not rasterExtractionApproach and not arcpy.Exists(DTMRaster) and not arcpy.Exists(DSMRaster):
                extractMultipatchFromPts(fullextent=fullextent, row=row)


# merge the MultiPatches into a single FC
outputMerge = os.path.join(outputWS, 'outputMergeMP')
arcpy.Merge_management(mp_list, outputMerge)

#TODO DJARRARD: delete all buildingMP* files that exist in the output workspace
# Delete Individual Multipatch Buildings
'''if arcpy.Exists(os.path.join(outputWS, "bldgMP_0")):
    for fc in arcpy.ListFeatureClasses("bldgMP*", "MULTIPATCH", outputWS):
        arcpy.Delete_management(fc)'''

if arcpy.Exists(DTMLASD):
    arcpy.Delete_management(DTMLASD)
if arcpy.Exists(LRDSMLASD):
    arcpy.Delete_management(LRDSMLASD)

print("Finished Process at ", str(datetime.now()))

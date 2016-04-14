__author__ = 'geof7015'

import arcpy
import os
import os.path
import tempfile
import glob
from arcpy.sa import *
from datetime import datetime
import gc

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('spatial')
arcpy.CheckOutExtension('3d')

inLAS = r'E:\3D_City_Data\United States\Georgia\Athens\LiDAR'
inLASD = r'' #r'E:\3D_City_Data\United States\Georgia\Athens\New LasDataset.lasd'
buildingFootprints = r'E:\3D_City_Data\United States\Georgia\Athens\Data.gdb\BuildingFootprints_1'
sr = "PROJCS['NAD_1983_StatePlane_Georgia_West_FIPS_1002_Feet',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',2296583.333333333],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-84.16666666666667],PARAMETER['Scale_Factor',0.9999],PARAMETER['Latitude_Of_Origin',30.0],UNIT['Foot_US',0.3048006096012192]]"
outputWS = r'E:\3D_City_Data\United States\Georgia\Athens\multipatch.gdb'
scratchGDB = arcpy.env.scratchGDB
tempFolder = tempfile.mkdtemp()

beginOnFeatureNumber = 0
pointSpacingCorrectionFactor = 0.5
interpolateBetweenPoints = False  # Currently Bugged...
reduceTesselations = True
rasterExtractionApproach = True

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


def interpolateBetweenLasPts(LrDSM):
    # Run raster interpolation algorithm on LiDAR derived rasters if interpolateAdditionalPoints is True and not Recursive
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

def reduceTesselationProcess(LrDSM,SlopedAreasPolygonBuffered):
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


def extractMultipatch(fullextent, row):
    try:
        arcpy.env.extent = fullextent
        # get raster extent
        geom = row[0]
        #ext = "{0} {1} {2} {3}".format(geom.extent.XMin, geom.extent.YMin, geom.extent.XMax, geom.extent.YMax)

        # copy the feature temporarily
        tp = os.path.join(scratchGDB, "tp{0}".format(i + beginOnFeatureNumber))
        tempGeom = arcpy.CopyFeatures_management(geom, tp)

        #extentgeom = arcpy.Describe(tp)
        #extent = "{0} {1} {2} {3}".format(extentgeom.extent.XMin, extentgeom.extent.YMin, extentgeom.extent.XMax, extentgeom.extent.YMax)
        #print("Building Footprint Extent = ", extent)
        extentgeom = arcpy.Describe(tp)
        arcpy.env.mask = tp
        print("extentgeom = ", extentgeom)
        extent = "{0} {1} {2} {3}".format(extentgeom.extent.XMin, extentgeom.extent.YMin, extentgeom.extent.XMax, extentgeom.extent.YMax)
        print("extent = ", extent)
        arcpy.env.extent = extent

        print("Begin Raster Creation Process")
        LrDSM = os.path.join(tempFolder, "LrDSM.tif")
        DTM = os.path.join(tempFolder, "DTM.tif")
        # Delete terrain rasters if existing.
        if arcpy.Exists(DTM):
            arcpy.Delete_management(DTM)
        if arcpy.Exists(LrDSM):
            arcpy.Delete_management(LrDSM)
        arcpy.LasDatasetToRaster_conversion("DTMLASD", DTM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                            valueField, "CELLSIZE", pointSpacing, heightValue)
        print("Created DTM Raster at location: " + DTM)
        arcpy.LasDatasetToRaster_conversion("LRDSMLASD", LrDSM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                            valueField, "CELLSIZE", pointSpacing, heightValue)
        print("Created  Last Return DSM Raster at location: " + LrDSM)

        if interpolateBetweenPoints:
            interpolateBetweenLasPts(LrDSM=LrDSM)

        if reduceTesselations:
            SlopedAreasPolygonBuffered = os.path.join(tempFolder, "SlopedAreasPolygonBuff.shp")
            reduceTesselationProcess(LrDSM, SlopedAreasPolygonBuffered)
            arcpy.env.mask = SlopedAreasPolygonBuffered

        # smooth the DSM
        print('smoothing DSM')
        # nbr = NbrRectangle(3, 3, "CELL")

        sm_clip = FilterRasterProcess = Filter(LrDSM, "LOW", "DATA")

        # clean up clipped raster
        #arcpy.Delete_management(DSMClipRast)

        # convert raster to points
        print('converting raster to points')

        out_points = os.path.join(tempFolder, "clipPoints.shp")
        if arcpy.Exists(out_points):
            arcpy.Delete_management(out_points)
        arcpy.RasterToPoint_conversion(sm_clip, out_points, "Value")

        arcpy.env.mask = tp

        # Create TIN with points
        print('making surface TIN')
        # feats_tin = "{} Shape.Z Mass_Points <None>;".format(out_points3d)
        feats_tin = "{0} grid_code Mass_Points <None>;".format(out_points)
        out_surf_tin = os.path.join(tempFolder, "surfTin")
        if arcpy.Exists(out_surf_tin):
            arcpy.Delete_management(out_surf_tin)
        arcpy.CreateTin_3d(out_surf_tin, sr, feats_tin, 'DELAUNAY')

        # clip the DTM
        print('clipping DTM')
        dtmClipRast = os.path.join(tempFolder, 'tempDEMclip{0}.tif'.format(i + beginOnFeatureNumber))
        arcpy.Clip_management(DTM, extent, dtmClipRast, tp, "true", "false")
        # convert DEM to Int
        #dtmClipRastInt = Int(dtmClipRast)

        # add Min Height to Building Footprints
        print('determining Minimum Building Elevation')
        arcpy.AddField_management(tp, "ID", "SHORT", None, None, None, "ID", "true", "true", None)
        arcpy.CalculateField_management(tp, "ID", 1, "PYTHON_9.3", None)
        minMaxElevTable = os.path.join(scratchGDB, "minMaxElevTable")
        arcpy.sa.ZonalStatisticsAsTable(tp, "ID", LrDSM, minMaxElevTable, "true", "MIN_MAX_MEAN")
        arcpy.JoinField_management(tp, "ID", minMaxElevTable, "ID", "MIN;MAX")


        # then, move building footprints to MIN Z Height
        out_poly3d = os.path.join(scratchGDB, "out_poly3d")
        arcpy.FeatureTo3DByAttribute_3d(tp, out_poly3d, "MIN", "")

        # make ground TIN
        gnd_feats_tin = "{} Shape.Z Hard_Clip <None>;".format(out_poly3d)
        out_gnd_tin = os.path.join(tempFolder, "gndTin")
        arcpy.CreateTin_3d(out_gnd_tin, sr, gnd_feats_tin, "DELAUNAY")

        # extrude polygon between TINs
        print('creating Multipatch')
        this_MP = os.path.join(outputWS, "bldgMP_{0}".format(i + beginOnFeatureNumber))
        arcpy.ExtrudeBetween_3d(out_surf_tin, out_gnd_tin, out_poly3d, this_MP)

        # add feature name to list
        mp_list.append(this_MP)

        # Delete Unnecessary files
        arcpy.Delete_management(tp)
        arcpy.Delete_management(out_points)
        arcpy.Delete_management(minMaxElevTable)
        arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(dtmClipRast)
        arcpy.Delete_management(out_gnd_tin)
        arcpy.Delete_management(out_surf_tin)
        arcpy.Delete_management(out_poly3d)
        arcpy.Delete_management(FilterRasterProcess)
        #del nbr
        del sm_clip
        del row, tp


        # TODO Geoff7015 Incorporate Cleanup Building CGA from Geof7015 rule into tool here:
        ''' every multipatch building must have LiDAR point spacing as a attribute and "Units: feet/meters
            will need to update CGA cleanup rules settings with a conditional calculator operation where
            it leverages these attributes and changes the cleanupGeometry operations optimally based on input features
            final output will be two file geodatabases. one with original buildings and other with cleaned.'''
        ''' Other Cleanup Utility tools/processes may be required to optimize building faces and roof geometries'''
    except:
        print("Unable to process feature {0}".format(i + beginOnFeatureNumber))

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
arcpy.MakeLasDatasetLayer_management(inLASD, DTMLASD, "2", "", "", "", "", "", "", "")
arcpy.MakeLasDatasetLayer_management(inLASD, LRDSMLASD, "1;2", "Last Return", "", "", "", "", "", "")
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
            extractMultipatch(fullextent=fullextent, row=row)

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

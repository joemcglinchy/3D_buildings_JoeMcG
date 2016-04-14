__author__ = 'geof7015'
#This Script Extracts and Cleans Building Footprints from Classified LiDAR

import arcpy
import os
from arcpy.sa import *

arcpy.env.overwriteOutput = True

global inLASD

ProductionMode = True  # Set to True for ArcGIS Pro GP Tool Use
if ProductionMode:
    inLASD = arcpy.GetParameterAsText(0)
    sr = arcpy.GetParameterAsText(1)
    OutputFootprints = arcpy.GetParameterAsText(2)
    scratchGDB = arcpy.env.scratchGDB
else:
    inLASD = r'C:\Users\geof7015\PycharmProjects\testData\Boulder\las.lasd'
    sr = "PROJCS['NAD_1983_HARN_StatePlane_Colorado_North_FIPS_0501_Feet',GEOGCS['GCS_North_American_1983_HARN',DATUM['D_North_American_1983_HARN',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',3000000.000316083],PARAMETER['False_Northing',999999.999996],PARAMETER['Central_Meridian',-105.5],PARAMETER['Standard_Parallel_1',39.71666666666667],PARAMETER['Standard_Parallel_2',40.78333333333333],PARAMETER['Latitude_Of_Origin',39.33333333333334],UNIT['Foot_US',0.3048006096012192]]"
    scratchGDB = arcpy.env.scratchGDB
    OutputFootprints = os.path.join(scratchGDB, "Footprints")

########################
# Set Global Variables #
########################
global lasList
global ptFileInfoFile
global ptFileInfoList
global ptSpacing
global avgPtSpacing
global mp_list


##################
# Define Modules #
##################

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
def obtainLiDARInfo(inLASD, lasList):
    if arcpy.Exists(inLASD):
        arcpy.AddMessage("Calculating Necessary Statistics for Feature Extraction Process")
        lasDatasetStatsText = os.path.join(scratchGDB, "lasDatasetStatsText.txt")
        if arcpy.Exists(lasDatasetStatsText):
            arcpy.Delete_management(lasDatasetStatsText)
        arcpy.LasDatasetStatistics_management(inLASD, "true", lasDatasetStatsText, "LAS_FILES", "COMMA",
                                              "DECIMAL_POINT")
        ptFileInfoFile = os.path.join(scratchGDB, 'ptFileInfoFile')
        if arcpy.Exists(ptFileInfoFile):
            arcpy.Delete_management(ptFileInfoFile)
        arcpy.PointFileInformation_3d(lasList, ptFileInfoFile, "LAS", None, sr, "false", "false", "DECIMAL_POINT",
                                      "false", "false")

        rows = arcpy.SearchCursor(ptFileInfoFile, fields="FileName; Pt_Spacing; Z_Min; Z_Max",
                                  sort_fields="FileName; Pt_Spacing; Z_Min; Z_Max")
        # Iterate through the rows in the cursor and store the "FileName; Pt_Spacing; Z_Min; Z_Max"
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
            ptspacinglist = float("{0}".format(
                row.getValue("Pt_Spacing")))
            PtSpacing.append(ptspacinglist)
        print(ptFileInfoList)
        print(PtSpacing)
        avgPtSpacing = sum(PtSpacing)/float(len(PtSpacing))
        print(avgPtSpacing)
        if arcpy.Exists(ptFileInfoFile):
            arcpy.Delete_management(ptFileInfoFile)
        return ptFileInfoFile, ptFileInfoList, PtSpacing, avgPtSpacing

################
# Begin Script #
################

if arcpy.Exists(inLASD):
    lasDatasetStatsText = os.path.join(scratchGDB, "lasDatasetStatsText.txt")
    if arcpy.Exists(lasDatasetStatsText):
        arcpy.Delete_management(lasDatasetStatsText)
    arcpy.LasDatasetStatistics_management(inLASD, "true", lasDatasetStatsText, "LAS_FILES", "COMMA", "DECIMAL_POINT")
    filenames = findLasDatasetStatisticsfilePaths(lasDatasetStatsText)

    if len(filenames) == 0:
        arcpy.AddMessage("1 LAS file detected in LASD DATASET")
    else:
        arcpy.AddMessage("{0} LAS files detected in LASD DATASET".format(len(filenames)))

    # Process lasList into Esri GP tool friendly input format
    newstr = str(filenames)[1:-1].replace("', ", ";")
    lasList = '"' + newstr.replace("'", "") + '"'

    avgPtSpacing = obtainLiDARInfo(inLASD, lasList)[3]

    if arcpy.Exists("BuildingLASD"):
        arcpy.Delete_management("BuildingLASD")
    arcpy.MakeLasDatasetLayer_management(inLASD, "BuildingLASD", 6, None, "true", "true", "true", "false", None,
                                         "false")
    arcpy.AddMessage("Created LASD Layer with Building Class Only")
else:
    arcpy.AddError("LASD Dataset Does not exist or is corrupt")
    exit()

arcpy.AddMessage("Converting Point-Cloud to Raster to Extract Footprint Outlines")
BuildingPointStatsRaster = os.path.join(scratchGDB, "PtStatsRaster")
if arcpy.Exists(BuildingPointStatsRaster):
    arcpy.Delete_management(BuildingPointStatsRaster)
arcpy.LasPointStatsAsRaster_management("BuildingLASD", BuildingPointStatsRaster, "INTENSITY_RANGE", "CELLSIZE",
                                       avgPtSpacing)
arcpy.AddMessage("Building Outline Raster Created from LASD Dataset")

BCRaster = os.path.join(scratchGDB, "BCRaster")
if arcpy.Exists(BCRaster):
    arcpy.Delete_management(BCRaster)
BCRasterProcess = arcpy.sa.BoundaryClean(BuildingPointStatsRaster, "ASCEND", "true")
BCRasterProcess.save(BCRaster)
arcpy.AddMessage("Filled No-Data Holes in Building Outline Raster")

ConRaster = os.path.join(scratchGDB, "ConRaster")
if arcpy.Exists(ConRaster):
    arcpy.Delete_management(ConRaster)
conRasterProcess = Con(BCRaster, 1, 0, "Value >= 0")
conRasterProcess.save(ConRaster)
arcpy.AddMessage("Converted Building Outline Raster to Raster Value of 1")

arcpy.AddMessage("Begin Converting Building Footprint Raster to Polygon Outline")
Raster2Poly = os.path.join(scratchGDB, "Raster2Poly")
if arcpy.Exists(Raster2Poly):
    arcpy.Delete_management(Raster2Poly)
arcpy.RasterToPolygon_conversion(ConRaster, Raster2Poly, "true", "Value")

arcpy.AddMessage("Begin Regularizing Building Footprints")
BuildingFootprints = OutputFootprints
arcpy.RegularizeBuildingFootprint_3d(Raster2Poly, BuildingFootprints, "ANY_ANGLE", avgPtSpacing, avgPtSpacing, 0.25,
                                     1.5, 0.1, 1000000)


arcpy.AddMessage("Begin Removing Erroneous Geometry in Building Footprints")
BuildingFootprintsCleaned = OutputFootprints + "Cleaned"
arcpy.EliminatePolygonPart_management(BuildingFootprints, BuildingFootprintsCleaned, "AREA",
                                      "600 SquareFeet", 0, "true")

arcpy.AddMessage("Extraction Complete. Removing Intermediate Data. DO NOT STOP PROCESS")

############################
# Delete Intermediate Data #
############################

if arcpy.Exists("BuildingLASD"):
    arcpy.Delete_management("BuildingLASD")
if arcpy.Exists(BuildingPointStatsRaster):
    arcpy.Delete_management(BuildingPointStatsRaster)
if arcpy.Exists(BCRaster):
    arcpy.Delete_management(BCRaster)
if arcpy.Exists(ConRaster):
    arcpy.Delete_management(ConRaster)
if arcpy.Exists(Raster2Poly):
    arcpy.Delete_management(Raster2Poly)

arcpy.AddMessage("Intermediate Data Deleted. You are Free to Close GP Tool")

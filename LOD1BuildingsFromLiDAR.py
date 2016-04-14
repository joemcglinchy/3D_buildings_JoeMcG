__author__ = 'geof7015'

import arcpy
import os


#inLASD = r'E:\3D_City_Data\United States\North Carolina\Charlotte\AIC\InnovationCorridorLiDAR.lasd'
#outputDerivatives = r'E:\3D_City_Data\United States\North Carolina\Charlotte\AIC\Workspace.gdb'
#buildingFootprint = r'E:\3D_City_Data\United States\North Carolina\Charlotte\AIC\Data.gdb\BuildingFootprints_test'

inLASD = r'E:\3D_City_Data\United States\North Carolina\Charlotte\Source\LiDAR\Charlotte.lasd'
outputDerivatives = r'E:\3D_City_Data\United States\North Carolina\Charlotte\AIC\Workspace.gdb'
buildingFootprint = r'E:\3D_City_Data\United States\North Carolina\Charlotte\Source\BuildingFootprints\Buildings.shp'

pointSpacing = 2
pointSpacingFactor = 3

buildingClassCode = 6
buildingReturn = "Last Return"

groundClassCode = 2
groundReturn = ""

DeleteIntermediateData = False


################
# Begin Script #
################

# Create Initial Raster Derivatives.

#Specify initial LiDAR Derivative Raster Settings
DTMLASD = "DTMLASD"
LRDSMLASD = "LRDSMLASD"
heightValue = 1
valueField = "FLOAT"
pointSpacingEquation = pointSpacing * pointSpacingFactor

# Create Las Dataset Layers for Ground and Buildings as inputs to raster tools.
if arcpy.Exists(DTMLASD):
    arcpy.Delete_management(DTMLASD)
if arcpy.Exists(LRDSMLASD):
    arcpy.Delete_management(LRDSMLASD)
arcpy.MakeLasDatasetLayer_management(inLASD, DTMLASD, str(groundClassCode), groundReturn, "", "", "", "", "", "")
arcpy.MakeLasDatasetLayer_management(inLASD, LRDSMLASD, str(buildingClassCode) + ";" + str(groundClassCode), buildingReturn, "", "", "", "", "", "")

# Create Building Height Raster
arcpy.AddMessage("Beginning processing of Building Height Raster. "
                 "This process may take a few hours for Large Datasets")
BuildingHeightRaster = os.path.join(outputDerivatives, "BuildingHeightRaster")
# Delete Building Height Raster if existing
if arcpy.Exists(BuildingHeightRaster):
    arcpy.Delete_management(BuildingHeightRaster)
# arcpy.LasPointStatsAsRaster_management(LRDSMLASD, BuildingHeightRaster, "Z_RANGE", "CELLSIZE", pointSpacingEquation)
arcpy.LasDatasetToRaster_conversion(LRDSMLASD, BuildingHeightRaster, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                    valueField, "CELLSIZE", pointSpacingEquation, heightValue)
print("Created DSM/Building Height Raster")

# Remove Artifacts from Building Height Raster
arcpy.AddMessage("Removing Potential Artifacts from Building Height Raster")
BuildingHtFilter = os.path.join(outputDerivatives, "BuildingHtFilter")
Filter_raster = arcpy.sa.Filter(BuildingHeightRaster, "LOW", "true")
Filter_raster.save(BuildingHtFilter)

# Create Digital Elevation Model from LiDAR
arcpy.AddMessage("Beginning processing of Terrain Raster. "
                 "This process may take a few hours for Large Datasets")
DTM = os.path.join(outputDerivatives, "DTM")
# Delete Terrain Raster if existing
if arcpy.Exists(DTM):
    arcpy.Delete_management(DTM)
arcpy.LasDatasetToRaster_conversion(DTMLASD, DTM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR", valueField,
                                    "CELLSIZE", pointSpacingEquation, heightValue)
print("Created DTM/Terrain Raster")

# Delete LASDataset Layers to clear up space system in_memory
if arcpy.Exists(DTMLASD):
    arcpy.Delete_management(DTMLASD)
if arcpy.Exists(LRDSMLASD):
    arcpy.Delete_management(LRDSMLASD)

# Begin Data Attribution Process

# Copy building footprints to GDB
''' attempting to keep all mathematical processes in table in_memory before joining to footprint at end of process'''
buildingFootprintsCopy = os.path.join(outputDerivatives, "BuildingFootprintsCopy")
if arcpy.Exists(buildingFootprintsCopy):
    arcpy.Delete_management(buildingFootprintsCopy)
arcpy.CopyFeatures_management(buildingFootprint, buildingFootprintsCopy)

def createID(buildingFootprints):
    # Determine if "ID" Field Exists

    #List of field names to be added
    to_add = ["ID"]

    #Create a list of existing field names
    fieldList = arcpy.ListFields(buildingFootprintsCopy)
    fieldName = [f.name for f in fieldList]

    for field in to_add:
      if field in fieldName:
          print("ID Already Exists!")
      else:
        # Add Height Attributes to Building Footprints
        arcpy.AddMessage("Calculating Minimum Building Elevations for Building Footprints")
        arcpy.AddField_management(buildingFootprintsCopy, "ID", "LONG", None, None, None, "ID", "true", "true", None)
        arcpy.CalculateField_management(buildingFootprintsCopy, "ID", "autoIncrement()", "PYTHON_9.3", r"rec=0\ndef autoIncrement():\n global rec\n pStart = 1 #adjust start value, if req'd \n pInterval = 1 #adjust interval value, if req'd\n if (rec == 0): \n  rec = pStart \n else: \n  rec = rec + pInterval \n return rec")


createID(buildingFootprints=buildingFootprintsCopy)

arcpy.AddMessage("Calculating Ground Elevation for Building Footprints")
# Calculate Min Elevation
ElevTable = os.path.join("in_memory", "minMaxElevTable")
arcpy.sa.ZonalStatisticsAsTable(buildingFootprintsCopy, "ID", DTM, ElevTable, "true", "MINIMUM")
arcpy.AddField_management(ElevTable, "baseElevation", "DOUBLE", None, None, None, None, "true", "false", None)
arcpy.CalculateField_management(ElevTable, "baseElevation", "!MIN!", "PYTHON_9.3", None)

arcpy.AddMessage("Calculating Building Roof Height for Building Footprints")

# Calculate MaxElevation
MaxElevTable = os.path.join("in_memory", "MaxElevTable")
arcpy.sa.ZonalStatisticsAsTable(buildingFootprintsCopy, "ID", BuildingHtFilter, MaxElevTable, "true", "MAXIMUM")
arcpy.JoinField_management(ElevTable, "ID", MaxElevTable, "ID", "MAX")
arcpy.AddField_management(ElevTable, "MaxElevation", "DOUBLE", None, None, None, None, "true", "false", None)
arcpy.CalculateField_management(ElevTable, "MaxElevation", "!MAX!", "PYTHON_9.3", None)

if arcpy.Exists(MaxElevTable):
    arcpy.Delete_management(MaxElevTable)

# Calculate totalHeight
arcpy.AddField_management(ElevTable, "totalHeight", "DOUBLE", None, None, None, None, "true", "false", None)
arcpy.CalculateField_management(ElevTable, "totalHeight", "!MaxElevation! - !baseElevation!", "PYTHON_9.3", None)

# Merge Calculations to Building Footprints
arcpy.JoinField_management(buildingFootprintsCopy, "ID", ElevTable, "ID", "baseElevation;MaxElevation;totalHeight")

if arcpy.Exists(ElevTable):
        arcpy.Delete_management(ElevTable)

arcpy.AddMessage("Z Enabling Building Footprints")
# then, move building footprints to MIN Z Height
BuildingFootprintsFinal = os.path.join(outputDerivatives, "buildingFootprints")
if arcpy.Exists(BuildingFootprintsFinal):
    arcpy.Delete_management(BuildingFootprintsFinal)
arcpy.FeatureTo3DByAttribute_3d(buildingFootprintsCopy, BuildingFootprintsFinal, "baseElevation", "")


# Delete Intermediate Data
if DeleteIntermediateData:
    arcpy.AddMessage("Deleting Intermediate Data")
    if arcpy.Exists(DTM):
        arcpy.Delete_management(DTM)
    if arcpy.Exists(BuildingHeightRaster):
        arcpy.Delete_management(BuildingHeightRaster)
    if arcpy.Exists(BuildingHtFilter):
        arcpy.Delete_management(BuildingHtFilter)
    if arcpy.Exists(buildingFootprintsCopy):
        arcpy.Delete_management(buildingFootprintsCopy)

arcpy.AddMessage("Process Complete")
arcpy.AddMessage("Building Footprints Saved to: " + BuildingFootprintsFinal)
if not DeleteIntermediateData:
    arcpy.AddMessage("Derivative Data Saved to: " + outputDerivatives)

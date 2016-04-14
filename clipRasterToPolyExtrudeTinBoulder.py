# this script will create MultiPatch feature classes from a polygon and DSM/DTM Raster.
# These rasters acn also be automatically derived from LAS or LASD in this process.
# Original Process & concept by Joseph McGlinchy & (intern name)
# LAS & LASD Integration, interpolateAdditionalPoints algorithm, & optimization by Geoff Taylor
# Enjoy :)

import arcpy
import os
import os.path
import tempfile
import glob
from arcpy.sa import *

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('spatial')
arcpy.CheckOutExtension('3d')

#################
# Define Inputs #
#################

ProductionMode = False  # Set to True for ArcGIS Pro GP Tool Use
if ProductionMode:
    ''' For ArcGIS PRO GP Tool '''
    inLASD = arcpy.GetParameterAsText(0)
    inLAS = arcpy.GetParameterAsText(1)
    DTMRaster = arcpy.GetParameterAsText(2)
    DSMRaster = arcpy.GetParameterAsText(3)
    buildingFootprints = arcpy.GetParameterAsText(4)
    sr = arcpy.GetParameterAsText(5)  # Spatial Reference
    outputWS = arcpy.GetParameterAsText(6)
    interpolateAdditionalPoints = arcpy.GetParameterAsText(7)
    recursivelyCreateAndClipRastersFromLasd = arcpy.GetParameterAsText(8)
    scratchGDB = arcpy.env.scratchGDB
    tempFolder = tempfile.mkdtemp()

else:
    ''' For Testing Purposes'''  # Comment out inputs
    inLASD = r'' #C:\Users\geof7015\PycharmProjects\LiDARTestData\testData\LiDAR\backup\LASD\test.lasd' #C:\Users\geof7015\PycharmProjects\LiDARTestData\testData\LiDAR\backup\LASD\test.lasd'  # C:\Users\geof7015\PycharmProjects\testData\Charlotte\LiDAR\LiDARAOI.lasd'
    inLAS = r'C:\Users\geof7015\PycharmProjects\testData\Boulder\LAS\329.las'     #C:\Users\geof7015\PycharmProjects\testData\Charlotte\LiDAR'  #C:\Users\geof7015\PycharmProjects\LiDARTestData\testData\LiDAR' # C:\Users\geof7015\PycharmProjects\LiDARTestData\testData\LiDAR13810E4750N.las # C:\Users\geof7015\PycharmProjects\LiDARTestData\testData #  C:\Users\geof7015\PycharmProjects\LiDARTestData\testData\LiDAR13810E4750N.las #  r'C:\Users\geof7015\PycharmProjects\testData\Charlotte\LiDAR'  # C:\Users\geof7015\PycharmProjects\testData\Charlotte\LiDAR
    DTMRaster = r''  # C:\workspace\data\testdata\bh12TVK1800084000.img
    DSMRaster = r''  # 'C:\workspace\data\testdata\hh12TVK1800084000.img'
    buildingFootprints = r'C:\Users\geof7015\PycharmProjects\testData\Boulder\Data.gdb\Building329'     # C:\Users\geof7015\PycharmProjects\testData\Charlotte\Data.gdb\BldgFootprints'
    sr = "PROJCS['NAD_1983_HARN_StatePlane_Colorado_North_FIPS_0501_Feet',GEOGCS['GCS_North_American_1983_HARN',DATUM['D_North_American_1983_HARN',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',3000000.000316083],PARAMETER['False_Northing',999999.999996],PARAMETER['Central_Meridian',-105.5],PARAMETER['Standard_Parallel_1',39.71666666666667],PARAMETER['Standard_Parallel_2',40.78333333333333],PARAMETER['Latitude_Of_Origin',39.33333333333334],UNIT['Foot_US',0.3048006096012192]]"
    # "PROJCS['NAD_1983_StatePlane_North_Carolina_FIPS_3200_Feet',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',2000000.002616666],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-79.0],PARAMETER['Standard_Parallel_1',34.33333333333334],PARAMETER['Standard_Parallel_2',36.16666666666666],PARAMETER['Latitude_Of_Origin',33.75],UNIT['Foot_US',0.3048006096012192]]"
    outputWS = r'C:\Users\geof7015\PycharmProjects\testData\Boulder\Workspace.gdb'
    scratchGDB = arcpy.env.scratchGDB
    tempFolder = tempfile.mkdtemp()
    # TODO Geof7015 resolve bug with Boundary Clean Raster Int to Float issue!
    # interpolateAdditionalPoints currently broken
    interpolateAdditionalPoints = False
    recursivelyCreateAndClipRastersFromLasd = False

########################
# Set Global Variables #
########################

global lasList
global DTM
global DSM
global LrDSM
global ptFileInfoFile
global ptFileInfoList
global ptSpacing
global avgPtSpacing

##################
# Define Modules #
##################
''' place all code modules here and notate each well...'''

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
        if recursivelyCreateAndClipRastersFromLasd:
            pass

        if not recursivelyCreateAndClipRastersFromLasd:
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
                ptspacinglist = float("{0}".format(
                    row.getValue("Pt_Spacing")))
                PtSpacing.append(ptspacinglist)
            print(ptFileInfoList)
            print(PtSpacing)
            avgPtSpacing = sum(PtSpacing)/float(len(PtSpacing))
            print(avgPtSpacing)
            return ptFileInfoFile, ptFileInfoList, PtSpacing, avgPtSpacing

# TODO Geof7015 & DJARRARD Integrate automated class code detection for feature extraction process for .las or .lasd use

def createlaslayers(inLASD):
    # Something like def createLasDatasetLayer(inLASD, #DTMLASDclasses, #DTMLASDreturns, #LrDSMLASDclasses, #LrDSMLASDreturns, interpolateAdditionalPoints):
    if arcpy.Exists("DTMLASD"):
        arcpy.Delete_management("DTMLASD")
    if arcpy.Exists("LRDSMLASD"):
        arcpy.Delete_management("LRDSMLASD")
    arcpy.MakeLasDatasetLayer_management(inLASD, "DTMLASD", "2", "", "", "", "", "", "", "")
    arcpy.MakeLasDatasetLayer_management(inLASD, "LRDSMLASD", "2;6", "", "", "", "", "", "", "")
    if arcpy.Exists("DTMLASD"):
        if arcpy.Exists("LRDSMLASD"):
            arcpy.AddMessage("LASD Layers Created")
    else:
        arcpy.AddMessage("Could Not Create LASD Layers")

def createSurfaceRasters(lasdLayerGround, lasdLayerSurface, outputDTM, outputDSM):
    # selector determining whether or not to interpolateAdditional points for tin creation. Helps with Terribe LiDAR.
    ''' if interpolateAdditionalPoints is enabled then input correct heightValue & valueField for raster processing '''
    ''' Determine the correct point spacing settings based on raster processing algorithm requirements '''
    if interpolateAdditionalPoints:
        if not recursivelyCreateAndClipRastersFromLasd:
            pointSpacing = obtainLiDARInfo(inLASD, lasList)[3]  # return 3 is Average LiDAR point Spacing
        else:
            # TODO make point spacing support recursions. obtainLiDARInfo(inLASD, lasList)[3] is currently placeholder
            ''' calculate average pt spacing of LiDAR tiles building footprints intersect'''
            pointSpace = obtainLiDARInfo(inLASD, lasList)[3]
            pointSpacing = pointSpace / 0.5
        heightValue = 100
        valueField = "INT"
    else:
        if not recursivelyCreateAndClipRastersFromLasd:
            pointSpacing = obtainLiDARInfo(inLASD, lasList)[3]
        else:
            # TODO make point spacing support recursions.  obtainLiDARInfo(inLASD, lasList)[3] is currently placeholder
            ''' calculate average pt spacing of LiDAR tiles building footprints intersect'''
            pointSpacing = obtainLiDARInfo(inLASD, lasList)[3]
        heightValue = 1
        valueField = "FLOAT"

    # Delete terrain rasters if existing.
    if arcpy.Exists(outputDTM):
        arcpy.Delete_management(outputDTM)
    if arcpy.Exists(outputDSM):
        arcpy.Delete_management(outputDSM)

    arcpy.LasDatasetToRaster_conversion(lasdLayerGround, outputDTM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                        valueField, "CELLSIZE", pointSpacing, heightValue)
    print("Created DTM Raster at location: " + outputDTM)
    arcpy.LasDatasetToRaster_conversion(lasdLayerSurface, outputDSM, "ELEVATION", "BINNING MAXIMUM NATURAL_NEIGHBOR",
                                        valueField, "CELLSIZE", pointSpacing, heightValue)
    print("Created  Last Return DSM Raster at location: " + outputDSM)
    return outputDTM, outputDSM


# TODO INT.  make sure that integer raster is being produced
def interpolateBetweenLasPts(DTM,LrDSM):
    # Raster Cleanup Process via Boundary Clean Operation (Fills & Smooths Voids based on nearest pixel value)
    DTMBC = os.path.join(scratchGDB, "DTMBC")
    DTMBoundaryClean = arcpy.sa.BoundaryClean(DTM, "ASCEND", "TWO_WAY")
    DTMBoundaryClean.save(DTMBC)
    LrDSMBC = os.path.join(scratchGDB, "LrDSMBC")
    LrDSMBoundaryClean = arcpy.sa.BoundaryClean(LrDSM, "ASCEND", "TWO_WAY")
    LrDSMBoundaryClean.save(LrDSMBC)

    # TODO INT.  check and see id produced raster is created correctly
    # Divide rasters by 1000 for accurate Height & Override the Original DEM & DSM Raster Inputs with the LiDAR Ext ones
    DTMFinal = arcpy.Raster(DTMBC) / 100  # May need to resolve output rasters format to true Float
    DTM = os.path.join(scratchGDB, "DTM")
    DTMFinal.save(DTM)
    LrDSMFinal = arcpy.Raster(LrDSMBC) / 100  # May need to resolve output rasters format to true Float
    DSM = os.path.join(scratchGDB, "DSM")
    LrDSMFinal.save(DSM)
    print("Complete Creation of Interpolation Data")
    return DTM, DSM

def createNDSM(DSM, DTM, nDSM):
    temp = arcpy.Raster(DSM) - arcpy.Raster(DTM)
    # TODO Geof7015 obtain feet or meters from SR as input for automated unit selector function
    ''' then pass conversion equation case units == meters: minBldgHeight * 0.3048 else: minBldgHeight '''
    minBldgHeight = 6   # Value is in Feet
    nDSMRaster = SetNull(temp < minBldgHeight, temp)
    nDSMRaster.save(nDSM)
    del temp
    return nDSM

def maskDSM(DSM,maskRaster,DSMMasked):
    # Mask DSM to nDSM height limiting value
    DSMMaskOperation = arcpy.sa.ExtractByMask(DSM, maskRaster)
    DSMMaskOperation.save(DSMMasked)
    print("DSMMasked File Location" + DSMMasked)
    print("DTM File Location" + DTM)
    return DSMMasked

##############
# Begin Code #
##############
''' commence Code/Script operation here and notate well...'''

# Conditional operation to allow for input of .lasd file, .las file folder or Rasters
''' Detects if  dataSet exists and runs the correct operation based on dataSet'''

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

if inLAS.lower().endswith('.las') or os.path.isdir(inLAS) or arcpy.Exists(inLASD):
    createlaslayers(inLASD=inLASD)
    if not recursivelyCreateAndClipRastersFromLasd:
        LrDSM = os.path.join(scratchGDB, "LrDSM")
        DTM = os.path.join(scratchGDB, "DTM")
        createSurfaceRasters(lasdLayerGround="DTMLASD", lasdLayerSurface="LRDSMLASD", outputDTM=DTM, outputDSM=LrDSM)
        arcpy.Delete_management("DTMLASD")
        arcpy.Delete_management("LRDSMLASD")
    else:
        pass

# TODO INT Make sure that INT raster is produced
# Run raster interpolation algorithm on LiDAR derived rasters if interpolateAdditionalPoints is True and not Recursive
if inLAS.lower().endswith('.las') or os.path.isdir(inLAS) or arcpy.Exists(inLASD):
    if interpolateAdditionalPoints and not recursivelyCreateAndClipRastersFromLasd:
            interpolateBetweenLasPts(DTM, LrDSM)
    else:
        DSM = os.path.join(outputWS, "LrDSM")

# Process Rasters prior to recursion process if Raster as user Input or recursivelyCreateAndClipRastersFromLasd deselect
if (arcpy.Exists(DTMRaster) and arcpy.Exists(DSMRaster)) or not recursivelyCreateAndClipRastersFromLasd:

    # Change names of initial input rasters to DTM and DSM to align with LiDAR derived rasters for less coding :)
    if arcpy.Exists(DTMRaster) and arcpy.Exists(DSMRaster):
        DTM = DTMRaster
        DSM = DSMRaster

    # Begin process of deriving rasters from created LAS DataSets
    if inLAS.lower().endswith('.las') or os.path.isdir(inLAS) or arcpy.Exists(inLASD):
        pass
        # get spatial reference
        sr = arcpy.Describe(DTM).spatialReference

        DSM = LrDSM
        nDSM = os.path.join(scratchGDB, "nDSM")
        createNDSM(DSM=DSM, DTM=DTM, nDSM=nDSM)

        DSMMaskedRaster = os.path.join(scratchGDB, "DSMMaskedRaster")
        maskDSM(DSM=DSM, maskRaster=nDSM, DSMMasked=DSMMaskedRaster)

        DSM = DSMMaskedRaster

# create list for multiPatch features
mp_list = []

# make search cursor for footprint polygons
fields = ["SHAPE@"]
with arcpy.da.SearchCursor(buildingFootprints, fields) as sc:
    for i, row in enumerate(sc):
        try:
            print("on feature {0}".format(i))
            # get raster extent
            geom = row[0]
            ext = "{0} {1} {2} {3}".format(geom.extent.XMin, geom.extent.YMin, geom.extent.XMax, geom.extent.YMax)

            # copy the feature temporarily
            tp = os.path.join(outputWS, "tp")
            tempGeom = arcpy.CopyFeatures_management(geom, tp)

            # clip the DSM
            print('clipping DSM')
            DSMClipRast = os.path.join(arcpy.env.scratchFolder, 'tempclip{0}.tif'.format(i))
            arcpy.Clip_management(DSM, ext, DSMClipRast, tp, "true", "false")

            if interpolateAdditionalPoints:
                # Int Raster
                DSMClipRast = Int(DSMClipRast)

            # smooth the DSM
            print('smoothing DSM')
            nbr = NbrRectangle(3, 3, "CELL")

            # sm_clip = FocalStatistics(out_raster, nbr, "MEAN", "DATA")
            sm_clip = Filter(DSMClipRast, "LOW", "DATA")

            # clean up clipped raster
            arcpy.Delete_management(DSMClipRast)

            # TODO INT Make sure that Raster Points are produced correctly and at precise elevations
            # convert raster to points
            print('converting raster to points')
            out_points = "in_memory/clipPoints"
            arcpy.RasterToPoint_conversion(sm_clip, out_points, "Value")

            # convert to 3d by values
            # out_points3d ="in_memory/clipPoints3d"
            # arcpy.FeatureTo3DByAttribute_3d(in_features=out_points, out_feature_class=out_points3d,
            # height_field="grid_code", to_height_field="")

            # TODO INT Make sure that TIN Surface is produced correctly and at precise elevation
            # Create TIN with points
            print('making surface TIN')
            # feats_tin = "{} Shape.Z Mass_Points <None>;".format(out_points3d)
            feats_tin = "{0} grid_code Mass_Points <None>;".format(out_points)
            out_surf_tin = os.path.join(arcpy.env.scratchFolder, "surfTin")
            arcpy.CreateTin_3d(out_surf_tin, sr, feats_tin, 'DELAUNAY')

            # clip the DTM
            print('clipping DTM')
            dtmClipRast = os.path.join(arcpy.env.scratchFolder, 'tempDEMclip{0}.tif'.format(i))
            arcpy.Clip_management(DTM, ext, dtmClipRast, tp, "true", "false")

            # TODO INT  Check to ensure whether or not this is necessary
            # convert DEM to Int
            if interpolateAdditionalPoints:
                dtmClipRast = Int(dtmClipRast)

            # add Min Height to Building Footprints
            print('determining Minimum Building Elevation')
            arcpy.AddField_management(tp, "ID", "SHORT", None, None, None, "ID", "true", "true", None)
            arcpy.CalculateField_management(tp, "ID", 1, "PYTHON_9.3", None)
            minMaxElevTable = "in_memory/minMaxElev"
            arcpy.sa.ZonalStatisticsAsTable(tp, "ID", dtmClipRast, minMaxElevTable, "true", "MIN_MAX_MEAN")
            arcpy.JoinField_management(tp, "ID", minMaxElevTable, "ID", "MIN;MAX")

            # then, move building footprints to MIN Z Height
            out_poly3d = "in_memory/out_poly3d"
            arcpy.FeatureTo3DByAttribute_3d(tp, out_poly3d, "MIN", "")

            # make ground TIN
            gnd_feats_tin = "{} Shape.Z Hard_Clip <None>;".format(out_poly3d)
            out_gnd_tin = os.path.join(arcpy.env.scratchFolder, "gndTin")
            arcpy.CreateTin_3d(out_gnd_tin, sr, gnd_feats_tin, "DELAUNAY")

            # extrude polygon between TINs
            print('creating Multipatch')
            this_MP = os.path.join(outputWS, "bldgMP_{0}".format(i))
            arcpy.ExtrudeBetween_3d(out_surf_tin, out_gnd_tin, out_poly3d, this_MP)

            # add feature name to list
            mp_list.append(this_MP)

            # Delete Unnecessary files
            arcpy.Delete_management(tp)
            arcpy.Delete_management(out_points)
            arcpy.Delete_management(minMaxElevTable)
            arcpy.Delete_management(out_poly3d)

            # TODO Geoff7015 Incorporate Cleanup Building CGA from Geof7015 rule into tool here:
            ''' every multipatch building must have LiDAR point spacing as a attribute and "Units: feet/meters
                will need to update CGA cleanup rules settings with a conditional calculator operation where
                it leverages these attributes and changes the cleanupGeometry operations optimally based on input features
                final output will be two file geodatabases. one with original buildings and other with cleaned.'''
            ''' Other Cleanup Utility tools/processes may be required to optimize building faces and roof geometries'''
        except:
            print("Unable to process feature {0}".format(i))

# merge the MultiPatches into a single FC
outputMerge = os.path.join(outputWS, 'outputMergeMP')
arcpy.Merge_management(mp_list, outputMerge)

#TODO DJARRARD: delete all buildingMP* files that exist in the output workspace
# Delete Individual Multipatch Buildings
'''if arcpy.Exists(os.path.join(outputWS, "bldgMP_0")):
    for fc in arcpy.ListFeatureClasses("bldgMP*", "MULTIPATCH", outputWS):
        arcpy.Delete_management(fc)'''

if arcpy.Exists("DTMLASD"):
    arcpy.Delete_management("DTMLASD")
if arcpy.Exists("LRDSMLASD"):
    arcpy.Delete_management("LRDSMLASD")

__author__ = 'geof7015 esri'

import arcpy
import os

arcpy.env.overwriteOutput = True

# TODO Geof7015 or DJARRARD integrate the following with clipRasterToPolyExtrudeTin gp tool once complete

buildingFootprints = \
        r"C:\Users\geof7015\PycharmProjects\FeatureExtraction\testData\buildingFootprints\buildings.gdb\Buildings_Subset"

scratchGDB = r'C:\Users\geof7015\PycharmProjects\FeatureExtraction\testData\buildingFootprints\scratch.gdb'

buildings = os.path.join(scratchGDB, "Buildings")
Field = 'OBJECTID'

arcpy.AddMessage("Copying Building Footprints")
arcpy.CopyFeatures_management(buildingFootprints, buildings)

arcpy.AddMessage("Adding Necessary Fields to procedurally model 3D Buildings from")
arcpy.AddField_management(buildings, "TEXT", "TEXT", None, None, 12, "Text", "true", "false", None)

buildingList = []
buildingCursor = arcpy.SearchCursor(buildings, Field)

result = arcpy.GetCount_management(buildings)
count = int(result.getOutput(0))

arcpy.AddMessage("Building Footprints detected..." + str(count))

# Set the progress Bar (progressor)
arcpy.SetProgressor("step", "Extracting Building Info From LiDAR",
                    0, count, 1)

# for row in buildingCursor:
#    buildingList.append(row[0])

# arcpy.AddMessage("{0} Building Footprints detected...".format(len(buildingList)))

for row in buildingCursor:
    buildingList.append(row[0])

    ######################
    # Begin GP Tool Here #
    ######################

    # TODO Geof7015 or DJARRARD integrate clipRasterToPolyExtrudeTin gp tool here:

    arcpy.AddMessage("Processing Building Footprint {0} ".format(max(buildingList)))
    arcpy.CalculateField_management(buildings, "TEXT", "!OBJECTID!", "PYTHON_9.3", None)

    # Update the progressor position
    arcpy.SetProgressorPosition()

import time
# Wait for 1 second
time.sleep(1)
arcpy.AddMessage("\n" + "Complete Feature Extraction of " + str(count) + " Buildings")
# Wait for 3 seconds
time.sleep(4)
arcpy.AddMessage("\n" + "You are now a 3D Jedi!")
# Wait for 4 seconds
time.sleep(4)
arcpy.AddMessage("\n" + "What Are You Waiting For?")
# Wait for 2 seconds
time.sleep(2)
arcpy.AddMessage("\n" + "Apply a CityEngine .rpk in ArcGIS Pro to the Building Footprints")
# Wait for 2 seconds
time.sleep(2)
arcpy.AddMessage("\n" + "Or export them to CityEngine to Design in 3D Context!")

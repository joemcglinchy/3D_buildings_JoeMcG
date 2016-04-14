__author__ = 'geof7015'

import arcpy
import os
import os.path
import tempfile
from datetime import datetime

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('spatial')
arcpy.CheckOutExtension('3d')

roofPoints = r"C:\Users\geof7015\PycharmProjects\testData\Charlotte\Workspace2.gdb\singlepts"  # This is the clipped roof Points
tp = r'C:\Users\geof7015\PycharmProjects\testData\Charlotte\Workspace2.gdb\tp'  # This is the building footprint
sr = "PROJCS['NAD_1983_HARN_StatePlane_Colorado_North_FIPS_0501_Feet',GEOGCS['GCS_North_American_1983_HARN',DATUM['D_North_American_1983_HARN',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Lambert_Conformal_Conic'],PARAMETER['False_Easting',3000000.000316083],PARAMETER['False_Northing',999999.999996],PARAMETER['Central_Meridian',-105.5],PARAMETER['Standard_Parallel_1',39.71666666666667],PARAMETER['Standard_Parallel_2',40.78333333333333],PARAMETER['Latitude_Of_Origin',39.33333333333334],UNIT['Foot_US',0.3048006096012192]]"
pointSpace = 1.9  # LiDAR PointSpacing **Auto-Determined in larger script
outputWS = r'C:\Users\geof7015\PycharmProjects\testData\Charlotte\Workspace2.gdb'  # Output
scratchGDB = arcpy.env.scratchGDB
tempFolder = tempfile.mkdtemp()

####################
# Define Functions #
####################

# Automatically removes Artifacts from point-clouds.
# Resolves issues where building sides were triangulated and other geometric flaws
def cleanupArtifacts(single_bldg_pts, single_bldg_pts_cleaned):

    arcpy.Near3D_3d(single_bldg_pts, single_bldg_pts, str(1.4 * pointSpace), "LOCATION", "ANGLE", "DELTA")

    bldgpts = os.path.join(outputWS, "bldgPoints")
    arcpy.MakeFeatureLayer_management(single_bldg_pts, bldgpts)

    arcpy.SelectLayerByAttribute_management(bldgpts, "NEW_SELECTION", "NEAR_DELTZ >= -1 Or NEAR_ANG_V <> 0")

    arcpy.CopyFeatures_management(bldgpts, single_bldg_pts_cleaned)
    print("Artifacts Removed")
    return single_bldg_pts_cleaned


# assigns nearest point value using planar distance to building footprint points from an input point feature class.
def interpolatePointsToBoundary(input_bldg_points, input_bldg_fp, output_bldg_points_with_border):

    # explode input multipoint FC to single part
    # it is understood the input point FC will be multipoint. Need to convert to single part features
    #single_bldg_pts = os.path.join(outputWS, "singlepts")
    single_bldg_pts = os.path.join("in_memory", "singlepts")
    arcpy.MultipartToSinglepart_management(input_bldg_points, single_bldg_pts)
    print("createdSinglePts")

    # Cleanup Artifacts
    # single_bldg_pts_cleaned = os.path.join(outputWS, "single_bldg_pts_cleaned")
    single_bldg_pts_cleaned = os.path.join("in_memory", "single_bldg_pts_cleaned")
    cleanupArtifacts(single_bldg_pts=single_bldg_pts, single_bldg_pts_cleaned=single_bldg_pts_cleaned)

    # add geometry attributes
    arcpy.AddGeometryAttributes_management(Input_Features=single_bldg_pts_cleaned, Geometry_Properties="POINT_X_Y_Z_M",
                                           Length_Unit="", Area_Unit="", Coordinate_System="")

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

################
# Begin Script #
################

buildingInsideAndBorderPoints = os.path.join(outputWS, "buildingBorderPoints")
interpolatePointsToBoundary(input_bldg_points=roofPoints, input_bldg_fp=tp,
                            output_bldg_points_with_border=buildingInsideAndBorderPoints)

if arcpy.Exists(roofPoints):
    arcpy.Delete_management(roofPoints)

roofTin = os.path.join(tempFolder, "roofTin")
roofPtsFormula = "{0} Shape.Z Mass_Points <None>;{1} <None> Soft_Clip <None>".format(buildingInsideAndBorderPoints, tp)
print(roofPtsFormula)
arcpy.CreateTin_3d(roofTin, sr, roofPtsFormula, "false")
print("roof Tin Created")

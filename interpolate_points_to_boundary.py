#-------------------------------------------------------------------------------
# Name:        interpolate_points_to_boundary.py
# Purpose:     this script assigns nearest point value using planar distance to
#              building footprint points from an input point feature class.
#
# Author:      Joe McGlinchy
#
# Created:     21/09/2015
# Copyright:   (c) jose6641 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import arcpy
#import time
import os


arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('3D')

#inputs
input_bldg_points = r"J:\Projects\ResearchGroup\3dcities\PointsToMultipatch\PointsToMultipatch.gdb\PointsClipped1"
input_bldg_fp = r"J:\Projects\ResearchGroup\3dcities\PointsToMultipatch\PointsToMultipatch.gdb\Footprint1"
output_bldg_points_with_border = r"J:\Projects\ResearchGroup\3dcities\PointsToMultipatch\PointsToMultipatch.gdb\border_incl1"

# define in memory workspace
mem = "in_memory"

## it is understood the input point FC will be multipoint. Need to convert to
## single part features
# explode input multipoint FC to single part
single_bldg_pts = os.path.join(mem, "singlepts")
arcpy.MultipartToSinglepart_management(input_bldg_points,single_bldg_pts)

# add geometry attributes
arcpy.AddGeometryAttributes_management(Input_Features=single_bldg_pts, Geometry_Properties="POINT_X_Y_Z_M", Length_Unit="", Area_Unit="", Coordinate_System="")


## process the building footprint
# convert to line
bldg_line = os.path.join(mem, "bldgline")
arcpy.FeatureToLine_management(in_features=input_bldg_fp, out_feature_class=bldg_line, cluster_tolerance="", attributes="NO_ATTRIBUTES")

# densify
arcpy.Densify_edit(in_features=bldg_line, densification_method="DISTANCE", distance="1 Feet", max_deviation="0.33 Feet", max_angle="10")

# convert to points
bldg_ln_pts = os.path.join(mem, "bldglinepts")
arcpy.FeatureVerticesToPoints_management(in_features=bldg_line, out_feature_class=bldg_ln_pts, point_location="ALL")

# use Near tool to identify point FID from building points to the boundary points
arcpy.Near_analysis(in_features=bldg_ln_pts, near_features=single_bldg_pts, search_radius="5 Feet", location="NO_LOCATION", angle="NO_ANGLE", method="PLANAR")


# now, grab the NEARI_FID field and assign that feature's z-value to the building footprint point z value
arcpy.AddField_management(bldg_ln_pts, "z_val", "DOUBLE")
tbl_fp = arcpy.da.FeatureClassToNumPyArray(bldg_ln_pts, ["NEAR_FID"])
tbl_pts = arcpy.da.FeatureClassToNumPyArray(single_bldg_pts, ["POINT_Z"])

# update the z_val attribute
with arcpy.da.UpdateCursor(bldg_ln_pts, ["z_val"]) as sc:
    for i,row in enumerate(sc):
        #print i
        fid = tbl_fp[i][0]
        row[0] = tbl_pts[fid-1][0]

        sc.updateRow(row)

# convert to 3D and copy
arcpy.FeatureTo3DByAttribute_3d(bldg_ln_pts, output_bldg_points_with_border, "z_val")





# using time, this took about 0.07 seconds
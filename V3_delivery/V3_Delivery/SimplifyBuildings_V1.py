# ---------------------------------------------------------------------------
# Name: SimplifyBuildings.py
# Purpose:     To Simplify building footprints
# Usage:
# Description:
# Author:       Roslyn Dunn
# Organization: Esri Inc.
#
# Created:     03/11/2016 Roslyn Dunn
# ---------------------------------------------------------------------------

import sys
import os
import inspect
import traceback
import time
import shutil
import arcpy


def main(input_fc, output_fgdb, output_fc_name, delete_intermediate_results):
    try:
        start = time.time()

        if not arcpy.Exists(input_fc):
            arcpy.AddMessage("\nExiting process...Input Feature Class does not exist: \n {0}".format(input_fc))
            sys.exit()

        if not arcpy.Exists(output_fgdb):
            arcpy.AddMessage("\nExiting process...Input File GDB does not exist: \n {0}".format(input_fc))
            sys.exit()

        output_fc = os.path.join(output_fgdb, output_fc_name)
        if arcpy.Exists(output_fc):
            arcpy.AddMessage("\nExiting process...Output Feature Class already exists: \n {0}".format(output_fc))
            sys.exit()

        # Send intermediate outputs to the same GDB as the final output
        scratch_gdb_name = output_fgdb
        arcpy.AddMessage("\nScratch GDB for intermediate products: \n {0}".format(scratch_gdb_name))

        buff1_output = os.path.join(scratch_gdb_name, r"BuffInward")
        arcpy.Buffer_analysis(in_features=input_fc, out_feature_class=buff1_output,
                              buffer_distance_or_field="-2 Meters", line_side="FULL",
                              line_end_type="ROUND", dissolve_option="NONE", dissolve_field="", method="GEODESIC")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from Buffer_analysis are: \n{0}\n".format(messages))

        single_part_features = os.path.join(scratch_gdb_name, r"SinglePart")
        arcpy.MultipartToSinglepart_management(in_features=buff1_output, out_feature_class=single_part_features)
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from MultipartToSinglepart_management are: \n{0}\n".format(messages))

        buff2_output = os.path.join(scratch_gdb_name, r"BuffOutward")
        arcpy.Buffer_analysis(in_features=single_part_features, out_feature_class=buff2_output,
                              buffer_distance_or_field="2 Meters", line_side="FULL",
                              line_end_type="ROUND", dissolve_option="NONE", dissolve_field="", method="GEODESIC")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from Buffer_analysis are: \n{0}\n".format(messages))

        large_buildings = r"LargeBuildings"
        arcpy.MakeFeatureLayer_management(in_features=buff2_output, out_layer=large_buildings,
                                          where_clause="Shape_Area > 500", workspace="", field_info="#")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from MakeFeatureLayer_management are: \n{0}\n".format(messages))

        small_buildings = r"SmallBuildings"
        arcpy.MakeFeatureLayer_management(in_features=buff2_output, out_layer=small_buildings,
                                          where_clause="Shape_Area <= 500 AND Shape_Area >= 10",
                                          workspace="", field_info="#")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from MakeFeatureLayer_management are: \n{0}\n".format(messages))

        large_reg = os.path.join(scratch_gdb_name, r"LargeReg")
        arcpy.RegularizeBuildingFootprint_3d(in_features=large_buildings, out_feature_class=large_reg,
                                             method="RIGHT_ANGLES_AND_DIAGONALS", tolerance="0.75",
                                             densification="0.75", precision="0.25", diagonal_penalty="1.5",
                                             min_radius="0.1", max_radius="1000000")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from RegularizeBuildingFootprint_3d are: \n{0}\n".format(messages))

        small_reg = os.path.join(scratch_gdb_name, r"SmallReg")
        arcpy.RegularizeBuildingFootprint_3d(in_features=small_buildings, out_feature_class=small_reg,
                                             method="RIGHT_ANGLES", tolerance="0.75",
                                             densification="0.75", precision="0.25", diagonal_penalty="1.5",
                                             min_radius="0.1", max_radius="1000000")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from RegularizeBuildingFootprint_3d are: \n{0}\n".format(messages))

        merge_list = []
        merge_list.append(large_reg)
        merge_list.append(small_reg)
        merge_reg = os.path.join(scratch_gdb_name, r"MergeReg")
        arcpy.Merge_management(inputs=merge_list, output=merge_reg,
                               field_mappings="#")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from Merge_management are: \n{0}\n".format(messages))

        simplify_once = os.path.join(scratch_gdb_name, r"SimplifyOnce")
        arcpy.SimplifyBuilding_cartography(in_features=merge_reg, out_feature_class=simplify_once,
                                           simplification_tolerance="2 Meters", minimum_area="0 SquareFeet",
                                           conflict_option="NO_CHECK")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from SimplifyBuilding_cartography are: \n{0}\n".format(messages))


        arcpy.SimplifyBuilding_cartography(in_features=simplify_once, out_feature_class=output_fc,
                                           simplification_tolerance="4 Meters", minimum_area="0 SquareFeet",
                                           conflict_option="NO_CHECK")
        messages = arcpy.GetMessages()
        arcpy.AddMessage("\nResults output from SimplifyBuilding_cartography are: \n{0}\n".format(messages))

        if (delete_intermediate_results == "true"):
            # Delete Mosaic Datasets in MosaicDatasets.gdb, then delete remaining fc's and then MosaicDatasets.gdb
            if arcpy.Exists(buff1_output):
                arcpy.Delete_management(buff1_output)
            if arcpy.Exists(single_part_features):
                arcpy.Delete_management(single_part_features)
            if arcpy.Exists(buff2_output):
                arcpy.Delete_management(buff2_output)
            if arcpy.Exists(large_buildings):
                arcpy.Delete_management(large_buildings)
            if arcpy.Exists(small_buildings):
                arcpy.Delete_management(small_buildings)
            if arcpy.Exists(large_reg):
                arcpy.Delete_management(large_reg)
            if arcpy.Exists(small_reg):
                arcpy.Delete_management(small_reg)
            if arcpy.Exists(merge_reg):
                arcpy.Delete_management(merge_reg)
            if arcpy.Exists(simplify_once):
                arcpy.Delete_management(simplify_once)

        arcpy.Delete_management(large_buildings)
        arcpy.Delete_management(small_buildings)

        end = time.time()
        delta = end - start
        # This is useful if the tool is run at the command line
        arcpy.AddMessage("***** Total elapsed time is {0} hours *****".format(delta/3600))

    except arcpy.ExecuteError:
        print(arcpy.GetMessages())
    except Exception:
        # Return any Python specific errors and any error returned by the geoprocessor
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "PYTHON ERRORS:\nTraceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        msgs = "GP ERRORS:\n" + arcpy.GetMessages() + "\n"
        arcpy.AddError(msgs)

if __name__ == '__main__':
    arcpy.AddMessage(inspect.getfile(inspect.currentframe()))
    arcpy.AddMessage(os.path.dirname(inspect.getfile(inspect.currentframe())))
    arcpy.AddMessage(sys.version)
    arcpy.AddMessage(sys.executable)
    executed_from = sys.executable.upper()

    arcpy.CheckOutExtension('3D')
    # if "ARCGISPRO" in executed_from:
    #     arcpy.AddMessage("Exiting...this tool does not yet run from ArcGIS Pro")
    #     sys.exit(0)

    if not ("ARCMAP" in executed_from or "ARCCATALOG" in executed_from or
            "RUNTIME" in executed_from):
        arcpy.AddMessage("Getting parameters from command line...")

        input_fc = sys.argv[1]
        output_fgdb = sys.argv[2]
        output_fc_name = sys.argv[3]
        delete_intermediate_results = sys.argv[4]
    else:
        arcpy.AddMessage("Getting parameters from GetParameterAsText...")
        input_fc = arcpy.GetParameterAsText(0)
        output_fgdb = arcpy.GetParameterAsText(1)
        output_fc_name = arcpy.GetParameterAsText(2)
        delete_intermediate_results = arcpy.GetParameterAsText(3)

    input_fc = input_fc.strip()
    output_fgdb = output_fgdb.strip()
    output_fc_name = output_fc_name.strip()
    arcpy.AddMessage("Input Feature Class:                      {0}".format(input_fc))
    arcpy.AddMessage("Output File GDB:                          {0}".format(output_fgdb))
    arcpy.AddMessage("Output Feature Class name:                {0}".format(output_fc_name))
    arcpy.AddMessage("Delete Intermediate Results:                {0}".format(delete_intermediate_results))

    main(input_fc, output_fgdb, output_fc_name, delete_intermediate_results)

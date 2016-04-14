#-------------------------------------------------------------------------------
# Name:        bldg_proc_cgal_pt_cleanup.py
# Purpose:     this is the final script for the first task in the point processing R&D project
#
# Author:      jose6641
#
# Created:     05/02/2016
# Copyright:   (c) jose6641 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Name:        bldg1pts_cgal.py
# Purpose:
#
# Author:      jose6641
#
# Created:     08/01/2016
# Copyright:   (c) jose6641 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os, sys, time

import numpy as np
import arcpy
import CGAL.CGAL_Point_set_processing_3 as cgal_pt_proc
from CGAL.CGAL_Kernel import Point_3

# plotting
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

arcpy.env.overwriteOutput = True

def numpy_cgal_point_conversion(array):
    ''' utility function for converting array of points to CGAL Point_3 '''

    points = []
    for pt in array:
        points.append( Point_3(pt[0], pt[1], pt[2]) )

    return points

def cgal_point_to_numpy_conversion(points, dt):
    ''' utility function for converting array of points to CGAL Point_3 '''

    npts = len(points)
    array = np.zeros((npts,3))
    for i,pt in enumerate(points):
        array[i,0] = pt.x()
        array[i,1] = pt.y()
        array[i,2] = pt.z()

    array.dtype = dt
    return array


def calc_similarity_dot_product(a, b, nbins):
    ''' this function will calculate the similarity of two equal sized arrays using the dot product '''

    # calculate histograms
    ha, abins = np.histogram(a['SHAPE@Z'], nbins)
    hb, abins = np.histogram(b['SHAPE@Z'], nbins)

    # calculate cosine similarity using dot product of the pmf arrays (divided by sum)
    sim =  1-np.dot(ha/float(ha.sum()), hb/float(hb.sum()))

    return sim


def iter_remove_outliers(cgal_pts, nbins, dt, nbrs=24, pct=1.5, sim=0.99, dif = 0.01):
    ''' this applies CGAL's remove outliers iteratively similiarity measure is reached'''

    # calculate first outlier removal
    first = cgal_pt_proc.remove_outliers(cgal_pts, nbrs, pct)
    new_pts = cgal_pts[0:first]
    first_dif = len(cgal_pts) - len(new_pts)

    if not first_dif:
        return 1, new_pts
    else:
        ctr = 0
        flag = True
        sim_list = [0.0]
        while flag:
            second = cgal_pt_proc.remove_outliers(new_pts, nbrs, pct)
            #a = len(new_pts)
            a = cgal_point_to_numpy_conversion(new_pts, dt)
            new_pts = new_pts[0:second]
            #b = len(new_pts)
            b = cgal_point_to_numpy_conversion(new_pts, dt)
            sim_it = calc_similarity_dot_product(a, b, nbins)
            sim_list.append(sim_it)
            ctr += 1

            print ('{}, {}, {} points'.format(ctr, sim_it, len(b)))
            arcpy.AddMessage('{}, {}, {} points'.format(ctr, sim_it, len(b)))

            # check for similarity threshold hit or no significant change between iterations
            if (sim_it > sim) or (abs(sim_list[-1] - sim_list[-2]) < dif):
                flag = False


        return ctr+1, new_pts


def plot_points(array):
    # show the points
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(array['SHAPE@X'], array['SHAPE@Y'], array['SHAPE@Z'])
    plt.show()
    plt.close()

## main stuff
# inputs for testing
#sourceGDB = r"V:\ResearchGroup\3dcities\research_project\data\DevTesting.gdb"
#pt_fcs = ("points1", "points4", "points5_singlePart", "points6_singlePart")
#out_gdb = r"V:\ResearchGroup\3dcities\research_project\data\test_pro.gdb"

# get inputs
fc = arcpy.GetParameterAsText(0)
out_gdb = arcpy.GetParameterAsText(1)

# get feature class name
pfc = os.path.basename(fc)


if os.path.exists(out_gdb):
    arcpy.Delete_management(out_gdb)

# check for existence of out_gdb
if not os.path.exists(out_gdb):
    arcpy.CreateFileGDB_management(os.path.dirname(out_gdb), os.path.basename(out_gdb))


# process iteratively using z-histogram dot product techinique #
# iterate through point feature classes and process            #
# convert to point FC to 3-tuple array and save dtype          #

# variables needed for conversion to and from numpy world
field_names = ["SHAPE@X", "SHAPE@Y", "SHAPE@Z"]
sr = arcpy.Describe(fc).SpatialReference
array = arcpy.da.FeatureClassToNumPyArray(fc, field_names)
array_dt = array.dtype

# number of height bins (0.5 meter resolution)
nbins = np.floor((array['SHAPE@Z'].max() - array['SHAPE@Z'].min()) * 2)

# convert to CGAL points
cgal_pts = numpy_cgal_point_conversion(array)

# get point spacing using 24 nearest neighbors
avg_space = cgal_pt_proc.compute_average_spacing(cgal_pts, 24)
arcpy.AddMessage("{} points. average point spacing is: {}".format(len(cgal_pts), avg_space))

# run simplification
grid_simp = True
if grid_simp:
    grid_size = avg_space/5
    new_pts = cgal_pt_proc.grid_simplify_point_set(cgal_pts, grid_size)
    cgal_pts = cgal_pts[0:new_pts]

    print ('grid_simplification: grid size {}, points remaining {}'.format(grid_size, new_pts))
    arcpy.AddMessage('grid_simplification: grid size {}, points remaining {}'.format(grid_size, new_pts))

    # convert back to numpy
    filt_array = cgal_point_to_numpy_conversion(cgal_pts, array_dt)

    # plot the points
    #plot_points(filt_array)

rand_simp = True
if rand_simp:
    new_pts = cgal_pt_proc.random_simplify_point_set(cgal_pts, 1)
    cgal_pts = cgal_pts[0:new_pts]

    print ('random_simplification: points remaining {}'.format(new_pts))
    arcpy.AddMessage('random_simplification: points remaining {}'.format(new_pts))

    # convert back to numpy
    filt_array1 = cgal_point_to_numpy_conversion(cgal_pts, array_dt)

    # plot the points
    #plot_points(filt_array1)


# iteratively remove outliers with large neighbor count and low percentage
iters,cgal_pts = iter_remove_outliers(cgal_pts, nbins, array_dt, 24, 1.5)
print ('remove_outliers: points remaining {} after {} iterations'.format(len(cgal_pts), iters))
arcpy.AddMessage('remove_outliers: points remaining {} after {} iterations'.format(len(cgal_pts), iters))

# convert back to numpy
filt_array = cgal_point_to_numpy_conversion(cgal_pts, array_dt)

# make out_fc name
ofc = "{}_iterative_outlier_removal".format(pfc)
out_fc = os.path.join(out_gdb, ofc)

# save the processed points to a feature class
arcpy.da.NumPyArrayToFeatureClass(filt_array, out_fc, field_names, sr)

print ('finished')

#!/usr/bin/env python
# Copyright 2011 Alex Zvoleff
#
# This file is part of the AccraABM agent-based model.
# 
# AccraABM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# AccraABM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# AccraABM.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact Alex Zvoleff in the Department of Geography at San Diego State 
# University with any comments or questions. See the README.txt file for 
# contact information.

"""
Wrapper to run a set of Accra ABM model runs: Reads in input parameters, then 
calls routines to initialize and run the model, and output model statistics.

NOTE: Borrows code from matplotlib, particularly for rcsetup functions.

Alex Zvoleff, azvoleff@mail.sdsu.edu
"""

import os
import sys
import getopt
import time
import pickle
import tempfile
import subprocess
import socket
import csv

import numpy as np

from PyABM.rcsetup import write_RC_file
from PyABM.file_io import write_single_band_raster

from AccraABM import rcParams
from AccraABM.initialize import generate_world
from AccraABM.modelloop import main_loop

if rcParams['model.use_psyco'] == True:
    import psyco
    psyco.full()

def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        rc_file = sys.argv[1]
        print "\nWARNING: using default rc params. Custom rc_file use is not yet implemented.\n"
    except IndexError:
        pass

    # Get machine hostname to print it in the results file and use in the 
    # run_ID_number.
    hostname = socket.gethostname()
    
    # The run_ID_number provides an ID number (built from the start time and 
    # machine name) to uniquely identify this model run.
    run_ID_number = time.strftime("%Y%m%d-%H%M%S") + '_' + hostname
    # First strip any trailing backslash from the model.resultspath value from 
    # rcparams, so that os.path.join-ing it to the scenario.name does not lead 
    # to having two backslashes in a row.
    model_results_path_root = str.strip(rcParams['model.resultspath'], "/\\")
    scenario_path = os.path.join(str(rcParams['model.resultspath']), rcParams['scenario.name'])
    results_path = os.path.join(scenario_path, run_ID_number)
    if not os.path.exists(scenario_path):
        try:
            os.mkdir(scenario_path)
        except OSError:
            raise OSError("error creating scenario directory %s"%(scenario_path))
    try:
        os.mkdir(results_path)
    except OSError:
        raise OSError("error creating results directory %s"%(results_path))
    

    if rcParams['model.reinitialize']:
        # Generate a new world (with new resampling, etc.)
        world = generate_world()
    else:
        # Load a pickled World for use in the model.
        input_data_file = rcParams['path.initialization_file']
        file = open(input_data_file, "r")
        try:
            world = pickle.load(file)
        except IOError:
            raise IOError('error loading world data from  %s'%input_data_file)

    # Run the model loop
    start_time = time.localtime()
    start_time_string = time.strftime("%m/%d/%Y %I:%M:%S %p", start_time)
    print """
************************************************************************************************
%s: started model run number %s.
************************************************************************************************
"""%(start_time_string, run_ID_number)
    time_strings = main_loop(world, results_path) # This line actually runs the model.
    end_time = time.localtime()
    end_time_string = time.strftime("%m/%d/%Y %I:%M:%S %p", end_time) 
    print """
************************************************************************************************
%s: finished model run number %s.
************************************************************************************************
"""%(end_time_string, run_ID_number)
    
    # Save the results
    print "Saving result files..."

    # Write out the world file and mask used to run the model. Update the 
    # rcparams to point to these files so they will be reused if this run is 
    # rerun.
    #DEM_data_file = os.path.join(results_path, "AccraABM_DEM.tif")
    #array, gt, prj = world.get_DEM_data()
    #write_single_band_raster(array, gt, prj, DEM_data_file)
    world_mask_data_file = os.path.join(results_path, "AccraABM_world_mask.tif")
    array, gt, prj = world.get_world_mask_data()
    write_single_band_raster(array, gt, prj, world_mask_data_file)

    lulc_data_file = os.path.join(results_path, "AccraABM_land_cover.tif")
    array, gt, prj = world.get_lulc_data()
    write_single_band_raster(array, gt, prj, lulc_data_file)

    # TODO: write a function to handle writing out a markov_matrix from the 
    # world/rcparams. Do this maybe after writing markov_matrix handling into 
    # the rcvalidation code.
    
    # Save the SHA-1 of the commit used to run the model, along with any diffs 
    # from the commit (the output of the git diff command). sys.path[0] gives 
    # the path of the currently running AccraABM code.
    git_diff_file = os.path.join(results_path, "git_diff.patch")
    commit_hash = save_git_diff(sys.path[0], git_diff_file)

    time_csv_file = os.path.join(results_path, "time.csv")
    write_time_csv(time_strings, time_csv_file)
    
    if rcParams['model.make_animations']:
        print "Plotting results..."
        Rscript_binary = rcParams['path.Rscript_binary']
        dev_null = open(os.devnull, 'w')
        try:
            subprocess.check_call([Rscript_binary, 'plot_results.R', results_path],
                    cwd=sys.path[0], stdout=dev_null, stderr=dev_null)
        except:
            print "WARNING: Error running plot_results.R."
        dev_null.close()

    # Calculate the number of seconds per month the model took to run (to 
    # simplify choosing what machine to do model runs on). This is equal to the 
    # length of time_strings divided by the timestep size (in months).
    speed = (time.mktime(end_time) - time.mktime(start_time)) / (len(time_strings['timestep']) / rcParams['model.timestep'])

    # After running model, save rcParams to a file, along with the SHA-1 of the 
    # code version used to run it, and the start and finish times of the model 
    # run. Save this file in the same folder as the model output.
    run_RC_file = os.path.join(results_path, "AccraABMrc")
    RC_file_header = """# This file contains the parameters used for a AccraABM model run.
# Model run ID:\t%s
# Start time:\t%s
# End time:\t\t%s
# Run speed:\t%.4f
# Code version:\t%s"""%(run_ID_number, start_time_string, end_time_string, 
        speed, commit_hash)
    write_RC_file(run_RC_file, RC_file_header, rcParams)

    print "\nFinished at", time.strftime("%m/%d/%Y %I:%M:%S %p") + "."

    return 0

def save_git_diff(code_path, git_diff_file):
    # First get commit hash from git show
    temp_file_fd, temp_file_path = tempfile.mkstemp()
    try:
        git_binary= rcParams['path.git_binary']
        subprocess.check_call([git_binary, 'show','--pretty=format:%H'], stdout=temp_file_fd, cwd=code_path)
    except:
        print "WARNING: Error running git. Skipping git-diff patch output."
        return "ERROR_RUNNING_GIT"
    os.close(temp_file_fd)
    temp_file = open(temp_file_path, 'r')
    commit_hash = temp_file.readline().strip('\n')
    temp_file.close()
    os.remove(temp_file_path)

    # Now write output of git diff to a file.
    try:
        out_file = open(git_diff_file, "w")
        git_binary = rcParams['path.git_binary']
        subprocess.check_call([git_binary, 'diff'], stdout=out_file, cwd=code_path)
        out_file.close()
    except IOError:
        print "WARNING: Error writing to git diff output file: %s"%(git_diff_file)
    return commit_hash

def write_time_csv(time_strings, time_csv_file):
    """
    Write a CSV file for conversion of timestep number, float, etc. to actual 
    year and month (for plotting).
    """
    out_file = open(time_csv_file, "w")
    csv_writer = csv.writer(out_file)
    col_headers = sorted(time_strings.keys())
    csv_writer.writerow(col_headers)
    columns = []
    for col_header in col_headers:
        # Subtract 1 as Python has zero indexing but the model uses 1 to denote 
        # the first timestep.
        if columns == []:
            columns = np.array((time_strings[col_header]))
        else:
            columns = np.vstack((columns, time_strings[col_header]))
    columns = np.transpose(columns)
    csv_writer.writerows(columns)
    out_file.close()

if __name__ == "__main__":
    sys.exit(main())

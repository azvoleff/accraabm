# Copyright 2011 Alex Zvoleff
#
# This file is part of the ChitwanABM agent-based model.
# 
# ChitwanABM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# ChitwanABM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with
# ChitwanABM.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact Alex Zvoleff in the Department of Geography at San Diego State 
# University with any comments or questions. See the README.txt file for 
# contact information.

"""
Contains main model loop: Contains the main loop for the model. Takes input 
parameters read from runModel.py, and passes results of model run back.

Alex Zvoleff, azvoleff@mail.sdsu.edu
"""

import os
import time
import copy

import numpy as np

from PyABM.file_io import write_single_band_raster
from PyABM.utility import TimeSteps
from ChitwanABM import rcParams

if rcParams['model.use_psyco'] == True:
    import psyco
    psyco.full()

timebounds = rcParams['model.timebounds']
timestep = rcParams['model.timestep']

model_time = TimeSteps(timebounds, timestep)

#TODO: Setup saved_data so it can store the raster results
saved_data=[]

def main_loop(world, results_path):
    """This function contains the main model loop. Passed to it is a list of 
    regions, which contains the person agents and land cover to be used in the 
    model, and the land-use parameters."""

    time_strings = {}
    # Store the date values (as timestep number (0),  float and date string) 
    # for time zero (T0) so that the initial values of the model (which are for 
    # time zero, the zeroth timestep) can be used in plotting model results.
    time_strings['timestep'] = [0]
    time_strings['time_float'] = [model_time.get_T0_date_float()]
    time_strings['time_date'] = [model_time.get_T0_date_string()]

    # Save the starting time of the model to use in printing elapsed time while 
    # it runs.
    modelrun_starttime = time.time()

    def write_results_CSV(world, results_path, timestep):
        """
        Function to periodically save model results to CSV (if this option is 
        selected in the rc file).
        """
        if rcParams['save_psn_data']:
            world.write_persons_to_csv(timestep, results_path)
        #TODO: Also print out land use

    # Write the results for timestep 0
    write_results_CSV(world, results_path, 0)

    while model_time.in_bounds():
        for region in world.iter_regions():
            # Save event, LULC, and population data for later output to CSV.
            #region.increment_age()
            pass

        #TODO: Transition land, calculate health
        world.lulc_markov_transition()

        # Write out LULC for this timestep
        output_file = os.path.join(results_path, "lulc_time_%s.tif"%model_time.get_cur_int_timestep())
        lulc, gt, prj = world.get_lulc_data()
        write_single_band_raster(lulc, gt, prj, output_file)
                
        # Print an information line to allow keeping tabs on the model while it 
        # is running.
        num_persons = region.num_persons()
        stats_string = "%s | P: %5s "%(model_time.get_cur_date_string().ljust(7), num_persons)
        print stats_string

        # Save timestep, year and month, and time_float values for use in 
        # storing results (to CSV) keyed to a particular timestep.
        time_strings['timestep'].append(model_time.get_cur_int_timestep())
        time_strings['time_float'].append(model_time.get_cur_date_float())
        time_strings['time_date'].append(model_time.get_cur_date_string())

        if num_persons == 0:
            print "End of model run: population is zero."
            break

        write_results_CSV(world, results_path, model_time.get_cur_int_timestep())

        model_time.increment()

    return time_strings

def elapsed_time(start_time):
    elapsed = int(time.time() - start_time)
    hours = elapsed / 3600
    minutes = (elapsed - hours * 3600) / 60
    seconds = elapsed - hours * 3600 - minutes * 60
    return "%ih %im %is" %(hours, minutes, seconds)

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
Sets up a AccraABM model run: Initializes person agents and land use using the 
original WHSA data.

Alex Zvoleff, azvoleff@mail.sdsu.edu
"""

import os
import sys
import copy

import numpy as np
import pickle
from subprocess import check_call, CalledProcessError

from PyABM.file_io import read_single_band_raster

from AccraABM import rcParams
from AccraABM.agents import World

def main():
    model_world = generate_world()
    processed_data_file = rcParams['path.initialization_file']

    try:
        save_world(model_world, processed_data_file)
    except:
        print "ERROR: while saving world file to %s"%(processed_data_file)

def read_WHSA_data(textfile, key_field):
    """
    Reads in WHSA data from a CSV file into a dictionary of dictionary objects, 
    where the first line of the file gives the column headings (used as keys 
    within the nested dictionary object). No conversion of the fields is done: 
    they are all stored as strings, EXCEPT for the key_field, which is 
    converted and stored as an int.
    """

    try:
        file = open(textfile, 'r')
        lines = file.readlines()
        file.close()
    except:
        raise IOError("error reading %s"%(textfile))
    
    # The first line of the data file gives the column names
    col_names = lines[0].split(',')
    for n in xrange(len(col_names)):
        col_names[n] = col_names[n].strip('\n \"')

    data = {}
    for line in lines[1:]:
        fields = line.split(',')
        for n in xrange(len(fields)):
            fields[n] = fields[n].strip('\n \"')

        new_data = {} 
        for field, column_name in zip(fields, col_names):
            new_data[column_name] = field

        data_key = int(new_data[key_field])
        if new_data.has_key(data_key):
            raise KeyError('key %s is already in use'%(data_key))
        data[data_key] = new_data

    return data

def assemble_persons(populationFile, model_world):
    """
    Reads data in from the WHSA, which was read in by the data_preprocess.R R 
    script. This function then assembles person agents from the data.
    """
    people = read_WHSA_data(populationFile, "id")

    persons = []
    for person_id in people.iterkeys():
        person = people[person_id]
        id = person['id']
        x = float(person['x_utm30'])
        y = float(person['y_utm30'])
        # Note that age is converted to months.
        age_months = int(person['age']) * 12
        sex = 'female'
        health = int(person['rand_sf36_pf'])
        education = int(person['education'])
        charcoal = bool(person['Charcoal'])
        yrs_in_house_cat = str(person['YrsInHouse_Cat'])
        own_toilet = bool(person['OwnToilet'])
        ethnicity = int(person['major_ethnic'])
        initial_agent = True
        person = model_world.new_person(None, id, age_months,
                sex, initial_agent, ethnicity, education, charcoal, own_toilet, yrs_in_house_cat, x, y, 
                health)
        persons.append(person)
        
    return persons

def assemble_world():
    """
    Puts together a single world (with, currently, only a single region) from 
    the WHSA data using the above functions to input WHSA data on sample 
    population.
    """
    model_world = World()

    raw_data_path = rcParams['path.raw_input_data']
    population_data_path = os.path.join(raw_data_path, rcParams['inputfile.sample_population'])

    persons = assemble_persons(population_data_path, model_world)

    # Add the DEM, VIS, NDVI and Study Area masks to the model_world instance.
    #DEM_file = os.path.join(raw_data_path, rcParams['path.DEM_file'])
    #model_world.set_DEM_data(read_single_band_raster(DEM_file))

    land_cover_file = os.path.join(raw_data_path, rcParams['inputfile.land_cover'])
    lulc, gt, prj = read_single_band_raster(land_cover_file)
    lulc[lulc < -1] = rcParams['lulc.NA_value'] # Assume anything coded <-1 is NaN
    model_world.set_lulc_data(lulc, gt, prj)

    world_mask_file = os.path.join(raw_data_path, rcParams['inputfile.world_mask'])
    world_mask, gt, prj = read_single_band_raster(land_cover_file)
    model_world.set_world_mask_data(world_mask, gt, prj)

    markov_matrix_path = os.path.join(raw_data_path, rcParams['inputfile.markov_matrix'])
    model_world.set_lulc_markov_matrix(markov_matrix_path)

    # TODO: Use the world mask to produce a square cutout of the data, 
    # eliminating as much of the NA area as possible.

    # Populate the Accra region.
    region = model_world.new_region()
    for person in persons:
        region.add_agent(person)

    print "\nPersons: %s, Region shape: %s"%(region.num_persons(), np.shape(model_world.get_lulc()))
    
    # Save a list of the original agents used to initialize the region. This 
    # list will be in re-sampling as the model runs, to repopulate the region 
    # as agents die.
    region._original_agents = []
    for person in persons:
        region._original_agents.append(copy.copy(person))
    return model_world

def save_world(world, filename):
    "Pickles a world for later reloading."
    file = open(filename, "w")
    pickle.dump(world, file)

def generate_world():
    """
    Performs the complete process necessary for initializing the model from    
    WHSA restricted data.
        1) Calls the necessary R script for preparing the necessary CSV 
        initialization files from the WHSA data. 

        2) Calls the assemble_world function to prepare an instance of the 
        World class to be used in the model.

        3) Saves this world instance in the standard location. NOTE: This must 
        be an encrypted directory that is not publically accessible to conform 
        to ICPSR and IRB requirements.
    """
    print
    if rcParams['model.reinitialize.run_R_script']:
        try:
            Rscript_binary = rcParams['path.Rscript_binary']
            print "Calling R to preprocess WHSA data..."
            raw_data_path = rcParams['path.raw_input_data']
            imagery_path = rcParams['inputfile.imagery']
            WHSA_path = rcParams['inputfile.WHSA']
            Accra_FMV_path = rcParams['inputpath.Accra_FMVs']
            buffer = str(rcParams['lulc.buffer'])
            windowed_Markov = str(rcParams['model.Markov_window'])
            Markov_window_size = str(rcParams['model.Markov_window_size'])
            check_call([Rscript_binary, "data_preprocess.R", raw_data_path, 
                imagery_path, WHSA_path, Accra_FMV_path, buffer, windowed_Markov,
                Markov_window_size])
        except CalledProcessError:
            print "ERROR: while running data_preprocess.R script"
    else:
        print "Skipping call of data_preprocess.R..."
    print "Generating world from preprocessed WHSA data..."
    model_world = assemble_world()

    return model_world

if __name__ == "__main__":
    sys.exit(main())

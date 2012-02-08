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
Contains the classes for Person and EA agents.  Person agents are subclasses of 
the Agent class, while EA and is a subclass of the Agent_set object.

Alex Zvoleff, azvoleff@mail.sdsu.edu
"""

import os
import csv

import numpy as np

from PyABM import IDGenerator, boolean_choice
from PyABM.agents import Agent, Agent_set, Agent_Store

from AccraABM import rcParams, random_state
#from AccraABM.statistics import *

if rcParams['model.use_psyco'] == True:
    import psyco
    psyco.full()

class Person(Agent):
    "Represents a single person agent"
    def __init__(self, world, birthdate, PID=None, age=0, sex=None, 
            initial_agent=False, ethnicity=None, religion=None, education=None, 
            EA=None, hweight08=None, ses=None):
        Agent.__init__(self, world, PID, initial_agent)

        # birthdate is the timestep of the birth of the agent. It is used to 
        # calculate the age of the agent. Agents have a birthdate of 0 if they 
        # were BORN in the first timestep of the model.  If they were used to 
        # initialize the model their birthdates will be negative.
        self._birthdate = birthdate

        # self._age is used as a convenience to avoid the need to calculate the 
        # agent's age from self._birthdate each time it is needed. It is         
        # important to remember though that all agent's ages must be 
        # incremented with each model timestep, and are expressed in months.
        # The age starts at 0 (it is zero for the entire first timestep of the 
        # model).
        self._age = age

        # deathdate is used for tracking agent deaths in the results, mainly 
        # for debugging.
        self._deathdate = None
        self._alive = True

        if sex==None:
            # Person agents are randomly assigned a sex if they were not given 
            # one.
            if boolean_choice():
                self._sex = 'female'
            else:
                self._sex = 'male'
        elif sex in ['male', 'female']:
            self._sex = sex
        else:
            raise ValueError("%s is not a valid gender"%(sex))

        # self._initial_agent is set to "True" for agents that were used to 
        # initialize the model.
        self._initial_agent = initial_agent

        self._ethnicity = ethnicity
        self._education = education
        self._religion = religion
        self._EA = EA
        self._hweight08 = hweight08
        self._ses = ses
        self._health = None

    def get_sex(self):
        return self._sex

    def get_age(self):
        return self._age

    def get_ethnicity(self):
        return self._ethnicity

    def is_initial_agent(self):
        return self._initial_agent

    def kill(self, time):
        self._alive = False
        self._deathdate = time

    def __str__(self):
        return "Person(PID: %s. EA: %s.)" %(self.get_ID(), self.get_parent_agent().get_ID())

class Region(Agent_set):
    """Represents a set of agents sharing a spatial area (and therefore 
    land use data), and demographic characteristics."""
    def __init__(self, world, ID=None, initial_agent=False):
        Agent_set.__init__(self, world, ID, initial_agent)

    def __repr__(self):
        use
        #TODO: Finish this
        return "__repr__ UNDEFINED"

    def __str__(self):
        return "Region(RID: %s, %s person(s))"%(self.get_ID(), 
                self.num_persons())

    def is_initial_agent(self):
        return self._initial_agent

    def iter_persons(self):
        "Returns an iterator over all the persons in the region"
        for person in self.iter_agents():
            yield person

    def increment_age(self):
        """Adds one to the age of each agent. The units of age are dependent on 
        the units of the input rc parameters."""
        for person in self.iter_persons():
            timestep = rcParams['model.timestep']
            person._age += timestep

    def num_persons(self):
        "Returns the number of persons in the population."
        return len(self._members)

class World():
    """The world class generates new agents, while tracking ID numbers to 
    ensure that they are always unique across each agent type. It also contains 
    a dictionary with all the regions in the model."""
    def __init__(self):
        # _members stores member regions in a dictionary keyed by RID
        self._members = {}

        # These IDGenerator instances generate unique ID numbers that are never 
        # reused, and always unique (once used an ID number cannot be 
        # reassigned to another agent). All instances of the Person class, for  
        # example, will have a unique ID number generated by the PIDGen 
        # IDGenerator instance.
        self._PIDGen = IDGenerator()
        self._HIDGen = IDGenerator()
        self._NIDGen = IDGenerator()
        self._RIDGen = IDGenerator()

    def set_DEM_data(world_mask, prj, gt):
        self._DEM_array = DEM
        self._DEM_prj = prj
        self._DEM_gt = gt
        return 0

    def get_DEM(self):
        return self._DEM_array

    def get_DEM_data(self):
        return self._DEM_array, self._DEM_prj, self._DEM_gt

    def set_world_mask_data(world_mask, prj, gt):
        self._world_mask_array = world_mask
        self._world_mask_prj = prj
        self._world_mask_gt = gt
        return 0

    def get_world_mask(self):
        return self._world_mask_array

    def get_world_mask_data(self):
        return self._world_mask_array, self._world_mask_prj, self._world_mask_gt

    def new_person(self, birthdate, PID=None, age=0, sex=None, 
            initial_agent=False, ethnicity=None, religion=None, education=None, 
            EA=None, hweight08=None, ses=None):
        "Returns a new person agent."
        if PID == None:
            PID = self._PIDGen.next()
        else:
            # Update the generator so the PID will not be reused
            self._PIDGen.use_ID(PID)
        return Person(self, birthdate, PID, age, sex, initial_agent, ethnicity, 
                religion, education, EA, hweight08, ses)

    def new_region(self, RID=None, initial_agent=False):
        "Returns a new region agent, and adds it to the world member list."
        if RID == None:
            RID = self._RIDGen.next()
        else:
            # Update the generator so the RID will not be reused
            self._RIDGen.use_ID(RID)
        region = Region(self, RID, initial_agent)
        self._members[region.get_ID()] = region
        return region

    def get_regions(self):
        return self._members.values()

    def iter_regions(self):
        "Convenience function for iteration over all regions in the world."
        for region in self._members.values():
            yield region

    def iter_persons(self):
        "Convenience function used for things like incrementing agent ages."
        for region in self.iter_regions():
            for person in region.iter_persons():
                yield person

    def write_persons_to_csv(self, timestep, results_path):
        """
        Writes a list of persons, with a header row, to CSV.
        """
        psn_csv_file = os.path.join(results_path, "psns_time_%s.csv"%timestep)
        out_file = open(psn_csv_file, "w")
        csv_writer = csv.writer(out_file)
        csv_writer.writerow(["pid", "rid", "gender", "age", "religion", "education", "ses", "health", "ethnicity"])
        for region in self.iter_regions():
            for person in region.iter_persons():
                new_row = []
                new_row.append(person.get_ID())
                new_row.append(person.get_parent_agent().get_ID())
                new_row.append(person.get_sex())
                new_row.append(person.get_age())
                new_row.append(person._religion)
                new_row.append(person._education)
                new_row.append(person._ses)
                new_row.append(person._health)
                new_row.append(person._ethnicity)
                csv_writer.writerow(new_row)
        out_file.close()

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
Contains statistical models to calulate probabilities (such as of death), and 
to run markov model (to calculate land use).
"""

import numpy as np

from ChitwanABM import rcParams

if rcParams['model.use_psyco'] == True:
    import psyco
    psyco.full()

probability_time_units = rcParams['probability.time_units']

class UnitsError(Exception):
    pass

def convert_probability_units(probability):
    """
    Converts probability so units match timestep used in the model, assuming probability 
    function is uniform across the interval.

    Conversions are made accordingly using conditional probability.
    """
    # If the probability time units don't match the model timestep units, then the 
    # probabilities need to be converted.
    if probability_time_units == 'months':
        pass
    elif probability_time_units == 'years':
        for key, value in probability.iteritems():
            probability[key] = 1 - (1 - value)**(1/12.)
    elif probability_time_units == 'decades':
        for key, value in probability.iteritems():
            probability[key] = 1 - (1 - value)**(1/120.)
    else:
        raise UnitsError("unhandled probability_time_units")
    return probability

#TODO: these probabilities should be derived from the region, not directly from rcParams
death_probabilities_male = convert_probability_units(rcParams['probability.death.male'])
death_probabilities_female = convert_probability_units(rcParams['probability.death.female'])

def __probability_index__(t):
    """
    Matches units of time in model to those the probability is expressed in. For 
    instance: if probabilities are specified for decades, whereas the model runs in 
    months, __probability_index__, when provided with an age in months, will convert 
    it to decades, rounding down. NOTE: all probabilities must be expressed with the 
    same time units.
    """
    if probability_time_units == 'months':
        return t
    elif probability_time_units == 'years':
        return int(round(t / 12.))
    elif probability_time_units == 'decades':
        return int(round(t / 120.))
    else:
        raise UnitsError("unhandled probability_time_units")

def calc_probability_death(person):
    "Calculates the probability of death for an agent."
    age = person.get_age()
    probability_index = __probability_index__(age)
    try:
        if person.get_sex() == 'female':
            return death_probabilities_female[probability_index]
        elif person.get_sex() == 'male':
            return death_probabilities_male[probability_index]
    except IndexError:
        raise IndexError("error calculating death probability (index %s)"%(probability_index))

def draw_from_prob_dist(prob_dist):
    """
    Draws a random number from a manually specified probability distribution,
    where the probability distribution is a tuple specified as:
        ([a, b, c, d], [1, 2, 3])
    where a, b, c, and d are bin limits, and 1, 2, and 3 are the probabilities 
    assigned to each bin. Notice one more bin limit must be specifed than the 
    number of probabilities given (to close the interval).
    """
    # First randomly choose the bin, with the bins chosen according to their 
    # probability.
    binlims, probs = prob_dist
    num = np.random.rand() * np.sum(probs)
    n = 0
    probcumsums = np.cumsum(probs)
    for upprob in probcumsums[1:]:
        if num < upprob:
            break
        n += 1
    upbinlim = binlims[n+1]
    lowbinlim = binlims[n]
    # Now we know the bin lims, so draw a random number evenly distributed 
    # between those two limits.
    return np.random.uniform(lowbinlim, upbinlim)

def markov_transition(region):
    # TODO: Right now this only works if there are only two classes. Rewrite 
    # this to work for a wider range of scenarios.
    current_lu = self._land_cover_array
    new_lu = current_lu
    random_values = np.random.random(np.shape(self._land_cover_array))
    for t1 in trans_matrix.keys():
        for t2 in trans_matrix[t1].keys():
            lower_prob_bound, upper_prob_bound = trans_matrix[t1][t2]
            new_lu[(current_lu == t1) & (random_values >= lower_prob_bound) & (random_values < upper_prob_bound)] <- t2
    self._land_cover_array = new_lu
    return self._land_cover_array

    

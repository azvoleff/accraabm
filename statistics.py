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
Contains statistical models to calculate probabilities (such as of death), and 
to run Markov model (to calculate land use).
"""

import numpy as np

from AccraABM import rcParams

if rcParams['model.use_psyco'] == True:
    import psyco
    psyco.full()

probability_time_units = rcParams['probability.time_units']

class UnitsError(Exception):
    pass

death_probabilities_female = rcParams['probability.death.female']

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
    age = person.get_age_months()
    probability_index = __probability_index__(age)
    try:
        return death_probabilities_female[probability_index]
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

def calculate_cover_fraction_NBH(person_ids, egocentric_nbhs, value, na_value):
    # note that areas are expressed in pixels.
    # in below line, np.invert is use for bitwise not (to select values that 
    # are not nan
    total_area = np.sum(np.sum(egocentric_nbhs != na_value, 1), 0)
    # convert total area to float or else it will be an integer/integer for the 
    # later division
    total_area = np.array(total_area, dtype="float")
    cover_area = np.sum(np.sum((egocentric_nbhs == value) & 
        (egocentric_nbhs != na_value), 1), 0)
    cover_fractions = (cover_area / total_area)
    cover_fractions_dict = {}
    for cover_fraction, person_id in zip(cover_fractions, person_ids):
        cover_fractions_dict[person_id] = cover_fraction
    return cover_fractions_dict

def calculate_cover_fraction_world(lulc, value, na_value):
    # note that areas are expressed in pixels.
    # in below line, np.invert is use for bitwise not (to select values that 
    # are not nan
    total_area = np.sum(np.sum(lulc != na_value, 1), 0)
    # convert total area to float or else it will be an integer/integer for the 
    # later division
    total_area = float(total_area)
    cover_area = np.sum(np.sum((lulc == value) & (lulc != na_value), 1), 0)
    cover_fraction = (cover_area / total_area)
    return cover_fraction

def predict_self_reported_health(person):
    """
    Calculates the self_reported_health of an agent, using the results of an 
    ordinal logistic regression (AKA proportional odds model). See Harell 
    (2001, 333) for formula.
    """
    levels = rcParams['srh.olr.depvar_levels']
    prob_y_gte_j = np.zeros(len(levels) - 1) # probability y >= j
    for n in np.arange(len(prob_y_gte_j)):
        intercept = rcParams['srh.olr.intercepts'][n]
        xb_sum = 0
        # Individual-level characteristics
        xb_sum += rcParams['srh.olr.coef.age'] * person.get_age_years()
        if person.get_ethnicity() == 1:
            xb_sum += rcParams['srh.olr.coef.major_ethnic_1']
        elif person.get_ethnicity() == 2:
            xb_sum += rcParams['srh.olr.coef.major_ethnic_2']
        elif person.get_ethnicity() == 3:
            xb_sum += rcParams['srh.olr.coef.major_ethnic_3']
        elif person.get_ethnicity() == 4:
            xb_sum += rcParams['srh.olr.coef.major_ethnic_4']
        else:
            raise StatisticsError("Ethnicity %s does not have a specified coefficient"%person.get_ethnicity())
        xb_sum += rcParams['srh.olr.coef.education'] * person._education
        # Neighborhood chararacteristics
        xb_sum += rcParams['srh.olr.coef.veg_fraction'] * person._veg_fraction
        prob_y_gte_j[n] = 1. / (1 + np.exp(-(intercept + xb_sum)))
    prob_y_eq_j = np.zeros(4) # probability y == j
    prob_y_eq_j[0] = 1 - prob_y_gte_j[0]
    # Loop over all but the first cell of prob_y_eq_j
    for j in np.arange(1, len(prob_y_gte_j)):
        prob_y_lt_j = np.sum(prob_y_eq_j[0:j])
        prob_y_eq_j[j] = 1 - prob_y_gte_j[j] - prob_y_lt_j
    prob_cutoffs = np.cumsum(prob_y_eq_j)
    prob_cutoffs[-1]  = 1
    rand = np.random.rand()
    for n in np.arange(len(prob_cutoffs)):
        if rand <= prob_cutoffs[n]:
            return levels[n]
    raise StatisticsError("Check level calculation - no class predicted")

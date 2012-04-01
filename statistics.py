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

def predict_physical_functioning(person):
    """
    Calculates the physical functioning score of an agent, using the results of 
    a spatial autoregressive lag model.
    """
    if rcParams['lulc.use_egocentric']:
        prefix = 'reg.pf.Ego'
    else:
        prefix = 'reg.pf.FMV'

    ##################################
    # Neighborhood characteristics
    ##################################
    # Note that in the regression model the coefficient for veg fraction is 
    # calculated on the log of the (percentage veg in NBH + 1)
    xb_sum = rcParams[prefix + '.coef.log_veg_fraction'] * np.log(person._veg_fraction*100 + 1)

    ##################################
    # Individual-level characteristics
    ##################################
    xb_sum += rcParams[prefix + '.coef.age'] * person.get_age_years()

    xb_sum += rcParams[prefix + '.coef.age_squared'] * (person.get_age_years()**2)

    if person.get_ethnicity() == 1:
        xb_sum += rcParams[prefix + '.coef.major_ethnic_1']
    elif person.get_ethnicity() == 2:
        xb_sum += rcParams[prefix + '.coef.major_ethnic_2']
    elif person.get_ethnicity() == 3:
        xb_sum += rcParams[prefix + '.coef.major_ethnic_3']
    elif person.get_ethnicity() == 4:
        xb_sum += rcParams[prefix + '.coef.major_ethnic_4']
    else:
        raise StatisticsError("Ethnicity %s does not have a specified coefficient"%person.get_ethnicity())

    if person._education == 0:
        xb_sum += rcParams[prefix + '.coef.education_0']
    elif person._education == 1:
        xb_sum += rcParams[prefix + '.coef.education_1']
    elif person._education == 2:
        xb_sum += rcParams[prefix + '.coef.education_2']
    elif person._education == 3:
        xb_sum += rcParams[prefix + '.coef.education_3']
    elif person._education == 4:
        xb_sum += rcParams[prefix + '.coef.education_4']
    else:
        raise StatisticsError("Education %s does not have a specified coefficient"%person._education)

    xb_sum += rcParams[prefix + '.coef.charcoal'] * person._charcoal

    xb_sum += rcParams[prefix + '.coef.own_toilet'] * person._own_toilet

    if person._yrs_in_house_cat == "5":
        xb_sum += rcParams[prefix + '.coef.yrs_in_house_cat_0']
    elif person._yrs_in_house_cat == "15":
        xb_sum += rcParams[prefix + '.coef.yrs_in_house_cat_5']
    elif person._yrs_in_house_cat == "30":
        xb_sum += rcParams[prefix + '.coef.yrs_in_house_cat_15']
    elif person._yrs_in_house_cat == ">30":
        xb_sum += rcParams[prefix + '.coef.yrs_in_house_cat_gt30']
    else:
        raise StatisticsError("Yrs_In_House_Cat %s does not have a specified coefficient"%person._yrs_in_house_cat)

    # Intercept
    xb_sum += rcParams[prefix + '.coef.intercept']

    return xb_sum

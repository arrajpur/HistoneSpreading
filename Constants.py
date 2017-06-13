## Constants.py
## Author: Aparna Rajpurkar

# imports
import math
import random
import numpy as np
import scipy.stats as stats
from MyEnum import States

# Set constants
TIMESTEPS_PER_CELLCYCLE = 0
# rates -- be sure to change these
# per nuc per cell cycle
CR_U_to_A = 1
CR_A_to_U = 1
CR_U_to_M = 1
CR_M_to_U = 1

POWER = 2 # arbitrary number

# set random seeds for reproducibility
SEED = 1
random.seed(SEED)
np.random.seed(SEED)

# set colors
GRAY = (0.662745,0.662745,0.662745)
RED = (0.545098,0,0)
BLUE = (0,0,0.545098)

# begin functions
def get_max_events():
    '''
    get_max_events()
    check which rate yields the maximum number of events
    return that number
    '''
    return max(CR_U_to_A, CR_A_to_U, CR_U_to_M, CR_M_to_U) / TIMESTEPS_PER_CELLCYCLE


def get_rate(old, new):
    ''' 
    given an old and new nucleosome state, return the rate of 
    conversion from old to new
    '''
    
    # currently, not using rates because could not find biological ratios
    # 0 means event happens instantaneously
    return 0

    # the below math may be wrong -- check. Dependent on what kind of rate is given
    if old == new:
        return 0
    else:
        if old == States.M_STATE:
            return 1 / CR_M_to_U * TIMESTEPS_PER_CELLCYCLE

        elif old == States.A_STATE:
            return 1 / CR_A_to_U * TIMESTEPS_PER_CELLCYCLE

        elif old == States.U_STATE:
            if new == States.M_STATE:
                return 1 / CR_U_to_M * TIMESTEPS_PER_CELLCYCLE

            if new == States.A_STATE:
                return 1 / CR_U_to_A * TIMESTEPS_PER_CELLCYCLE

def truncated_power_law(power, limit_neg, limit_pos):
    '''
    truncated_power_law(power_constant, negative_limit, positive_limit)
    use to calculate random probability based on a double truncated power law
    not currently used
    '''
    # pick a random number
    prob_direction = random.random()

    # get the probability of going negative
    neg_prob = limit_neg / (limit_neg + limit_pos)

    if prob_direction < neg_prob:
        # negative direction
        # limit_neg must be a positive number of nucs in neg direction
        x = np.arange(1, limit_neg+1, dtype='float')
        pmf = 1/x**power
        pmf /= pmf.sum()
        dist = stats.rv_discrete(values=(range(1, limit_neg+1), pmf))
        return 0 - dist.rvs(size=1)
    else:
        # positive direction
        x = np.arange(1, limit_pos+1, dtype='float')
        pmf = 1/x**power
        pmf /= pmf.sum()
        dist = stats.rv_discrete(values=(range(1, limit_pos+1), pmf))
        return dist.rvs(size=1)

def state_to_color(state):
    '''
    state_to_color(STATE)
    convert a MyEnum State class constant to a color
    '''
    if state == States.A_STATE:
        return RED

    if state == States.U_STATE:
        return GRAY

    if state == States.M_STATE:
        return BLUE

def exponential(mean):
    '''
    exponential(mean)
    calculate a random number from an exponential distribution
    not currently used
    '''
    return - mean * math.log(random.random())


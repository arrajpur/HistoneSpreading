## Chromatin.py
## Author: Aparna Rajpurkar

# imports
import sys
import math
import random
import numpy as np
import operator
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.stats as stats
import Constants
from MyEnum import ProbSpread, States, Domain, DomainBleed, Divisions, ProbConv, Recruit

## begin class definitions ##

class Nucleosome:
    '''
    Nucleosome class
    store information about a single nucleosome
    and perform functions on that nucleosome
    '''

    def __init__(self, init_state, left, right):
        ''' initialization function '''
        # check if we initialize the state randomly
        if init_state == States.INIT_STATE:
            self.state = int(random.choice(
                [States.U_STATE, States.M_STATE, States.A_STATE]
                ))
        else:
            # else set all states to the init state
            self.state = init_state

        # domain limits
        self.left_limit = left
        self.right_limit = right

    def change_state(self, new):
        ''' 
        change state 
        this was used with the old version where states were not MyEnums
        '''

        if new > self.state:
            self.state += 1
        elif new < self.state:
            self.state -= 1

class Chromatin:
    '''
    Chromatin class
    Major class of the simulation. This class handles all chromatin dynamics
    across the string of nucleosomes. Modify with caution.
    '''

    def __init__(self, input_dat):
        '''
        initialization function
        '''
        # set the input data as a class variable
        self.dat = input_dat

        # initialize class variables
        self.events = []
        self.totals = {
                States.M_STATE:0, 
                States.A_STATE:0, 
                States.U_STATE:0
                }

        self.colors = []
        self.nucleosomes = []

        self.TIME = 0
        self.time_array = []

        self.prob_mat = np.zeros(shape=(input_dat['n'],input_dat['n']))
        self.M_mat = np.zeros(shape=(input_dat['n']))
        self.A_mat = np.zeros(shape=(input_dat['n']))

        # run initialization functions
        self.init_nucs(input_dat['adv']['domain'], input_dat['n'], input_dat['i'], input_dat['data']['domains'])

        self.init_prob_mat(input_dat['n'], input_dat['adv']['domainbleed'], input_dat['data']['domainbleed'])

        self.init_colors_and_state_mats(input_dat['n'])


    ## init functions ##
    def init_colors_and_state_mats(self, n):
        '''
        init_colors_and_state_mats()
        initialize colors and state matricies for entire string of 
        nucleosomes
        '''
        # iterate over all nucleosomes
        for i in range(n):
            # handle colors and state totals
            curr_state = self.nucleosomes[i].state
            self.colors.append(Constants.state_to_color(curr_state))
            self.totals[curr_state] += 1

            # handle boolean arrays of states
            if curr_state == States.M_STATE:
                self.M_mat[i] = 1
            if curr_state == States.A_STATE:
                self.A_mat[i] = 1

    def init_prob_mat(self, n_nucs, db_enum, db_val):
        '''
        init_prob_mat()
        initialize probability matrix based on input options
        '''

        # iterate over all nucs
        for col in range(n_nucs):
            # get the limit of each nuc
            left = self.nucleosomes[col].left_limit
            right = self.nucleosomes[col].right_limit
            
            # for every other nuc:
            for row in range(n_nucs):
                # if we are on the same index, skip
                if col == row:
                    continue
                
                # otherwise, calculate probability
                prob = self.calc_prob(col, row, right, left, db_enum, db_val, n_nucs)

                if prob == -1:
                    continue
                elif prob == -2:
                    break

                # add calculated probability to correct cell of matrix
                self.prob_mat[col,row] = prob

    def init_nucs(self, domain_enum, n, initstate, num_domains):
        '''
        init_nucs()
        initialize string of nucleosomes and set limits
        '''

        # check domain option
        if domain_enum == Domain.NONE:
            # if no domains, just make the limits the start and end of the string
            self.nucleosomes = [ 
                    Nucleosome(initstate, 0, n) for x in range(n) 
                    ]
        elif domain_enum == Domain.EQUAL_DEFAULT:
            # else, calculate limits
            # NOTE: domain bleedthrough is NOT currently considered
            domain_sizes = math.floor(n / num_domains)
            num_d = num_domains

            for i in range(n):
                # calculate indicies for each nucleosome
                low_index = math.floor(i / domain_sizes) * domain_sizes
                high_index = low_index + domain_sizes - 1

                if high_index >= n :
                    high_index = n - 1

                self.nucleosomes.append(
                        Nucleosome(initstate, low_index, high_index)
                        )


    ## helper functions ##
    def calc_prob(self, col, row, right, left, domainbleed_enum, domainbleed_val, n):
        '''
        calc_prob()
        calculate probability given domain options and indicies for a nucleosome
        '''
        # if no domain bleed
        if domainbleed_enum == DomainBleed.NONE:
            # hard limits on probability
            if row < left:
                return -1 # continue
            elif row > right:
                return -2 # break
            
            if row < col:
                return stats.powerlaw.ppf( (col - row) / (col - left), Constants.POWER)
            else:
                return stats.powerlaw.ppf( (row - col) / (right - col), Constants.POWER)
            
        else:
            # calculate probabilities if domain bleedthrough is possible
            extreme_left = left
            extreme_right = right

            if left > 0:
                extreme_left = self.nucleosomes[left - 1].left_limit

            if right < n - 1:
                extreme_right = self.nucleosomes[right + 1].right_limit

            if row < extreme_left:
                return -1 # continue

            if row > extreme_right:
                return -2 # break

            prob = 0

            if row < left:
                prob = domainbleed_val * stats.powerlaw.ppf( 
                        (col - row) / (col - extreme_left), Constants.POWER
                        )
            elif row > right:
                prob = domainbleed_val * stats.powerlaw.ppf(
                        (row - col) / (extreme_right - col), Constants.POWER
                        )
            elif row < col:
                prob = stats.powerlaw.ppf( 
                        (col - row) / (col - left), Constants.POWER
                        )
            elif row > col:
                prob = stats.powerlaw.ppf( 
                        (row - col) / (right - col), Constants.POWER
                        )

            return prob

    def handle_timers(self, index, old, new, timers, nuc_seq, map_seq, lim):
        '''
        handle_timers()
        calculate t_next, add new timer if greater than this timestep, else update
        '''
        # calculate t_next from an exponential distribution based on the 
        # rate of conversion
        t_next = int(np.random.exponential(Constants.get_rate(old, new)))

        if t_next > 0:
            # add new timer
            timers[index] = {
                'timer' : t_next,
                'old' : old,
                'new' : new
                    }
            return self.fake_del(map_seq, nuc_seq, lim, index)

        else:
            # update if t_next fals in this timespan
            self.update(old, new, index)
            return lim

    def fake_del(self, map_seq, nuc_seq, lim, del_elem):
        ''' 
        fake_del()
        lazy deletion function
        '''

        if lim >= 0:
            # perform lazy deletion by rearranging indicies and returning
            # a pseudo-limit on the array
            lim -= 1

            tmp = nuc_seq[lim]

            del_elem_index = map_seq[del_elem]
            nuc_seq[lim] = nuc_seq[del_elem_index]
            map_seq[del_elem] = lim

            nuc_seq[del_elem_index] = tmp
            map_seq[tmp] = del_elem_index

            return lim
        else:
            print("ERROR trying to delete but everything is deleted!")

    def fake_readd(self, map_seq, nuc_seq, lim, add_elem):
        '''
        fake_readd()
        lazy re-addition function to pair with fake_del()
        '''

        if lim < len(nuc_seq):
            # re-add an element by rearranging indicies and updating a fake limit
            add_elem_index = map_seq[add_elem]
            tmp = nuc_seq[lim]

            nuc_seq[lim] = nuc_seq[add_elem_index]
            map_seq[add_elem] = lim

            nuc_seq[add_elem_index] = tmp
            map_seq[tmp] = add_elem_index

            return lim + 1
        else:
            print("ERROR trying to re-add but nothing deleted!")

    def print_nucs(self, fp):
        '''
        print_nucs()
        print out current state of all nucleosomes to a file
        '''
        for n in self.nucleosomes:
            fp.write(States.enum_to_string(n.state))

        fp.write("\n")
        
    ## Timestep simulation
    def timesim(self, n_nucs, sim_num):
        '''
        timesim()
        the major simulation function. Simulate chromatin spreading through time
        using parallel event simulation
        '''

        # calculate the number of events per timestep
        EVENTS_PER_TIMESTEP = int(Constants.get_max_events() * self.dat['n'])
        print("Events_per_timestep:", EVENTS_PER_TIMESTEP)
        TOT_TIMESTEPS = self.dat['t'] 
        prob_event = EVENTS_PER_TIMESTEP / n_nucs
    
        # initialize data structures
        timers = {}
        map_to_seq = [ x for x in range(n_nucs) ]
        nuc_index_seq = [ x for x in range(n_nucs) ]
        lim = n_nucs

        # open outfile
        fp = open(self.dat['o'] + "_" + str(sim_num) + ".txt", "w")

        # iterate over all timesteps
        for t in range(TOT_TIMESTEPS):
            # handle timers first
            self.print_nucs(fp)

            # initialize array to hold timer indicies to delete
            delete_timers = []
            # check each timer
            for nuc_index in timers:
                if timers[nuc_index]['timer'] == 0:
                    # add back to pool & update
                    lim = self.fake_readd(map_to_seq, nuc_index_seq, lim, nuc_index)

                    old = timers[nuc_index]['old']
                    new = timers[nuc_index]['new']
                    
                    self.update(old, new, nuc_index)
                    delete_timers.append(nuc_index)
                else:
                    # decrement timer's count
                    timers[nuc_index]['timer'] -= 1

            # remove all timers to be deleted
            for timer in delete_timers:
                del timers[timer]

            # handle divisions
            if self.dat['d'] != Divisions.NONE and t != 0:
                div_rate = self.dat['data']['divisions']

                # check if we are dividing now
                if t % div_rate == 0:
                    # decide a random number of nucleosomes to be replaced
                    # centered around a poisson of half of available nucleosomes
                    num_nucs_replaced = int(np.random.poisson(lim / 2))

                    # check if we exceeded the limit
                    # unlikely but may happen bc poisson unbounded
                    if num_nucs_replaced > lim:
                        num_nucs_replaced = lim

                    # randomly sample which indicies to replace
                    nucs_replaced = random.sample(nuc_index_seq[:lim], num_nucs_replaced)

                    # replace each of the chosen indicies with a U-state
                    for nuc in nucs_replaced:
                        old = self.nucleosomes[nuc].state
                        self.update(old, States.U_STATE, nuc)

                    # skip everything else for this timestep--just go to next one
                    continue
    
            # handle recruitment
            if self.dat['r'] != Recruit.NONE and \
                    t >= self.dat['data']['recruit_time_init'] and \
                    t <= (self.dat['data']['recruit_time_init'] + 
                            self.dat['data']['recruit_time']):

                # get indicies where recruitment will happen
                start_nuc = int(self.dat['n']/2 - self.dat['data']['recruit_n']/2)
                end_nuc = start_nuc + self.dat['data']['recruit_n']

                # iterate over these indicies and simulate recruitment
                # by moving each nuc 1 step towards M
                for i in range(start_nuc, end_nuc):
                    nuc_seq_i = map_to_seq[i]
                    if nuc_seq_i < lim:
                        lim = self.handle_timers(i, self.nucleosomes[i].state, States.M_STATE, timers, nuc_index_seq, map_to_seq, lim)
                
            # choose number of events to happen in this timeslice
            num_events = int(np.random.poisson(EVENTS_PER_TIMESTEP * (lim / n_nucs)))

            # handle if poisson overshoots limit
            if num_events >= lim:
                num_events = lim 

            # select indicies to have an event
            nucs_w_event = random.sample(nuc_index_seq[:lim], num_events)

            # calculate alpha: probability of random events
            a = 1/(self.dat['f'] + 1)
            
            # choose number of random events
            num_rand_events = np.random.poisson(EVENTS_PER_TIMESTEP * (lim / n_nucs) * a)

            # handle if poisson overshoots 
            if num_rand_events > len(nucs_w_event):
                num_rand_events = len(nucs_w_event)

            # choose indicies to have a random event
            nucs_w_rand_event = random.sample(nucs_w_event, num_rand_events)
            # remainder of nucs with event which were not chosen for random
            # will have a feedback event
            nucs_w_feedback_event = list( set(nucs_w_event) - set(nucs_w_rand_event) )

            # handle all random events
            for nuc in nucs_w_rand_event:
                # get old state
                old = self.nucleosomes[nuc].state

                # if old == U-state, then we have equal chance of getting M or A, given that we
                # have a CR floating around which allows that conversion
                if old == States.U_STATE:
                    if random.random() < 2/3:
                        if random.random() < 0.5:
                            lim = self.handle_timers(nuc, old, States.A_STATE, timers, nuc_index_seq, map_to_seq, lim)
                        else:
                            lim = self.handle_timers(nuc, old, States.M_STATE, timers, nuc_index_seq, map_to_seq, lim)
                elif random.random() < 1/3:
                    lim = self.handle_timers(nuc, old, States.U_STATE, timers, nuc_index_seq, map_to_seq, lim)
            
            ## VERY IMPORTANT that M_mat and A_mat are np ARRAYS, not matricies! WE DO NOT WANT DOT PRODUCT

            # handle feedback events
            # calculate probability matricies for feedback events for M and A
            prob_conv_M = self.prob_mat[nucs_w_feedback_event,:] * self.M_mat
            prob_conv_A = self.prob_mat[nucs_w_feedback_event,:] * self.A_mat

            # get the total probability
            tot_prob_per_nuc_M = np.sum(prob_conv_M, axis = 1) / self.dat['n'] 
            tot_prob_per_nuc_A = np.sum(prob_conv_A, axis = 1) / self.dat['n']

            # iterate over all nucs with feedback events
            for nuc in range(len(nucs_w_feedback_event)):
                # get current state
                curr_state = self.nucleosomes[nucs_w_feedback_event[nuc]].state

                if curr_state == States.M_STATE:
                    # if current state is M, we can only move towards A
                    # check probability of moving to A
                    if random.random() < tot_prob_per_nuc_A[nuc] :
                        lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.U_STATE, timers, nuc_index_seq, map_to_seq, lim)
                elif curr_state == States.A_STATE:
                    # if current state is A, we can only move towards M
                    # check probability of moving towards M
                    if random.random() < tot_prob_per_nuc_M[nuc] :
                        lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.U_STATE, timers, nuc_index_seq, map_to_seq, lim)
                else: # U state
                    # if we're in U state, we can move towards M or A
                    if tot_prob_per_nuc_A[nuc] != 0 and \
                            tot_prob_per_nuc_M[nuc] != 0:

                        # get the probabilities of going to A or M
                        A_prob = tot_prob_per_nuc_A[nuc]
                        M_prob = tot_prob_per_nuc_M[nuc]
                        
                        added_prob = A_prob + M_prob

                        if added_prob > 1:
                            # normalize to 1
                            scaling = 1 / added_prob 
                            A_prob *= scaling
                            M_prob *= scaling

                        # use cumulative sum trick to pick whether we go to M or A or do nothing
                        cumsum = np.cumsum([1 - A_prob - M_prob, A_prob, M_prob])
                        int_sums = cumsum < random.random()
                        index = np.sum(int_sums.astype(int))

                        if index == 1:
                            lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.A_STATE, timers, nuc_index_seq, map_to_seq, lim)
                        elif index == 2:
                            lim = self.handle_timers(nucs_w_feedback_event[nuc], curr_state, States.M_STATE, timers, nuc_index_seq, map_to_seq, lim)
                        # else nothing

        fp.close()
    ##

    def divide(self):
        '''
        divide()
        handle divisions for each nucleosome
        '''
        for i in range(self.dat['n']):
            if random.random() <= 0.5:
                prev = self.nucleosomes[i].state
                self.nucleosomes[i].state = States.U_STATE
                self.update(prev, States.U_STATE)
                self.colors[i] = Constants.state_to_color(States.U_STATE)

    def update(self, old, new, i):
        '''
        update()
        convenient function which updates all Class-wide datastructures
        if we convert nucleosomes
        '''
        if old == new:
            return

        self.totals[old] -= 1
        self.totals[new] += 1
        self.colors[i] = Constants.state_to_color(new)
    
        self.nucleosomes[i].state = new

        if old == States.M_STATE:
            self.M_mat[i] = 0
        elif old == States.A_STATE:
            self.A_mat[i] = 0

        if new == States.M_STATE:
            self.M_mat[i] = 1
        elif new == States.A_STATE:
            self.A_mat[i] = 1


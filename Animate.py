## Animate.py
## Author: Aparna Rajpurkar

# imports
import random
import math
import operator
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.stats as stats
import Constants
from MyEnum import ProbSpread, States, Domain, DomainBleed, Divisions, ProbConv

def animate_from_file(filename, n, f, t, outfile, div_count):
    '''
    animate_from_file(simulation_filename, N_nucs, F_val, T_timesteps, basefilename, num_divisions)
    function which will animate the output from the simulation function (a file) into a
    animated plot using matplotlib
    '''
    # set some globals to use later
    # these will be used in the frame animation function--must be global
    global colors
    global totals
    global file_sim

    # initialize globals 
    colors = [Constants.GRAY]*n
    totals = {
        States.A_STATE : 0,
        States.M_STATE : 0,
        States.U_STATE : 0
            }

    # base figure
    fig = plt.figure()

    # axes
    # figure 1: nucleosomes
    ax1 = fig.add_subplot(2,2,1)
    ax1.set_xticks([])
    ax1.set_yticks([])

    titlestr = "N = " + str(n) + " F = " + str(f)
    ax1.set_title(titlestr)

    # figure 2: proportion vs events
    ax2 = fig.add_subplot(2,2,2)
    ax2.set_ylim([-5,105])
    ax2.set_xlim(t / 40, t + t / 40)
    ax2.set_xlabel("Timesteps")
    ax2.set_ylabel("% Nucleosomes")

    # adjust spacing between figures
    fig.tight_layout()

    # calculate size of square for nucleosomes in fig 1
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    # set x and y coordinates for each nucleosome
    x_int = 1 / cols
    y_int = 1 / rows

    x_vals = []
    y_vals = []

    count = 0
    for i in range(rows):
        for j in range(cols):
            if count >= n :
                break
            
            x_vals.append(j * x_int)
            y_vals.append(i * y_int)
            count += 1

    # initialize fig1 as a scatterplot
    scat = ax1.scatter(x_vals, y_vals, facecolors = colors)

    # initialize fig2 as a lineplot
    lines = []
    lines.append(ax2.plot([], [], lw = 2, color = "red", label = "A")[0])
    lines.append(ax2.plot([], [], lw = 2, color = "blue", label = "M")[0])
    ax2.legend(loc = "upper right")

    # initialize to empty data
    for line in lines:
        line.set_data([],[])

    linex = []
    liney = [[], []]

    # open the input file
    file_sim = open(filename, "r")

    # initialize
    def init_an():
        '''initialize animation: this is needed'''
        return scat, (*lines)

    # update function
    def update_an(i):
        '''update function for animation'''
        # check if this is first loop, if yes then go to next
        if i == 0:
            return scat, (*lines)

        # get a line from the file
        line = file_sim.readline()

        # initialize all totals to 0
        totals[States.A_STATE] = 0
        totals[States.U_STATE] = 0
        totals[States.M_STATE] = 0
        # init index
        index = 0

        for nuc in list(line.rstrip()):
            # go through each nucleosome and set the color & update totals
            totals[States.string_to_enum(nuc)] += 1
            colors[index] = Constants.state_to_color(States.string_to_enum(nuc))
            index += 1

            # check if we're dividing in this frame
            if div_count != 0 and i != 0:
                if i % div_count == 0:
                    # if we are dividing, add a vertical line to the plot to indicate division
                    lines.append(ax2.plot([], [], lw = 1, ls = "dotted", color = "black")[0])
                    lines[len(lines) - 1].set_data([i, i], [-5, 105])

        # append the new x,y coordinates of the plot
        linex.append(i)
        liney[0].append(totals[States.A_STATE] / n * 100)
        liney[1].append(totals[States.M_STATE] / n * 100)

        # set the data
        lines[0].set_data(linex, liney[0])
        lines[1].set_data(linex, liney[1])

        # set the colors
        scat.set_facecolors(colors)

        # return updated data
        return scat, (*lines)

    # run animation and store output in a variable
    anim = animation.FuncAnimation(fig, update_an, init_func = init_an, frames = t, interval = 1, repeat = False, blit=True)

    # write the animation to an mp4
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps = 200, metadata = dict(artist = 'Me'), bitrate = 1800)
    anim.save(outfile + ".mp4", writer = writer)

    # close the input file
    file_sim.close()


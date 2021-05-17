# adapted from Johannes Asplund Samuelsson, Ph. D. at KTH Royal Institute of Technology

# Import libraries
import numpy as np
import time
import re
import pandas as pd
import random

start_time = time.time()

from SMatrix_from_reactions import read_equations
from SMatrix_from_reactions import create_Smatrix

reactions_infile = open('reactions.txt', 'r')

reaction_dict, seen_met = read_equations(reactions_infile)
S_pandas, S = create_Smatrix(reaction_dict, seen_met) # S matrix as pandas data frame as well as as numpy array

for i in S_pandas.index:
    print(i)
print(S_pandas)

from read_concentration_ranges import read_ranges
from read_concentration_ranges import ranges_to_array
from read_concentration_ranges import ratios_to_array

conc_infile = open('concentration_ranges.txt', 'r')

conc_dict, ratio_dict, met_ratio_dict = read_ranges(conc_infile)
conc_np = ranges_to_array(conc_dict, S_pandas)
ratio_limits, ratio_mat = ratios_to_array(ratio_dict, met_ratio_dict, S_pandas)

c_lim = np.log(conc_np) # log is natural logarithm
print(c_lim)
ratio_lim = np.log(ratio_limits)

from read_equilibrator_dg import read_dg

equilibrator_infile = open('equilibrator_results.tsv', 'r')

g = read_dg(equilibrator_infile)

T=303.15
R=8.31e-3
RT = R*T

# MDF solution
start_c = np.log(np.array([7.00E-05, 0.00015, 0.01, 0.000417131, 0.004, 4.36E-07, 0.000340349, 0.0012, 0.0011, 2.80E-05, 2.50E-05, 2.00E-05, 0.0047, 1.45E-05, 0.000120239, 0.001571813, 7.47E-05, 1, 2.00E-07, 2.88E-05, 0.0008, 0.0029, 0.0099, 0.0006, 0.0099, 0.00095, 0.000209648, 1.36E-05, 0.0099, 0.007, 0.0099, 0.01625, 2.00E-07, 0.00244, 2.00E-07, 3.315E-06, 0.000407499, 0.0099, 0.0099, 0.000273, 8.00E-06, 0.0099, 0.0099]))

# Define function for random sampling of concentrations
def random_c(c_lim):
    sample = np.array([np.random.random() for n in range(0, c_lim.shape[0])])
    return sample * (c_lim[:,1] - c_lim[:,0]) + c_lim[:,0]


# Define function for checking if set is thermodynamically feasible
def df_ok(c):
    # Calculate delta G prime
    df = -(g + RT * np.sum(np.transpose(S) * c, 1)) # c is already log data
    return sum(df > 0) == df.shape[0]

# Define function for checking if set has acceptable ratios
def ratios_ok(c): 
    ratios = np.sum(ratio_mat.T * c, 1).reshape([ratio_lim.shape[0], 1]) 
    min = np.sum(np.subtract(ratios, ratio_lim) >= 0, 0)[0] == ratios.shape[0]
    max = np.sum(np.subtract(ratios, ratio_lim) <= 0, 0)[1] == ratios.shape[0]
    return min and max

# Define function for checking that sum of concentrations is not too high
def sum_ok(c, max_tot_c = 1.15): # max concentration is 1150mM = 1.15M (c(H2O) = 1M)
    t_f = np.sum(np.exp(c)) <= max_tot_c
    return np.sum(np.exp(c)) <= max_tot_c

# Define function that checks concentrations are within limits
def limits_ok(c):
    c_l = c.reshape([c.shape[0],1])
    min = np.sum(np.subtract(c_l, c_lim) >= 0, 0)[0] == c.shape[0]
    max = np.sum(np.subtract(c_l, c_lim) <= 0, 0)[1] == c.shape[0]
    return min and max

# Define function for checking feasibility, ratios, sum, and limits in one go
def is_feasible(c):
    return df_ok(c) and ratios_ok(c) and sum_ok(c[2:]) and limits_ok(c)

# Modify direction in order to get unstuck from concentration limits, a.k.a. The Unsticking Function TM
def unstick_direction(c, direction, c_lim):
    # Determine what metabolites are stuck at limits
    stuck = c.reshape((c.size,1)) == c_lim
    # Determine current signs of direction vector
    dirsign = np.sign(direction)
    # Pick a random sign for metabolites stuck at max
    max_sign = np.random.choice([-1,1], 1)
    # All directions for metabolites stuck at max must be the same sign
    dirsign[stuck[:,1] * dirsign != 0] = max_sign
    # All directions for metabolites stuck at min must be the opposite sign
    dirsign[stuck[:,0] * dirsign != 0] = -max_sign
    # Determine the directions that must change sign
    change_sign = dirsign != np.sign(direction)
    # Change the sign of directions that must change sign
    direction[change_sign] = direction[change_sign] * -1
    # Return the compatibility-modified "unstuck" direction vector
    return direction

# Define function for selecting a random direction
def random_direction(c):
    # Create a random vector of the same length as c
    direction = np.array([np.random.random() for n in range(0, c.shape[0])])
    # Subtract 0.5 to introduce negative directions
    direction = direction - 0.5
    # Set fixed concentration direction to zero
    direction[c_lim[:,1] - c_lim[:,0] == 0] = 0
    # Normalize length of direction vector
    normalized_direction = direction / np.linalg.norm(direction)
    return normalized_direction

# Define function to generate one feasible metabolite concentration set
def generate_feasible_c(c_lim):
    c = start_c
    while not is_feasible(c):
        c = random_c(c_lim) # Generate new c until 
    return c

# Determine minimum and maximum possible theta given concentration limits
def calculate_theta_hard_limit(c, direction, c_lim):
    # Find smallest fraction of direction that hits limit if added
    theta_max = np.hstack([
        (c_lim[:,1] - c)[direction != 0] / direction[direction != 0],
        (c_lim[:,0] - c)[direction != 0] / direction[direction != 0]
    ])
    theta_max = min(theta_max[theta_max >= 0])
    # Find smallest fraction of direction that hits limit if subtracted
    theta_min = np.hstack([
        (c - c_lim[:,1])[direction != 0] / direction[direction != 0],
        (c - c_lim[:,0])[direction != 0] / direction[direction != 0]
    ])
    theta_min = -min(theta_min[theta_min >= 0])
    return (theta_min, theta_max)

# Define function for determining minimum and maximum step length (theta)
def theta_range(c, direction, precision=1e-3):
    # Define function for honing in on a theta limit
    def hone_theta(theta_outer, theta_inner=0):
        if is_feasible(c + theta_outer * direction):
            # If the outer theta is feasible, accept that solution
            theta_inner = theta_outer
        else:
            while abs(theta_outer - theta_inner) > precision:
                # Calculate a theta between outer and inner limits
                theta_cur = (theta_outer + theta_inner) / 2
                if is_feasible(c + theta_cur * direction):
                    # Move outwards, set inner limit to current theta
                    theta_inner = theta_cur
                else:
                    # Move inwards, set outer limit to current theta
                    theta_outer = theta_cur
        # Return inner theta
        return theta_inner
    # Get hard limits on theta from concentrations
    theta_lim = calculate_theta_hard_limit(c, direction, c_lim)
    # Hone in on upper theta
    theta_upper = hone_theta(theta_lim[1])
    # Hone in on lower theta
    theta_lower = hone_theta(theta_lim[0])
    # Return results
    return [theta_lower, theta_upper]

# Define function for performing hit-and-run sampling within the solution space
def hit_and_run(S, g, c_lim, ratio_lim, ratio_mat, n_samples, precision=1e-3):
    # Generate starting point
    c = generate_feasible_c(c_lim)
    # Set up concentration storage list
    fMCSs = [c]
    # Perform n steps
    for i in range(0, n_samples - 1):
        # Generate random direction
        direction = random_direction(c)
        # Make sure that the algorithm doesn't get stuck at the boundaries of the solution space
        direction_unstuck = unstick_direction(c, direction,c_lim)
        # Determine minimum and maximum step length
        theta = theta_range(c, direction_unstuck, precision=precision)
        # Perform a random sampling of the step length
        theta = theta[0] + np.random.random() * (theta[1] - theta[0])
        # Perform step
        c = c + theta * direction
        # Ensure feasibility
        if not is_feasible(c):
            print("Warning: Infeasible point reached.")
            break
        # Store concentration
        fMCSs.append(c)
    # Return list of concentrations
    return fMCSs

# Define function for performing rejection sampling
def rejection_sampling(S, g, c_lim, ratio_lim, ratio_mat, n_samples):
    return [generate_feasible_c(c_lim) for i in range(0, n_samples)]

import sys

outfile = sys.argv[1] # file for further analysis in matlab in .tsv file format
outfile = open(outfile, "w")

# list of compounds (same order as in S matrix)
comp_list = seen_met.keys()

outfile.write("\t".join(comp_list) + "\n")

# start sampling 10x from same starting concentration with 500,000 steps each
# afterwards randomly select approximately 1 out of 990 until 5,000 sets of metabolite concentrations are written to output file
counter = 0
for z in range(0,10):
    print(z)
    for c in hit_and_run(S, g, c_lim, ratio_lim, ratio_mat, 500000):
        # Print CSV in mM
            if counter < 5000:
                if random.uniform(0,990) <= 1:
                    outfile.write("\t".join([str(np.exp(x)*1000) for x in c]) + "\n")
                    counter += 1
                else:
                    pass
            else:
                print("--- %s seconds ---" % (time.time() - start_time))
                outfile.close()
                exit()

#!/usr/bin/python
import os
import sys
import argparse
import json
import numpy as np

sys.path.insert(1, "./utils/")
from initialization_functions import initial_conditions_generation

# arguments
parser = argparse.ArgumentParser()

parser.add_argument("-o_dir", help="output directory")
parser.add_argument("-sim_length", type=int, help="length of histories simulation in ms")
parser.add_argument("-stim_interval", type=int, help="stimulus interval poisson parameter in ms")
parser.add_argument("-num_histories", type=int, help="number of histories to generate")
parser.add_argument("-num_preceding_stimuli", type=int, help="number of preceding stimuli to walk")

args = parser.parse_args()

# write args to log file
with open('log.txt', 'a') as log:
    log.write('args\n-o_dir: {}\n-sim_length: {}\n-stim_interval: {}\n-num_histories: {}\n-num_preceding_stimuli: {}\n'.format(args.o_dir, args.sim_length, args.stim_interval, args.num_histories, args.num_preceding_stimuli))

# generate the initial conditions for the histories
icg_histories = initial_conditions_generation(
    sim_length=args.sim_length,
    stim_interval=args.stim_interval,
    num_histories=args.num_histories,
    num_preceding_stimuli=0,
)

# generate histories
with open('log.txt', 'a') as log:
    log.write('generating histories\n')
icg_histories.gen_histories()

# calculate the time constants
with open('log.txt', 'a') as log:
    log.write('calculating time constants\n')
icg_histories.calc_t_constants()

# save histories and thier initial conditions

# save generated simulation data to files
with open('log.txt', 'a') as log:
    log.write('saving initial conditions files')
# sim_df
np.save(args.o_dir + 'initial_conditions_histories/' + 'sim_df', icg_histories.sim_df)

# stimuli
with open(args.o_dir + 'initial_conditions_histories/' + 'stimuli.txt', 'w') as writer:
    for stim in icg_histories.stimuli:
        writer.write("%s\n" % stim)

# spike_times        
with open(args.o_dir + 'initial_conditions_histories/' + 'spike_times.txt', 'w') as writer:
    for spike in icg_histories.spike_times:
        writer.write("%s\n" % spike)
        
# taus
np.save(args.o_dir + 'initial_conditions_histories/' + 'taus', icg_histories.taus)

# infs
np.save(args.o_dir + 'initial_conditions_histories/' + 'infs', icg_histories.infs)

# histories
np.save(args.o_dir + 'histories/' + 'histories', icg_histories.histories)



# generate preceding_stims
with open('log.txt', 'a') as log:
    log.write('generating preceding stimuli\n')
# generate new simulation to sample preceding stims from
icg_preceding_stimuli = initial_conditions_generation(
    sim_length=args.sim_length,
    stim_interval=args.stim_interval,
    num_histories=args.num_histories,
    num_preceding_stimuli=args.num_preceding_stimuli,
)

icg_preceding_stimuli.gen_histories()
icg_preceding_stimuli.gen_non_overlapping_preceding_stim_sets()

# save generated simulation data to files
with open('log.txt', 'a') as log:
    log.write('saving preceding stimuli to files\n')
# sim_df
np.save(args.o_dir + 'initial_conditions_preceding_stimuli/' + 'sim_df', icg_preceding_stimuli.sim_df)

# stimuli
with open(args.o_dir + 'initial_conditions_preceding_stimuli/' + 'stimuli.txt', 'w') as writer:
    for stim in icg_preceding_stimuli.stimuli:
        writer.write("%s\n" % stim)

# spike_times        
with open(args.o_dir + 'initial_conditions_preceding_stimuli/' + 'spike_times.txt', 'w') as writer:
    for spike in icg_preceding_stimuli.spike_times:
        writer.write("%s\n" % spike)
  
        
# save experiment setup trains (preceding stimuli) in json format
for i in range(len(icg_preceding_stimuli.centered_spikes)):
    
    with open(args.o_dir + 'preceding_stimuli/' + 'experiment_train_' + str(i) + '.json', 'w') as json_out:
        json.dump({
            'ind':i,
            'sampled_spike':icg_preceding_stimuli.centered_spikes[i],
            'preceding_stimuli':icg_preceding_stimuli.centered_stim_sets[i],
            'preceding_spikes':icg_preceding_stimuli.centered_spike_sets[i]
        }, json_out)

        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

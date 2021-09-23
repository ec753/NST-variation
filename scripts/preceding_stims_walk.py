#!/usr/bin/python

import os
import sys
import argparse
import json
import pickle as pkl
import numpy as np

sys.path.insert(1, "./utils/")
from experiment_functions import walking_stims_experiment

# arguments
parser = argparse.ArgumentParser()

parser.add_argument("-start", type=int, help="starting index for experiments to run") #this helps with parallelization
parser.add_argument("-stop",  type=int, help="ending index for experiments to run")

parser.add_argument("-o_dir", help="out directory")
parser.add_argument("-extend_dur", type=int, 
    help="extended duration: the amount of time to continue each simulation after the pivot stimulus (ms")

args = parser.parse_args()

# write args to log file
with open('log.txt', 'a') as log:log.write('args\n-o_dir: {}\n-extend_dur: {}\n-'.format(args.o_dir, args.extend_dur))

# load initialization data
with open('log.txt', 'a') as log:log.write('loading initialization data\n')
# histories

histories = np.load(args.o_dir + 'histories/histories.npy')

#preceding_stimuli_files = os.listdir(args.o_dir + 'preceding_stimuli')
preceding_stimuli_files = ['experiment_train_' + str(ind) + '.json' for ind in range(args.start, args.stop + 1)]

with open('log.txt', 'a') as log:log.write('beginning experiments\n')

for i, preceding_stimuli_file in enumerate(preceding_stimuli_files):
  
    with open('log.txt', 'a') as log: log.write(str(i+1)+'/'+str(len(preceding_stimuli_files))+'\n')
    
    with open(args.o_dir + 'preceding_stimuli/' + preceding_stimuli_file) as json_file:
        experiment = json.load(json_file)
    sample_spike = experiment['sampled_spike']
    preceding_stimuli = experiment['preceding_stimuli']
    ind = experiment['ind']
    
    # run the preceding stims walk
    wse = walking_stims_experiment(
        extended_duration=args.extend_dur,
        preceding_stims=preceding_stimuli,
        histories=histories
    )    
    wse.preceding_stims_walk()
    
    # save to file
    np.save(args.o_dir + 'output/' + 'rNSTs_' + str(ind) + '.txt', np.array(wse.resulting_NSTs))

    with open(args.o_dir + 'output/' + 'rPSTs_' + str(ind) + '.pkl', 'wb') as writer:
            pkl.dump(wse.resulting_PSTs, writer)
            
    
with open('log.txt', 'a') as log:
    log.write('experiments done\n')












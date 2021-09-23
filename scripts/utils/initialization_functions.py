# functions for the NST variation and AP event horizon experiments
from neuron import h
from neuron.units import mV, ms

h.load_file("stdrun.hoc")

import numpy as np
import random
import math

# simulation functions
class initial_conditions_generation:
    # functions to generate initial conditions for the walking stims experiment

    def __init__(self, sim_length, stim_interval, num_histories, num_preceding_stimuli):
        # gen_histories parameters
        self.sim_length = sim_length
        self.stim_interval = stim_interval
        self.num_histories = num_histories

        # gen_histories output
        self.sim_df = None
        self.histories = None
        self.spike_times = None
        self.stimuli = None
        self.taus = None
        self.infs = None

        # gen_preceding_stims parameters
        self.num_preceding_stimuli = num_preceding_stimuli

        # gen_preceding_stims outut
        self.centered_spikes = [] # the sampled spikes for which these trains are generated
        self.centered_stim_sets = [] # the preceding stimuli corresponding to the centered spikes
        self.centered_spike_sets = [] # the spikes that occur within the train

    def center_preceding_stims(self, spike, stims): 
        spike = (math.floor((spike - stims[-1]) * 40))/40
        stims = [(math.floor((stim - stims[-1]) * 40))/40 for stim in stims]

        return spike, stims

    def gen_non_overlapping_preceding_stim_sets(self):
        possible_sampled_spikes = {i:[] for i, spike in enumerate(self.spike_times) if spike > self.stimuli[self.num_preceding_stimuli-1]}

    	# get the preceding stims for each possible sampled spike
        for spike in possible_sampled_spikes:
            preceding_stims = [stim for stim in self.stimuli if stim < self.spike_times[spike]]
            preceding_stims = preceding_stims[-self.num_preceding_stimuli:]
            possible_sampled_spikes[spike] = preceding_stims
        
        used_spike_inds = [] # hold the sampled spikes indices that do not overlap
        used_windows = [] # hold the used windows (time spans of sampled spike and stim sets, to troubleshoot) 
        used_stim_inds = [] # hold the stimuli used for sampled spikes

        while len(possible_sampled_spikes) > 0:
            current_spike = random.choice(list(possible_sampled_spikes.keys()))

            if not bool(set(possible_sampled_spikes[current_spike]) & set(used_stim_inds)):
                # if this sampled spike's stimuli are already used, then toss it
                # otherwise we keep it
                used_stim_inds += possible_sampled_spikes[current_spike]
                used_spike_inds.append(current_spike)

                # output this set's window to check if it worked correctly
                window_start = possible_sampled_spikes[current_spike][0]
                window_end = self.spike_times[current_spike]
                used_windows.append([window_start, window_end])

                # center spikes and thier preceding stims sets around the pivot
                centered_spike, centered_stim_set = self.center_preceding_stims(self.spike_times[current_spike], possible_sampled_spikes[current_spike])
                self.centered_spikes.append(centered_spike)
                self.centered_stim_sets.append(centered_stim_set)

                # gather spikes withing this window and center them as well
                centered_spike_set = [spike for spike in self.spike_times if spike > window_start and spike < window_end]
                centered_spike_set = [(math.floor((spike - possible_sampled_spikes[current_spike][-1]) * 40))/40 for spike in centered_spike_set]
                self.centered_spike_sets.append(centered_spike_set)

            del possible_sampled_spikes[current_spike]
        
        
        with open('log.txt', 'a') as log:
            log.write('fraction of sampled spikes: {}/{} - {}\n'.format(len(used_spike_inds), len(self.spike_times), len(used_spike_inds)/len(self.spike_times)))

    def gen_histories(self):
        # Set up biophysical model
        axon = h.Section(name="axon")
        axon.insert(h.hh)

        # add a synapse
        syn = h.ExpSyn(axon(0))
        syn.tau = 1 * ms
        syn.e = 0 * mV
        syn_current = h.Vector().record(syn._ref_i)

        # add a stimulus
        stim = h.NetStim()
        stim.number = 9999999
        stim.interval = self.stim_interval * ms
        stim.noise = True
        stim.start = 0 * ms

        self.stimuli = h.Vector()

        # connect stimulus to synapse
        nc = h.NetCon(stim, syn)
        nc.delay = 0 * ms
        nc.weight[0] = 0.2
        nc.record(self.stimuli)

        # setup recording
        _t = h.Vector().record(h._ref_t)
        _v = h.Vector().record(axon(0.5)._ref_v)
        _m = h.Vector().record(axon(0.5).hh._ref_m)
        _n = h.Vector().record(axon(0.5).hh._ref_n)
        _h = h.Vector().record(axon(0.5).hh._ref_h)
        self.spike_times = h.Vector()
        nc_self = h.NetCon(axon(0.5)._ref_v, None, sec=axon)
        nc_self.record(self.spike_times)

        # run simulation
        h.finitialize(-65 * mV)
        h.continuerun(self.sim_length * ms)

        self.histories = np.array((_v, _m, _h, _n))[:, random.sample(range(0, len(_v)), self.num_histories)]

        self.sim_df = np.array((_t, _v, _m, _h, _n, syn_current))
        


    # equations for the time constants
    def vtrap(self, x, y):
        if abs(x / y) < 1e-6:
            return y * (1 - x / y / 2)
        else:
            return x / (math.exp(x / y) - 1)

    def calc_mtau(self, v, q10):
        alpha = 0.1 * self.vtrap(-(v + 40), 10)
        beta = 4 * math.exp(-(v + 65) / 18)
        _sum = alpha + beta
        mtau = 1 / (q10 * _sum)
        minf = alpha / _sum
        return mtau, minf

    def calc_htau(self, v, q10):
        alpha = 0.07 * math.exp(-(v + 65) / 20)
        beta = 1 / (math.exp(-(v + 35) / 10) + 1)
        _sum = alpha + beta
        htau = 1 / (q10 * _sum)
        hinf = alpha / _sum
        return htau, hinf

    def calc_ntau(self, v, q10):
        alpha = 0.01 * self.vtrap(-(v + 55), 10)
        beta = 0.125 * math.exp(-(v + 65) / 80)
        _sum = alpha + beta
        ntau = 1 / (q10 * _sum)
        ninf = alpha / _sum
        return ntau, ninf

    def calc_t_constants(self):
        q10 = math.pow(3, ((h.celsius - 6.3) / 10))  # some important variable
        taus = []
        infs = []
        for i in range(self.sim_df.shape[1]):
            mtau, minf = self.calc_mtau(self.sim_df[1, i], q10)
            htau, hinf = self.calc_htau(self.sim_df[1, i], q10)
            ntau, ninf = self.calc_ntau(self.sim_df[1, i], q10)

            taus.append([mtau, htau, ntau])
            infs.append([minf, hinf, ninf])

        self.taus = np.array(taus).T
        self.infs = np.array(infs).T




# functions for the NST variation and AP event horizon experiments
from neuron import h
from neuron.units import mV, ms

h.load_file("stdrun.hoc")

import numpy as np
import random
import math

# simulation functions
class walking_stims_experiment:
    # functions to perform a single preceding stims walk experiment

    def __init__(self, preceding_stims, histories, extended_duration):
        # hyperparameters
        self.extended_duration = extended_duration #amount of time after the last stimulus to run the experiment
        
        # input
        self.preceding_stims = preceding_stims
        self.histories = histories

        # outut
        self.resulting_NSTs = [] # next-spike-times
        self.resulting_PSTs = [] # pre-spike-times (those spikes that occur before the last stimulus)
        
    def sim_init(self, axon, _v_init, _m_init, _h_init, _n_init):
        # initialize simulations with new initial conditions and custom stimuli
        for seg in axon:
            seg.v = _v_init
            seg.hh.m = _m_init
            seg.hh.h = _h_init
            seg.hh.n = _n_init
        
    def seed_sim(self, history, stimuli, sim_length):

        # Set up biophysical model
        axon = h.Section(name="axon")
        axon.insert(h.hh)

        # add a synapse
        syn = h.ExpSyn(axon(0))
        syn.tau = 1 * ms
        syn.e = 0 * mV
        syn_current = h.Vector().record(syn._ref_i)

        # add custum stimuli
        netstims = [h.NetStim() for stimulus in stimuli]
        netcons = []  # store netcons here as they need to exist in memory to parallelize
        for netstim, stimulus in zip(netstims, stimuli):
            netstim.number = 1
            netstim.start = stimulus
            netcon = h.NetCon(netstim, syn)
            netcon.weight[0] = 0.2
            netcons.append(netcon)

            netcon.delay = 0 * ms

        stim_times = h.Vector()

        # setup recording
        _t = h.Vector().record(h._ref_t)
        _v = h.Vector().record(axon(0.5)._ref_v)
        _m = h.Vector().record(axon(0.5).hh._ref_m)
        _n = h.Vector().record(axon(0.5).hh._ref_n)
        _h = h.Vector().record(axon(0.5).hh._ref_h)

        spike_times = h.Vector()
        nc_self = h.NetCon(axon(0.5)._ref_v, None, sec=axon)
        nc_self.record(spike_times)

        # run simulation

        # initialize simulation
        _v_init = history[0]
        _m_init = history[1]
        _h_init = history[2]
        _n_init = history[3]

        fih = h.FInitializeHandler((self.sim_init, (axon, _v_init, _m_init, _h_init, _n_init)))
        h.finitialize()
        h.continuerun(sim_length * ms)

        return np.array((_t, _v, _m, _h, _n, syn_current)), list(spike_times)

    def run_simulation(self, stimuli):
        # function to run a simulation with a given set of histories, a set of predetermined stimuli, and an extended duration        
        exp_out = []
        
        for i in range(self.histories.shape[1]):
            history = self.histories[:,i]

            sim_out, spike_times = self.seed_sim(history, stimuli, sim_length = self.extended_duration + stimuli[-1])
            exp_out.append(spike_times)
            
        return exp_out
    
    def preceding_stims_walk(self):
        # function to walk through each of the preceding stimuli

        for num_stims_included in range(1, len(self.preceding_stims) + 1):
            current_stims = self.preceding_stims[len(self.preceding_stims) - num_stims_included:]

            # set the first stim to t0
            current_stims = [cs - current_stims[0] for cs in current_stims]

            #print('number of preceding stims:', len(current_stims))
            exp_out = self.run_simulation(current_stims)
            NSTs = [] # next-spike-times
            PSTs = [] # pre-spike-times (those spikes that occur before the last stimulus)

            # process each exp_out
            for spike_times in exp_out:
                # separate spike_times before and after the last stimulus
                pre_spike_times = [sp for sp in spike_times if sp <= current_stims[-1]]
                post_spike_times = [sp for sp in spike_times if sp > current_stims[-1]]
                
                #print(len(spike_times), len(pre_spike_times), len(post_spike_times))

                if len(post_spike_times) > 0:
                    NSTs.append(post_spike_times[0])
                else:
                    NSTs.append('na') # did not spike
                    
                PSTs += pre_spike_times 

            # set the last stimulus to t0 (re-center)
            _NSTs = []
            for nst in NSTs:
                if nst != 'na':
                    nst = nst - current_stims[-1]
                _NSTs.append(nst)
            
            PSTs = [pst - current_stims[-1] for pst in PSTs]

            self.resulting_NSTs.append(_NSTs)
            self.resulting_PSTs.append(PSTs)
        


# visualization functions
def view_preceding_stimuli(preceding_stims, pivot):
    # view the preceding stimuli
    plt.figure(figsize=(15, 3))
    plt.vlines(preceding_stims, 0, 1, color="black")
    plt.vlines(pivot, 0, 1, color="red")
    plt.xlabel("time (ms)")
    plt.yticks([])
    plt.show()


def view_preceding_stims_walk(preceding_stims, resulting_NSTs, rPSTs, pivot):
    fig, axes = plt.subplots(len(preceding_stims) - 1, 1, figsize=(15, 10), sharex=True)

    # get the latest spike time to set the xlim
    latest_spike_time = np.nanmax(resulting_NSTs)

    for i in range(1, len(preceding_stims)):
        i_ax = i - 1
        # xlims
        left_buffer = 2
        axes[i_ax].set_xlim(
            [math.floor(preceding_stims[0]) - left_buffer, latest_spike_time + 2]
        )

        # uncertain history window
        left = preceding_stims[0] - left_buffer
        bottom = 0
        width = (
            left_buffer
            + preceding_stims[len(preceding_stims) - 1 - i]
            - preceding_stims[0]
        )
        height = 1

        # add an extra 10 buffer so there isn't any white space
        axes[i_ax].add_patch(
            matplotlib.patches.Rectangle(
                (left - 10, bottom), width + 10, height, color="black", alpha=0.5
            )
        )

        # certain history window
        left = left + width
        width = 0 - left
        axes[i_ax].add_patch(
            matplotlib.patches.Rectangle(
                (left, bottom), width, height, color="red", alpha=0.2
            )
        )

        # output window
        left = left + width
        width = 100
        axes[i_ax].add_patch(
            matplotlib.patches.Rectangle(
                (left, bottom), width, height, color="gold", alpha=0.3
            )
        )

        # output spikes
        axes[i_ax].vlines(resulting_NSTs[i, :], 0, 1, color="green")

        # pre-last-stimulus spikes
        # axes[i_ax].vlines(rPSTs[i], 0, 1, color = 'limegreen')

        # stimuli
        axes[i_ax].vlines(preceding_stims[-i - 1 :], 0, 1, color="blue")
        axes[i_ax].vlines(preceding_stims[-i - 1 :][-1], 0, 1, color="cyan")

        # pivot spike
        axes[i_ax].vlines(pivot, 0, 1, color="fuchsia")

        # set axes ticks and such
        if i != len(preceding_stims) - 1:
            axes[i_ax].set_xticks([])
        else:
            axes[i_ax].set_xlabel("time (ms)")

        axes[i_ax].set_ylabel(str(i + 1) + "     ", rotation=0)
        axes[i_ax].set_yticks([])

        axes[i_ax].set_ylim(0.1, 0.2)
    last_x_tick = math.floor(preceding_stims[0] / 10) * 10
    axes[-1].set_xticks(range(last_x_tick, 20, 10))
    plt.show()

    return


def view_NST_hists(preceding_stims, resulting_NSTs, pivot):
    all_resulting_NSTs = np.array([st for sts in resulting_NSTs[1:] for st in sts])
    # remove nans
    all_resulting_NSTs = all_resulting_NSTs[~np.isnan(all_resulting_NSTs)]

    # run histogram on all output spikes to generate the bins for the specific histograms
    # n, bins, patches = plt.hist([st for sts in resulting_NSTs[1:] for st in sts], bins = 100)
    count, bins = np.histogram(all_resulting_NSTs, 100)

    # NST hist
    fig, axes = plt.subplots(len(preceding_stims) - 1, 1, figsize=(10, 15), sharex=True)

    # get the latest spike time to set the xlim
    latest_spike_time = max([st for sts in resulting_NSTs for st in sts])

    for i in range(1, len(preceding_stims)):
        i_ax = i - 1
        n, bins, patches = axes[i_ax].hist(resulting_NSTs[i], bins, color="darkgreen")
        axes[i_ax].vlines(pivot, -1000000, 100000, color="fuchsia")
        axes[i_ax].set_ylabel(str(i + 1) + "     ", rotation=0)
        axes[i_ax].set_yticks([])

        # window
        left = 0
        bottom = -1000
        width = 100
        height = 100000
        axes[i_ax].add_patch(
            matplotlib.patches.Rectangle(
                (left, bottom), width, height, color="gold", alpha=0.3
            )
        )
        axes[i_ax].set_ylim(0, max(patches.datavalues) * 1.1)

    axes[i_ax].set_xlim(min(bins), max(bins))

    axes[i_ax].set_xlabel("time (ms)")
    plt.show()

    return patches

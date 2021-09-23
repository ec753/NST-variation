#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=nst
#SBATCH --cpus-per-task 20

module load miniconda

conda activate nst_experiments

# create output directory

data_dir='./data'

# experiment setup hyperparameters
sim_length=600000
stim_interval=10
num_histories=10000
extend_dur=20
num_preceding_stimuli=20

printf "setting up data directory" > log.txt
rm -rf $data_dir

mkdir $data_dir
mkdir $data_dir/initial_conditions_histories
mkdir $data_dir/initial_conditions_preceding_stimuli
mkdir $data_dir/histories
mkdir $data_dir/preceding_stimuli
mkdir $data_dir/output

printf "\nSTEP 1: initial conditions" >> log.txt
python3 gen_initial_conditions.py -o_dir $data_dir/ -sim_length $sim_length -stim_interval $stim_interval -num_histories $num_histories -num_preceding_stimuli $num_preceding_stimuli

printf "\nSTEP 2: running experiments" >> log.txt
#python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur
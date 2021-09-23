#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=nst
#SBATCH --cpus-per-task 20

module load miniconda

conda activate nst_experiments

data_dir='./data'

# experiment setup hyperparameters
sim_length=600000
stim_interval=10
num_histories=10000
extend_dur=20
num_preceding_stimuli=20

printf "\nSTEP 2: running experiments" >> log.txt
echo "start: $1"
echo "stop: $2"

python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 500 -stop 510 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 510 -stop 520 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 520 -stop 530 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 530 -stop 540 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 540 -stop 550 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 550 -stop 560 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 560 -stop 570 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 570 -stop 580 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 580 -stop 590 &
python3 preceding_stims_walk.py -o_dir $data_dir/ -extend_dur $extend_dur -start 590 -stop 600 &

wait;echo "oink"

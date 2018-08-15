#!/bin/bash

RUN_ID=2
PACKING_FRACTION=0.6
PARTICLE_NUMBER=504
SIMULATION_TIME=4085
ASPECT_RATIO=2.8
BROWNIAN_TIMESTEP=1
SAVING_TIMESTEP=0.001
DATA="/data/scc/aleena/check_k_v4/equilibrated/mw-run${RUN_ID}-n${PARTICLE_NUMBER}-t${SIMULATION_TIME}-phi${PACKING_FRACTION}-k${ASPECT_RATIO}-tb${BROWNIAN_TIMESTEP}-ts${SAVING_TIMESTEP}"


#$ -N proc_msd # the name of the job
#$ -cwd          # run in current directory
#$ -m eba        # send a mail at end ...
#$ -M aleena.laganapan@uni.kn # specify mail address
#$ -e data/logs/  # specify directory for error output
#$ -o data/logs/  # specify directory for standard output
#$ -q scc
##$ -q grendel # TODO: check why this cannot be included in above -q option
##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
#$ -l h_vmem=1G
#$ -l h_rt=10:00:00
##$ -notify         

# local job
SGE_TASK_ID=1

### uncomment when you need modules
# . /etc/profile.d/modules.sh
module purge
#module load python/3.3.3
#module load python/3
#module load numpy/1.9.2-python33-intel13
#module load scipy/0.16.0-python33-intel13
#module load h5py/2.5.0-python33-intel13
#module load simpy/3.0.8-python3

#cmd=python3
cmd="python3 msd_cumulative_average_v2.py"
args="--source_directory $DATA --natoms $PARTICLE_NUMBER"
echo ${cmd} ${args}
exec ${cmd} ${args}

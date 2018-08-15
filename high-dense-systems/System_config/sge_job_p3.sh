#!/bin/bash
#$ -N pyBDsimP3n4T300phi074 # the name of the job
#$ -cwd          # run in current directory
#$ -m eba        # send a mail at end ...
#$ -M markus.heinrich@uni.kn # specify mail address
#$ -e data/logs/  # specify directory for error output
#$ -o data/logs/  # specify directory for standard output
#$ -q scc,long
##$ -q grendel # TODO: check why this cannot be included in above -q option
##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
#$ -l h=!scc013&!scc051
#$ -l h_vmem=3G   
#$ -l h_rt=20:00:00
##$ -notify         
#$ -t 1-100

# local job
#SGE_TASK_ID=14

### uncomment when you need modules
# . /etc/profile.d/modules.sh
module purge
#module load python/3.3.3
#module load python/3
#module load numpy/1.9.2-python33-intel13
module load scipy/0.16.0-python33-intel13
module load h5py/2.5.0-python33-intel13
module load simpy/3.0.8-python3

#cmd=python3
cmd="python3 p3pin.py"
DATA="/data/scc/mahe/p3-pinned/p3-n4-t300-phi0.74-tb0.1-ts0.1"
mkdir -p $DATA/logs
args="--path $DATA --brownian-timestep 0.1 --saving-timestep 0.1 --system-lifetime 300 --packing-fraction 0.74 --p3-iterations 4 --simulation-id $SGE_TASK_ID --verbose"
echo ${cmd} ${args}
exec ${cmd} ${args}

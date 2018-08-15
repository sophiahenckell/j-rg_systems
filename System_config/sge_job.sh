#!/bin/bash
#$ -N swell1000 # the name of the job
#$ -cwd          # run in current directory
#$ -m eba        # send a mail at end ...
#$ -M aleena.laganapan@uni.kn # specify mail address
#$ -e data/logs/  # specify directory for error output
#$ -o data/logs/  # specify directory for standard output
#$ -q scc,long
##$ -q grendel # TODO: check why this cannot be included in above -q option
##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
#$ -l h_vmem=12G   
#$ -l h_rt=72:00:00
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
module load simpy/3.0.10-python3
module load cython
module list

#cmd=python3
cmd="python3 mw.py"
#mkdir -p ./data/logs
#args="--path ./data --brownian-timestep 0.1 --saving-timestep 0.02 --system-lifetime 100 --packing-fraction 0.02 --particle-number 512 --aspect-ratio 1.71 --simulation-id $SGE_TASK_ID --verbose"
echo ${cmd} ${args}
exec ${cmd} ${args}

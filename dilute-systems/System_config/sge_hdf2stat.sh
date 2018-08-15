#!/bin/bash
##$ -N pyBDsimH2SRun1T300N1024k18Phi05 # the name of the job
#$ -N pyBDsimH2SRun1T300N1024k18Phi06 # the name of the job
#$ -cwd          # run in current directory
#$ -m eba        # send a mail at end ...
#$ -M sophia.henckell@uni.kn # specify mail address
#$ -e data/logs/testrun/  # specify directory for error output
#$ -o data/logs/testrun/  # specify directory for standard output
#$ -q scc
##$ -q grendel # TODO: check why this cannot be included in above -q option
##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
#$ -l h=!scc013&!scc051
#$ -l h_vmem=1G   
#$ -l h_rt=1:00:00
##$ -notify         
##$ -t 1-100

# local job
#SGE_TASK_ID=1

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
cmd="python3 hdf2stat.py"
#DATA="/data/scc/sophia/sophia/hello"
DATA="/data/scc/sophia/hello/"
mkdir -p $DATA/logs
for f in $DATA/*.h5
do
	args="--filename $f"
	echo ${cmd} ${args}
	`${cmd} ${args}`
done

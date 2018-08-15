#!/bin/bash

#ACHTUNG! LOG SCALE MUSS für den submit EINGESTELLT SEIN!!!!
RUN_ID=1
TASKS=1
PACKING_FRACTION=0.4
PARTICLE_NUMBER=1372
SIMULATION_TIME=10
ASPECT_RATIO=0
RADIUS=1
BROWNIAN_TIMESTEP=0.01
SAVING_TIMESTEP=0.0001 
SWELLING_RATE=0
NAME="Brown_Test"
DATA="/data/scc/sophia/AA_Brown_test_free_particle"
mkdir -p $DATA/logs3

TMPFILE="tmp-${RUN_ID}.sh"

cat > $TMPFILE <<- EOF
	#$ -N ${NAME} # the name of the job
	#$ -cwd          # run in current directory
	#$ -m eba        # send a mail at end ...
	#$ -M sophia.henckell@uni.kn # specify mail address
	#$ -e ${DATA}/logs3/  # specify directory for error output
	#$ -o ${DATA}/logs3/  # specify directory for standard output
	#$ -q scc,long,longer  
	##$ -q grendel # TODO: check why this cannot be included in above -q option
	##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
	##$ -l h=!scc013&!scc051
	#$ -l h_vmem=10G   
	#$ -l h_rt=10:00:00
	##$ -notify         
	#$ -t 1-${TASKS}

	# local job
	#SGE_TASK_ID=1

	### uncomment when you need modules
	module purge
	module load simpy/3.0.10-python3
	module load cython
	module list
	
	##echo "$JOB_ID.$SGE_TASK_ID $DATA" >> /data/scc/sophia/globallogs/clustersubmit.log

	cmd="python3 mw_log.py"

EOF

qsub $TMPFILE
rm $TMPFILE

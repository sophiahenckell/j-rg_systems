#!/bin/bash

#ACHTUNG! LOG SCALE MUSS für den submit EINGESTELLT SEIN!!!!
RUN_ID=1
TASKS=250
BROWNIAN_TIMESTEP=0.01
SAVING_TIMESTEP=0.001
SWELLING_RATE=0
SIMULATION_TIME=2000
PACKING_FRACTION=0.19
PARTICLE_NUMBER=1372
ASPECT_RATIO=0
RADIUS=1
 

NAME="SYSTEM3${RUN_ID}T${SIMULATION_TIME}N${PARTICLE_NUMBER}k${ASPECT_RATIO//.}Phi${PACKING_FRACTION//.}"
DATA="/data/scc/sophia/PaperSystems/SYSTEM3/${RUN_ID}-n${PARTICLE_NUMBER}-t${SIMULATION_TIME}-phi${PACKING_FRACTION}-k${ASPECT_RATIO}-tb${BROWNIAN_TIMESTEP}-ts${SAVING_TIMESTEP}-r${RADIUS}"
mkdir -p $DATA/logs

TMPFILE="tmp-${RUN_ID}.sh"

cat > $TMPFILE <<- EOF
	#$ -N ${NAME} # the name of the job
	#$ -cwd          # run in current directory
	#$ -m eba        # send a mail at end ...
	#$ -M sophia.henckell@uni.kn # specify mail address
	#$ -e ${DATA}/logs/  # specify directory for error output
	#$ -o ${DATA}/logs/  # specify directory for standard output
	#$ -q scc,long,longer  
	##$ -q grendel # TODO: check why this cannot be included in above -q option
	##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
	##$ -l h=!scc013&!scc051
	#$ -l h_vmem=15G   
	#$ -l h_rt=719:00:00
	##$ -notify         
	#$ -t 1-${TASKS}

	# local job
	#SGE_TASK_ID=1

	### uncomment when you need modules
	module purge
	module load simpy/3.0.10-python3
	module load cython
	module list
	
	echo "Joerg Systems $JOB_ID.\$SGE_TASK_ID $DATA" >> /data/scc/sophia/globallogs/clustersubmit.log

	cmd="python3 mw_log.py"
	args="--path $DATA --brownian-timestep ${BROWNIAN_TIMESTEP} --saving-timestep ${SAVING_TIMESTEP} --swelling-rate ${SWELLING_RATE} --system-lifetime ${SIMULATION_TIME} --packing-fraction ${PACKING_FRACTION} --particle-number ${PARTICLE_NUMBER} --aspect-ratio ${ASPECT_RATIO} --radius ${RADIUS} --simulation_id \$SGE_TASK_ID --verbose"
	echo \${cmd} \${args}
	exec \${cmd} \${args}
EOF

qsub $TMPFILE
rm $TMPFILE

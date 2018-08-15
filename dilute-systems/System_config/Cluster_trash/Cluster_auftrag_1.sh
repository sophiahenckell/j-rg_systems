#!/bin/bash
RUN_ID=3
TASKS=1
PACKING_FRACTION=0.222
PARTICLE_NUMBER=500
SIMULATION_TIME=4000
ASPECT_RATIO=1.2
BROWNIAN_TIMESTEP=0.01
SAVING_TIMESTEP=10
SWELLING_RATE=10
NAME="ZweiteMSDBerechnungen${RUN_ID}T${SIMULATION_TIME}N${PARTICLE_NUMBER}k${ASPECT_RATIO//.}Phi${PACKING_FRACTION//.}"
DATA="/data/scc/sophia/Second_MSDCalculations${RUN_ID}-n${PARTICLE_NUMBER}-t${SIMULATION_TIME}-phi${PACKING_FRACTION}-k${ASPECT_RATIO}-tb${BROWNIAN_TIMESTEP}-ts${SAVING_TIMESTEP}"
mkdir -p $DATA/logs

TMPFILE="tmp-${RUN_ID}.sh"

cat > $TMPFILE <<- EOF
	#$ -N ${NAME} # the name of the job
	#$ -cwd          # run in current directory
	#$ -m eba        # send a mail at end ...
	#$ -M sophia.henckell@uni.kn # specify mail address
	#$ -e ${DATA}/logs/  # specify directory for error output
	#$ -o ${DATA}/logs/  # specify directory for standard output
	#$ -q scc  
	##$ -q grendel # TODO: check why this cannot be included in above -q option
	##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
	##$ -l h=!scc013&!scc051
	#$ -l h_vmem=5G   
	#$ -l h_rt=167:00:00
	##$ -notify         
	#$ -t 1-${TASKS}

	# local job
	#SGE_TASK_ID=1

	### uncomment when you need modules
	module purge
	module load simpy/3.0.10-python3
	module load cython
	module list
	
	echo "$JOB_ID.$SGE_TASK_ID $DATA" >> /data/scc/sophia/globallogs/clustersubmit.log

	cmd="python3 mw.py"
	args="--path $DATA --brownian-timestep ${BROWNIAN_TIMESTEP} --saving-timestep ${SAVING_TIMESTEP} --swelling-rate ${SWELLING_RATE} --system-lifetime ${SIMULATION_TIME} --packing-fraction ${PACKING_FRACTION} --particle-number ${PARTICLE_NUMBER} --aspect-ratio ${ASPECT_RATIO} --simulation-id \$SGE_TASK_ID --verbose"
	echo \${cmd} \${args}
	exec \${cmd} \${args}
EOF

qsub $TMPFILE
rm $TMPFILE

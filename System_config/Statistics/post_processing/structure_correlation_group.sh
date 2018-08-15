#!/bin/bash
RUN_ID=2
TASKS=5
PACKING_FRACTION=0.4
PARTICLE_NUMBER=504
SIMULATION_TIME=4085
ASPECT_RATIO=2.8
BROWNIAN_TIMESTEP=0.001
SAVING_TIMESTEP=0.001
NAME="sk_calc"
DATA="/data/scc/aleena/check_k_v4/equilibrated/mw-run${RUN_ID}-n${PARTICLE_NUMBER}-t${SIMULATION_TIME}-phi${PACKING_FRACTION}-k${ASPECT_RATIO}-tb${BROWNIAN_TIMESTEP}-ts${SAVING_TIMESTEP}"
mkdir -p $DATA/logs2

TMPFILE="tmp-${RUN_ID}.sh"

cat > $TMPFILE <<- EOF
	#$ -N ${NAME} # the name of the job
	#$ -cwd          # run in current directory
	##$ -m eba        # send a mail at end ...
	##$ -M aleena.laganapan@uni.kn # specify mail address
	#$ -e ${DATA}/logs2/  # specify directory for error output
	#$ -o ${DATA}/logs2/  # specify directory for standard output
	#$ -q scc
	##$ -q grendel # TODO: check why this cannot be included in above -q option
	##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
	#$ -l h=!scc013&!scc051
	#$ -l h_vmem=1G   
	#$ -l h_rt=10:00:00
	##$ -notify         
	#$ -t 1-${TASKS}

	# local job
	#SGE_TASK_ID=1

	### uncomment when you need modules
	module purge

	cmd="python3 structure_correlation_group.py"
	args="--path $DATA --brownian-timestep ${BROWNIAN_TIMESTEP} --saving-timestep ${SAVING_TIMESTEP} --system-lifetime ${SIMULATION_TIME} --packing-fraction ${PACKING_FRACTION} --particle-number ${PARTICLE_NUMBER} --aspect-ratio ${ASPECT_RATIO} --simulation-id \$SGE_TASK_ID --verbose"
	echo \${cmd} \${args}
	exec \${cmd} \${args}
EOF

qsub $TMPFILE
rm $TMPFILE

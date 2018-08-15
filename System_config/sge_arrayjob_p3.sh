#!/bin/bash
RUN_ID=106
TASKS=1
PACKING_FRACTION=0.785398
P3_ITERATIONS=6
SIMULATION_TIME=50
BROWNIAN_TIMESTEP=0.1
SAVING_TIMESTEP=0.1
NAME="pyBDsimP3Run${RUN_ID}T${SIMULATION_TIME}N${P3_ITERATIONS}Phi${PACKING_FRACTION//.}"
DATA="/data/scc/mahe/p3-pinned/p3-run${RUN_ID}-n${P3_ITERATIONS}-t${SIMULATION_TIME}-phi${PACKING_FRACTION}-k${ASPECT_RATIO}-tb${BROWNIAN_TIMESTEP}-ts${SAVING_TIMESTEP}"
mkdir -p $DATA/logs

TMPFILE="tmp-${RUN_ID}.sh"

cat > $TMPFILE <<- EOF
	#$ -N ${NAME} # the name of the job
	#$ -cwd          # run in current directory
	#$ -m eba        # send a mail at end ...
	#$ -M markus.heinrich@uni.kn # specify mail address
	#$ -e ${DATA}/logs/  # specify directory for error output
	#$ -o ${DATA}/logs/  # specify directory for standard output
	#$ -q scc,long
	##$ -q grendel # TODO: check why this cannot be included in above -q option
	##$ -l h=!ramsey&!hooke&!gilbert&!pauli&!fermat&!hund&!fano&!onsager&!geiger&!huygens&!kerr&!hertz&!brown
	#$ -l h=!scc013&!scc051
	#$ -l h_vmem=2G   
	#$ -l h_rt=80:00:00
	##$ -notify         
	#$ -t 1-${TASKS}

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

	cmd="python3 p3pin.py"
	args="--path $DATA --brownian-timestep ${BROWNIAN_TIMESTEP} --saving-timestep ${SAVING_TIMESTEP} --system-lifetime ${SIMULATION_TIME} --packing-fraction ${PACKING_FRACTION} --p3-iterations ${P3_ITERATIONS} --simulation-id \$SGE_TASK_ID --verbose"
	echo \${cmd} \${args}
	exec \${cmd} \${args}
EOF

qsub $TMPFILE
rm $TMPFILE

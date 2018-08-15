 for i in {1..500}
 do

 #!/bin/bash
ID=$i
RUN_ID=1
TASKS=1
PACKING_FRACTION=0.494
PARTICLE_NUMBER=1372
SIMULATION_TIME=2000
ASPECT_RATIO=0
BROWNIAN_TIMESTEP=0.01
SAVING_TIMESTEP=0.001
NAME="sk_calc"
DATA="/data/scc/sophia/check_incoherent/equilibrated/mw-run${RUN_ID}-n${PARTICLE_NUMBER}-t${SIMULATION_TIME}-phi${PACKING_FRACTION}-k${ASPECT_RATIO}-tb${BROWNIAN_TIMESTEP}-ts${SAVING_TIMESTEP}"
mkdir -p $DATA/logs

TMPFILE="tmp-${ID}.sh"

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
	#$ -l h=!scc013&!scc051
	#$ -l h_vmem=1G   
	#$ -l h_rt=2:00:00
	##$ -notify         
	#$ -t 1-${TASKS}

# 	# local job
	#SGE_TASK_ID=1

	### uncomment when you need modules
	module purge

	cmd="python3 structure_correlation_group_NEW.py --id $ID"

	echo \${cmd} 
	exec \${cmd} 
EOF

qsub $TMPFILE
rm $TMPFILE

done

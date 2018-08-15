#!/bin/bash

#ACHTUNG! LINEARER SCALE MUSS für den submit EINGESTELLT SEIN!!!!
RUN_ID=1
TASKS=1

NAME="stat4.2_II${RUN_ID}T"
DATA="/data/scc/sophia/StaticStructureFig4_2_II_neu${RUN_ID}"
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
	#$ -l h_vmem=5G   
	#$ -l h_rt=168:00:00
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

	cmd="python3 ./Statistics/post_processing/structure_correlation_group.py"
	
	echo \${cmd} \${args}
	exec \${cmd} \${args}
EOF

qsub $TMPFILE
rm $TMPFILE

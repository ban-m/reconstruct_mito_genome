#!/bin/bash
#$ -S /bin/bash
#$ -N annotation
#$ -cwd
#$ -pe smp 12
#$ -e ./logfiles/annotate.log
#$ -o ./logfiles/annotate.out
#$ -V
#$ -m e
# -M banmasutani@gmail.com ## To send a mail when ended
# -t 1:n ## For array job
## "-S" specifies which shell do the script run on.
## "-N" specifies the name of the job(should NOT be TEMPLATE).
## "-cwd" means the script run on the current directly, not on $HOME.
## "-pe smp N" specifies the number of the core the script requires.
## with "-e [path]", the stderr would be written in [path], not ${HOME}/$JOB_NUMBER
## with "-o [path]" is the same as "-e", except for stdout.
## with "-V", all envirnmental variables would be inherited.
## "with -m e", a mail would be sent if the job halt with an error.
## For more detail, see `man qsub`
OUTPUT_DIR=/work/ban-m/arabidopsis_mitochondria/annotate_mitochondrial_sequence/result/test
PREFIX=test
DATABASE=./data/plant_mitochondria_from_genbank.fa
ANNOATION_FILE=./data/plant_mitochondria_from_genbank.tab

bash ./script/annotate.sh $1 ${OUTPUT_DIR} ${PREFIX} ${DATABASE} ${ANNOATION_FILE}

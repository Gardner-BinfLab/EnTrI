#$ -S /bin/sh
#$ -e /home/fas31/logs
#$ -o /home/fas31/logs

esl-reformat -o /data/people/fas31/homology/cmalign.out.phylip phylip /data/people/fas31/homology/cmalign.out

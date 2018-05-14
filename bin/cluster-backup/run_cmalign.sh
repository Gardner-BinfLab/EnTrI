#$ -S /bin/sh
#$ -e /home/fas31/logs
#$ -o /home/fas31/logs

cmalign -o /data/people/fas31/homology/cmalign.out /data/people/fas31/homology/RF00177.cm /data/people/fas31/homology/cmsearch-sequences.fa

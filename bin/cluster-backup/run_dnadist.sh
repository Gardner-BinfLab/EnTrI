#$ -S /bin/sh
#$ -e /home/fas31/logs
#$ -o /home/fas31/logs

dnadist < /data/people/fas31/homology/dnadist.options
mv outfile /data/people/fas31/homology/dnadist.out

#$ -S /bin/sh
#$ -e /home/fas31/logs
#$ -o /home/fas31/logs

export PATH=/home/fas31/anaconda3/bin:/home/fas31/local/bin:$PATH

for filename in /data/genomes/bacteria2016-04-14/*.embl; do /data/people/fas31/homology/embl2peptides.py $filename /data/genomes/bacteria2016-04-14/fasta-prot/"$(basename $filename | cut -d. -f1).fasta"; done

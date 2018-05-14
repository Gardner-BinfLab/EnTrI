#$ -S /bin/sh
#$ -e /home/fas31/logs
#$ -o /home/fas31/logs

for filename in /data/genomes/bacteria2016-04-14/fasta-prot/*.fasta; do /home/fas31/program-bank/bin/hmmsearch -o /dev/null --noali --domtblout /data/people/fas31/homology/hmmsearch/$(basename $filename | cut -d. -f1).domtblout /data/people/fas31/homology/corefam.hmm $filename; done

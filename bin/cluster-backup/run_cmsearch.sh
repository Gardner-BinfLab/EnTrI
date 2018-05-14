#$ -S /bin/sh
#$ -e /home/fas31/logs
#$ -o /home/fas31/logs

for filename in /data/genomes/bacteria2016-04-14/fasta/*.fasta; do cmsearch -o /dev/null --noali --tblout /data/people/fas31/homology/cmsearch/$(basename $filename | cut -d. -f1).tblout /data/people/fas31/homology/RF00177.cm $filename; done

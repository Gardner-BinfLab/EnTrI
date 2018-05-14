1- qsub -v PATH=$PATH -l h_vmem=32g embl2peptide.sh
2- qsub -v PATH=$PATH -l h_vmem=32g run_hmmsearch.sh
3- qsub -v PATH=$PATH -l h_vmem=32g run_con-ess.sh
4- Find SSU rRNA from rfam and save it in RF00177.cm
5- qsub -v PATH=$PATH -l h_vmem=32g run_cmsearch.sh
6- qsub -v PATH=$PATH -l h_vmem=32g run_read_cmsearch_output.sh
7- qsub -v PATH=$PATH -l h_vmem=32g run_cmalign.sh
8- qsub -v PATH=$PATH -l h_vmem=32g run_esl-reformat.sh
9- qsub -v PATH=$PATH -l h_vmem=32g run_dnadist.sh
10- qsub -v PATH=$PATH -l h_vmem=32g run_con-ess-areba2.sh

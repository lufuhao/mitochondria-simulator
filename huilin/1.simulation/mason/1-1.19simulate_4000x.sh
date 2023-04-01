#!/bin/bash
#SBATCH -o %x-%J.20210119.log
#SBATCH -e %x-%J.20210119.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J mason_simulate.4000x
#SBATCH -p fat
#SBATCH --mail-user=huilin23@outlook.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1




module add  mason/v0.1.2

#Mt:READNUM=452528*4000/(2*150)=6033706
#Ct:READNUM=135885*500/（2*150）=226475
#Gene:READNUM=14300719029*500/(2*150)=23834531715
##Mt模拟4000x的
mason illumina -N 6033706 -aNg -sq  -n 150 -ll 500 -le 200 -mp --read-naming 2 --read-name-prefix mt.4000x.6033706 -o /home/hpcusers/master/huilin_hu/simulate/Second_simulate/mt_ct_gn/Triticum.chrM.4000x.fastq  /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.mason.500bp.fa
##叶绿体和基因组模拟500x的
mason illumina -N 226475 -aNg -sq -n 150 -ll 500 -le 200 -mp --read-naming 2 --read-name-prefix chloroplast.500x.226475 -o /home/hpcusers/master/huilin_hu/simulate/Second_simulate/500x_ct_gn/Triticum_cs.chrC.500x.fastq /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrC.fasta
##gene
mason illumina -N 23834531715 -aNg -sq -n 150 -ll 500 -le 200 -mp --read-naming 2 --read-name-prefix gene.500x.23834531715 -o /home/hpcusers/master/huilin_hu/simulate/Second_simulate/500x_ct_gn/Triticum_cs.gene.500x.fastq /home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta

#!/bin/bash
#SBATCH -o %x-%J.20201121.log
#SBATCH -e %x-%J.20201121.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J bwa_mason
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add bwa/v0.7.17-13b5637
module add samtools/v0.1.20

COV=$1
mkdir /home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x/1.bwa${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x/1.bwa${COV}x




for i in 1 2 3
do
	bwa mem -t 6 /home/hpcusers/master/huilin_hu/simulate/bwa_index/Triticum_cs.chrM.mason.500bp.fa  ${HOMEDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-10.trim.fastq.gz  ${HOMEDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-10.trim.fastq.gz | samtools view -Sb - > ${OUTDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-10.trim.bam
	samtools sort ${OUTDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-10.trim.bam ${OUTDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-10.trim.sort
done

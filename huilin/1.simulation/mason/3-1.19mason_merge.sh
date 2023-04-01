#!/bin/bash
#SBATCH -o %x-%J.20210119.log
#SBATCH -e %x-%J.20210119.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J merge-fastqc.sh
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1

module add seqkit/v0.13.2
module add fastqc/v0.11.9

COV=$1
RATIO=$2
COV1=$(echo "scale=0;$COV/20" |bc -l)
mkdir /home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x
mkdir /home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/0.fastqc_merge${COV}x

#such as:sbatch 11.06merge_mason.sh 500 10

for i in 1 2 3 4 5
do
	cat ${HOMEDIR}/Triticum.chrC.${COV1}x.${i}_1-${RATIO}_1.fastq ${HOMEDIR}/Triticum.gene.${COV1}x.${i}_1-${RATIO}_1.fastq >${OUTDIR}/starter_${i}_1-${RATIO}.${COV1}x_1.fastq
	cat ${HOMEDIR}/Triticum.chrC.${COV1}x.${i}_1-${RATIO}_2.fastq ${HOMEDIR}/Triticum.gene.${COV1}x.${i}_1-${RATIO}_2.fastq >${OUTDIR}/starter_${i}_1-${RATIO}.${COV1}x_2.fastq

	cat ${OUTDIR}/starter_${i}_1-${RATIO}.${COV1}x_1.fastq  ${HOMEDIR}/Triticum.chrM.${COV}x.${i}_1-${RATIO}_1.fastq  | seqkit shuffle --threads 8 --rand-seed 9000 --out-file  ${OUTDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.fastq
	cat ${OUTDIR}/starter_${i}_1-${RATIO}.${COV1}x_2.fastq  ${HOMEDIR}/Triticum.chrM.${COV}x.${i}_1-${RATIO}_2.fastq  | seqkit shuffle --threads 8 --rand-seed 9000 --out-file  ${OUTDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.fastq
	fastqc --noextract --nogroup --format fastq ${OUTDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.fastq  ${OUTDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.fastq  -o ${OUTDIR}/0.fastqc_merge${COV}x
done
rm ${OUTDIR}/starter*

#!/bin/bash
#SBATCH -o %x-%J.20201107.log
#SBATCH -e %x-%J.20201107.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J trim_merge.sh
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add fastp/v0.9.1-424900e
module add fastqc/v0.11.9

COV=$1
RATIO=$2

mkdir /home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x
mkdir /home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x/0.fastqc_trim${COV}x

for i in 1 2 3 4 5
do 
  fastp -i ${HOMEDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.fastq  -o ${OUTDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.trim.fastq  -I ${HOMEDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.fastq -O ${OUTDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.trim.fastq 
  fastqc --noextract --nogroup --format fastq ${OUTDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.trim.fastq ${OUTDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.trim.fastq  -o ${OUTDIR}/0.fastqc_trim${COV}x
done 

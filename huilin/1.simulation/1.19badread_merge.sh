#!/bin/bash
#SBATCH -o %x-%J.20210119.log
#SBATCH -e %x-%J.20210119.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J badread_merge.sh
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

mkdir /home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/1.pollution${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/1.pollution${COV}x


#such as:sbatch 10.28badread_merge.sh 500 10


for i in 1 2 3 4 5
do
	cat ${HOMEDIR}/Triticum.chrC_badread${COV1}x.1_1-${RATIO}.fastq ${HOMEDIR}/Triticum.gene_badread${COV1}x.1_1-${RATIO}.fastq >${OUTDIR}/starter_${i}_1-${RATIO}.${COV}x_1.fastq
	cat ${OUTDIR}/starter_${i}_1-${RATIO}.${COV}x_1.fastq  ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.fastq  | seqkit shuffle --threads 8 --rand-seed 9000 --out-file  ${OUTDIR}/Triticum.merge.${i}_${COV}x_1-${RATIO}.fastq
done
rm ${OUTDIR}/starter*
rm ${HOMEDIR}/*gz


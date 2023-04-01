#!/bin/bash
#SBATCH -o %x-%J.20201109.log
#SBATCH -e %x-%J.20201109.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J flye_badread
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add flye/v2.8.1-3ee5b33

COV=$1
RATIO=$2
mkdir /home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/flye/badread${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/1.pollution${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/flye/badread${COV}x

for i in 1 2 3
do 
	mkdir ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}
	time flye --pacbio-raw ${HOMEDIR}/Triticum.merge.${i}_${COV}x_1-${RATIO}.fastq  --out-dir ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/
done

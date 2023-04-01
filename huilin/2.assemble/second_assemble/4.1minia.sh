#!/bin/bash
#SBATCH -o %x-%J.20210401.log
#SBATCH -e %x-%J.20210401.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J minia_mason
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1



module add minia/v3.2.4-bfd91a7

COV=$1
RATIO=$2
k=$3
mkdir /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/minia/mason${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/minia/mason${COV}x

for i in 1 2 3 4 5
do
	#kmer41
	mkdir ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}_k${k}
	time minia -kmer-size ${k}  -in ${HOMEDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.trim.fastq -in ${HOMEDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.trim.fastq  -out ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}_k${k}/Triticum.sim_${i}_1-${RATIO}_k${k}
done


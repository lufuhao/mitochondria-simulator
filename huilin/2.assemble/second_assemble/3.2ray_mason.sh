#!/bin/bash
#SBATCH -o %x-%J.20210302.log
#SBATCH -e %x-%J.20210302.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J ray_mason
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add ray/v2.3.1-rc1-621dfe84

COV=$1
RATIO=$2
k=$3

mkdir /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/ray/mason${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/ray/mason${COV}x

for i in 1 2 3 4 5
do
	#kmer41
	time mpiexec Ray -k ${k} -p ${HOMEDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.trim.fastq  ${HOMEDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.trim.fastq -o ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}_k${k}
done


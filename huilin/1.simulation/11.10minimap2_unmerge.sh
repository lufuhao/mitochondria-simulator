#!/bin/bash
#SBATCH -o %x-%J.20201118.log
#SBATCH -e %x-%J.20201118.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J minimap2_unmerge.sh
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add minimap2/v2.17-c9874e2
module add samtools/v0.1.20

COV=$1
RATIO=$2
mkdir /home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/2.minimap${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/2.minimap${COV}x

for i in 1 2 3
do
	minimap2 -ax map-pb /home/hpcusers/master/huilin_hu/simulate/Triticum.badread.15000bp.fa ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.fastq >${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.sam
	samtools view -Sb ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.sam > ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.bam
	samtools sort ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.bam ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.sort
done


#!/bin/bash
#SBATCH -o %x-%J.20201128.log
#SBATCH -e %x-%J.20201128.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J mummer_rmdup
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add mummer4/v4.0.0rc1

COV=$1
RATIO=$2
BIOSOFT=$3
CONTINGS=$4
HOMEDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/badread${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/badread${COV}x

mkdir ${OUTDIR}/mummer_rmdup

for i in 1 2 3
do
nucmer --maxmatch -c 65 -g 90 -p  ${OUTDIR}/mummer_rmdup/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}  /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa  ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/${CONTINGS}
mummerplot -f -l -p ${OUTDIR}/mummer_rmdup/Triticum.chrM_badread${COV}x.${i}_1-${RATIO} -s large  --png ${OUTDIR}/mummer_rmdup/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.delta
show-coords -T -r -l -c  ${OUTDIR}/mummer_rmdup/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.delta  > ${OUTDIR}/mummer_rmdup/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.coord
done

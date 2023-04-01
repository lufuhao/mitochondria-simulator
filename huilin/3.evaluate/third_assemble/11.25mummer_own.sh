#!/bin/bash
#SBATCH -o %x-%J.20201125.log
#SBATCH -e %x-%J.20201125.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J mummer_own
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

mkdir ${OUTDIR}/mummer_own

for i in 1 2 3 4 5
do 
nucmer --maxmatch -c 65 -g 90 -p  ${OUTDIR}/mummer_own/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/${CONTINGS}  ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/${CONTINGS}
mummerplot -p ${OUTDIR}/mummer_own/Triticum.chrM_badread${COV}x.${i}_1-${RATIO} -s large  --png ${OUTDIR}/mummer_own/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.delta
show-coords -T -r -l -c  ${OUTDIR}/mummer_own/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.delta > ${OUTDIR}/mummer_own/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.coord
done

###去重复片段
####得到${OUTDIR}/mummer_own/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.coord，找到最长的片段，坐标位置减1，得到去除重复片段的fasta
samtools faidx contigs.fasta tig00000001:1-451745 > conting.rmdup.fa
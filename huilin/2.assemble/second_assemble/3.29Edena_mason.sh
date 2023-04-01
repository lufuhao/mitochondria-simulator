#!/bin/bash
#SBATCH -o %x-%J.20210329.log
#SBATCH -e %x-%J.20210329.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J Edena_mason
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add edena/v3.131028

## 50x-500x:minoverlap:100 1000x:120 1500x:130 2000x:140

COV=$1
RATIO=$2
minOverlap=$3
mkdir /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/Edena/mason${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/Edena/mason${COV}x


for i in 1 2 3 4 5
do
        mkdir ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}
        time edena -nThreads 24 -minOverlap ${minOverlap} -DRpairs ${HOMEDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.fastq  ${HOMEDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.fastq  -p  ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}/
        time edena -e  ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}/.ovl -p ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}/
done
~      

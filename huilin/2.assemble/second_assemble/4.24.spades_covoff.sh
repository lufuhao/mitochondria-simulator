#!/bin/bash
#SBATCH -o %x-%J.20210424.log
#SBATCH -e %x-%J.20210424.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J spades_cov
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add  spades/v3.14.1

COV=$1
RATIO=$2
COV2=$(echo "scale=0;$COV/2" | bc -l)

mkdir /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/12.18spades/mason${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x/1.pollution${COV}x/1.trim${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/12.18spades/mason${COV}x

for i in 1 2 3 4 5
do
        #multiple kmer
        mkdir ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}
        time spades.py -m 32 --careful --cov-cutoff ${COV2} -t 16 -k 81,91,101,111   -1 ${HOMEDIR}/Triticum.merge_${i}_1.${COV}x.150bp.1-${RATIO}.trim.fastq  -2  ${HOMEDIR}/Triticum.merge_${i}_2.${COV}x.150bp.1-${RATIO}.trim.fastq  -o ${OUTDIR}/Triticum.sim_${i}_1-${RATIO}
done


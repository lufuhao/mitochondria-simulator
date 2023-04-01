#!/bin/bash
#SBATCH -o %x-%J.20201031.log
#SBATCH -e %x-%J.20201031.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J Nextdenovo_badread
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


COV=$1
RATIO=$2
RunDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/Nextdenovo/badread${COV}x/

#sbatch 10.31next_denovo_badread.sh 500 10
for i in 1 2 3
do
time /home/hpcusers/master/huilin_hu/biosoft/NextDenovo/nextDenovo  ${RunDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/run.cfg
done



#/home/hpcusers/master/huilin_hu/biosoft/NextDenovo  /home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/Nextdenovo/badread500x//Triticum.chrM_badread500x.1_1-10/run.cfg

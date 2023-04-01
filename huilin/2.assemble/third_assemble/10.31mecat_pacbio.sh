#!/bin/bash
#SBATCH -o %x-%J.20201109.log
#SBATCH -e %x-%J.20201109.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J Mecat2_badread
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add MECAT2/v20190314-f54c542


COV=$1
RATIO=$2
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/1.pollution${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/Mecat2/badread${COV}x



for i in 1 2 3
do
  time mecat.pl correct ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}_config_file.txt
  time mecat.pl trim ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}_config_file.txt
  time mecat.pl assemble ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}_config_file.txt
done

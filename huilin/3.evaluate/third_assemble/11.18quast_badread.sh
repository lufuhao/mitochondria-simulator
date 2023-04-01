#!/bin/bash
#SBATCH -o %x-%J.20201118.log
#SBATCH -e %x-%J.20201118.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J quast_badread
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1

###########注意，运行命令的时候加上倍数和运行的哪个软件以及k值##########
###########sbatch 10.14quast_badread 500 spades contings.fasta#######

module add quast/v5.1.0rc1-e010ca4


COV=$1
BIOSOFT=$2
CONTINGS=$3

mkdir /home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/badread${COV}x/quast
HOMEDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/badread${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/badread${COV}x/quast


for i in 1 2 3 4 5
do
	mkdir ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-10
	quast.py -r /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/${CONTINGS} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-10
done

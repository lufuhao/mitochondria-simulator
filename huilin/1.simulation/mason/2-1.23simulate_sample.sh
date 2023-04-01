#!/bin/bash
#SBATCH -o %x-%J.20210123.log
#SBATCH -e %x-%J.20210123.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J mason_sample.sh
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1

# Mitochondria 500X: readnumber=(452528*500)/(2*150) = 754210

module add  mason/v0.1.2
module add fastqtools/v0.8 fastqc/v0.11.9

COV=$1
GENESIZE1=$2
RATIO=$3
#Chinese spring 452528bp

#sbatch  1.19simulate_sample.sh 500 452528 10

COV1=$(echo "scale=0;$COV/20" |bc -l)
READNUM1=$(echo "scale=0;($GENESIZE1*$COV)/(2*150)" | bc -l)

###chrC和gene随机取样的readnum,为500x的一半再除以2
READNUM2=$(echo "scale=0;$READNUM1/($RATIO*2)" | bc -l)

mkdir /home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mt_ct_gn
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x 
for i in 1 2 
do
	##Mt
	fastq-sample -s 10 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM.${COV}x.1_1-${RATIO}_${i}  ${HOMEDIR}/Triticum.chrM.4000x_${i}.fastq
	fastq-sample -s 11 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM.${COV}x.2_1-${RATIO}_${i}  ${HOMEDIR}/Triticum.chrM.4000x_${i}.fastq
	fastq-sample -s 12 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM.${COV}x.3_1-${RATIO}_${i}  ${HOMEDIR}/Triticum.chrM.4000x_${i}.fastq
	fastq-sample -s 13 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM.${COV}x.4_1-${RATIO}_${i}  ${HOMEDIR}/Triticum.chrM.4000x_${i}.fastq
	fastq-sample -s 14 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM.${COV}x.5_1-${RATIO}_${i}  ${HOMEDIR}/Triticum.chrM.4000x_${i}.fastq
	
	fastq-sample -s 10 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC.${COV1}x.1_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.chrC.500x_${i}.fastq
	fastq-sample -s 11 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC.${COV1}x.2_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.chrC.500x_${i}.fastq
	fastq-sample -s 12 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC.${COV1}x.3_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.chrC.500x_${i}.fastq
	fastq-sample -s 13 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC.${COV1}x.4_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.chrC.500x_${i}.fastq
	fastq-sample -s 14 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC.${COV1}x.5_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.chrC.500x_${i}.fastq
	
	fastq-sample -s 10 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene.${COV1}x.1_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.gene.ct_500x_${i}.fastq
	fastq-sample -s 11 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene.${COV1}x.2_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.gene.ct_500x_${i}.fastq
	fastq-sample -s 12 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene.${COV1}x.3_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.gene.ct_500x_${i}.fastq
	fastq-sample -s 13 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene.${COV1}x.4_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.gene.ct_500x_${i}.fastq
	fastq-sample -s 14 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene.${COV1}x.5_1-${RATIO}_${i} ${HOMEDIR}/Triticum_cs.gene.ct_500x_${i}.fastq
done
mkdir ${OUTDIR}/0.fastqc${COV}x
for i in 1 2 3 4 5
do
	fastqc --noextract --nogroup --format fastq ${OUTDIR}/Triticum.chrM.${COV}x.${i}_1-${RATIO}_1.fastq  ${OUTDIR}/Triticum.chrM.${COV}x.${i}_1-${RATIO}_2.fastq -o ${OUTDIR}/0.fastqc${COV}x
	fastqc --noextract --nogroup --format fastq ${OUTDIR}/Triticum.chrC.${COV1}x.${i}_1-${RATIO}_1.fastq  ${OUTDIR}/Triticum.chrC.${COV1}x.${i}_1-${RATIO}_2.fastq -o ${OUTDIR}/0.fastqc${COV}x
	fastqc --noextract --nogroup --format fastq ${OUTDIR}/Triticum.gene.${COV1}x.${i}_1-${RATIO}_1.fastq  ${OUTDIR}/Triticum.gene.${COV1}x.${i}_1-${RATIO}_2.fastq -o ${OUTDIR}/0.fastqc${COV}x
done

#!/bin/bash
#SBATCH -o %x-%J.20210123.log
#SBATCH -e %x-%J.20210123.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J badread_sample.sh
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1

# Mitochondria 500X: readnumber=(452528*500)/(2*150) = 754210


module add fastqtools/v0.8 fastqc/v0.11.9

COV=$1
GENESIZE1=$2
RATIO=$3
#Chinese spring 452528bp

sbatch  1.19badread_sample.sh 50 452528 10

COV1=$(echo "scale=0;$COV/20" |bc -l)
READNUM1=$(echo "scale=0;($GENESIZE1*$COV)/(15000)" | bc -l)

###chrC和gene随机取样的readnum,为500x的一半再除以2
READNUM2=$(echo "scale=0;$READNUM1/($RATIO*2)" | bc -l)

mkdir /home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/mt_ct_gn
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x

gunzip -k ${HOMEDIR}/Triticum.chrC_badread500x.fastq.gz

#取六次重复，取三次相近的值
	#Mt
	fastq-sample -s 10 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.1_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 11 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.2_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 13 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.3_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 14 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.4_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 15 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.5_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 15 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.6_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	#ct
	fastq-sample -s 10 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.1_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 11 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.2_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 13 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.3_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 14 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.4_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 15 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.5_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 15 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.6_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	#gn
	fastq-sample -s 10 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.1_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 11 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.2_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 13 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.3_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 14 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.4_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 15 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.5_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 15 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.6_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq


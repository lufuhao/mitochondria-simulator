#!/bin/bash
#SBATCH -o %x-%J.20210313.log
#SBATCH -e %x-%J.20210313.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J minimap2_bedtools
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1

module add bedtools/v2.29.2-4ebba703 htslib/v1.10.2 samtools/v1.10
#sbatch 3.11spades_minimap2.sh 50 12.18spades  3.11count.tsv
COV=$1
BIOSOFT=$2
RESULT=$3
HOMEDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/badread${COV}x
for i in 1 2 3 4 5
do
	#echo在单引号里加入单引号可以写入变量
	echo -ne ''$BIOSOFT'\t'$COV'\t' >> /home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/${RESULT}
	minimap2 -a -N 5 -x asm5 --cs -c /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/contigs.ct.gn.filter.fa | samtools view -Sb | samtools sort > ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/biosoft_contigs_vs_ref.bam
	bedtools genomecov -ibam ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/biosoft_contigs_vs_ref.bam -bga | perl -ane '$l=$F[2]-$F[1];$t+=$l;$c+=$l if($F[3]>0);END{print "$t\t$c\t".($c/$t)."\t".($t-$c)}'>> /home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/${RESULT}
	minimap2 -a -N 5 -x asm5 --cs -c ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/contigs.ct.gn.filter.fa /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa | samtools view -Sb | samtools sort >${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/ref_vs_biosoft_contigs.bam
	bedtools genomecov -ibam ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/ref_vs_biosoft_contigs.bam -bga | perl -ane '$l=$F[2]-$F[1];$t+=$l;$c+=$l if($F[3]>0);END{print "\t$t\t$c\t".($c/$t)."\t"}'>>/home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/${RESULT}
	grep -c  "^>" ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-10/contigs.ct.gn.filter.fa >> /home/hpcusers/master/huilin_hu/Mt_assemble/Third_biosoft/${BIOSOFT}/${RESULT}
done

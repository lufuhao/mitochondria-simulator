#!/bin/bash
#SBATCH -o %x-%J.20210311.log
#SBATCH -e %x-%J.20210311.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J minimap2_bedtools
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1

#统计score脚本
module add bedtools/v2.29.2-4ebba703 htslib/v1.10.2 samtools/v1.10
#sbatch 3.10minimap2.sh 50 abyss 41 3.10count.tsv
COV=$1
BIOSOFT=$2
k=$3
#CONTIGS=$4
RESULT=$4

HOMEDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/${BIOSOFT}/mason${COV}x
#echo -e 'biosoft\tcov\tkmer\tref_tot\tref_cov\tref_cov_frac\t\tsnp\tqry_tot\tqry_cov\tqry_cov_frac\tcontig_num' >> /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/${BIOSOFT}/${RESULT}
for i in 1 2 3 4 5
do
	#echo在单引号里加入单引号可以写入变量
	echo -ne ''$BIOSOFT'\t'$COV'\t'$k'\t' >> /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/${BIOSOFT}/${RESULT}
	minimap2 -a -N 5 -x asm5 --cs -c /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa ${HOMEDIR}/Triticum.sim_${i}_1-10_k${k}/scaffold.samtools.fa | samtools view -Sb | samtools sort > ${HOMEDIR}/Triticum.sim_${i}_1-10_k${k}/biosoft_contigs_vs_ref.bam
	bedtools genomecov -ibam ${HOMEDIR}/Triticum.sim_${i}_1-10_k${k}/biosoft_contigs_vs_ref.bam -bga | perl -ane '$l=$F[2]-$F[1];$t+=$l;$c+=$l if($F[3]>0);END{print "$t\t$c\t".($c/$t)}'>> /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/${BIOSOFT}/${RESULT}
	minimap2 -a -N 5 -x asm5 --cs -c ${HOMEDIR}/Triticum.sim_${i}_1-10_k${k}/scaffold.samtools.fa /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa | samtools view -Sb | samtools sort >${HOMEDIR}/Triticum.sim_${i}_1-10_k${k}/ref_vs_biosoft_contigs.bam
	bedtools genomecov -ibam ${HOMEDIR}/Triticum.sim_${i}_1-10_k${k}/ref_vs_biosoft_contigs.bam -bga | perl -ane '$l=$F[2]-$F[1];$t+=$l;$c+=$l if($F[3]>0);END{print "\t$t\t$c\t".($c/$t)."\t"}'>>/home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/${BIOSOFT}/${RESULT}
	grep -c  "^>" ${HOMEDIR}/Triticum.sim_${i}_1-10_k${k}/scaffold.samtools.fa >> /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/${BIOSOFT}/${RESULT}
done

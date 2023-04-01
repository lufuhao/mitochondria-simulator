#!/bin/bash
#SBATCH -o %x-%J.20210326.log
#SBATCH -e %x-%J.20210326.err
#SBATCH -c 1
#SBATCH --ntasks 1
#SBATCH -J samtools_depth
#SBATCH -p fat
#SBATCH --mail-user=2721356963@qq.com
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=30000
#SBATCH -N 1


module add bwa/v0.7.17-13b5637 htslib/v1.10.2 samtools/v1.10 seqtk/v1.3-r106-c91a6af quast/v5.1.0rc1-e010ca4
#sbatch 3.26samtools_depth.sh 50 ray 25 Scaffolds.fasta 41

#samtools版本1.10才不会出错
COV=$1
BIOSOFT=$2
#depth=$(echo "scale=0;$COV/2" |bc -l);export depth
depth=$3;export depth
CONTIGS=$4
k=$5
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Second_simulate/mason${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/${BIOSOFT}/mason${COV}x
mkdir ${OUTDIR}/mummer_depth

for i in 1 2 3 4 5
do
#ct比对:
mkdir ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index
bwa index ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/${CONTIGS}
mv ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/${CONTIGS}.*  ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/
#Mt比对
time bwa mem -t 6  ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/${CONTIGS}  ${HOMEDIR}/Triticum.chrM.${COV}x.${i}_1-10_1.fastq.gz ${HOMEDIR}/Triticum.chrM.${COV}x.${i}_1-10_2.fastq.gz | samtools view -Sb | samtools sort > ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.scaffold.mt.sort.bam
#samtools得到深度列表
samtools depth ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.scaffold.mt.sort.bam > ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.mt.coverage
#python脚本求平均深度
python3 /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/py/3.25depth_max_avg.py ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.mt.coverage ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.mt.avgdepth
Rscript /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/py/3.30depth.R ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.mt.avgdepth ${OUTDIR}/mummer_depth/Triticum.sim_${i}_1-10_k${k}.scaffold.tiff
#Rscript /home/hpcusers/master/huilin_hu/Mt_assemble/Second_biosoft/py/4.14meddepth.R ${HOMEDIR}/Triticum.sim_${i}_1-10__k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.mt.avgdepth ${HOMEDIR}/mummer_depth/Triticum.sim_${i}_1-10_k${k}.med.depth.tiff
#取大于低depth的reads,并进行mummer比对和quast得到contigN50等
cat ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/bwa_index/Triticum.sim_${i}_1-10_k${k}.mt.avgdepth | awk -F"\t" '$3> ENVIRON["depth"]{print $1}' >${OUTDIR}/mummer_depth/Triticum.sim_${i}_1-10_k${k}.mt.list
seqtk subseq ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/${CONTIGS} ${OUTDIR}/mummer_depth/Triticum.sim_${i}_1-10_k${k}.mt.list > ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/scaffold.samtools.fa
nucmer --maxmatch -c 65 -g 90 -p ${OUTDIR}/mummer_depth/Triticum.sim_${i}_1-10_k${k}.scaffold.fiter  '/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa' ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/scaffold.samtools.fa
mummerplot -f -l -p ${OUTDIR}/mummer_depth/Triticum.sim_${i}_1-10_k${k}.scaffold.fiter --large  --png ${OUTDIR}/mummer_depth/Triticum.sim_${i}_1-10_k${k}.scaffold.fiter.delta
quast.py -r /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.old.fa  ${OUTDIR}/Triticum.sim_${i}_1-10_k${k}/${CONTIGS} -o ${OUTDIR}/quast_ct.gn/quast${k}/Triticum.sim_${i}_1-10_k${k}
done


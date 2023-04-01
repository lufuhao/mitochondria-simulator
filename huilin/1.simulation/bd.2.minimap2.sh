module add minimap2/v2.17-c9874e2
module add samtools/v0.1.20

COV=$1
RATIO=$2
mkdir /home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/2.minimap${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/2.minimap${COV}x

for i in 1 2 3
do
	minimap2 -ax map-pb /home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.badread.15000bp.fa ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.fastq.gz >${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.sam
	samtools view -Sb ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.sam > ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.bam
	samtools sort ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.bam ${OUTDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.sort
done


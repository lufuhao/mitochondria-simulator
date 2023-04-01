

module  add  badread/v0.1.5-9e030e8
####模拟4000x的线粒体和500x的叶绿体和1x的基因组



HOMEDIR=/home/hpcusers/master/huilin_hu/simulate
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/mt_ct_gn


###Mitochondria,$COV,15000bp
        badread simulate --reference  ${HOMEDIR}/Triticum_cs.chrM.badread.15000bp.fa --quantity 4000x \
                --error_model pacbio --qscore_model pacbio --identity 85,95,3  \
                | gzip > ${OUTDIR}/Triticum.chrM_badread4000x.fastq.gz
                

###Chloroplast,$COV1,15000bp
        badread simulate --reference  ${HOMEDIR}/Triticum_cs.chrC.fasta --quantity 500x \
                --error_model pacbio --qscore_model pacbio --identity 85,95,3  \
                | gzip > ${OUTDIR}/Triticum.chrC_badread500x.fastq.gz

#Gene
        badread simulate --reference  /home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta --quantity 1M \
                --error_model pacbio --qscore_model pacbio --identity 85,95,3  \
                >Triticum_genebadread.fastq




module add fastqtools/v0.8 fastqc/v0.11.9

COV=$1
GENESIZE1=$2
RATIO=$3
#Chinese spring 452528bp

#such as :sbatch 1.19badread_sample.sh 50 452528 10

COV1=$(echo "scale=0;$COV/20" |bc -l)
READNUM1=$(echo "scale=0;($GENESIZE1*$COV)/(15000)" | bc -l)

###chrC和gene随机取样的readnum,为500x的一半再除以2
READNUM2=$(echo "scale=0;$READNUM1/($RATIO*2)" | bc -l)

mkdir /home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/mt_ct_gn
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x

gunzip -k ${OUTDIR}/Triticum.chrM_badread4000x.fastq.gz
gunzip -k ${HOMEDIR}/Triticum.chrC_badread500x.fastq.gz

#取五次重复，然后去掉一个大的，一个小的
	#Mt
	fastq-sample -s 10 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.1_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 11 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.2_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 13 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.3_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 14 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.4_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	fastq-sample -s 15 -n ${READNUM1} -o ${OUTDIR}/Triticum.chrM_badread${COV}x.5_1-${RATIO}  ${HOMEDIR}/Triticum.chrM_badread4000x.fastq
	#ct
	fastq-sample -s 10 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.1_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 11 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.2_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 13 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.3_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 14 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.4_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	fastq-sample -s 15 -n ${READNUM2} -o ${OUTDIR}/Triticum.chrC_badread${COV1}x.5_1-${RATIO} ${HOMEDIR}/Triticum.chrC_badread500x.fastq
	#gn
	fastq-sample -s 10 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.1_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 11 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.2_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 13 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.3_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 14 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.4_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	fastq-sample -s 15 -n ${READNUM2} -o ${OUTDIR}/Triticum.gene_badread${COV1}x.5_1-${RATIO} ${HOMEDIR}/Triticum_badread1x.fastq
	gzip -9 ${OUTDIR}/


module add seqkit/v0.13.2  
module add fastqc/v0.11.9

COV=$1
RATIO=$2
COV1=$(echo "scale=0;$COV/20" |bc -l)

mkdir /home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/1.pollution${COV}x
HOMEDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x
OUTDIR=/home/hpcusers/master/huilin_hu/simulate/Third_simulate/badread${COV}x/1.pollution${COV}x


#such as:sbatch 10.28badread_merge.sh 500 10

for i in 1 2 3 4 5 
do
	zcat ${HOMEDIR}/Triticum.chrC_badread${COV1}x.1_1-${RATIO}.fastq.gz ${HOMEDIR}/Triticum.gene_badread${COV1}x.1_1-${RATIO}.fastq.gz >${OUTDIR}/starter_${i}_1-${RATIO}.${COV}x_1.fastq.gz
	zcat ${OUTDIR}/starter_${i}_1-${RATIO}.${COV}x_1.fastq.gz  ${HOMEDIR}/Triticum.chrM_badread${COV}x.${i}_1-${RATIO}.fastq.gz  | seqkit shuffle --threads 8 --rand-seed 9000 --out-file  ${OUTDIR}/Triticum.merge.${i}_${COV}x_1-${RATIO}.fastq.gz
done
rm ${OUTDIR}/starter*








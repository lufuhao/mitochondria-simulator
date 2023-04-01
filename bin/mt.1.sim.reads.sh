#!/bin/bash
### Exit if command fails
#set -o errexit
### Set readonly variable
#readonly passwd_file=”/etc/passwd”
### exit when variable undefined
#set -o nounset
### Script Root
RootDir=$(cd `dirname $(readlink -f $0)`; pwd)
### MachType
if [ ! -z $(uname -m) ]; then
	machtype=$(uname -m)
elif [ ! -z "$MACHTYPE" ]; then
	machtype=$MACHTYPE
else
	echo "Warnings: unknown MACHTYPE" >&2
fi

#export NUM_THREADS=`grep -c '^processor' /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`;
ProgramName=${0##*/}
echo "MachType: $machtype"
echo "RootPath: $RootDir"
echo "ProgName: $ProgramName"
RunPath=$PWD
echo "RunDir: $RunPath"

################# help message ######################################
help() {
cat<<HELP

$0 --- Simulate Illumina Paired-End reads using Mason

Version: v20210828

Requirements:
    Linux: gzip, mkdir, cd, perl
    1. 
    2. mason/badread, samtools
    3. fastqtools
    4. fastqc, seqkit, fastq_checkid.pl[FuhaoBin]
    5. fastp, fastqc
    6. BWA/minimap2, samtools, deepTools
    7. abyss,edena,ray+mpiexec,minia,spades
       canu,flye,mecat,nextdenovo,smartdenovo+make

Descriptions:
    1. Calculate the number of read pairs required based on the DNAsize and coverage
    2. Separately simulate mtDNA/ctDNA/gDNA paired reads using mason / badread
    3. Sample required number of reads from Step2
    4. Merge mtDNA + ctDNA + gDNA reads to output
    5. Trim fastq using Fastp
    6. Mapping simulated read back to mitochondria
    7. Assemble fastq to contigs/scaffolds
    8. To be defined

Options:
  -h    Print this help message
Fasta
  -mt  <FILE>  Mitochondria fasta; TIPS: for circular Mt, copy a read length of sequence from start to the end
  -gn  <FILE>  Genome fasta
  -ct  <FILE>  Chloroplast fasta

Simulation
  -pf  <STR>   Sequencing platform: illumina, pacbio
  -rl  <INT>   Read length, default: Illumina PE: 150, PacBio: 15000,13000
  -cp  <FLOAT> contamination proportion for gDNA and ctDNA, default: 0.05,0.05
  -mc  <INT>   Expected read depth for mitochrondria, multiple numbers delimited by comma, default: 1000,500
  -rp  <INT>   Replicates required, default: 5
  -fl  <INT>   Fragment length for Illumina, default: 350

Quality Control


Trimming [Illumina only]
  -tl  <INT>   Fastp: reads shorter than length_required will be discarded, default is 15

Assemble
  -as  <STR>   Assembler list, default: 
                   [For Illumina]: abyss,edena,ray,minia,spades
                   [for PacBio]  : canu,flye,mecat,nextdenovo,smartdenovo
  -km  <STR>   Kmer list for Illumina ONLY (abyss, ray, minia, spades)
                   *SPAdes and Edena are multi-Kmer assemblers
SLURM options
  -sa  <STR>   Run assembly manually by generated commands if you want to
                  evaluate the CPU and memory usage
  -em          Your Email address which slurm send messege [optional]

General
  -t    Number of threads, default: 1
  -s    Run Steps, default: 1,2,3,4,5,6,7
  -d    Output directory/path
  -k    Keep current existing file
  -mm <INT>    Max memory allowed in Gb, default 20
                   * SPAdes will terminates if exceeded
  -dg   Debug mode: Print more running details


Example:
  $0 -pf illumina -rl 150 -cp 0.05,0.05 -mc 2000,1000,500,250,100,50 -rp 5 -fl 350 -km 71,91,101 \
     -tl 50 -d \$PWD \
     -sa -em 'lufuhao@henu.edu.cn'
     -gn "/home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta" \
     -mt "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.mason.500bp.fa" \
     -ct "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrC.fasta"

  $0 -pf pacbio -cp 0.05,0.05 -mc 2000,1000,500,250,100,50 -rp 5 \
     -d \$PWD \
     -gn "/home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta" \
     -mt "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.mason.15000bp.fa" \
     -ct "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrC.fasta"

Author:

  Developer: Hui-Lin Hu (胡惠林)
  Master Students
  Grade 2019
  
  In the Laboratory of
  
  Fu-Hao Lu
  Professor, PhD
  State Key Labortory of Crop Stress Adaptation and Improvement
  College of Life Science
  Jinming Campus, Henan University
  Kaifeng 475004, P.R.China
  E-mail: lufuhao@henu.edu.cn
HELP
exit 2
}
[ $# -lt 1 ] && help
[ "$1" = "-h" ] || [ "$1" = "--help" ] && help
#################### Environments ###################################
echo -e "\n######################\nProgram $ProgramName initializing ...\n######################\n"


#################### Initializing ###################################
verbose=1
### genome
opt_gn=""
opt_mt=""
opt_ct=""
### read length
opt_fl=350
opt_rl=150
opt_rd=0
masonOptions=" -aNg -sq -ll 500 -le 200 -rn 2 -mp "
opt_rp=5
### trim
opt_tl=15
### mapping
binsize=10
### assembler
opt_as=""
opt_km=""
opt_sa=0;
opt_em=""
### 
###
declare -a opt_contamination_proportion=()
opt_t=1
outdir=$PWD
declare -a opt_mc=()
declare -a opt_s=(1 2 3 4 5 6 7)
keep_files=0
opt_mm=20
opt_debug=0
opt_genome_size=500000



#################### Parameters #####################################
while [ -n "$1" ]; do
  case "$1" in
    -h) help;shift 1;;
    -gn) opt_gn=$2;shift 2;;
    -mt) opt_mt=$2;shift 2;;
    -ct) opt_ct=$2;shift 2;;
    -rl) opt_rl=$2;shift 2;;
    -cp) opt_contamination_proportion=($(echo $2 | tr ',' "\n"));shift 2;;
    -mc) opt_mc=($(echo $2 | tr ',' "\n"));shift 2;;
    -rp) opt_rp=$2;shift 2;;
    -fl) opt_fl=$2;shift 2;;
    -tl) opt_tl=$2;shift 2;;
    -pf) opt_pf=$2;shift 2;;
    -as) opt_as=$2;shift 2;;
    -km) opt_km=$2;shift 2;;
    -sa) opt_sa=1;shift;;
    -em) opt_em=$2;shift 2;;
    -t) opt_t=$2;shift 2;;
    -s) opt_s=($(echo $2 | tr ',' "\n"));shift 2;;
    -dg) opt_debug=1; shift;;
    -d) outdir=$2;shift 2;;
    -k) keep_files=1; shift 1;;
    -mm) opt_mm=$2; shift 2;;
    --) shift;break;;
    -*) echo "error: no such option $1. -h for help" > /dev/stderr;exit 1;;
    *) break;;
  esac
done
#    -i) FastQR1Arr=($(echo $2 | tr  "\n"));shift 2;;
#    -1) seq_rfn=(${seq_rfn[@]} "$2");shift 2;;

#################### Subfuctions ####################################
###Detect command existence
CmdExit () {
  if command -v $1 >/dev/null 2>&1; then
    return 0
  else
    echo "Error: Program $1 not found" >&2
    exit 100
  fi
}
headerSlurm () {
	local HSjobname=$1
	local HSoutfile=$2
	
	if [ ! -s "$RootDir/../aux/slurm_header" ]; then
		echo "Error: SLURM header not found" >&2
		exit 100
	fi
	
	cat $RootDir/../aux/slurm_header > $HSoutfile
	echo "#SBATCH -J $HSjobname" >> $HSoutfile
	if [ ! -z "$opt_em" ]; then
		echo "#SBATCH --mail-type=END,FAIL" >> $HSoutfile
		echo "#SBATCH --mail-user=$opt_em" >> $HSoutfile
	fi
	echo "#SBATCH -o %J.%x.running.log" >> $HSoutfile
	echo "#SBATCH -o %J.%x.running.err" >> $HSoutfile
	echo "#SBATCH --mem=$opt_mm" >> $HSoutfile
	echo -e "\n\n\n" >> $HSoutfile
	
	return 0
}
getFastaSize () {
	local FSfa=$1
	
	local FSfaSize=0;
	
	if [ ! -s "$FSfa" ]; then
		echo "Error: invalid fasta file: $FSfa" >&2
		exit 100
	fi
	
	if [[ "$FSfa" =~ \.[gG][zZ]$ ]]; then
		echo "Info: genome in gzipped files" >&2
		FSfaSize=$(zcat $FSfa | grep -v ^'>' | perl -lne 'BEGIN{$sum=0;} $sum+=length($_); END {print $sum;}')
	else
		FSfaSize=$(grep -v ^'>' $FSfa | perl -lne 'BEGIN{$sum=0;} $sum+=length($_); END {print $sum;}')
	fi
	echo "Info: fasta: $FSfa ; Size: $FSfaSize" >&2
	
	if [[ "$FSfaSize" =~ ^[0-9]+$ ]] && [ $FSfaSize -gt 0 ]; then
		echo $FSfaSize
	else
		echo "Error: invalid fasta size: $FSfaSize for fasta $FSfa" >&2
		exit 100
	fi
}
printSeperator() {
	echo -e "\n\n\n"
	echo -e "\n\n\n" >&2
	echo -e "##################################################"
	echo -e "##################################################" >&2
	return 0
}
fileGunzip () {
	local DGFfile_in=$1
	local DGFfile_out=$2
	
	gzip -d -c $DGFfile_in > $DGFfile_out
	if [ $? -ne 0 ] || [ ! -s $DGFfile_out ]; then
		echo "(Func:fileGunzip)Error: decompress FASTA file failed" >&2
		echo "(Func:fileGunzip)    CMD used: gzip -d -c $DGFfile_in > $DGFfile_out" >&2
		exit 100
	fi
	
	return 0
}
sam2Bam () {
	local SBbam_in=$1
	local SBbam_out=$2
	
	if [ ! -s "$SBbam_in" ]; then
		echo "Error: sam2Bam invalid SAM input: $SBbam_in" >&2
		exit 100
	fi
	samtools view -S -b -h $SBbam_in > $SBbam_out
	if [ $? -ne 0 ] || [ ! -s $SBbam_out ]; then
		echo "Error: sam2Bam runninf failed" >&2
		echo "    CMD used: samtools view -S -b -h $SBbam_in > $SBbam_out" >&2
		exit 100
	fi
	
	return 0
}
bamSort () {
	local BSbam_in=$1
	local BSbam_out=$2
	
	if [ ! -s "$BSbam_in" ]; then
		echo "Error: invalid BAM input: $BSbam_in" >&2
		exit 100
	fi
	samtools sort -o $BSbam_out $BSbam_in
	if [ $? -ne 0 ] || [ ! -s "$BSbam_out" ]; then
		echo "Error: samtools sort error" >&2
		echo "    CMD used: samtools sort -o $BSbam_out $BSbam_in" >&2
		exit 100
	fi
	
	return 0
}



### Global: $samtoolsVers
bamIndex () {
	local BIbam_in=$1
	
	if [ ! -s "$BIbam_in" ]; then
		echo "Error: invalid BAM input: $BIbam_in" >&2
		exit 100
	fi
	
	local BIlongest=$(samtools view -H $BIbam_in | perl -lne 'BEGIN{$max=0;}if(/^\@SQ/){s/^.*\tLN://;s/\t.*$//;if ($_=~/^\d+$/ and $_>$max){$max=$_;}}END{print $max;}')
	
	if [ $BIlongest -ge 536870912 ]; then
		if [ $samtoolsVers -eq 1 ]; then
			samtools index -c $BIbam_in
			if [ $? -ne 0 ] || [ ! -s "${BIbam_in}.csi" ]; then
				echo "Error: samtools -c index error" >&2
				echo "    CMD used: samtools index $BIbam_in" >&2
				exit 100
			fi
		else
			echo "Error: please use SAMtools v1+ for longer chromosomes" >&2
			exit 100
		fi
	else
		samtools index $BIbam_in
		if [ $? -ne 0 ] || [ ! -s "${BIbam_in}.bai" ]; then
			echo "Error: samtools index error" >&2
			echo "    CMD used: samtools index $BIbam_in" >&2
			exit 100
		fi
	fi
	
	return 0
}



###Global: samtools
bamFlagstat () {
	local BFbam_in=$1
	
	if [ ! -s "$BFbam_in" ]; then
		echo "Error: invalid BAM input: $BFbam_in" >&2
		exit 100
	fi
	samtools flagstat $BFbam_in > $BFbam_in.flagstat
	if [ $? -ne 0 ] || [ ! -s "$BFbam_in.flagstat" ]; then
		echo "Error: samtools flagstat error" >&2
		echo "    CMD used: samtools flagstat $BFbam_in > $BFbam_in.flagstat" >&2
		exit 100
	fi
	
	return 0
}
### global: deepTools
bam2Bigwig () {
	local BBbam_in=$1
	local BBbwout=$2
	local BBbinsize=$3
	
	if [ ! -s "$BBbam_in" ]; then
		echo "Error: invalid BAM input: $BBbam_in" >&2
		exit 100
	fi
	if [ ! -z "$BBbinsize" ] && [[ ! "$BBbinsize" =~ ^[0-9]+$ ]] && [ $BBbinsize -ge 0 ]; then
		echo "Warnings: bamCoverage reset binsize to 10" >&2
		BBbinsize=10
	fi
	
	bamCoverage -b $BBbam_in -o $BBbwout --binSize $BBbinsize
	if [ $? -ne 0 ] || [ ! -s "$BBbwout" ]; then
		echo "Error: bamCoverage error" >&2
		echo "    CMD used: bamCoverage -b $BBbam_in -o $BBbwout --binSize $BBbinsize" >&2
		exit 100
	fi
	
	return 0
}
### $readnum=$(fqReadNum xxx.fq[.gz])
### Global: zcat, cat, wc
fqReadNum () {
	local FRNfq=$1
	
	local FRNlines=0;
	
	if [[ "$FRNfq" =~ \.[gG][zZ]$ ]]; then
		FRNlines=$(zcat $FRNfq | wc -l)
	else
		FRNlines=$(cat $FRNfq | wc -l)
	fi
	
	if [ $FRNlines -eq 0 ] || [ $(($FRNlines%4)) -ne 0 ]; then
		echo "Error: invalid FASTQ: $FRNfq" >&2
		exit 100
	fi
	
	echo $(($FRNlines/4))
}
### 
faSeqNum () {
	local FSNfa=$1
	
	local FSNlines=0;
	
	if [[ "$FSNfa" =~ \.[gG][zZ]$ ]]; then
		FSNlines=$(zcat $FSNfa | grep '^>' | wc -l)
	else
		FSNlines=$(cat $FSNfa | grep '^>' | wc -l)
	fi
	
	if [ $FSNlines -eq 0 ]; then
		echo "Error: invalid or empty FASTA: $FSNfa" >&2
		exit 100
	fi
	
	echo $FSNlines
}
### Global: fqReadNum
fqVerifyNumber () {
	local -a FVNargs=("$@")
	
	local FVNexp_num=${FVNargs[0]};
	unset FVNargs[0];
	local FVNcal_num=0;
	
	if [ ${#FVNargs[@]} -eq 0 ]; then
		echo "Info: fqVerifyNumber: empty FastQ list" >&2
		exit 100
	fi
	
	for FVNindfq in ${FVNargs[@]}; do
		FVNcal_num=$(fqReadNum $FVNindfq)
		if [ $FVNexp_num -eq $FVNcal_num ]; then
			echo "Info: Verified read number: $FVNexp_num for FASTQ: $FVNindfq"
		else
			echo "Info: NOT verified read number: expected $FVNexp_num but get $FVNcal_num for FASTQ: $FVNindfq" >&2
			exit 100
		fi
	done
	
	return 0
}
fqGzip () {
	local -a FGarr=("$@")
	
	for FGindfq in ${FGarr[@]}; do
		echo "Info: compressing fastq: $FGindfq"
		gzip -9 -c $FGindfq > "$FGindfq.gz"
		if [ $? -ne 0 ] || [ ! -s "$FGindfq.gz" ]; then
			echo "Error: compress FastQ failed: $FGindfq" >&2
			echo "    CMD used: gzip -9 $FGindfq" >&2
			exit 100
		fi
	done
	
	return 0
}
### Global: $verbose
fqGzipDel () {
	local -a FGarr=("$@")
	
	for FGindfq in ${FGarr[@]}; do
		if [ $verbose -eq 1 ]; then
			echo "Info: compressing fastq: $FGindfq"
		fi
		gzip -9 $FGindfq
		if [ $? -ne 0 ] || [ ! -s "$FGindfq.gz" ]; then
			echo "Error: compress FastQ failed: $FGindfq" >&2
			echo "    CMD used: gzip -9 $FGindfq" >&2
			exit 100
		fi
	done
	
	return 0
}
### fqGunzip in.fq.gz [in.fq] [1/0]
### arg3: 1=keep in.fq.gz; 0= delete in.fq.gz
### global: $verbose
fqGunzip () {
	local FGin=$1
	local FGout=$2
	local FGkeep=$3
	
	if [ -z "$FGout" ]; then
		FGout=${FGin%.*}
		if [ "$FGin" == "$FGout" ]; then
			echo "Error: FastQ out needs to be specified for: $FGin" >&2
			exit 100
		fi
	fi
	if [ -z "$FGkeep" ]; then
		FGkeep=0
	fi
	if [ $verbose -eq 1 ]; then
		echo "Info: decompressing fastq: $FGin"
	fi
	if [ $FGkeep -eq 0 ]; then
		gzip -d $FGin
	elif [ $FGkeep -eq 1 ]; then
		gzip -dc $FGin > $FGout
		if [ $? -ne 0 ] || [ ! -s "$FGout" ]; then
			echo "Error: compress FastQ failed: $FGin" >&2
			echo "    CMD used: gzip -dc $FGin > $FGout" >&2
			exit 100
		fi
	else
		echo "Error: unknown KEEP code: $FGkeep" >&2
		exit 100
	fi
	
	return 0
}
fq2Fa () {
	local FFfq_in=$1
	local FFfa_out=$2
	
	if [ ! -s "$FFfq_in" ]; then
		echo "Error: fq2Fa: invalid fastq input: $FFfq_in" >&2
		exit 100
	fi
	
	if [[ "$FFfq_in" =~ \.[gG][zZ]$ ]]; then
		if [[ "$FFfa_out" =~ \.[gG][zZ]$ ]]; then
			if [ $verbose -eq 1 ]; then
				echo "Info: Converting gzipped FastQ to gzipped Fasta"
			fi
			zcat $FFfq_in | perl -lne 'BEGIN{$num=0;}$num++;if ($num%4==1) {s/^\@/>/;print;}elsif($num%4==2){print;}' | gzip -9 > $FFfa_out
		else
			if [ $verbose -eq 1 ]; then
				echo "Info: Converting gzipped FastQ to flat Fasta"
			fi
			zcat $FFfq_in | perl -lne 'BEGIN{$num=0;}$num++;if ($num%4==1) {s/^\@/>/;print;}elsif($num%4==2){print;}' > $FFfa_out
		fi
	else
		if [[ "$FFfa_out" =~ \.[gG][zZ]$ ]]; then
			if [ $verbose -eq 1 ]; then
				echo "Info: Converting flat FastQ to gzipped Fasta"
			fi
			cat $FFfq_in | perl -lne 'BEGIN{$num=0;}$num++;if ($num%4==1) {s/^\@/>/;print;}elsif($num%4==2){print;}' | gzip -9 > $FFfa_out
		else
			if [ $verbose -eq 1 ]; then
				echo "Info: Converting flat FastQ to flat Fasta"
			fi
			cat $FFfq_in | perl -lne 'BEGIN{$num=0;}$num++;if ($num%4==1) {s/^\@/>/;print;}elsif($num%4==2){print;}' > $FFfa_out
		fi
	fi
	
	local FFnum1=$(fqReadNum $FFfq_in)
	local FFnum2=$(faSeqNum $FFfa_out)
	if [ ! -s $FFfa_out ] || [ -z "$FFnum1" ] || [ -z "$FFnum2" ] || [ $FFnum1 -ne $FFnum2 ]; then
		echo "Error: unequal seq number" >&2
		echo "    FastQ: $FFnum1    $FFfq_in" >&2
		echo "    FastA: $FFnum2    $FFfa_out" >&2
		exit 100
	fi
	
	return 0
}
### Global: mason/0.1.2
mason2reads () {
	local MRfa=$1
	local MRnum=$2
	local MRout=$3
	local MRreadpfx=$4
	
	MRfastqR1="${MRout%.fastq}_1.fastq"
	MRfastqR2="${MRout%.fastq}_2.fastq"
	
	mason illumina $masonOptions -N $MRnum -n $opt_rl --read-name-prefix $MRreadpfx -o $MRout $MRfa > $MRout.mason.log 2>&1
	if [ $? -ne 0 ]; then
		echo "Error: mason illumina running error" >&2
		echo "CMD used: mason illumina $masonOptions -N $MRnum -n $opt_rl --read-name-prefix $MRreadpfx -o $MRout $MRfa > $MRout.mason.log 2>&1" >&2
		exit 100
	fi
	
	if [ -s "$MRout.sam" ]; then
		sam2Bam "$MRout.sam" "$MRout.bam"
	fi
	if [ -s "$MRout.bam" ]; then
		bamSort "$MRout.bam" "$MRout.sort.bam"
		rm "$MRout.bam" > /dev/null 2>&1
	fi
	if [ -s "$MRout.sort.bam" ]; then
		bamIndex "$MRout.sort.bam"
		bamFlagstat "$MRout.sort.bam"
	fi
	
	if [ -s "$MRfastqR1" ] && [ -s "$MRfastqR2" ]; then
		mv $MRfastqR1 "${MRout%.*}.R1.fastq";
		# gzip -9 "${MRout%.*}.R1.fastq"
		mv $MRfastqR2 "${MRout%.*}.R2.fastq";
		# gzip -9 "${MRout%.*}.R2.fastq"
		return 0
	else
		echo "Error: mason read simution output failed: $MRout" >&2
		exit 100
	fi
	
	fqVerifyNumber $MRnum "${MRout%.*}.R1.fastq" "${MRout%.*}.R2.fastq"
	
	return 0
}
badread2reads () {
	local BRfa=$1
	local BRnum=$2
	local BRout=$3
	
	BRbdOptions=" --error_model pacbio --qscore_model pacbio --identity 85,95,3 "
	
	if [ $opt_rd -gt 0 ]; then
		BRbdOptions="$BRbdOptions --length $opt_rl,$opt_rd "
	fi

#   --quantity QUANTITY                   Either an absolute value (e.g. 250M) or a relative depth (e.g. 25x)
#   250M indicate bases, NOT read number
	echo "badread simulate --reference  ${BRfa} --quantity $(($BRnum*$opt_rl)) $BRbdOptions > ${BRout} 2> ${BRout}.badread.error"
	badread simulate --reference  ${BRfa} --quantity $(($BRnum*$opt_rl)) $BRbdOptions > ${BRout} 2> ${BRout}.badread.error
	if [ $? -ne 0 ]; then
		echo "Error: badread simulate running error" >&2
		echo "CMD used: badread simulate --reference  ${BRfa} --quantity $(($BRnum*$opt_rl)) $BRbdOptions > ${BRout} 2> ${BRout}.badread.error" >&2
		exit 100
	fi
	if [ -e "${BRout}.badread.error" ]; then
		rm "${BRout}.badread.error" >/dev/null 2>&1
	fi
	BRnumfq=$(fqReadNum $BRout)
	echo "Info: expected $BRnum and get $BRnumfq in FastQ: $BRout"
#	fqVerifyNumber $BRnum $BRout
	
	return 0
}
fqSampler () {
	local FSfq_in=$1
	local FSfq_out=$2
	local FSfq_num=$3
	local FSrandom=$4
	
	local FStmp_out=${FSfq_out%.fastq.gz}
	
	if [ ! -s $FSfq_in ]; then
		echo "Error: fqSampler: Fastq not found: $FSfq_in" >&2
		exit 100
	fi
	if [ -e "${FSfq_out%.gz}" ]; then
		rm "${FSfq_out%.gz}" >/dev/null 2>&1
	fi
	
	if [ $verbose -eq 1 ]; then
		echo "fastq-sample -s $FSrandom -n $FSfq_num -o $FStmp_out $FSfq_in"
	fi
	fastq-sample -s $FSrandom -n $FSfq_num -o $FStmp_out $FSfq_in
	if [ $? -ne 0 ] || [ ! -s "${FSfq_out%.gz}" ]; then
		echo "Error: fastq-sample running error" >&2
		echo "    CMD used: fastq-sample -s $FSrandom -n $FSfq_num -o $FStmp_out $FSfq_in" >&2
		exit 100
	fi
	
	if [ $verbose -eq 1 ]; then
		echo "fqVerifyNumber $FSfq_num ${FSfq_out%.gz}"
	fi
	fqVerifyNumber $FSfq_num ${FSfq_out%.gz}
	if [ $verbose -eq 1 ]; then
		echo "fqGzipDel ${FSfq_out%.gz}"
	fi
	fqGzipDel ${FSfq_out%.gz}
	
	return 0
}
fqMerge () {
	local FMfq_out=$1;
	local FMfq_in1=$2;
	local FMfq_in2=$3;
	local FMfq_in3=$4;
	local FMrandseed=$5

	if [ ! -s "$FMfq_in1" ]; then
		echo "Error: fqMerge file 1 not found: $FMfq_in1" >&2
		exit 100
	fi
	if [ ! -s "$FMfq_in2" ]; then
		echo "Error: fqMerge file 2 not found: $FMfq_in2" >&2
		exit 100
	fi
	if [ ! -s "$FMfq_in3" ]; then
		echo "Error: fqMerge file 3 not found: $FMfq_in3" >&2
		exit 100
	fi
	if [ -z "$FMrandseed" ] || [[ ! "$FMrandseed" =~ ^[0-9]+$ ]]; then
		echo "Warnings: seqkit shuffle --rand-seed parameter is invalid; reset to 9000" >&2
		FMrandseed=9000
	fi
	echo "Info: merging $FMfq_in1 $FMfq_in2 $FMfq_in3 into $FMfq_out"
	echo "Info: merging $FMfq_in1 $FMfq_in2 $FMfq_in3 into $FMfq_out" >&2
	zcat $FMfq_in1 $FMfq_in2 $FMfq_in3 | seqkit shuffle --threads $opt_t --rand-seed $FMrandseed --out-file $FMfq_out
	if [ $? -ne 0 ] || [ ! -s $FMfq_out ]; then
		echo "Error: Merge FastQ failed: $FMfq_out" >&2
		exit 100
	fi
	
	local RMall=$(fqReadNum $FMfq_out)
	local RMnum1=$(fqReadNum $FMfq_in1)
	local RMnum2=$(fqReadNum $FMfq_in2)
	local RMnum3=$(fqReadNum $FMfq_in3)
	if [ $(($RMnum1+$RMnum2+$RMnum3)) -ne $RMall ] || [ $RMall -eq 0 ]; then
		echo "Error: read lines are not equal: expected $(($RMnum1+$RMnum2+$RMnum3)) but get $RMall"
		exit 100
	fi
	
	return 0
}
### Global: kmers=(41, 51, 61, 71)
runAbyss () {
	local RAr1=$1
	local RAr2=$2
	local RAdir=$3
	
	local RAslurmfile=""
	
	if [ -z "$RAdir" ];then
		RAdir=$PWD
	fi
	if [ ! -d $RAdir ]; then
		mkdir -p $RAdir
	fi
	
	for RAindk in ${kmers[@]}; do
		cd $RAdir
		if [ $keep_files -eq 1 ] && [ -d $RAdir/K$RAindk ] && [ -L "$RAdir/K$RAindk/${RAdir##*/}-contigs.fa" ]; then
			echo "Info: using existing ABySS out at K $RAindk: $RAdir/K$RAindk"
		else
			if [ -d "$RAdir/K$RAindk" ]; then
				rm -rf "$RAdir/K$RAindk" >/dev/null 2>&1
			fi
			RAslurmfile="$RAdir/${RAdir##*/}.K$RAindk.command"
			headerSlurm "${RAdir##*/}.K$RAindk" $RAslurmfile
			echo "mkdir -p $RAdir/K$RAindk" >> $RAslurmfile
			echo "cd $RAdir/K$RAindk" >> $RAslurmfile
			echo "time abyss-pe k=$RAindk name=${RAdir##*/} in=\"$RAr1 $RAr2\" > ${RAdir##*/}.abyss.log 2>&1" >> $RAslurmfile
			
			if [ $opt_sa -eq 0 ]; then
				mkdir -p $RAdir/K$RAindk
				cd $RAdir/K$RAindk
				echo "time abyss-pe k=$RAindk name=${RAdir##*/} in=\"$RAr1 $RAr2\" > ${RAdir##*/}.abyss.log 2>&1"
				time abyss-pe k=$RAindk name=${RAdir##*/} in="$RAr1 $RAr2" > ${RAdir##*/}.abyss.log 2>&1
				if [ $? -ne 0 ] || [ ! -L "$RAdir/K$RAindk/${RAdir##*/}-contigs.fa" ]; then
					echo "Error: abyss-pe running error" >&2
					echo "    CMD used: time abyss-pe k=$RAindk name=${RAdir##*/} in=\"$RAr1 $RAr2\" > ${RAdir##*/}.abyss.log 2>&1" >&2
					exit 100
				fi
			else
				echo "Info: SLURM file: $RAslurmfile"
			fi
		fi
	done
	
	return 0
}
# Edena program requires the reads to be all the same length
runEdena () {
	local REr1=$1
	local REr2=$2
	local REdir=$3
	local REmin_overlap=$4
	
	local REslurmfile=""
	
	if [ -z "$REdir" ];then
		REdir=$PWD
	fi
	if [ ! -d $REdir ]; then
		mkdir -p $REdir
	fi
	if [[ "$REr1" =~ \.[gG][zZ]$ ]]; then
		echo "Warnings: Edena Not support gzipped R1 fastq: $REr1" >&2
		if [ -s "${REr1%.gz}" ]; then
			echo "Info: Edena using existing R1 fastq: ${REr1%.gz}"
		else
			fqGunzip $REr1 ${REr1%.gz} 1
		fi
		REr1=${REr1%.gz}
	elif [[ "$REr1" =~ \.[fF][qQ]$ ]] || [[ "$REr1" =~ \.[fF][aA][sS][tT][qQ]$ ]]; then
		echo "Info: Edena use flat FASTQ R1: check"
	else
		echo "Error: Edena unknown R1 fastq format: $REr1" >&2
		exit 100
	fi
	if [[ "$REr2" =~ \.[gG][zZ]$ ]]; then
		echo "Warnings: Edena Not support gzipped R2 fastq: $REr2" >&2
		if [ -s "${REr2%.gz}" ]; then
			echo "Info: Edena using existing R2 fastq: ${REr2%.gz}"
		else
			fqGunzip $REr2 ${REr2%.gz} 1
		fi
		REr2=${REr2%.gz}
	elif [[ "$REr2" =~ \.[fF][qQ]$ ]] || [[ "$REr2" =~ \.[fF][aA][sS][tT][qQ]$ ]]; then
		echo "Info: Edena use flat FASTQ R2: check"
	else
		echo "Error: Edena unknown R2 fastq format: $REr2" >&2
		exit 100
	fi
	REslurmfile="$REdir/${REdir##*/}.command"
	headerSlurm "${REdir##*/}" $REslurmfile
	echo "cd $REdir" >> $REslurmfile
	echo "time edena -minOverlap ${REmin_overlap} -DRpairs $REr1  $REr2 -p ${REdir##*/} > ${REdir##*/}.EDENA.overlap.log 2>&1" >> $REslurmfile
	echo "time edena -e $REdir/${REdir##*/}.ovl -p ${REdir##*/} > ${REdir##*/}.EDENA.assembler.log 2>&1" >> $REslurmfile
	
	if [ $opt_sa -eq 0 ]; then
		cd $REdir
		echo "time edena -minOverlap ${REmin_overlap} -DRpairs $REr1  $REr2 -p ${REdir##*/} > ${REdir##*/}.EDENA.overlap.log 2>&1"
		time edena -minOverlap ${REmin_overlap} -DRpairs $REr1  $REr2 -p ${REdir##*/} > ${REdir##*/}.EDENA.overlap.log 2>&1
		if [ $? -ne 0 ] || [ ! -s "$REdir/${REdir##*/}.ovl" ]; then
			echo "Error: Edena overlapping mode failed" >&2
			echo "    CMD used: time edena -minOverlap ${REmin_overlap} -DRpairs $REr1  $REr2 -p ${REdir##*/} > ${REdir##*/}.EDENA.overlap.log 2>&1" >&2
			exit 100
		fi
		echo "time edena -e $REdir/${REdir##*/}.ovl -p ${REdir##*/} > ${REdir##*/}.EDENA.assembler.log 2>&1"
		time edena -e $REdir/${REdir##*/}.ovl -p ${REdir##*/} > ${REdir##*/}.EDENA.assembler.log 2>&1
		if [ $? -ne 0 ] || [ ! -s "$REdir/${REdir##*/}_contigs.fasta" ]; then
			echo "Error: Edena assembler mode failed" >&2
			echo "    CMD used: time edena -e $REdir/${REdir##*/}.ovl -p ${REdir##*/} > ${REdir##*/}.EDENA.assembler.log 2>&1" >&2
			exit 100
		fi
	else
		echo "Info: SLURM file: $REslurmfile"
	fi
	
	return 0
}
### Global: kmers=(41, 51, 61, 71)
### discards pairing information
runMinia () {
	local RMr1=$1
	local RMr2=$2
	local RMdir=$3
	
	local RMslurmfile=""
	
	if [ -z "$RMdir" ];then
		RMdir=$PWD
	fi
	if [ ! -d $RMdir ]; then
		mkdir -p $RMdir
	fi
	
	for RMindk in ${kmers[@]}; do
		cd $RMdir
		if [ $keep_files -eq 1 ] && [ -d $RRdir/K$RRindk ] && [ -s "$RRdir/K$RRindk/Contigs.fasta" ]; then
			echo "Info: using existing Ray out at K $RRindk: $RRdir/K$RRindk"
		else
			RMslurmfile="$RMdir/${RMdir##*/}.K$RMindk.command"
			headerSlurm "${RMdir##*/}.K$RMindk" "$RMslurmfile"
			echo "cd $RMdir/K$RMindk" >> "$RMslurmfile"
			echo "time minia -kmer-size ${RMindk}  -in $RMr1 -in $RMr2  -out ${RMdir##*/} > ${RMdir##*/}.minia.K$RMindk.log 2>&1" >> "$RMslurmfile"
			
			if [ $opt_sa -eq 0 ]; then
				mkdir $RMdir/K$RMindk
				cd $RMdir/K$RMindk
				echo "time minia -kmer-size ${RMindk}  -in $RMr1 -in $RMr2  -out ${RMdir##*/} > ${RMdir##*/}.minia.K$RMindk.log 2>&1"
				time minia -kmer-size ${RMindk}  -in $RMr1 -in $RMr2  -out ${RMdir##*/} > ${RMdir##*/}.minia.K$RMindk.log 2>&1
				if [ $? -ne 0 ] || [ ! -s "$RMdir/K$RMindk/${RMdir##*/}.contigs.fa" ]; then
					echo "Error: minia running error" >&2
					echo "    CMD used: time minia -kmer-size ${RMindk}  -in $RMr1 -in $RMr2  -out ${RMdir##*/} > ${RMdir##*/}.minia.K$RMindk.log 2>&1" >&2
					exit 100
				fi
			else
				echo "Info: SLURM file: $RMslurmfile"
			fi
		fi
	done
	
	return 0
}
### Global: kmers=(41, 51, 61, 71)
### Ray -o outdir will auto-create, do NOT create it before running
runRay () {
	local RRr1=$1
	local RRr2=$2
	local RRdir=$3
	
	local RRslurmfile=""
	
	if [ -z "$RRdir" ];then
		RRdir=$PWD
	fi
	if [ ! -d $RRdir ]; then
		mkdir -p $RRdir
	fi
	
	for RRindk in ${kmers[@]}; do
		cd $RRdir
		if [ $keep_files -eq 1 ] && [ -d $RRdir/K$RRindk ] && [ -s "$RRdir/K$RRindk/Contigs.fasta" ]; then
			echo "Info: using existing Ray out at K $RRindk: $RRdir/K$RRindk"
		else
			RRslurmfile="$RRdir/${RRdir##*/}.K$RRindk.command"
			headerSlurm "${RRdir##*/}.K$RRindk" "$RRslurmfile"
			echo "cd $RRdir/" >> "$RRslurmfile"
			echo "time mpiexec -n 1 Ray -k ${RRindk} -p $RRr1 $RRr2 -o $RRdir/K$RRindk" >> "$RRslurmfile"
			if [ $opt_sa -eq 0 ]; then
				echo "time mpiexec -n 1 Ray -k ${RRindk} -p $RRr1 $RRr2 -o $RRdir/K$RRindk"
#				time mpiexec Ray -k ${RRindk} -p $RRr1 $RRr2 -o $RRdir/K$RRindk
				time mpiexec -n 1 Ray -k ${RRindk} -p $RRr1 $RRr2 -o $RRdir/K$RRindk > $RRdir/${RRdir##*/}.Ray.K$RRindk.log 2>&1
				if [ $? -ne 0 ] || [ ! -s "$RRdir/K$RRindk/Contigs.fasta" ] || [ ! -s "$RRdir/K$RRindk/Scaffolds.fasta" ]; then
					echo "Error: Ray running error" >&2
					echo "    CMD used: time mpiexec -n 1 Ray -k ${RRindk} -p $RRr1 $RRr2 -o $RRdir/K$RRindk > $RRdir/${RRdir##*/}.Ray.K$RRindk.log 2>&1" >&2
					exit 100
				fi
			else
				echo "Info: SLURM file: $RRslurmfile"
			fi
		fi
	done
	
	return 0
}



###global: $opt_mm
runSpades () {
	local RSr1=$1
	local RSr2=$2
	local RSdir=$3
	
	local RSslurmfile="$RSdir/${RSdir##*/}.command"
	
	if [ -z "$RSdir" ];then
		RSdir=$PWD
	fi
	if [ ! -d $RSdir ]; then
		mkdir -p $RSdir
	fi
	
	headerSlurm "${RSdir##*/}" "$RSslurmfile"
	echo "cd $RSdir" >> "$RSslurmfile"
	echo "time spades.py -m $opt_mm --careful --cov-cutoff auto -k $opt_km -1 $RSr1 -2 $RSr2 -o $RSdir > ${RSdir##*/}.spades.log 2>&1" >> "$RSslurmfile"
# -m RAM limit for SPAdes in Gb (terminates if exceeded). [default: 250]
#--careful tries to reduce number of mismatches and short indels *****
#-t <int>, --threads <int>   number of threads. [default: 16]
#	time spades.py -m $opt_mm --careful --cov-cutoff auto -1 $RSr1 -2 $RSr2 -o $RSdir
	if [ $opt_sa -eq 0 ]; then
		cd $RSdir
		echo "time spades.py -m $opt_mm --careful --cov-cutoff auto -k $opt_km -1 $RSr1 -2 $RSr2 -o $RSdir > ${RSdir##*/}.spades.log 2>&1"
		time spades.py -m $opt_mm --careful --cov-cutoff auto -k $opt_km -1 $RSr1 -2 $RSr2 -o $RSdir > ${RSdir##*/}.spades.log 2>&1
		if [ $? -ne 0 ] || [ ! -s "$RSdir/contigs.fasta" ] || [ ! -s "$RSdir/scaffolds.fasta" ]; then
			echo "Error: SPAdes running failed" >&2
			echo "    CMD used: time spades.py -m $opt_mm --careful --cov-cutoff auto -k $opt_km -1 $RSr1 -2 $RSr2 -o $RSdir > ${RSdir##*/}.spades.log 2>&1" >&2
			exit 100
		fi
	else
		echo "Info: SLURM file: $RSslurmfile"
	fi
	
	return 0
}



### Global: $verbose
runCanu () {
	local RCfq=$1
	local RCdir=$2
	
	local RCslurmfile="$RCdir/${RCdir##*/}.command"
	
#	local RCspecfile="$RCdir/myspec.txt"
	
	if [ -z "$RCdir" ];then
		RCdir=$PWD
	fi
	if [ ! -d $RCdir ]; then
		mkdir -p $RCdir
	fi
	
	headerSlurm "${RCdir##*/}" "$RCslurmfile"
	echo "cd $RCdir" >> $RCslurmfile
	echo "time canu -p ${RCdir##*/} -d $RCdir genomeSize=$opt_genome_size  -pacbio $RCfq > ${RCdir##*/}.canu.log 2>&1" >> $RCslurmfile
	
	if [ $opt_sa -eq 0 ]; then
		cd $RCdir
#	echo "useGrid=0" > $RCspecfile
#	 -s $RCspecfile
		if [ $verbose -eq 1 ]; then
			echo "time canu -p ${RCdir##*/} -d $RCdir genomeSize=$opt_genome_size  -pacbio $RCfq"
		fi
		time canu -p ${RCdir##*/} -d $RCdir genomeSize=$opt_genome_size  -pacbio $RCfq > ${RCdir##*/}.canu.log 2>&1
		if [ $? -ne 0 ] && [ -s "$RCdir/${RCdir##*/}.contigs.fasta" ]; then
			echo "Error: CANU running failed" >&2
			echo "    CMD used: time canu -p ${RCdir##*/} -d $RCdir genomeSize=$opt_genome_size  -pacbio $RCfq > ${RCdir##*/}.canu.log 2>&1"
			exit 100
		fi
	else
		echo "Info: SLURM file: $RCslurmfile"
	fi
#out: "$RCdir/${RCdir##*/}.contigs.fasta"
#lufuhao  1259558 1259246 99 16:28 ?        00:02:09 /home/hpcsoft/ProdSoft/jdk/v11.0.8/x64/bin/java -XX:ParallelGCThreads=8 -server -Xms5530m -Xmx5530m -jar /home/hpcsoft/TestSoft/canu/v2.1/x86_64/bin/../share/java/classes/mhap-2.1.3.jar --repeat-weight 0.9 --repeat-idf-scale 10 -k 16 --store-full-id --num-hashes 256 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --threshold 0.8 --filter-threshold 0.0000001 --min-olap-length 500 --num-threads 8 -f ../../0-mercounts/ALLpacbio.mt500x.Rep1.canu.ms16.ignore.gz -p ./000001.input.fasta -q .
#lufuhao  1255059 1255030 99 15:56 ?        00:01:27 /home/hpcsoft/TestSoft/canu/v2.1/x86_64/bin/overlapInCore -t 8 -k 22 -k ../0-mercounts/ALLpacbio.mt500x.Rep1.canu.ms22.dump --hashbits 22 --hashload 0.8 --maxerate 0.045 --minlength 500 -h 1-16572 -r 1-16572 --hashdatalen 14263149 -o ./001/000001.ovb.WORKING -s ./001/000001.stats ../../ALLpacbio.mt500x.Rep1.canu.seqStore

	return 0
}
runFlye () {
	local RFfq=$1
	local RFdir=$2
	
	local RFslurmfile="$RFdir/${RFdir##*/}.command"
	
	if [ -z "$RFdir" ];then
		RFdir=$PWD
	fi
	if [ ! -d $RFdir ]; then
		mkdir -p $RFdir
	fi
	
	headerSlurm "${RFdir##*/}" "$RFslurmfile"
	echo "cd $RFdir" >> $RFslurmfile
	echo "time flye --pacbio-raw $RFfq  --out-dir $RFdir > ${RFdir##*/}.log 2>&1" >> $RFslurmfile
	
	if [ $opt_sa -eq 0 ]; then
		cd $RFdir
		if [ $verbose -eq 1 ]; then
			echo "time flye --pacbio-raw $RFfq  --out-dir $RFdir > ${RFdir##*/}.log 2>&1"
		fi
		time flye --pacbio-raw $RFfq  --out-dir $RFdir > ${RFdir##*/}.flye.log 2>&1
		if [ $? -ne 0 ] && [ -s "$RFdir/assembly.fasta" ]; then
			echo "Error: Flye running failed" >&2
			echo "    CMD used: time flye --pacbio-raw $RFfq  --out-dir $RFdir > ${RFdir##*/}.log 2>&1"
			exit 100
		else
			if [ -e "${RFdir##*/}.flye.log" ]; then
				rm "${RFdir##*/}.flye.log" >/dev/null 2>&1
			fi
		fi
	else
		echo "Info: SLURM file: $RFslurmfile"
	fi
	
	return 0
}



runMecat () {
	local RMfq=$1
	local RMdir=$2
	
	local RMcfg="$PWD/${RMdir##*/}.cfg"
	local RMslurmfile="$PWD/${RMdir##*/}.command"
	
	if [ -z "$RMdir" ];then
		RMdir=$PWD
	fi

cat > $RMcfg <<EOMECAT
PROJECT=${RMdir##*/}
RAWREADS=$RMfq
GENOME_SIZE=$opt_genome_size
MIN_READ_LENGTH=3000
CNS_OVLP_OPTIONS="-kmer_size 13"
CNS_PCAN_OPTIONS="-p 100000 -k 100"
CNS_OPTIONS=""
CNS_OUTPUT_COVERAGE=30
TRIM_OVLP_OPTIONS="-skip_overhang"
TRIM_PM4_OPTIONS="-p 100000 -k 100"
TRIM_LCR_OPTIONS=""
TRIM_SR_OPTIONS=""
ASM_OVLP_OPTIONS=""
FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1"
FSA_ASSEMBLE_OPTIONS=""
CLEANUP=0
EOMECAT
#THREADS=4
#PROJECT=ecoli, the name of the project. In this example, a directory ecoli will be created in the current directory, and then everything will take place in the directory ecoli.
#RAWREADS=, the raw reads (with full path) to be processed by MECAT2. See Input Format.
#GENOME_SIZE=, the size (in bp) of the underlying genome.
#THREADS=, number of CPU threads used by MECAT2.
#MIN_READ_LENGTH=, minimal length of corrected reads and trimmed reads.
#CNS_OVLP_OPTIONS="", options for detecting overlap candidates in the correction stage. Run mecat2map -help for details. Note that the output format is seqidx (-outfmt seqidx), which is set internally by mecat.pl.
#CNS_PCAN_OPTIONS: mecat2pcan -help: 
#    -p <Integer, >0>    Number of reads in each part,     Default = '100000'
#    -k <Integer, >0>    Number of parts processed at each round,     (If <0, then it will be set to system limit value),     Default = '100'
#CNS_OPTIONS="", options for correcting raw reads. Run mecat2cns -help for details.
#CNS_OUTPUT_COVERAGE=30, number of coverage of the longest corrected reads are extracted to be trimed and then assembled. In this example, 30x (specifically, 30 * 4800000 = 144 MB) of the longest corrected reads will be extracted.
#TRIM_OVLP_OPTIONS="", options for detecting overlaps in the trimming stage. Run mecat2map for details. Note that output format is m4x (-outfmt m4x), which is set internally by mecat.pl.
#TRIM_PM4_OPTIONS="-p 100000 -k 100"
#TRIM_LCR_OPTIONS=""
#TRIM_SR_OPTIONS=""
#ASM_OVLP_OPTIONS="", options for detecting overlaps in the assemble stage. Run mecat2map -help for details. The output format is m4 (-outfmt m4), which is set internally by mecat.pl.
#FSA_OL_FILTER_OPTIONS="", options for filtering overlaps. See below for details.
#FSA_ASSEMBLE_OPTIONS="", options for assembling trimmed reads. See below for details.
#USE_GRID=false, using multiple computing nodes (true) or not (false).
#CLEANUP=0, delete intermediate date genrated by MECAT2 (1) or not (0). Please note the in assemblying large genomes, the intermediate data can be very large.

	headerSlurm "${RMdir##*/}" "$RMslurmfile"
	echo "cd $PWD" >> $RMslurmfile
	echo "time mecat.pl correct $RMcfg > $PWD/${RMdir##*/}.correct.log  2>&1" >> $RMslurmfile
	echo "time mecat.pl trim $RMcfg > $PWD/${RMdir##*/}.trim.log  2>&1" >> $RMslurmfile
	echo "time mecat.pl assemble $RMcfg > $PWD/${RMdir##*/}.assemble.log  2>&1" >> $RMslurmfile
	if [ $opt_sa -eq 0 ]; then
		if [ $verbose -eq 1 ]; then
			echo "time mecat.pl correct $RMcfg > $PWD/${RMdir##*/}.correct.log  2>&1"
		fi
		time mecat.pl correct $RMcfg > $PWD/${RMdir##*/}.correct.log  2>&1
		if [ $? -ne 0 ]; then
			echo "Error: MECAT correct failed" >&2
			echo "    CMD: time mecat.pl correct $RMcfg > $PWD/${RMdir##*/}.correct.log  2>&1" >&2
			exit 100
		fi
		if [ $verbose -eq 1 ]; then
			echo "time mecat.pl trim $RMcfg > $PWD/${RMdir##*/}.trim.log  2>&1"
		fi
		time mecat.pl trim $RMcfg > $PWD/${RMdir##*/}.trim.log  2>&1
		if [ $? -ne 0 ]; then
			echo "Error: MECAT trim failed" >&2
			echo "    CMD: time mecat.pl trim $RMcfg > $PWD/${RMdir##*/}.trim.log  2>&1" >&2
			exit 100
		fi
		if [ $verbose -eq 1 ]; then
			echo "time mecat.pl assemble $RMcfg > $PWD/${RMdir##*/}.assemble.log  2>&1"
		fi
		
		time mecat.pl assemble $RMcfg > $PWD/${RMdir##*/}.assemble.log  2>&1
		if [ $? -ne 0 ] || [ ! -s "$RMdir/4-fsa/contigs.fasta" ]; then
			echo "Error: MECAT assemble failed" >&2
			echo "    CMD: time mecat.pl assemble $RMcfg > $PWD/${RMdir##*/}.assemble.log  2>&1" >&2
			exit 100
		fi
	else
		echo "Info: SLURM file: $RMslurmfile"
	fi

	return 0
}



runNextdenovo () {
	local RNfq=$1
	local RNdir=$2
	
	local RNcfg="$RNdir/${RNdir##*/}.setting.cfg"
	local RNslurmfile="$RNdir/${RNdir##*/}.command"
	
	if [ -z "$RNdir" ];then
		RNdir=$PWD
	fi
	if [ ! -d $RNdir ]; then
		mkdir -p $RNdir
	fi
	
	headerSlurm "${RNdir##*/}" "$RNslurmfile"
	echo "cd $RNdir" >> "$RNslurmfile"
	echo "time nextDenovo $RNcfg > $RNdir/${RNdir##*/}.log 2>&1" >> "$RNslurmfile"
	
	
	echo "$RNfq" > "$RNdir/run.fofn"

cat > $RNcfg <<EONEXTDENOVO
[General]
job_type = local
job_prefix = ${RNdir##*/}
task = all
rewrite = yes
deltmp = yes
rerun = 3
parallel_jobs = 1
input_type = raw
read_type = hifi
input_fofn = run.fofn
workdir = ${RNdir}

[correct_option]
read_cutoff = 1k
seed_cutoff = 3k
blocksize = 3g
pa_correction = 20
seed_cutfiles = 20
sort_options = -m ${opt_mm}g -k 40
minimap2_options_raw = -x ava-pb
correction_options = -p 8

[assemble_option]
random_round = 20
minimap2_options_cns = -x ava-pb -k17 -w17
nextgraph_options = -a 1
EONEXTDENOVO

#task = all # 'all', 'correct', 'assemble'
#rewrite = yes # yes/no
	if [ $opt_sa -eq 0 ]; then
		cd $RNdir
		if [ $verbose -eq 1 ]; then
			echo "time nextDenovo $RNcfg > $RNdir/${RNdir##*/}.log 2>&1"
		fi
		
		time nextDenovo $RNcfg > $RNdir/${RNdir##*/}.log 2>&1
		if [ $? -ne 0 ] || [ ! -d "$RNdir/03.ctg_graph" ] || [ ! -s "$RNdir/03.ctg_graph/nd.asm.fasta" ]; then
			echo "Error: nextDenovo running failed" >&2
			echo "    CMD: time nextDenovo $RNcfg > $RNdir/${RNdir##*/}.log 2>&1" >&2
			exit 100
		else
			if [ -e "$RNdir/${RNdir##*/}.log" ]; then
				rm "$RNdir/${RNdir##*/}.log" >/dev/null 2>&1
			fi
		fi
	else
		echo "Info: SLURM file: $RNslurmfile"
	fi
	
	return 0
}
#[General]
#job_type = local #job_type 设置运行环境，可以使用（local， sge， pbs等）
#job_prefix = nextDenovo #运行工作的前缀. (default: nextDenovo)
#task = all # 'all', 'correct', 'assemble' #可以进行针对性选择
#rewrite = yes # yes/no # 再次运行是否覆盖之前结果
#deltmp = yes  #删除中间结果（默认：是）
#rerun = 3 #重新运行未完成的作业，直到完成或达到重新运行循环，0=否。 （默认：3）
#parallel_jobs = 5 #用于并行运行的任务数。 （默认：10）
#input_type = raw  #  raw, corrected； 输入数据是否是raw或corrected
#read_type = clr # 输入数据类型, clr=PacBio 连续长读取，hifi=PacBio 高精度长读取，ont=NanoPore 读取。 （必需的）
#input_fofn = run.fofn # reads 文件（必需）
#workdir = 01_rundir # 工作目录
# cluster_options =auto # 用于定义每个作业的资源需求的模板
#
#[correct_option]
#read_cutoff = 1k #过滤读取长度 < read_cutoff。 （默认：1k）
#genome_size = 0.4M  #估计基因组大小，识别后缀K/M/G，用于计算seed_cutoff/seed_cutfiles/blocksize和平均深度，手动设置seed_cutoff时可以省略。
#seed_cutoff = 0 #最小种子长度，<=0 表示使用 bin/seq_stat 自动计算
#blocksize = 3g #并行运行的块大小，将非种子读取拆分为小文件，每个文件的最大大小为块大小。 （默认：10g）
#pa_correction = 3 #用于并行运行的更正任务数，每个更正任务需要 ~TOTAL_INPUT_BASES/4 字节的内存使用量，仅为此步骤覆盖 parallel_jobs。 （默认：3）
#seed_cutfiles = 5  #拆分种子读入 seed_cutfiles 子文件。 （默认：pa_correction）
#sort_options = -m 20g -t 8 -k 40 #排序选项，-m指设置最大可用缓冲区大小，后缀 K/M/G [40G]，-t:要使用的线程数 [8],-k指max depth of each overlap, should <= average #sequencing depth [40]
#minimap2_options_raw = -x ava-pb -t 8 #minimap2 选项，用于查找原始读取之间的重叠
#correction_options = -p 8 #设置用于校正的进程数。 （默认：10）
#
#[assemble_option]
#random_round = 20 #建议设置20-100. 该参数是设置随机组装参数的数量，会基于每一套随机参数做一次组装， 避免默认参数效果不好
#minimap2_options_cns = -x ava-pb -t 8 -k17 -w17 #用于查找更正读取之间的重叠,-k指k-mer 大小（不大于 28），用于重新对齐 [17]，-w指minizer window size, used to re-align [10]
#nextgraph_options = -a 1 #输出格式，0=无，1=fasta，2=graphml，3=gfa2，4=path [1]



runSmartdenovo () {
	local RSfq=$1
	local RSdir=$2
	
	local RSslurmfile="${RSdir}/${RSdir##*/}.command"
	
	if [ -z "$RSdir" ];then
		RSdir=$PWD
	fi
	if [ ! -d $RSdir ]; then
		mkdir -p $RSdir
	fi
	
	fq2Fa $RSfq "${RSdir}/${RSdir##*/}.fq2fa.fa"
	headerSlurm "${RSdir##*/}" "$RSslurmfile"
	echo "cd $RSdir" >> "$RSslurmfile"
	echo "time smartdenovo.pl -p ${RSdir##*/} -c 1 ${RSdir}/${RSdir##*/}.fq2fa.fa > ${RSdir}/${RSdir##*/}.mak 2> ${RSdir}/${RSdir##*/}.fq2fa.err" >> "$RSslurmfile"
	echo "time make -f ${RSdir}/${RSdir##*/}.mak > ${RSdir}/${RSdir##*/}.mak.log 2>&1" >> "$RSslurmfile"
	if [ $opt_sa -eq 0 ]; then
		cd $RSdir
		if [ $verbose -eq 1 ]; then
			echo "time smartdenovo.pl -p ${RSdir##*/} -c 1 ${RSdir}/${RSdir##*/}.fq2fa.fa > ${RSdir}/${RSdir##*/}.mak 2> ${RSdir}/${RSdir##*/}.fq2fa.err"
		fi
	
		time smartdenovo.pl -p ${RSdir##*/} -c 1 ${RSdir}/${RSdir##*/}.fq2fa.fa > ${RSdir}/${RSdir##*/}.mak 2> ${RSdir}/${RSdir##*/}.fq2fa.err
		if [ $? -ne 0 ] || [ ! -s "${RSdir}/${RSdir##*/}.mak" ]; then
			echo "Error: smartdenovo.pl running error" >&2
			echo "    CMD: time smartdenovo.pl -p ${RSdir##*/} -c 1 ${RSdir}/${RSdir##*/}.fq2fa.fa > ${RSdir}/${RSdir##*/}.mak 2> ${RSdir}/${RSdir##*/}.fq2fa.err" >&2
			exit 100
		else
			if [ -e "${RSdir}/${RSdir##*/}.fq2fa.err" ]; then
				rm "${RSdir}/${RSdir##*/}.fq2fa.err" >/dev/null 2>&1
			fi
		fi
		if [ $verbose -eq 1 ]; then
			echo "time make -f ${RSdir}/${RSdir##*/}.mak > ${RSdir}/${RSdir##*/}.mak.log 2>&1"
		fi
		time make -f ${RSdir}/${RSdir##*/}.mak > ${RSdir}/${RSdir##*/}.mak.log 2>&1
		if [ $? -ne 0 ] || [ ! -s "$RSdir/${RSdir##*/}.dmo.cns" ]; then
			echo "Error: smartdenovo make running error" >&2
			echo "    CMD: time make -f ${RSdir}/${RSdir##*/}.mak > ${RSdir}/${RSdir##*/}.mak 2>&1" >&2
			exit 100
		else
			if [ -e "${RSdir}/${RSdir##*/}.mak.log" ]; then
				rm "${RSdir}/${RSdir##*/}.mak.log" >/dev/null 2>&1
			fi
		fi
	else
		echo "Info: SLURM file: $RSslurmfile"
	fi
	
	return 0
}



#################### Command test ###################################
runStep1=0;
runStep2=0;
runStep3=0;
runStep4=0;
runStep5=0;
runStep6=0;
runStep7=0;
for inds in "${opt_s[@]}"; do
	case $inds in
		1)  runStep1=1;;
		2)  runStep2=1;;
		3)  runStep3=1;;
		4)  runStep4=1;;
		5)  runStep5=1;;
		6)  runStep6=1;;
		7)  runStep7=1;;
		*)  echo "Error: invalid step: $inds" >&2; exit 100;;
	esac
done


if [ $runStep2 -eq 1 ]; then
	CmdExit 'mason'
fi
if [ $runStep3 -eq 1 ]; then
	CmdExit 'fastq-sample'
fi
if [ $runStep4 -eq 1 ] || [ $runStep5 -eq 1 ]; then
	CmdExit 'fastqc'
fi
if [ $runStep4 -eq 1 ]; then
	CmdExit 'seqkit'
	CmdExit 'fastq_checkid.pl'
fi
if [ $runStep5 -eq 1 ]; then
	CmdExit 'fastp'
fi
if [ $runStep6 -eq 1 ]; then
	CmdExit 'bwa'
	CmdExit 'samtools'
	CmdExit 'bamCoverage'
fi



#################### Defaults #######################################
test2run_abyss=0;
test2run_edena=0;
test2run_minia=0;
test2run_ray=0;
test2run_spades=0;
test2run_canu=0;
test2run_flye=0;
test2run_mecat=0;
test2run_nextdenovo=0;
test2run_smartdenovo=0;
declare -a kmers=()
printSeperator
if [ ! -s $opt_mt ]; then
	echo "    Mt: invalid file: $opt_mt" >&2
	exit 100
fi
if [ ! -s $opt_ct ]; then
	echo "    Ct: invalid file: $opt_ct" >&2
	exit 100
fi
if [ ! -s $opt_gn ]; then
	echo "    Gn: invalid file: $opt_gn" >&2
	exit 100
fi
echo "Genomes:"
echo "    Mt: $opt_mt"
echo "    Ct: $opt_ct"
echo "    Gn: $opt_gn"
echo "Info: Platform: $opt_pf"
if [ "$opt_pf" == "illumina" ]; then
### check fragment length
	if [[ "$opt_fl" =~ ^[0-9]+$ ]] && [ $opt_fl -gt 0 ]; then
		fragment_length_stdev=$(($opt_fl*2/5))
		echo "    Mate-pair mean library fragment length: $opt_rl"
		echo "    Mate-pair library tolerance: $fragment_length_stdev"
	else
		echo "Error: invalid fragment length --fl: $opt_fl"
		exit 100
	fi

	if [[ "$opt_rl" =~ ^[0-9]+$ ]] && [ $opt_rl -gt 0 ]; then
		echo "    Read length: $opt_rl"
	else
		opt_rl=150
		echo "    Read length: $opt_rl; reset"
	fi
	echo "    Minimum read length to keep a read during trimming: $opt_tl"
	if [ -z "$opt_as" ]; then
		opt_as="abyss,edena,minia,spades,ray"
	fi
	Assembler2arr=($(echo $opt_as | tr ',' "\n"))
	for indas in ${Assembler2arr[@]}; do
		case $indas in
			abyss) test2run_abyss=1;CmdExit 'abyss-pe'; echo "    Assemblers: abyss";;
			edena) test2run_edena=1;CmdExit 'edena';echo "    Assemblers: edena";;
			minia) test2run_minia=1;CmdExit 'minia';echo "    Assemblers: minia";;
			spades) test2run_spades=1;CmdExit 'spades.py';echo "    Assemblers: spades";;
			ray) test2run_ray=1;CmdExit 'Ray';echo "    Assemblers: ray";;
		esac
	done
	### check Kmers
	if [ $test2run_abyss -eq 1 ] || [ $test2run_minia -eq 1 ] || [ $test2run_ray -eq 1 ] || [ $test2run_spades -eq 1 ]; then
		if [ -z "$opt_km" ] ; then
			echo "Error: please specify assembling kmers with --km for K-mer dependent assemblers" >&2
			exit 100
		fi
		declare -a kmer_tmp=($(echo $opt_km | tr ',' "\n"))
		for indk in ${kmer_tmp[@]}; do
			if [[ "$indk" =~ ^[0-9]+$ ]] && [ $(($indk%2)) -eq 1 ] && [ $indk -gt 0 ] && [ $indk -lt $opt_rl ]; then
				kmers=(${kmers[@]} $indk)
			else
				echo "Error: invalid kmer setting: $indk" >&2
				exit 100
			fi
		done
		if [ ${#kmers[@]} -lt 1 ]; then
			echo "Error: empty usable kmers" >&2
			exit 100
		else
			echo "    Kmers: ${kmers[@]}"
		fi
	fi
elif [ "$opt_pf" == "pacbio" ]; then
	if [[ "$opt_rl" =~ ^[0-9]+,[0-9]+$ ]]; then
		opt_rd=${opt_rl##*,}
		opt_rl=${opt_rl%%,*}
		echo "    Read length: $opt_rl"
		echo "    Read length sd: $opt_rd"
	else
		opt_rd=15000
		opt_rl=13000
		echo "    Read length: $opt_rl; reset"
		echo "    Read length sd: $opt_rd; reset"
	fi
	if [ -z "$opt_as" ]; then
		opt_as="canu,flye,mecat,nextdenovo,smartdenovo"
	fi
	Assembler3arr=($(echo $opt_as | tr ',' "\n"))
	for indas in ${Assembler3arr[@]}; do
		case $indas in
			canu) test2run_canu=1;CmdExit 'canu'; echo "    Assemblers: canu";;
			flye) test2run_flye=1;CmdExit 'flye'; echo "    Assemblers: flye";;
			mecat) test2run_mecat=1;CmdExit 'mecat.pl'; echo "    Assemblers: mecat";;
			nextdenovo) test2run_nextdenovo=1;CmdExit 'nextDenovo'; echo "    Assemblers: nextdenovo";;
			smartdenovo) test2run_smartdenovo=1;CmdExit 'smartdenovo.pl'; CmdExit 'make';echo "    Assemblers: smartdenovo";;
		esac
	done
fi
if [[ "$opt_mm" =~ ^[0-9]$ ]]; then
	echo "Error: max memory should be an integer: $opt_mm" >&2
	exit 100
else
	echo "Max memory: $opt_mm"
fi
echo "Contamination rate: ${opt_contamination_proportion[@]}"
echo "Read depth: $opt_mc"
echo "Replicate: $opt_rp"
echo "Steps: ${opt_s[@]}"
echo "Threads: $opt_t"
echo "Out path: $outdir"
### mason options
readonly masonOptions=" -aNg -sq -n $opt_rl -ll $opt_fl -le $fragment_length_stdev -mp --read-naming 2 "
### step number
step=0
readonly randSeed=9000
sizeMt=0
samtoolsVers=$(samtools 2>&1 | grep ^'Version' | sed 's/^Version:\s\+//;s/\..*$//;')



#################### Input and Output ###############################
#opt_gn="/home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta"
#opt_mt="/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.mason.500bp.fa"
#opt_ct="/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrC.fasta"
if [ -z "$opt_gn" ] && [ ! -s $opt_gn ]; then
	echo "Error: invalid genome fasta" >&2
	exit 100
fi
if [ -z "$opt_mt" ] && [ ! -s $opt_mt ]; then
	echo "Error: invalid mitochondria fasta" >&2
	exit 100
fi
if [ -z "$opt_ct" ] && [ ! -s $opt_ct ]; then
	echo "Error: invalid chloroplast fasta" >&2
	exit 100
fi



#################### Main ###########################################



### Step1: calculate genome sizes
((step++))
printSeperator
echo "(Step$step)Info: calculate genome sizes"
echo "(Step$step)Info: calculate genome sizes" >&2
runDir="$outdir/1.genome"
if [ ! -d "$runDir" ]; then
	mkdir -p $runDir
fi
cd $runDir
if [ $runStep1 -eq 1 ] || [ $runStep2 -eq 1 ] || [ $runStep3 -eq 1 ]; then
	if [[ "$opt_mt" =~ \.[gG][zZ]$ ]]; then
		opt_mt_decompress=${opt_mt##*/};
		opt_mt_decompress=${opt_mt_decompress%.*};
		opt_mt_decompress="$runDir/${opt_mt_decompress%.*}.fa"
		if [ $keep_files -eq 1 ] && [ -s "$opt_mt_decompress" ];then
			echo "(Step$step)Info: using existing mitchondria genome sequences: $opt_mt_decompress"
		else
			echo "(Step$step)Info: decompress mitchondria genome sequences"
			fileGunzip $opt_mt $opt_mt_decompress
		fi
		opt_mt=$opt_mt_decompress
	else
		echo "(Step$step)Info: mitchondria genome sequences in flat text: $opt_mt"
	fi
	if [[ "$opt_ct" =~ \.[gG][zZ]$ ]]; then
		opt_ct_decompress=${opt_ct##*/};
		opt_ct_decompress=${opt_ct_decompress%.*};
		opt_ct_decompress="$runDir/${opt_ct_decompress%.*}.fa"
		if [ $keep_files -eq 1 ] && [ -s "$opt_ct_decompress" ];then
			echo "(Step$step)Info: using existing chloroplast genome sequences: $opt_ct_decompress"
		else
			echo "(Step$step)Info: decompress chloroplast genome sequences"
			fileGunzip $opt_ct $opt_ct_decompress
		fi
		opt_ct=$opt_ct_decompress
	else
		echo "(Step$step)Info: chloroplast genome sequences in flat text: $opt_ct"
	fi
	if [[ "$opt_gn" =~ \.[gG][zZ]$ ]]; then
		opt_gn_decompress=${opt_gn##*/};
		opt_gn_decompress=${opt_gn_decompress%.*};
		opt_gn_decompress="$runDir/${opt_gn_decompress%.*}.fa"
		if [ $keep_files -eq 1 ] && [ -s "$opt_gn_decompress" ];then
			echo "(Step$step)Info: using existing nuclear genome sequences: $opt_gn_decompress"
		else
			echo "(Step$step)Info: decompress nuclear genome sequences"
			fileGunzip $opt_gn $opt_gn_decompress
		fi
		opt_gn=$opt_gn_decompress
	else
		echo "(Step$step)Info: nuclear genome sequences in flat text: $opt_gn"
	fi
	sizeMt=$(getFastaSize $opt_mt)
	if [[ "$sizeMt" =~ ^[0-9]+$ ]] && [ $sizeMt -gt 0 ]; then
		echo "(Step$step)     Mitochondria seqence size: $sizeMt"
	else
		echo "(Step$step)Error: invalid mitochondria size: $sizeMt" >&2
		exit 100
	fi
else
	echo "(Step$step)Warnings: Step1 skipped"
fi



### Step2: simulate total reads
((step++))
printSeperator
echo "(Step$step)Info: read simulation using mason"
echo "(Step$step)Info: read simulation using mason" >&2
runDir="$outdir/2.simulation"
if [ ! -d "$runDir" ]; then
	mkdir -p $runDir
fi
cd $runDir
covMt=0
for indn in "${opt_mc[@]}"; do
	((indn > covMt)) && covMt=$indn
done
if [[ "$covMt" =~ ^[0-9]+$ ]] && [ $covMt -gt 0 ]; then
	echo "(Step$step)      Mitochondria max read depth: $covMt"
else
	echo "(Step$step)Error: invalid max mtDNA read depth: $covMt, please check -mc parameter" >&2
	exit 100
fi
((covMt=covMt*2))
if [ "$opt_pf" == "illumina" ]; then
	readnumMt=$(echo -e "$sizeMt\t$covMt\t$opt_rl" | awk '{print int($1*$2/(2*$3))}')
elif [ "$opt_pf" == "pacbio" ]; then
	readnumMt=$(echo -e "$sizeMt\t$covMt\t$opt_rl" | awk '{print int($1*$2/($3))}')
else
	echo "(Step$step)Error: unknown platform: $opt_pf" >&2
	exit 100
fi
echo "(Step$step)      Number of read pairs for Mitochondria to be simulated: $readnumMt, depth: $covMt, size: $sizeMt"
readnumCt=$(echo -e "$readnumMt\t${opt_contamination_proportion[1]}" | awk '{print int($1*$2)}')
echo "(Step$step)      Number of read pairs for chloroplast  to be simulated: $readnumCt, Proportion: ${opt_contamination_proportion[1]}"
readnumGn=$(echo -e "$readnumMt\t${opt_contamination_proportion[0]}" | awk '{print int($1*$2)}')
echo "(Step$step)      Number of read pairs for genome       to be simulated: $readnumGn, Proportion: ${opt_contamination_proportion[0]}"
mtFqPfx="mtDNA.${opt_pf}.mt${covMt}x"
ctFqPfx="ctDNA.${opt_pf}.mt${covMt}x"
gnFqPfx="gnDNA.${opt_pf}.mt${covMt}x"
if [ "$opt_pf" == "illumina" ]; then
	masonMtR1="$runDir/$mtFqPfx.R1.fastq"
	masonMtR2="$runDir/$mtFqPfx.R2.fastq"
	masonCtR1="$runDir/$ctFqPfx.R1.fastq"
	masonCtR2="$runDir/$ctFqPfx.R2.fastq"
	masonGnR1="$runDir/$gnFqPfx.R1.fastq"
	masonGnR2="$runDir/$gnFqPfx.R2.fastq"
	if [ $runStep2 -eq 1 ]; then
		if [ $keep_files -eq 1 ] && [ -s $masonMtR1 ] && [ -s $masonMtR2 ]; then
			echo "(Step$step)Info: using existing mitochondia FastQ: $masonMtR1 $masonMtR2"
		else
			if [ $verbose -eq 1 ]; then
				echo "mason2reads $opt_mt $readnumMt $runDir/$mtFqPfx.fastq mtDNA.mt${covMt}x.${readnumMt}"
			fi
			mason2reads $opt_mt $readnumMt "$runDir/$mtFqPfx.fastq" "mtDNA.mt${covMt}x.${readnumMt}"
		fi
		if [ $keep_files -eq 1 ] && [ -s $masonCtR1 ] && [ -s $masonCtR2 ]; then
			echo "(Step$step)Info: using existing chloroplast FastQ: $masonCtR1 $masonCtR2"
		else
			if [ $verbose -eq 1 ]; then
				echo "mason2reads $opt_ct $readnumCt $runDir/$ctFqPfx.fastq ctDNA.mt${covMt}x.${readnumCt}"
			fi
			mason2reads $opt_ct $readnumCt "$runDir/$ctFqPfx.fastq" "ctDNA.mt${covMt}x.${readnumCt}"
		fi
		if [ $keep_files -eq 1 ] && [ -s $masonGnR1 ] && [ -s $masonGnR2 ]; then
			echo "(Step$step)Info: using existing nuclear genome FastQ: $masonGnR1 $masonGnR2"
		else
			if [ $verbose -eq 1 ]; then
				echo "mason2reads $opt_gn $readnumGn $runDir/$gnFqPfx.fastq gnDNA.mt${covMt}x.${readnumGn}"
			fi
			mason2reads $opt_gn $readnumGn "$runDir/$gnFqPfx.fastq" "gnDNA.mt${covMt}x.${readnumGn}"
		fi
	else
		echo "(Step$step)Warnings: Step2 skipped"
	fi
elif [ "$opt_pf" == "pacbio" ]; then
	badreadMt="$runDir/$mtFqPfx.fastq"
	badreadCt="$runDir/$ctFqPfx.fastq"
	badreadGn="$runDir/$gnFqPfx.fastq"
	if [ $runStep2 -eq 1 ]; then
		if [ $keep_files -eq 1 ] && [ -s $badreadMt ]; then
			echo "(Step$step)Info: using existing mitochondia PacBio FastQ: $badreadMt"
		else
			if [ $opt_debug -eq 1 ]; then
				echo "(Step$step)Info: badread2reads $opt_mt $readnumMt $badreadMt"
			fi
			badread2reads $opt_mt $readnumMt $badreadMt
		fi
		if [ $keep_files -eq 1 ] && [ -s $badreadCt ]; then
			echo "(Step$step)Info: using existing chloroplast PacBio FastQ: $badreadCt"
		else
			if [ $opt_debug -eq 1 ]; then
				echo "(Step$step)Info: badread2reads $opt_ct $readnumCt $badreadCt"
			fi
			badread2reads $opt_ct $readnumCt $badreadCt
		fi
		if [ $keep_files -eq 1 ] && [ -s $badreadGn ]; then
			echo "(Step$step)Info: using existing nuclear genome PacBio FastQ: $badreadGn"
		else
			if [ $opt_debug -eq 1 ]; then
				echo "(Step$step)Info: badread2reads $opt_gn $readnumGn $badreadGn"
			fi
			badread2reads $opt_gn $readnumGn $badreadGn
		fi
	else
		echo "(Step$step)Warnings: Step2 skipped"
	fi
fi



((step++))
printSeperator
runDir="$outdir/3.sampling"
echo "(Step$step)Info: read sampling using fastq-sample"
echo "(Step$step)Info: read sampling using fastq-sample" >&2
if [ ! -d "$runDir" ]; then
	mkdir -p $runDir
fi
cd $runDir
if [ $runStep3 -eq 1 ]; then
	for indcov in "${opt_mc[@]}"; do
		if [ "$opt_pf" == "illumina" ]; then
			samplenumMt=$(echo -e "$sizeMt\t$indcov\t$opt_rl" | awk '{print int($1*$2/(2*$3))}')
		elif [ "$opt_pf" == "pacbio" ]; then
			samplenumMt=$(echo -e "$sizeMt\t$indcov\t$opt_rl" | awk '{print int($1*$2/($3))}')
		fi
		echo "(Step$step)      Number of read pairs for Mitochondria to be simulated: $samplenumMt, depth: $indcov, size: $sizeMt"
		samplenumCt=$(echo -e "$samplenumMt\t${opt_contamination_proportion[1]}" | awk '{print int($1*$2)}')
		echo "(Step$step)      Number of read pairs for chloroplast  to be simulated: $samplenumCt, Proportion: ${opt_contamination_proportion[1]}"
		samplenumGn=$(echo -e "$samplenumMt\t${opt_contamination_proportion[0]}" | awk '{print int($1*$2)}')
		echo "(Step$step)      Number of read pairs for genome       to be simulated: $samplenumGn, Proportion: ${opt_contamination_proportion[0]}"
		for ((rep=1;rep<=opt_rp;rep++)); do
			((randNum=rep+9))
			if [ "$opt_pf" == "illumina" ]; then
				masonMtR1="$outdir/2.simulation/$mtFqPfx.R1.fastq"
				masonMtR2="$outdir/2.simulation/$mtFqPfx.R2.fastq"
				masonCtR1="$outdir/2.simulation/$ctFqPfx.R1.fastq"
				masonCtR2="$outdir/2.simulation/$ctFqPfx.R2.fastq"
				masonGnR1="$outdir/2.simulation/$gnFqPfx.R1.fastq"
				masonGnR2="$outdir/2.simulation/$gnFqPfx.R2.fastq"
				
				sampleMtR1="${runDir}/mtDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				sampleMtR2="${runDir}/mtDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $sampleMtR1 ] && [ -s $sampleMtR2 ]; then
					echo "(Step$step)Info: using existing mitchondria Illumina FastQ sample: $sampleMtR1 $sampleMtR2"
				else
					fqSampler $masonMtR1 $sampleMtR1 ${samplenumMt} ${randNum}
					fqSampler $masonMtR2 $sampleMtR2 ${samplenumMt} ${randNum}
#					fqVerifyNumber $samplenumMt $sampleMtR1 $sampleMtR2
				fi
				sampleCtR1="${runDir}/ctDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				sampleCtR2="${runDir}/ctDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $sampleCtR1 ] && [ -s $sampleCtR2 ]; then
					echo "(Step$step)Info: using existing mitchondria Illumina FastQ sample: $sampleCtR1 $sampleCtR2"
				else
					fqSampler $masonCtR1 $sampleCtR1 ${samplenumCt} ${randNum}
					fqSampler $masonCtR2 $sampleCtR2 ${samplenumCt} ${randNum}
#					fqVerifyNumber ${samplenumCt} $sampleCtR1 $sampleCtR2
				fi
				sampleGnR1="${runDir}/gnDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				sampleGnR2="${runDir}/gnDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $sampleGnR1 ] && [ -s $sampleGnR2 ]; then
					echo "(Step$step)Info: using existing mitchondria Illumina FastQ sample: $sampleGnR1 $sampleGnR2"
				else
					fqSampler $masonGnR1 $sampleGnR1 ${samplenumGn} ${randNum}
					fqSampler $masonGnR2 $sampleGnR2 ${samplenumGn} ${randNum}
#					fqVerifyNumber ${samplenumGn} $sampleGnR1 $sampleGnR2
				fi
			elif [ "$opt_pf" == "pacbio" ]; then
				badreadMt="$outdir/2.simulation/$mtFqPfx.fastq"
				badreadCt="$outdir/2.simulation/$ctFqPfx.fastq"
				badreadGn="$outdir/2.simulation/$gnFqPfx.fastq"
				
				sampleMtBd="${runDir}/mtDNA.${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $sampleMtBd ]; then
					echo "(Step$step)Info: using existing mitchondria PacBio FastQ sample: $sampleMtBd"
				else
					if [ $opt_debug -eq 1 ]; then
						echo "fqSampler $badreadMt $sampleMtBd ${samplenumMt} ${randNum}"
					fi
					fqSampler $badreadMt $sampleMtBd ${samplenumMt} ${randNum}
				fi
#				fqVerifyNumber ${samplenumMt} $sampleMtBd
				sampleCtBd="${runDir}/ctDNA.${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $sampleCtBd ]; then
					echo "(Step$step)Info: using existing mitchondria PacBio FastQ sample: $sampleCtBd"
				else
					if [ $opt_debug -eq 1 ]; then
						echo "fqSampler $badreadCt $sampleCtBd ${samplenumCt} ${randNum}"
					fi
					fqSampler $badreadCt $sampleCtBd ${samplenumCt} ${randNum}
				fi
#				fqVerifyNumber ${samplenumCt} $sampleCtBd
				sampleGnBd="${runDir}/gnDNA.${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $sampleGnBd ]; then
					echo "(Step$step)Info: using existing mitchondria PacBio FastQ sample: $sampleGnBd"
				else
					if [ $opt_debug -eq 1 ]; then
						echo "fqSampler $badreadGn $sampleGnBd ${samplenumGn} ${randNum}"
					fi
					fqSampler $badreadGn $sampleGnBd ${samplenumGn} ${randNum}
				fi
#				fqVerifyNumber ${samplenumGn} $sampleGnBd
			fi
		done
	done
else
	echo "(Step$step)Warnings: Step3 skipped"
fi



((step++))
printSeperator
runDir="$outdir/4.merge"
echo "(Step$step)Info: merge and QC"
echo "(Step$step)Info: merge and QC" >&2
if [ ! -d "$runDir" ]; then
	mkdir -p $runDir
fi
cd $runDir
if [ $runStep4 -eq 1 ]; then
	for indcov in "${opt_mc[@]}"; do
		for ((rep=1;rep<=$opt_rp;rep++)); do
			if [ "$opt_pf" == "illumina" ]; then
				mergeR1="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				mergeR2="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				mtR1="$outdir/3.sampling/mtDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				mtR2="$outdir/3.sampling/mtDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				ctR1="$outdir/3.sampling/ctDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				ctR2="$outdir/3.sampling/ctDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				gnR1="$outdir/3.sampling/gnDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				gnR2="$outdir/3.sampling/gnDNA.${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $mergeR1 ]; then
					echo "(Step$step)Info: using existing merged Fastq R1: $mergeR1"
				else
					fqMerge $mergeR1 $mtR1 $ctR1 $gnR1 $randSeed
				fi
				if [ $keep_files -eq 1 ] && [ -s $mergeR2 ]; then
					echo "(Step$step)Info: using existing merged Fastq R2: $mergeR2"
				else
					fqMerge $mergeR2 $mtR2 $ctR2 $gnR2 $randSeed
				fi
				if [ $keep_files -eq 1 ] && [ -s "${mergeR1%.fastq.gz}_fastqc.zip" ] && [ -s "${mergeR2%.fastq.gz}_fastqc.zip" ]; then
					echo "(Step$step)Info: using existing FastQC report: $mergeR1  $mergeR2"
				else
					fastq_checkid.pl  $mergeR1  $mergeR2 '\@(\S+)\/[12]\s*\S*'
					if [ $? -ne 0 ]; then
						echo "(Step$step)Error: FastQ ID not paired: $mergeR1 $mergeR2" >&2
						exit 100
					fi
					fastqc --noextract --nogroup --format fastq -o $runDir --threads $opt_t --quiet $mergeR1  $mergeR2
					if [ $? -ne 0 ]; then
						echo "(Step$step)Error: FastQC failed: $mergeR1 $mergeR2" >&2
						exit 100
					fi
				fi
			elif [ "$opt_pf" == "pacbio" ]; then
				mergeFq="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				mtFq="$outdir/3.sampling/mtDNA.${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				ctFq="$outdir/3.sampling/ctDNA.${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				gnFq="$outdir/3.sampling/gnDNA.${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $mergeFq ]; then
					echo "(Step$step)Info: using existing merged Fastq: $mergeFq"
				else
					fqMerge $mergeFq $mtFq $ctFq $gnFq $randSeed
				fi
			fi
		done
	done
	echo "(Step$step)Info: simulation fastqs are in $runDir"
else
	echo "(Step$step)Warnings: Step4 skipped"
fi



((step++))
printSeperator
runDir="$outdir/5.fastp"
echo "(Step$step)Info: FastQ trim and QC"
echo "(Step$step)Info: FastQ trim and QC" >&2
if [ ! -d "$runDir" ]; then
	mkdir -p $runDir
fi
cd $runDir
if [ $runStep5 -eq 1 ]; then
	for indcov in "${opt_mc[@]}"; do
		for ((rep=1;rep<=$opt_rp;rep++)); do
			if [ "$opt_pf" == "illumina" ]; then
				trimIn1="$outdir/4.merge/ALL${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				trimIn2="$outdir/4.merge/ALL${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				trimOut1="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.R1.fastq.gz"
				trimOut2="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.R2.fastq.gz"
				
				if [ $keep_files -eq 1 ] && [ -s $trimOut1 ] && [ -s $trimOut2 ]; then
					echo "(Step$step)Info: using existing trim FastQ: $trimOut1 $trimOut2"
				else
					echo "(Step$step)Info: trim $trimIn1";
					echo "(Step$step)Info: trim $trimIn1" >&2;
					echo "(Step$step)Info: trim $trimIn2"
					echo "(Step$step)Info: trim $trimIn2" >&2
					fastp -i $trimIn1  -o $trimOut1  -I $trimIn2 -O $trimOut2 -l $opt_tl --thread $opt_t
					if [ $? -ne 0 ] || [ ! -s $trimOut1 ] || [ ! -s $trimOut2 ]; then
						echo "(Step$step)Error: fastp runing error" >&2
						echo "(Step$step)    CMD used: fastp -i $trimIn1  -o $trimOut1  -I $trimIn2 -O $trimOut2 -l $opt_tl --thread $opt_t" >&2
						exit 100
					fi
				fi
				if [ $keep_files -eq 1 ] && [ -s "${trimOut1%.fastq.gz}_fastqc.zip" ] && [ -s "${trimOut2%.fastq.gz}_fastqc.zip" ]; then
					echo "(Step$step)Info: using existing trim FastQC: $trimOut1 $trimOut2"
				else
					fastqc --noextract --nogroup --format fastq -o $runDir --threads $opt_t --quiet $trimOut1  $trimOut2
					if [ $? -ne 0 ]; then
						echo "(Step$step)Error: FastQC failed" >&2
						echo "(Step$step)    CMD used: fastqc --noextract --nogroup --format fastq -o $runDir --threads $opt_t $trimOut1  $trimOut2" >&2
					fi
				fi
			fi
		done
	done
else
	echo "(Step$step)Warnings: Step5 skipped"
fi



((step++))
printSeperator
runDir="$outdir/6.bwa"
echo "(Step$step)Info: Step$step: Mapping with BWA"
echo "(Step$step)Info: Step$step: Mapping with BWA" >&2
if [ ! -d "$runDir" ]; then
	mkdir -p $runDir
fi
cd $runDir
if [ $runStep6 -eq 1 ]; then
	if [ "$opt_pf" == "illumina" ]; then
		bwaIndexPfx="$runDir/Mt"
		if [ $keep_files -eq 1 ] && [ -s "$bwaIndexPfx.amb" ] && [ -s "$bwaIndexPfx.ann" ] && [ -s "$bwaIndexPfx.pac" ] && [ -s "$bwaIndexPfx.sa" ] && [ -s "$bwaIndexPfx.bwt" ]; then
			echo "(Step$step)Info: using existing BWA index: $bwaIndexPfx"
		else
			bwa index -p "Mt" $opt_mt > $runDir/bwa.index.log 2>&1
			if [ $? -ne 0 ]; then
				echo "(Step$step)Error: BWA index error" >&2
				echo "(Step$step)    CMD used: bwa index -p Mt $opt_mt > $runDir/bwa.index.log 2>&1" >&2
				exit 100
			else
				if [ -e "$runDir/bwa.index.log" ]; then
					rm -rf "$runDir/bwa.index.log" >/dev/null 2>&1
				fi
			fi
		fi
	fi
	for indcov in "${opt_mc[@]}"; do
		for ((rep=1;rep<=$opt_rp;rep++)); do
			Bam1Out="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.bam"
			Bam2Sort="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.sort.bam"
			Bam3depth="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.sort.bam.depth"
			Bam4bigwig="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.sort.bam.bw"
			if [ "$opt_pf" == "illumina" ]; then
				trimOut1="$outdir/5.fastp/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.R1.fastq.gz"
				trimOut2="$outdir/5.fastp/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.R2.fastq.gz"
				if [ $keep_files -eq 1 ] && [ -s $Bam2Sort -o -s $Bam1Out ]; then
					if [ -s $Bam2Sort ]; then
						echo "(Step$step)Info: existing sorted BAM: $Bam2Sort"
					elif [ -s $Bam1Out ]; then
						echo "(Step$step)Info: existing BAM: $Bam1Out"
					fi
				else
					bwa mem -t $opt_t $bwaIndexPfx  $trimOut1  $trimOut2 2> $Bam1Out.BWA.bwamem.error | samtools view -h -S -b - > $Bam1Out 2> $Bam1Out.BWA.samtoolsview.error
					if [ $? -ne 0 ] || [ ! -s $Bam1Out ]; then
						echo "(Step$step)Error: BWA mapping error" >&2
						echo "(Step$step)    CMD used: bwa mem -t $opt_t $bwaIndexPfx  $trimOut1  $trimOut2 2> $Bam1Out.BWA.bwamem.error | samtools view -h -S -b - > $Bam1Out 2> $Bam1Out.BWA.samtoolsview.error" >&2
						exit 100
					else
						if [ -e "$Bam1Out.BWA.bwamem.error" ]; then
							rm -rf "$Bam1Out.BWA.bwamem.error" > /dev/null 2>&1
						fi
						if [ -e "$Bam1Out.BWA.samtoolsview.error" ]; then
							rm -rf "$Bam1Out.BWA.samtoolsview.error" > /dev/null 2>&1
						fi
					fi
				fi
			elif [ "$opt_pf" == "pacbio" ]; then
				if [ $keep_files -eq 1 ] && [ -s "$Bam1Out" -o -s "$Bam2Sort" ]; then
					if [ -s "$Bam2Sort" ]; then
						echo "(Step$step)Info: existing sorted BAM: $Bam2Sort"
					elif [ -s "$Bam1Out" ]; then
						echo "(Step$step)Info: existing BAM: $Bam1Out"
					fi
				else
					minimap2 -ax map-pb $opt_mt "$outdir/4.merge/ALL${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz" 2> $Bam1Out.minimap2.err | samtools view -h -S -b - > $Bam1Out 2> $Bam1Out.samtoolsview.err
					if [ $? -ne 0 ] || [ ! -s $Bam1Out ]; then
						echo "(Step$step)Error: minimap2 mapping error" >&2
						echo "(Step$step)    CMD used: minimap2 -ax map-pb $opt_mt "$outdir/4.merge/ALL${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz" 2> $Bam1Out.minimap2.err | samtools view -h -S -b - > $Bam1Out 2> $Bam1Out.samtoolsview.err" >&2
						exit 100
					else
						if [ -e "$Bam1Out.err" ]; then
							rm "$Bam1Out.err" >/dev/null 2>&1
						fi
					fi
				fi
			fi
			if [ $keep_files -eq 1 ] && [ -s $Bam2Sort ]; then
				echo "(Step$step)Info: using existing sorted BAM: $Bam2Sort"
			else
				bamSort $Bam1Out $Bam2Sort
				rm $Bam1Out > /dev/null 2>&1
			fi
			if [ $keep_files -eq 1 ] && [ -s "$Bam2Sort.bai" ]; then
				echo "(Step$step)Info: using existing sorted BAM index: $Bam2Sort.bai"
			else
				bamIndex $Bam2Sort
			fi
			if [ $keep_files -eq 1 ] && [ -s $Bam3depth ]; then
				echo "(Step$step)Info: using existing BAM depth: $Bam3depth"
			else
				samtools depth -aa -o $Bam3depth $Bam2Sort
				if [ $? -ne 0 ] || [ ! -s $Bam3depth ]; then
					echo "(Step$step)Error: samtools depth error" >&2
					echo "(Step$step)    CMD used: samtools depth -aa -o $Bam3depth $Bam2Sort" >&2
					exit 100
				fi
			fi
			perl -lne 'BEGIN{$sum=0;$count=0;}@F=split(/\t/);$count++;$sum+=$F[2];END{$avg=sprintf "%.2f",$sum/$count; print "Bases: $count\nAverage Depth: $avg";}' $Bam3depth
			if [ $keep_files -eq 1 ] && [ -s $Bam4bigwig ]; then
				echo "(Step$step)Info: using existing BAM BIGWIG: $Bam4bigwig"
			else
				genomeSize=$(samtools view -H $Bam2Sort | grep ^'@SQ' | perl -lne 'BEGIN{print STDERR "Sequences add to total length"; $sum=0;} $line=$_; @arr=split(/\t/, $line); print STDERR "$arr[1]\i$arr[2]";$line=~s/^\@SQ\s+SN\S+\s+LN://; $sum+=$line;END{print $sum;}')
				echo "(Step$step)Info: Genome size: $genomeSize"
				bam2Bigwig $Bam2Sort $Bam4bigwig $binsize
			fi
		done
	done
else
	echo "(Step$step)Warnings: Step6 skipped"
fi



### Assemble
((step++))
printSeperator
runDir="$outdir/7.assembly"
echo "(Step$step)Info: Assembly"
echo "(Step$step)Info: Assembly" >&2
if [ ! -d "$runDir" ]; then
	mkdir -p $runDir
fi
cd $runDir
if [ $runStep7 -eq 1 ]; then
	if [ "$opt_pf" == "illumina" ]; then
		for indcov in "${opt_mc[@]}"; do
			for ((rep=1;rep<=$opt_rp;rep++)); do
				trimIn1="$outdir/4.merge/ALL${opt_pf}.mt${indcov}x.Rep${rep}.R1.fastq.gz"
				trimIn2="$outdir/4.merge/ALL${opt_pf}.mt${indcov}x.Rep${rep}.R2.fastq.gz"
				trimOut1="$outdir/5.fastp/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.R1.fastq.gz"
				trimOut2="$outdir/5.fastp/ALL${opt_pf}.mt${indcov}x.Rep${rep}.trim.R2.fastq.gz"
				if [ ! -s $trimOut1 ] || [ ! -s $trimOut2 ]; then
					echo "(Step$step)Error: missing Illumina fastQ: $trimOut1 $trimOut2" >&2
					exit 100
				fi
				if [ $test2run_abyss -eq 1 ]; then
					abyssOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.abyss"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $abyssOutDir ]; then
					# && [ -L $abyssOutDir/${abyssOutDir##*/}-contigs.fa ]
						echo "(Step$step)Info: existing ABySS out: $abyssOutDir/; skipping"
					else
						runAbyss $trimOut1 $trimOut2 $abyssOutDir
						### out: symlink: $abyssOutDir/K$k/${abyssOutDir##*/}-contigs.fa
						### out: symlink: $abyssOutDir/K$k/${abyssOutDir##*/}-scaffolds.fa
					fi
				fi
				if [ $test2run_edena -eq 1 ]; then
					edenaOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.edena"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $edenaOutDir ] && [ -s "$edenaOutDir/${edenaOutDir##*/}_contigs.fasta" ]; then
						echo "(Step$step)Info: existing Edena out: $edenaOutDir/; skipping"
					else
						min_ol=$(echo "scale=0;($indcov/50)+(2*$opt_rl/3)" |bc -l)
						### Edena require un-trimmed fastq
						runEdena $trimIn1 $trimIn2 $edenaOutDir $min_ol
						### out: $edenaOutDir/${edenaOutDir##*/}_contigs.fasta
					fi
				fi
				if [ $test2run_minia -eq 1 ]; then
					miniaOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.minia"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $miniaOutDir ]; then
						echo "(Step$step)Info: existing Minia out: $miniaOutDir/; skipping"
					else
						runMinia $trimOut1 $trimOut2 $miniaOutDir
#						$miniaOutDir/K$k/${miniaOutDir##*/}.contigs.fa
					fi
				fi
				if [ $test2run_ray -eq 1 ]; then
					rayOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.ray"
					cd $runDir
					runRay $trimOut1 $trimOut2 $rayOutDir
#					$rayOutDir/K$k/Contigs.fasta
#					$rayOutDir/K$k/Scaffolds.fasta
				fi
				if [ $test2run_spades -eq 1 ]; then
					spadesOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.spades"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $spadesOutDir ] && [ -s "$spadesOutDir/contigs.fasta" ]; then
						echo "(Step$step)Info: existing SPAdes out: $spadesOutDir/; skipping"
					else
						runSpades $trimOut1 $trimOut2 $spadesOutDir
#						$spadesOutDir/contigs.fasta
#						$spadesOutDir/scaffolds.fasta
					fi
				fi
			done
		done
	elif [ "$opt_pf" == "pacbio" ]; then
		for indcov in "${opt_mc[@]}"; do
			for ((rep=1;rep<=$opt_rp;rep++)); do
				mixFq="$outdir/4.merge/ALL${opt_pf}.mt${indcov}x.Rep${rep}.fastq.gz"
				if [ ! -s $mixFq ]; then
					echo "(Step$step)Error: missing PacBio fastQ: $trimOut1" >&2
					exit 100
				fi
				if [ $test2run_canu -eq 1 ]; then
					canuOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.canu"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $canuOutDir ] && [ -s "$canuOutDir/${canuOutDir##*/}.contigs.fasta" ]; then
						echo "(Step$step)Info: existing Canu out: $canuOutDir/${canuOutDir##*/}.contigs.fasta; skipping"
					else
						runCanu $mixFq $canuOutDir
#						$canuOutDir/${canuOutDir##*/}.contigs.fasta
					fi
				fi
				if [ $test2run_flye -eq 1 ]; then
					flyeOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.flye"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $flyeOutDir ] && [ -s "$flyeOutDir/assembly.fasta" ]; then
						echo "(Step$step)Info: existing Flye out: $flyeOutDir/; skipping"
					else
						runFlye $mixFq $flyeOutDir
#						"$flyeOutDir/assembly.fasta"
					fi
				fi
				if [ $test2run_mecat -eq 1 ]; then
					mecatOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.mecat"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $mecatOutDir ] && [ -d "$mecatOutDir/4-fsa" ] && [ -s "$mecatOutDir/4-fsa/contigs.fasta" ]; then
						echo "(Step$step)Info: existing mecat out: $mecatOutDir/; skipping"
					else
						runMecat $mixFq $mecatOutDir
#						$mecatOutDir/4-fsa/contigs.fasta
					fi
				fi
				if [ $test2run_nextdenovo -eq 1 ]; then
					nextdenovoOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.nextdenovo"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $nextdenovoOutDir ]; then
						echo "(Step$step)Info: existing nextdenovo out: $nextdenovoOutDir/; skipping"
					else
						runNextdenovo $mixFq $nextdenovoOutDir
#						$nextdenovoOutDir/03.ctg_graph/nd.asm.fasta
					fi
				fi
				if [ $test2run_smartdenovo -eq 1 ]; then
					smartdenovoOutDir="$runDir/ALL${opt_pf}.mt${indcov}x.Rep${rep}.smartdenovo"
					cd $runDir
					if [ $keep_files -eq 1 ] && [ -d $smartdenovoOutDir ] && [ -s "$smartdenovoOutDir/${smartdenovoOutDir##*/}.dmo.cns" ]; then
						echo "(Step$step)Info: existing smartdenovo out: $smartdenovoOutDir/; skipping"
					else
						runSmartdenovo $mixFq $smartdenovoOutDir
#						"$smartdenovoOutDir/${smartdenovoOutDir##*/}.dmo.cns"
					fi
				fi
			done
		done
	fi
else
	echo "(Step$step)Warnings: Step7 skipped"
fi



exit 0

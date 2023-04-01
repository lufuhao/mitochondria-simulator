# mitochondia-simulator

These scripts are designed for the mtDNA assembly project

> bin

- [x] mtDNA sequencing simulation to optimize parameters for Illumina

- [x] assembly programs for 2nd /3rd - Generation sequencing

- [x] SPAdes was selected for its good performace, and can generate GFA file

> huilin

- [x] Contains some simulation scripts for excise before real data

- [x] Hui-Lin Hu was a Master student starting 2019

- [x] Hard working to learn bioinformatics from the very beginning

- [x] So be patient for some errors or traps. Aha...

> Illumina

- [x] Contains cmds for Illumina data assembly, annotation, and analysis

- [x] We wrote many custom scripts for data mining and plotting in 4 years

- [x] Scripts would be available on request

- [x] Be nice when you ask, or Huilin may deny your request upon her mood. Aha...

> Be gentle before its publication (Will let you know here)

- [x] Developer: Hui-Lin Hu (Master student, 2019)

- [x] Group Leader: Fu-Hao Lu





## Requirements

### ms.1.mason.reads.sh
    [x] Linux: gzip, mkdir, cd, perl

    [x] BWA

    [x] deepTools

    [x] Fastp

    [x] FastqTools

    [x] FastQC

    [x] FuhaoBin: fastq_checkid.pl

    [x] Mason

    [x] SAMtools

    [x] SeqKit



### 

## Descriptions


### 2nd-Generation simulation

```
ms.1.mason.reads.sh -rl 150 -cp 0.05,0.05 -mc 2000,1000,500,250,100,50 -rp 5 -d "/home/hpcusers/admin/lufuhao/mit/test" \
     -gn "/home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta" \
     -mt "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.mason.500bp.fa" \
     -ct "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrC.fasta"

mitochondria-simulator/bin/mt.1.sim.reads.sh -pf "illumina" -cp "0.05,0.05" -mc 500 -rp 3 -k -rl 150 -tl 51 -mm 10 -km 71,91,101 \
  -s "1,2,3,4,5,6,7" -sa -em 'lufuhao@henu.edu.cn' -fl 350 \
  -as "abyss,edena,ray,minia,spades" -dg -d "/home/hpcusers/admin/lufuhao/mit/test" \
  -gn "/home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta" \
  -mt "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.badread.15000bp.fa" \
  -ct "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrC.fasta" 

mitochondria-simulator/bin/mt.1.sim.reads.sh -pf "pacbio" -cp "0.05,0.05" -mc 500 -rp 3 -k -rl 15000,13000 -mm 10 \
   -s "1,2,3,4,5,6,7" -sa -em 'lufuhao@henu.edu.cn' \
  -as "canu,flye,mecat,nextdenovo,smartdenovo" -dg -d "/home/hpcusers/admin/lufuhao/mit/test" \
  -gn "/home/hpcshared/Databases/genomes/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta" \
  -mt "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrM.badread.15000bp.fa" \
  -ct "/home/hpcusers/master/huilin_hu/simulate/Triticum_cs.chrC.fasta"


```


## Author:

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

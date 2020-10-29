## Steps in producing the final fasta format paired end reads for use in alignment.
## A. Execute the primary filtering script, ffq.TS_PE_SE_RNA.pl
##   1. Produces a LOT of output files, having various uses.
##   2. We care about 4
##     a. r1.pe.fq.gz and r2.pe.fq.gz
##        these are "normal" read-pairs where no overlap was detected
##     b. r1.se.fq.gz
##        these are read-pairs that were interpreted to have an overlap.
##        each such pair has been merged into a single read, based on the overlap.
##     c. failedQuality.r1.fq.gz
##        these are the R1s for all pairs that failed due to low quality in either read
##
## B. Convert r1.pe.fq.gz and r2.pe.fq.gz to fasta files
##   1. new names: r1.pe.fa.gz, r2.pe.fa.gz
##
## C. Convert r1.se.fq.gz into "pseudo-pair" fasta files
##   1. r1.se.fa.gz is r1.se.fq.gz stripped of quality scores
##   2. r2.se.fa.gz is reverse-complements of r1.seq.fa.gz
##
## D. Rescue reads in failedQuality.r1.fq.gz that by themselves have sufficient quality
##   1. fQ_rescue.r1.seq.fa.gz holds the rescued reads, trimmed to include only high-quality basecalls
##   2. fQ_rescue.r2.seq.fa.gz holds the reverse complements of fQ_rescue.r1.seq.fa.gz
##
## E. Concatenate the 3 sets of paired reads and pseudo-paired reads into the final fasta files
##    that will be input to the aligner
##   1. r1.sepe.fa.gz
##   2. r2.sepe.fa.gz
##
## IMPORTANT NOTE!!!
## ffq.TS_PE_SE_RNA.pl assumes that the raw fastq.gz files will
##   match the format in the command below.
## if yours do not you MUST edit ffq.TS_PE_SE_RNA.pl accordingly!!
##
## The arguments to scripts/ffq.TS_PE_SE_RNA.pl must be given in the exact order shown
## A. ../../Samples/$1/*_R1_001.fastq.gz: the list of all R1 raw fastq
## B. ../../Samples/$1/*_R2_001.fastq.gz: the list of all R2 raw fastq
## C. ../filter_refdata/chrM.50mers.txt: the file containing mitochondrial N-mers for pre-filtering of mtRNA
## D. 50: the length N of mitochondrial N-mers to be used for pre-filtering
## E. 15: the minimum allowable Quality score for basecalls
## F. 22: the minimum allowable length of the remaining read after trimming low-Quality basecalls
## G. 0: the number of reads to skip at the beginning of each raw fastq file
## H. 1000000000: the number of reads to process from each raw fastq file
mkdir $1
cd $1

../scripts/ffq.TS_PE_SE_RNA.pl \
../../Samples/$1/*_R1_001.fastq.gz \
../../Samples/$1/*_R2_001.fastq.gz \
../filter_refdata/chrM.50mers.txt 50 \
15 22 \
0 1000000000 \
> ffq.TS_PE_SE_RNA.pl.out 2> ffq.TS_PE_SE_RNA.pl.log

pigz -d -p 3 --stdout < ./r1.pe.fq.gz | perl -nle '$km4 = $k++ % 4; if ($km4 == 0){substr($_,0,1)=">"} if (2 > $km4){print}' | pigz > ./r1.pe.fa.gz
pigz -d -p 3 --stdout < ./r2.pe.fq.gz | perl -nle '$km4 = $k++ % 4; if ($km4 == 0){substr($_,0,1)=">"} if (2 > $km4){print}' | pigz > ./r2.pe.fa.gz

pigz -d -p 3 --stdout < ./r1.se.fq.gz | perl -nle '$km4 = $k++ % 4; if ($km4 == 0){substr($_,0,1)=">"} if (2 > $km4){print}' | pigz > ./r1.se.fa.gz
pigz -d -p 3 --stdout < ./r1.se.fa.gz | perl -nle 'if (m/^>/){print}else{tr/ACGT/TGCA/; $rs=reverse $_; print $rs}' | pigz > ./r2.se.fa.gz

pigz -d -p 3 --stdout failedQuality.r1.fq.gz \
| perl -nle \
'$kmodz=($k++ % 4);
 if (3 == $kmodz){
   if (($good,$bad) = m/([F:]{40}[F:]+)(.*)$/){
     substr($h,0,1)=">";
     print $h;
     $lenGood = length $good;
     $lenBad = length $bad;
     $lenSeq = length $_;
     $goodStart = $lenSeq - $lenGood - $lenBad;
     print substr($s, $goodStart, $lenGood);
   }
 }
 elsif (0 == $kmodz) {
   $h=$_;
 }
 elsif (1 == $kmodz) {
   $s=$_;
 }
' \
| pigz -p 3 --stdout \
> fQ_rescue.r1.se.fa.gz

pigz -d -p 3 --stdout < ./fQ_rescue.r1.se.fa.gz | perl -nle 'if (m/^>/){print}else{tr/ACGT/TGCA/; $rs=reverse $_; print $rs}' | pigz > ./fQ_rescue.r2.se.fa.gz
pigz -d -p 3 --stdout ./r1.se.fa.gz ./r1.pe.fa.gz ./fQ_rescue.r1.se.fa.gz | pigz -p 3 > ./r1.sepe.fa.gz
pigz -d -p 3 --stdout ./r2.se.fa.gz ./r2.pe.fa.gz ./fQ_rescue.r2.se.fa.gz | pigz -p 3 > ./r2.sepe.fa.gz

cd ..

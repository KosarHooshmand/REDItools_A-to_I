# REDItools-Installation
git clone [REDItools](https://github.com/BioinfoUNIBA/REDItools)
cd REDItools/
python setup.py install

# Installatoion of associated packages:
## Pblat installation. need to be version 0.6 and 0.7

git clone [Pblat](https://github.com/icebert/pblat.git)

cd pblat/

make

## pysam version 0.17 installation:

pip install pysam==0.17

## samtools version 1.9 installation:

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2

tar -vxjf samtools-1.9.tar.bz2

cd samtools-1.9

make

## bcftools version 1.9 installation:

wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2

tar -vxjf bcftools-1.9.tar.bz2

cd bcftools-1.9

make

## bedtools version 2.28.0 installation

wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz


tar -zxvf bedtools-2.28.0.tar.gz

cd bedtools2

make

##########################################################################################################

# Downloading and organizing required data

##  ● Index the reference genome for REDItools

~/Apps/samtools-1.9/samtools faidx ~/Ref/GRCh38.p13.genome.fa

##  ● Create the nochr file for REDItools

grep ">" ~/Ref/GRCh38.p13.genome.fa  | awk '{if (substr($1,1,3)==">GL") print $2}' > nochr1


##  ● Download and unzip hg38 RefSeq annotations in .bed format for strand detection

wget https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_RefSeq.bed.gz

gunzip hg38_RefSeq.bed.gz

##  ● Download and unzip hg38 RefSeq annotations in .txt format for Gene symbols

wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

gunzip hg38.refGene.txt.gz

conda install -c bioconda ucsc-genepredtogtf

cut -f 2- hg38.refGene.txt | genePredToGtf -utr -source=hg38_refseq file stdin hg38-RefGene.gtf

sort -k1,1 -k4,4n hg38-RefGene.gtf > Sorted-hg38-RefGene.gtf 

##  ● Prepare RepeatMasker annotations for REDItools

http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

awk '{OFS="\t"} {print $6,"rmsk_hg38",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$13"\";"}' rmsk.txt > rmsk38.gtf

sort -k1,1 -k4,4n rmsk38.gtf > rmsk38.sorted.gtf

bgzip rmsk38.sorted.gtf

tabix -p gff rmsk38.sorted.gtf.gz

##  ● Prepare dbSNP annotations for REDItools

http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz

awk '{OFS="\t"} {if ($11=="genomic" && $12=="single") print $2,"ucsc_snp151_hg38","snp",$4,$4,".",$7,".","gene_id \""$5"\"; transcript_id \""$5"\";"}' snp151.txt > snp151.gtf

sort -k1,1 -k4,4n snp151.gtf > snp151.sorted.gtf

bgzip snp151.sorted.gtf

tabix -p gff snp151.sorted.gtf.gz

##  ● Prepare REDIportal annotations for REDItools and extract recoding events

http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz

### remove the header:

tail -n +2 TABLE1_hg38.txt > TABLE1_hg38_without.txt

awk '{OFS="\t"} {sum+=1; print $1,"rediportal","ed",$2,$2,".",$5,".","gene_id \""sum"\"; transcript_id\""sum"\";"}' TABLE1_hg38_without.txt > atlas38.gtf

sort -V -k1,1 -k4,4n atlas38.gtf > sorted_atlas38.gtf

bgzip sorted_atlas38.gtf

tabix -p gff sorted_atlas38.gtf.gz

python2 ~/Apps/REDItools/accessory/rediportal2recoding.py TABLE1_hg38_without.txt > atlas38_recoding.gff

sort -V -k1,1 -k4,4n atlas38_recoding.gff > sorted_atlas38_recoding.gff

bgzip sorted_atlas38_recoding.gff

tabix -p gff sorted_atlas38_recoding.gff.gz

##  ● Prepare REDIportal annotations For REDItoolKnown.py

sort -k1,1 -k2,2n TABLE1_hg38_without.txt > sorted-TABLE1_hg38_without.txt

bgzip sorted-TABLE1_hg38_without.txt

python2 ~/Apps/REDItools/main/REDItoolKnown.py -i STAR_Alignment/*.bam -f ~/Ref/GRCh38.p13.genome.fa -o ./known-edits -l ./sorted-TABLE1_hg38_without.txt.gz

##  ● Prepare splice sites annotations for REDItools

https://github.com/juliangehring/GMAP-GSNAP

In case GMAT is a typo, GMAP and GSNAP can be found here:

http://research-pub.gene.com/gmap/

wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2021-08-25.tar.gz

tar -xvzf gmap-gsnap-2021-08-25.tar.gz

cd gmap-gsnap-2021-08-25

./configure --prefix /home/kosar/Apps 

make 

make check

make install

export PATH=/home/kosar/Apps/gmap-gsnap-2021-08-25:$PATh


### If you have a GTF file, you can use the included programs gtf_splicesites and gtf_introns like this:

cat <gtf file> | /gmap-2021-08-25/util/gtf_splicesites > foo.splicesites

cat <gtf file> | /gmap-2021-08-25/util/gtf_introns > foo.introns


### Second, if you retrieve an alignment tracks from UCSC:

 ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/ene.txt.gz

ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

### you can process this track using the included program psl_splicesites or psl_introns, like this:

gunzip -c refGene.txt.gz | /gmap-2021-08-25/util/psl_splicesites -s 1 > foo.splicesites

gunzip -c refGene.txt.gz | /gmap-2021-08-25/util/psl_introns -s 1 > foo.introns

mawk -F" " '{split($2,a,":"); split(a[2],b,"."); if (b[1]>b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}' foo.splicesites > mysplicesites.ss
  
  
  
  

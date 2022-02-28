# ⚙️Technologies & Tools
![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

## Installatoion of REDItools and the associated packages:

* REDItools installation:

git clone [REDItools](https://github.com/BioinfoUNIBA/REDItools)

`cd REDItools/`

`python setup.py install`

* Pblat version 0.6 and 0.7 installation:

git clone [Pblat](https://github.com/icebert/pblat.git)

`cd pblat/`

`make`

* pysam version 0.17 installation:

`pip install pysam==0.17`

* samtools version 1.9 installation:

wget [samtools](https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2) 

`tar -vxjf samtools-1.9.tar.bz2`

`cd samtools-1.9`

`make`

* bcftools version 1.9 installation:

wget [bcftools](https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2) 

`tar -vxjf bcftools-1.9.tar.bz2`

`cd bcftools-1.9`

`make`

* bedtools version 2.28.0 installation

wget [bedtools](https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz)

`tar -zxvf bedtools-2.28.0.tar.gz`

`cd bedtools2`

`make`

# :pushpin: Downloading and organizing required data

* Index the reference genome for REDItools

`~/Apps/samtools-1.9/samtools faidx ~/Ref/GRCh38.p13.genome.fa`

* Create the nochr file for REDItools

`grep ">" ~/Ref/GRCh38.p13.genome.fa  | awk '{if (substr($1,1,3)==">GL") print $2}' > nochr1`


* Download and unzip hg38 RefSeq annotations in .bed format for strand detection

wget -c -O [RefSeq](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)

`gunzip hg38.refGene.txt.gz`

`conda install -c bioconda ucsc-genepredtogtf`

`cut -f 2- hg38.refGene.txt | genePredToGtf -utr -source=hg38_refseq file stdin hg38-RefGene.gtf`

`sort -k1,1 -k4,4n hg38-RefGene.gtf > Sorted-hg38-RefGene.gtf`

* Prepare RepeatMasker annotations for REDItools

wget [rmsk](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz)

`awk '{OFS="\t"} {print $6,"rmsk_hg38",$12,$7+1,$8,".",$10,".","gene_id \""$11"\"; transcript_id \""$13"\";"}' rmsk.txt > rmsk38.gtf`

`sort -k1,1 -k4,4n rmsk38.gtf > rmsk38.sorted.gtf`

`bgzip rmsk38.sorted.gtf`

`tabix -p gff rmsk38.sorted.gtf.gz`

* Prepare dbSNP annotations for REDItools

wget [dbSNP](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz)

`awk '{OFS="\t"} {if ($11=="genomic" && $12=="single") print $2,"ucsc_snp151_hg38","snp",$4,$4,".",$7,".","gene_id \""$5"\"; transcript_id \""$5"\";"}' snp151.txt > snp151.gtf`

`sort -k1,1 -k4,4n snp151.gtf > snp151.sorted.gtf`

`bgzip snp151.sorted.gtf`

`tabix -p gff snp151.sorted.gtf.gz`

* Prepare REDIportal annotations for REDItools and extract recoding events

wget [REDIprotal](http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz)

`tail -n +2 TABLE1_hg38.txt > TABLE1_hg38_without.txt`

`awk '{OFS="\t"} {sum+=1; print $1,"rediportal","ed",$2,$2,".",$5,".","gene_id \""sum"\"; transcript_id\""sum"\";"}' TABLE1_hg38_without.txt > atlas38.gtf`

`sort -V -k1,1 -k4,4n atlas38.gtf > sorted_atlas38.gtf`

`bgzip sorted_atlas38.gtf`

`tabix -p gff sorted_atlas38.gtf.gz`

* Prepare REDIportal annotations For REDItoolKnown.py

`sort -k1,1 -k2,2n TABLE1_hg38_without.txt > sorted-TABLE1_hg38_without.txt`

`bgzip sorted-TABLE1_hg38_without.txt`

* Prepare splice sites annotations for REDItools

[GMAP-GSNAP]( https://github.com/juliangehring/GMAP-GSNAP)

[gmap](http://research-pub.gene.com/gmap/) 

wget [GMAP-GSNAP](http://research-pub.gene.com/gmap/src/gmap-gsnap-2021-08-25.tar.gz)

`tar -xvzf gmap-gsnap-2021-08-25.tar.gz`

`cd gmap-gsnap-2021-08-25`

`./configure --prefix ~/Apps`

`make`

`make check`

`make install`

`export PATH=~/Apps/gmap-gsnap-2021-08-25:$PATh`


### If you have a GTF file, you can use the included programs gtf_splicesites and gtf_introns like this:

`cat <gtf file> | /gmap-2021-08-25/util/gtf_splicesites > foo.splicesites`

`cat <gtf file> | /gmap-2021-08-25/util/gtf_introns > foo.introns`

`mawk -F" " '{split($2,a,":"); split(a[2],b,"."); if (b[1]>b[3]) print a[1],b[3],b[1],toupper(substr($3,1,1)),"-"; else print a[1],b[1],b[3],toupper(substr($3,1,1)),"+"}' foo.splicesites > mysplicesites.ss`

# :writing_hand: RNA Editing Detection 

## 1a.Detect known RNA edit sites

`python2 ~/Apps/REDItools/main/REDItoolKnown.py -i ../STAR_Alignment/*.bam -f ~/kosar/Ref/GRCh38.p13.genome.fa -o ./STARcheck -l Sorted-TABLE1_hg38.tab`
  
## 1b. Detect all potential DNA–RNA variants

`python2 ~/Apps/REDItools/main/REDItoolDnaRna.py -i ../STAR_Alignment/*.bam -j ../BWA-mem2_Alignment/*.bam -o NEW_RNAEdits_picard -f ~/Ref/GRCh38.p13.genome.fa -t10 -c1,1 -m30,255 -v1 -q30,30 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S`

## 2. Exclude invariant positions as well as positions not supported by ≥10 WGS reads:

  `awk '{FS="\t"} {if ($8!="-" && $10>=10 && $13=="-") print}' outTable_443931662 > outTable_443931662.out`

## 3. Annotate positions using RepeatMasker and dbSNP annotations:

`python2 ~/Apps/REDItools/accessory/AnnotateTable.py -a ../../RMSK/rmsk38.sorted.gtf.gz -n rmsk -i outTable_443931662.out -o Table_443931662.out.out.rmsk -u`

`python2 ~/Apps/REDItools/accessory/AnnotateTable.py -a ../../SNP151/snp151.sorted.gtf.gz -n snp151 -i Table_443931662.out.out.rmsk -o Table_443931662.out.rmsk.snp -u`

## 4. Create a first set of positions selecting sites supported by at least five RNAseq reads and a single mismatch:

`python2 ~/Apps/REDItools/accessory/selectPositions.py -i Table_443931662.out.rmsk.snp -c 5 -v 1 -f 0.0 -o Table_443931662.out.rmsk.snp.sel1`

## 5. Create a second set of positions selecting sites supported by ≥10 RNAseq reads, three mismatches and minimum editing frequency of 0.1:

`python2 ~/Apps/REDItools/accessory/selectPositions.py -i Table_443931662.out.rmsk.snp -c 10 -v 3 -f 0.1 -o Table_443931662.out.rmsk.snp.sel2`

## 6. Select ALU sites from the first set of positions:
  
`awk '{FS="\t"} {if ($1!="chrM" && substr($16,1,3)=="Alu" && $17=="-"&& $8!="-" && $10>=10 && $13=="-") print}' Table_443931662.out.rmsk.snp.sel1 > Table_443931662.out.rmsk.snp.alu`

## 7. Select REP NON ALU sites from the second set of positions, excluding sites in Simple repeats or Low complexity regions:

`awk '{FS="\t"} {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15!="-" && $15!="Simple_repeat" && $15!="Low_complexity" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' Table_443931662.out.rmsk.snp.sel2 > Table_443931662.out.rmsk.snp.nonalu`

## 8. Select NON REP sites from the second set of positions:

`awk '{FS="\t"} {if ($1!="chrM" && substr($16,1,3)!="Alu" && $15=="-" && $17=="-" && $8!="-" && $10>=10 && $14<=0.05 && $9>=0.1) print}' Table_443931662.out.rmsk.snp.sel2 > Table_443931662.out.rmsk.snp.nonrep`

## 9. Annotate ALU, REP NON ALU and NON REP sites using known editing events from REDIportal:
  
`python2 ~/Apps/REDItools/accessory/AnnotateTable.py -a ../../REDIPortal/sorted_atlas38.gtf.gz -n ed -k R -c 1 -i Table_443931662.out.rmsk.snp.alu -o Table_443931662.out.rmsk.snp.alu.ed -u`
  
`python2 ~/Apps/REDItools/accessory/AnnotateTable.py -a ../../REDIPortal/sorted_atlas38.gtf.gz -n ed -k R -c 1 -i Table_443931662.out.rmsk.snp.nonalu -o Table_443931662.out.rmsk.snp.nonalu.ed -u`
  
`python2 ~/Apps/REDItools/accessory/AnnotateTable.py -a ../../REDIPortal/sorted_atlas38.gtf.gz -n ed -k R -c 1 -i Table_443931662.out.rmsk.snp.nonrep -o Table_443931662.out.rmsk.snp.nonrep.ed -u`

## 10. Extract known editing events from ALU, REP NON ALU and NON REP sites:
  
`mv Table_443931662.out.rmsk.snp.alu.ed alu`
  
`mv Table_443931662.out.rmsk.snp.nonalu.ed nonalu`
  
`mv Table_443931662.out.rmsk.snp.nonrep.ed nonrep`
  
`cat alu nonalu nonrep > alu-nonalu-nonrep`
  
`awk '{FS="\t"} {if ($19=="ed") print}' alu-nonalu-nonrep > knownEditing`
 
## 11. Convert editing candidates in REP NON ALU and NON REP sites in GFF format for further filtering:
  
`cat nonalu nonrep > nonalu-nonrep`
  
`awk '{FS="\t"} {if ($19!="ed") print}' nonalu-nonrep > pos.txt`
  
`python2 ~/Apps/REDItools/accessory/TableToGFF.py -i pos.txt -s -t -o pos.gff`

## 12. Convert editing candidates in ALU sites in GFF format for further filtering:
  
`awk '{FS="\t"} {if ($19!="ed") print}' alu > posalu.txt`

`python2 ~/Apps/REDItools/accessory/TableToGFF.py -i posalu.txt -s -t -o posalu.gff`

## 13. Launch REDItoolDnaRna.py on ALU sites using stringent criteria to recover potential editing candidates:

`python2.7 ~/Apps/REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 4 -i ../STAR_Alignment/CGND-HRA-00072_AGCACCTC_ACB20PANXX_L005_001_marked_duplicates.bam -f ~/Ref/GRCh38.p13.genome.fa -c 5,5 -q 30,30 -m255,255 -O 5,5 -p -u -a 11-6 -l -v 1 -n 0.0 -e -T RNA72-DNA191-posalu.sorted.gff.gz -w  ../REDItools2/Splicesites/mysplicesites.ss -k ../REDItools2/Genome_hg38/nochr -R -o RNA72-DNA191-firstalu`

## 14. Launch REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria to recover RNAseq reads harboring reference mismatches:
  
`python2.7 ~/Apps/REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 4 -i ../STAR_Alignment/CGND-HRA-00072_AGCACCTC_ACB20PANXX_L005_001_marked_duplicates.bam -f ~/Ref/GRCh38.p13.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T RNA72-DNA191-pos.sorted.gff.gz -w ../REDItools2/Splicesites/mysplicesites.sss -k ../REDItools2/Genome_hg38/nochr --reads -R --addP -o RNA72-DNA191-first`

## 15. Launch pblat on RNAseq reads harboring reference mismatches from Step 22 and select multi-mapping reads:
  
`~/Apps/pblat/pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 ~/Ref/GRCh38.p13.genome.fa first/DnaRna_51144481/outReads_51144481 reads.psl
python2 ~/Apps/REDItools/accessory/readPsl.py reads.psl badreads.txt`

## 16. Extract RNAseq reads harboring reference mismatches from Step 14 and remove duplicates:
  
`sort -k1,1 -k2,2n -k3,3n first/DnaRna_51144481/outPosReads_51144481 | mergeBed > bed`
  
`~/Apps/samtools-1.9/samtools view -@ 4 -L bed -h -b ../STAR_Alignment/*.bam > *_bed.bam`
  
`~/Apps/samtools-1.9/samtools sort -@ 4 -n *_bed.bam -o *_bed_ns.bam`
  
`~/Apps/samtools-1.9/samtools fixmate -@ 4 -m *_bed_ns.bam *_bed_ns_fx.bam`
  
`~/Apps/samtools-1.9/samtools sort -@ 4 *_bed_ns_fx.bam -o *_bed_ns_fx_st.bam`
  
`~/Apps/samtools-1.9/samtools markdup -r -@ 4 *_bed_ns_fx_st.bam *_bed_dedup.bam`
  
`~/Apps/samtools-1.9/samtools index *_bed_dedup.bam`

## 17. Re-run REDItoolDnaRna.py on REP NON ALU and NON REP sites using stringent criteria, deduplicated reads and mis-mapping info:
:heavy_exclamation_mark: Our bam files does not have the indels info so we need to remove --rmIndels from the command line
  
`python2 ~/Apps/REDItools/main/REDItoolDnaRna.py -s 2 -g 2 -S -t 4 -i *_bed_dedup.bam -f .~/Ref/GRCh38.p13.genome.fa -c 10,10 -q 30,30 -m 255,255 -O 5,5 -p -u -a 11-6 -l -v 3 -n 0.1 -e -T pos.sorted.gff.gz -w ../Splicesites/mysplicesites.ss -R -k ../Genome_hg38/nochr -b badreads.txt --rmIndels -o second python2 ~/Apps/REDItools/main/REDItoolDnaRna.py -i ../STAR_Alignment/*.bam -j ../BWA-mem2_Alignment/*.bam -o NEW_RNAEdits_picard -f ~/Ref/GRCh38.p13.genome.fa -t10 -c1,1 -m30,255 -v1 -q30,30 -e -n0.0 -N0.0 -u -l -p -s2 -g2 -S` 

## 18. Collect filtered ALU, REP NON ALU and NON REP sites:
:heavy_exclamation_mark: Run the command in the directory which include the outputs related to Second, firstalu and knownEditing.txt per individual (sample).

`python2 ~/Apps/REDItools/NPscripts/collect_editing_candidates.py`

`sort -k1,1 -k2,2n editing.txt > editing_sorted.txt`

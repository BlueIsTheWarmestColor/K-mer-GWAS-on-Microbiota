The pipeline of k-mer GWAS on maize leaf microbiota
==========================================================

These scripts show how to analyse the microbiome-realated k-mers in k-mer GWAS results.

1.Extract unmapped reads of maize WGS data
---------------------------------------------
/kmerGWAS/gwas_test/meta_genome/extract_unmapped_reads.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate david
export LD_LIBRARY_PATH=/data05/bxin/kmerGWAS/my_lib:$LD_LIBRARY_PATH
library_dir="/data05/bxin/softwares/kmerGWAS/" # The path to the directory with the library
base_dir="/data05/bxin/kmerGWAS" # This will be our working directory
operation="/data05/bxin/softwares/fetch_reads_with_kmers/fetch_reads"
bowtie2="/data05/bxin/softwares/bowtie2-2.4.4-linux-x86_64/bowtie2"
DTT_contig="/data05/bxin/kmerGWAS/gwas_test/DTT_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/contigs.fasta"
DTA_contig="/data05/bxin/kmerGWAS/gwas_test/DTA_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/contigs.fasta"
DTS_contig="/data05/bxin/kmerGWAS/gwas_test/DTS_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/contigs.fasta"
EH_contig="/data05/bxin/kmerGWAS/gwas_test/EH_redo_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge/contigs.fasta"
KRN_contig="/data05/bxin/kmerGWAS/gwas_test/KRN_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/contigs.fasta"
KT_contig="/data05/bxin/kmerGWAS/gwas_test/KT_maf0.01/from_kmer_get_loci/step5_merge_reads_assembly/5per_merge_redo/contigs.fasta"
META_genome="/data05/bxin/kmerGWAS/gwas_test/meta_genome/meta_genome_reference/ref_file_origin/all_genome.fna"
THREADS="4"
 
for name in `cat ../accessions_id|tail -n +16`
do 
echo "Working on sample  : $name"
 
        # 1. Download sequence files
        bam=obs://mem-result/$name/$name\_mem_rd.sort.bam
        prefix=$name\_mem_rd
        cd /evs/ &&/bin/rm -rf /evs/* &&  mkdir $name && cd $name 
        date
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >  obs_bam_download.log  2> speed1  
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 2>> speed1 
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 2>> speed1 
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 
        obsutil cp $bam /evs/$name/$prefix.sort.bam  -f  -u  -j=3  -vmd5 -vlength >>  obs_bam_download.log 
        date
        date
        samtools view -@ $THREADS -bh -f 4 $prefix.sort.bam \
            -o $prefix.sort.unmapped.bam  && echo 'unmapped bam finished'
        date
        mkdir /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name
        date
        samtools fastq -@ $THREADS $prefix.sort.unmapped.bam \
            -1 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$prefix\_unmapped_1.fq \
            -2 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$prefix\_unmapped_2.fq \
            -0 /dev/null -s /dev/null -n  && echo 'fastq finished'
        date
        /bin/rm $prefix.sort.bam
        /bin/rm $prefix.sort.unmapped.bam
done
```


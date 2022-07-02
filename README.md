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

2.Map unmapped reads to contigs signicifant k-mer-realated
-------------------------------------------------------------
/kmerGWAS/gwas_test/meta_genome/unmapped_reads_remapping_contigs_3367.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate py_37_simon
bowtie2="/data05/bxin/softwares/bowtie2-2.4.4-linux-x86_64/bowtie2"
DTT_contig="/data05/bxin/kmerGWAS/gwas_test/meta_genome/DTT/needed_contigs_final.fa"
DTA_contig="/data05/bxin/kmerGWAS/gwas_test/meta_genome/DTA/needed_contigs_final.fa"
DTS_contig="/data05/bxin/kmerGWAS/gwas_test/meta_genome/DTS/needed_contigs_final.fa"
EH_contig="/data05/bxin/kmerGWAS/gwas_test/meta_genome/EH/needed_contigs_final.fa"
KRN_contig="/data05/bxin/kmerGWAS/gwas_test/meta_genome/KRN/needed_contigs_final.fa"
KT_contig="/data05/bxin/kmerGWAS/gwas_test/meta_genome/KT/needed_contigs_final.fa"
THREADS="4"
 
for name in `cat ../accessions_id`
do 
echo "Working on sample  : $name"
mkdir /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/mapping_contig/${name}
$bowtie2 --threads $THREADS --very-sensitive-local -x $DTT_contig -1 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_1.fq.gz -2 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_2.fq.gz |samtools view -bS -q 20 - | samtools sort - -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/mapping_contig/${name}/DTT_remapping.q20.sorted.bam && echo 'DTT mapping finished'
$bowtie2 --threads $THREADS --very-sensitive-local -x $DTA_contig -1 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_1.fq.gz -2 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_2.fq.gz |samtools view -bS -q 20 - | samtools sort - -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/mapping_contig/${name}/DTA_remapping.q20.sorted.bam && echo 'DTA mapping finished'
$bowtie2 --threads $THREADS --very-sensitive-local -x $DTS_contig -1 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_1.fq.gz -2 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_2.fq.gz |samtools view -bS -q 20 - | samtools sort - -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/mapping_contig/${name}/DTS_remapping.q20.sorted.bam && echo 'DTS mapping finished'
$bowtie2 --threads $THREADS --very-sensitive-local -x $EH_contig -1 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_1.fq.gz -2 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_2.fq.gz |samtools view -bS -q 20 - | samtools sort - -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/mapping_contig/${name}/EH_remapping.q20.sorted.bam && echo 'EH mapping finished'
$bowtie2 --threads $THREADS --very-sensitive-local -x $KRN_contig -1 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_1.fq.gz -2 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_2.fq.gz |samtools view -bS -q 20 - | samtools sort - -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/mapping_contig/${name}/KRN_remapping.q20.sorted.bam && echo 'KRN mapping finished'
$bowtie2 --threads $THREADS --very-sensitive-local -x $KT_contig -1 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_1.fq.gz -2 /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/unmapped_reads/$name/$name\_unmapped_2.fq.gz |samtools view -bS -q 20 - | samtools sort - -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/3367/mapping_contig/${name}/KT_remapping.q20.sorted.bam && echo 'KT mapping finished'
done
```

3.Mutation calling
------------------------
In this step, the different genotypes of specific microbiome in maize leaf are identified

#### Step1 all mutation on the contigs
/kmerGWAS/gwas_test/meta_genome/bcftools_3367.sh
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
conda activate py_37_simon
KRN_all_bam_list=KRN_3367_bam
KT_all_bam_list=KT_3367_bam
DTA_all_bam_list=DTA_3367_bam
DTS_all_bam_list=DTS_3367_bam
DTT_all_bam_list=DTT_3367_bam
EH_all_bam_list=EH_3367_bam
DTT_contigs="/data05/bxin/kmerGWAS/gwas_test/meta_genome/DTT/needed_contigs_final.fa"
DTA_contigs="/data05/bxin/kmerGWAS/gwas_test/meta_genome/DTA/needed_contigs_final.fa"
DTS_contigs="/data05/bxin/kmerGWAS/gwas_test/meta_genome/DTS/needed_contigs_final.fa"
EH_contigs="/data05/bxin/kmerGWAS/gwas_test/meta_genome/EH/needed_contigs_final.fa"
KRN_contigs="/data05/bxin/kmerGWAS/gwas_test/meta_genome/KRN/needed_contigs_final.fa"
KT_contigs="/data05/bxin/kmerGWAS/gwas_test/meta_genome/KT/needed_contigs_final.fa"
bcftools mpileup --threads 32 -b $KRN_all_bam_list --fasta-ref $KRN_contigs > /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/KRN/KRN.vcf
bcftools call --threads 32 /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/KRN/KRN.vcf -mv -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/KRN/KRN_variants.vcf
bcftools mpileup --threads 32 -b $EH_all_bam_list --fasta-ref $EH_contigs > /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH.vcf
bcftools call --threads 32 /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH.vcf -mv -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH_variants.vcf
bcftools mpileup --threads 32 -b $KT_all_bam_list --fasta-ref $KT_contigs > /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/KT/KT.vcf
bcftools call --threads 32 /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/KT/KT.vcf -mv -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/KT/KT_variants.vcf
bcftools mpileup --threads 32 -b $DTA_all_bam_list --fasta-ref $DTA_contigs > /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTA/DTA.vcf
bcftools call --threads 32 /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTA/DTA.vcf -mv -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTA/DTA_variants.vcf
bcftools mpileup --threads 32 -b $DTS_all_bam_list --fasta-ref $DTS_contigs > /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTS/DTS.vcf
bcftools call --threads 32 /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTS/DTS.vcf -mv -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTS/DTS_variants.vcf
bcftools mpileup --threads 32 -b $DTT_all_bam_list --fasta-ref $DTT_contigs > /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTT/DTT.vcf
bcftools call --threads 32 /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTT/DTT.vcf -mv -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTT/DTT_variants.vcf
```
#### Step2 filter the mutation we interested
```Bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
 
bgzip /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTS/DTS_remapping_variants.vcf
bgzip /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH_remapping_variants.vcf
bcftools filter -T kmer_region_all_needed /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH_remapping_variants.vcf.gz -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH_remapping_variants_kmer_region.vcf
bcftools filter -T kmer_region_all_needed /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTS/DTS_remapping_variants.vcf.gz -o /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTS/DTS_remapping_variants_kmer_region.vcf

conda activate Huang_population
 
cd /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/DTS/
vcftools --minQ 50 --vcf DTS_remapping_variants_kmer_region.vcf --out DTS_remapping_variants_kmer_region_filter_QUAL.vcf --recode --recode-INFO-all
cd /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/
vcftools --minQ 50 --vcf EH_remapping_variants_kmer_region.vcf --out EH_remapping_variants_kmer_region_filter_QUAL.vcf --recode --recode-INFO-all

gatk SelectVariants \
    -V /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH_remapping_variants_kmer_region_filter_QUAL.vcf.recode.vcf \
    -restrict-alleles-to BIALLELIC \
    -O /data05/bxin/kmerGWAS/gwas_test/meta_genome/all_accessions_mutation/EH/EH_remapping_variants_kmer_region_filter_QUAL_nomultiallelic.vcf.recode.vcf
```


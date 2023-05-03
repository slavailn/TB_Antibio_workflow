## Introduction 
Workflow that attempts to recreate the detection of the variants conferring antibiotic resistance in M. tuberculosis implemented as a web application at https://bioinf.fz-borstel.de/mchips/phyresse/

Supporting publication for the web app: https://journals.asm.org/doi/10.1128/JCM.00025-15

## Installing the necessary software with Conda
1. Install Miniconda for Linux
Download the installer
```
 wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
 bash Miniconda3-latest-Linux-x86_64.sh
```
Run the installer and follow the prompts.
Restart the session for the changes to take effect.
Create new conda environment names 'varcall'.

```
conda create -n varcall
```
Activate the newly created environment.

```
conda activate varcall
```
2. Install SRA toolkit
```
 conda install -c bioconda sra-tools
```

3. Install FastQC.
```
conda install -c bioconda fastqc
```
4. Install QualiMap
```
conda install -c bioconda qualimap
```

5. Install TrimGalore!
```
conda install -c bioconda trim-galore```

```
6. Install BWA (aligner).
```
conda install -c bioconda bwa
```
7. Install the variant caller (freebayes).
```
conda install -c bioconda freebayes
```
8. Install samtools
```
conda install -c bioconda samtools
```
Conda install sometimes results in a bug: Samtools shared library libcrypto.so.1.0.0 not found

In this case, update samtools to the latest build:
```
conda install -c bioconda samtools=1.9 --force-reinstall
```

## Obtaining sample data from GEO
Download a sample M.tuberculosis sequencing file from GEO.
Link to the file: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR17227596
Related project: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL31090

1) First, download SRA file with Aspera connect from SRAtoolkit. 
```
 prefetch -v prefetch -v SRR17227596
```
This command with create ncbi/ directory with SRA file downloaded from GEO.

Now we have to convert SRA file to fastq.

```
 fastq-dump ncbi/public/sra/SRR17227596.sra
```
This will create a fastq file *SRR17227596.fastq* in the working directory. 
I'm going to take top 200000 reads to reduce computational time.

```
head -n 800000 SRR17227596.fastq > sample.fastq
mkdir fastq # create directory to hold fastq files
mv sample.fastq fastq # move this file to fastq/ directory 
```
This command will take top 800,000 lines from SRR17578862.fastq and direct them into a new *sample.fastq* file. Since one read corresponds to 4 lines in еру fastq file this number translated into 200,000 reads.

2) Run initial quality control with FastQC.
Related links.
- Software: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Documentation: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/
- Training videos: https://www.youtube.com/watch?v=lUk5Ju3vCDM
  https://www.youtube.com/watch?v=bz93ReOv87Y
- Exampled of good and bad data: https://sequencing.qcfail.com/

Run FastQC with *sample.fastq* file:
```
fastqc sample.fastq
```
Examine the results and make the decision regarding trimming. Trimming is required when when the reads are of poor quality and/or in the presence of adapter contamination.

Trim the reads with TrimGalore!
```
 trim_galore -j 4 --length 18 --illumina --fastqc -q 30 sample.fastq
```
TrimGalore options explained:

  * -j 4 - use 4 processor cores in parallel
  * --length 18 - remove sequences that became shorter than 18 bp after trimming
  * --illumina - remove illumina adapters
  * --fastqc - generate FastQC report
  * -q 30 - remove low quality bases below Phred score of 30

With these options TrimGalore will generate trimmed fastq file, trimming report and fastqc report files. Examine FastQC report again after trimming. 

3) Mapping the reads to the reference genome.

BWA is a popular choice of an aligner for variant calling http://bio-bwa.sourceforge.net/. The reference genome used in PhyReSE is Mycobacterium tuberculosis H37Rv NC_000962.3 available from GenBank https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3?report=fasta

Downloaded M. tuberculosis genome and renamed the file to *Mtubersulosis.fasta*.

Genome needs to be indexed before the alignment:
```
 bwa index Mtuberculosis.fasta
```
Indexing will generate the following files that comprise the index:
Mtuberculosis.fasta.amb
Mtuberculosis.fasta.ann
Mtuberculosis.fasta.bwt
Mtuberculosis.fasta.pac
Mtuberculosis.fasta.sa

Map the reads to the genome:
```
 bwa mem -t 4 Mtuberculosis.fasta ../fastq/sample_trimmed.fq > sample.sam
```
Explanation:

  * -t 4 - use 4 processor cores
  * Mtuberculosis.fasta - index prefix, the working directory must contain all of the index files
  * ../fastq/sample_trimmed.fq - single-end reads (trimmed)
  * \> sample.sam - redirect the the output in SAM format to a new file
  
Convert SAM to BAM with samtools, sort and index:
```
samtools view -bS sample.sam > sample.bam
samtools sort sample.bam > sample.sorted.bam
samtools index sample.sorted.bam
```
Analyze quality of the alignment with *qualimap* http://qualimap.conesalab.org/doc_html/index.html

```
qualimap bamqc -bam sample.sorted.bam -outdir qualimap_results
```
View *qualimapReport.html* in qualimap_results/ directory to see various alignment results. In this case 87.54% of reads were mapped to the genome. Duplication rate was 13.57%. Mean coverage was 3.22 with the standard deviation of 36.17. Low coverage is understandable since we used only about 200,000 for the alignment. 

4) Call variants.

Mark duplicate reads and index marked BAM file:
```
picard MarkDuplicates I=sample.sorted.bam O=sample.dedup.bam M=marked_dup_metrix.txt

samtools index sample.dedup.bam
```
This command takes sorted BAM file and outputs another BAM file where PCR and optical duplicates are marked. This command does not remove the duplicate reads, only flags them as such within the newly created BAM file.

Call variants with freebayes:
```
freebayes -f ../genome/Mtuberculosis.fasta sample.dedup.bam --ploidy 1 > var.vcf

```
Explanation:

* -f ../genome/Mtuberculosis.fasta - point freebayes to the reference genome
* sample.dedup.bam - input should be sorted, indexed, deduplicated BAM file
* --ploidy 1 - ploidy is 1 for bacteria
* var.vcf - the output is a file containing variants in VCF format, this is a tab delimited text file containing variants. VCF specification can be found here: https://en.wikipedia.org/wiki/Variant_Call_Format

Filter resulting variant by quality:

```
vcffilter -f "QUAL > 20" var.vcf > filtered.vcf
```
This command will retain variants with quality over 20 on Phred scale, which translates in 0.01 probability of not being polymorphic or 0.99 probability of true polymorphism at this position. 

Find variants present the master file with antibiotic resistance SNVs and those obtained in our analysis. We use an ad-hoc perl script to match the genomic positions, reference and alternative bases and output matching variants as a tab delimited file *matching_variants.txt*.

```
./compare_to_master.pl nph-phyresse.csv filtered.vcf
```

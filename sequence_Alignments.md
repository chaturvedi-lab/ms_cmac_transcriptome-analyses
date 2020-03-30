This is the pipeline for doing sequence alignment for paired end reads generated for RNA sequencing data. The steps here show quality checks of raw reads, sequence alignment, and creating a count matrix for the differential expression analyses.

**1. FastQC**

I used FastQC to check quality of sequences and just visualized the quality metrics.

Here is the bash script I used:

```bash
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J fastqc

module load fastqc

##subsitute path to directory for your fastq files ###
for f in /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/RNAseq20151015/*.fastq; do fastqc –o  "./" $f; done

```

Fastqc outputs two sets of files for each sample: 1). html file 2). Zip file with all the statistics

I used the following bash code to compile the summaries and statistics from the zip folders for each sample. The output summaries are saved in the folder fq_aggregate.

```bash
# Run this script in a directory containing zip files from fastqc. It aggregates images of each type in individual folders
# So looking across data is quick. Output will be saved in the folder: fq_aggregate

zips=`ls *.zip`

for i in $zips; do
    unzip -o $i &>/dev/null;
done

fastq_folders=${zips/.zip/}

rm -rf fq_aggregated # Remove aggregate folder if present
mkdir fq_aggregated

# Rename Files within each using folder name.
for folder in $fastq_folders; do
    folder=${folder%.*}
    img_files=`ls ${folder}/Images/*png`;
    for img in $img_files; do
        img_name=$(basename "$img");
        img_name=${img_name%.*}
        new_name=${folder};
        mkdir -p fq_aggregated/${img_name};
        mv $img fq_aggregated/${img_name}/${folder/_fastqc/}.png;
    done;
done;


# Concatenate Summaries
for folder in $fastq_folders; do
    folder=${folder%.*}
    cat ${folder}/summary.txt >> fq_aggregated/summary.txt
done;

# Concatenate Statistics
for folder in $fastq_folders; do
    folder=${folder%.*}
    head -n 10 ${folder}/fastqc_data.txt | tail -n 7 | awk -v f=${folder/_fastqc/} '{ print $0 "\t" f }' >> fq_aggregated/statistics.txt
    rm -rf ${folder}
done
```

The results show some summary stats as failed. I did the next step and then reran fastqc.

**2.rcorrector**

Since rcorrector is not present as a module on the UofU CHPC, I installed it locally in my directory: /uufs/chpc.utah.edu/common/home/u6007910/projects/sam

Install rcorrector:

```bash
module load git
git clone https://github.com/mourisl/rcorrector.git #downloads rcorrector repository
cd rcorrector/
ls
make #installs and downloads jellyfish2 if it is not available in the rcorrector path
perl run_rcorrector.pl #test for installation
```

Running rcorrector:

```bash
mkdir rcorrector

#bash script for rcorrector:

#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 50:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J rcorrect
#SBATCH --mail-type=ALL  # Type of email notification- BEGIN,END,FAIL,ALL  
#SBATCH --mail-user=<samridhi.chaturvedi@gmail.com>  # Email to which notifications will be sent

module load perl

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/rcorrector/


for prefix in $(ls /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/RNAseq20151015/*.fastq | sed -r 's/_R[12]_001[.]fastq//' | uniq)

do

perl /uufs/chpc.utah.edu/common/home/u6007910/projects/sam/rcorrector/run_rcorrector.pl -t 12 -1 "${prefix}_R1_001.fastq" -2 "${prefix}_R2_001.fastq"

done
```

After rcorrector has run, we will get outfiles with extensions "fq". I used the following steps to run some statistics on these files and get numbers for how many files are tagged with "cor" or "unfixable". I saved these numbers in this file: https://docs.google.com/spreadsheets/d/1vU6Fj-MqAYnps0ksBL_gx5hS5eSfsp1y8X0TnuJgpDk/edit#gid=0

```bash
#get file names
ls *.fq | cat

#get the total number of reads in each file:
for i in *fq; do grep ^@ $i | wc -l; done

#get reads with "cor" tag
for i in *fq; do grep "cor" $i | wc -l; done

#get reads with "unfixable" tag
for i in *fq; do grep "unfixable" $i | wc -l; done
```

After this I ran the script *FilterUncorrectabledPEfastq.py* to remove the unfixable reads. This outputs files with names: unfixrm*cor.fq. This script is present in the repository.

**3. Trim galore**

After running rcorrector and removing unfixable reads, I ran Trim_galore to trim adapter sequences and erroneous k-mers from the sequences. Here is the bash script I used for running trim galore on the cluster:

```bash
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 150:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J trimgalore
#SBATCH --mail-type=ALL  # Type of email notification- BEGIN,END,FAIL,ALL  
#SBATCH --mail-user=<samridhi.chaturvedi@gmail.com>  # Email to which notifications will be sent

module load cutadapt
module load fastqc

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/trim_galore/


for prefix in $(ls /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/rcorrector/unfixrm_*_001.cor.fq | sed -r 's/_R[12]_001[.]cor.fq//' | uniq) 

do 

/uufs/chpc.utah.edu/sys/pkg/trim_galore/trim_galore --dont_gzip --paired --phred33 --length 36 -q 5 --stringency 1 -e 0.1 "${prefix}_R1_001.cor.fq" "${prefix}_R2_001.cor.fq"

done 
```

I then ran fastqc on the output files and rechecked the statistics and summary (see step 1 above). The files are now ready for alignments.

**5. Running STAR for sequence alignments**

I downloaded the data from Sahayadi et al's paper (https://www.nature.com/articles/s41559-019-1041-9) to get the gff files which are in this folder: /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/main_gff_ref.

I used this genome assembly (fasta file) and gff file for the the STAR aligment.

I then started a trial run for STAR.

Directory for STAR: /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/star

Results of first run are saved in star_results_run1
Results of second run are saved in star_results_run2 (not very different just output bam file)

Here is the bash script for running the STAR index:

```bash
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J star-index
#SBATCH --mail-type=ALL  # Type of email notification- BEGIN,END,FAIL,ALL  
#SBATCH --mail-user=<samridhi.chaturvedi@gmail.com>  # Email to which notifications will be sent

module load star

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/star/star_results_run1/

#with gff file
STAR --runMode genomeGenerate --runThreadN 12 --limitSjdbInsertNsj 300000 --genomeDir cmacgenome_starindex --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/main_gff_ref/GCA_900659725.1_ASM90065972v1_genomic.gff --genomeFastaFiles /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/main_gff_ref/GCA_900659725.1_ASM90065972v1_genomic.fasta
```

Here is the bash script to run the first run of STAR alignment:

```bash
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J star-align
#SBATCH --mail-type=ALL  # Type of email notification- BEGIN,END,FAIL,ALL  
#SBATCH --mail-user=<samridhi.chaturvedi@gmail.com>  # Email to which notifications will be sent

module load star

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/star/

for prefix in $(ls ../trimgalore/unfixrm_*_001.cor_val_*.fq | sed -r 's/_R[12]_001[.]cor_val_[12].fq//' | uniq)

do

STAR --runThreadN 12 --runMode alignReads --genomeDir ./cmacgenome_starindex --sjdbGTFfile /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/main_gff_ref/GCA_900659725.1_ASM90065972v1_genomic.gtf --readFilesIn  "${prefix}_R1_001.cor_val_1.fq" "${prefix}_R2_001.cor_val_2.fq" --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

done
```

Here is the bash script to run the first run of STAR alignment:

```
#!/bin/bash
#SBATCH -n 12 
#SBATCH -N 1
#SBATCH -t 300:00:00
#SBATCH -p usubio-kp
#SBATCH -A usubio-kp
#SBATCH -J star-align
#SBATCH --mail-type=ALL  # Type of email notification- BEGIN,END,FAIL,ALL  
#SBATCH --mail-user=<samridhi.chaturvedi@gmail.com>  # Email to which notifications will be sent

module load star

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/star/

for prefix in $(ls ../trimgalore/unfixrm_*_001.cor_val_*.fq | sed -r 's/_R[12]_001[.]cor_val_[12].fq//' | uniq)

do

STAR --runThreadN 12 --runMode alignReads --sjdbFileChrStartEnd ./star_results_run1/SJ.out.tab --genomeDir ./cmacgenome_starindex --readFilesIn  "${prefix}_R1_001.cor_val_1.fq" "${prefix}_R2_001.cor_val_2.fq" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts

done
```

The output files generated from second run are:
Aligned.sortedByCoord.out.bam  
Log.final.out  --> has the final run statistics
Log.out  
Log.progress.out  
ReadsPerGene.out.tab  
SJ.out.tab

I then checked the alignment statistics as follows:

```bash
module load samtools
samtools flagstat Aligned.sortedByCoord.out.bam 
```

The results of this are:
64079449 + 0 in total (QC-passed reads + QC-failed reads)
26409432 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
64079449 + 0 mapped (100.00% : N/A)
37670017 + 0 paired in sequencing
18835451 + 0 read1
18834566 + 0 read2
37668530 + 0 properly paired (100.00% : N/A)
37668530 + 0 with itself and mate mapped
1487 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

These look fine for downstream analyses.

**6. featureCount **

I used the output from STAR to create the final count matrix using the program featureCount. 

For running this script we need to specify the path to STAR alignments (bam files), annotated gff file for the genome. Here is the bash script I used to run featureCount on alignment outputs from STAR:

```bash
#!/bin/bash
#SBATCH --job-name=featcount
#SBATCH --time=300:00:00 #walltime
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=usubio-kp #PI account
#SBATCH --partition=usubio-kp #specify computer cluster, other option is kinspeak

cd /uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/featurecount/

DIR=/uufs/chpc.utah.edu/common/home/gompert-group1/data/callosobruchus/Annotation/transcriptome/star/alignments/

for file in ${DIR}*.sortedByCoord.out.bam
do
        sample=$(echo $file| cut -d'/' -f 13 | cut -d'.' -f 1)
        subread-2.0.0-source/bin/featureCounts -t exon -g Parent -a ../../main_gff_ref/GCA_900659725.1_ASM90065972v1_genomic.gff -o ${sample}.featureCounts $file

done
```

The output from this script are two files per sample: one with extension ".featurecounts" and one with extension ".summary".

I created the final matrix of counts by combining the data from each sample output file (".featurecounts") by using the R script *prepare_countmatrix.R*. I then annotated these transcripts with gene IDs using the python script *transcript_annot.py*. Both these scripts are in the repository. I used the annotated count matrix for differential expression analyses as described in the "differential_expression_analyses" file in the repository.

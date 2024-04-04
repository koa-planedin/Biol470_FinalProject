# Notes from Mar 27 Lab

Getting started

```console
#open a screen to save progress and keep programs running
screen -S biol470

#navigate to directory
cd /project/ctb-grego/biol470/koa_fola_lars/biol470_project

#load packages
module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9 fastqc/0.12.1
```

Looking over the quality - I ran a small subset, but the data is simulated so they're all the same.

```console
#fastQC on a subset
fastqc Biol470.p1.*.500000_R2.fastq.gz

#clean up directory
mkdir output_html
chmod ugo+rwx output_html/
mv *.html output_html/
```

Time to trim! Again, the reads are perfect, so lets not mess with that.
```console
#load the program
module load trimmomatic

#an example run
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE Biol470.p3.i2.500000_R1.fastq.gz Biol470.p3.i2.500000_R2.fastq.gz Biol470.p3.i2.500000_R1.trim.fastq Biol470.p3.i2.500000_R1.unpaired.fastq Biol470.p3.i2.500000_R2.trim.fastq Biol470.p3.i2.500000_R2.unpaired.fastq ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

```

Looking at the reference genome
```console
module load samtools/1.18
module load bwa-mem2/2.2.1

samtools faidx SalmonReference.fasta
bwa-mem2 index SalmonReference.fasta
```

Aligning to the reference
```console
#make new directory
mkdir alignment
chmod ugo+rwx alignment

#making a script for alignment
cd fastq/

#making the actual script - the contents are in another block below
nano align.sh

#make it executable
chmod +x align.sh

#run it
./align.sh
```

Script for alignments. Running this in the fastq folder and moving to alignment.
```console
for file_r1 in *R1.fastq.gz; do
    file_r2="${file_r1/R1/R2}"  # Get corresponding R2 file

    # Check if R2 file exists
    if [ -f "$file_r2" ]; then
        output_base="${file_r1%_*}"  # Remove R1.fastq extension
        bwa-mem2 mem -t 2 ../ref/SalmonReference.fasta $file_r1 $file_r2 > $output_base.sam
    else
        echo "Warning: No corresponding R2 file found for $file_r1"
    fi
mv *.sam ../alignment/

done
```

Now we need to sort these cuties
```console
cd alignment
mkdir bam
chmod 777 bam
nano sort_sam.sh
chmod +x sort_sam.sh
./sort_sam.sh
```
sort_sam.sh file:
```console
for file in *.sam; do

    # Check if R2 file exists
    if [ -f "$file" ]; then
        output_base="${file%.*}"  # Remove R1.fastq extension
       echo "starting: $output_base"
         samtools sort $file > $output_base.sort.bam
    else
        echo "Not found in directory: $file"
    fi

mv $output_base.sort.bam bam/

done
```
#marking duplicates

```console
module load picard/3.1.0

mkdir duplicates
nano mark_dups.sh
chmod +x mark_dups.sh
./mark_dups.sh
```
Script for mark_dups.sh:
```console
for file in *.bam; do

    # Check if R2 file exists
    if [ -f "$file" ]; then
        output_base="${file%.*}"  # Remove R1.fastq extension
       echo "starting: $output_base"
        java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$file O=$output_base.markdup.bam M=$o>    else
        echo "Not found in directory: $file"
    fi

mv $output_base.markdup.bam duplicates/
mv $output_base.dupmetrics.txt duplicates/

done
```
Output: bam files


Now indexing and coverage, in a file called index.sh that is running in duplicates/:
```console
for file in *.markdup.bam; do

    #check if file is there
    if [ -f "$file" ]; then
        output_base="${file%.*}"  #base name
       echo "starting: $output_base"
         samtools index $file
        samtools coverage $file > $output_base.coverage.txt
    else
        echo "Not found in directory: $file"
    fi

done
```
Output: indexed file (.bai) and coverage file (.coverage.txt)

Variant calling:
```console
#adding read groups
mkdir bam_rg
nano name_groups.sh
chmod +x name_groups.sh
./name_groups.sh

#make list of bamfiles
ls bam_rg/*.bam > bamlist.txt

######################
## bcftools version ##
######################

#running mpileup in alignment/bam/duplicates/
bcftools mpileup -q 20 -f ../../../ref/SalmonReference.fasta -b bamlist.txt > salmon_int.g.vcf

#call genotypes with bcftools
bcftools call -mv salmon_int.g.vcf | bcftools +fill-tags> salmon_int.bcftools.vcf

#filter for maf, two versions
bcftools view -q 0.1:minor salmon_int.bcftools.vcf > salmon_int.bcftools.maf10.vcf
bcftools view -q 0.01:minor salmon_int.bcftools.vcf > salmon_int.bcftools.maf1.vcf

#######################
## freebayes version ##
#######################

#calling variants
freebayes -L bamlist.txt -f ../../../ref/SalmonReference.fasta > salmon_int.fb.vcf
```

name_groups.sh file:
```console
for file in *.bam; do

    #check if file is there
    if [ -f "$file" ]; then
        output_base="${file%.*}"  #base name
       echo "starting: $output_base"

         java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
  I=$file \
  O=bam_rg/$output_base.rg.bam \
  RGSM=$output_base \
  RGLB=$output_base \
  RGID=$output_base \
  RGPL=Illumina \
  RGPU=NULL
 
  samtools index bam_rg/$output_base.rg.bam

    else
        echo "Not found in directory: $file"
    fi

done
```






Note: Greg mentioned there was an issue with renaming the chromosomes using his old code, but this would work
```console
bcftools annotate --rename_chrs
```


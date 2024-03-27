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
mv *.html output_html/
```

Time to trim!
```console
#load the program
module load trimmomatic

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE DRR053219_1.fastq DRR053219_2.fastq DRR053219_1.trim.fastq DRR053219_1.unpaired.fastq DRR053219_2.trim.fastq DRR053219_2.unpaired.fastq ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

```


Note: Greg mentioned there was an issue with renaming the chromosomes using his old code, but this would work
```console
bcftools annotate --rename_chrs
```

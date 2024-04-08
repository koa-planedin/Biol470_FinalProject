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

Admixture
(Don't trust the names to be consectutive here, but the steps are)
```console
module load plink

#renaming chrs - repeat with all vcfs, or repeat with each unfiltered vcf and make new filtered versions
bcftools annotate --rename-chrs chr_map.txt salmon_int.fb.vcf > salmon_int.fb.rename.vcf

#likewise, repeat
#I'll make a nice loop for this soon
plink --vcf vcf/salmon_int.fb.maf10.rename.vcf --out salmon_int.fb.maf10 --make-bed --allow-extra-chr --double-id --autosome-num 95

#another repeat
for K in `seq 12` /; do /; admixture --cv bed/salmon_int.fb.maf1.bed $K -j3 | tee log${K}.out; /; done

#view errors for more recent set
#when making a script I'd rename these and automatically grep the output into txt file
grep -h CV log*.out
#for saving the errors
grep -h CV log*.out > cv_error_fb.maf1.txt
```
R scripts
```console
#Get sample names in order from the .fam 
sample_names <- read_table("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/vcf/bed/salmon_int.fb.fam",col_names = F)
names(sample_names)[2]<-"sample.id"
all_data <- tibble()

#Loop for each K value
for (k in 1:12){
  #Read the Q file 
  data <- read_delim(paste0("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/vcf/salmon_int.fb.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  #Add in sample names
  data$sample.id <- sample_names$sample.id
  #Add in the K value
  data$k <- k
  
  #Convert the wide table to a long table, which is easier to plot
  data %>% gather(Q, value, -sample.id,-k) -> data
  #Bind together all the outputs of each K value
  all_data <- rbind(all_data,data)
}

#sep samples names
all_data <- all_data %>%
  separate_wider_delim(sample.id, ".", names = c(NA,"pop","indiv", NA,NA,NA)) %>%
  mutate(name = paste0(pop,".",indiv))

#remove p from population to order them numerically
all_data$pop<-sub("^p", "", all_data$pop)
all_data<-all_data %>% mutate(pop = as.numeric(pop)) %>% 
  arrange(pop)

#graph of all k values by row
#to just look at one k value, filter data and remove facet_grid
all_data %>%
  ggplot(.,aes(x=fct_reorder(name,pop),y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_grid(rows=vars(k),scales = "free", space = "free")
```


PCA and FST
```console
#PCA
plink --vcf vcf/salmon_int.fb.maf10.rename.vcf --out salmon_int.fb.maf10 --pca --allow-extra-chr --double-id --autosome-num 95


#sample names
bcftools query -l salmon_int.fb.maf10.rename.vcf > sample_names.txt
#I needed up doing some stuff with sed here
$but importantly the new name files are made:
#such as samples_pop5.rn.txt

#another repeat step
#A cute little script coming to you soon
vcftools --vcf  salmon_int.fb.maf10.rename.vcf \
--weir-fst-pop samples_pop1.rn.txt \
--weir-fst-pop samples_pop2.rn.txt \
--out salmon_int.fb.maf10
```

```console

#########
## pca ##
#########

#meta data for testing correlations
meta<-read.table("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/info/body_size.txt",header = T)
loc<-read.xlsx("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/info/Final_project_locations.xlsx")
head(meta)

#pca data
pca_data <- read_table("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/salmon_int.bcftools.maf1.eigenvec",
                       col_names = c("sample.id","spacer",paste0("PC",1:20)))  
pca_data<-pca_data %>% 
  separate_wider_delim(spacer, ".", names = c(NA,"pop","indiv", NA,NA,NA))

#remove p and/or i as needed
pca_data$pop<-sub("^p", "", pca_data$pop)
pca_data$indiv<-sub("^i", "", pca_data$indiv)

#also, make population a factor or numeric as needed
#to make numeric, just change "as.factor" to "as.numeric"
pca_data<-pca_data %>% mutate(pop = as.factor(pop),
         name = paste0(pop,".",indiv),
         indiv = as.numeric(indiv))

#adding name to meta so that they join
meta<-meta %>% mutate(name = paste0(pop,".",individual))

#join info
pca_meta <- left_join(pca_data,meta)
pca_meta$pop<-as.factor(pca_meta$pop)

#plot data
pca_data %>%
  ggplot(.,aes(x=PC1,y=PC3,col=pop)) +
  geom_point()


#########
## fst ##
#########

#load fst
fst <- read_tsv("/project/ctb-grego/biol470/koa_fola_lars/biol470_project/salmon_int.bcftools.maf1.pop1.pop2.weir.fst")

#plot it
#we only got this far, but we should probably make a grid of fst plots?
#maybe a table of SNPs too?
fst %>%
  ggplot(.,aes(x=POS,y=WEIR_AND_COCKERHAM_FST)) +
  geom_point() + facet_wrap(~CHROM)
```



# Data for the final project is available now at "/project/ctb-grego/sharing/biol470_project". It includes fastq files for the samples, a reference genome, a phenotype file (for the GWAS) and excel file telling you where each population came from (which you can use creatively).

# 0. Set up for the project and logging in

# Log in to the remote server after turning on the UVic VPN (Cisco)
ssh your_username@indri.rcs.uvic.ca

# Source the bash.sh file so that you have access to the modules
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

# Navigate to our project directory
cd /project/ctb-grego/biol470/koa_fola_lars

# Gave everyone permission to use the directory. Might need to change this and update future files and folders. U = users, G = group. 
chmod ugo+rwx koa_ola_lars/

# Copy project files from Greg's directory
cp -r  /project/ctb-grego/sharing/biol470_project /project/ctb-grego/biol470/koa_fola_lars

# Make a symbolic link on your area of the server (because it's hard to scp off Greg's)
ln -s /project/ctb-grego/biol470/koa_fola_lars

# 1. Filtering reads for adapter sequences and base quality.

# 1a. Unzip the .gz fastq files
cd fastq

gunzip *fastq.gz

# 1b. Examine the quality of the files and look for red flags using FastQC (maybe better to do a subset because this takes a long time for each file)
module load  fastqc/0.12.1

fastqc *.fastq

# Put the html files in a folder to clean up the workspace
mkdir html

mv *.html /qc_html

rm *.zip

# Open a new terminal window and run this to look at resulting .html files on your own desktop
scp username@indri.rcs.uvic.ca:/home/username/koa_fola_lars/biol470_project/fastq/html /path/to/your/desktop

# Now you can look at the quality reports

# 1c. Running Trimmomatic on the pairs of reads for each individual (Trimmomatric: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

# This is the script from our lab for unpaired reads
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE DRR053219_1.fastq DRR053219_2.fastq DRR053219_1.trim.fastq DRR053219_1.unpaired.fastq DRR053219_2.trim.fastq DRR053219_2.unpaired.fastq ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

# What i think the command will be for each paired read for each individual
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE Biol470.p1.i1.500000_R1.fastq Biol470.p1.i1.500000_R2.fastq Biol470.p1.i1.500000_R1.trim.fastq Biol470.p1.i1.500000_R1.paired.fastq Biol470.p1.i1.500000_R1.unpaired.fastq Biol470.p1.i1.500000_R2.trim.fastq Biol470.p1.i1.500000_R2.paired.fastq Biol470.p1.i1.500000_R2.unpaired.fastq ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

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

#Open a new terminal window and run this to look at resulting .html files on your own desktop
scp 
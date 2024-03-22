# Data for the final project is available now at "/project/ctb-grego/sharing/biol470_project". It includes fastq files for the samples, a reference genome, a phenotype file (for the GWAS) and excel file telling you where each population came from (which you can use creatively).

# 1. Set up for the project and logging in

# Log in to the remote server after turning on the UVic VPN (Cisco)
ssh your_username@indri.rcs.uvic.ca

# Source the bash.sh file so that you have access to the modules
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

# Navigate to our project directory
cd /project/ctb-grego/biol470/koa_ola_lars

# Gave everyone permission to use the directory. Might need to change this and update future files and folders. U = users, G = group. 
chmod ug+rwx koa_ola_lars/

# Copy project files from Greg's directory
cp -r  /project/ctb-grego/sharing/biol470_project /project/ctb-grego/biol470/koa_ola_lars

# Talk to the group about how you want to work together on code. Right now you have set up a github repository and it has been working alright to install git on your computer, as well as the GitHub desktop and Visual Code. Does Koa have a better work around? 

# How Larissa set up her coding party

# Install git
https://git-scm.com/book/en/v2/Getting-Started-Installing-Git

# Set up git basics
$ git config --global user.name "Human Name"
$ git config --global user.email humanname@example.com

# Install GitHub desktop for pushing updates to code
https://desktop.github.com/ # This might be redundant...I was just having a hard time getting VisualCode to push updates to GitHub before I downloaded GitHub desktop

# Install visual code to have a nice environment to work on your scripts within
https://code.visualstudio.com/ 

# Clone the repository for our project
https://github.com/larissaissabron/Biol470_FinalProject

#2. Assignment instructions

#For this project, you will be given a set of read data for multiple samples from multiple populatioons. Each sample will have a phenotype trait value and a population ID. The goal is to understand how the populations are related, and what that says about their history or biology, and determine the genetic basis for a phenotypic trait. Your analysis should include at least the following steps:

#1) Filtering reads for adapter sequences and base quality.
#### Lab 5: https://owensgl.github.io/intro-genomics/labs/lab_5/

#2) Alignment to a reference genome.
#### Lab 6: https://owensgl.github.io/intro-genomics/labs/lab_6/

#3) Marking or removal of duplicate reads.
#### Lab 6 (marking duplicates, could do more to remove?): https://owensgl.github.io/intro-genomics/labs/lab_6/ 

#4) Variant calling.
##### Lab 7: https://owensgl.github.io/intro-genomics/labs/lab_7/

#5) Filtering of the VCF file.
##### Lab 7: https://owensgl.github.io/intro-genomics/labs/lab_7/ 

#6) Principal Component Analysis.
#### Lab 8: https://owensgl.github.io/intro-genomics/labs/lab_8/ 

#7) Popula6on structure analysis (e.g. STRUCTURE, admixture, dapc).
#### Lab 8: https://owensgl.github.io/intro-genomics/labs/lab_8/ 

#8) Genome-wide associa6on for the trait.
##### Lab 9: https://owensgl.github.io/intro-genomics/labs/lab_9/ 


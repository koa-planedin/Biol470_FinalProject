# Data for the final project is available now at "/project/ctb-grego/sharing/biol470_project". It includes fastq files for the samples, a reference genome, a phenotype file (for the GWAS) and excel file telling you where each population came from (which you can use creatively).

# 1. Set up for the project and logging in

# Log in to the remote server after turning on the UVic VPN (Cisco)
ssh larissabron@indri.rcs.uvic.ca

# Source the bash.sh file so that you have access to the modules
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

# Navigate to our project directory
/project/ctb-grego/biol470/koa_ola_lars

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




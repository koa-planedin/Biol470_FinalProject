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

# What a good time to be alive


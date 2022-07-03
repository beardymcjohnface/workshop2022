---
layout: default
permalink: /patric
---

# Metagenome assembled genomes (MAG) Analysis

Once the MAGs are generated, here are some downstream analysis that you can perform

## Step 1: Download the MAGs generated to laptop
Goto - https://cloudstor.aarnet.edu.au/plus/s/2FJKmMoF8bOoGZd \
Select "bin_20_seqs.fasta" and download 

## Step 2: Upload files to PATRIC
PATRIC is web service available for microbial genome annotations, https://patricbrc.org

- First register for an account, if you dont already have one. Requires your name, and email. Should be quick
- Once you have an account set up, login with your credentials
- Click on "Workspaces" on the top bar
- On the right corner of the window, there is a button to create a new workspace
- Type in "metagenomics-workshop 2022", and save 
- This should open a new workspace, click on the link "Upload" on the top right corner of the workspace
- This will open a new pop up window, with the option "Select Files". Click this option 
- Navigate your laptop for the bin downloaded, and upload this file
- Change Upload Type on the pop up window from "Unspecified" to "Contigs"
- Then select "Start Upload" 
- Once all the files are uploaded, they will be listed in the workspace
  - Check the file sizes to confirm they were uploaded, shouldn't be 0 
  - Check the symbol next to the filetypes, it should look like two lines if not, then the filetype needs to be changed
    - To change filetype, click on file and click on "Edit type" to the right bar
    - Select contigs, and save
    - The symbal should be updated to the lines now

## Annotating bins/MAGs
- In PATRIC, click on "Services" on the top bar, then "Annotation"
- This should opena new tab. Next to contigs, clik on the "folder icon"
- This will open a new pop up window that allows you to navigate through your worksapce
  - find the metagenomics-workshop 2022 folder
  - select "bin_20_seqs.fasta" and click on "OK" 
- This will take you back to the Annotation window, fil the form with the below inforamtion
  - Domain: bacteria
  - Taxonomy ID 2
  - My label : Workshop bin20
  - Output folder: metagenomics-workshop 20022
  - Click on Annotate

This may take a while, so have to wait 

## Downloading the protein sequences
- The in the annotation directory
- Click on "Bacteria workshop bin20.feature_protein.fasta", and click on 
the option "Download" on the right bar


## The rest of the steps will shown as demo 

The next few steps are dependent on your research questions. So I will just demonostrate the some of the data PATRIC generates that could be useful


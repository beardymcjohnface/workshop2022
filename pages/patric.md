---
layout: default
permalink: /patric
---

# Metagenome assembled genomes (MAG) Analysis

[Link to presentation](https://flinders-my.sharepoint.com/:p:/g/personal/nala0006_flinders_edu_au/ES4Mk3jakl1GhygE67e_7GgBCj808SMLwzDpracec2xqxA?e=cRhPtc)

Once the MAGs are generated, here are some downstream analysis that you can perform

## Step 1: Download the MAGs generated to laptop
- Goto - [Click on this link to get MAGs](https://cloudstor.aarnet.edu.au/plus/s/2FJKmMoF8bOoGZd) 
- Select "bin_20_seqs.fasta" and download 

## Step 2: Upload files to PATRIC
[PATRIC](https://patricbrc.org) is web service available for microbial genome 
annotations

**Register**
- First register for an account, if you dont already have one. Requires your name, and email. Should be quick

**Setup workspace**
- Once you have an account set up, login with your credentials
- Click on "Workspaces" on the top bar
- On the right corner of the window, there is a button to create a new workspace
- Type in "metagenomics-workshop 2022", and save 

**Upload files**
- This should open a new workspace, click on the link "Upload" on the top right corner of the workspace
- This will open a new pop up window, with the option "Select Files". Click this option 
- Navigate your laptop for the bin downloaded, and upload this file
- Change Upload Type on the pop up window from "Unspecified" to "Contigs"
- Then select "Start Upload"

**Check the uploads** 
- Once all the files are uploaded, they will be listed in the workspace
  - Check the file sizes to confirm they were uploaded, shouldn't be 0 
  - Check the symbol next to the filetypes, it should look like two lines if not, then the filetype needs to be changed
    - To change filetype, click on file and click on "Edit type" to the right bar
    - Select contigs, and save
    - The symbol should be updated to the lines now

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

**Patience**
The command was submit to the cluster, waiting for the magic to happen


## Downloading the protein sequences
- In the Annotation directory
- Click on "Bacteria workshop bin20.feature_protein.fasta", and click on 
the option "Download" on the right bar

## Next steps... 
The next few steps are dependent on your research questions. So I will just demonostrate the some of the data PATRIC generates that could be useful

Here is a [Youtube video](https://www.youtube.com/watch?v=uixhLFj-L8U) recorded by Amber Skye on some 
downstream analysis that you can perform.

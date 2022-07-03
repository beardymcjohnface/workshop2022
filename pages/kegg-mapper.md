---
layout: default
permalink: /kegg-mapper
---

# KEGG mapper 

KEGG is a database resource for understanding high level functions and 
utilities in biological systems. 

[Click here to go to KEGG website](https://www.genome.jp/kegg/)

You can use this databse to identify pathway maps, biochemical reactions 
in a cell. 

## Starting file

The protein sequnces "Bacteria workshop 
bin20.feature_protein.fasta" from PATRIC. 

This was the last step in PATRIC tutorial

## Annotating the proteins to get KEGG ID
- Goto [GhostKoala](https://www.kegg.jp/ghostkoala/)
- Click on choose file, select "Bacteria workshop 
bin20.feature_protein.fasta" 
	- Select the option "genus_prokaryotes"
	- add your email address 
	- Select "request email confirmation"
- Once the form is submitted, there will be an email sent to the email 
address provided
- Check for email from "noreply@scl.kyoto-u.ac.jp via kegg.jp", and clik 
on the link to submit the job
	- There will be two links provided, first one is to submit the 
commannd and the second to delete the command 

- Wait.... Depending on how busy the cluster is on their end, this may 
take a while

- Once the command is done running, they send you an email with the 
results. So check email 
	- Email from "noreply@scl.kyoto-u.ac.jp via kegg.jp" with a link 
to the results 
	- the results are only provided for 7 days so make sure to 
download them 

 
## Results from KEGG
- The results link takes you to a new page with the annotations identified in the MAG
- Click on preview 100, next to Annotation data
	- The KO IDs are saved in this file

- Back on the main results page
	- Click on "Reconstrt Pathway" 
	- This takes you to the list of the pathways in KEGG
		- Can select for a specific pathway and have the genes found highlighted

## Comparing genes in two genomes
[NCGAS blog on KEGG analysis to compare 
two genomes](https://blogs.iu.edu/ncgas/2019/02/21/visualizing-kegg-pathways/)

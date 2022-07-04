---
layout: default
permalink: /qc-qa
---

# QC/QA Hands on tutorial

[Link to presentation](https://flinders-my.sharepoint.com/:p:/g/personal/nala0006_flinders_edu_au/Edk3RVzYyiVIgn1VRONqsXUBc_OC-J_ZAQjAcI5nRwPrCg?e=557EOV)

## Data
Login
- Using putty/mobaxterm 
- Command line using 
	
	ssh username@115.146.84.253

List all the files, using command 
	
	ls
  
There should hopefully be a folder labelled "subsample"

If you have only "subsample.tar" folder, this means you have to decompress this file. To do this run the command 
  
  `tar -xvf subsample.tar`

### [If you need to download the files, follow the instructions here](/workshop2022/files/fastq/README.md)

List out the files in the subsample folder, 

  `ls subsample`
  
This should list the 8 metagenomes in the folder 

    SRR1237780_1.fastq.gz  SRR1237781_1.fastq.gz  SRR1237782_1.fastq.gz  SRR1237783_1.fastq.gz
    SRR1237780_2.fastq.gz  SRR1237781_2.fastq.gz  SRR1237782_2.fastq.gz  SRR1237783_2.fastq.gz

Decompress the files by running the command 

    gzip -d SRR1237780_1.fastq.gz
    gzip -d SRR1237780_2.fastq.gz
    gzip -d SRR1237781_1.fastq.gz
    gzip -d SRR1237781_2.fastq.gz
    gzip -d SRR1237782_1.fastq.gz
    gzip -d SRR1237782_2.fastq.gz
    gzip -d SRR1237783_1.fastq.gz
    gzip -d SRR1237783_2.fastq.gz

If the files are decompressed they should no longer have the .gz extension

## FastQC 
A quality control tool for high throughput sequence data - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Instead of downloading the application locally, we can run the application as commands. 
Run the below commands, 

    cd
    mkdir fastqc
    fastqc subsample/SRR1237780_1.fastq -o fastqc/
 
The FastQC report for this file is saved to the directory fastqc. Check this diretcory after the command is run

    cd 
    cd fastqc/
    ls

This should list the below files

    (base) nala0006@workshoppy:~/fastqc$ ls
    SRR1237780_1_fastqc.html  SRR1237780_1_fastqc.zip

To visualize the output, copy the file to your laptop. 
- This can be done using WinsSCP (https://winscp.net/eng/index.php)
- Command line 

      scp -r nala0006@115.146.84.253:/home/nala0006/fastqc/SRR1237780_1_fastqc.html .

Once copied to laptop, click on the file and it should open in a browser

**Repeat the same steps for reverse read, SRR1237780_2**
Commands to run 

    fastqc subsample/SRR1237780_2.fastq -o fastqc/
    
Now there should be 4 files in fastqc

    (base) nala0006@workshoppy:~$ ls fastqc/
    SRR1237780_1_fastqc.html  SRR1237780_1_fastqc.zip  SRR1237780_2_fastqc.html  SRR1237780_2_fastqc.zip
   
Download SRR1237780_2_fastqc.html, and open in browser

## Prinseq 
Prinseq++ can be used to trim, filter and reformat sequences.

    cd 
    mkdir prinseq 
    
    prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 0 -derep 1  -out_format 0 -trim_tail_left 5 -trim_tail_right 5  -trim_qual_type min -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10-threads 2 -out_name prinseq/SRR1237780 -fastq subsample/SRR1237780_1.fastq -fastq2 subsample/SRR1237780_2.fastq 

    prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 0 -derep 1  -out_format 0 -trim_tail_left 5 -trim_tail_right 5  -trim_qual_type min -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10-threads 2 -out_name prinseq/SRR1237781 -fastq subsample/SRR1237781_1.fastq -fastq2 subsample/SRR1237781_2.fastq 

    prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 0 -derep 1  -out_format 0 -trim_tail_left 5 -trim_tail_right 5  -trim_qual_type min -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10-threads 2 -out_name prinseq/SRR1237782 -fastq subsample/SRR1237782_1.fastq -fastq2 subsample/SRR1237782_2.fastq 

    prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 0 -derep 1  -out_format 0 -trim_tail_left 5 -trim_tail_right 5  -trim_qual_type min -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10-threads 2 -out_name prinseq/SRR1237783 -fastq subsample/SRR1237783_1.fastq -fastq2 subsample/SRR1237783_2.fastq 

Output from prinseq would contain, 
  - good_out: filtered sequences for forward and reverse reads
  - bad_out: removed sequences for forward and reverse reads
  - single_out: singletons for forward and reverse reads

        (base) nala0006@workshoppy:~$ ls prinseq/
        SRR1237780_bad_out_R1.fastq     SRR1237781_good_out_R1.fastq    SRR1237782_single_out_R1.fastq
        SRR1237780_bad_out_R2.fastq     SRR1237781_good_out_R2.fastq    SRR1237782_single_out_R2.fastq
        SRR1237780_good_out_R1.fastq    SRR1237781_single_out_R1.fastq  SRR1237783_bad_out_R1.fastq
        SRR1237780_good_out_R2.fastq    SRR1237781_single_out_R2.fastq  SRR1237783_bad_out_R2.fastq
        SRR1237780_single_out_R1.fastq  SRR1237782_bad_out_R1.fastq     SRR1237783_good_out_R1.fastq
        SRR1237780_single_out_R2.fastq  SRR1237782_bad_out_R2.fastq     SRR1237783_good_out_R2.fastq
        SRR1237781_bad_out_R1.fastq     SRR1237782_good_out_R1.fastq    SRR1237783_single_out_R1.fastq
        SRR1237781_bad_out_R2.fastq     SRR1237782_good_out_R2.fastq    SRR1237783_single_out_R2.fastq
        
 Yay!! All done with this tuorial too

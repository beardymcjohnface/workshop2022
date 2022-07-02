# QC/QA Hands on tutorial

## Data 
Login
- Using putty/mobaxterm 
- Command line using `ssh username@115.146.84.253`

List all the files, using command `ls`
  
  A folder labelled "subsample"

If you have only "subsample.tar" folder, this means you have to decompress this file. To do this run the command 
  
  `tar -xvf subsample.tar`

List out the files in the subsample folder, 

  `ls subsample`
  
This should list the 8 metagenomes in the folder 

`SRR1237780_1.fastq.gz  SRR1237781_1.fastq.gz  SRR1237782_1.fastq.gz  SRR1237783_1.fastq.gz
 SRR1237780_2.fastq.gz  SRR1237781_2.fastq.gz  SRR1237782_2.fastq.gz  SRR1237783_2.fastq.gz`

Decompress the files by running the command 

`gzip -d SRR1237780_1.fastq.gz
 gzip -d SRR1237780_2.fastq.gz
 gzip -d SRR1237781_1.fastq.gz
 gzip -d SRR1237781_2.fastq.gz
 gzip -d SRR1237782_1.fastq.gz
 gzip -d SRR1237782_2.fastq.gz
 gzip -d SRR1237783_1.fastq.gz
 gzip -d SRR1237783_2.fastq.gz`

If the files are decompressed they should no longer have the .gz extension

## FastQC 
A quality control tool for high throughput sequence data - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Instead of downloading the application locally, we can run the application as commands. 
Run the below commands, 

`cd
 mkdir fastqc
 fastqc subsample/SRR1237780_1.fastq -o fastqc/`
 
The FastQC report for this file is saved to the directory fastqc. Check this diretcory after the command is run

`cd 
 cd fastqc/
 ls`

This should list the below files

`(base) nala0006@workshoppy:~/fastqc$ ls
SRR1237780_1_fastqc.html  SRR1237780_1_fastqc.zip`

To visualize the output, copy the file to your laptop. 
- This can be done using WinsSCP (https://winscp.net/eng/index.php)
- Command line 
`scp -r nala0006@115.146.84.253:/home/nala0006/fastqc/SRR1237780_1_fastqc.html .`

Once copied to laptop, click on the file and it should open in a browser

## Prinseq 



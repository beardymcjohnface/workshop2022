---
layout: default
permalink: /metagenomic-binning
---

# Metagenomic Binning

Metagenomic binning is the process of separating sequences into bins representing different taxonomic groups to form draft genomes, generally known as **metagenome assembled genomes**. The most popular tools bin contigs as they are longer than reads and contain more information that can characterise the underlying organisms in the metagenomic sample.

In this session, we will use [MetaBAT](https://bitbucket.org/berkeleylab/metabat/src/master/), a popular programs for binning metagenomes to bin our contigs.

# Before binning...

We have already assembled the reads and provided the necessary assembly files to do binning. Let's download these files and transfer them to the server.

## Download the files to your laptop

Go to [https://cloudstor.aarnet.edu.au/plus/s/2oFEKQv68UZoFS9](https://cloudstor.aarnet.edu.au/plus/s/2oFEKQv68UZoFS9).

Select all the files and click on `Download`.

The files will download to your laptop as `assembly_binning.tar`.

## Transfer the files to the server provided for workshop

On your laptop, use [WinSCP](https://winscp.net/eng/index.php) to transfer files.

If you are using the command line on your laptop, enter the following command and your password when prompted.

```
scp assembly_binning.tar username@115.146.84.253:~/.
```

You should see the progress of your upload.

## Confirm files were transferred to the server 

Log in to the server using `ssh` command and enter your password when prompted.

```
ssh username@115.146.84.253
```
    
OR use putty or MobaXterm to login instead.

Once logged in, type the command `ls`. You should see `assembly_binning.tar` in the list.

## Decompress `assembly_binning.tar`

Type in the following command to extract the files from `assembly_binning.tar`.

```
tar -xvf assembly_binning.tar && rm assembly_binning.tar
ls
```
  
If you type the command `ls`, you will see the folder `assembly_binning`. 

Go into the folder `assembly_binning` and list the content using the following command.

```
cd assembly_binning
ls
```

You will see the following main files.

* `assembly_graph_with_scaffolds.gfa`
* `contigs.fasta`
* `contigs.paths`
* `depth.txt`
* `metabat2_contig_bins.csv`
* `prep_result.py`

Now we are set to bin our data using MetaBAT.


# Binning using MetaBAT

MetaBAT does the profiling based on tetranucleotide frequencies in the samples. It counts the 4-mers, (i.e. `AAAA`, `AAAT`, `AAAG`, `AAAC`, `AATA`, `AATT`, …) in the sequences and uses those to suggest which samples should go together. The advantage of MetaBAT is that it will return your contigs for you.

MetaBAT also uses coverage information of the contigs which is calculated from BAM files of every input sample. We need to map the individual reads to the contigs and get separate BAM files.

Since it takes a lot of time to generate the BAM files and create the depth file, we have already provided you with the depth file (`depth.txt`).

## Generate `depth.txt` file for your own data

If you need to create BAM files for your own data, you can use the following command to create BAM files for each  sample name listed in a file named `file.txt`.

```
mkdir bam_files
for FASTA in `cat file.txt`; do minimap2 -ax sr contigs.fasta $f_good_out_R1.fastq $f_good_out_R2.fastq | samtools sort > bam_files/$f.bam done
```

Now we can generate a depth file from the BAM files using the command `jgi_summarize_bam_contig_depths` provided from MetaBAT. This goes through the BAM files and calculates the average number of reads covering each base of the contigs. You can run the following command to generate the depth file.

```
jgi_summarize_bam_contig_depths --outputDepth bam_files/depth.txt *.bam
```

If you enter `head depth.txt` you will see the following output.

```
contigName	contigLen	totalAvgDepth	SRR1237780.bam	SRR1237780.bam-var	SRR1237781.bam	SRR1237781.bam-var	SRR1237782.bam	SRR1237782.bam-var	SRR1237783.bam	SRR1237783.bam-var
NODE_1_length_110371_cov_8.539967	110371	20.851	0.535061	3.59679	0.0114134	0.0112765	18.813	41.6837	1.49151	2.26044
NODE_2_length_60492_cov_10.103049	60492	24.3014	0.0307083	0.0435139	0.00884956	0.00876961	22.3116	22.6173	1.95033	2.2689
NODE_3_length_57035_cov_7.080186	57035	21.9248	2.48215	2.94401	0.804096	0.81977	15.9555	16.8353	2.68306	3.23609
NODE_4_length_54641_cov_7.481717	54641	23.0121	2.11929	2.316	0.465783	0.547674	16.9546	28.0057	3.47245	6.37836
NODE_5_length_54014_cov_16.257900	54014	44.0697	3.03223	3.37595	2.11247	2.12222	35.9221	44.55	3.00295	3.1778
NODE_6_length_48409_cov_10.290069	48409	24.9984	0.0504362	0.0953332	0.0121014	0.015298	22.8362	22.5818	2.09973	2.32797
NODE_7_length_47877_cov_7.628414	47877	23.3906	2.46823	2.47848	0.97712	1.0791	17.2922	22.9389	2.65305	3.35134
NODE_8_length_46824_cov_7.158524	46824	22.1051	2.49087	2.48935	0.862172	0.872394	16.2341	18.6561	2.51791	3.07479
NODE_9_length_46478_cov_5.442970	46478	20.1159	0.621849	0.710255	0.233876	0.256314	12.0661	13.6028	7.19403	7.6685
```

The rows represent contigs and columns represent further information about the contigs and coverage in each sample.

## Run MetaBAT2

Let's bin our contigs using MetaBAT. We will use **MetaBAT2** which is the latest version of MetaBAT.

```
metabat2 -i contigs.fasta -a depth.txt -o metabat_bins/bin
```

You can find the bins identified by MetaBAT2 in the folder `metabat_bins`.

# Refining bins using GraphBin

[GraphBin](https://github.com/metagentools/GraphBin) is a metagenomic contig bin refinement tool that makes use of the contig connectivity information from the assembly graph to refine bins and recover contigs discarded. 

Once we have our bins from MetaBAT2, we can use GraphBin to refine the MetaBAT2 binning result.

## Preparing the binning result

We have to format our binning result so that GraphBin can read the binning result. We have provided a script named [prep_result.py](https://github.com/beardymcjohnface/workshop2022/blob/gh-pages/scripts/prep_result.py) which you can find in the bundled data. You can run it as follows.

```
python prep_result.py --binned metabat_bins --output ./
```

## Run GraphBin

Since, we used metaSPAdes to assembly out data, we can run the metaSPAdes version of GraphBin as follows.

```
graphbin --assembler spades --graph assembly_graph_with_scaffolds.gfa --contigs contigs.fasta --paths contigs.paths --binned initial_contig_bins.csv --output ./
```

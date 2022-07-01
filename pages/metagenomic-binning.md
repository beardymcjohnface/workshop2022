---
layout: default
---

# Metagenomic Binning

Metagenomic binning is the process of separating sequences into bins representing different taxonomic groups. These bins are generally known as *metagenome assembled genomes*. The most popular tools bin contigs as they are longer than reads and contain more information that can characterise the underlying organisms in the metagenomic sample.

## Before binning...

We have already assembled the reads and provided the necessary assembly files at https://cloudstor.aarnet.edu.au/plus/s/2oFEKQv68UZoFS9. 

## Binning using MetaBAT

One of the most popular programs for binning metagenomes is [MetaBAT](https://bitbucket.org/berkeleylab/metabat/src/master/).

MetaBAT does the profiling based on tetranucleotide frequencies in the samples. It counts the 4-mers, (i.e. `AAAA`, `AAAT`, `AAAG`, `AAAC`, `AATA`, `AATT`, …) in the sequences and uses those to suggest which samples should go together. The advantage of MetaBAT is that it will return your contigs for you.

MetaBAT also uses coverage information of the contigs which is calculated from BAM files for every input sample. We need to map the individual reads to the contigs in separate files.

You can use the following command to create BAM files for each input sample:

```
mkdir bam_files
for FASTA in *_good_out_R1.fastq; do minimap2 -ax sr contigs.fasta $f_good_out_R1.fastq $f_good_out_R2.fastq | samtools sort > bam_files/$f.bam done
```

Now we can generate a depth file from the BAM files using the following command:

```
jgi_summarize_bam_contig_depths --outputDepth bam_files/depth.txt *.bam
```

Since it takes a lot of time to generate the BAM files and create the depth file, we have already provided you with the depth file (`depth.txt`) which can be found in https://cloudstor.aarnet.edu.au/plus/s/2oFEKQv68UZoFS9.

Now we can run MetaBAT. We will use MetaBAT2 which is the latest version of of MetaBAT.

```
metabat2 -i contigs.fasta -a depth.txt -o metabat_bins/bin
```

You can find the identified bins in the folder `metabat_bins`.

## Refining bins using GraphBin

[GraphBin](https://github.com/metagentools/GraphBin) is a metagenomic contig bin refinement tool that makes use of the contig connectivity information from the assembly graph to refine bins and recover contigs discarded. 

Once we have our bins from MetaBAT2, we can use GraphBin to refine the MetaBAT2 binning result.

### Preparing the binning result

We have to format our binning result so that GraphBin can read the binning result. We have provided a script named [prep_result.py](https://github.com/beardymcjohnface/workshop2022/blob/gh-pages/scripts/prep_result.py) which you can download.

```
python prep_result.py --binned metabat_bins --output ./
```

We can run the metaSPAdes version of GraphBin as follows.

```
graphbin --assembler spades --graph assembly_graph_with_scaffolds.gfa --contigs contigs.fasta --paths contigs.paths --binned initial_contig_bins.csv --output ./
```

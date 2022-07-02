# Kraken 
Kraken is a tool which provides DNA sequences with taxonomic labels by examining the k-mers within a query sequence. In this tutorial we will use Kraken to determine the lowest common ancestor (LCA) for each sequence.

## Getting started 
To run kraken, we require a database which allow us to determine the taxonomy of our sequences based on k-mers. A kraken2 database has already been downloaded for us and is located at `/data/kraken` 

## Running Kraken 
Kraken has already been installed for us. We can run kraken on one of our samples (SRR1237782) by running the command

```
kraken2 --threads 4 --memory-mapping --db /data/kraken2/ --report SRR1237782_kraken_report.txt --use-names --report-zero-counts  --output SRR1237782_kraken_output.txt prinseq/SRR1237782_good_out_R2.fastq 
```

In this command: 
- `SRR1237782_kraken_report.txt` is the name of the output file which contains the abundance of each taxon. 
- `SRR1237782_kraken_output.txt` is the name of the output file which lists what taxon each read was mapped to. 

The `--memory-mapping` flag is important to ensure that the entire database isn't loaded into RAM and the `--use-names` flag is also important to make sure that our output contains the scientific name of each species rather than their taxonomy IDs. We don't need to include the `-report-zero-counts` flag, this just adds rows to our output for species which aren't present. Sometimes scripts used to analyse kraken output require this rows to be included. 

## Visualising with Krona 

Krona can be used to generate pretty figures which tell us about the composition of our sample. 

#TODO get rob to install taxonomy for krona 

We can create a Krona plot for the Kraken output of one of our samples by running the command 

```
ktImportTaxonomy -q 2 -t 3  -o  SRR1237782_krona.html  SRR1237782_kraken_report.txt
```

Next download the Krona html file to your desktop using WinSCP. You can open this html  browser to reveal a plot of the distribution of functions in the samples. You can zoom in and zoom out to see the different levels of annotations. 

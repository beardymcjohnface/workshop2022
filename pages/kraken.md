# Kraken 
Kraken is a tool which provides DNA sequences with taxonomic labels by examining the k-mers within a query sequence. In this tutorial we will use Kraken to determine the lowest common ancestor (LCA) for each sequence.

## Getting started 
To run kraken, we require a database which allow us to determine the taxonomy of our sequences based on our k-mers. A kraken2 database has already been downloaded for us and is located at `/data/kraken` 

## Running Kraken 
Kraken has already been installed for us. We can run kraken on our sequences by running the command 

```
kraken2 --threads 4 --memory-mapping --db /data/kraken2/ --report kraken_report.txt --use-names --quick --output kraken_output.txt good_out_R1/*  
```

In this command: 
- `kraken_report.txt` is the name of the output file which contains the number of reads mapped to each taxon. 
- `kraken_output.txt` refers to the #TODO 

the `--memory-mapping` flag is important to ensure that the entire database isn't loaded into RAM and the `--use-names` flag is also important to make sure that our output contains the scientific name of each species rather than their taxonomy IDs. 

## Visualising with Krona 

Krona can be used to generate pretty figures which tell us about the composition of our sample. Krona has already been installed for us as well.

```
ktImportTaxonomy -q 2 -t 3 Sample1.txt Sample2.txt -o krona.html
```

#TODO do steps to figure out how we can run the taxonomy 
Next download the Krona html file to your desktop using WinSCP. You can open this html  browser to reveal a plot of the distribution of functions in the samples. You can zoom in and zoom out to see the different levels of annotations. 

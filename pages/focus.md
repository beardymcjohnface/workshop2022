---
layout: default
permalink: /focus
---

# FOCUS 
FOCUS (Find Organisms by Composition Usage) is a computational tool which allows us to determine the taxa present in a metagenomic sequencing using an alignment free approach. 

In this tutorial you will learn how to determine which taxa are present in a metagenome. 

## Before we begin ... 

FOCUS requires that your input fastq files be in a directory. We will create a directory to store your fastq files. 
As we will be running SUPER-FOCUS only on the R1 reads we will copy just the R1 fastq files into our directory. This is commonly done as the R2 reads should have a similar (virtually identical) profile. 

```
mkdir good_out_R1
cp prinseq/*good_out_R1.fastq good_out_R1
```

## Running FOCUS 

Lucky for us, FOCUS has already been downloaded and configured with a database. This means are ready to run FOCUS!\
(If you need to install FOCUS in the future you should refer to the [FOCUS github](https://github.com/metageni/FOCUS) for instructions)

Run the command 
```
focus -q good_out_R1 -o focus_out
```

When we run this command, a new directory will be created named `focus_out` which will contain files generated by FOCUS. 

## Great, now what? 

We can start to look at the output which FOCUS by taking looking in the output directory
```
cd focus_out
```
You'll notice that there is separate table for each taxonomic level. 

To look at the output for the genus level, run the command

```
column -t -s $'\t' -n string output_Genus_tabular.xls | less
```

You will notice that this file contains four columns, one for each sample. The values in the table correspond to the relative abundance of each taxa. 

When you are done looking at the output, press the letter 'q' on your keyboard. 

## Visualising FOCUS with Krona 
     
We can build a Krona plot on our SUPER-FOCUS output just like we did for our Kraken output earlier today.

The *slightly* annoying part of building a Krona plot on our FOCUS output is that we need to rearrange it into a format which Krona can understand. The good news is that we can rearrange it using this bash command. 

```
tail -n +5 output_All_levels.xls | awk -F '\t' '{n=$9+$10; print n"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > focus_out_krona.tsv
``` 

This creates a file, `focus_out_krona.tsv` which can be read by Krona.

Next we can generate our krona plot by running the command

```
ktImportText focus_out_krona.tsv -o focusKronaPlot.html
```

Next download the Krona html file to your desktop using WinSCP. 
```
scp -r grig0076@115.146.84.253:/home/grig0076/focus_out/focus_out_krona.tsv .
```

scp -r grig0076@115.146.84.253:/home/grig0076/focus_out/focus_out_krona.tsv .
You can open this html file in your favourite browser to reveal a plot of the taxonomic composition of the samples. You can zoom in and zoom out to see the different levels of annotations. 

Congratulations on making it to the end of the tutorial! Stay tuned for the next tutorial. 

__[Feeling lazy? Here's the final product!](/workshop2022/files/focusKronaPlot.html)__

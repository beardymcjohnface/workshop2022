---
layout: default
---

# Data visualisation Part 2

In this section, we will work through some different plots.
Files to download:
- [kraken.taxon.tsv](/files/kraken/kraken.taxon.tsv)

# Preparation

You can do these steps yourself if you like, 
but we'll be skipping to the next step in the interest of time.

<details style="background-color:lightyellow;">
<summary><b>Click to expand</b></summary>

Files I used:
- [SRR1237780_kraken_output.txt](/files/kraken/SRR1237780_kraken_output.txt)
- [SRR1237781_kraken_output.txt](/files/kraken/SRR1237781_kraken_output.txt)
- [SRR1237782_kraken_output.txt](/files/kraken/SRR1237782_kraken_output.txt)
- [SRR1237783_kraken_output.txt](/files/kraken/SRR1237783_kraken_output.txt)

Concatonate and convert using taxonkit.

```shell
cat *output.txt > kraken.all.tsv
head kraken.all.tsv
```

```text
C       SRR1237780.2151832.1    Comamonas piscis (taxid 1562974)        101     1562974:1 0:54 80840:1 0:11
C       SRR1237780.4926005.1    Bifidobacterium adolescentis (taxid 1680)       101     1680:67
C       SRR1237780.501486.1     Bacteria (taxid 2)      101     2:2 1:65
U       SRR1237780.1331889.1    unclassified (taxid 0)  89      0:55
U       SRR1237780.3407698.1    unclassified (taxid 0)  101     0:67
U       SRR1237780.973930.1     unclassified (taxid 0)  101     0:67
U       SRR1237780.3542719.1    unclassified (taxid 0)  101     0:67
U       SRR1237780.4639012.1    unclassified (taxid 0)  101     0:67
C       SRR1237780.1088785.1    Arcobacter cryaerophilus ATCC 43158 (taxid 1032070)     101     0:15 1032070:3 0:37 1032070:2 0:10
U       SRR1237780.1066259.1    unclassified (taxid 0)  101     0:67
```

This is one line per read, we just want the counts for each taxid per sample.

```shell
cut -f2,3 kraken.all.tsv \
  | sed 's/\.\S\+//' \
  | sed 's/\t.\+taxid \([0-9]\+\))/\t\1/' \
  | awk '{if ($2==0){$2=1};print $1"\t"$2}' \
  | sort \
  | uniq -c \
  | sed 's/^\s\+//' \
  | sed 's/ /\t/' \
  > kraken.clean.tsv

head kraken.clean.tsv
```

This leaves us with a 3-column file: count <tab> sampleID <tab> taxid.

```text
14382   SRR1237780      0
353     SRR1237780      1
1       SRR1237780      1005058
1       SRR1237780      1028416
1       SRR1237780      1028801
1       SRR1237780      1028989
1       SRR1237780      1029756
14      SRR1237780      1032070
2       SRR1237780      1032071
2       SRR1237780      1032239
```

We convert taxid to a taxon path with taxonkit.

```shell
taxonkit lineage -i 3 kraken.clean.tsv \
  | taxonkit reformat -i 4 -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -F --fill-miss-rank \
  | cut --complement -f3,4 \
  > kraken.taxon.tsv

head kraken.taxon.tsv
```

```text
14735   SRR1237780      unclassified root superkingdom  unclassified root phylum        unclassified root class unclassified root order unclassified root family        unclassified root genus unclassified root species
1       SRR1237780      Bacteria        Proteobacteria  Gammaproteobacteria     Pasteurellales  Pasteurellaceae Gallibacterium  Gallibacterium anatis
1       SRR1237780      Bacteria        Proteobacteria  Gammaproteobacteria     Pseudomonadales Moraxellaceae   Psychrobacter   Psychrobacter sp. DAB_AL43B
1       SRR1237780      Bacteria        Proteobacteria  Alphaproteobacteria     Hyphomicrobiales        Rhizobiaceae    Neorhizobium    Neorhizobium galegae
1       SRR1237780      Bacteria        Proteobacteria  Gammaproteobacteria     Pseudomonadales Pseudomonadaceae        Pseudomonas     Pseudomonas sp. StFLB209
1       SRR1237780      Bacteria        Proteobacteria  Alphaproteobacteria     Hyphomicrobiales        Hyphomicrobiaceae       Hyphomicrobium  Hyphomicrobium nitrativorans
14      SRR1237780      Bacteria        Proteobacteria  Epsilonproteobacteria   Campylobacterales       Campylobacteraceae      Aliarcobacter   Aliarcobacter cryaerophilus
2       SRR1237780      Bacteria        Proteobacteria  Epsilonproteobacteria   Campylobacterales       Campylobacteraceae      Aliarcobacter   Aliarcobacter cryaerophilus
2       SRR1237780      Bacteria        Proteobacteria  Epsilonproteobacteria   Campylobacterales       Campylobacteraceae      Aliarcobacter   Aliarcobacter skirrowii
1       SRR1237780      Bacteria        Proteobacteria  Epsilonproteobacteria   Campylobacterales       Campylobacteraceae      Aliarcobacter   Aliarcobacter thereius
```

Let's add a header.

```shell
sed  -i '1i Count\tSampleID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' kraken.taxon.tsv

head kraken.taxon.tsv
```

```text
Count   SampleID        Kingdom Phylum  Class   Order   Family  Genus   Species
14735   SRR1237780      unclassified root superkingdom  unclassified root phylum        unclassified root class unclassified root order unclassified root family        unclassified root genus unclassified root species
1       SRR1237780      Bacteria        Proteobacteria  Gammaproteobacteria     Pasteurellales  Pasteurellaceae Gallibacterium  Gallibacterium anatis
1       SRR1237780      Bacteria        Proteobacteria  Gammaproteobacteria     Pseudomonadales Moraxellaceae   Psychrobacter   Psychrobacter sp. DAB_AL43B
1       SRR1237780      Bacteria        Proteobacteria  Alphaproteobacteria     Hyphomicrobiales        Rhizobiaceae    Neorhizobium    Neorhizobium galegae
1       SRR1237780      Bacteria        Proteobacteria  Gammaproteobacteria     Pseudomonadales Pseudomonadaceae        Pseudomonas     Pseudomonas sp. StFLB209
1       SRR1237780      Bacteria        Proteobacteria  Alphaproteobacteria     Hyphomicrobiales        Hyphomicrobiaceae       Hyphomicrobium  Hyphomicrobium nitrativorans
14      SRR1237780      Bacteria        Proteobacteria  Epsilonproteobacteria   Campylobacterales       Campylobacteraceae      Aliarcobacter   Aliarcobacter cryaerophilus
2       SRR1237780      Bacteria        Proteobacteria  Epsilonproteobacteria   Campylobacterales       Campylobacteraceae      Aliarcobacter   Aliarcobacter cryaerophilus
2       SRR1237780      Bacteria        Proteobacteria  Epsilonproteobacteria   Campylobacterales       Campylobacteraceae      Aliarcobacter   Aliarcobacter skirrowii
```

</details>

# Relative abundance comparisons (Rstudio)

Download the kraken.taxon.tsv file above and start up Rstudio.
Set your working directory to the files location.
Load the libraries and read in the file.

```r
library(tidyr)
library(dplyr)
library(ggplot2)

data = read.csv('kraken.taxon.tsv',header=T,sep='\t')
View(data)
```

We can plot this directly if we want to look at the species compositions,
but, there are a ton of species and any plots would have way too many points.
Instead, we'll look at a higher taxon level and filter for the most abundant taxa.

```r
phylumCounts = data %>% 
  group_by(SampleID, Phylum) %>% 
  summarise(n = sum(Count))
View(phylumCounts)
```

The most commonly seen visualisation for this would be a stacked bar chart of % abundances:

```r
ggplot(phylumCounts, aes(x=SampleID, y=n, fill=Phylum)) +
  geom_bar(stat='identity',position='fill') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

There are a bunch of very low abundance phyla that we dont care about.
Let's get a list of the most abundant phyla and lump the rest into an 'other' category.

```r
topPhyla = data %>% group_by(Phylum) %>% 
  summarise(n=sum(Count)) %>%
  filter(n>50) %>% 
  pull(Phylum)

phylumOtherCounts = phylumCounts
phylumOtherCounts[!(phylumOtherCounts$Phylum %in% topPhyla),]$Phylum = 'Other'
phylumOtherCounts = phylumOtherCounts %>% group_by(SampleID, Phylum) %>% summarise(n = sum(n))
View(phylumOtherCounts)
```

Let's make sure the Phyla are plotted in order of total abundance (across all samples).
We have to recalculate the summaries as we now have a new 'other' category.
Then, refactor _Phylum_ in _phylumOtherCounts_:

```r
PhylaRank = phylumOtherCounts %>% 
  group_by(Phylum) %>% 
  summarise(n=sum(n)) %>% 
  arrange(n) %>% 
  pull(Phylum)

phylumOtherCounts$Phylum = factor(phylumOtherCounts$Phylum, levels = PhylaRank)
```

Remake the plot:

```r
ggplot(phylumOtherCounts, aes(x=SampleID, y=n, fill=Phylum)) +
  geom_bar(stat='identity',position='fill') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](/files/kraken/stackedBarPlot.png)

A much better version of this is a bubble plot.
Bubble plots are not only clearer, but allow an extra dimension of information as well.
Let's collect the top 20ish families, and plot them coloured by their phyla, for each sample.

```r
# grab the top 20 families
top20Fams = data %>% 
  group_by(Family) %>% 
  summarise(n=sum(Count)) %>% 
  arrange(desc(n)) %>% 
  pull(Family)
top20Fams = top20Fams[1:20]

# get the family counts per sample
FamSums = data %>% 
  group_by(SampleID,Phylum,Family) %>% 
  summarise(n=sum(Count))

# rename non top20 families
FamSums[!(FamSums$Family %in% top20Fams),]$Family = 'Other'
FamSums[!(FamSums$Family %in% top20Fams),]$Phylum = 'NA'

# work out the plotting order for families
FamRank = FamSums %>%
  group_by(Phylum,Family) %>%
  summarise(n=sum(n)) %>%
  arrange(Phylum, n) %>%
  pull(Family)

FamSums$Family = factor(FamSums$Family, levels=(FamRank))
```

To plot, we simply use geom_point().
Plot the SampleID against the Family instead of the abundance values.
The abundance values are mapped to the points' sizes.
The extra dimension of information now becomes the coloring of the points.
I've sorted and colored the points by Phylum, but you could choose some other rank,
or get creating and split the data another way!

```r
ggplot(FamSums, aes(x=SampleID,y=Family,fill=Phylum,size=n)) + 
  geom_point(shape=21,color='black') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_size(range=c(0,10))
```

![](/files/kraken/bubblePlot.png)

# todo: PCA

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

[File preparation steps](/pages/data-viz-2-prep.html)

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
topPhyla = data %>% 
  group_by(Phylum) %>% 
  summarise(n=sum(Count)) %>%
  filter(n>50) %>% 
  pull(Phylum)

phylumOtherCounts = phylumCounts
phylumOtherCounts[!(phylumOtherCounts$Phylum %in% topPhyla),]$Phylum = 'Other'
phylumOtherCounts = phylumOtherCounts %>% 
  group_by(SampleID, Phylum) %>% 
  summarise(n = sum(n))
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

![](/workshop2022/files/kraken/stackedBarPlot.png)

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

![](/workshop2022/files/kraken/bubblePlot.png)

# todo: PCA

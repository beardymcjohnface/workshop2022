---
layout: default
---

# Preparation of Kraken output for Data visualation part 2

You don't have to do these steps, they're just a handy reference for messing with Kraken output in bash.

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


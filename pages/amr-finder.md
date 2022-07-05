---
layout: default
permalink: /amrfinder
---

# Finding AMRs with AMRFinder

The NCBI has [a very nice tool](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) for identifying Anti-Micriobial Resistance (AMR) genes,
as well as virulence factors and toxins.

## Installation

Installing AMRFinder is easily done with conda.

```shell
# create a new environemnt
conda create -n amrfinder

# activate the environment
conda activate amrfinder

# install AMRFinder (specify bioconda and conda-forge channels)
conda install -c bioconda -c conda-forge ncbi-amrfinderplus
```

## Running

Running AMRFinder is relatively easy.
Before the first time you run AMRFinder, you need to download its databases.

```shell
amrfinder --update
```

This will download the database files for you.
Now you can start looking for AMRs!
If you have an assembly, or some sequencing reads, you would use the nucleotide option.

```shell
# run on an assembly
amrfinder -n assembly.fasta > amrfinder.output.tsv

# run on some reads
amrfinder -n someReads.fasta > amrfinder.output.tsv

# oh no, I have fastq.gz files, not fasta files
zcat reads.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' > reads.fasta
amrfinder -n reads.fasta > amrfinder.output.tsv
```

If you've already annotated your genome and have a bunch of predicted protein sequences, 
you can feed these into amrfinder instead.

```shell
amrfinder -p proteins.fasta > amrfinder.output.tsv
```

If you have an annotated genome, you can feed in both the protein and nucleotide sequences.
HOWEVER, amrfinder will want to know where the protein sequences belong on the genome.
You feed it this information with a GFF file, which specifies the genome coordinates of all proteins.
The problem you will probably encounter with GFF files in general is that they are often not formatted correctly,
or different programs will want them formatted in slightly different ways.
Assuming you have a GFF file that AMRFinder is happy with:

```shell
amrfinder -p proteins.fasta -n assembly.fasta -g annotations.gff > amrfinder.output.tsv
```

## Output

The output of AMRFinder is a tab separated file (TSV).
Everything should be TSV.
The output looks like this:

```text
Protein identifier  Contig id            Start  Stop  Strand  Gene symbol  Sequence name                                             Scope  Element type  Element subtype  Class                Subclass               Method               Target length  Reference sequence length  % Coverage of reference sequence  % Identity to reference sequence  Alignment length  Accession of closest sequence  Name of closest sequence                                  HMM id  HMM description
NA                  ERR1356692:1:261704  1      123   -       dfrB5        trimethoprim-resistant dihydrofolate reductase DfrB5      core   AMR           AMR              TRIMETHOPRIM         TRIMETHOPRIM           PARTIAL_CONTIG_ENDX  41             78                         52.56                             92.68                             41                WP_063844482.1                 trimethoprim-resistant dihydrofolate reductase DfrB5      NA      NA
NA                  ERR1543789:16:20590  3      251   -       lnu(C)       lincosamide nucleotidyltransferase Lnu(C)                 core   AMR           AMR              LINCOSAMIDE          LINCOSAMIDE            PARTIAL_CONTIG_ENDX  83             164                        50.61                             96.39                             83                WP_063851341.1                 lincosamide nucleotidyltransferase Lnu(C)                 NA      NA
NA                  ERR1301981:5:3301    2      250   -       lnu(C)       lincosamide nucleotidyltransferase Lnu(C)                 core   AMR           AMR              LINCOSAMIDE          LINCOSAMIDE            PARTIAL_CONTIG_ENDX  83             164                        50.61                             98.80                             83                WP_063851341.1                 lincosamide nucleotidyltransferase Lnu(C)                 NA      NA
NA                  ERR1664649:16:12246  3      299   -       erm(37)      23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(37)   core   AMR           AMR              MACROLIDE            MACROLIDE              PARTIAL_CONTIG_ENDX  99             179                        55.31                             100.00                            99                WP_003900446.1                 23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(37)   NA      NA
NA                  ERR1664649:21:8743   3      299   -       erm(37)      23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(37)   core   AMR           AMR              MACROLIDE            MACROLIDE              PARTIAL_CONTIG_ENDX  99             179                        55.31                             100.00                            99                WP_003900446.1                 23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(37)   NA      NA
NA                  ERR1664649:25:49979  1      300   -       aac(2')-Ic   aminoglycoside N-acetyltransferase AAC(2')-Ic             core   AMR           AMR              AMINOGLYCOSIDE       GENTAMICIN/TOBRAMCYIN  PARTIAL_CONTIG_ENDX  100            181                        55.25                             100.00                            100               WP_003899880.1                 aminoglycoside N-acetyltransferase AAC(2')-Ic             NA      NA
NA                  ERR1356738:1:295758  8      124   -       dfrB         trimethoprim-resistant dihydrofolate reductase DfrB       core   AMR           AMR              TRIMETHOPRIM         TRIMETHOPRIM           PARTIAL_CONTIG_ENDX  39             78                         50.00                             94.87                             39                WP_063844479.1                 trimethoprim-resistant dihydrofolate reductase DfrB6      NA      NA
NA                  ERR1467153:1:136128  13     243   -       qacH         quaternary ammonium compound efflux SMR transporter QacH  core   STRESS        BIOCIDE          QUATERNARY AMMONIUM  QUATERNARY AMMONIUM    PARTIAL_CONTIG_ENDX  77             125                        61.60                             89.61                             77                WP_000121134.1                 quaternary ammonium compound efflux SMR transporter QacH  NA      NA
NA                  ERR1467153:2:333502  2      214   +       qacL         quaternary ammonium compound efflux SMR transporter QacL  core   STRESS        BIOCIDE          QUATERNARY AMMONIUM  QUATERNARY AMMONIUM    PARTIAL_CONTIG_ENDX  71             110                        64.55                             97.18                             71                WP_015060824.1                 quaternary ammonium compound efflux SMR transporter QacL  NA      NA
```

As you can see this shows you what predicted AMR genes etc are identified in your samples,
and importantly the sequence IDs and coordinates.
You can look at the alignment metrics to help you decide if you believe the gene is there and in tact etc.

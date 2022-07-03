---
layout: default
permalink: /data_viz_part_1
---

# Data visualisation Part 1

__R tutorial tested with R version 4.1.0 (2021-05-18) and Rstudio version 2022.02.3 Build 492 on Windows 10__

__Python tutorial tested with PyCharm Pro 2020.3 with Python version 3.9 on Windows 10__

In this session you will learn how to import and manage your data tables into R using the tidyverse packages tidyr and dplyr.
You will find [this tidyr/dplyr cheatsheet](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf) handy when using these packages.
We will also walk through some equivalent operations in Python using Pandas.
Here is a [handy Pandas cheat sheet](https://pandas.pydata.org/Pandas_Cheat_Sheet.pdf) for future reference.

# Download some data

Lets start with a superfocus run.
Download the files:
- superfocus output file [output_all_levels_and_function.tsv.gz](/workshop2022/files/output_all_levels_and_function.tsv.gz)
- metadata file [metadata.tsv](/workshop2022/files/metadata.tsv)

# Open RStudio

It's time to load the data.
Set your working directory to wherever your downloaded tables are under _Session_ > _Set working directory_ > _Choose directory..._
and install (if needed) and load the libraries that you will be using.

```r
install.packages(c('tidyr','dplyr', 'ggplot2'))
library(tidyr)
library(dplyr)
library(ggplot2)
```

Load the data into new tables with read.csv().
View your dataframes with View().

```r
df = read.csv('output_all_levels_and_function.tsv.gz',
              header=T,
              sep='\t')
meta = read.csv('metadata.tsv',
                header=T,
                sep='\t')
View(df)
View(meta)
```

# Wrangling your tables with Tidyr and Dplyr

Let's remove the raw counts and just use the % columns.
We can do this with either base R or dplyr select() with a helper function.

```r
df_perc = df[c(1:4,9:12)]
# or
df_perc = select(df, -matches("fastq$"))
```

We want to convert the table into a long format table, where one row is one observation.
We do this with tidyr's gather().
The last four columns will be gathered into two columns, one for the sample name, and one for the perc value.
The first four columns will be duplicated for each sample.

```r
df_long = gather(df_perc, "SampleName", "PercentReads", 5:8)
```

This long format is the format that most Tidyverse packages will prefer.

We have some metadata for our samples that we want to use with our analysis.
Typically, we would merge this into our table.
I could write a whole tutorial on tables merging, but, 
as long as the index IDs are exactly the same in both tables it's fine.

Our problem is that our samples are named like "good_out_R1.SRR1237780_good_out_R1.fastq" in the superfocus output,
whereas the metadata only has the SRA accessions like "SRR1237780".
Rename the samples in the long table to match the metadata table.
We can do this easily with mutate to create a new column, then drop the old sample name column.

```r
df_long_clean = df_long %>% 
  mutate(SampleID = gsub('good_out_R1.|_good_out_R1.fastq..','',SampleName)) %>% 
  select(-SampleName)
```

Now we can merge in the metadata.

```r
# if the column heading are different
df_meta = merge(df_long_clean, meta, by.x = 'SampleID', by.y= 'SampleID')
# if they're the same you can just do this
df_meta = merge(df_long_clean, meta, by= 'SampleID')
```

# Collect some summaries etc.

Dplyr has some excellent functions for grouping data in the data frames and generating useful information.
We can pool the upstream and downstream samples, find the average, sort by abundance, filter out low abundance hits etc.

```r
top_group_perc = df_meta %>%
  group_by(Subsystem.Level.1, Group) %>% 
  summarise(mean = mean(PercentReads)) %>%
  arrange(desc(mean)) %>%
  filter(mean > 0.01)
```

We can then plot the result.

```r
ggplot(top_group_perc,
       aes(x=Subsystem.Level.1,
           y=mean,
           fill=Group)) +
  geom_bar(stat='identity', 
           position = 'dodge') +
  coord_flip()
```

# Open PyCharm and redo everything in Python

Create a new project for the workshop, and initialise the project with a virtual environment (venv).

Install that packages we need (at the moment just Pandas).
_File_ > _Settings_ > _Project: yourProjectName_ > _Python interpreter_ > click the plus sign bottom left to add packages.
Search for _Pandas_ and install it.
Now do the same for _Seaborn_.

You can work out of the Python Console, 
but it's probably a better idea to create a python file to work from as this will save all your commands.
_Select your base directory_, _File_ > _New_ > _Python File_ > give it a name.

Open the python file and load pandas (in PyCharm, execute a line of code with ALT + SHIFT + E)
Read in with read_csv() and you can view with .head(), .tail(), 
or in PyCharm you can right click and view as dataframe for the _df_ and _meta_ tables directly.

```python
import pandas as pd
import seaborn as sns
import re
import numpy as np

df = pd.read_csv('output_all_levels_and_function.tsv.gz',
                 compression='gzip',
                 header=0,
                 sep='\t')
meta = pd.read_csv('metadata.tsv',
                   header=0,
                   sep='\t')
df.head()
meta.head()
```

To drop the counts and keep the percentage columns:

```python
df_perc = df.iloc[0,1,2,3,8,9,10,11]

# you can make this look neater with numpy
df_perc = df.iloc[np.r_[0:4,8:12]]
```

convert to long format; you have to specify columns by name, 
so using df_perc.columns[0:4].to_list() is easier than typing it out.

```python
df_long = df_perc.melt(
    id_vars=df_perc.columns[0:4].to_list(),
    value_vars=df_perc.columns[4:8].to_list(), 
    var_name='SampleID', 
    value_name="PercentReads")
```

Clean the sample IDs:

```python
df_long['SampleID'] = [re.sub('good_out_R1/|_good_out.*', '', i) for i in df_long['SampleID']]
```

merge in the metadata:

```python
# because they both have columns named "SampleID" we can just do this
df_meta = pd.merge(df_long, meta)

# if the column names are different you need to specify them
df_meta = pd.merge(df_long, meta, left_on='SampleID', right_on='SampleID')
```

collect some summaries:

```python
top_group_perc = df_meta.groupby(by=["Subsystem Level 3", "Group"], as_index=False).agg('mean').query('PercentReads > 0.05').sort_values('PercentReads')

# same as above but easier to read
top_group_perc = df_meta.groupby(by=["Subsystem Level 3", "Group"], as_index=False) \
                        .agg('mean') \
                        .query('PercentReads > 0.05') \
                        .sort_values('PercentReads')

# you can also filter like so
top_group_perc = top_group_perc[top_group_perc['PercentReads']>0.05]
```

Make a nice figure:

```python
sns.barplot(data=top_group_perc, 
            y='Subsystem Level 3', 
            x='PercentReads', 
            hue='Group')
```


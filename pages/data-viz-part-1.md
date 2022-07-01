---
layout: default
---

# Data visualisation Part 1

In this session you will learn how to import and manage your data tables into R using the tidyverse packages tidyr and dplyr.
You will find [this tidyr/dplyr cheatsheet](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf) handy when using these packages.
We will also walk through some equivalent operations in Python using Pandas.
Here is a [handy Pandas cheat sheet](https://pandas.pydata.org/Pandas_Cheat_Sheet.pdf) for future reference.

# Download some data

Lets start with a superfocus run.
Download the files:
- superfocus output file [output_all_levels_and_function.tsv.gz](/files/output_all_levels_and_function.tsv.gz)
- metadata file [metadata.tsv](/files/metadata.tsv)

# Open RStudio

It's time to load the data.
Set your working directory to wherever your downloaded tables are under _Session_ > _Set working directory_ > _Choose directory..._
and install (if needed) and load the libraries that you will be using.

```r
install.packages(c('tidyr','dplyr'))
library(tidyr)
library(dplyr)
```

Load the data into new tables with read.csv().
View your dataframes with View().

```r
df = read.csv('output_all_levels_and_function.tsv.gz',header=T,sep='\t')
meta = read.csv('metadata.tsv',header=T,sep='\t')
View(df)
View(meta)
```

# Open PyCharm and do everything again in Python

Create a new project for the workshop, and initialise the project with a virtual environment (venv).

Install that packages we need (at the moment just Pandas).
_File_ > _Settings_ > _Project: yourProjectName_ > _Python interpreter_ > click the plus sign bottom left to add packages.
Search for _Pandas_ and install it.

You can work out of the Python Console, 
but it's probably a better idea to create a python file to work from as this will save all your commands.
_Select your base directory_, _File_ > _New_ > _Python File_ > give it a name.

Open the python file and load pandas (in PyCharm, execute a line of code with ALT + SHIFT + E)
Read in with read_csv() and you can view with .head(), .tail(), or in PyCharm you can open and view _df_ and _meta_ tables directly.

```python
import pandas as pd

df = pd.read_csv('output_all_levels_and_function.tsv.gz',compression='gzip',header=0,sep='\t')
meta = pd.read_csv('metadata.tsv',header=0,sep='\t')
df.head()
meta.head()
```


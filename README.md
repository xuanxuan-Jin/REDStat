# REDStat
REDStat(RNA Editing Detection and Statistics) is a tool developed to investigate RNA editing patterns using R scripts.
In this package, we provide tools including RNA_seq data analysis,statistical analysis of RNA editing events.

# Requirements
In order to run these scripts, you will need to install R on your computer. We also suggest using RStudio as an IDE.

  Downloading and Installing R - http://cran.us.r-project.org/  
  Downloading and Installing RStudio - http://www.rstudio.com/products/rstudio/download/  
  
  R >= 4.2  
  
Once you have your installation up and running, you will need to make sure you have some packages installed. The required packages are:
ShortRead  
Biostrings  
parallel  
dplyr  
data.table  

# Workflow  
FASTQ files  
↓  
Barcode extraction  
↓  
Read pairing  
↓  
Construct reconstruction  
↓  
Editing rate calculation  
# How to use?
## Scripts

| Script | Function |
|---------|---------|
| 01_fastq_process.R | FASTQ preprocessing |
| 02_editing_rate_count.R | Editing rate quantification |

## Usage

```R
process_fastq_files(
    file1,
    file2,
    output_path
)
```

## Output

```text
df_A_POS_part.csv
```

Contains A-to-G editing rates for all target positions.

这是一个符合生信项目规范、适合放在博士论文代码仓库里的 README 模板。

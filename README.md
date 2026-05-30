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
L1<-'TGAAGACTAATAGTAGGACCACGAGAAGTACTATTACGACCAACAGAAACAGTATCACGAGTAATACGAGCAATAGGAGCAGGACTACTAGCATAAACAGGACGATCACAACGACAAGAATGATTAGTACGAGGAGGATCATGATTAGGAT'   ###mRNA靶向序列
  
    fixed_start= "GAATTC"  # 固定起始序列
    rev= "CAACGC"  # 反向序列标记
    umi1= "CGACT"  # UMI1 条件
    umi2= "GCTCAAC" # UMI2 条件
    linker_left= "CCAGGC"  # 左侧 linker 标记

file1 <- 'sample_R1.fastq'  # Constructing file1 name
file2 <- 'sample_R2.fastq'  # Constructing file2 name
output_path <- 'sample'  # Assigning output_path as a character string
design_seq_const <- L1

rna_stat_editing(file1,file2,output_path,design_seq_const)

## Output

freq.csv   ###每种barcode对应的序列频率
seq.csv    ###靶向序列
df_A_POS_part.csv  ###每个barcode对应的双链中靶向序列上靶位点的编辑效率
```

Contains A-to-G editing rates for all target positions.
# Contact
Xuanxuan Jin
jean_jxx@163.com

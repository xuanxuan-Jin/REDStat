# 1. FASTQ预处理
# Load required libraries
library(ShortRead)
library(Biostrings)
library(parallel)
library(data.table)

process_fastq_read <- function(fastq_file, output_path, output_file_name1,  output_file_name2,chunk_size = 1e5) {
  #  ####处理1
  streamer <- FastqStreamer(fastq_file, n = chunk_size)
  max_read1 <- 0
  exclude_reads1<-c()
  # 确保输出文件夹存在
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
   }
  
    # 定义输出文件的完整路径
    output_file1 <- file.path(output_path, output_file_name1)
  
    repeat {
    batch <- yield(streamer)
    if (length(batch) == 0) break  # 如果没有更多数据，退出循环
    
    # 转换为数据框并添加序号
    R1_sequence <- as.data.frame(sread(batch))
    R1_sequence$number_read1 <- max_read1 + (1:nrow(R1_sequence))
    max_read1 <- max(R1_sequence$number_read1)
    
    # 筛选包含 fixed_start 的序列
    R1_sequence_1 <- R1_sequence[grepl(fixed_start, R1_sequence$x), ]
    R1_sequence_1$constructR1 <- gsub(paste0(".*", fixed_start), "", R1_sequence_1$x)

    # 重命名列
    colnames(R1_sequence_1) <- c("non_processed_R1", "number_read1", "constructR1")
    
    # 过滤长度大于等于 110 的序列
    R1_sequence_1$length <- nchar(R1_sequence_1$constructR1)
    R1_sequence_1 <- R1_sequence_1[which(nchar(R1_sequence_1$constructR1) >= 110), ]
    
    # 提取桥接序列
    R1_sequence_1$brige <- substr(R1_sequence_1$constructR1, start = 115, stop = 125)
    R1_sequence_1 <- R1_sequence_1 %>% select(non_processed_R1, number_read1, constructR1, brige)

    exclude_reads1<-c(exclude_reads1,R1_sequence_1$number_read1)
    # 将结果写入 CSV 文件，支持追加
    write.table(R1_sequence_1, output_file1, sep = ',', col.names = !file.exists(output_file1), append = TRUE, row.names = FALSE)
  }
  
  # 关闭 FastqStreamer
  close(streamer)
  print('Processing done for read1')

  ####处理2
  streamer <- FastqStreamer(fastq_file, n = chunk_size)
  
  max_read2 <- 0

  # 确保输出文件夹存在
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # 定义输出文件的完整路径
  output_file2 <- file.path(output_path, output_file_name2)
  
  repeat {
    batch <- yield(streamer)
    if (length(batch) == 0) break  # 如果没有更多数据，退出循环

    # 反向互补序列，并转换为数据框
    R1_sequence_reverse_comp <- as.data.frame(sread(reverseComplement(batch)))
    R1_sequence_reverse_comp$number_read2 <- max_read2 + (1:nrow(R1_sequence_reverse_comp))
    max_read2 <- max(R1_sequence_reverse_comp$number_read2)

    # 筛选出不在 exclude_reads1 中的序列
    R1_sequence_2 <- R1_sequence_reverse_comp[!R1_sequence_reverse_comp$number_read2 %in% exclude_reads1, ]

    # 筛选包含特定序列的行
    R1_sequence_2 <- R1_sequence_2[grepl(rev, R1_sequence_2$x), ]  # 存在M13 rev

    # 处理barcode
    R1_sequence_2$barcode <- gsub(paste0(".*",umi1), "", R1_sequence_2$x)    # 前面删去
    R1_sequence_2$barcode <- gsub(paste0("\\", umi2, ".*"), "", R1_sequence_2$barcode)    # 后面删去
    R1_sequence_2 <- R1_sequence_2[nchar(R1_sequence_2$barcode) == 20, ]  # 长度为20的条形码
    R1_sequence_2$barcode <- paste(substr(R1_sequence_2$barcode, 1, 9), 
                                     substr(R1_sequence_2$barcode, 12, 20), sep = '')

    R1_sequence_2$barcode1<-gsub(paste0(".*",barcode1), "", R1_sequence_2$x)
    R1_sequence_2$barcode1<-substr(R1_sequence_2$barcode1, 1, 5)

    R1_sequence_2$barcode  <-paste0(R1_sequence_2$barcode,'_',R1_sequence_2$barcode1)
    # 筛选包含CCAGGC的序列
    R1_sequence_2 <- R1_sequence_2[grep(linker_left, R1_sequence_2$x), ]
    
    # 提取constructR2
    R1_sequence_2$constructR2 <- gsub(paste0("\\", linker_left, ".*"), "", R1_sequence_2$x)
    R1_sequence_2 <- R1_sequence_2[nchar(R1_sequence_2$constructR2) >= 35, ]

    # 将结果写入 CSV 文件，支持追加
    write.table(R1_sequence_2, output_file2, sep = ',', col.names = !file.exists(output_file2), append = TRUE, row.names = FALSE)
  }
  
  # 关闭 FastqStreamer
  close(streamer)
  print('Processing done for read2')
}
process_fastq_files <- function(file1, file2, file_path) {
  # 处理第一个文件
  output_file_name1 <- "seq1_1.csv"             # 输出文件名
  output_file_name2 <- "seq1_2.csv" 
  process_fastq_read(file1,output_path,output_file_name1,output_file_name2)
  
  # 处理第二个文件
  output_file_name1 <- "seq2_2.csv"             # 输出文件名
  output_file_name2 <- "seq2_1.csv" 
  process_fastq_read(file2,output_path,output_file_name1,output_file_name2)
}
# 定义处理函数
process_merge <- function(file1, file2, output_file, output_path,chunk_size = 1e6) {
  # 读取第二个文件并设置索引
  file1 <- file.path(output_path, file1)
  file2 <- file.path(output_path, file2)
  columns_to_select <- c('barcode', 'number_read2', 'constructR2')
  seq2 <- fread(file2, sep = ',', header = TRUE, select = columns_to_select)
  setDT(seq2)
  setkey(seq2, number_read2)
  
  seq2[, number_read2 := as.character(number_read2)]
  
  # 打开第一个文件进行分块读取
  con <- file(file1, open = "r")
  
  # 读取文件的列名
  col_names <- c("non_processed_R1", "number_read1", "constructR1", "brige")
  # 初始化结果数据框
  result <- NULL
  output_file <- file.path(output_path, output_file)
  repeat {
    # 读取一块数据
    chunk <- read.table(con, sep = ",", header = FALSE, nrows = chunk_size, col.names = col_names, fill = TRUE)
    chunk <- chunk[, c("number_read1", "constructR1", "brige")]
    
    # 检查是否读取到数据
    if (nrow(chunk) == 0) {
      break
    }
    
    chunk$number_read1 <- as.character(chunk$number_read1)
    
    # 过滤数据，合并两个数据集
    setDT(chunk)
    result_chunk <- merge(chunk, seq2, by.x = "number_read1", by.y = "number_read2", all.x = TRUE)

    # 写入结果到文件
    fwrite(result_chunk, output_file, sep = ',', row.names = FALSE, col.names = !file.exists(output_file), append = TRUE)
  }
  
  # 关闭文件连接
  close(con)
  
  # 清理内存
  rm(result, chunk, result_chunk)
  gc()
  
  print('merge data')
}
process_fastq_merge1 <- function(file_path) {
  # 处理第一个文件
  file1 <- "seq1_1.csv"             # 输出文件名
  file2 <- "seq2_1.csv" 
  output_file <- 'merge_lib1.csv'  # 输出文件路径
  process_merge(file1,file2,output_file,output_path=output_path)

  # 处理第二个文件
  file1 <- "seq2_2.csv"             # 输出文件名
  file2 <- "seq1_2.csv"
  output_file <- 'merge_lib2.csv'  # 输出文件路径
  process_merge(file1,file2,output_file,output_path=output_path)

}

process_whole <- function(input_file, output_file,output_path, chunk_size = 1e4, max_cores = 8) {
  # 打开文件进行分块读取
  input_file <- file.path(output_path, input_file)
  output_file <- file.path(output_path, output_file)
  con <- file(input_file, open = "r")
  
  # 读取文件的列名
  col_names <- c('number_read1', 'constructR1', 'brige', 'barcode', 'constructR2')
  
  repeat {
    # 读取一块数据
    chunk <- read.table(con, sep = ",", header = FALSE, nrows = chunk_size, col.names = col_names, fill = TRUE)
    chunk <- chunk[, c('constructR1', 'brige', 'barcode', 'constructR2')]
    chunk$constructR1 <- substr(chunk$constructR1, start = 1, stop = 114)
    
    # 检查是否读取到数据
    if (nrow(chunk) == 0) {
      break
    }
    
    # 设置并行处理的核心数目
    cores <- min(max_cores, detectCores())
    
    # 使用 mclapply 进行并行处理
    constructR2_1 <- mclapply(
      1:nrow(chunk), mc.cores = cores, function(x) {
        tryCatch({
          # 检查 brige 是否在 constructR2 中
          if (grepl(chunk$brige[x], chunk$constructR2[x])) {
            # 删除包含 brige 及其之前的序列
            construct <- gsub(paste0(".*", chunk$brige[x]), "", chunk$constructR2[x])
            return(construct)
          } else {
            return(NA)  # 如果 brige 不在 constructR2 中，返回 NA
          }
        }, error = function(e) {
          cat("Error in processing row", x, ": ", conditionMessage(e), "\n")
          return(NA)  # 在发生错误时返回 NA
        })
      }
    )
    
    # 将结果添加到数据框中
    chunk$constructR2_1 <- unlist(constructR2_1)

    # 创建 Whole_Construct 时处理 NA 情况
    chunk$Whole_Construct <- paste0(
      chunk$constructR1, 
      chunk$brige, 
      ifelse(is.na(chunk$constructR2_1), "", chunk$constructR2_1)
    )
    
    # 根据 Whole_Construct 的长度进行过滤
    chunk <- chunk[nchar(chunk$Whole_Construct) == 151, ]
    chunk <- chunk[, c('barcode', 'Whole_Construct')]
    
    # 写入处理后的数据到文件
    write.table(chunk, output_file, sep = ',', col.names = !file.exists(output_file), row.names = FALSE, append = TRUE)
  }
  
  # 关闭文件连接
  close(con)

  # 清理内存
  rm(chunk)
  gc()
  
  print('merge data')
}
###删去中间文件
process_fastq_merge2 <- function(output_path) {
  # 处理第一个文件
  # 示例用法
  input_file <- "merge_lib1.csv"  # 输入文件路径
  output_file <- 'whole_1.csv'  # 输出文件路径
  process_whole(input_file,output_file,output_path=output_path)
  # 处理第二个文件
  # 示例用法
  input_file <- "merge_lib2.csv"  # 输入文件路径
  output_file <- 'whole_2.csv'  # 输出文件路径
  process_whole(input_file,output_file,output_path=output_path)
  delete1<-file.path(output_path,'seq1_1.csv')
  delete2<-file.path(output_path,'seq1_2.csv')
  delete3<-file.path(output_path,'seq2_1.csv')
  delete4<-file.path(output_path,'seq2_2.csv')
  delete5<-file.path(output_path,'merge_lib1.csv')
  delete6<-file.path(output_path,'merge_lib2.csv')
    file.remove(delete1)
    file.remove(delete2)
    file.remove(delete3)
    file.remove(delete4)
    file.remove(delete5)
    file.remove(delete6)
}
#### 合并拆分好的序列
process_data1 <- function(output_path,output_file='seq.csv'){
    output_path <- output_path  # 假设输出路径
    file1 <- "whole_1.csv"  # 第一个文件名
    file2 <- "whole_2.csv"  # 第二个文件名
    # 构建文件路径
    file1 <- file.path(output_path, file1)
    file2 <- file.path(output_path, file2)
    output_file <- file.path(output_path, output_file)

  # 读取两个文件并设置索引
  columns_to_select <- c('barcode', 'Whole_Construct')
  seq1 <- fread(file1, sep = ',', header = TRUE, select = columns_to_select)
  seq2 <- fread(file2, sep = ',', header = TRUE, select = columns_to_select)
  
  # 合并数据
  seq <- rbind(seq1, seq2)
  aa <- as.data.frame(table(seq$barcode))
  #aa <- aa[which(aa$Freq >= min_freq), ]
  aa1<-file.path(output_path, 'freq.csv')
  write.table(aa,aa1,sep=',',col.names=T,row.names=F)
  # 写入文件
  fwrite(seq, output_file, sep = ',', col.names = T, row.names = FALSE)
  }

rna_stat_seq<-function(file1,file2,output_path){
    #处理reads1,reads2
    process_fastq_files(file1,file2, output_path)
    # 合并数据
    process_fastq_merge1(output_path)
    process_fastq_merge2(output_path)
    # 计算编辑率
    process_data1(output_path)
}

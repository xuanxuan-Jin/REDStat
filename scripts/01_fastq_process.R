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

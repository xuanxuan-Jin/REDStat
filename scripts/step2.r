####
process_data2 <- function(design_seq_const,output_path,barcode,output_file){
  file1<-'seq.csv'
  file1<-file.path(output_path, file1)
  output_file <- file.path(output_path, output_file)
   # 读取两个文件并设置索引
   columns_to_select <- c('barcode', 'Whole_Construct')
   seq1 <- fread(file1, sep = ',', header = TRUE, select = columns_to_select)
   # 生成频率表
    
    barcode<-file.path(output_path, barcode)
    barcode_freq <- read.table(barcode,header=T,sep=',')
    barcode_freq<-barcode_freq[which(barcode_freq$Freq>=50),]
    process_barcode <- function(barcode_name) {
    tryCatch({
  
    data_matrix <- seq1[seq1$barcode == barcode_name, ]
    colnames(data_matrix) <- c("barcode", "Whole_Construct")
    
    # 创建 DNAStringSet 对象
    Sequence_barcode_StringSet <- DNAStringSet(data_matrix$Whole_Construct, use.names = TRUE)
    
    # 计算共识矩阵
    consensus_matrix_barcode <- t(consensusMatrix(Sequence_barcode_StringSet))[, c(1:4, 15:16)]
    
    consensus_matrix_barcode <- as.data.frame(consensus_matrix_barcode)
    consensus_matrix_barcode$position <- 1:nrow(consensus_matrix_barcode)
    consensus_matrix_barcode$A_to_G_rate <- (consensus_matrix_barcode$G / (consensus_matrix_barcode$A +      consensus_matrix_barcode$G)) * 100
    
    # 提取 A-to-G 率
    df_A_POS <- as.data.frame(t(as.data.frame(consensus_matrix_barcode[A_G_positions, "A_to_G_rate"])))
    df_A_POS$barcode_names <- barcode_name
    colnames(df_A_POS) <- c(paste0("A", A_G_positions), "barcode")
    fwrite(df_A_POS, output_file, sep = ',', col.names = !file.exists(output_file), append = TRUE, row.names = FALSE)
    return(df_A_POS)
    
  }, error = function(e) {
    message("Error processing barcode ", barcode_name, ": ", e$message)
    return(NULL)  # 如果发生错误，返回 NULL
  })
}

            

A_G_positions <- unlist(gregexpr('A', design_seq_const))
for (j in 1:dim(barcode_freq)[1]) {
  batch_results <- process_barcode(barcode_freq$Var1[j])
}
    
  print('Editing count done!'))
}

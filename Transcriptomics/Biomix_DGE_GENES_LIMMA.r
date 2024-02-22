#### INPUT PROMPT ----
# 
# 
# MANUAL INPUT
# library(vroom)
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"C1"
# args[2] <-"Control"
# args[3] <-"/home/henry/Desktop/BiomiX2.2"
# #
# directory <- args[3]
# iterations = 2
# i=2
# selection_samples = "YES"
# purity_filter = "NO"
# Cell_type = "Whole_RNA"
# setwd(directory[[1]])
# COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
# COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
# DIR_METADATA <- readLines("directory.txt")



if (length(args) ==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  if (args == "help") {
    print("CELL TYPE available: BLymphocytes / Monocytes / TLymphocytes / Neutrophils ....... DISEASES AVAILABLE: SjS / RA / SSc / PAPs / MCTD / UCTD / SLE ")
    stop()
  } else {
    stop("The number of argument is not enough, did you miss the cell type or the disease?")
  }
  
} else if (length(args)==3) {
  print("Correct number of argument :)")
  paste("Disease:", args[1])
}else if (length(args)> 3) {
  stop("Too many argument.. are you typing random words?") 
}

#### DATASET REARRANGEMENT ----

print("starting data upload and preparation")

library(vroom)
library(dplyr)
library(sjmisc)

directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
setwd(directory2)

print(iterations)

COMMAND$DIRECTORIES[i]
Matrix <-vroom(COMMAND$DIRECTORIES[i] , delim = "\t", col_names = TRUE)
Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)
COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
log2FC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[2])
Gene_panel <- COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3]
n_genes_heat<-as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])


genes <- Matrix$ID


if (selection_samples == "YES") {
  num <- grep(Cell_type, Metadata_total$SAMPLE_TYPE)
  Metadata_Bcell <- Metadata_total[num,]
  Metadata_total <- Metadata_total[num,]
  
  
  if (length(grep("[0-9].*_[A-Z].*",Metadata_Bcell$ID)) == length(Metadata_Bcell$ID)){
    IDs <- sapply(strsplit(Metadata_Bcell$ID, "_"), function(x) x[1])
    num <-which(colnames(Matrix) %in% Metadata_Bcell$ID)
  }else{num <-which(colnames(Matrix) %in% Metadata_Bcell$ID)}
  
  
  Matrix <- as.matrix(Matrix[,num])
  print(Metadata_Bcell$ID)
  
  if (length(grep("[0-9].*_[A-Z].*",Metadata_Bcell$ID)) == length(Metadata_Bcell$ID)){
    num <-which( Metadata_Bcell$ID %in% colnames(Matrix))
  }else{num <-which( Metadata_Bcell$ID %in% colnames(Matrix) )}
  
  
  Metadata_Bcell <- Metadata_Bcell[num,]
  #If multiple select only the first column
  Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
  
  Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
  Matrix<-Matrix[,order(colnames(Matrix))]
  
  if (length(grep("[0-9].*_[A-Z].*",Metadata_Bcell$ID)) == length(Metadata_Bcell$ID)){
    IDs <- sapply(strsplit(as.character(Metadata_Bcell$ID), "_"), function(x) x[1])
    IDs<- as.numeric(IDs)
    Metadata_Bcell$ID <- IDs
    colnames(Matrix) <- IDs
  }
  
} else {
  print("No samples selection")
  Metadata_Bcell <- Metadata_total
  Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
  Matrix <- as.matrix(Matrix[,-1])
  Matrix<-Matrix[,order(colnames(Matrix))]
  
  Metadata_Bcell <- Metadata_total
  num <-which(colnames(Matrix) %in% Metadata_Bcell$ID)
  Matrix <- as.matrix(Matrix[,num])
  #print(Metadata_Bcell$ID_CELL_TYPE)
  num <-which(Metadata_Bcell$ID  %in% colnames(Matrix) )
  Metadata_Bcell <- Metadata_Bcell[num,]
  #If multiple select only the first column
  Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
  
  Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
  Matrix<-Matrix[,order(colnames(Matrix))]
}

rownames(Matrix) <- genes
summary(as.factor(Metadata_Bcell$CONDITION))

Metadata_individual=NULL
Metadata_reads=NULL
Metadata=NULL
Metadata_tot = NULL
Metadata_total = NULL






#STATEMENT IF DATA ARE NORMALIZED
if (is_float(Matrix[-1])){
  NORMALIZATION <- "YES"
  
  library(DESeq2)
  setwd(directory2)
  
  if(sum(grep("ENSG", rownames(Matrix)[1])) == 1){
    Mart<-read.table(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt", sep=""), sep = ",", header=TRUE)
    Matrix <- as.data.frame(Matrix)
    Matrix$X <- genes
    DGE <- merge(Matrix, Mart, by.x = "X", by.y = "Gene.stable.ID")
    genes_name<-DGE$X
    genes_name_C <-DGE$Gene.name
    DGE2<-DGE[c(-1,-ncol(DGE))]
    rownames(DGE2) <- genes_name
    GENE_ANNOTATION <- "ENSEMBL"
  }else{
    Mart<-read.table(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt", sep=""), sep = ",", header=TRUE)
    Matrix <- as.data.frame(Matrix)
    Matrix$X <- genes
    DGE <- Matrix
    genes_name<-DGE$X
    DGE$Gene.name <- DGE$X
    genes_name_C <-DGE$Gene.name
    DGE2<-as.matrix(DGE[c(-(ncol(DGE)-1),-ncol(DGE))])
    rownames(DGE2) <- genes_name
    DGE2 <- as.data.frame(DGE2)
    GENE_ANNOTATION <- "GENE_NAME"
  }
  
  print(paste("Automatic detection of normalized data. Are the data normalized? ",NORMALIZATION))
  DGE2 <- DGE2[,order(colnames(DGE2))]
  Metadata_Bcell <- arrange(Metadata_Bcell,ID)
  
  num <- Metadata_Bcell$CONDITION == args[1]
  #Choise disease + HC
  numero <- Metadata_Bcell$CONDITION == args[2]
  
  Metadata_Bcells <- Metadata_Bcell
  Metadata_Bcell <- Metadata_Bcell[num | numero,]
  DGE2 <- DGE2[, num | numero]
  
  print("X")
  
  
  #Purity from MCP
  library(dplyr)
  
  setwd(paste(directory,"/Transcriptomics/INPUT",sep=""))
  
  METADATA_FILT <- !is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])
  METADATA_FILT_INDEX <-grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))
  
  repetition = 0
  for (meta_filter in METADATA_FILT_INDEX){
    repetition <- repetition + 1 
    if (!is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])[repetition]){
      COLNAME<-as.character(COMMAND_ADVANCED[1,meta_filter])
      if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Numerical"){
        To_filter<-as.numeric(unlist(Metadata_Bcell[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,1)
        characters_to_remove <- c(">", "<", "=", " ")
        value_threshold <- as.numeric(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "<" = function(a, b) a < b,
                                      ">" = function(a, b) a > b,
                                      "=" = function(a, b) a == b,
                                      ">=" = function(a, b) a >= b,
                                      "<=" = function(a, b) a <= b,
                                      NA)
        
        Metadata_Bcell <- Metadata_Bcell[comparison_operator(To_filter, value_threshold),]
        DGE2 <- DGE2[,c(comparison_operator(To_filter, value_threshold))]
        
      }else if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Factors"){
        To_filter<- as.character(unlist(Metadata_Bcell[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,2)
        characters_to_remove <- c("!=", "==", " ")
        value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "==" = function(a, b) a == b,
                                      "!=" = function(a, b) a != b,
                                      NA)
        
        Metadata_Bcell <- Metadata_Bcell[comparison_operator(To_filter, value_threshold),]
        DGE2 <- DGE2[,c(comparison_operator(To_filter, value_threshold),TRUE)]
      }
      
    }
  }
  
  Metadata_Bcell$CONDITION
  summary(as.factor(Metadata_Bcell$CONDITION))
  
  ###
  
  DGE3 <- DGE2
  DGE3$Gene.name <- rownames(DGE3)
  
  dir.create(path = paste(directory,"/MOFA/INPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  directory2 <- paste(directory, "/MOFA/INPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="")
  
  setwd(directory2)
  write.table(x=Metadata_Bcell   , file= paste("Metadata_",Cell_type,"_",args[1],".tsv",sep = "")  ,sep= "\t", row.names = FALSE, col.names = T,  quote = FALSE)
  
  colnames(DGE3) <- Metadata_Bcell$ID
  DGE3 <- as.data.frame(DGE3)
  colnames(DGE3)[ncol(DGE3)] <-"Gene.name"
  
  #DGE2$Gene.name <- rownames(DGE2)
  setwd(directory2)
  write.table(x=DGE3  , file= paste(Cell_type,"_",args[1],"_vs_", args[2],"_normalized_vst.tsv",sep = "") ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
  
}else{
  NORMALIZATION <- "NO"
  print(paste("Automatic detection of normalized data. Are the data normalized? ",NORMALIZATION))
  
  
  #### NORMALIZATION ----
  
  
  #Normalization
  library(DESeq2)
  setwd(directory2)
  
  if(sum(grep("ENSG", rownames(Matrix)[1])) == 1){
    Mart<-read.table(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt", sep=""), sep = ",", header=TRUE)
    Matrix <- as.data.frame(Matrix)
    Matrix$X <- genes
    DGE <- merge(Matrix, Mart, by.x = "X", by.y = "Gene.stable.ID")
    genes_name<-DGE$X
    genes_name_C <-DGE$Gene.name
    DGE2<-DGE[c(-1,-ncol(DGE))]
    rownames(DGE2) <- genes_name
    GENE_ANNOTATION <- "ENSEMBL"
  }else{
    Mart<-read.table(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt", sep=""), sep = ",", header=TRUE)
    Matrix <- as.data.frame(Matrix)
    Matrix$X <- genes
    DGE <- Matrix
    genes_name<-DGE$X
    DGE$Gene.name <- DGE$X
    genes_name_C <-DGE$Gene.name
    DGE2<-as.matrix(DGE[c(-(ncol(DGE)-1),-ncol(DGE))])
    rownames(DGE2) <- genes_name
    DGE2 <- as.data.frame(DGE2)
    GENE_ANNOTATION <- "GENE_NAME"
  }
  
  DGE2 <- DGE2[,order(colnames(DGE2))]
  Metadata_Bcell <- arrange(Metadata_Bcell,ID)
  
  num <- Metadata_Bcell$CONDITION == args[1]
  #Choise disease + HC
  numero <- Metadata_Bcell$CONDITION == args[2]
  
  Metadata_Bcells <- Metadata_Bcell
  Metadata_Bcell <- Metadata_Bcell[num | numero,]
  DGE2 <- DGE2[, num | numero]
  
  print("X")
  
  
  #Purity from MCP
  library(dplyr)
  
  setwd(paste(directory,"/Transcriptomics/INPUT",sep=""))
  
  METADATA_FILT <- !is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])
  METADATA_FILT_INDEX <-grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))
  
  repetition = 0
  for (meta_filter in METADATA_FILT_INDEX){
    repetition <- repetition + 1 
    if (!is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])[repetition]){
      COLNAME<-as.character(COMMAND_ADVANCED[1,meta_filter])
      if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Numerical"){
        To_filter<-as.numeric(unlist(Metadata_Bcell[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,1)
        characters_to_remove <- c(">", "<", "=", " ")
        value_threshold <- as.numeric(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "<" = function(a, b) a < b,
                                      ">" = function(a, b) a > b,
                                      "=" = function(a, b) a == b,
                                      ">=" = function(a, b) a >= b,
                                      "<=" = function(a, b) a <= b,
                                      NA)
        
        Metadata_Bcell <- Metadata_Bcell[comparison_operator(To_filter, value_threshold),]
        DGE2 <- DGE2[,c(comparison_operator(To_filter, value_threshold))]
        
      }else if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Factors"){
        To_filter<- as.character(unlist(Metadata_Bcell[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,2)
        characters_to_remove <- c("!=", "==", " ")
        value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "==" = function(a, b) a == b,
                                      "!=" = function(a, b) a != b,
                                      NA)
        
        Metadata_Bcell <- Metadata_Bcell[comparison_operator(To_filter, value_threshold),]
        DGE2 <- DGE2[,c(comparison_operator(To_filter, value_threshold),TRUE)]
      }
      
    }
  }
  
  
  
  DGE2<- apply(DGE2,2,round)
  Metadata_Bcell$CONDITION
  summary(as.factor(Metadata_Bcell$CONDITION))
  
  
  dds <- DESeqDataSetFromMatrix(countData = DGE2,
                                colData = Metadata_Bcell,
                                design = ~CONDITION)
  dds
  
  y<-unlist(args[2])
  dds$CONDITION <- relevel(dds$CONDITION, ref = y)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  
  normalized_counts <- counts(dds, normalized=TRUE) #Normalization size factor
  normalized_counts2 <- vst(assay(dds)) #Normalization VST
  
  
  DGE2 <- normalized_counts
  DGE2 <- as.data.frame(DGE2)
  DGE2$Gene.name <- rownames(DGE2)
  DGE3 <- assay(dds)
  DGE3 <- as.data.frame(DGE3)
  DGE3$Gene.name <- rownames(DGE3)
  DGE4 <- normalized_counts2
  DGE4 <- as.data.frame(DGE4)
  DGE4$Gene.name <- rownames(DGE4)
  
  
  #### DATA INTEGRATION MATRIX ----
  
  #create directory
  dir.create(path = paste(directory,"/MOFA/INPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  directory2 <- paste(directory, "/MOFA/INPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="")
  
  setwd(directory2)
  write.table(x=Metadata_Bcell   , file= paste("Metadata_",Cell_type,"_",args[1],".tsv", sep = "")  ,sep= "\t", row.names = FALSE, col.names = T,  quote = FALSE)
  
  colnames(DGE2) <- Metadata_Bcell$ID
  DGE2 <- as.data.frame(DGE2)
  colnames(DGE2)[ncol(DGE2)] <-"Gene.name"
  
  colnames(DGE3) <- Metadata_Bcell$ID
  DGE3 <- as.data.frame(DGE3)
  colnames(DGE3)[ncol(DGE3)] <-"Gene.name"
  
  colnames(DGE4) <- Metadata_Bcell$ID
  DGE4 <- as.data.frame(DGE4)
  colnames(DGE4)[ncol(DGE4)] <-"Gene.name"
  
  
  #DGE2$Gene.name <- rownames(DGE2)
  setwd(directory2)
  write.table(x=DGE2  , file= paste(Cell_type,args[1],"vs", args[2],"normalized.tsv",sep = "_") ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
  write.table(x=DGE3  , file= paste(Cell_type,args[1],"vs", args[2], "unormalized.tsv",sep = "_")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
  write.table(x=DGE4  , file= paste(Cell_type,args[1],"vs", args[2], "normalized_vst.tsv",sep = "_")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
  
  
}

#### HEATMAP ----

#Genes interferon signature
directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
setwd(directory2)

if(COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] != "X"){
  genes <-vroom(Gene_panel, delim = "\t", col_names = TRUE)
  genes <- genes$GENES_FOR_SUBPOPULATION}

directory2 <- paste(directory, "/MOFA/INPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="")
setwd(directory2)

if (NORMALIZATION=="NO"){
  DGE <- counts(dds,normalized=TRUE)
}else{
  DGE <- DGE2     
}

print("OK")
DGE <- as.data.frame(DGE)
DGE$X <- rownames(DGE)

if(GENE_ANNOTATION == "GENE_NAME"){
  DGE$Gene.name <- DGE$X
  DGE<-DGE[,-grep("X", colnames(DGE))]
}else{
  DGE <- merge(DGE, Mart, by.x = "X", by.y = "Gene.stable.ID")        
}

print("OK2")
if (COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] != "X") {
  
  Panel=NULL
  for (i in genes){
    Panel <-append(Panel, which(DGE$Gene.name == i))
  }
  
  DGE <- DGE[Panel,]
  rownames(DGE)<- DGE$Gene.name
  if(GENE_ANNOTATION == "ENSEMBL"){
    DGE<- DGE[,c(-1,-ncol(DGE))]
    rownames(DGE) <- genes
  }else{
    DGE<- DGE[,c(-ncol(DGE))]
    rownames(DGE) <- genes  
  }
  
  IFN_score= as.data.frame(matrix(0, ncol = ncol(DGE), nrow = nrow(DGE)))
  colnames(IFN_score) <- colnames(DGE)
  rownames(IFN_score) <- genes
  iteration <- 0
  
  x <- Metadata_Bcell$CONDITION == args[2]
  HC <- DGE[,x]
  
  tmp=NULL
  for (gene in 1:nrow(DGE)){
    if (iteration > 0){ 
      IFN_score[iteration,] <- tmp}
    Mean_hd <- as.numeric(apply((HC[gene,]),1,mean))
    SD_hd <- as.numeric(apply(HC[gene,],1,sd))
    tmp=NULL
    iteration <- iteration + 1
    for (sample in 1:ncol(DGE)){
      value <- ((as.numeric(DGE[gene,sample]) - Mean_hd))/ SD_hd
      tmp<-append(tmp,value)
    }}
  IFN_score[length(Panel),] <- tmp
  
  #colnames(IFN_score)<- Metadata_Bcells$CONDITION
  
  # library(ComplexHeatmap)
  # Heatmap(IFN_score)
  # 
  # library(circlize)
  # col_fun = colorRamp2(c(-0.3, 0, 0.3), c("blue", "black", "yellow"))
  # col_fun(seq(-3, 3))
  # Heatmap(IFN_score,km =2, name = "IFN_score", col = col_fun, clustering_distance_rows = "pearson", clustering_method_rows= "complete", row_dend_width = unit(0.5, "cm"), column_names_gp = grid::gpar(fontsize = 6),
  #         row_names_gp = grid::gpar(fontsize = 8))
  
  
  #HEATMAP
  num<-Metadata_Bcell$CONDITION == args[2]
  num2<-Metadata_Bcell$CONDITION == args[1]
  
  IFN_score<-IFN_score[, num | num2]
  
  
  # #### ADD HEATMAP AFTER IFNPOS PREDICION ----
  
  
  # #### AUTOMATIC IDENTIFICATION POSITIVE - NEGATIVE ----
  library(stringr)
  
  z1 <- IFN_score > 2 
  
  colSums(z1, na.rm = TRUE) 
  positivity1 <- colSums(z1, na.rm = TRUE) >= as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING[1])
  
  z2 <- IFN_score > 1
  
  colSums(z2, na.rm = TRUE) 
  positivity2 <- colSums(z2, na.rm = TRUE) >= as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING[2])
  
  
  # z3 <- IFN_score > 0.5
  # 
  # colSums(z3, na.rm = TRUE) 
  # positivity3 <- colSums(z3, na.rm = TRUE) > 10
  # 
  positivity<- positivity1|positivity2
  
  
  
  Metadata_Bcell$condition <- positivity
  
  Metadata_Bcell$condition <- str_replace(Metadata_Bcell$condition, "TRUE", "pos")
  Metadata_Bcell$condition <- str_replace(Metadata_Bcell$condition, "FALSE", "neg")
  
  
  library(ComplexHeatmap)
  setwd(directory2)
  pdf(file= paste("Gene_panel_subgroups_", Cell_type,"_", args[1],"_vs_", args[2], ".pdf",sep=""))
  #dev.new()
  
  
  library(circlize)
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow"))
  col_fun(seq(-3, 3))
  p<- Heatmap(IFN_score,km =2, name = "SD_score", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],  clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
              row_names_gp = grid::gpar(fontsize = 8))
  print(p)
  
  #If you want to visualize the group
  colnames(IFN_score)<- Metadata_Bcell$CONDITION[ num |num2]
  
  library(circlize)
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow"))
  col_fun(seq(-3, 3))
  
  
  t<-c("CTRL" = "blue", "SLE" = "red")
  attr(t, "names")[1]<- args[2]
  attr(t, "names")[2]<- args[1]
  attr(t, "names")
  
  ha = HeatmapAnnotation(condition = Metadata_Bcell$CONDITION,
                         col = list(condition = t))
  p<- Heatmap(IFN_score,km =2, name = "SD_score", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1], clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
              row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
  
  print(p)
  
  library(circlize)
  col_fun = colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow"))
  col_fun(seq(-3, 3))
  
  ha = HeatmapAnnotation(condition = Metadata_Bcell$condition,
                         col = list(condition = c("neg" = "blue", "pos" = "red" )))
  p<- Heatmap(IFN_score,km =2, name = "SD_score", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1], clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
              row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
  
  print(p)
  
  colnames(IFN_score) <- Metadata$ID
  #Saved for annotation
  setwd(directory2)
  #write.table(Metadata_Bcell, file="Metadata_Bcells", sep = "\t", quote = FALSE, row.names = F)
  dev.off()
  
  
  # #### VALIDATION SIGNATURE ----
  
  library(dplyr)
  
  if ("MARKER" %in% colnames(Metadata_Bcell)){
    
    sink("Validation_MARKER_vs_GENE_PANEL.tsv")
    
    #Results classification
    print(paste("total number of", args[1], "samples =", nrow(Metadata_Bcell %>% filter(CONDITION == args[1])) , sep = " "))
    print(paste("total number of", args[1], "samples with IFN measurement", nrow(Metadata_Bcell %>% filter(CONDITION == args[1]) %>% filter(MARKER > 0 | MARKER < 0 )) , sep = " "))
    
    
    #Results classification
    print(paste("total number of ", args[2], " samples = ", nrow(Metadata_Bcell %>% filter(CONDITION == args[2])) , sep = ""))
    print(paste("total number of ", args[2], " samples with IFN measurement = ", nrow(Metadata_Bcell %>% filter(CONDITION == args[2]) %>% filter(MARKER > 0 | MARKER < 0 )) , sep = ""))
    
    
    #To the metadata table is important to add the column called condition,
    #here we must define the classification pos/neg based on the heatmap
    
    setwd(directory2)
    #write.table(x=Metadata  , file= "Metadata_condition"  ,sep= ",", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    # Metadata_cond <- vroom("Metadata_condition" , delim = ",", col_names = TRUE)
    # Metadata$condition <- Metadata_cond$condition
    
    #Setting the IFN value to 1 this is the neg/pos condition accuracy
    print(paste("Number of samples IFN pos correctly classified = ", nrow(Metadata_Bcell %>% filter(MARKER > 0 | MARKER < 0 ) %>% filter(CONDITION == args[1]) %>% filter(condition == "pos")) , "/",  nrow(Metadata_Bcell %>% filter(CONDITION == args[1]) %>% filter(MARKER > 1.25 )), sep = ""))
    
    #Setting the IFN value to 1 this is the neg/pos condition accuracy
    print(paste("Number of samples IFN neg correctly classified = ", nrow(Metadata_Bcell %>% filter(MARKER > 0 | MARKER < 0 ) %>% filter(CONDITION == args[1]) %>% filter(condition == "neg")) , "/",  nrow(Metadata_Bcell %>% filter(CONDITION == args[1]) %>% filter(MARKER < 1.25 )), sep = ""))
    
    #Setting the IFN value to 1 this is the neg/pos condition accuracy
    print(paste("Number of samples CTRL pos correctly classified = ", nrow(Metadata_Bcell %>% filter(MARKER > 0 | MARKER < 0 ) %>% filter(CONDITION == args[2]) %>% filter(condition == "pos")) , "/",  nrow(Metadata_Bcell %>% filter(CONDITION == args[2]) %>% filter(MARKER > 1.25 )), sep = ""))
    
    #Setting the IFN value to 1 this is the neg/pos condition accuracy
    print(paste("Number of samples CTRL neg correctly classified = ", nrow(Metadata_Bcell %>% filter(MARKER > 0 | MARKER < 0 ) %>% filter(CONDITION == args[2]) %>% filter(condition == "neg")) , "/",  nrow(Metadata_Bcell %>% filter(CONDITION == args[2]) %>% filter(MARKER < 1.25 )), sep = ""))
    
    sink()
    
  }else{
    print("No validation marker added")
  }
  
  #### SAVING NORMALIZED DATA ----
  
}else{
  print("No panel selected")
}

if (COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] != "X"){
if (COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING[3] != "YES") {
  
  Metadata<-Metadata_Bcell
  num <- which(Metadata$CONDITION == args[2] & Metadata$condition == "neg")
  Metadata$condition[num] <- "healthy"
  DGE2 <- as.data.frame(DGE2)
  DGE3 <- as.data.frame(DGE3)
  Metadata_all <- Metadata
  
  #colnames(DGE2) <- c(Metadata$OMICID, "gene_name")
  DGE2<-as.data.frame(t(DGE2))
  
  if (NORMALIZATION=="NO"){
          DGE2<-DGE2[-nrow(DGE2),]}else{print("")}
  #DGE2<-DGE2[-nrow(DGE2),]
  #ADD CENTER AND BATCH TO THE MODEL
  #DGE2$batch <- Metadata$Batch
  #DGE2$center <- Metadata$Center
  #colnames(DGE2) <- c(genes_name,"condition","Batch", "Center")
  DGE2$condition <- Metadata$condition
  
  colnames(DGE2) <- c(genes_name,"condition")
  rownames(DGE2) <- Metadata$ID
  DGE2$samples <- Metadata$ID  
  
}else{
  
  Metadata<-Metadata_Bcell
  num <- which(Metadata$CONDITION == args[2] & Metadata$condition == "neg")
  Metadata$condition[num] <- "healthy"
  num <- which(Metadata$CONDITION == args[2] & Metadata$condition == "pos")
  DGE2 <- as.data.frame(DGE2)
  DGE3 <- as.data.frame(DGE3)
  DGE2 <-DGE2[,-num]
  DGE3 <-DGE3[,-num]
  Metadata <- Metadata[-num,]
  Metadata_all <- Metadata
  
  #colnames(DGE2) <- c(Metadata$OMICID, "gene_name")
  DGE2<-as.data.frame(t(DGE2))
  
  if (NORMALIZATION=="NO"){
    DGE2<-DGE2[-nrow(DGE2),]}else{print("")}
  
  DGE2$condition <- Metadata$condition
  
  #ADD CENTER AND BATCH TO THE MODEL
  #DGE2$batch <- Metadata$Batch
  #DGE2$center <- Metadata$Center
  #colnames(DGE2) <- c(genes_name,"condition","Batch", "Center")
  colnames(DGE2) <- c(genes_name,"condition")
  rownames(DGE2) <- Metadata$ID
  DGE2$samples <- Metadata$ID
  
}}else{
        Metadata<-Metadata_Bcell
        DGE2 <- as.data.frame(DGE2)
        DGE3 <- as.data.frame(DGE3)
        Metadata_all <- Metadata
        
        #colnames(DGE2) <- c(Metadata$OMICID, "gene_name")
        DGE2<-as.data.frame(t(DGE2))
        
        if (NORMALIZATION=="NO"){
                DGE2<-DGE2[-nrow(DGE2),]}else{print("")}
        
        #ADD CENTER AND BATCH TO THE MODEL
        #DGE2$batch <- Metadata$Batch
        #DGE2$center <- Metadata$Center
        #colnames(DGE2) <- c(genes_name,"condition","Batch", "Center")
        colnames(DGE2) <- c(genes_name)
        rownames(DGE2) <- Metadata$ID
        DGE2$samples <- Metadata$ID
}

setwd(directory2)
write.table(x=Metadata   , file= paste("Metadata_",Cell_type ,"_", args[1],"_vs_",args[2],"_final.tsv",sep ="")  ,sep= "\t", row.names = FALSE, col.names = T,  quote = FALSE)
write.table(x=DGE3  , file= paste("gene_count_",Cell_type, "_", args[1], "_vs_", args[2],"_unnormalized.tsv", sep ="")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)



counts<- paste("gene_count_",Cell_type, "_", args[1], "_vs_", args[2],"_unnormalized.tsv", sep ="") 
Meta<- paste("Metadata_", Cell_type, "_", args[1],"_vs_",args[2],"_final.tsv",sep ="")


if (COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] != "X") {
  
  #### DGE ANALYSIS IFN POS vs CONTROL ----
  
  library(vroom)
  library(dplyr)
  
  setwd(directory2)
  
  DGE3 <-vroom(counts , delim = "\t", col_names = TRUE)
  Metadata <-vroom(Meta , delim = "\t", col_names = TRUE)
  
  rownames(DGE3) <- DGE3$Gene.name
  
  
  #colnames(DGE3) <- c(Metadata$OMICID, "gene")
  gene_name <- DGE3$Gene.name
  
  
  
  #Purity from MCP
  library(dplyr)
  
  setwd(paste(directory,"/Transcriptomics/INPUT",sep=""))
  
  
  METADATA_FILT <- !is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])
  METADATA_FILT_INDEX <-grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))
  
  repetition = 0
  for (meta_filter in METADATA_FILT_INDEX){
    repetition <- repetition + 1 
    if (!is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])[repetition]){
      COLNAME<-as.character(COMMAND_ADVANCED[1,meta_filter])
      if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Numerical"){
        To_filter<-as.numeric(unlist(Metadata[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,1)
        characters_to_remove <- c(">", "<", "=", " ")
        value_threshold <- as.numeric(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "<" = function(a, b) a < b,
                                      ">" = function(a, b) a > b,
                                      "=" = function(a, b) a == b,
                                      ">=" = function(a, b) a >= b,
                                      "<=" = function(a, b) a <= b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
        DGE3 <- DGE3[,c(comparison_operator(To_filter, value_threshold),TRUE)]
        
      }else if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Factors"){
        To_filter<- as.character(unlist(Metadata[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,2)
        characters_to_remove <- c("!=", "==", " ")
        value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "==" = function(a, b) a == b,
                                      "!=" = function(a, b) a != b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
        DGE3 <- DGE3[,c(comparison_operator(To_filter, value_threshold),TRUE)]
      }
      
    }
  }
  
  
  #FROM HERE
  DGE3 <- DGE3[,-ncol(DGE3)]
  
  #Metadata<- Metadata %>% arrange(condition)
  which(colnames(DGE3) %in% Metadata$ID )
  
  
  
  
  tmp = NULL
  for (i in 1:length(Metadata$ID)){
    x <-grep(Metadata$ID[i], colnames(DGE3))
    tmp<-append(tmp,x)
  }
  
  DGE3 <- DGE3[,tmp]
  
  colnames(DGE3) <- c(Metadata$ID)
  
  num <-Metadata$condition == "neg"
  
  DGE3 <- DGE3[,!num]
  Metadata <- Metadata[!num,]
  
  if (sum(Metadata$condition == "pos") > 2){
    
    if(NORMALIZATION=="YES"){
      
      library(limma)
      library(edgeR)
      
      rownames(DGE3) <- gene_name
      
      if ("AGE" %in% colnames(Metadata) & length(unique(Metadata$AGE)) > 1){
        
        
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          
          #=========================================
          
          CONDITION_<-as.factor(Metadata$CONDITION)
          AGE <- Metadata$AGE
          GENDER <- Metadata$GENDER
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_ + AGE + GENDER, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
          
        }else{
          CONDITION_<-as.factor(Metadata$CONDITION)
          AGE <- Metadata$AGE
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_ + AGE, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
          
        }
      }else{
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          
          CONDITION_<-as.factor(Metadata$CONDITION)
          GENDER <- Metadata$GENDER
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_ + GENDER, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
          
        }else{
          
          CONDITION_<-as.factor(Metadata$CONDITION)
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
        }}
      
      
      #create directory
      dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_pos_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      
      directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_pos_vs_",args[2], sep ="")
      
      setwd(directory2)
      
      dds<-dds[,c(2,1,6,3,4,5)]
      colnames(dds) <- c("AveExpr", "log2FoldChange", "B" , "t", "pvalue", "padj"  )
      
      write.table(dds,
                  file= paste(args[1], "_pos_", args[2],"_",Cell_type, ".tsv" ,sep =""),sep= "\t")
      
      Normalized_heatmap <- DGE3
      
      
    }else{
      
      
      
      #colnames(DGE3) <- c(Metadata$OMICID, "gene")
      rownames(DGE3) <- gene_name
      #DGE3$GENE <- gene_name
      #DGE3$GENE
      
      library("DESeq2")
      
      if ("AGE" %in% colnames(Metadata) & length(unique(Metadata$AGE)) > 1){
        
        Metadata$AGE[Metadata$AGE <= 35] <- "less30"
        Metadata$AGE[Metadata$AGE <= 50 & Metadata$AGE > 35] <- "35_50"
        Metadata$AGE[Metadata$AGE > 50] <- "more50"
        
        
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          
          #=========================================
          
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ GENDER + AGE + condition)
          
        }else{
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ AGE + condition)}
      }else{
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ GENDER + condition)
          
        }else{
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ condition)}}
      
      
      
      dds$condition <- relevel(dds$condition, ref = "healthy")
      
      #dds$condition <- relevel(dds$condition, ref = "neg")
      
      #crea l'oggetto di classe DESeqDataSet, quello utilizzato dal pacchetto
      ddsHTSeq2 <- dds
      
      
      dds <- DESeq(dds)
      
      res <- results(dds)
      res
      #genera l'espressione differeziale sulla base del log2 fold change
      #aggiungendo il p-value e il p-value adjusted
      
      resultsNames(dds)
      #Number of comparison available
      
      
      resOrdered <- res[order(res$padj),]
      #ordina i geni in base alla significatività del p_value adjusted
      
      summary(res)
      
      sum(res$padj < 0.1, na.rm=TRUE)
      #quanti p_value adjusted sono meno di 0.1?
      sum(res$padj > 0.1, na.rm=TRUE)
      #e maggiori di 0.1?
      
      #informazioni su variabili e test utilizzati
      mcols(res)$description
      
      
      #create directory
      dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_pos_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      
      directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_pos_vs_",args[2], sep ="")
      
      setwd(directory2)
      
      write.table(as.data.frame(resOrdered),
                  file= paste(args[1], "_pos_", args[2],"_",Cell_type,".tsv" ,sep =""),sep = "\t",row.names = TRUE)
      
      Normalized_heatmap <- as.data.frame(vst(assay(dds))) #Normalization VST
      
    }
    
    #permette di scrivere un file csv contenente tutti i risultati (in questo caso con pvalue adj ordinato)
    setwd(directory2)
    DGE <- read.table(file = paste(args[1], "_pos_", args[2],"_", Cell_type,".tsv", sep =""), sep="\t", header = TRUE)
    DGE$X <- rownames(DGE)
    DGE <- DGE %>% select(X, everything())
    
    dir.create(path = paste(directory2,"/Subpopulation_input/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    DGE3_modified <- DGE3
    DGE3_modified$ID <- rownames(DGE3_modified)
    DGE3_modified <- DGE3_modified %>% select(ID, everything())
    write.table(x=DGE3_modified   , file= paste(directory2,"/Subpopulation_input/",args[1], "_pos_and_", args[2],"_", Cell_type,"_counts.tsv", sep ="")  ,sep= "\t", row.names = FALSE, col.names = T,  quote = FALSE)
    write.table(x=Metadata   , file= paste(directory2,"/Subpopulation_input/",args[1], "_pos_and_", args[2],"_", Cell_type,"_metadata.tsv", sep ="")  ,sep= "\t", row.names = FALSE, col.names = T,  quote = FALSE)
    DGE3_modified = NULL
    gc()
    
    if(GENE_ANNOTATION == "GENE_NAME"){
      print(" ")
      colnames(DGE)[1] <- "Gene.name"
    }else{
      setwd(paste(directory,"/MOFA/x_BiomiX_DATABASE",sep=""))
      Mart<-read.table("mart_export_37.txt", sep = ",", header=TRUE)
      DGE <- merge(DGE, Mart, by.x = "X", by.y = "Gene.stable.ID")}
    
    #colnames(DGE)[7] <- "Gene.name"
    
    setwd(directory2)
    write.table(as.data.frame(DGE),
                file= paste(args[1], "_pos_",args[2],"_",Cell_type,".tsv", sep =""), sep = "\t", row.names = FALSE)
    
    
    setwd(directory2)
    
    if(NORMALIZATION=="NO"){
      tab<- as.data.frame(counts(dds,normalized=TRUE))
    }else{tab <- DGE3}
    
    tab$X<- rownames(tab)
    
    if(GENE_ANNOTATION == "GENE_NAME"){
      tab <- tab %>% select(X, everything())
    }else{
      tab <- merge(tab, Mart, by.x = "X", by.y = "Gene.stable.ID")
      #tab<-tab[,c(ncol(tab),3:ncol(tab)-1)]
      tab<-tab[,-ncol(tab)]}
    
    colnames(tab)[1] <- "NAME"
    tab$Description <- "NA"
    
    tab<-tab[,c(1,ncol(tab),3:ncol(tab)-1)]
    
    
    file.remove(paste(args[1], "pos_",args[2],"_",Cell_type, "_RAW.gct", sep =""))
    tab<- as.data.frame(tab)
    line="#1.2"
    write(line,file=paste(args[1], "pos_",args[2],"_",Cell_type, "_RAW.gct", sep =""),append=TRUE)
    line=c((nrow(tab)),(ncol(tab)-2))
    write(line,file=paste(args[1], "pos_",args[2],"_",Cell_type, "_RAW.gct", sep =""),append=TRUE,sep = "\t")
    write.table(tab, append =TRUE,
                file= paste(args[1], "pos_",args[2],"_",Cell_type, "_RAW.gct", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
    write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE),],
                file= paste(args[1], "pos_",args[2],"_",Cell_type, "_RAW_meta.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
    
    dis<-Metadata[order(Metadata$CONDITION),]
    
    write.table(dis[,c("ID","CONDITION")],
                file= paste(args[1], "pos_",args[2],"_",Cell_type, "_RAW_meta_groups.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
    #Ho preso i nomi dei diversi geni attraverso bio.mart da questo ho scaricato
    #la tabella e l'ho aggiunta a quella DGE.
    #
    #
    #
    ALTI <- subset(DGE, padj < padju)
    ALTI <-subset(ALTI, log2FoldChange > log2FC)
    ALTI <- ALTI[order(ALTI$padj), ]
    #
    BASSI <- subset(DGE, padj < padju)
    BASSI <- subset(BASSI, log2FoldChange < -log2FC)
    BASSI <- BASSI[order(BASSI$padj), ]
    NO <- subset(DGE,padj > padju)
    
    
    setwd(directory2)
    
    # write.csv(arrange(ALTI,padj),
    #           file="HC_vs_SLEpos_upregulated")
    #
    #
    # write.csv(arrange(BASSI,padj),
    #           file="HC_vs_SLEpos_downregulated")
    
    
    library(ggplot2)
    library(ggrepel)
    setwd(directory2)
    pdf(file=paste("plot_DGE_", args[1],"_pos_",args[2],"_", Cell_type,".pdf"))
    
    p <-ggplot() +
      ggtitle( paste("DGE", Cell_type,"_", args[1], "_POS_vs_", args[2],sep="")) + theme(
        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
      geom_point(data = ALTI, aes(x = log2FoldChange, y = -log10(padj)), color = "red") +
      geom_text_repel(data = ALTI[c(1:25), ], aes(x = log2FoldChange , y = -log10(padj), label=Gene.name),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1)+
      geom_point(data = BASSI, aes(x = log2FoldChange, y = -log10(padj)), color = "blue") +
      geom_text_repel(data = BASSI[1:25, ], aes(x = log2FoldChange , y = -log10(padj), label=Gene.name),hjust=0, vjust=0,size=3, force = 20, force_pull = 1)+
      geom_point(data = NO, aes(x = log2FoldChange, y = -log10(padj)), color = "black")
    
    p2 <- p + geom_vline(xintercept=c(-log2FC,log2FC), col="red") +
      geom_hline(yintercept=-log10(padju), col="red")
    
    
    print(p2)
    dev.off()
    
    
    setwd(directory2)
    
    write.table(x=ALTI$Gene.name   , file= paste(Cell_type,"_", args[1],"pos_",args[2],"_","UP_ENRICHR.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x=BASSI$Gene.name  , file= paste(Cell_type,"_",args[1],"pos_",args[2],"_","DOWN_ENRICHR.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x=ALTI   , file= paste(Cell_type,"_", args[1],"pos_",args[2],"_","GENESUP.tsv",sep="") ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    write.table(x=BASSI  , file= paste(Cell_type,"_", args[1],"pos_",args[2],"_","GENESDOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    
    
    if ((nrow(ALTI) + nrow(BASSI)) > 2){
      
      #ADD HEATMAP CREATION HERE
      if (nrow(ALTI) > n_genes_heat){
        ALTI <- ALTI[1:n_genes_heat,]}
      if (nrow(BASSI) > n_genes_heat){
        BASSI <- BASSI[1:n_genes_heat,]}
      
      #HEATMAP INPUT ARE NORMALIZED
      
      Heat<-Normalized_heatmap[rownames(Normalized_heatmap) %in% c(ALTI$X, BASSI$X),]
      to_add<-rownames(Normalized_heatmap)[rownames(Normalized_heatmap) %in% c(ALTI$X, BASSI$X)]
      rownames(Heat) <-to_add
      
      if (GENE_ANNOTATION == "ENSEMBL"){
        
        Heat$GENE <-rownames(Heat)
        Heat <- merge(Heat, Mart, by.x = "GENE", by.y = "Gene.stable.ID")
        Heat <- Heat[!duplicated( Heat[,ncol(Heat)]),]
        rownames(Heat) <- Heat[,ncol(Heat)]
        Heat <- Heat[,c(-1,-ncol(Heat))]
      }
      
      library(ComplexHeatmap)
      setwd(directory2)
      pdf(file= paste("Heatmap_top_genes", Cell_type, args[2],"vs", args[1],"pos.pdf",sep="_"))
      #dev.new()
      
      library(circlize)
      col_fun = colorRamp2(c(min(Heat) - 3, mean(colMeans(Heat)), max(Heat) - 3), c("blue", "black", "yellow"))
      col_fun(seq(-3, 3))
      
      t<-c("CTRL" = "blue", "SLE" = "red")
      attr(t, "names")[1]<- args[2]
      attr(t, "names")[2]<- args[1]
      attr(t, "names")
      
      ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                             col = list(condition = t))
      
      p<- Heatmap(Heat,km =2, name = "vst_counts", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],  clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                  row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
      
      print(p)
      
      #If you want to visualize the group
      colnames(Heat)<- Metadata$CONDITION
      
      ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                             col = list(condition = t))
      p<- Heatmap(Heat,km =2, name = "vst_counts", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],  clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                  row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
      
      
      print(p)
      
      
      dev.off()
    }
    
    #ENRICHR ANALYSIS
    
    dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_pos_vs_",args[2],"/Pathway_analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    
    directory_path <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_pos_vs_",args[2],"/Pathway_analysis", sep ="")
    
    pdf(file=paste(directory_path,"/Pathways_DGE_EnrichR.pdf", sep=""), width = 20, height = 9)
    
    library(enrichR)
    dbs <- c("Reactome_2022", "GO_Biological_Process_2023", "CODE_and_ChEA_Consensus_TFs_from_ChIP-X")
    websiteLive <- getOption("enrichR.live")
    
    if(length(ALTI$Gene.name) != 0){
      if (websiteLive) {
        enriched <- enrichr(ALTI$Gene.name, dbs)}
      for(paths in 1:3){
        if (nrow(enriched[[paths]]) != 0){
          if (websiteLive) {
            dir.create(paste(directory_path,"/TABLES", sep=""))
            write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",names(enriched)[paths],"_upregulated.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
            xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_upregulated_genes","_",dbs[paths], sep=""))
            print(xx)
          }}}}
    
    if(length(BASSI$Gene.name) != 0){
      if (websiteLive) {
        enriched <- enrichr(BASSI$Gene.name, dbs)}
      for(paths in 1:3){
        if (nrow(enriched[[paths]]) != 0){
          if (websiteLive) {
            dir.create(paste(directory_path,"/TABLES", sep=""))
            write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",names(enriched)[paths],"_downregulated.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
            xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_downregulated_genes","_",dbs[paths], sep=""))
            print(xx)
          }}}}
    
    dev.off()
    
    
    
    
    
  }else{
    print("Positive samples lower than 3, the DGE analysis is not possible")
  }
  
  
  
  #### DGE ANALYSIS IFN NEG vs CONTROL ----
  
  library(vroom)
  library(dplyr)
  
  directory2 <- paste(directory, "/MOFA/INPUT/",Cell_type, "_",args[1], "_vs_", args[2], sep ="")
  setwd(directory2)
  
  
  DGE3 <-vroom(counts , delim = "\t", col_names = TRUE)
  Metadata <-vroom(Meta , delim = "\t", col_names = TRUE)
  
  rownames(DGE3) <- DGE3$Gene.name
  
  
  #colnames(DGE3) <- c(Metadata$OMICID, "gene")
  gene_name <- DGE3$Gene.name
  
  
  
  #Purity from MCP
  library(dplyr)
  
  setwd(paste(directory,"/Transcriptomics/INPUT",sep=""))
  
  #purity_filter= "YES"
  
  METADATA_FILT <- !is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])
  METADATA_FILT_INDEX <-grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))
  
  repetition = 0
  for (meta_filter in METADATA_FILT_INDEX){
    repetition <- repetition + 1 
    if (!is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])[repetition]){
      COLNAME<-as.character(COMMAND_ADVANCED[1,meta_filter])
      if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Numerical"){
        To_filter<-as.numeric(unlist(Metadata[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,1)
        characters_to_remove <- c(">", "<", "=", " ")
        value_threshold <- as.numeric(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "<" = function(a, b) a < b,
                                      ">" = function(a, b) a > b,
                                      "=" = function(a, b) a == b,
                                      ">=" = function(a, b) a >= b,
                                      "<=" = function(a, b) a <= b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
        DGE3 <- DGE3[,c(comparison_operator(To_filter, value_threshold),TRUE)]
        
      }else if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Factors"){
        To_filter<- as.character(unlist(Metadata[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,2)
        characters_to_remove <- c("!=", "==", " ")
        value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "==" = function(a, b) a == b,
                                      "!=" = function(a, b) a != b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
        DGE3 <- DGE3[,c(comparison_operator(To_filter, value_threshold),TRUE)]
      }
      
    }
  }
  
  DGE3 <- DGE3[,-ncol(DGE3)]
  
  #Metadata<- Metadata %>% arrange(condition)
  which(colnames(DGE3) %in% Metadata$ID )
  
  
  
  
  tmp = NULL
  for (i in 1:length(Metadata$ID)){
    x <-grep(Metadata$ID[i], colnames(DGE3))
    tmp<-append(tmp,x)
  }
  
  DGE3 <- DGE3[,tmp]
  
  colnames(DGE3) <- c(Metadata$ID)
  
  num <-Metadata$condition == "pos"
  
  DGE3 <- DGE3[,!num]
  Metadata <- Metadata[!num,]
  
  if (sum(Metadata$condition == "neg") > 2){
    
    if(NORMALIZATION=="YES"){
      
      library(limma)
      library(edgeR)
      
      rownames(DGE3) <- gene_name
      
      if ("AGE" %in% colnames(Metadata) & length(unique(Metadata$AGE)) > 1){
        
        
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          
          #=========================================
          
          CONDITION_<-as.factor(Metadata$CONDITION)
          AGE <- Metadata$AGE
          GENDER <- Metadata$GENDER
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_ + AGE + GENDER, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
          
        }else{
          CONDITION_<-as.factor(Metadata$CONDITION)
          AGE <- Metadata$AGE
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_ + AGE, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
          
        }
      }else{
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          
          CONDITION_<-as.factor(Metadata$CONDITION)
          GENDER <- Metadata$GENDER
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_ + GENDER, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
          
        }else{
          
          CONDITION_<-as.factor(Metadata$CONDITION)
          rownames(DGE3) <- gene_name
          
          design <- model.matrix(~ 0 + CONDITION_, data = DGE3)
          colnames(design)[1] <- args[[2]]
          colnames(design)[2] <- args[[1]]
          design
          
          contr.matrix <- makeContrasts(
            CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
            levels = design)
          contr.matrix
          
          colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
          rownames(contr.matrix)[1] <- args[2]
          rownames(contr.matrix)[1] <- args[1]
          
          vfit <- lmFit(DGE3, design)
          vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
          efit <- eBayes(vfit,trend=TRUE)
          plotSA(efit, main="Final model: Mean-variance trend")
          dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
        }}
      
      
      #create directory
      dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_pos_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      
      directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_pos_vs_",args[2], sep ="")
      
      setwd(directory2)
      
      dds<-dds[,c(2,1,6,3,4,5)]
      colnames(dds) <- c("AveExpr", "log2FoldChange", "B" , "t", "pvalue", "padj"  )
      
      write.table(dds,
                  file= paste(args[1], "_neg_", args[2],"_",Cell_type ,".tsv",sep =""), sep="\t")
      
      Normalized_heatmap <- DGE3
      
      
    }else{
      
      #colnames(DGE3) <- c(Metadata$OMICID, "gene")
      rownames(DGE3) <- gene_name
      #DGE3$GENE <- gene_name
      #DGE3$GENE
      
      library("DESeq2")
      if ("AGE" %in% colnames(Metadata) & length(unique(Metadata$AGE)) > 1){
        
        Metadata$AGE[Metadata$AGE <= 35] <- "less30"
        Metadata$AGE[Metadata$AGE <= 50 & Metadata$AGE > 35] <- "35_50"
        Metadata$AGE[Metadata$AGE > 50] <- "more50"
        
        
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          
          #=========================================
          
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ GENDER + AGE + condition)
          
        }else{
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ AGE + condition)}
      }else{
        if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ GENDER + condition)
          
        }else{
          dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                        colData = Metadata,
                                        design = ~ condition)}}
      
      
      
      
      dds$condition <- relevel(dds$condition, ref = "healthy")
      
      #dds$condition <- relevel(dds$condition, ref = "neg")
      
      #crea l'oggetto di classe DESeqDataSet, quello utilizzato dal pacchetto
      ddsHTSeq2 <- dds
      
      
      dds <- DESeq(dds)
      
      res <- results(dds)
      res
      #genera l'espressione differeziale sulla base del log2 fold change
      #aggiungendo il p-value e il p-value adjusted
      
      resultsNames(dds)
      #Number of comparison available
      
      
      resOrdered <- res[order(res$padj),]
      #ordina i geni in base alla significatività del p_value adjusted
      
      summary(res)
      
      sum(res$padj < 0.1, na.rm=TRUE)
      #quanti p_value adjusted sono meno di 0.1?
      sum(res$padj > 0.1, na.rm=TRUE)
      #e maggiori di 0.1?
      
      #informazioni su variabili e test utilizzati
      mcols(res)$description
      
      
      #create directory
      dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_neg_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      
      directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_neg_vs_",args[2], sep ="")
      
      setwd(directory2)
      
      write.table(as.data.frame(resOrdered),
                  file= paste(args[1], "_neg_", args[2],"_",Cell_type ,".tsv",sep =""), sep = "\t")
      
      Normalized_heatmap <- as.data.frame(vst(assay(dds))) #Normalization VST
      
    }
    
    #permette di scrivere un file csv contenente tutti i risultati (in questo caso con pvalue adj ordinato)
    
    setwd(directory2)
    DGE <- read.table(file = paste(args[1], "_neg_", args[2],"_", Cell_type,".tsv", sep =""), sep="\t", header = TRUE)
    DGE$X <- rownames(DGE)
    DGE <- DGE %>% select(X, everything())
    
    
    dir.create(path = paste(directory2,"/Subpopulation_input/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    DGE3_modified <- DGE3
    DGE3_modified$ID <- rownames(DGE3_modified)
    DGE3_modified <- DGE3_modified %>% select(ID, everything())
    write.table(x=DGE3_modified   , file= paste(directory2,"/Subpopulation_input/",args[1], "_neg_and_", args[2],"_", Cell_type,"_counts.tsv", sep ="")  ,sep= "\t", row.names = FALSE, col.names = T,  quote = FALSE)
    write.table(x=Metadata   , file= paste(directory2,"/Subpopulation_input/",args[1], "_neg_and_", args[2],"_", Cell_type,"_metadata.tsv", sep ="")  ,sep= "\t", row.names = FALSE, col.names = T,  quote = FALSE)
    DGE3_modified = NULL
    gc()
    
    
    if(GENE_ANNOTATION == "GENE_NAME"){
      print(" ")
      colnames(DGE)[1] <- "Gene.name"
    }else{
      setwd(paste(directory,"/MOFA/x_BiomiX_DATABASE",sep=""))
      Mart<-read.table("mart_export_37.txt", sep = ",", header=TRUE)
      DGE <- merge(DGE, Mart, by.x = "X", by.y = "Gene.stable.ID")}
    
    #colnames(DGE)[7] <- "Gene.name"
    
    setwd(directory2)
    write.table(as.data.frame(DGE),
                file= paste(args[1], "_neg_",args[2],"_",Cell_type,".tsv", sep =""),sep = "\t", row.names = FALSE)
    
    
    setwd(directory2)
    
    if(NORMALIZATION=="NO"){
      tab<- as.data.frame(counts(dds,normalized=TRUE))
    }else{tab <- DGE3}
    tab$X<- rownames(tab)
    
    
    if(GENE_ANNOTATION == "GENE_NAME"){
      tab <- tab %>% select(X, everything())
    }else{
      tab <- merge(tab, Mart, by.x = "X", by.y = "Gene.stable.ID")
      #tab<-tab[,c(ncol(tab),3:ncol(tab)-1)]
      tab<-tab[,-ncol(tab)]}
    
    
    colnames(tab)[1] <- "NAME"
    tab$Description <- "NA"
    
    tab<-tab[,c(1,ncol(tab),3:ncol(tab)-1)]
    
    
    file.remove(paste(args[1], "neg_",args[2],"_",Cell_type, "_RAW.gct", sep =""))
    tab<- as.data.frame(tab)
    line="#1.2"
    write(line,file=paste(args[1], "neg_",args[2],"_",Cell_type, "_RAW.gct", sep =""),append=TRUE)
    line=c((nrow(tab)),(ncol(tab)-2))
    write(line,file=paste(args[1], "neg_",args[2],"_",Cell_type, "_RAW.gct", sep =""),append=TRUE,sep = "\t")
    write.table(tab, append =TRUE,
                file= paste(args[1], "neg_",args[2],"_",Cell_type, "_RAW.gct", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
    write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE),],
                file= paste(args[1], "neg_",args[2],"_",Cell_type, "_RAW_meta.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
    
    dis<-Metadata[order(Metadata$CONDITION),]
    
    write.table(dis[,c("ID","CONDITION")],
                file= paste(args[1], "neg_",args[2],"_",Cell_type, "_RAW_meta_groups.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
    #Ho preso i nomi dei diversi geni attraverso bio.mart da questo ho scaricato
    #la tabella e l'ho aggiunta a quella DGE.
    #
    #
    #
    ALTI <- subset(DGE, padj < padju)
    ALTI <-subset(ALTI, log2FoldChange > log2FC)
    ALTI <- ALTI[order(ALTI$padj), ]
    #
    BASSI <- subset(DGE, padj < padju)
    BASSI <- subset(BASSI, log2FoldChange < -log2FC)
    BASSI <- BASSI[order(BASSI$padj), ]
    NO <- subset(DGE,padj > padju)
    
    
    setwd(directory2)
    
    # write.csv(arrange(ALTI,padj),
    #           file="HC_vs_SLEpos_upregulated")
    #
    #
    # write.csv(arrange(BASSI,padj),
    #           file="HC_vs_SLEpos_downregulated")
    
    
    library(ggplot2)
    library(ggrepel)
    setwd(directory2)
    pdf(file=paste("plot_DGE_", args[1],"_neg_",args[2],"_", Cell_type,".pdf"))
    
    p <-ggplot() +
      ggtitle( paste("DGE", Cell_type,"_", args[1], "_NEG_vs_", args[2],sep="")) + theme(
        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold")) +
      geom_point(data = ALTI, aes(x = log2FoldChange, y = -log10(padj)), color = "red") +
      geom_text_repel(data = ALTI[c(1:25), ], aes(x = log2FoldChange , y = -log10(padj), label=Gene.name),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1)+
      geom_point(data = BASSI, aes(x = log2FoldChange, y = -log10(padj)), color = "blue") +
      geom_text_repel(data = BASSI[1:25, ], aes(x = log2FoldChange , y = -log10(padj), label=Gene.name),hjust=0, vjust=0,size=3, force = 20, force_pull = 1)+
      geom_point(data = NO, aes(x = log2FoldChange, y = -log10(padj)), color = "black")
    
    p2 <- p + geom_vline(xintercept=c(-log2FC,log2FC), col="red") +
      geom_hline(yintercept=-log10(padju), col="red")
    
    
    print(p2)
    dev.off()
    
    
    setwd(directory2)
    
    write.table(x=ALTI$Gene.name   , file= paste(Cell_type,"_", args[1],"neg_",args[2],"_","UP_ENRICHR.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x=BASSI$Gene.name  , file= paste(Cell_type,"_",args[1],"neg_",args[2],"_","DOWN_ENRICHR.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x=ALTI   , file= paste(Cell_type,"_", args[1],"neg_",args[2],"_","GENESUP.tsv",sep="") ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    write.table(x=BASSI  , file= paste(Cell_type,"_", args[1],"neg_",args[2],"_","GENESDOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    
    
    if ((nrow(ALTI) + nrow(BASSI)) > 2){                
      
      #ADD HEATMAP CREATION HERE
      if (nrow(ALTI) > n_genes_heat){
        ALTI <- ALTI[1:n_genes_heat,]}
      if (nrow(BASSI) > n_genes_heat){
        BASSI <- BASSI[1:n_genes_heat,]}
      
      #HEATMAP INPUT ARE NORMALIZED
      
      Heat<-Normalized_heatmap[rownames(Normalized_heatmap) %in% c(ALTI$X, BASSI$X),]
      to_add<-rownames(Normalized_heatmap)[rownames(Normalized_heatmap) %in% c(ALTI$X, BASSI$X)]
      rownames(Heat) <-to_add
      
      if (GENE_ANNOTATION == "ENSEMBL"){
        
        Heat$GENE <-rownames(Heat)
        Heat <- merge(Heat, Mart, by.x = "GENE", by.y = "Gene.stable.ID")
        Heat <- Heat[!duplicated( Heat[,ncol(Heat)]),]
        rownames(Heat) <- Heat[,ncol(Heat)]
        Heat <- Heat[,c(-1,-ncol(Heat))]
      }
      
      library(ComplexHeatmap)
      setwd(directory2)
      pdf(file= paste("Heatmap_top_genes", Cell_type, args[2],"vs", args[1],"neg.pdf",sep="_"))
      #dev.new()
      
      library(circlize)
      col_fun = colorRamp2(c(min(Heat) - 3, mean(colMeans(Heat)), max(Heat) - 3), c("blue", "black", "yellow"))
      col_fun(seq(-3, 3))
      
      t<-c("CTRL" = "blue", "SLE" = "red")
      attr(t, "names")[1]<- args[2]
      attr(t, "names")[2]<- args[1]
      attr(t, "names")
      
      ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                             col = list(condition = t))
      p<- Heatmap(Heat,km =2, name = "vst_counts", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],  clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                  row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
      
      
      print(p)
      
      
      #If you want to visualize the group
      colnames(Heat)<- Metadata$CONDITION
      
      ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                             col = list(condition = t))
      p<- Heatmap(Heat,km =2, name = "vst_counts", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],  clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                  row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
      
      print(p)
      
      dev.off()
      
    }
    
    #ENRICHR ANALYSIS
    
    dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_neg_vs_",args[2],"/Pathway_analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    
    directory_path <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_neg_vs_",args[2],"/Pathway_analysis", sep ="")
    
    pdf(file=paste(directory_path,"/Pathways_DGE_EnrichR.pdf", sep=""), width = 20, height = 9)
    
    library(enrichR)
    dbs <- c("Reactome_2022", "GO_Biological_Process_2023", "CODE_and_ChEA_Consensus_TFs_from_ChIP-X")
    websiteLive <- getOption("enrichR.live")
    
    if(length(ALTI$Gene.name) != 0){
      if (websiteLive) {
        enriched <- enrichr(ALTI$Gene.name, dbs)}
      for(paths in 1:3){
        if (nrow(enriched[[paths]]) != 0){
          if (websiteLive) {
            dir.create(paste(directory_path,"/TABLES", sep=""))
            write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",names(enriched)[paths],"_upregulated.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
            xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_upregulated_genes","_",dbs[paths], sep=""))
            print(xx)
          }}}}
    
    if(length(BASSI$Gene.name) != 0){
      if (websiteLive) {
        enriched <- enrichr(BASSI$Gene.name, dbs)}
      for(paths in 1:3){
        if (nrow(enriched[[paths]]) != 0){
          if (websiteLive) {
            dir.create(paste(directory_path,"/TABLES", sep=""))
            write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",names(enriched)[paths],"_downregulated.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
            xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_downregulated_genes","_",dbs[paths], sep=""))
            print(xx)
          }
        }
      }
    }
    
    
    dev.off()
    
    
    
  }else{
    print("Negative samples lower than 3, the DGE analysis is not possible")
  }
  
}else{
  print("NO DGE on subpopulations")
}

#### DGE DISEASE VS CONTROL ----


directory2 <- paste(directory, "/MOFA/INPUT/",Cell_type, "_",args[1], "_vs_", args[2], sep ="")
setwd(directory2)


Metadata <- Metadata_all
DGE3 <-vroom(counts , delim = "\t", col_names = TRUE)
gene_name <- DGE3$Gene.name



#Purity from MCP
library(dplyr)

setwd(paste(directory,"/Transcriptomics/INPUT",sep = ""))


METADATA_FILT <- !is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])
METADATA_FILT_INDEX <-grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))

repetition = 0
for (meta_filter in METADATA_FILT_INDEX){
  repetition <- repetition + 1 
  if (!is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])[repetition]){
    COLNAME<-as.character(COMMAND_ADVANCED[1,meta_filter])
    if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Numerical"){
      To_filter<-as.numeric(unlist(Metadata[,COLNAME]))
      simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,1)
      characters_to_remove <- c(">", "<", "=", " ")
      value_threshold <- as.numeric(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
      
      comparison_operator <- switch(simbol,
                                    "<" = function(a, b) a < b,
                                    ">" = function(a, b) a > b,
                                    "=" = function(a, b) a == b,
                                    ">=" = function(a, b) a >= b,
                                    "<=" = function(a, b) a <= b,
                                    NA)
      
      Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
      DGE3 <- DGE3[,c(comparison_operator(To_filter, value_threshold),TRUE)]
      
    }else if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Factors"){
      To_filter<- as.character(unlist(Metadata[,COLNAME]))
      simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,2)
      characters_to_remove <- c("!=", "==", " ")
      value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
      
      comparison_operator <- switch(simbol,
                                    "==" = function(a, b) a == b,
                                    "!=" = function(a, b) a != b,
                                    NA)
      
      Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
      DGE3 <- DGE3[,c(comparison_operator(To_filter, value_threshold),TRUE)]
    }
    
  }
}



DGE3 <- DGE3[,-ncol(DGE3)]

#Metadata<- Metadata %>% arrange(CONDITION)
which(colnames(DGE3) %in% Metadata$ID )

tmp = NULL
for (i in 1:length(Metadata$ID)){
  x <-grep(Metadata$ID[i], colnames(DGE3))
  tmp<-append(tmp,x)
}

DGE3 <- DGE3[,tmp]
rownames(DGE3) <- gene_name

if(NORMALIZATION=="YES"){
  
  library(limma)
  library(edgeR)
  
  rownames(DGE3) <- gene_name
  
  if ("AGE" %in% colnames(Metadata) & length(unique(Metadata$AGE)) > 1){
    
    
    if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
      
      #=========================================
      
      CONDITION_<-as.factor(Metadata$CONDITION)
      AGE <- Metadata$AGE
      GENDER <- Metadata$GENDER
      rownames(DGE3) <- gene_name
      
      design <- model.matrix(~ 0 + CONDITION_ + AGE + GENDER, data = DGE3)
      colnames(design)[1] <- args[[2]]
      colnames(design)[2] <- args[[1]]
      design
      
      contr.matrix <- makeContrasts(
        CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
        levels = design)
      contr.matrix
      
      colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
      rownames(contr.matrix)[1] <- args[2]
      rownames(contr.matrix)[1] <- args[1]
      
      vfit <- lmFit(DGE3, design)
      vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
      efit <- eBayes(vfit,trend=TRUE)
      plotSA(efit, main="Final model: Mean-variance trend")
      dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
      
    }else{
      CONDITION_<-as.factor(Metadata$CONDITION)
      AGE <- Metadata$AGE
      rownames(DGE3) <- gene_name
      
      design <- model.matrix(~ 0 + CONDITION_ + AGE, data = DGE3)
      colnames(design)[1] <- args[[2]]
      colnames(design)[2] <- args[[1]]
      design
      
      contr.matrix <- makeContrasts(
        CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
        levels = design)
      contr.matrix
      
      colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
      rownames(contr.matrix)[1] <- args[2]
      rownames(contr.matrix)[1] <- args[1]
      
      vfit <- lmFit(DGE3, design)
      vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
      efit <- eBayes(vfit,trend=TRUE)
      plotSA(efit, main="Final model: Mean-variance trend")
      dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
      
    }
  }else{
    if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
      
      CONDITION_<-as.factor(Metadata$CONDITION)
      GENDER <- Metadata$GENDER
      rownames(DGE3) <- gene_name
      
      design <- model.matrix(~ 0 + CONDITION_ + GENDER, data = DGE3)
      colnames(design)[1] <- args[[2]]
      colnames(design)[2] <- args[[1]]
      design
      
      contr.matrix <- makeContrasts(
        CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
        levels = design)
      contr.matrix
      
      colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
      rownames(contr.matrix)[1] <- args[2]
      rownames(contr.matrix)[1] <- args[1]
      
      vfit <- lmFit(DGE3, design)
      vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
      efit <- eBayes(vfit,trend=TRUE)
      plotSA(efit, main="Final model: Mean-variance trend")
      dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
      
    }else{
      
      CONDITION_<-as.factor(Metadata$CONDITION)
      rownames(DGE3) <- gene_name
      
      design <- model.matrix(~ 0 + CONDITION_, data = DGE3)
      colnames(design)[1] <- args[[2]]
      colnames(design)[2] <- args[[1]]
      design
      
      contr.matrix <- makeContrasts(
        CTRLvsCONDITION = paste(colnames(design)[1],"-",colnames(design)[2]),
        levels = design)
      contr.matrix
      
      colnames(contr.matrix) <- paste(args[2],"_vs_",args[1], sep="")
      rownames(contr.matrix)[1] <- args[2]
      rownames(contr.matrix)[1] <- args[1]
      
      vfit <- lmFit(DGE3, design)
      vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
      efit <- eBayes(vfit,trend=TRUE)
      plotSA(efit, main="Final model: Mean-variance trend")
      dds <- as.data.frame(topTreat(efit, coef=1, n=Inf, adjust="BH"))
    }}
  
  
  #create directory
  dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
  dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_vs_",args[2], sep ="")
  
  setwd(directory2)
  
  dds<-dds[,c(2,1,6,3,4,5)]
  colnames(dds) <- c("AveExpr", "log2FoldChange", "B" , "t", "pvalue", "padj"  )
  
  write.table(dds,
              file= paste(args[1], "_", args[2],"_",Cell_type ,".tsv",sep =""), sep="\t")
  
  Normalized_heatmap <- DGE3
  
  
}else{
  
  library("DESeq2")
  if ("AGE" %in% colnames(Metadata) & length(unique(Metadata$AGE)) > 1){
    
    Metadata$AGE[Metadata$AGE <= 35] <- "less30"
    Metadata$AGE[Metadata$AGE <= 50 & Metadata$AGE > 35] <- "35_50"
    Metadata$AGE[Metadata$AGE > 50] <- "more50"
    
    
    if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
      
      #=========================================
      
      dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                    colData = Metadata,
                                    design = ~ GENDER + AGE + CONDITION)
      
    }else{
      dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                    colData = Metadata,
                                    design = ~ AGE + CONDITION)}
  }else{
    if ("GENDER" %in% colnames(Metadata) & length(unique(Metadata$GENDER)) > 1) {
      dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                    colData = Metadata,
                                    design = ~ GENDER + CONDITION)
      
    }else{
      dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                    colData = Metadata,
                                    design = ~ CONDITION)}}
  
  
  
  
  
  
  y<-unlist(args[2])
  dds$CONDITION<- as.factor(dds$CONDITION)
  dds$CONDITION <- relevel(dds$CONDITION, ref = y)
  
  #dds$condition <- relevel(dds$condition, ref = "neg")
  
  #crea l'oggetto di classe DESeqDataSet, quello utilizzato dal pacchetto
  ddsHTSeq2 <- dds
  
  
  dds <- DESeq(dds)
  
  
  res <- results(dds)
  res
  #genera l'espressione differeziale sulla base del log2 fold change
  #aggiungendo il p-value e il p-value adjusted
  
  resultsNames(dds)
  #Number of comparison available
  
  
  resOrdered <- res[order(res$padj),]
  #ordina i geni in base alla significatività del p_value adjusted
  
  summary(res)
  
  sum(res$padj < 0.1, na.rm=TRUE)
  #quanti p_value adjusted sono meno di 0.1?
  sum(res$padj > 0.1, na.rm=TRUE)
  #e maggiori di 0.1?
  
  #informazioni su variabili e test utilizzati
  mcols(res)$description
  
  #create directory
  dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_vs_",args[2], sep ="")
  
  setwd(directory2)
  
  write.table(as.data.frame(resOrdered),
              file= paste(args[1], "_",args[2],"_",Cell_type ,".tsv",sep =""), sep = "\t") #
  
  Normalized_heatmap <- as.data.frame(vst(assay(dds))) #Normalization VST
  
}

#permette di scrivere un file csv contenente tutti i risultati (in questo caso con pvalue adj ordinato)

DGE <- read.table(file = paste(args[1], "_",args[2],"_",Cell_type ,".tsv",sep ="") , sep="\t", header = TRUE)
DGE$X <- rownames(DGE)
DGE <- DGE %>% select(X, everything())


if(GENE_ANNOTATION == "GENE_NAME"){
  print(" ")
  colnames(DGE)[1] <- "Gene.name"
}else{
  setwd(paste(directory,"/MOFA/x_BiomiX_DATABASE",sep=""))
  Mart<-read.table("mart_export_37.txt", sep = ",", header=TRUE)
  DGE <- merge(DGE, Mart, by.x = "X", by.y = "Gene.stable.ID")}
#colnames(DGE)[7] <- "Gene.name"


setwd(directory2)
write.table(as.data.frame(DGE),
            file= paste(args[1], "_",args[2],"_",Cell_type,".tsv", sep =""),sep = "\t", row.names = FALSE)

if(NORMALIZATION=="NO"){
  tab<- as.data.frame(counts(dds,normalized=TRUE))
}else{tab <- DGE3}
tab$X<- rownames(tab)

if(GENE_ANNOTATION == "GENE_NAME"){
  tab <- tab %>% select(X, everything())
}else{
  tab <- merge(tab, Mart, by.x = "X", by.y = "Gene.stable.ID")
  #tab<-tab[,c(ncol(tab),3:ncol(tab)-1)]
  tab<-tab[,-ncol(tab)]}


colnames(tab)[1] <- "NAME"
tab$Description <- "NA"

tab<-tab[,c(1,ncol(tab),3:ncol(tab)-1)]


setwd(directory2)


tab<- as.data.frame(tab)
file.remove(paste(args[1], "_",args[2],"_",Cell_type, "_RAW.gct", sep =""))
line="#1.2"
write(line,file=paste(args[1], "_",args[2],"_",Cell_type, "_RAW.gct", sep =""),append=TRUE)
line=c((nrow(tab)),(ncol(tab)-2))
write(line,file=paste(args[1], "_",args[2],"_",Cell_type, "_RAW.gct", sep =""),append=TRUE,sep = "\t")
write.table(tab, append =TRUE,
            file= paste(args[1], "_",args[2],"_",Cell_type, "_RAW.gct", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE),],
            file= paste(args[1], "_",args[2],"_",Cell_type, "_RAW_meta.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)

dis<-Metadata[order(Metadata$CONDITION),]

write.table(dis[,c("ID","CONDITION")],
            file= paste(args[1], "_", args[2],"_",Cell_type, "_RAW_meta_groups.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)#Ho preso i nomi dei diversi geni attraverso bio.mart da questo ho scaricato
#la tabella e l'ho aggiunta a quella DGE.
#
#
#
ALTI <- subset(DGE, padj < padju)
ALTI <-subset(ALTI, log2FoldChange > log2FC)
ALTI <- ALTI[order(ALTI$padj), ]
#ALTI <- ALTI %>% filter(baseMean > 3)
#
BASSI <- subset(DGE, padj < padju)
BASSI <- subset(BASSI, log2FoldChange < -log2FC)
BASSI <- BASSI[order(BASSI$padj), ]
#BASSI <- BASSI %>% filter(baseMean > 3)
NO <- subset(DGE,log2FoldChange < log2FC & log2FoldChange > -log2FC)

library(ggplot2)
library(ggrepel)
setwd(directory2)
pdf(file=paste("plot_DGE_", args[1],"_pos_",args[2],"_", Cell_type,".pdf"))

p <-ggplot() +
  ggtitle( paste("DGE", Cell_type,"_", args[1], "_vs_", args[2],sep="")) + theme(
    plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")) +
  geom_point(data = ALTI, aes(x = log2FoldChange, y = -log10(padj)), color = "red") +
  geom_text_repel(data = ALTI[c(1:25), ], aes(x = log2FoldChange , y = -log10(padj), label=Gene.name),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1,max.overlaps = getOption("ggrepel.max.overlaps", default = 25))+
  geom_point(data = BASSI, aes(x = log2FoldChange, y = -log10(padj)), color = "blue") +
  geom_text_repel(data = BASSI[1:25, ], aes(x = log2FoldChange , y = -log10(padj), label=Gene.name),hjust=0, vjust=0,size=3, force = 20, force_pull = 1,max.overlaps = getOption("ggrepel.max.overlaps", default = 25))+
  geom_point(data = NO, aes(x = log2FoldChange, y = -log10(padj)), color = "black")

p2 <- p + geom_vline(xintercept=c(-log2FC,log2FC), col="red") +
  geom_hline(yintercept=-log10(padju), col="red")


print(p2)
dev.off()





DGE3$Gene.name <- gene_name

setwd(directory2)

write.table(x=ALTI$Gene.name   , file= paste(Cell_type,"_",args[1],"_",args[2],"_UP_ENRICHR.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
write.table(x=BASSI$Gene.name  , file= paste(Cell_type,"_",args[1],"_",args[2],"_DOWN_ENRICHR.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
write.table(x=ALTI   , file= paste(Cell_type,"_",args[1],"_",args[2],"_GENESUP.tsv",sep="") ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(x=BASSI  , file= paste(Cell_type,"_",args[1],"_",args[2],"_GENESDOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)


if ((nrow(ALTI) + nrow(BASSI)) > 2){
  
  #ADD HEATMAP CREATION HERE
  if (nrow(ALTI) > n_genes_heat){
    ALTI <- ALTI[1:n_genes_heat,]}
  if (nrow(BASSI) > n_genes_heat){
    BASSI <- BASSI[1:n_genes_heat,]}
  
  #HEATMAP INPUT ARE NORMALIZED
  
  Heat<-Normalized_heatmap[rownames(Normalized_heatmap) %in% c(ALTI$X, BASSI$X),]
  to_add<-rownames(Normalized_heatmap)[rownames(Normalized_heatmap) %in% c(ALTI$X, BASSI$X)]
  rownames(Heat) <-to_add
  
  if (GENE_ANNOTATION == "ENSEMBL"){
    
    Heat$GENE <-rownames(Heat)
    Heat <- merge(Heat, Mart, by.x = "GENE", by.y = "Gene.stable.ID")
    Heat <- Heat[!duplicated( Heat[,ncol(Heat)]),]
    rownames(Heat) <- Heat[,ncol(Heat)]
    Heat <- Heat[,c(-1,-ncol(Heat))]
  }
  
  library(ComplexHeatmap)
  setwd(directory2)
  pdf(file= paste("Heatmap_top_genes_", Cell_type,"_", args[2],"_vs_", args[1],".pdf",sep=""))
  #dev.new()
  
  library(circlize)
  col_fun = colorRamp2(c(min(Heat) - 3, mean(colMeans(Heat)), max(Heat) - 3), c("blue", "black", "yellow"))
  col_fun(seq(-3, 3))
  
  t<-c("CTRL" = "blue", "SLE" = "red")
  attr(t, "names")[1]<- args[2]
  attr(t, "names")[2]<- args[1]
  attr(t, "names")
  
  ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                         col = list(condition = t))
  p<- Heatmap(Heat,km =2, name = "vst_counts", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],  clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
              row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
  
  print(p)
  
  
  #If you want to visualize the group
  colnames(Heat)<- Metadata$CONDITION
  
  ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                         col = list(condition = t))
  p<- Heatmap(Heat,km =2, name = "vst_counts", col = col_fun, clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],  clustering_method_rows= "complete", clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2], row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
              row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
  
  print(p)
  
  dev.off()
  
}

#ENRICHR ANALYSIS

dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_vs_",args[2],"/Pathway_analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory_path <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_vs_",args[2],"/Pathway_analysis", sep ="")

pdf(file=paste(directory_path,"/Pathways_DGE_EnrichR.pdf", sep=""), width = 20, height = 9)

library(enrichR)
dbs <- c("Reactome_2022", "GO_Biological_Process_2023", "CODE_and_ChEA_Consensus_TFs_from_ChIP-X")
websiteLive <- getOption("enrichR.live")

if(length(ALTI$Gene.name) != 0){
  if (websiteLive) {
    enriched <- enrichr(ALTI$Gene.name, dbs)}
  for(paths in 1:3){
    if (nrow(enriched[[paths]]) != 0){
      if (websiteLive) {
        dir.create(paste(directory_path,"/TABLES", sep=""))
        write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",names(enriched)[paths],"_upregulated.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_upregulated_genes","_",dbs[paths], sep=""))
        print(xx)
      }}}}

if(length(BASSI$Gene.name) != 0){
  if (websiteLive) {
    enriched <- enrichr(BASSI$Gene.name, dbs)}
  for(paths in 1:3){
    if (nrow(enriched[[paths]]) != 0){
      if (websiteLive) {
        dir.create(paste(directory_path,"/TABLES", sep=""))
        write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",names(enriched)[paths],"_downregulated.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_downregulated_genes","_",dbs[paths], sep=""))
        print(xx)
      }}}}

dev.off()

#### GENERATION MATRIX FOR DATA INTEGRATION ----

if (NORMALIZATION=="NO"){
  print("PATH NO NORMALIZED 1")
  
  directory2 <- paste(directory,"/MOFA/INPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="")
  setwd(directory2)
  
  normalizzati <-vroom(file = paste(Cell_type, "_",args[1],"_vs_", args[2],"_normalized_vst.tsv", sep="") , delim = "\t", col_names = TRUE)
  normalizzati_size <-vroom(file = paste(Cell_type, "_",args[1],"_vs_", args[2],"_normalized.tsv", sep="") , delim = "\t", col_names = TRUE)
  
  genes <- normalizzati$Gene.name
  normalizzati <- normalizzati[,-ncol(normalizzati)]
  rownames(normalizzati) <- genes
  
  genes <- normalizzati_size$Gene.name
  normalizzati_size <- normalizzati_size[,-ncol(normalizzati_size)]
  rownames(normalizzati_size) <- genes
  
  
  x <- normalizzati_size > 0
  righe_da_eliminare=NULL
  for (i in 1:nrow(normalizzati_size)){
    if ((sum(x[i,]) <= (ncol(normalizzati_size)/2))){
      righe_da_eliminare<-append(righe_da_eliminare,i)
    }
  }
  
  
  print("numero di geni prima del filtraggio")
  print(nrow(normalizzati_size))
  
  if (length(righe_da_eliminare) > 0){
    normalizzati_size<-normalizzati_size[-righe_da_eliminare,]
    normalizzati<-normalizzati[-righe_da_eliminare,]
    gene <- genes[-righe_da_eliminare]}
  
  print("numero di geni dopo il filtraggio")
  print(nrow(normalizzati_size))
  print(nrow(normalizzati))
  
}else{
  print("PATH NORMALIZED 1")
  directory2 <- paste(directory,"/MOFA/INPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="")
  setwd(directory2)
  
  normalizzati <-vroom(file = paste(Cell_type, "_",args[1],"_vs_", args[2],"_normalized_vst.tsv", sep="") , delim = "\t", col_names = TRUE)
  
  genes <- normalizzati$Gene.name
  normalizzati <- normalizzati[,-ncol(normalizzati)]
  rownames(normalizzati) <- genes
  
  x <- normalizzati > 0
  righe_da_eliminare=NULL
  for (i in 1:nrow(normalizzati)){
    if ((sum(x[i,]) <= (ncol(normalizzati)/2))){
      righe_da_eliminare<-append(righe_da_eliminare,i)
    }
  }
  
  print("numero di geni prima del filtraggio")
  print(nrow(normalizzati))
  
  if (length(righe_da_eliminare) > 0){
    normalizzati<-normalizzati[-righe_da_eliminare,]
    gene <- genes[-righe_da_eliminare]}else{gene <- genes}
  
  print("numero di geni dopo il filtraggio")
  print(nrow(normalizzati))
  
}



#=====================GENE FILTER ON VARIANCE==============================

if (NORMALIZATION=="NO"){
  print("PATH NO NORMALIZED 2")
  
  varianze = NULL
  #Vari <- read.table("tab_finaleDESeq2", header=TRUE, sep = ",")
  Vari <- as.matrix(normalizzati_size)
  #rownames(Vari) <- Vari[,8]
  Vari <- as.data.frame(Vari)
  #Vari<- Vari[,c(-1,-ncol(Vari))]
  genes <- Vari$Gene.name
  Vari<- Vari[,-ncol(Vari)]
  #generazione e filtraggio per varianza
  for (i in 1:nrow(Vari)){
    varianze <- append(varianze, var(as.numeric(Vari[i,])))
  }
  
  normalizzati<- as.data.frame(normalizzati)
  normalizzati$Gene.name <- gene
  normalizzati$variance <- log10(varianze)
  
  normalizzati_size<- as.data.frame(normalizzati_size)
  normalizzati_size$Gene.name <- gene
  normalizzati_size$variance <- log10(varianze)
  
  if(GENE_ANNOTATION == "GENE_NAME"){
    MART<-read.table(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt", sep=""), sep = ",", header=TRUE)
    colnames(MART)[2] <-"Gene"
    normalizzati <- merge(normalizzati, MART, by.x = "Gene.name", by.y = "Gene")
    normalizzati <- normalizzati[,-1]
    colnames(normalizzati)[ncol(normalizzati)] <- "Gene.name"
    normalizzati <- normalizzati[, c(1:(ncol(normalizzati) - 2), ncol(normalizzati), ncol(normalizzati) -1 )]
  }
  
  #OVERWRITE OR REPLACE A FINAL NORMALIZED FILE WITH VARIANCE IN IT.
write.table(x=normalizzati  , file= paste(Cell_type,args[1], "vs", args[2], "normalized_vst_variance.tsv",sep = "_")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
  
}else{
  print("PATH YES NORMALIZED 2")
  varianze = NULL
  #Vari <- read.table("tab_finaleDESeq2", header=TRUE, sep = ",")
  Vari <- as.matrix(normalizzati)
  #rownames(Vari) <- Vari[,8]
  Vari <- as.data.frame(Vari)
  #Vari<- Vari[,c(-1,-ncol(Vari))]
  
  genes <- rownames(Vari)
  
  #generazione e filtraggio per varianza
  for (i in 1:nrow(Vari)){
    varianze <- append(varianze, var(as.numeric(Vari[i,])))
  }
  
  normalizzati<- as.data.frame(normalizzati)
  normalizzati$Gene.name <- gene
  normalizzati$variance <- log10(varianze)
  
  
  if(GENE_ANNOTATION == "GENE_NAME"){
    MART<-read.table(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt", sep=""), sep = ",", header=TRUE)
    colnames(MART)[2] <-"Gene"
    normalizzati <- merge(normalizzati, MART, by.x = "Gene.name", by.y = "Gene")
    normalizzati <- normalizzati[,-1]
    colnames(normalizzati)[ncol(normalizzati)] <- "Gene.name"
    normalizzati <- normalizzati[, c(1:(ncol(normalizzati) - 2), ncol(normalizzati), ncol(normalizzati) -1 )]
  }
  
  write.table(x=normalizzati  , file= paste(Cell_type,args[1], "vs", args[2], "normalized_vst_variance.tsv",sep = "_")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
  normalizzati_size <- normalizzati
}


#PLOT VARIANCE
library(ggplot2)
pdf(file = paste("Gene_expression_Variance_", Cell_type,"_",args[1],"_vs_", args[2],".pdf" ,sep =""))

t<- ggplot(normalizzati_size, aes(x=variance)) +
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.7, bins=1000) + stat_bin(bins = 1000)
print(t)

p <-ggplot(normalizzati_size, aes(x=variance)) +
  geom_density(color="darkblue", fill="lightblue")
p <- p + geom_vline(aes(xintercept=2.5))
print(p)

dev.off()


# OUT<- normalizzati$varianze > 1.5
# normalizzati<-normalizzati[OUT,-ncol(normalizzati)]

directory2 <- paste(directory, "/MOFA/INPUT/",Cell_type, "_",args[1], "_vs_", args[2], sep ="")
setwd(directory2)

if(file.exists(paste("Gene_panel_subgroups_", Cell_type,"_", args[1],"_vs_", args[2], ".pdf",sep=""))){
  file.copy(paste("Gene_panel_subgroups_", Cell_type,"_", args[1],"_vs_", args[2], ".pdf",sep=""),paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_vs_",args[2], sep =""), overwrite=TRUE )
  
}

if(file.exists(paste("Gene_expression_Variance_", Cell_type,"_",args[1],"_vs_", args[2],".pdf" ,sep =""))){
  file.copy(paste("Gene_expression_Variance_", Cell_type,"_",args[1],"_vs_", args[2],".pdf" ,sep =""),paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_vs_",args[2], sep =""), overwrite=TRUE )
  
}

if(file.exists("Validation_MARKER_vs_GENE_PANEL.tsv")){
  file.copy("Validation_MARKER_vs_GENE_PANEL.tsv",paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/", args[1],"_vs_",args[2], sep =""), overwrite=TRUE )
  
}

gc()


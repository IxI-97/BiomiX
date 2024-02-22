library(vroom)
library(dplyr)

#### INPUT PROMPT ----

# args = commandArgs(trailingOnly=TRUE)
# 
# args = as.list(c("Wholeblood","SLE"))
# print(args)
# 
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"C1"
# args[2] <-"Control"
# args[3] <-"/home/henry/Desktop/BiomiX2.2"
# 
# directory <- args[3]
# iterations = 1
# i=5
# selection_samples = "NO"
# Cell_type = "Whole_methyl"
# DIR_METADATA <- readLines("/home/henry/Desktop/BiomiX2.2/directory.txt")

# Matrix <- vroom("WHOLEBLOOD_METHYLOMICS_SSC_vs_CTRL.tsv" , delim = "\t", col_names = TRUE)
# Metadata_total <- vroom("/home/henry/Desktop/BiomiX2.0/Metadata/Metadata_PRECISESADS.tsv" , delim = "\t", col_names = TRUE)


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

directory2 <- paste(directory,"/Methylomics/INPUT",sep="")
setwd(directory2)

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
LogFC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS[2])
array <- COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS[3]
n_genes_heat<-as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])




Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)
Matrix <- vroom(COMMAND$DIRECTORIES[i] , delim = "\t", col_names = TRUE)



if (selection_samples == "YES") {
        num <- grep(Cell_type, Metadata_total$ID_CELL_TYPE)
        Metadata_Bcell <- Metadata_total[num,]
        Metadata_total <- Metadata_total[num,]
        num <-which(colnames(Matrix) %in% Metadata_Bcell$ID_CELL_TYPE )
        Matrix <- as.matrix(Matrix[,num])
        print(Metadata_Bcell$ID_CELL_TYPE)
        num <-which(Metadata_Bcell$ID_CELL_TYPE  %in% colnames(Matrix) )
        Metadata_Bcell <- Metadata_Bcell[num,]
        #If multiple select only the first column
        Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
        
        Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
        Matrix<-Matrix[,order(colnames(Matrix))]
        
} else {
        print("No samples selection")
        Metadata_Bcell <- Metadata_total
        Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
        Matrix <- as.matrix(Matrix)
        Matrix<-Matrix[,order(colnames(Matrix))]
        
        Metadata_Bcell <- Metadata_total
        num <-which(colnames(Matrix) %in% Metadata_Bcell$ID)
        Matrix <- as.data.frame(Matrix)
        Identifier<-Matrix$ID
        Matrix <- as.matrix(Matrix[,num])
        print(Metadata_Bcell$ID)
        num <-which(Metadata_Bcell$ID  %in% colnames(Matrix) )
        Metadata_Bcell <- Metadata_Bcell[num,]
        #If multiple select only the first column
        Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
        
        Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
        Matrix<-Matrix[,order(colnames(Matrix))]
}

Metadata <- Metadata_Bcell




#FILTERING BASED ON METADATA

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
                        Matrix <- Matrix[,comparison_operator(To_filter, value_threshold)]
                        
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
                        Matrix <- Matrix[,comparison_operator(To_filter, value_threshold)]
                }
                
        }
}




Metadata_individual <- Metadata

num <- Metadata_individual$CONDITION == args[1]
#Choise disease + HC
numero <- Metadata_individual$CONDITION == args[2]
#Metadata_individual <- Metadata_individual[num | numero,]
Matrix <- as.data.frame(Matrix)
Matrix$ID <- Identifier

Matrix <- Matrix[, c(ncol(Matrix), 1:(ncol(Matrix)-1))]

vector<- num | numero
vector2 <- c(TRUE,vector)

Matrix <- Matrix[,vector2]

dim(Matrix)

#### METHYLATION FILTER ON VARIANCE + INTEGRATION ----


varianze = NULL

Vari <- as.matrix(Matrix[,-1])

varianze <-apply(Vari,1,var)
# for (i in 1:nrow(Vari)){
#         varianze <- append(varianze, var(as.numeric(Vari[i,])))
# }
gc()

Matrix2 <- Matrix
Matrix2$variance <- varianze

dir.create(path = paste(directory,"/MOFA/INPUT/", "Methylome_",Cell_type,"_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory,"/MOFA/INPUT/", "Methylome_",Cell_type, "_",args[1],"_vs_", args[2], sep ="") 

setwd(directory2)

#IN THE METADATA SELECT ONLY THE SAMPLES IN THE MATRIX
Metadata_individual<-Metadata_individual[vector,]

#SELECT JUST THE TOP 10000 FEATURES WITH HIGHER VIABILITY
Matrix2<- Matrix2 %>% arrange(desc(variance))


write.table(x=Matrix2[1:10000,]  , file= paste(Cell_type,"_matrix_MOFA.tsv",sep = "")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
write.table(x=Metadata_individual  , file= paste(Cell_type, "_metadata_MOFA.tsv", sep = "")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE,)


#### DMP ANALYSIS ----

dir.create(path = paste(directory,"/Methylomics/OUTPUT/", sep =""))
dir.create(path = paste(directory,"/Methylomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep =""))

directory2 <- paste(directory,"/Methylomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="")

setwd(directory2)

library(vroom)
library(dplyr)
library(ChAMP)

spot <- Matrix$ID
Matrix<- as.data.frame(Matrix[, -1])
rownames(Matrix) <- spot

Metadata_individual$CONDITION <- as.factor(Metadata_individual$CONDITION)

Matrix <- as.matrix(Matrix)
Matrix <-apply(Matrix,2, as.numeric)
Matrix<- as.data.frame(Matrix)
rownames(Matrix) <- spot

myDMP <- champ.DMP(beta = Matrix,pheno=Metadata_individual$CONDITION, adjPVal = 1, arraytype = array)

results <- myDMP[[1]]
results$CpG_island <-rownames(results)

directory2 <- paste(directory,"/Methylomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], sep ="")
setwd(directory2)

write.table(results , file= paste("DMP_", Cell_type, "_Methylome_",args[1],"_vs_", args[2],".tsv",sep = "")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)


ALTI <- subset(results, adj.P.Val < padju)
ALTI <-subset(ALTI, deltaBeta > LogFC)
ALTI <- ALTI[order(ALTI$adj.P.Val), ]
#
BASSI <- subset(results, adj.P.Val < padju)
BASSI <- subset(BASSI, deltaBeta < -LogFC)
BASSI <- BASSI[order(BASSI$adj.P.Val), ]
NO <- subset(results,deltaBeta < LogFC & deltaBeta > -LogFC)

pdf(file=paste("plot_DMA_", args[1],"_",args[2],"_", Cell_type, ".pdf",sep=""))

library(ggplot2)
library(ggrepel)

p <-ggplot() +
        ggtitle( paste("DMA_", Cell_type,"_", args[1], "_vs_", args[2],sep="")) + theme(
                plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
                axis.title.x = element_text(color="black", size=14, face="bold"),
                axis.title.y = element_text(color="black", size=14, face="bold")) +
        geom_point(data = ALTI, aes(x = deltaBeta, y = -log10(adj.P.Val)), color = "red") +
        geom_text_repel(data = ALTI[c(1:25), ], aes(x = deltaBeta , y = -log10(adj.P.Val), label=gene),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1, max.overlaps = 25)+
        geom_point(data = BASSI, aes(x = deltaBeta, y = -log10(adj.P.Val)), color = "blue") +
        geom_text_repel(data = BASSI[1:25, ], aes(x = deltaBeta , y = -log10(adj.P.Val), label=gene),hjust=0, vjust=0,size=3, force = 20, force_pull = 1, max.overlaps = 25)+
        geom_point(data = NO, aes(x = deltaBeta, y = -log10(adj.P.Val)), color = "black")

p2 <- p + geom_vline(xintercept=c(-LogFC,LogFC), col="red") +
        geom_hline(yintercept=-log10(padju), col="red")


print(p2)
dev.off()



write.table(x=ALTI$gene   , file= paste(Cell_type,"_",args[1],"_",args[2],"_UP.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
write.table(x=BASSI$gene  , file= paste(Cell_type,"_",args[1],"_",args[2],"_DOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
write.table(x=ALTI   , file= paste(Cell_type,"_",args[1],"_",args[2],"_GENESUP.tsv",sep="") ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(x=BASSI  , file= paste(Cell_type,"_",args[1],"_",args[2],"_GENESDOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
write.table(Metadata_individual[order(Metadata_individual$CONDITION, decreasing = TRUE),],
            file= paste(args[1], "_",args[2],"_",Cell_type, "_metadata.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)


dir.create(path = paste(directory,"/Methylomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/Pathway_analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory_path <- paste(directory,"/Methylomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2], "/Pathway_analysis", sep ="")


#HEATMAP

if (length(c(ALTI$gene, BASSI$gene)) > 2){
        
        #ADD HEATMAP CREATION HERE
        if (nrow(ALTI) > n_genes_heat){
                ALTI <- ALTI[1:n_genes_heat,]}
        if (nrow(BASSI) > n_genes_heat){
                BASSI <- BASSI[1:n_genes_heat,]}
        
        #HEATMAP INPUT ARE NORMALIZED
        
        
        Heat<-Matrix[rownames(Matrix) %in% c(ALTI$CpG_island, BASSI$CpG_island),]
        Cpg_isl <- rownames(Heat)
        Heat <-apply(Heat,2,as.numeric)
        rownames(Heat) <-Cpg_isl
        
        SAVED <- NULL
        for (ix in Cpg_isl){
        iu<-grep(ix, results$CpG_island)
        Gene <- as.character(results$gene[iu])
        SAVED<-append(Gene, SAVED)}
        SAVED <- as.character(SAVED)
        SAVED<-rev(SAVED)
        samples <-nrow(Heat)
        
        for (ix in 1:samples){
                if(SAVED[ix] !=""){
                        rownames(Heat)[ix] <- SAVED[ix]
                }
        }
        
        library(ComplexHeatmap)
        setwd(directory2)
        pdf(file= paste("Heatmap_top_genes_", Cell_type,"_", args[1],"_vs_", args[2],sep=""))
        #dev.new()
        
        library(circlize)
        col_fun = colorRamp2(c(min(Heat), mean(colMeans(Heat)), max(Heat)), c("blue", "black", "yellow"))
        col_fun(seq(-3, 3))
        
        t<-c("CTRL" = "blue", "SLE" = "red")
        attr(t, "names")[1]<- args[2]
        attr(t, "names")[2]<- args[1]
        attr(t, "names")
        
        ha = HeatmapAnnotation(condition = Metadata_individual$CONDITION,
                               col = list(condition = t))
        p<- Heatmap(Heat,km =2, name = "SD_score", col = col_fun, clustering_distance_rows = "pearson", clustering_method_rows= "complete", clustering_method_columns ="ward.D", row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                    row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
        
        print(p)
        
        #If you want to visualize the group
        colnames(Heat)<- Metadata_individual$CONDITION
        
        ha = HeatmapAnnotation(condition = Metadata_individual$CONDITION,
                               col = list(condition = t))
        p<- Heatmap(Heat,km =2, name = "SD_score", col = col_fun, clustering_distance_rows = "pearson", clustering_method_rows= "complete", clustering_method_columns ="ward.D", row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                    row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
        
        print(p)
        
        dev.off()
        
}




pdf(file=paste(directory_path,"/Pathways_DMG_EnrichR.tsv", sep=""), width = 20, height = 9)


#PATHWAY ANALYSIS
library(enrichR)
dbs <- c("Reactome_2022", "GO_Biological_Process_2023", "CODE_and_ChEA_Consensus_TFs_from_ChIP-X")
websiteLive <- getOption("enrichR.live")

if(length(ALTI$gene) != 0){
        if (websiteLive) {
                enriched <- enrichr(as.character(ALTI$gene), dbs)}
        for(paths in 1:3){
                if (nrow(enriched[[paths]]) != 0){
                        if (websiteLive) {
                          print(enriched[[paths]])
                                dir.create(paste(directory_path,"/TABLES", sep=""))
                                write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",names(enriched)[paths],".tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
                                xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_upregulated_genes","_",dbs[paths], sep=""))
                                print(xx)
                        }}}}

if(length(BASSI$gene) != 0){
        if (websiteLive) {
                enriched <- enrichr(as.character(BASSI$gene), dbs)}
        for(paths in 1:3){
                if (nrow(enriched[[paths]]) != 0){
                        if (websiteLive) {
                                xx<-plotEnrich(enriched[[paths]], showTerms = 25 , numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_downregulated_genes","_",dbs[paths], sep=""))
                                print(xx)
                        }}}}

dev.off()


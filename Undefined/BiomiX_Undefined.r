library(vroom)
library(dplyr)
library(stringr)
library(rlist)
library(tibble)
library(readxl)

# MANUAL INPUT
# # # #
# library(vroom)
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"mutated"
# args[2] <-"unmutated"
# args[3] <-"/home/cristia/BiomiX2.2"
# 
# directory <- args[3]
# selection_samples = "NO"
# Cell_type = "populations"
# DIR_METADATA <- readLines("/home/cristia/BiomiX2.2/directory.txt")
# i = 2

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")


COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
Heatmap_genes <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])

LogFC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[2])


#loading matrix
directory2 <- paste(directory,"/Undefined/INPUT",sep="")
setwd(directory2)

print(i)
print(COMMAND$LABEL[i])

if (grepl("\\.xlsx$|\\.xls$", COMMAND$DIRECTORIES[i])) {
        matrix <- read_excel(COMMAND$DIRECTORIES[i])
        print("Matrix Excel File read successfully!")
}else{
        matrix <- vroom(COMMAND$DIRECTORIES[i],delim="\t",col_names = TRUE, comment = "#")
}


sam <- colnames(matrix)
pea <-matrix$ID
matrix <- t(matrix[,-1])
matrix <- add_column(as.data.frame(matrix), sam[-1], .after = 0) 
colnames(matrix) <-c("ID",pea)


#create directory
dir.create(path = paste(directory,"/MOFA/INPUT/", "Undefined_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory, "/MOFA/INPUT/", "Undefined_", Cell_type,  "_",args[1],"_vs_", args[2], sep ="")

setwd(directory2)


if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
        Metadata_total <- read_excel(DIR_METADATA)
        print("Metadata Excel File read successfully!")
}else{
        Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)}


#If statement in case of sample selection

if (selection_samples == "YES") {
        num <- grep(Cell_type, Metadata_total$ID_CELL_TYPE)
        Metadata_Bcell <- Metadata_total[num,]
        Metadata_total <- Metadata_total[num,]
        Identifier<-matrix$ID
        num <-which(colnames(matrix) %in% Metadata_Bcell$ID_CELL_TYPE )
        matrix <- as.matrix(matrix[,num])
        print(Metadata_Bcell$ID_CELL_TYPE)
        num <-which(Metadata_Bcell$ID_CELL_TYPE  %in% colnames(matrix) )
        Metadata_Bcell <- Metadata_Bcell[num,]
        #If multiple select only the first column
        Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
        
        Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
        matrix<-matrix[,order(colnames(matrix))]
        
} else {
        print("No samples selection")
        Metadata_Bcell <- Metadata_total
        Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
        matrixi <- matrix[,2:ncol(matrix)]
        matrix[,2:ncol(matrix)]<-matrixi[,order(colnames(matrixi))]
        matrix<-as.data.frame(matrix)
        
        Metadata_Bcell <- Metadata_total
        num <-which(matrix$ID %in% Metadata_Bcell$ID)
        Identifier<-matrix$ID
        matrix <- as.matrix(matrix[num,])
        print(Metadata_Bcell$ID)
        num <-which(Metadata_Bcell$ID  %in% matrix[,1] )
        Metadata_Bcell <- Metadata_Bcell[num,]
        #If multiple select only the first column
        Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
        
        Metadata_Bcell<-Metadata_Bcell[order(Metadata_Bcell$ID),]
        matrix<-matrix[order(matrix[,1]),]
}

#The metadata file contains the metadata of only the samples within the omics
#Excluding the ones not included in the original metadata file.
Metadata <- Metadata_Bcell
Metadata_individual=NULL
Metadata_reads=NULL
Metadata_Bcell=NULL
Metadata_tot = NULL
Metadata_total = NULL

#Check if metadata and matrix contain the same samples.
print("Same samples in matrix and metadata?")
print(all(Metadata$ID == matrix[,1]))
outs <- Metadata$CONDITION  == args[2]|Metadata$CONDITION  == args[1]
Metadata<- Metadata[outs,]
matrix<- matrix[outs,]


#FILTERING BASED ON METADATA
#Based on the criteria defined in the advance options.

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
                        matrix <- matrix[,comparison_operator(To_filter, value_threshold)]
                        
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
                        DGE3 <- DGE3[,comparison_operator(To_filter, value_threshold)]
                }
                
        }else{
                print("skip_filtering")
        }
}

###
matrix<- as.data.frame(matrix)
rownames(matrix) <- matrix$ID
###


#This part of code is to refill the column having an NA name and rename the duplicated ones
na_cols <- which(is.na(colnames(matrix)))
new_names <- paste("UNKNOWN",1:sum(is.na(colnames(matrix))),sep="") 
colnames(matrix)[na_cols] <- new_names

colnames(matrix) <- make.unique(names(matrix), sep = "_")

#====================================

matrixs <- matrix 


matrixs <- add_column(matrixs, Metadata$CONDITION, .after = 1)
colnames(matrixs)[2] <- "CONDITION"


tryCatch({
        # If there is a missmatch between metadata and matrix this line will not work
        matrixs[,3:ncol(matrixs)]<-apply(matrixs[,3:ncol(matrixs)],2,as.numeric)
        
}, error = function(e) {
        # Personalized error message
        print(paste("Oops! An error occurred:", e$message))
        print("Are you sure the Undefined matrix samples are the same of the metadata file?")
})



#Calculus variance
varianze = NULL
Vari <- as.matrix(matrixs[,c(-1:-2)])
Vari <- as.data.frame(Vari)

#generazione e filtraggio per varianza
for (i in 1:ncol(Vari)){
        varianze <- append(varianze, var(as.numeric(Vari[,i])))
}


normalizzati<- as.data.frame(matrixs)
variance <- log10(varianze)

ordering<-data.frame(colnames(Vari), variance)
ordering$variance<-abs(ordering$variance)
ordering <- ordering[order(ordering$variance, decreasing = TRUE), ]
limit<-as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA[1])

#SAVING THE INPUT FOR THE MOFA ANALYSIS
#These data should not be normalized by log.

if(length(rownames(ordering)) > limit){
        ordering<-ordering[!ordering$variance == Inf,]
        normalizzati_MOFA <- ordering[1:limit,]
        ond<-colnames(matrixs) %in% normalizzati_MOFA[,1] 
        ond[1:2] <- TRUE
        normalizzati_MOFA<-matrixs[,ond]
}else{
normalizzati_MOFA <- matrixs}

write.table(normalizzati_MOFA,paste(directory2,"/Undefined_",Cell_type, "_MOFA.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")



# #### FUNCTION DEFINITION FOR CONTROL AND TESTED SAMPLES----


DMS_REFERENCE <- function(Condition1){
        print(Condition1)
        CON <- matrix[Metadata$CONDITION == Condition1,] #Example CONTROL vs SLE patients
        samp_CON<-row.names(CON) 
        CON <-CON[,c(-1,-2)] 
        CON <-apply(CON,2,as.numeric)
        

        
        y=NULL
        for (i in seq(1:length(colnames(CON)))) {
                if (all_same(CON[,i])){
                        y <-append(y, 1)
                }else{
                res <- (shapiro.test(CON[,i]))
                y <-append(y, res[["p.value"]])
        }}
        
        CON<-as.data.frame(CON)
        shapiro_test_CON <- y
        shapiro <-shapiro_test_CON < 0.05
        print(paste(length(shapiro[shapiro== TRUE]), "variables out of", length(shapiro), "not having a normal distribution in", args[2]))
        if (length(shapiro[shapiro== TRUE]) > length(shapiro)/2 ) {
                print(paste("Suggested Log normalization before the MOFA analysis for", args[2]))
        }
        
        rownames(CON) <-samp_CON
        return(CON)
}





DMS_TESTED <- function(Condition2){
        print(Condition2)
        SLE <- matrix[Metadata$CONDITION == Condition2,] 
        samp_SLE<-row.names(SLE) 
        SLE <-SLE[,c(-1,-2)]
        SLE <-apply(SLE,2,as.numeric) 
        
        y=NULL
        for (i in seq(1:length(colnames(SLE)))) {
                if (all_same(SLE[,i])){
                        y <-append(y, 1)
                }else{
                        res <- (shapiro.test(SLE[,i]))
                        y <-append(y, res[["p.value"]])
                }}
        
        SLE<-as.data.frame(SLE)
        shapiro_test_SLE <- y
        shapiro <-shapiro_test_SLE < 0.05
        length(shapiro[shapiro== TRUE]) 
        print(paste(length(shapiro[shapiro== TRUE]), "variables out of", length(shapiro), "not having a normal distribution in",args[1]))
        if (length(shapiro[shapiro== TRUE]) > length(shapiro)/2 ) {
                print(paste("Suggested Log normalization before the MOFA analysis for", args[1]))
        }
        #If the p.value is bigger than 0.05 we don't have difference
        #between our and a normal distribution (generated randomly)
        
        rownames(SLE) <-samp_SLE
        
        # Wilcox test CON vs SLE
        pval=NULL
        for (i in 1:ncol(SLE)) {
                res <-wilcox.test(CON[,i],SLE[,i], alternative = "two.sided") 
                pval <-append(pval, res[["p.value"]])
        }
        
        # T test CON vs SLE
        pval_t=NULL
        for (i in 1:ncol(SLE)) {
                res_t <-t.test(CON[,i],SLE[,i], alternative = "two.sided") 
                pval_t <-append(pval_t, res[["p.value"]])
        }
        
        
        
        fold=NULL
        for (i in 1:ncol(CON)) {
                up <-CON[,i]
                up <-up[!is.na(up)]
                down <-SLE[,i]
                down <-down[!is.na(down)]
                FC <- log2(abs(median(down) / median(up)))
                fold <-append(fold, FC)}
        
        
        SLE<-as.data.frame(t(SLE)) #RUN THESE LINES FOR THE CTRL VS SLE COMPARISON
        colnames(SLE) <- samp_SLE
        SLE$shapiro_pvalue <- shapiro_test_SLE
        SLE$p_val_wilcox <- pval
        SLE$p_val_t_test <- pval_t
        SLE$log2FC <- fold
        
        return(SLE)
}


#Function to check if all the variable elements are the same
all_same <- function(x) {
        all(x == x[1])
}



# #### DMS DIFFERENTIAL METABOLITE SIGNALS ANALYSIS----

CON <-DMS_REFERENCE(args[2])        
TEST <- DMS_TESTED(args[1])        
TEST$NAME <- row.names(TEST)
gc()

#Saving results
        
        dir.create(path = paste(directory,"/Undefined/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        dir.create(path = paste(directory,"/Undefined/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        directory2 <- paste(directory,"/Undefined/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="")
        
        setwd(directory2)
        
        total <- TEST
        total$padj_wilcox = p.adjust(total$p_val_wilcox, method = "fdr")
        total$padj_t_test = p.adjust(total$p_val_t_test, method = "fdr")
        totalshOK <- total %>% filter(padj_t_test < padju)
        total <- total %>% arrange(padj_t_test)
        x <- colnames(total) %in% matrix$ID
        total <- total[,!x]
        total_3 <- total
        total_3 <- total_3[,-2:-3]
        write.table(total_3,paste(directory2,"/",Cell_type,"_",args[1],"_vs_",args[2],"_results.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
        
        
        total$padj_wilcox<-as.numeric(total$padj_wilcox)
        total$log2FC<-as.numeric(total$log2FC)
        
        
        #PLOTS SECTIONS
        
        ALTI <- subset(total, padj_wilcox < padju)
        ALTI <-subset(ALTI, log2FC > LogFC)
        ALTI <- ALTI[order(ALTI$padj_wilcox), ]
        #
        BASSI <- subset(total, padj_wilcox < padju)
        BASSI <- subset(BASSI, log2FC < -LogFC)
        BASSI <- BASSI[order(BASSI$padj_wilcox), ]
        x <- total$NAME %in% ALTI$NAME
        NO<-total[!x,]
        x <- NO$NAME %in% BASSI$NAME
        NO<-NO[!x,]
        
        
        write.table(x=ALTI   , file= paste(Cell_type,"_",args[1],"_",args[2],"_VARIABLES_UP.tsv",sep="") ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        write.table(x=BASSI  , file= paste(Cell_type,"_",args[1],"_",args[2],"_VARIABLES_DOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE),c("ID","CONDITION")],
                    file= paste(args[1], "_",args[2],"_",Cell_type, "_metadata.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
        
        
        
        if (nrow(ALTI) == 0){
                SKIP_ALTI <- FALSE}else{SKIP_ALTI <- TRUE}
        if (nrow(BASSI) == 0){
                SKIP_BASSI <- FALSE}else{SKIP_BASSI <- TRUE}
        
        pdf(file=paste("plot_Undefined_", args[1],"_",args[2],"_", Cell_type, ".pdf", sep=""))
        
        library(ggplot2)
        library(ggrepel)
        
        
        p <-ggplot() +
                ggtitle( paste("Variables analysis", Cell_type,"_", args[1], "_vs_", args[2],sep="")) + theme(
                        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
                        axis.title.x = element_text(color="black", size=14, face="bold"),
                        axis.title.y = element_text(color="black", size=14, face="bold")) +
                ylab("-log10(padj_wilcox)") + xlab("log2FC")
        if (SKIP_ALTI){
                p <- p + geom_point(data = ALTI, aes(x = log2FC, y = -log10(padj_wilcox)), color = "red") +
                        geom_text_repel(data = ALTI[c(1:25), ], aes(x = log2FC , y = -log10(padj_wilcox), label=NAME),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1, max.overlaps = 5)}
        if (SKIP_BASSI){
                p <- p + geom_point(data = BASSI, aes(x = log2FC, y = -log10(padj_wilcox)), color = "blue") +
                        geom_text_repel(data = BASSI[1:25, ], aes(x = log2FC , y = -log10(padj_wilcox), label=NAME),hjust=0, vjust=0,size=3, force = 20, force_pull = 1, max.overlaps = 5)}
        p <- p + geom_point(data = NO, aes(x = log2FC, y = -log10(padj_wilcox)), color = "black")
        
        p2 <- p + geom_vline(xintercept=c(-LogFC,LogFC), col="red") +
                geom_hline(yintercept=-log10(padju), col="red")
        
        print(p2)
        
        
        
        ALTI_2 <- subset(total, padj_t_test < padju)
        ALTI_2 <-subset(ALTI_2, log2FC > LogFC)
        ALTI_2 <- ALTI_2[order(ALTI_2$padj_t_test), ]
        #
        BASSI_2 <- subset(total, padj_t_test < padju)
        BASSI_2 <- subset(BASSI_2, log2FC < -LogFC)
        BASSI_2 <- BASSI_2[order(BASSI_2$padj_t_test), ]
        x <- total$NAME %in% ALTI_2$NAME
        NO_2<-total[!x,]
        x <- NO_2$NAME %in% BASSI_2$NAME
        NO_2<-NO_2[!x,]

        if (nrow(ALTI_2) == 0){
                SKIP_ALTI_2 <- FALSE}else{SKIP_ALTI_2 <- TRUE}
        if (nrow(BASSI_2) == 0){
                SKIP_BASSI_2 <- FALSE}else{SKIP_BASSI_2 <- TRUE}        
        
        p <-ggplot() +
                ggtitle( paste("Variables analysis", Cell_type,"_", args[1], "_vs_", args[2],sep="")) + theme(
                        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
                        axis.title.x = element_text(color="black", size=14, face="bold"),
                        axis.title.y = element_text(color="black", size=14, face="bold")) +
                ylab("-log10(padj_t_test)") + xlab("log2FC")
        if (SKIP_ALTI_2){
                p <- p + geom_point(data = ALTI_2, aes(x = log2FC, y = -log10(padj_t_test)), color = "red") +
                        geom_text_repel(data = ALTI_2[c(1:25), ], aes(x = log2FC , y = -log10(padj_t_test), label=NAME),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1, max.overlaps = 5)}
        if (SKIP_BASSI_2){
                p <- p + geom_point(data = BASSI_2, aes(x = log2FC, y = -log10(padj_t_test)), color = "blue") +
                        geom_text_repel(data = BASSI_2[1:25, ], aes(x = log2FC , y = -log10(padj_t_test), label=NAME),hjust=0, vjust=0,size=3, force = 20, force_pull = 1, max.overlaps = 5)}
        p <- p + geom_point(data = NO_2, aes(x = log2FC, y = -log10(padj_t_test)), color = "black")
        
        p2 <- p + geom_vline(xintercept=c(-LogFC,LogFC), col="red") +
                geom_hline(yintercept=-log10(padju), col="red")
        
        print(p2)
        
        
        
        dev.off()
        
        
        ###  HEATMAP ###
        
        if (length(c(ALTI$NAME, BASSI$NAME)) > 2){
                
                #ADD HEATMAP CREATION HERE
                if (nrow(ALTI) > Heatmap_genes){
                        ALTI <- ALTI[1:Heatmap_genes,]}
                if (nrow(BASSI) > Heatmap_genes){
                        BASSI <- BASSI[1:Heatmap_genes,]}
                
                #HEATMAP INPUT ARE NORMALIZED
                
                
                Heat<-matrix[,colnames(matrix) %in% c(ALTI$NAME, BASSI$NAME)]
                samples <- rownames(Heat)
                Heat <-apply(Heat,2,as.numeric)
                rownames(Heat) <-samples
                Heat<-t(Heat)
                Heat[Heat == 0] <- 1
                Heat <-apply(Heat,2,log10)
                
                library(ComplexHeatmap)
                setwd(directory2)
                pdf(file= paste("Heatmap_top_variables_", Cell_type,"_", args[2],"_vs_", args[1],".pdf",sep=""))
                #dev.new()
                
                library(circlize)
                col_fun = colorRamp2(c(min(Heat), mean(colMeans(Heat)), max(Heat)), c("blue", "black", "yellow"))
                col_fun(seq(-3, 3))
                
                t<-c("CTRL" = "blue", "SLE" = "red")
                attr(t, "names")[1]<- args[2]
                attr(t, "names")[2]<- args[1]
                attr(t, "names")
                
                ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                                       col = list(condition = t))
                p<- Heatmap(Heat,km =2, name = "SD_score", col = col_fun, clustering_distance_rows = "pearson", clustering_method_rows= "complete", clustering_method_columns ="ward.D", row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                            row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
                
                print(p)
                
                #If you want to visualize the group
                colnames(Heat)<- Metadata$CONDITION
                
                ha = HeatmapAnnotation(condition = Metadata$CONDITION,
                                       col = list(condition = t))
                p<- Heatmap(Heat,km =2, name = "SD_score", col = col_fun, clustering_distance_rows = "pearson", clustering_method_rows= "complete", clustering_method_columns ="ward.D", row_dend_width = unit(0.5, "cm"),column_dend_height = unit(60, "mm"), column_names_gp = grid::gpar(fontsize = 6),
                            row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha)
                
                print(p)
                
                dev.off()
                
}

        
                
        



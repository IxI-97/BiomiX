
cat("\n\n\n\n\n          /////////      ///    /////////    ///       ///     ///    ///      ///\n         ///    ///     ///    ///   ///    /////  //////     ///      ///  ///  \n        /////////      ///    ///   ///    ///  ///  ///     ///        ///    \n       ///    ///     ///    ///   ///    ///  ///  ///     ///      ///  ///  \n      ///////////    ///    /////////    ///  ///  ///     ///    ///      ///\n                                                               ///            ///\n                                                                                 /////\n\n\n\n\n\n")

#directory <- "/media/henry/My_Passport/Article_MULTI_BETA"

# #MANUAL INPUT
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"SJS"
# args[2] <-"CTRL"
# args[3] <-"/home/cristia/BiomiX2.2"
library(vroom)


#Sys.sleep(5)
print("Welcome to BiomiX toolkit")

#Sys.sleep(5)
Cell_type = readline(prompt = "Which sorted cell type do you want to analyze? ^-^ : ")
Cell_Disease = readline(prompt = "Which disease type do you want to analyze? ^-^ : ")

#Sys.sleep(5)

print("Information acquired, running the analysis")

#Sys.sleep(5)

#if(args[1] == ""){
        # taking input with showing the message
        args = commandArgs(trailingOnly=TRUE)
#}

#directory <- "/home/cristia/Scrivania/PhD/Bioinfo/Article_MULTI_BETA"
directory <- args[[3]]
setwd(directory)
#setwd(directory)


#COMMAND <- vroom(paste(directory,"Commands",sep="/"), delim = "\t")
COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
DIR_METADATA <- readLines("directory.txt")
DIR_METADATA_output <- readLines("directory_out.txt")

#args = as.list(c("BLymphocytes","SJS"))
print(args)
#print(args[1])
#print(args[2])
#print(args[3])


#FIND A WAY TO SELECT ANY POSSIBLE INPUT WITHOUT CONSIDERING THE ORDER

n_iteration<- sum(COMMAND$DATA_TYPE == "Transcriptomics")
position <- which(COMMAND$DATA_TYPE %in% "Transcriptomics")
iterations <- 0
for (i in position){
if(COMMAND$ANALYSIS[i] == "YES"){
        
        
        Cell_type <- COMMAND$LABEL[i]
        cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
        iterations = iterations + 1
        selection_samples = COMMAND$SELECTION[i]
        #purity_filter =COMMAND$PURITY[i]
        directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
        #files <- grep("MATRIX*",list.files(directory2),value=TRUE)
        #files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
        source(paste(directory,"/Transcriptomics/Biomix_DGE_GENES_LIMMA.r",sep=""))
        cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
        gc()
}
}




#===============================================================================


#Sys.sleep(5)

n_iteration<- sum(COMMAND$DATA_TYPE == "Metabolomics")
position <- which(COMMAND$DATA_TYPE %in% "Metabolomics")
iterations <- 0
for (i in position){
        
        Cell_type <- COMMAND$LABEL[i]
        cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
        iterations = iterations + 1
        selection_samples = COMMAND$SELECTION[i] # TO ADD
        #purity_filter =COMMAND$PURITY[i] #TO ADD
        directory2 <- paste(directory,"/Transcriptomics/INPUT_PRECISESADS",sep="")
        files <- grep("MATRIX*",list.files(directory2),value=TRUE)
        files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
        if(COMMAND$ANALYSIS[i] == "YES"){
        source(paste(directory,"/Metabolomics/BiomiX_DMA.r",sep = ""))
        cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
        gc()
}
}

#================================================================================

n_iteration<- sum(COMMAND$DATA_TYPE == "Methylomics")
position <- which(COMMAND$DATA_TYPE %in% "Methylomics")
iterations <- 0
for (i in position){
        
        Cell_type <- COMMAND$LABEL[i]
        cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
        iterations = iterations + 1
        selection_samples = COMMAND$SELECTION[i] # TO ADD
        #purity_filter =COMMAND$PURITY[i] #TO ADD
        directory2 <- paste(directory,"/Transcriptomics/INPUT_PRECISESADS",sep="")
        files <- grep("MATRIX*",list.files(directory2),value=TRUE)
        files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
        if(COMMAND$ANALYSIS[i] == "YES"){
                source(paste(directory,"/Methylomics/BiomiX_DMA.r",sep=""))
                cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
                gc()
        }
}

#================================================================================



if(COMMAND_MOFA[2,2] == "YES") {
        
        if(COMMAND_MOFA[3,2] == 0) {
                
                cat(paste( "\n\n\n\n\nStarting the AUTOMATIC MOFA analysis \n\n\n\n\n", sep =""))
                Sys.sleep(5)
                source(paste(directory,"/MOFA/MOFA_MULTI2_AUTO.R",sep=""))
                cat("\n\n\n\n\n AUTOMATIC MOFA analysis complete ^-^")
                
        }else{
                
                cat(paste( "\n\n\n\n\nStarting the MOFA analysis \n\n\n\n\n", sep =""))
                Sys.sleep(5)
                source(paste(directory,"/MOFA/MOFA_MULTI2.R",sep=""))
                cat("\n\n\n\n\n MOFA analysis complete ^-^")
}

}


#================================================================================

cat(paste( "\n\n\n\n\nStarting the MOFA factors correlation vs clinical factors \n\n\n\n\n", sep =""))

#Sys.sleep(5)

if(COMMAND_MOFA[2,2] == "YES") {
        source(paste(directory,"/Clinical_data/BiomiX_Clinical.R",sep=""))
        cat("\n\n\n\n\n MOFA factors correlation vs clinical factors completed ^-^")
}



#================================================================================

cat(paste( "\n\n\n\n\nStarting the MOFA factors articles matching \n\n\n\n\n", sep =""))

#Sys.sleep(5)

if(COMMAND_MOFA[2,2] == "YES") {
        source(paste(directory,"/MOFA/BiomiX_PUBMED.R",sep=""))
        cat("\n\n\n\n\n MOFA factors articles matching completed ^-^")
}


#================================================================================

cat(paste( "\n\n\n\n\nStarting the MOFA factors pathway analysis \n\n\n\n\n", sep =""))

#Sys.sleep(5)

if(COMMAND_MOFA[2,2] == "YES") {
        source(paste(directory,"/MOFA/BiomiX_MOFA_MINER.R",sep=""))
        cat("\n\n\n\n\n MOFA factors pathway analysis completed ^-^")
}


#================================================================================
#args[1] = "ICPLUS"

cat(paste( "\n\n\n\n\n Saving the files in the directory \n\n\n\n\n", sep =""))

Sys.sleep(5)

#directory2 <- paste(directory, "/",COMMAND$DATA_TYPE[1], "/OUTPUT/",COMMAND$LABEL[1], "_",args[1], "_vs_", args[2], sep ="")
#print(directory2)
#setwd(directory2)

#COPY OMICS DIRECTORY
for(output in 1:length(na.omit(COMMAND$LABEL))) {
        print(output)
        print(paste(directory,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))
        if(file.exists(paste(directory,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))){
                dir.create(path =  paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                file.copy(paste(directory,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""), paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
                print(paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))
}
}

if(COMMAND_MOFA[2,2] == "YES") {

#COPY MOFA FILES
output<-list.files(path = paste(directory,"/","MOFA", "/OUTPUT/",sep =""), pattern = paste("^","AUTOMATIC_SEARCH_BEST_FACTORS_",args[1], "_vs_", args[2], sep=""))
if(file.exists(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""))){
        dir.create(path =  paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        file.copy(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep =""), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE) 
        print(paste(output,"copied", sep = " "))
}

#COPY MOFA DIRECTORY
file_to_transfer<-dir(path = paste(directory,"/","MOFA", "/OUTPUT/",sep =""), pattern = paste("^","MOFA_",args[1], "_vs_", args[2], sep="") )
for(output in file_to_transfer){
        if(file.exists(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""))){
                dir.create(path =  paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                file.copy(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
                #unlink(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""), recursive=TRUE)
                print(paste(output,"copied", sep = " "))
        }}

output<-dir(path = paste(directory,"/","Clinical_data", "/OUTPUT/",sep =""), pattern = paste("^","Clinical_",args[1], "_vs_", args[2], sep="") )
        if(file.exists(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""))){
                dir.create(path =  paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                file.copy(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
                #unlink(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""), recursive=TRUE)
                print(paste(output,"copied", sep = " "))
        }

}
cat("\n\n\n\n\n file saving complete ^-^")

matr<-data.frame(c(17:2))
matr[1:3,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS
matr[1:3,1]<-c("LOG2FC_TRANSCRIPTOMICS", "P.ADJ_TRANSCRIPTOMICS", "GENE_PANEL_FILE")
matr[4:6,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS
matr[4:6,1]<-c("LOG2FC_METABOLOMICS", "P.ADJ_METABOLOMICS", "CPU_USED")
matr[7:9,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS
matr[7:9,1]<-c("LOG2FC_METHYLOMICS", "P.ADJ_METHYLOMICS", "ARRAY_TYPE")
matr[10:12,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING
matr[10:12,1]<-c("PANEL POSITIVITY WITH N°GENES WITHIN THE PANEL WITH Z-SCORE > 2", "PANEL POSITIVITY WITH N°GENES WITHIN THE PANEL WITH Z-SCORE > 1", "REMOVE CONTROL POSITIVE FOR GENE PANEL?")
matr[13:15,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS
matr[13:15,1]<-c("CLUSTERING DISTANCE", "CLUSTERING METHOD", "N° GENES TO VISUALIZE IN THE HEATMAP?")
matr[16,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3]
matr[16,1]<-c("N° MOFA INPUT FEATURES")


matr_2<-data.frame(c(30:2))
matr_2[1:2,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL[1:2]
matr_2[1:2,1]<-c("LEVEL_ANNOTATION_METABOLOMICS_DATASET?", "IF_ANNOTATED_WHICH_TYPE_OF ANNOTATION?")
matr_2[3:5,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1
matr_2[3:5,1]<-c("ION_MODE_ONLY_MS1_MODE", "NEURAL_MODE? _ONLY_MS1_MODE", "TOLERANCE_IN_PPM_ONLY_MS1_MODE")
matr_2[6:8,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2
matr_2[6:8,1]<-c("POSITIVE_ADDUCTS_ONLY_MS1_MODE", "NEGATIVE_ADDUCTS_ONLY_MS1_MODE", "DATABASE_FOR_MS1_ANNOTATION_ONLY_MS1_MODE")
matr_2[9:11,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX
matr_2[9:11,1]<-c("INDEXING_1", "INDEXING_2", "INDEXING_3")
matr_2[12:14,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES
matr_2[12:14,1]<-c("FILE_MS1_ANNOTATION_INDEXING_1_ONLY_MS1_MODE", "FILE_MS1_ANNOTATION_INDEXING_2_ONLY_MS1_MODE", "FILE_MS1_ANNOTATION_INDEXING_3_ONLY_MS1_MODE")
matr_2[15:17,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2
matr_2[15:17,1]<-c("TOLERANCE PPM MATCH MS1 MS2 PEAKS", "RETENTION TIME MATCH MS1 MS2", "MS2 DATABASE ANNOTATION PRIORITY")
matr_2[18,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3[1]
matr_2[18,1]<-c("ION_MODE_ONLY_MS1-2_MODE")
matr_2[19:21,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2
matr_2[19:21,1]<-c("POSITIVE_ADDUCTS__MS1-2_MODE", "NEGATIVE_ADDUCTS__MS1-2_MODE", "TYPE LC COLUMN")
matr_2[22:23,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4[1:2]
matr_2[22:23,1]<-c("DATABASE_FOR_MS1_ANNOTATION_ONLY_MS1-2_MODE", "TOLERANCE_IN_PPM_ONLY_MS1-2_MODE")
matr_2[24:26,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX
matr_2[24:26,1]<-c("FILE_MS1_ANNOTATION_INDEXING_1_ONLY_MS1-2_MODE", "FILE_MS1_ANNOTATION_INDEXING_2_ONLY_MS1-2_MODE", "FILE_MS1_ANNOTATION_INDEXING_3_ONLY_MS1-2_MODE")
matr_2[27:29,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2
matr_2[27:29,1]<-c("FILE_MS1_ANNOTATION_INDEXING_1_ONLY_MS1-2_MODE", "FILE_MS1_ANNOTATION_INDEXING_2_ONLY_MS1-2_MODE", "FILE_MS1_ANNOTATION_INDEXING_3_ONLY_MS1-2_MODE")
matr_2[30,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1]
matr_2[30,1]<-c("DIRECTORY_TO_MS2_FILES")


matr_3<-data.frame(c(12:2))
matr_3[1:3,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_1
matr_3[1:3,1]<-c("FILTERING_METADATA_1", "TYPE_OF_DATA_1", "THRESHOLD_OR_NAME_1")
matr_3[4:6,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_2
matr_3[4:6,1]<-c("FILTERING_METADATA_2", "TYPE_OF_DATA_2", "THRESHOLD_OR_NAME_2")
matr_3[7:9,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_3
matr_3[7:9,1]<-c("FILTERING_METADATA_3", "TYPE_OF_DATA_3", "THRESHOLD_OR_NAME_3")
matr_3[10:12,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_4
matr_3[10:12,1]<-c("FILTERING_METADATA_4", "TYPE_OF_DATA_4", "THRESHOLD_OR_NAME_4")


matr_4<-data.frame(c(13:2))
matr_4[1:3,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_MOFA
matr_4[1:3,1]<-c("MAX_ITERATION", "CONVERGENCE_SPEED", "THRESHOLD_CONTRIBUTION_WEIGHT")
matr_4[4:5,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[1:2]
matr_4[4:5,1]<-c("TYPE_OF_RESEARCH", "N°_OF_ARTICLES")
matr_4[6:7,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[1:2]
matr_4[6:7,1]<-c("N°_TOP_CONTRIBUTORS_TO_SEARCH_IN_PUBMED", "N°_KEYWORDS_EXTRACTED_BY_THE_KEYWORD_GENERATOR")
matr_4[8:9,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY[1:2]
matr_4[8:9,1]<-c("P.ADJ_THRESHOLD_PATHWAY", "N°_OF_PATHWAY_TO_VISUALIZE")
matr_4[10:11,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL[1:2]
matr_4[10:11,1]<-c("NUMERIC_CLINICAL_DATA_UPLOADED?", "BINARY_CLINICAL_DATA_UPLOADED")
matr_4[12:13,2]<-COMMAND_ADVANCED$ADVANCED_OPTION_CLINIC_DATA_DIRECTORY[1:2]
matr_4[12:13,1]<-c("NUMERIC_CLINICAL_DATA_FILE?", "BINARY_CLINICAL_DATA_FILE")

dir.create(path =  paste(directory,"/","Report_parameters/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
line="#PARAMETERS OMICS INPUT#"
write(line,file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE)
write.table(COMMAND[,-1], file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE , row.names = F)
line="#PARAMETERS INTEGRATION#"
write(line,file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE)
COMMAND_MOFA[,1]<-c("TYPE_INTEGRATION", "MULTIOMICS INTEGRATION?", "N° MOFA FACTOR CALCULATED? 0 = AUTOMATIC MODE", "MOFA FACTOR TO EXPLORE WITH VISUALIZATION", "MIN OMICS SHARED PER SAMPLE")
write.table(COMMAND_MOFA, file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE, col.names = F, row.names = F)
line="#ADVANCED OPTION PARAMETER (GENERAL)#"
write(line,file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE)
write.table(matr, file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE , col.names = F, row.names = F)
line="#ADVANCED OPTION PARAMETER (METABOLOMICS)#"
write(line,file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE)
write.table(matr_2, file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE , col.names = F, row.names = F)
line="#ADVANCED OPTION PARAMETER (METADATA)#"
write(line,file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE)
write.table(matr_3, file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE , col.names = F, row.names = F)
line="#ADVANCED OPTION PARAMETER (MOFA)#"
write(line,file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE)
write.table(matr_3, file= paste(directory,"/","Report_parameters/","Report", args[1],"_vs_", args[2], substr(Sys.time(), 1,16),sep=""),append=TRUE, col.names = F, row.names = F)


cat("\n\n\n\n\n Report parameters saved ^-^")

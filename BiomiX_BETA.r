
cat("\n\n\n\n\n          /////////      ///    /////////    ///       ///     ///    ///      ///\n         ///    ///     ///    ///   ///    /////  //////     ///      ///  ///  \n        /////////      ///    ///   ///    ///  ///  ///     ///        ///    \n       ///    ///     ///    ///   ///    ///  ///  ///     ///      ///  ///  \n      ///////////    ///    /////////    ///  ///  ///     ///    ///      ///\n                                                               ///            ///\n                                                                                 /////\n\n\n\n\n\n")

#directory <- "/media/henry/My_Passport/Article_MULTI_BETA"

# #MANUAL INPUT
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"SJS"
# args[2] <-"CTRL"
# args[3] <-"/home/cristia/Scrivania/BiomiX2.2"
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
                unlink(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""), recursive=TRUE)
                print(paste(output,"copied", sep = " "))
        }}

output<-dir(path = paste(directory,"/","Clinical_data", "/OUTPUT/",sep =""), pattern = paste("^","Clinical_",args[1], "_vs_", args[2], sep="") )
        if(file.exists(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""))){
                dir.create(path =  paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                file.copy(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
                unlink(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""), recursive=TRUE)
                print(paste(output,"copied", sep = " "))
        }

}
cat("\n\n\n\n\n file saving complete ^-^")

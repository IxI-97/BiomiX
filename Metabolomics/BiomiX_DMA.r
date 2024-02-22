library(vroom)
library(dplyr)
library(cmmr)
library(stringr)
library(rlist)
library(tibble)

#Testing if the number of argument added in the interface is correct (disease,control and directory)

if (length(args) == 0) {
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

# MANUAL INPUT
# # #
# library(vroom)
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"cachexic"
# args[2] <-"control"
# args[3] <-"/home/cristia/Scrivania/BiomiX2.2"
# 
# directory <- args[3]
# iterations = 1
# selection_samples = "NO"
# Cell_type = "Cachexia"
# i = 1
# ANNOTATION = "Annotated"
# DIR_METADATA <- readLines("/home/cristia/Scrivania/BiomiX2.2/directory.txt")


# library(vroom)
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"cachexic"
# args[2] <-"control"
# args[3] <-"/home/cristia/Scrivania/BiomiX2.1"
# #
# directory <- args[3]
# iterations = 1
# selection_samples = "NO"
# Cell_type = "CACHEXIA"
# i = 1
# ANNOTATION = "Annotated"

# DIR_METADATA <- readLines("/home/cristia/Scrivania/BiomiX2.1/directory.txt")
# matrix <-vroom("/home/henry/Desktop/BiomiX2.0/Metabolomics/INPUT/ORIGINAL/PLASMA_STRD.tsv" , delim = "\t", col_names = TRUE)
# matrix <-vroom("/home/cristia/Scrivania/BiomiX1.9/Metabolomics/INPUT/PLASMA_SLE.tsv" , delim = "\t", col_names = TRUE)
# Metadata_total <- vroom("/home/cristia/Scrivania/BiomiX1.9/Metadata/Metadata_PRECISESADS.tsv", delim = "\t", col_names = TRUE)
# annotation <- vroom("/home/cristia/Scrivania/BiomiX1.9/Metabolomics/INPUT/ANNOTATION_PLASMA_SLE.tsv",delim="\t",col_names = TRUE)

# MS2_annotation = "YES"   #ADD TO INTERFACE
MS2_databases = c("HMDB","MONA","MASSBANK")  #ADD TO INTERFACE

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")


COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
Heatmap_genes <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])

LogFC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[2])

#MS2_DATABASES
Input_ms2<- COMMAND_ADVANCED[3,8]
Input_ms2 <- unlist(strsplit(as.character(Input_ms2), "/"))
Input_ms2<-as.numeric(substr(Input_ms2, 1, 1))
if(sum(is.na(Input_ms2)) > 0){
        MS2_databases<-MS2_databases[!is.na(Input_ms2)]
        Input_ms2<-Input_ms2[!is.na(Input_ms2)]
}
MS2_databases <-MS2_databases[order(Input_ms2)]


#loading annotation and matrix
directory2 <- paste(directory,"/Metabolomics/INPUT",sep="")
setwd(directory2)

print(i)
print(COMMAND$LABEL[i])

#LOADING ANNOTATION AND PEAKS
matrix <- vroom(COMMAND$DIRECTORIES[i],delim="\t",col_names = TRUE, comment = "#")

sam <- colnames(matrix)
pea <-matrix$ID
matrix <- t(matrix[,-1])
matrix <- add_column(as.data.frame(matrix), sam[-1], .after = 0) 
colnames(matrix) <-c("ID",pea)

ANNOTATION = COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL[1]


if(ANNOTATION == "MS1"){  
        input_ms1<-which(substr(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX,1,1) %in% i)
        if (length(input_ms1) != 0){
                if (file.exists(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES[input_ms1])){
                        annotation <- vroom(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES[input_ms1],delim="\t",col_names = TRUE)
                        if (ncol(annotation) == 3){
                                if (colnames(annotation)[3] == "RT_sec"){
                                        annotation$RT_min <- annotation$RT_sec/60
                                }
                                if (colnames(annotation)[3] == "RT_min"){
                                        annotation$RT_sec <- annotation$RT_min*60
                                }
                                
                        } 
                }
        }}


if(ANNOTATION == "MS2"){  
        input_ms1<-which(substr(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX,1,1) %in% i)
        if (length(input_ms1) != 0){
                if (file.exists(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2[input_ms1])){
                        annotation <- vroom(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2[input_ms1],delim="\t",col_names = TRUE)
                        if (ncol(annotation) == 3){
                                if (colnames(annotation)[3] == "RT_sec"){
                                        annotation$RT_min <- annotation$RT_sec/60
                                }
                                if (colnames(annotation)[3] == "RT_min"){
                                        annotation$RT_sec <- annotation$RT_min*60
                                }
                                
                        } 
                }
        }}

#create directory
dir.create(path = paste(directory,"/MOFA/INPUT/", "Metabolomics_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory, "/MOFA/INPUT/", "Metabolomics_", Cell_type,  "_",args[1],"_vs_", args[2], sep ="")

setwd(directory2)

#ANNOTATION <- "Annotated"
#UPLOAD SEUM URINE ANNOTATIONS FROM HMDB AND FILTERING (IF THE SAMPLE IS FROM SERUM OR URINE)

if (str_detect(COMMAND$LABEL[i], fixed("Serum", ignore_case=TRUE))| str_detect(COMMAND$LABEL[i], fixed("Plasma", ignore_case=TRUE))){
        if (file.exists("serum_metabolite_annotated.tsv") == TRUE){
                print("File available locally, using the local version")
                serum_metabolite <- vroom(paste(directory2,"serum_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                serum_metabolite <- vroom(url("https://hmdb.ca/metabolites.csv?action=index&blood=1&c=hmdb_id&controller=metabolites&d=up&detected=1&expected=1&filter=true&predicted=1&quantified=1&utf8=%E2%9C%93"),delim=",",col_names = TRUE )
                serum_metabolite <- as.data.frame(serum_metabolite)
                write.table(serum_metabolite,paste(directory2,"serum_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}

if (str_detect(COMMAND$LABEL[i], fixed("Urine", ignore_case=TRUE))){
        if (file.exists("urine_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                urine_metabolite <- vroom(paste(directory2,"urine_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                urine_metabolite <- vroom(url("https://hmdb.ca/metabolites.csv?action=index&c=hmdb_id&controller=metabolites&d=up&detected=1&expected=1&filter=true&predicted=1&quantified=1&urine=1&utf8=%E2%9C%93"),delim=",",col_names = TRUE )
                urine_metabolite <- as.data.frame(urine_metabolite)
                write.table(urine_metabolite,paste(directory2,"urine_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }
}

if (str_detect(COMMAND$LABEL[i], fixed("Saliva", ignore_case=TRUE))){
        if (file.exists("saliva_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                saliva_metabolite <- vroom(paste(directory2,"saliva_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                saliva_metabolite <- vroom(url("https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&saliva=1&filter=true"),delim=",",col_names = TRUE )
                saliva_metabolite <- as.data.frame(saliva_metabolite)
                write.table(saliva_metabolite,paste(directory2,"saliva_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}


if (str_detect(COMMAND$LABEL[i], fixed("Cerebrospinal Fluid", ignore_case=TRUE))){
        if (file.exists("cerebrospinal_fluid_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                CSF_metabolite <- vroom(paste(directory2,"cerebrospinal_fluid_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                CSF_metabolite <- vroom(url("https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&csf=1&filter=true"),delim=",",col_names = TRUE )
                CSF_metabolite <- as.data.frame(CSF_metabolite)
                write.table(CSF_metabolite,paste(directory2,"cerebrospinal_fluid_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}


if (str_detect(COMMAND$LABEL[i], fixed("Feces", ignore_case=TRUE))){
        if (file.exists("feces_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                feces_metabolite <- vroom(paste(directory2,"feces_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                feces_metabolite <- vroom(url("https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&feces=1&filter=true"),delim=",",col_names = TRUE )
                feces_metabolite <- as.data.frame(feces_metabolite)
                write.table(feces_metabolite,paste(directory2,"feces_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}

if (str_detect(COMMAND$LABEL[i], fixed("Sweat", ignore_case=TRUE))){
        if (file.exists("sweat_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                sweat_metabolite <- vroom(paste(directory2,"sweat_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                sweat_metabolite <- vroom(url("https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&sweat=1&filter=true"),delim=",",col_names = TRUE )
                sweat_metabolite <- as.data.frame(sweat_metabolite)
                write.table(sweat_metabolite,paste(directory2,"sweat_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}

if (str_detect(COMMAND$LABEL[i], fixed("Breast milk", ignore_case=TRUE))){
        if (file.exists("breast_milk_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                breast_milk_metabolite <- vroom(paste(directory2,"breast_milk_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                breast_milk_metabolite <- vroom(url("https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&breast_milk=1&filter=true"),delim=",",col_names = TRUE )
                breast_milk_metabolite <- as.data.frame(breast_milk_metabolite)
                write.table(breast_milk_metabolite,paste(directory2,"breast_milk_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}


if (str_detect(COMMAND$LABEL[i], fixed("Bile", ignore_case=TRUE))){
        if (file.exists("bile_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                bile_metabolite <- vroom(paste(directory2,"bile_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                bile_metabolite <- vroom(url("https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&bile=1&filter=true"),delim=",",col_names = TRUE )
                bile_metabolite <- as.data.frame(bile_metabolite)
                write.table(bile_metabolite,paste(directory2,"bile_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}

if (str_detect(COMMAND$LABEL[i], fixed("Amniotic Fluid", ignore_case=TRUE))){
        if (file.exists("AF_metabolite_annotated") == TRUE){
                print("File available locally, using the local version")
                AF_metabolite <- vroom(paste(directory2,"AF_metabolite_annotated.tsv",sep = "/"),delim="\t",col_names = TRUE)
        }else{
                print("File unavailable locally, downloading it from HMDB database")
                AF_metabolite <- vroom(url("https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&amniotic_fluid=1&filter=true"),delim=",",col_names = TRUE )
                AF_metabolite <- as.data.frame(AF_metabolite)
                write.table(AF_metabolite,paste(directory2,"AF_metabolite_annotated.tsv", sep = "/"),quote = FALSE, row.names = F, sep = "\t")
        }}

#CHECK OF THE CTRL AND CONDITION LABEL (ex SLE or CTRL)

Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)
# Identifier<-matrix$ID
# matrix<-matrix[,c(-1,-2)]
# matrix<-as.data.frame(t(matrix))
# colnames(matrix) <- Identifier
# matrix$ID <- row.names(matrix)
# matrix <- matrix[, c(ncol(matrix), 1:(ncol(matrix)-1))]
# rownames(matrix) = NULL
# write.table(x=matrix , file= "Urine.tsv" ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)


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

Metadata <- Metadata_Bcell
Metadata_individual=NULL
Metadata_reads=NULL
Metadata_Bcell=NULL
Metadata_tot = NULL
Metadata_total = NULL

print("Same samples in matrix and metadata?")
print(all(Metadata$ID == matrix[,1]))
outs <- Metadata$CONDITION  == args[2]|Metadata$CONDITION  == args[1]
Metadata<- Metadata[outs,]
matrix<- matrix[outs,]


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

# 
# 
# summary(as.factor(Metadata$CONDITION))
# matrix <- t(matrix)
# colnames(matrix) <- Identifier
# matrix<-as.data.frame(matrix)
# matrix$ID <- rownames(matrix)
# matrix$CONDITION <- Metadata$CONDITION
# samples <- matrix$ID
# matrix <- matrix[, c(ncol(matrix) -1,ncol(matrix), 1:(ncol(matrix)-2))]

###
matrix<- as.data.frame(matrix)
rownames(matrix) <- matrix$ID
###


#matrixs <- matrix %>% filter(CONDITION == args[2]| CONDITION ==args[1])
matrixs <- matrix 
matrixs <- add_column(matrixs, Metadata$CONDITION, .after = 1)
colnames(matrixs)[2] <- "CONDITION"


#SAVING THE INPUT FOR THE MOFA ANALYSIS
write.table(matrixs,paste(directory2,"/Metabolomics_",Cell_type, "_MOFA.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")



# #### FUNCTION DEFINITION FOR CONTROL AND TESTED SAMPLES----


DMS_REFERENCE <- function(Condition1){
        print(Condition1)
        CON <- matrix[Metadata$CONDITION == Condition1,] #Example CONTROL vs SLE patients
        samp_CON<-row.names(CON) #RUN IT ALWAYS + ONE DISEASE
        CON <-CON[,c(-1,-2)] #RUN IT ALWAYS + ONE DISEASE
        CON <-apply(CON,2,as.numeric) #RUN IT ALWAYS + ONE DISEASE
        
        y=NULL
        for (i in seq(1:length(colnames(CON)))) {
                res <- (shapiro.test(CON[,i]))
                y <-append(y, res[["p.value"]])
        }
        
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
        SLE <- matrix[Metadata$CONDITION == Condition2,] #RUN LINE 58 and 59 for CTRL vs SLE
        samp_SLE<-row.names(SLE) #RUN LINE 80 and 81 for CTRL vs SLE
        SLE <-SLE[,c(-1,-2)] #RUN LINE 109 and 110 for CTRL vs SLE
        SLE <-apply(SLE,2,as.numeric) #RUN LINE 130 and 131 for CTRL vs SLE
        
        y=NULL
        for (i in seq(1:length(colnames(SLE)))) {
                res <- (shapiro.test(SLE[,i]))
                y <-append(y, res[["p.value"]])
        }
        
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
        SLE$p_val <- pval
        SLE$log2FC <- fold
        
        return(SLE)
}






# #### DMS DIFFERENTIAL METABOLITE SIGNALS ANALYSIS----

CON <-DMS_REFERENCE(args[2])        
TEST <- DMS_TESTED(args[1])        
TEST$NAME <- row.names(TEST)
gc()

#CHECK IF THE PEAKS ARE ALREADY ANNOTATED OR NOT, IF NOT IT WILL COMPARE THE PEAK
#SIGNALS DIRECTLY, OTHERWISE IT WILL START THE ANNOTATION AT FIRST.

if(ANNOTATION == "Annotated"){
        
        dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        directory2 <- paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="")
        
        setwd(directory2)
        
        total <- TEST
        total$padj = p.adjust(total$p_val, method = "fdr")
        totalshOK <- total %>% filter(padj < padju)
        total <- total %>% arrange(padj)
        x <- colnames(total) %in% matrix$ID
        total <- total[,!x]
        write.table(total,paste(directory2,"/",Cell_type,"_",args[1],"_vs_",args[2],"_results.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
        
        
        total$padj<-as.numeric(total$padj)
        total$log2FC<-as.numeric(total$log2FC)
        
        total_min <- total[,c("NAME", "log2FC", "p_val", "padj")]
        total_min <- total_min %>% arrange(padj)
        write.table(total_min,paste(directory2,"/",Cell_type,"_",args[1],"_vs_",args[2],"_peak_statistics.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")

        
        ALTI <- subset(total, padj < padju)
        ALTI <-subset(ALTI, log2FC > LogFC)
        ALTI <- ALTI[order(ALTI$padj), ]
        #
        BASSI <- subset(total, padj < padju)
        BASSI <- subset(BASSI, log2FC < -LogFC)
        BASSI <- BASSI[order(BASSI$padj), ]
        x <- total$NAME %in% ALTI$NAME
        NO<-total[!x,]
        x <- NO$NAME %in% BASSI$NAME
        NO<-NO[!x,]
        
        
        write.table(x=ALTI$Name   , file= paste(Cell_type,"_",args[1],"_",args[2],"_UP.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x=BASSI$Name  , file= paste(Cell_type,"_",args[1],"_",args[2],"_DOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x=ALTI   , file= paste(Cell_type,"_",args[1],"_",args[2],"_METABUP.tsv",sep="") ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        write.table(x=BASSI  , file= paste(Cell_type,"_",args[1],"_",args[2],"_METABDOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE),c("ID","CONDITION")],
                    file= paste(args[1], "_",args[2],"_",Cell_type, "_metadata.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
        
        
        
        if (nrow(ALTI) == 0){
                SKIP_ALTI <- FALSE}else{SKIP_ALTI <- TRUE}
        if (nrow(BASSI) == 0){
                SKIP_BASSI <- FALSE}else{SKIP_BASSI <- TRUE}
        
        pdf(file=paste("plot_DMA_", args[1],"_",args[2],"_", Cell_type, ".pdf", sep=""))
        
        library(ggplot2)
        library(ggrepel)
        
        
        p <-ggplot() +
                ggtitle( paste("DGE", Cell_type,"_", args[1], "_vs_", args[2],sep="")) + theme(
                        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
                        axis.title.x = element_text(color="black", size=14, face="bold"),
                        axis.title.y = element_text(color="black", size=14, face="bold")) +
                ylab("-log10(padj)") + xlab("log2FC")
        if (SKIP_ALTI){
                p <- p + geom_point(data = ALTI, aes(x = log2FC, y = -log10(padj)), color = "red") +
                        geom_text_repel(data = ALTI[c(1:25), ], aes(x = log2FC , y = -log10(padj), label=NAME),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1, max.overlaps = 5)}
        if (SKIP_BASSI){
                p <- p + geom_point(data = BASSI, aes(x = log2FC, y = -log10(padj)), color = "blue") +
                        geom_text_repel(data = BASSI[1:25, ], aes(x = log2FC , y = -log10(padj), label=NAME),hjust=0, vjust=0,size=3, force = 20, force_pull = 1, max.overlaps = 5)}
        p <- p + geom_point(data = NO, aes(x = log2FC, y = -log10(padj)), color = "black")
        
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
                Heat <-apply(Heat,2,log10)
                
                library(ComplexHeatmap)
                setwd(directory2)
                pdf(file= paste("Heatmap_top_genes_", Cell_type,"_", args[2],"_vs_", args[1],".pdf",sep=""))
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
        
        
        
        
        ##### MetPath  #####
        library(metpath)
        library(tidyverse)
        library(dplyr)
        
        query_id<-as.data.frame(total[,c("NAME","log2FC","p_val")])
        
        
        
        if (nrow(query_id) != 0){
                
                if (COMMAND_ADVANCED[2,9] == "HMDB" ){
                        colnames(query_id)[1] <- "HMDB"
                        
                        dir.create(path = paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                        setwd(paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
                        
                        
                        pathway_class_HMDB = 
                                metpath::pathway_class(hmdb_pathway)
                        pathway_class_KEGG = 
                                metpath::pathway_class(kegg_hsa_pathway)
                        
                        
                        pdf(file=paste("Pathway_analysis_HMDB_", args[1],"_",args[2],"_", Cell_type,".pdf", sep=""))
                        
                        
                        gc()
                        remain_idx = which(unlist(pathway_class_HMDB) == "Metabolic;primary_pathway")
                        hmdb_pathway = hmdb_pathway[remain_idx]
                        hmdb_pathway
                        
                        
                        result = 
                                enrich_hmdb(query_id = unique(query_id$HMDB), 
                                            query_type = "compound", 
                                            id_type = "HMDB",
                                            pathway_database = hmdb_pathway,
                                            only_primary_pathway = TRUE,
                                            p_cutoff = 0.05, 
                                            p_adjust_method = "BH", 
                                            threads = as.numeric(COMMAND_ADVANCED[3,3]))
                        
                        result
                        
                        if (length(result) != 0){
                                
                                x<-enrich_bar_plot(
                                        object = result,
                                        x_axis = "p_value_adjust",
                                        cutoff = 1.1,
                                        top = 10
                                )
                                
                                print(x)
                                
                                x <-enrich_scatter_plot(object = result)
                                
                                print(x)
                                
                                write.table(result@result,paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], "/HMDB_table_results", sep =""),quote = FALSE, row.names = F, sep = "\t")
                                
                                dev.off()
                                
                                gc()
                        }
                        
                }
                
                
                if (COMMAND_ADVANCED[2,9] == "KEGG" ){
                        colnames(query_id)[1] <- "KEGG"
                        
                        dir.create(path = paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                        setwd(paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
                        
                        
                        pdf(file=paste("Pathway_analysis_KEGG_", args[1],"_",args[2],"_", Cell_type),".pdf", sep="")
                        
                        
                        head(pathway_class_KEGG)
                        remain_idx =
                                pathway_class_KEGG %>%
                                unlist() %>%
                                stringr::str_detect("Disease") %>%
                                `!`() %>%
                                which()
                        
                        remain_idx
                        
                        pathway_database =
                                kegg_hsa_pathway[remain_idx]
                        
                        pathway_database
                        
                        
                        result = 
                                enrich_kegg(query_id = unique(query_id$Kegg), 
                                            query_type = "compound", 
                                            id_type = "KEGG",
                                            pathway_database = pathway_database, 
                                            p_cutoff = 0.05, 
                                            p_adjust_method = "BH", 
                                            threads = as.numeric(COMMAND_ADVANCED[3,3]))
                        
                        
                        result
                        
                        if (length(result) != 0){
                                
                                x <-enrich_bar_plot(
                                        object = result,
                                        x_axis = "p_value_adjust",
                                        cutoff = 1.1,
                                        top = 10
                                )
                                print(x)
                                
                                x<-enrich_scatter_plot(object = result)
                                print(x)
                                
                                write.table(result@result,paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], "/KEGG_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
                                
                                dev.off()
                        }
                }
                
        }else{
                print("NO STATISTICAL SIGNIFICANT PEAK FOR METABOLITE PATHWAY ANALYSIS")
        }
        
        
        ##### MetaboAnalistR  #####
        
        
        #Enrichment analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep =""))
        
        if (COMMAND_ADVANCED[2,9] == "HMDB" ){
                write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "KEGG" ){
                write.table(x= query_id$Kegg[!is.na(query_id$Kegg )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "compound_name" ){
                write.table(x= query_id$Name[!is.na(query_id$Name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        
        
        #Pathway analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep =""))
        
        if (COMMAND_ADVANCED[2,9] == "HMDB" ){
                write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "KEGG" ){
                write.table(x= query_id$Kegg[!is.na(query_id$Kegg )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "compound_name" ){
                write.table(x= query_id$Name[!is.na(query_id$Name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        
        
        #Joint-Pathway analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep =""))
        
        query_id_select <- query_id[,c(1,2)]
        
        if (COMMAND_ADVANCED[2,9] == "HMDB" ){
                write.table(x= query_id_select, file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "KEGG" ){
                write.table(x= query_id_select , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "compound_name" ){
                write.table(x= query_id_select, file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        
        
        regeX <- paste("*",args[1],"_vs_",args[2], sep="")
        regeX2 <- paste(args[1],"_vs_",args[2], sep="")
        files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)
        
        if (length(files) != 0){
                for (fil in 1:length(files)){
                        print(fil)
                        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
                        LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
                        LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
                        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
                        
                        Sources<-str_split(files2[1], "_")[[1]][1]
                        print(Sources)
                        
                        dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                        setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep =""))
                        
                        write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
                        
                }}
        
        
        #Network analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep =""))
        
        query_id_select <- query_id[,c(1,2)]
        
        if (COMMAND_ADVANCED[2,9] == "HMDB" ){
                write.table(x= query_id_select, file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "KEGG" ){
                write.table(x= query_id_select , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],sep="") ,".tsv" ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        if (COMMAND_ADVANCED[2,9] == "compound_name" ){
                write.table(x= query_id_select, file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
        
        
        regeX <- paste("*",args[1],"_vs_",args[2], sep="")
        regeX2 <- paste(args[1],"_vs_",args[2], sep="")
        files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)
        
        if (length(files) != 0){
                for (fil in 1:length(files)){
                        print(fil)
                        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
                        LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
                        LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
                        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
                        
                        Sources<-str_split(files2[1], "_")[[1]][1]
                        print(Sources)
                        
                        dir.create(path = paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                        setwd(paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep =""))
                        
                        write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
                        
                }
        }
        
} else {
        
        
        
        if(ANNOTATION == "MS2"){  
                
                
                #### MS/MS FILTERING ####
                library(metid)
                
                directory3 <- paste(directory,"/Metabolomics",sep="")
                setwd(directory3)
                
                #path <- file.path(".", "INPUT")
                path <- as.character(COMMAND_ADVANCED[1,19])
                rt_threshold = as.numeric(COMMAND_ADVANCED[2,8])
                
                #LOADING MS2 DATABASES
                param = list()
                if (sum(MS2_databases == "HMDB") == 1){
                        HMDA <- load(paste(directory,"/MOFA/x_BiomiX_DATABASE/hmdb_database0.0.3.rda",sep=""))
                        param1<-identify_metabolites_params(
                                ms1.ms2.match.mz.tol = as.numeric(COMMAND_ADVANCED[1,8]),
                                #rt.match.tol = rt_threshold,
                                ms1.ms2.match.rt.tol = rt_threshold,
                                polarity = as.character(COMMAND_ADVANCED[1,24]),
                                ce = "all",
                                column = as.character(COMMAND_ADVANCED[3,23]),
                                total.score.tol = 0.5,
                                candidate.num = 3,
                                threads = as.numeric(COMMAND_ADVANCED[3,3]),
                                database = hmdb_database0.0.3)
                        param<-append(param,param1)
                }
                if (sum(MS2_databases == "MASSBANK") == 1){
                        MASSBANK <- load(paste(directory,"/MOFA/x_BiomiX_DATABASE/massbank_database0.0.3.rda",sep=""))
                        param1<-identify_metabolites_params(
                                ms1.ms2.match.mz.tol = as.numeric(COMMAND_ADVANCED[1,8]),
                                #rt.match.tol = rt_threshold,
                                ms1.ms2.match.rt.tol = rt_threshold,
                                polarity = as.character(COMMAND_ADVANCED[1,24]),
                                ce = "all",
                                column = as.character(COMMAND_ADVANCED[3,23]),
                                total.score.tol = 0.5,
                                candidate.num = 3,
                                threads = as.numeric(COMMAND_ADVANCED[3,3]),
                                database = massbank_database0.0.3)
                        param<-append(param,param1)
                }
                if (sum(MS2_databases == "MONA") == 1){
                        MONA <- load(paste(directory,"/MOFA/x_BiomiX_DATABASE/mona_database0.0.3.rda",sep=""))
                        param1<-identify_metabolites_params(
                                ms1.ms2.match.mz.tol = as.numeric(COMMAND_ADVANCED[1,8]),
                                #rt.match.tol = rt_threshold,
                                ms1.ms2.match.rt.tol = rt_threshold,
                                polarity = as.character(COMMAND_ADVANCED[1,24]),
                                ce = "all",
                                column = as.character(COMMAND_ADVANCED[3,23]),
                                total.score.tol = 0.5,
                                candidate.num = 3,
                                threads = as.numeric(COMMAND_ADVANCED[3,3]),
                                database = mona_database0.0.3)
                        param<-append(param,param1)
                }
                
                param2 <- param
                
                for (it in 1:length(param)){
                        if( MS2_databases[it] == toupper(param[[1]]$database@database.info$Source)){
                                param2[[it]] <- param[[1]]
                        }else if(MS2_databases[it] == toupper(param[[2]]$database@database.info$Source)){
                                param2[[it]] <- param[[2]]
                        }else{
                                param2[[it]] <- param[[3]]  
                        }}
                
                
                param <- param2
                
                gc()
                setwd(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1])
                annotation2 <- annotation[,-4]
                colnames(annotation2) <- c("\"name\"","\"mz\"","\"rt\"")
                write.table(annotation2,"TEMP.csv" ,quote = FALSE, row.names = F, sep = ",")
                
                #FIND THE LIST OF mzML file belonging to each matrix and add all of them in a vector
                files <- grep("MS2*",list.files(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1]),value=TRUE)
                files_MS2<-grep(paste("MS2_",COMMAND$LABEL[i],"*",sep=""),files,value=TRUE,ignore.case =TRUE)
                
                Sys.time()
                setwd(directory3)
                
                
                
                
                annotate_result5 <- 
                        identify_metabolite_all(ms1.data = "TEMP.csv", 
                                                ms2.data = files_MS2, #USING ONLY ONE FILE
                                                parameter.list = param,
                                                path = path)
                
                Sys.time()
                
                setwd(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1])
                
                if (file.exists("TEMP.csv")) {
                        file.remove("TEMP.csv")
                }
                
                ## RESTORE ####
                if (dir.exists("intermediate_data")) {
                        unlink("intermediate_data", recursive = TRUE)
                }
                
                ## RESTORE ####
                if (dir.exists("Result")) {
                        unlink("Result", recursive = TRUE)
                }
                
                
                
                
                
                #This part of the script is made to select all the peaks found in all the
                #databases. The final file (annot) will be used to replace the metabolites
                #with a level of annotation of 3 or 4.
                
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                directory2 <- paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="")
                
                
                setwd(directory2)
                
                gc()        
                unione <- list()
                unione_all <- list()
                iter=0
                for (dat in 1:length(param)){
                        iter = iter + 1
                        xx<-names(annotate_result5[[dat]]@identification.result)
                        saved<- annotate_result5[[dat]]@match.result[annotate_result5[[dat]]@match.result$MS2.spectra.name %in% xx,]
                        
                        y=NULL
                        
                        for (pe in seq(1:length(rownames(saved)))) {
                                
                                col_names<-colnames(annotate_result5[[dat]]@identification.result[[pe]])
                                
                                for(c in seq(1:nrow(annotate_result5[[dat]]@identification.result[[pe]]))){
                                        
                                        line<-annotate_result5[[dat]]@identification.result[[pe]][c,]
                                        x<-names(annotate_result5[[dat]]@identification.result[pe])
                                        x<-saved$MS2.spectra.name %in% x
                                        peak_number <- saved$MS1.peak.name[x]
                                        t<-unlist(c(peak_number, line))
                                        y <-rbind(y, t)
                                }
                        }
                        
                        if(iter==1){
                                unione <- as.data.frame(y)
                                colnames(unione)[1] <- "peak_number"
                                #write.table(unione,paste(directory2,"/",Cell_type,"_", "MS2_",names(annotate_result5)[dat],"_results", sep = ""),quote = FALSE, row.names = F, sep = "\t")
                                unione_all <- rbind(unione_all,unione)
                        }else{
                                unione2 <- as.data.frame(y)
                                colnames(unione2)[1] <- "peak_number"
                                #write.table(unione2,paste(directory2,"/",Cell_type,"_", "MS2_",names(annotate_result5)[dat],"_results", sep = ""),quote = FALSE, row.names = F, sep = "\t")
                                unione<-rbind(unione2[which(!unione2$peak_number %in% unione$peak_number),],unione)
                                unione_all <- rbind(unione_all,unione2)
                                
                        }
                        
                }
                
                unione_all
                unione_all <- unione_all[!duplicated(unione_all$Compound.name) & !duplicated(unione_all$Lab.ID ), ]
                unione_all$DATABASE <- "NA"
                unione_all$DATABASE[grep("*HMDB*",unione_all$Lab.ID )] <- "HMDB"
                unione_all$DATABASE[grep("MassBank*",unione_all$Lab.ID )] <- "MassBank"
                unione_all$DATABASE[grep("MONA*",unione_all$Lab.ID )] <- "MONA"
                
                unione <- unione[!duplicated(unione$Compound.name) & !duplicated(unione$Lab.ID ), ]  #?
                annot<-unione
                unione$DATABASE <- "FINAL_SELECTED"
                unione_all <- rbind(unione_all,unione)
                numbe<-as.numeric(gsub("peak", "", unione_all$peak_number))
                unione_all<-unione_all[order(numbe),]
                write.table(unione_all,paste(directory2,"/",Cell_type,"_", "MS2_results.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
                
                
                
                #ADD PART WHERE THERE ARE THE PLOTS FOR EACH PEAK2 IDENTIFIED
                dir.create(path = paste(directory2,"/MS2_SPECTRA", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                setwd(paste(directory2,"/MS2_SPECTRA", sep =""))
                
                if (sum(MS2_databases == "HMDB") == 1){
                        Y<- grep("HMDB*", names(annotate_result5), ignore.case=TRUE)
                        HMDA <- load(paste(directory,"/MOFA/x_BiomiX_DATABASE/hmdb_database0.0.3.rda",sep=""))
                        peak_to_show<-names(annotate_result5[[Y]]@identification.result)
                        matching<-annotate_result5[[Y]]@match.result
                        show<-matching[matching$MS2.spectra.name %in% peak_to_show,]
                        ms2.plot1 <- ms2plot(object = annotate_result5[[Y]],
                                             database = hmdb_database0.0.3,
                                             which.peak = "all")
                        file.rename("ms2_match_plot", "ms2_match_plot_HMDB")}
                
                if (sum(MS2_databases == "MASSBANK") == 1){
                        Y<- grep("MASSBANK*", names(annotate_result5),ignore.case=TRUE)
                        HMDA <- load(paste(directory,"/MOFA/x_BiomiX_DATABASE/massbank_database0.0.3.rda",sep=""))
                        peak_to_show<-names(annotate_result5[[Y]]@identification.result)
                        matching<-annotate_result5[[Y]]@match.result
                        show<-matching[matching$MS2.spectra.name %in% peak_to_show,]
                        ms2.plot1 <- ms2plot(object = annotate_result5[[Y]],
                                             database = massbank_database0.0.3,
                                             which.peak = "all")
                        file.rename("ms2_match_plot", "ms2_match_plot_MASSBANK")}
                
                if (sum(MS2_databases == "MONA") == 1){
                        Y<- grep("MONA*", names(annotate_result5),ignore.case=TRUE)
                        HMDA <- load(paste(directory,"/MOFA/x_BiomiX_DATABASE/mona_database0.0.3.rda",sep=""))
                        peak_to_show<-names(annotate_result5[[Y]]@identification.result)
                        matching<-annotate_result5[[Y]]@match.result
                        show<-matching[matching$MS2.spectra.name %in% peak_to_show,]
                        ms2.plot1 <- ms2plot(object = annotate_result5[[Y]],
                                             database = mona_database0.0.3,
                                             which.peak = "all")
                        file.rename("ms2_match_plot", "ms2_match_plot_MONA")}
                
                
                
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                directory2 <- paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="")
                
                
                setwd(directory2)      
        }else{
                
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                directory2 <- paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="")
                
                
                setwd(directory2)        
        }
        
        gc()
        
        total <- merge(TEST,annotation,by.x = "NAME", by.y = "name",all.x = TRUE)
        
        total$padj = p.adjust(total$p_val, method = "fdr") 
        
        total_min <- total[,c("NAME", "log2FC", "p_val", "padj")]
        total_min <- total_min %>% arrange(padj)
        write.table(total_min,paste(directory2,"/",Cell_type,"_",args[1],"_vs_",args[2],"_peak_statistics.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
        
        
        
        library(stringr)
        
        
        #Retrieve the annotation using the advanced research in ceu database
        
        if(ANNOTATION == "MS1"){ 
                adduct_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2[1]
                mode_ions <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1[1]
                adduct_list_2 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2[2]
                tolerance_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1[3]
                list_database_ms1 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2[3]}
        
        if(ANNOTATION == "MS2"){ 
                adduct_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2[1]
                mode_ions <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3[1]
                adduct_list_2 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2[2]
                tolerance_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4[2]
                list_database_ms1 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4[1]}
        
        #==============================================================================
        if(mode_ions == "positive" | mode_ions == "negative"){
                
                Positive_adduct <- as.array(adduct_list)
                Positive_adduct <- strsplit(Positive_adduct, "/")
                Positive_adduct <- Positive_adduct[[1]]
                
                Negative_adduct <- as.array(adduct_list_2)
                Negative_adduct <- strsplit(Negative_adduct, "/")
                Negative_adduct <- Negative_adduct[[1]]
                
                if (sum(is.na(Negative_adduct)) == 0){
                        ducts <- Positive_adduct %>% str_c(collapse = ",")
                        ducts <- paste("[",ducts,"]",sep = "")
                } else if (sum(is.na(Positive_adduct)) == 0){
                        ducts <- Positive_adduct  %>% str_c(collapse = ",") 
                        ducts <- paste("[",ducts,"]",sep = "")
                } else{
                        #ducts <- c(Positive_adduct,Negative_adduct)
                        print("ERROR YOU SELECTED POSITIVE AND NEGATIVE ADDUCTS AT THE SAME TIME!!")
                }
                
                masses_mode <- as.array("mz")
                Ion_mode <- as.array(mode_ions)
                
        } else if(mode_ions == "neutral"){
                Ion_mode <- as.array("neutral")
                masses_mode <- as.array("neutral")
                ducts <- "[\"all\"]"
        }
        
        tolerance <- as.array(tolerance_list)
        
        
        databases <- as.array(list_database_ms1)
        databases <- strsplit(databases, "/")
        databases <- databases[[1]]
        databases <- databases %>% str_c('"', ., '"') %>% str_c(collapse = ",") 
        databases<-paste("[",databases,"]",sep = "")
        
        
        
        MASS <- as.array(round(as.numeric(total$`m/z`),4))
        RT <- as.array(total$RT_min)
        
        cat(MASS, fill = getOption("width"), sep = ",")
        xx<-paste(MASS, collapse = " ")
        xx<-gsub(" ", ",", xx, fixed=TRUE)
        MASS<-paste("[",xx,"]",sep = "")
        
        cat(RT, fill = getOption("width"), sep = ",")
        xx<-paste(RT, collapse = " ")
        xx<-gsub(" ", ",", xx, fixed=TRUE)
        RT<-paste("[",xx,"]",sep = "")
        
        
        #THIS IF ELSE ALLOWS TO UPLOAD A PREVIOUS ANNOTATION MADE ON THE SAME DATASET.
        #IF THIS ANNOTATION IS NOT AVAILABLE OR IF IT IS THE FIRST ANALYSIS WILL AS TO
        #CEU MASS MEDIATOR TO DO IT.
        if (file.exists(paste("Annotation_metabolites_", Cell_type,"_", args[1], "_vs_", args[2],".tsv", sep = "")) == TRUE){
                print("File available locally, using the local version")
                advanced_batch_df <- vroom(paste(directory2,"/Annotation_metabolites_",Cell_type,"_", args[1], "_vs_", args[2],".tsv", sep = ""),delim="\t",col_names = TRUE)
        }else{
                
                
                
                
                library(cmmr)
                
                #Allow to user to change the option based on their request  
                
                #  '["hmdb","metlin", "kegg", "lipidmaps"]'
                #databases3 <-"[\"hmdb\",\"lipidmaps\",\"metlin\",\"kegg\"]"
                # '["M+H","M+Na","M+NH4","M+H-H2O"]'
                #database3 <-"[M+H,M+Na,M+NH4,M+H-H2O]"
                
                advanced_batch_df <- advanced_batch_search(
                        cmm_url             = paste0(
                                'http://ceumass.eps.uspceu.es/mediator/api/v3/',
                                'advancedbatch'),
                        chemical_alphabet   = 'all',
                        modifiers_type      = 'none',
                        metabolites_type    = 'all-except-peptides',
                        databases           = databases,
                        masses_mode         = masses_mode,
                        ion_mode            = Ion_mode,
                        adducts             = ducts,
                        deuterium           = 'false',
                        tolerance           = tolerance,
                        tolerance_mode      = "ppm",
                        masses              = MASS,
                        all_masses          = '[]',
                        retention_times     = RT,
                        all_retention_times = '[]'
                )
                
                #NOTE ADD THE OPTIONS WHEN CEU MASS MEDIATOR WILL WORK
                
                
                
                
                
                print("File unavailable locally, generating the annotation from Ceu Mass Mediator database")
                write.table(advanced_batch_df,paste(directory2,"/Annotation_metabolites_",Cell_type,"_", args[1], "_vs_", args[2] ,".tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
                
        }
        
        
        #save the annotation to avoid to run again the annotation search
        
        #===================================UPLOAD SAVED FILE==============================================
        
        head(advanced_batch_df)
        str(advanced_batch_df)
        gc()
        
        total$`m/z` %in% advanced_batch_df$Experimental.mass
        
        #Match annotation file and peaks file (total) using the m/z to match
        
        advanced_batch_df$Experimental.mass <- as.character(round(advanced_batch_df$Experimental.mass,4))
        total$`m/z` <- as.character(round(total$`m/z`,4))
        total$`m/z`<-as.character(total$`m/z`)
        
        
        
        tat<-merge(total,advanced_batch_df, by.x="m/z", by.y="Experimental.mass", all.x=TRUE)
        total$`m/z` %in% advanced_batch_df$Experimental.mass
        
        
        x <- colnames(tat) %in% matrix$ID
        tat <- tat[,!x]
        #merge the metabolites found with the table containing the results
        
        
        print(i)
        print(COMMAND$LABEL[i])
        filtering= 0
        
        
        #### FILTER OF TISSUES SPECIFIC METABOLITES IN HMDB #####
        #IF THE ANNOTATION FOR SERUM OR URINE IS AVAILABLE IT IS USED
        #TO FILTER THE RESULTS OBTAINED USING CEU MASS MEDIATOR
        if (str_detect(COMMAND$LABEL[i], fixed("Serum", ignore_case=TRUE)) | str_detect(COMMAND$LABEL[i], fixed("Plasma", ignore_case=TRUE))){
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% serum_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,serum_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Urine", ignore_case=TRUE))) {
                
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% urine_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,urine_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Saliva", ignore_case=TRUE))) {
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% saliva_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,saliva_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Cerebrospinal Fluid", ignore_case=TRUE))) {
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% CSF_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,CSF_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Feces", ignore_case=TRUE))) {
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% Feces_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,Feces_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
                
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Sweat", ignore_case=TRUE))) {
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% sweat_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,sweat_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Brest Milk", ignore_case=TRUE))) {
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% breast_milk_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,breast_milk_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Bile", ignore_case=TRUE))) {
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% bile_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,bile_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
                
                
        } else if (str_detect(COMMAND$LABEL[i], fixed("Amniotic Fluid", ignore_case=TRUE))) {
                
                #IF the peak annotation in ms1 doesn't have a match in the tissue selected
                #its annotation is eliminated and degradated to level4.
                
                to_eliminate2 =NULL
                for (pe in unique(tat$NAME)){
                        sat<-which(tat$NAME == pe)
                        if(sum(tat$HMDB[sat] %in% AF_metabolite$HMDB_ID) == 0){
                                tat[sat,9:ncol(tat)] <- NA
                                to_eliminate2<-append(to_eliminate2,sat[-1])
                                
                        }
                        
                }
                tat <- tat[-to_eliminate2,]
                
                total <-merge(tat,AF_metabolite, by.x="HMDB", by.y="HMDB_ID",all.x=TRUE )
                exit<-!is.na(total$Identifier) & is.na(total$SMILES)
                total <- total[!exit,]
                filtering <- "YES"
                gc()
                
                
        } else {
                print("NO BIOSPECIMENS filtering")
                total <- tat
                filtering <- "NO"
                colnames(total)[which(colnames(total) == "NAME")] <- "NAME.x"
        }
        
        
        totalshOK <- total %>% filter(padj <padju)
        total <- total %>% arrange(padj)
        print(colnames(total))
        
        # x <- colnames(total) %in% matrix$ID
        # total <- total[,!x]
        
        total$MS2_annot <- "NA"
        total$annot_level <- "Level_3"
        
        to_eliminate=NULL
        to_eliminate2=NULL
        
        if(ANNOTATION == "MS2"){  
                
                for(yy in 1:length(unique(annot$peak_number)) ){
                        t<-annot$peak_number %in% unique(annot$peak_number)[yy]
                        mix<-paste(annot$Compound.name[t], collapse = "/")
                        mix2<-paste(annot$Adduct[t], collapse = "/")
                        mix3<-paste(annot$mz.error[t], collapse = "/")
                        
                        if (sum(is.na(annot$HMDB.ID[t])) == 0){
                                x<-substr(total$HMDB,7,11) %in% substr(annot$HMDB.ID[t], 5,9)
                        }else{
                                x <- FALSE
                        }
                        if(sum(x) == 0){
                                print(paste(unique(annot$peak_number)[yy], "not_annotated",sep="_"))
                                z <-which(total$NAME.x %in% unique(annot$peak_number)[yy])
                                total[z,c(11:(ncol(total)-2))] <- "NA"
                                total[z[1],"Name"] <- mix
                                total[z[1],"Adduct"] <- mix2
                                total[z[1],"PPM.Error"] <- mix3
                                total[z[1],"MS2_annot"] <- mix
                                total[z[1],"annot_level"] <- "Level_2"
                                to_eliminate2<-append(to_eliminate2,z[-1])
                                
                        }else{
                                t<-total$HMDB[x]
                                t<-which(total$HMDB %in% t & total$NAME.x == unique(annot$peak_number)[yy])
                                
                                for(c in t){
                                        total$MS2_annot[c] <- mix
                                        total$annot_level[c]  <- "Level_2"
                                        
                                }
                                out<-which(total$NAME.x == unique(annot$peak_number)[yy] & total$annot_level == "Level_3")
                                to_eliminate <-append(to_eliminate,out)
                        }
                }
                
                if(length(to_eliminate) != 0){
                        total<-total[-to_eliminate,]
                }
                
                total$annot_level[is.na(total$Name)] <- "Level_4"
                
        }
        
        
        write.table(total,paste(directory2,"/",Cell_type,"_",args[1],"_vs_",args[2],"_results.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
        
        
        
        
        pdf(file=paste("plot_DMA_", args[1],"_",args[2],"_", Cell_type,".pdf",sep=""))
        
        library(ggplot2)
        library(ggrepel)
        
        #IF THE SERUM/URINE FILTERING HAS BEEN MADE IT WILL WILL DO THE PLOT AND 
        #SAVE THE UP AND DOWN REGULATED METABOLITES
        
        total$padj<-as.numeric(total$padj)
        total$log2FC<-as.numeric(total$log2FC)
        
        ALTI <- subset(total, padj < padju)
        ALTI <-subset(ALTI, log2FC > LogFC)
        ALTI <- ALTI[order(ALTI$padj), ]
        #
        BASSI <- subset(total, padj < padju)
        BASSI <- subset(BASSI, log2FC < -LogFC)
        BASSI <- BASSI[order(BASSI$padj), ]
        x <- total$NAME.x %in% ALTI$NAME.x
        NO<-total[!x,]
        x <- NO$NAME.x %in% BASSI$NAME.x
        NO<-NO[!x,]
        
        
        write.table(x=ALTI$Name   , file= paste(Cell_type,"_",args[1],"_",args[2],"_UP.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x=BASSI$Name  , file= paste(Cell_type,"_",args[1],"_",args[2],"_DOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x=ALTI   , file= paste(Cell_type,"_",args[1],"_",args[2],"_METABUP.tsv",sep="") ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        write.table(x=BASSI  , file= paste(Cell_type,"_",args[1],"_",args[2],"_METABDOWN.tsv",sep=""),sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        
        write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE),c("ID","CONDITION")],
                    file= paste(args[1], "_",args[2],"_",Cell_type, "_metadata.tsv", sep =""),sep = "\t", quote= FALSE, row.names = FALSE)
        
        
        
        ALTI <- ALTI %>% distinct(NAME.x, .keep_all = TRUE)
        BASSI <- BASSI %>% distinct(NAME.x, .keep_all = TRUE)
        NO <- NO %>% distinct(NAME.x, .keep_all = TRUE)
        
        
        if (nrow(ALTI) == 0){
                SKIP_ALTI <- FALSE}else{SKIP_ALTI <- TRUE}
        if (nrow(BASSI) == 0){
                SKIP_BASSI <- FALSE}else{SKIP_BASSI <- TRUE}
        
        
        p <-ggplot() +
                ggtitle( paste("DGE", Cell_type,"_", args[1], "_vs_", args[2],sep="")) + theme(
                        plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),
                        axis.title.x = element_text(color="black", size=14, face="bold"),
                        axis.title.y = element_text(color="black", size=14, face="bold")) +
                ylab("-log10(padj)") + xlab("log2FC")
        if (SKIP_ALTI){
                p <- p + geom_point(data = ALTI, aes(x = log2FC, y = -log10(padj)), color = "red") +
                        geom_text_repel(data = ALTI[c(1:25), ], aes(x = log2FC , y = -log10(padj), label=Name),hjust=0, vjust=0,size=3, arrow = NULL, force = 10, force_pull = 1, max.overlaps = 10)}
        if (SKIP_BASSI){
                p <- p + geom_point(data = BASSI, aes(x = log2FC, y = -log10(padj)), color = "blue") +
                        geom_text_repel(data = BASSI[1:25, ], aes(x = log2FC , y = -log10(padj), label=Name),hjust=0, vjust=0,size=3, force = 20, force_pull = 1, max.overlaps = 10)}
        
        p <- p + geom_point(data = NO, aes(x = log2FC, y = -log10(padj)), color = "black")
        
        p2 <- p + geom_vline(xintercept=c(-LogFC,LogFC), col="red") +
                geom_hline(yintercept=-log10(padju), col="red")
        
        
        print(p2)
        dev.off()
        
        
        ###  HEATMAP ###
        
        if (length(c(ALTI$NAME.x, BASSI$NAME.x)) > 2){
                
                #ADD HEATMAP CREATION HERE
                if (nrow(ALTI) > Heatmap_genes){
                        ALTI <- ALTI[1:Heatmap_genes,]}
                if (nrow(BASSI) > Heatmap_genes){
                        BASSI <- BASSI[1:Heatmap_genes,]}
                
                #HEATMAP INPUT ARE NORMALIZED
                
                Heat<-matrix[,colnames(matrix) %in% c(ALTI$NAME.x, BASSI$NAME.x),drop = FALSE]
                samples <- rownames(Heat)
                Heat <- as.data.frame(Heat)
                Heat <-apply(Heat,2,as.numeric)
                rownames(Heat) <-samples
                Heat<-t(Heat)
                Heat <-apply(Heat,2,log10)
                annot_vis<-rownames(Heat)
                pool<-rbind(ALTI,BASSI)
                
                SAVED <- NULL
                for (ix in annot_vis){
                        iu<- pool$NAME.x %in% ix 
                        Gene <- as.character(pool$Name[which(iu)])
                        SAVED<-append(Gene, SAVED)}
                SAVED <- as.character(SAVED)
                SAVED<-rev(SAVED)
                
                for (ix in 1:length(samples)){
                        if(!is.na(SAVED[ix])){
                                rownames(Heat)[ix] <- SAVED[ix]
                        }
                }
                
                
                library(ComplexHeatmap)
                setwd(directory2)
                pdf(file= paste("Heatmap_top_genes_", Cell_type,"_", args[2],"_vs_", args[1],".pdf",sep=""))
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
        
        ##### MetPath  #####
        library(metpath)
        library(tidyverse)
        library(dplyr)
        
        dir.create(path = paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
        
        saved<-total
        #Filtering of peak with one or max 2 annotation
        peak_repeated <-total %>% count(NAME.x)
        total<- total[!total$NAME.x %in% peak_repeated$NAME.x[peak_repeated$n > 2],]
        total <- total %>% arrange(NAME.x)
        
        #Eliminate level 3 peaks with level 2 annotation found
        peak_repeated<-total$NAME.x[total$annot_level == "Level_2"]
        outs<-total$NAME.x %in% peak_repeated & total$annot_level == "Level_3"
        total <-total[!outs,]
        
        #Selection of best annotation in peak with 2 annotation based on minimun ppm error
        iter<-which(duplicated(total$NAME.x))
        to_elim = NULL
        iterat=NULL
        for (u in iter){
                n1<-u
                s1<-total$PPM.Error[u]
                n2<-u-1
                s2<-total$PPM.Error[u-1]
                
                elim <- max(as.numeric(s1), as.numeric(s2))
                if(elim == as.numeric(s1)){
                        print = as.numeric(s1)
                        to_elim<-append(to_elim, n1)
                }
                if(elim == as.numeric(s2)){
                        to_elim<-append(to_elim, n2)
                }
                
        }
        
        if(length(to_elim) != 0){
                total <- total[-to_elim, ]}
        
        query_id <- total[which(total$p_val < 0.05),] #REPLACE IT WITH P.ADJ
        query_id_saved <- saved[which(saved$p_val < 0.05),]
        
        #SAVE PATHWAY ANALYSIS RESULTS
        if (nrow(query_id) != 0){
                
                dir.create(path = paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                setwd(paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
                
                
                pathway_class_HMDB = 
                        metpath::pathway_class(hmdb_pathway)
                pathway_class_KEGG = 
                        metpath::pathway_class(kegg_hsa_pathway)
                
                
                pdf(file=paste("Pathway_analysis_HMDB_", args[1],"_",args[2],"_", Cell_type,".pdf", sep=""))
                
                
                gc()
                remain_idx = which(unlist(pathway_class_HMDB) == "Metabolic;primary_pathway")
                hmdb_pathway = hmdb_pathway[remain_idx]
                hmdb_pathway
                
                
                result = 
                        enrich_hmdb(query_id = unique(query_id$HMDB), 
                                    query_type = "compound", 
                                    id_type = "HMDB",
                                    pathway_database = hmdb_pathway,
                                    only_primary_pathway = TRUE,
                                    p_cutoff = 0.05, 
                                    p_adjust_method = "BH", 
                                    threads = as.numeric(COMMAND_ADVANCED[3,3]))
                
                result
                
                
                if (length(result) != 0){
                        
                        x<-enrich_bar_plot(
                                object = result,
                                x_axis = "p_value_adjust",
                                cutoff = 1.1,
                                top = 10
                        )
                        
                        print(x)
                        
                        x <-enrich_scatter_plot(object = result)
                        
                        print(x)
                        
                        write.table(result@result,paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], "/HMDB_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
                        
                        dev.off()
                        
                }
                
                gc()
                dir.create(path = paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                setwd(paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
                
                
                pdf(file=paste("Pathway_analysis_KEGG_", args[1],"_",args[2],"_", Cell_type,".pdf", sep=""))
                
                
                head(pathway_class_KEGG)
                remain_idx =
                        pathway_class_KEGG %>%
                        unlist() %>%
                        stringr::str_detect("Disease") %>%
                        `!`() %>%
                        which()
                
                remain_idx
                
                pathway_database =
                        kegg_hsa_pathway[remain_idx]
                
                pathway_database
                
                
                result = 
                        enrich_kegg(query_id = unique(query_id$Kegg), 
                                    query_type = "compound", 
                                    id_type = "KEGG",
                                    pathway_database = pathway_database, 
                                    p_cutoff = 0.05, 
                                    p_adjust_method = "BH", 
                                    threads = as.numeric(COMMAND_ADVANCED[3,3]))
                
                
                result
                
                if (length(result) != 0){
                        
                        x <-enrich_bar_plot(
                                object = result,
                                x_axis = "p_value_adjust",
                                cutoff = 1.1,
                                top = 10
                        )
                        print(x)
                        
                        x<-enrich_scatter_plot(object = result)
                        print(x)
                        
                        write.table(result@result,paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], "/KEGG_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
                        
                        dev.off()
                        
                }
                
        }else{
                print("NO STATISTICAL SIGNIFICANT PEAK FOR METABOLITE PATHWAY ANALYSIS")
        }
        
        
        ##### MetaboAnalistR  #####
        
        
        #Enrichment analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep =""))
        
        numbe<-as.numeric(gsub("peak", "", query_id$NAME.x))
        query_id<-query_id[order(numbe),]
        
        write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv", sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id$Kegg[!is.na(query_id$Kegg )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id$Name[!is.na(query_id$Name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        
        numbe<-as.numeric(gsub("peak", "", query_id_saved$NAME.x))
        query_id_saved<-query_id_saved[order(numbe),]
        write.table(x= query_id_saved[,c(3,1,21,14)] , file= paste("HMDB_KEGG_Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],"_complete",".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        
        
        #Pathway analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep =""))
        
        write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id$Kegg[!is.na(query_id$Kegg )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id$Name[!is.na(query_id$Name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id_saved[,c(3,1,21,14)] , file= paste("HMDB_KEGG_Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],"_complete.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        
        
        #Joint-Pathway analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep =""))
        
        query_id_select <- query_id[,c("HMDB", "Kegg", "Name", "log2FC")]
        
        write.table(x= query_id_select[,c(1,4)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id_select[,c(2,4)] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id_select[,c(3,4)], file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        
        
        regeX <- paste("*",args[1],"_vs_",args[2], sep="")
        regeX2 <- paste(args[1],"_vs_",args[2], sep="")
        files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)
        
        if (length(files) != 0){
                for (fil in 1:length(files)){
                        print(fil)
                        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
                        LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
                        LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
                        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
                        
                        Sources<-str_split(files2[1], "_")[[1]][1]
                        print(Sources)
                        
                        dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                        setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep =""))
                        
                        write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
                        
                }}
        
        
        #Network analysis
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep =""))
        
        query_id_select <- query_id[,c("HMDB", "Kegg", "Name", "log2FC")]
        
        write.table(x= query_id_select[,c(1,4)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id_select[,c(2,4)] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        write.table(x= query_id_select[,c(3,4)], file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
        
        
        regeX <- paste("*",args[1],"_vs_",args[2], sep="")
        regeX2 <- paste(args[1],"_vs_",args[2], sep="")
        files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)
        
        if (length(files) != 0){
                for (fil in 1:length(files)){
                        print(fil)
                        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
                        LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
                        LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
                        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
                        
                        Sources<-str_split(files2[1], "_")[[1]][1]
                        print(Sources)
                        
                        dir.create(path = paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                        setwd(paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep =""))
                        
                        write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
                        
                }}
        
        
        
        #MSEA ANALYSIS
        
        dir.create(path = paste(directory2,"/MetaboAnalyst/", "Functional_analysis_MSEA", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/", "Functional_analysis_MSEA", sep =""))
        
        MSEA <- total[,c("m/z","p_val", "log2FC", "RT_min" )]
        colnames(MSEA) <- c("m.z","p.value", "t.score", "r.t" )
        
        
        write.table(x= MSEA, file= paste("MSEA_",args[1],"_vs_",args[2],".txt",sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
        gc()
        
}


advanced_batch_df = NULL
param=NULL
param1=NULL
annotation = NULL
annotation2 =NULL
massbank_database0.0.3 =NULL
mona_database0.0.3= NULL
pathway_class_KEGG =NULL
pathway_class_HMDB =NULL
hmdb_database0.0.3=NULL
tat= NULL
NO=NULL
serum_metabolite =NULL
annotate_result5=NULL
results=NULL
gc()


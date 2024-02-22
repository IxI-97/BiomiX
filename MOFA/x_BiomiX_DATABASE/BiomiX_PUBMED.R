library(vroom)
library(dplyr)
library(XML)
library(xml2)
library(stringr)
library(rentrez)


# MANUAL INPUT
# # #
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"SLE"
# args[2] <-"CTRL"
# args[3] <-"/home/cristia/Scrivania/BiomiX1.5"
# #
# directory <- args[3]

COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED",sep="/"), delim = "\t")

TYPE_OF_RESEARCH <- as.character(COMMAND_ADVANCED[1,10])   #Text Word  #"Title/Abstract"
retmaxi <- as.numeric(COMMAND_ADVANCED[2,10])  


# Function to retrieve keywords from PubMed given a PMID
get_keywords_from_pubmed <- function(pmid) {
        # Fetch the PubMed record in XML format
        record <- entrez_fetch(db = "pubmed", id = pmid, rettype = "xml")
        
        # Parse the XML record
        xml_doc <- xmlParse(record)
        
        # Extract the keyword elements from the XML
        keyword_nodes <- getNodeSet(xml_doc, "//KeywordList/Keyword")
        
        # Extract the keyword values from the nodes
        keywords <- sapply(keyword_nodes, xmlValue)
        
        return(keywords)
}



get_doi_from_pubmed <- function(pmid) {
        # Fetch the PubMed record in XML format
        record <- entrez_fetch(db = "pubmed", id = pmid, rettype = "xml")
        
        # Parse the XML record
        xml_doc <- xmlParse(record)
        
        # Extract the DOI from the XML
        doi <- xpathSApply(xml_doc, "//ArticleId[@IdType='doi']", xmlValue)
        
        return(doi)
}


# Function to count the occurrences of a list of words in a text string
count_word_occurrences <- function(text, word_list) {
        counts <- vector("numeric", length(word_list))
        for (i in seq_along(word_list)) {
                counts[i] <- sum(str_count(text, regex(word_list[i], ignore_case = TRUE)))
        }
        names(counts) <- word_list
        return(counts)
}

count_match_in_top_article <- function(nami) {
        
        iter= 0
        
        for (omik in (nami)){
                iter <- iter + 1
                print(omik)
                val<-grep(paste(omik,".*_pos_","*", sep=""), files)
                LIST_GENE<- vroom(paste(directory2,"/",files[val], sep=""), delim="\t")
                LIST_GENE<-LIST_GENE %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
                #IF there is only peak the query is removed, peak is removed from the total name
                LIST_GENE$features <-gsub("peak[0-9]*/", "", LIST_GENE$features)
                #Removal peak name
                sav<-grep("peak[0-9]*", LIST_GENE$features)
                if(length(sav) == 0){
                        print("selection...")}else{
                                LIST_GENE<- LIST_GENE[-sav,]}
                
                val<-grep(paste(omik,".*_neg_","*", sep=""), files)
                LIST_GENE2<-vroom(paste(directory2,"/",files[val], sep=""), delim="\t")
                LIST_GENE2<-LIST_GENE2 %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
                LIST_GENE2$features <-gsub("peak[0-9]*/", "", LIST_GENE2$features)
                #Removal peak name
                sav<-grep("peak[0-9]*", LIST_GENE2$features)
                if(length(sav) == 0){
                        print("selection...")}else{
                                LIST_GENE2<- LIST_GENE2[-sav,]}
                
                LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
                LIST_GENE3<-LIST_GENE3 %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
                
                #Retrieve gene name from cpg islands
                LIST_GENE3$features <-gsub("cg[0-9]*/", "", LIST_GENE3$features)
                sav<-grep("NA", LIST_GENE3$features)
                if(length(sav) != 0){
                LIST_GENE3<- LIST_GENE3[-sav,]}
                
                #TEST TO REMOVE COIL GENE (ONLY TO TEST)
                sav<-grep("COIL", LIST_GENE3$features)
                if(length(sav) != 0){
                        LIST_GENE3<- LIST_GENE3[-sav,]}
                
                if(length(LIST_GENE3$features) < 10){
                        x<-paste(LIST_GENE3$features[1:length(LIST_GENE3$features)-1], "OR ", sep= paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""))
                        xx<-c(x, paste(LIST_GENE3$features[length(LIST_GENE3$features)],paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""),sep=""))
                        columns<-LIST_GENE3$features[1:length(LIST_GENE3$features)]
                }else{
                x<-paste(LIST_GENE3$features[1:9], "OR ", sep= paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""))
                xx<-c(x, paste(LIST_GENE3$features[10],paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""),sep=""))
                columns<-LIST_GENE3$features[1:10]
                }
                
                xxx<-paste(xx, collapse="")
                xxx<-paste("(", xxx, ")", sep="")
                #print(xxx)
                
                if (iter > 1){
                        yyy<-paste(yyy, xxx, sep=" AND ")
                        columns_tot <-c(columns_tot, columns)
                }else{
                        yyy <- xxx
                        columns_tot <-columns
                }
        }
        
        print(yyy)
        esearch <- entrez_search(db = "pubmed", term = yyy, retmax = retmaxi, use_history = TRUE) #, use_history = TRUE
        ids <- esearch$ids
        
        if(length(ids) != 0){
                
                columns_tot<-unique(columns_tot[!is.na(columns_tot)])
                word_match<-matrix(data=NA,nrow=1,ncol=length(columns_tot) + 3)
                
                
                # Iterate over each article and check for a match with the words

                efetch <- entrez_fetch(db = "pubmed", web_history = esearch$web_history, rettype = "abstract", retmode = "xml")        
                parsed_xml <- read_xml(efetch)
                
                #Remove comments/errantum
                nodes_to_remove <- xml_find_all(parsed_xml, "//CommentsCorrectionsList/CommentsCorrections/PMID")
                xml_remove(nodes_to_remove)
                nodes <- xml_text(nodes_to_remove)
                nodes
                
                articles <- xml_find_all(parsed_xml, "//PubmedArticle")
                
                
                fu<-function(article) {
                        pmid_node <- xml_find_first(article, ".//PMID")
                        pmid_abstract <- xml_find_first(article, ".//Abstract")
                        pmid_doi <- xml_find_first(article, ".//ArticleId[@IdType='doi']")
                        pmid_key <- xml_find_all(article, ".//Keyword")
                        id <-xml_text(pmid_node)
                        doi <-xml_text(pmid_doi)
                        article_abstract<- xml_text(pmid_abstract)
                        article_abstract <- gsub("\\(", "", article_abstract)
                        article_abstract <- gsub("\\)", "", article_abstract)
                        text<-unlist(strsplit(article_abstract, " "))
                        #print(article_abstract)
                        
                        keywords<-xml_text(pmid_key)
                        keywords <- paste(keywords, collapse = "/")
                        #print(keywords)
                        
                        # Define the list of words you want to count
                        word_list <- columns_tot
                        # Count the occurrences of the words in the text string
                        occurrences <- count_word_occurrences(text, word_list)
                        occurrences <- c(occurrences, id, keywords, doi)
                        #print(occurrences)
                        #print(rbind(occurrences,word_match))
                        #word_match<-rbind(occurrences,word_match)
                        return(occurrences)
                }
                
                # Extract PMID for each article
                pmids <- sapply(articles,fu)        
                word_match<-as.data.frame(t(pmids))
                colnames(word_match) <- c(columns_tot,"Pubmed_ID","Keywords", "DOI")
                        
        }else{
                word_match = matrix()
}
                
                
                if (ncol(word_match) > 1) {
                if (nrow(word_match) > 1) {
                        
                        #SELETION BEST ARTICLE
                        df<-apply(word_match[,-c((ncol(word_match)-2):ncol(word_match)),drop=FALSE],2, as.numeric) > 0
                        somma<-rowSums(apply(word_match[,-c((ncol(word_match)-2):ncol(word_match)),drop=FALSE],2, as.numeric))
                        row_sums <- rowSums(df)
                        word_match <- as.data.frame(word_match)
                        word_match$tot_match <- somma
                        word_match$tot_word <- row_sums
                        print(word_match)
                        
                        #WITH TOTAL NUMBER OF MATCH SELECTION BASED ON THE TOTAL NUMBER OF MATCH
                        top_3_rows <- order(-word_match$tot_word, -word_match$tot_match)[1:5]
                        word_match <- word_match[top_3_rows,]
                        word_match <- word_match[!is.na(word_match$tot_word),] 
                        #print(word_match)
                        return(word_match)}else if(nrow(word_match) == 1) {
                                print("one_results")
                                return(word_match)}else{
                                        print("no_results")
                                }
}else{
        print("no_results")  
}
        
}
        
        




# # #MANUAL INPUT
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"SLE"
# args[2] <-"CTRL"
# args[3] <- "/home/cristia/Scrivania/BiomiX1.4"
# #
# directory <-args[3]


files_tot <- grep(paste("MOFA", "_", args[2] ,"_vs_", args[1],"*", sep=""),list.files(paste(directory,"/MOFA/OUTPUT/",sep="")),value=TRUE)
directory <- args[[3]]

for (fil in files_tot){


MART <- vroom(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

COMMAND <- vroom(paste(directory,"COMMANDS",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA",sep="/"), delim = "\t")

dir.create(path = paste(directory,"/MOFA/OUTPUT/",fil, "/Factors_Associated_Articles" ,sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory,"/MOFA/OUTPUT/",fil,sep="") 
directory3 <- paste(directory,"/MOFA/OUTPUT/",fil, "/Factors_Associated_Articles" ,sep="")
setwd(directory2)


#detection total significant factors
files <- grep("*factor_",list.files(directory2),value=TRUE)
files<- files[grep("\\Metabolomics|\\RNAseq|\\Methylomics", files)]
factors<-strsplit(files, "_")
factors<-unlist(factors)
factors<-unique(as.numeric(factors[!is.na(as.numeric(factors))]))

for (numb in factors){
        
#####BLOCK RESEARCH QUERY ALL OMICS TOGETHER ####

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Metabolomics|\\RNAseq|\\Methylomics", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

iter= 0

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)

line="#ALL OMICS (TRASCRIPTOMICS-METABOLOMICS-METHYLOMICS)#"
write(line,file=paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)
if (length(uu) == 0){
write("NULL",file=paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)} else{
        write.table(uu,file=paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
}


#####BLOCK RESEARCH QUERY OMICS PAIRS TOGETHER ####
#TRASCRIPTOMICS-METABOLOMICS

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Metabolomics|\\RNAseq", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

iter= 0

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)

line="#PAIRS OMICS (TRASCRIPTOMICS-METABOLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)
if (length(uu) == 0){
        write("NULL",paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
}

#TRASCRIPTOMICS-METHYLOMICS

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Methylomics|\\RNAseq", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

iter= 0

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)

line="#PAIRS OMICS (TRASCRIPTOMICS-METHYLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)
if (length(uu) == 0){
        write("NULL",file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
}




#METABOLOMICS-METHYLOMICS

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Methylomics|\\Metabolomics", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

iter= 0

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)

line="#PAIRS OMICS (METABOLOMICS-METHYLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)
if (length(uu) == 0){
        write("NULL",file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
}

#BLOCK SINGLE OMICS ANALYSIS
#TRANSCRIPTOMICS

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\RNAseq", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

iter= 0

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)

line="#SINGLE OMICS (TRASCRIPTOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)
if (length(uu) == 0){
        write("NULL",file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")
        }
}

#METABOLOMICS

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Metabolomics", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

iter= 0

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)


line="#SINGLE OMICS (METABOLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)
if (length(uu) == 0){
        write("NULL",file=paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)} else{
                write.table(uu,file=paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
}

#BLOCK SINGLE OMICS ANALYSIS
#METHYLOMICS

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Methylomics", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)


line="#SINGLE OMICS (METHYLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)
if (length(uu) == 0){
        write("NULL",file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }

}
}

}


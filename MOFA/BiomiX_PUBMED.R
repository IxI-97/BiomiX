library(vroom)
library(dplyr)
library(XML)
library(xml2)
library(stringr)
library(rentrez)
library(litsearchr)


# MANUAL INPUT
# # #
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"C4"
# args[2] <-"Control"
# args[3] <-"/home/henry/Desktop/BiomiX2.2"
# #
# directory <- args[3]

COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
ANNOTATION_LIPIDS <- vroom(paste(directory,"MOFA/x_BiomiX_DATABASE/METABOLOMICS_ANNOTATION_PUBMED.txt",sep="/"), delim = "\t",col_names = FALSE )

VOCABOLARY <- vroom(paste(directory,"/MOFA/x_BiomiX_DATABASE/Factor_names_selection/REA_HPO_vocabolary.txt",sep=""), delim = "\t", col_names = FALSE)


TYPE_OF_RESEARCH <- as.character(COMMAND_ADVANCED[1,10])   #Text Word  #"Title/Abstract"
retmaxi <- as.numeric(COMMAND_ADVANCED[2,10])  
n_top_contributors_pubmed = as.numeric(COMMAND_ADVANCED[1,22])
n_keywords_generated = as.numeric(COMMAND_ADVANCED[2,22])


#DEBUGGING SECTION
#########################################################################################
#TO IMPROVE CODE STABILITY IN CASE OF PUBMED LACK OF CONTINUOUS CONNECTIOR OR INTERNET
#NOT STABLE, THIS COSE WILL PROVIDE MORE STABILITY RE-RUNNING THE CODE THAT WAS AFFECTED BY INTERNET LACKING.

fetch_abstract <- function(pubIDs, max_retries = 5) {
        retries <- 0
        
        while (retries < max_retries) {
                retries <- retries + 1
                try_result <- tryCatch({
                        abstract <- entrez_fetch(db = "pubmed", id = pubIDs, rettype = "abstract", retmode = "text")
                        return(abstract)
                }, error = function(err) {
                        message(paste("PubMed ID:", pubIDs, "- Attempt", retries, "failed. Error:", conditionMessage(err)))
                        Sys.sleep(5)  # Add a delay before retrying
                        return(NULL)  # Return NULL if an error occurs
                })
                
                if (!inherits(try_result, "try-error")) {
                        return(try_result)  # Return the result if successful
                }
        }
        
        message(paste("Max retries reached for PubMed ID:", pubmed_id))
        return(NULL)  # Return NULL if max retries reached
}


#########################################################################################

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



get_text_mined_keywords <- function(uu) {
for (pubID in uu$Pubmed_ID) {
        #pubmed_record <- entrez_fetch(db = "pubmed", id = pubID, rettype = "abstract", retmode = "text")
        pubmed_record<- fetch_abstract(pubID)
        print(pubID)
        print(pubmed_record)
        pubmed_record <- gsub("\n", "", pubmed_record)
        #print(pubmed_record)
        
        terms_articles <- extract_terms(
                text = pubmed_record,
                method = "fakerake",
                min_freq=1,
                min_n=2
        )
        
        terms_articles_sum <- append(terms_articles_sum, terms_articles)
}
        
frequency <- table(terms_articles_sum)
frequency <-frequency[order(-frequency)]
#Filtering text word not interesting..
terms_to_keep <- VOCABOLARY$X1
# Remove the terms from the frequency table using grep
updated_frequency <- frequency[frequency > 1]

frequency_split<-unlist(str_split(names(updated_frequency), " "))
frequency_split<-unlist(str_split(frequency_split, "-"))
frequency_split <-toupper(frequency_split)

frequency_split<- frequency_split[frequency_split %in% terms_to_keep]

if(length(frequency_split) != 0){
frequency_split <-tolower(unique(frequency_split))

total<-grepl("VERY_UNLIKELY_TEXT", names(updated_frequency))
for (motif in 1:length(frequency_split)){
        single<- grepl(frequency_split[motif], names(updated_frequency))
        print(frequency_split[motif])
        total <-single | total
}

updated_frequency<-updated_frequency[total]

#print(updated_frequency)
return(updated_frequency)
}
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
        #HERE INSERT A VOCABULARY TO CONVERT THE NAMES THAT ARE COMPLEX
        
        # reduction_words<-str_count(word_list, "-") >= 3
        # reduction_words[is.na(reduction_words)] <- FALSE
        # if (sum(reduction_words) != 0){
        #         #word_list[reduction_words]<-unlist(lapply(lapply(strsplit(as.character(word_list[reduction_words]), split = '-'), function(x) x[c(length(x))]), paste, collapse = '-'))
        #         #word_list[reduction_words]<-unlist(lapply(lapply(strsplit(as.character(word_list[reduction_words]), split = '-'), function(x) x[c(length(x)-1,length(x))]), paste, collapse = '-')) #TWO LAST ONES
        #         word_list<-unlist(lapply(lapply(strsplit(as.character(word_list), split = '-'), function(x) x[c(length(x))]), paste, collapse = '-'))
        #         
        #         }
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
                position_list <- which(COMMAND$LABEL %in% substr(omik, 1, nchar(omik)-1))
                val<-grep(paste(omik,".*_pos_","*", sep=""), files)
                LIST_GENE<- vroom(paste(directory2,"/",files[val], sep=""), delim="\t")
                LIST_GENE<-LIST_GENE %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
                #IF there is only peak the query is removed, peak is removed from the total name
                LIST_GENE$features <-gsub("peak[0-9]*/", "", LIST_GENE$features)
                
                #Removal peak name
                if (COMMAND$DATA_TYPE[position_list] == "Metabolomics"){
                sav<-grep("peak[0-9]*", LIST_GENE$features)
                if(length(sav) == 0){
                        print("selection...")}else{
                                LIST_GENE<- LIST_GENE[-sav,]}}
                
                val<-grep(paste(omik,".*_neg_","*", sep=""), files)
                LIST_GENE2<-vroom(paste(directory2,"/",files[val], sep=""), delim="\t")
                LIST_GENE2<-LIST_GENE2 %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
                LIST_GENE2$features <-gsub("peak[0-9]*/", "", LIST_GENE2$features)
                
                #Removal peak name
                if (COMMAND$DATA_TYPE[position_list] == "Metabolomics"){
                sav<-grep("peak[0-9]*", LIST_GENE2$features)
                if(length(sav) == 0){
                        print("selection...")}else{
                                LIST_GENE2<- LIST_GENE2[-sav,]}}
                
                
                LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
                LIST_GENE3<-LIST_GENE3 %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
         
                #Retrieve gene name from cpg islands
                if (COMMAND$DATA_TYPE[position_list] == "Transcriptomics" | COMMAND$DATA_TYPE[position_list] == "Methylomics" ){
                LIST_GENE3$features <-gsub("cg[0-9]*/", "", LIST_GENE3$features)
                sav<-grep("NA", LIST_GENE3$features)
                if(length(sav) != 0){
                LIST_GENE3<- LIST_GENE3[-sav,]}
                
                #TEST TO REMOVE COIL GENE (ONLY TO TEST)
                sav<-grep("COIL", LIST_GENE3$features)
                if(length(sav) != 0){
                        LIST_GENE3<- LIST_GENE3[-sav,]}}
                
                #REMOVING NA FEATURES (ALL OMICS)
                sav<- LIST_GENE3$features == "NA"
                LIST_GENE3<- LIST_GENE3[!sav,]
                
                #Renaming metabolites to improve research
                if (COMMAND$DATA_TYPE[position_list] == "Metabolomics"){
                        LIST_GENE3$features<-gsub("\\(.*\\)","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("\\[.*\\]","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("\\{.*\\}","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("^-[0-9]-","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("\\'","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("\\,","-",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub(";.*","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("--","-",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("^-","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("^ ","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub(" $","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("-$","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("/.*","",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub(" -","-",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("- ","-",LIST_GENE3$features)
                        LIST_GENE3$features<-gsub("- ","-",LIST_GENE3$features)
                        for (value in 1:length(LIST_GENE3$features)){
                                out<-ANNOTATION_LIPIDS$X1 %in% LIST_GENE3$features[value]
                                if (sum(out) == 1){
                                        LIST_GENE3$features[value] <- ANNOTATION_LIPIDS$X2[grep(TRUE, out)] 
                                }
                        }
                        LIST_GENE3_out<- unlist(lapply(lapply(strsplit(as.character(LIST_GENE3$features), split = '-'), function(x) x[c(length(x))]), paste, collapse = '-'))
                        
                }else{
                        LIST_GENE3_out <-  LIST_GENE3$features     
                        }
                
                
                if(length(LIST_GENE3$features) < n_top_contributors_pubmed){
                        x<-paste(LIST_GENE3$features[1:length(LIST_GENE3$features)-1], "OR ", sep= paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""))
                        xx<-c(x, paste(LIST_GENE3$features[length(LIST_GENE3$features)],paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""),sep=""))
                        if (COMMAND$DATA_TYPE[position_list] == "Metabolomics"){
                                columns<-LIST_GENE3_out[1:length(LIST_GENE3$features)]
                        }else{columns<-LIST_GENE3$features[1:length(LIST_GENE3$features)]}
                }else{
                x<-paste(LIST_GENE3$features[1:(n_top_contributors_pubmed -1)], "OR ", sep= paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""))
                xx<-c(x, paste(LIST_GENE3$features[n_top_contributors_pubmed],paste("[", TYPE_OF_RESEARCH, "]"," ", sep=""),sep=""))
                if (COMMAND$DATA_TYPE[position_list] == "Metabolomics"){
                        columns<-LIST_GENE3_out[1:n_top_contributors_pubmed]
                }else{columns<-LIST_GENE3$features[1:n_top_contributors_pubmed]}
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
        
        #print(yyy)
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
                #nodes
                
                articles <- xml_find_all(parsed_xml, "//PubmedArticle")
                
                if (retmaxi > length(articles)) {
                        articles<-articles
                }else{
                        articles<-articles[1:retmaxi]}
                
                
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
                        print(doi)
                        
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
                        top_3_rows <- order(-word_match$tot_word, -word_match$tot_match)[1:15]  
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
        
        

files_tot <- grep(paste("MOFA", "_", args[1] ,"_vs_", args[2],"*", sep=""),list.files(paste(directory,"/MOFA/OUTPUT/",sep="")),value=TRUE)
directory <- args[[3]]

######
# fil <- files_tot[1]
# retmaxi <- 50
#####

for (fil in files_tot){
print(paste("folder analyzed:", fil))

MART <- vroom(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")

dir.create(path = paste(directory,"/MOFA/OUTPUT/",fil, "/Factors_Associated_Articles" ,sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory,"/MOFA/OUTPUT/",fil,sep="") 
directory3 <- paste(directory,"/MOFA/OUTPUT/",fil, "/Factors_Associated_Articles" ,sep="")
setwd(directory2)


#detection total significant factors
files <- grep("*factor_",list.files(directory2),value=TRUE)
files<- files[grep("\\Metabolomics|\\RNAseq|\\Methylomics", files)]
factors<-strsplit(files, "_")
factors<-unlist(factors)
factors<-str_remove(factors,".tsv")
factors<-unique(as.numeric(factors[!is.na(as.numeric(factors))]))

####
# numb <- factors[1]
####

if (length(factors) != 0){

for (numb in factors){
print(paste("factor analyzed:", numb, "*********************************"))        
#####BLOCK RESEARCH QUERY ALL OMICS TOGETHER ####

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Metabolomics|\\RNAseq|\\Methylomics", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])

iter= 0

if(length(nam) != 0){
uu<-count_match_in_top_article(nam)

if (length(uu) > 1){
uu$Keywords <- gsub("\n", "", uu$Keywords)
terms_articles_sum =NULL
it<-as.data.frame(get_text_mined_keywords(uu))
if (nrow(it) > n_keywords_generated){
        it<-it[1:n_keywords_generated,]   
}

mined<-paste(it$terms_articles_sum,it$Freq, collapse = "/")
uu$Keywords_text_mined <- mined
}


line="#ALL OMICS (TRASCRIPTOMICS-METABOLOMICS-METHYLOMICS)#"
write(line,file=paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""))
if (length(uu) == 1){
write("no_results",file=paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)} else{
        write.table(uu,file=paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
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

if (length(uu) > 1){
uu$Keywords <- gsub("\n", "", uu$Keywords)
terms_articles_sum =NULL
it<-as.data.frame(get_text_mined_keywords(uu))
if (nrow(it) > n_keywords_generated){
        it<-it[1:n_keywords_generated,]   
}

mined<-paste(it$terms_articles_sum,it$Freq, collapse = "/")
uu$Keywords_text_mined <- mined
}

line="#PAIRS OMICS (TRASCRIPTOMICS-METABOLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)
if (length(uu) == 1){
        write("no_results",paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
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

if (length(uu) > 1){
uu$Keywords <- gsub("\n", "", uu$Keywords)
terms_articles_sum =NULL
it<-as.data.frame(get_text_mined_keywords(uu))
if (nrow(it) > n_keywords_generated){
        it<-it[1:n_keywords_generated,]   
}

mined<-paste(it$terms_articles_sum,it$Freq, collapse = "/")
uu$Keywords_text_mined <- mined
}

line="#PAIRS OMICS (TRASCRIPTOMICS-METHYLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)
if (length(uu) == 1){
        write("no_results",file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
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

if (length(uu) > 1){
uu$Keywords <- gsub("\n", "", uu$Keywords)
terms_articles_sum =NULL
it<-as.data.frame(get_text_mined_keywords(uu))
if (nrow(it) > n_keywords_generated){
        it<-it[1:n_keywords_generated,]   
}

mined<-paste(it$terms_articles_sum,it$Freq, collapse = "/")
uu$Keywords_text_mined <- mined
}

line="#PAIRS OMICS (METABOLOMICS-METHYLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)
if (length(uu) == 1){
        write("no_results",file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
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

if (length(uu) > 1){
uu$Keywords <- gsub("\n", "", uu$Keywords)
terms_articles_sum = NULL
it<-as.data.frame(get_text_mined_keywords(uu))
if (nrow(it) > n_keywords_generated){
        it<-it[1:n_keywords_generated,]   
}

mined<-paste(it$terms_articles_sum,it$Freq, collapse = "/")
uu$Keywords_text_mined <- mined
}

line="#SINGLE OMICS (TRASCRIPTOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)
if (length(uu) == 0){
        write("no_results",file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")
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

if (length(uu) > 1){
uu$Keywords <- gsub("\n", "", uu$Keywords)
terms_articles_sum =NULL
it<-as.data.frame(get_text_mined_keywords(uu))
if (nrow(it) > n_keywords_generated){
        it<-it[1:n_keywords_generated,]   
}

mined<-paste(it$terms_articles_sum,it$Freq, collapse = "/")
uu$Keywords_text_mined <- mined
}

line="#SINGLE OMICS (METABOLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)
if (length(uu) == 0){
        write("no_results",file=paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)} else{
                write.table(uu,file=paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }
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

if (length(uu) > 1){
uu$Keywords <- gsub("\n", "", uu$Keywords)
terms_articles_sum =NULL
it<-as.data.frame(get_text_mined_keywords(uu))
if (nrow(it) > n_keywords_generated){
        it<-it[1:n_keywords_generated,]   
}

mined<-paste(it$terms_articles_sum,it$Freq, collapse = "/")
uu$Keywords_text_mined <- mined
}

line="#SINGLE OMICS (METHYLOMICS)#"
write(line,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)
if (length(uu) == 0){
        write("no_results",file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE)} else{
                write.table(uu,file= paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        }

}



terms_articles_tot = NULL
pubmed_codes =NULL
terms_articles_single=NULL

Total_articles <- vroom(paste(directory3, "/", "Factor_", numb, "_articles.tsv", sep=""), delim = "/t", col_names=FALSE)
Total_articles<- Total_articles[!grepl("#.*", Total_articles$X1),]
Total_articles<- Total_articles[!(Total_articles$X1 == "x"),]
Total_articles<- Total_articles[!grepl("no_results", Total_articles$X1),]
frame_identification<- which(grepl("Pubmed", Total_articles$X1))
frame_identification<- c(frame_identification,(length(Total_articles$X1) + 1))

for (number in 1:(length(frame_identification)-1) ){
        print(paste(number,"*************************"))
frame_identifications <- frame_identification[number + 1] - frame_identification[number]
next_spot <- frame_identification[number] + frame_identifications - 1
df_split <- data.frame(do.call(rbind, strsplit(Total_articles$X1[frame_identification[number]:next_spot], "\t")))
colnames(df_split) <- df_split[1,]
df_split<- df_split[-1,]
pubmed_codes <- append(pubmed_codes, df_split$Pubmed_ID)
print(pubmed_codes)

for (pubID in pubmed_codes) {
#pubmed_record <- entrez_fetch(db = "pubmed", id = pubID, rettype = "abstract", retmode = "text")
        pubmed_record<- fetch_abstract(pubID)
        pubmed_record <- gsub("\n", "", pubmed_record)
#print(pubmed_record)

terms_articles <- extract_terms(
        text = pubmed_record,
        method = "fakerake",
        min_freq=1,
        min_n=2
)

terms_articles_tot <- append(terms_articles_tot, terms_articles)
}

pubmed_codes =NULL
terms_articles_single=NULL

}



}

frequency <- table(terms_articles_tot)
frequency <-frequency[order(-frequency)]
#Filtering text word not interesting..
terms_to_keep <- VOCABOLARY$X1
# Remove the terms from the frequency table using grep
updated_frequency <- frequency[frequency > 1]

frequency_split<-unlist(str_split(names(updated_frequency), " "))
frequency_split<-unlist(str_split(frequency_split, "-"))
frequency_split <-toupper(frequency_split)

frequency_split<- frequency_split[frequency_split %in% terms_to_keep]

if(length(frequency_split) != 0){
        frequency_split <-tolower(unique(frequency_split))
        
        total<-grepl("VERY_UNLIKELY_TEXT", names(updated_frequency))
        for (motif in 1:length(frequency_split)){
                single<- grepl(frequency_split[motif], names(updated_frequency))
                print(frequency_split[motif])
                total <-single | total
        }
        
        updated_frequency<-updated_frequency[total]
        
        write.table(updated_frequency,file= paste(directory3, "/", "Keywords_mined_Factor_", numb,"_all_omics.tsv", sep=""),append=TRUE, quote = FALSE, row.names = FALSE, sep="\t")        
        }


}

}



        

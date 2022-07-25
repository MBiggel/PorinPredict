#!/usr/bin/env Rscript --vanilla

suppressMessages(library(dplyr))
suppressMessages(library(Biostrings))

# Retrieve arguments from Python
args <- commandArgs(trailingOnly = TRUE)

# Import Diamond and Blastn results
df = read.table(args[2],sep="\t",header=FALSE)
dfp = read.table(args[1],sep="\t",header=FALSE)

# Define output
outfile = paste(args[4],args[3],"_PorinPredict.tsv",sep="")

### Classification OprD CDS
# Diamond: no hit
if(df[1,1] == "missing") {
  names(df) = c("result_CDS","Assembly")
  df$"Subject ID" = NA
  df$Mutation = NA
  df$Allele = NA

# Diamond: all others   
} else {
names(df) = c("Query ID","Subject ID","Subject length","Percentage of identical matches",
              "Subject coverage per HSP","Alignment length","Number of mismatches",
              "Number of gap openings","Start of alignment in query",
              "End of alignment in query","Start of alignment in subject",
              "End of alignment in subject","Expected value","Bit score",
              "Query strand","Aligned part of query sequence (translated)",
              "Subject sequence","Assembly")


# Function: Identify amino acid substitutions
mismatches <- function(query, ref) {
  pairwiseAlignment(ref, query, substitutionMatrix = "BLOSUM50",
                    gapOpening = 3, gapExtension = 1) %>%
    mismatchTable() %>%
    mutate(ID=names(query), 
           Pos=PatternStart, 
           Reference_AA=as.character(PatternSubstring),
           Sample_AA=as.character(SubjectSubstring)) %>% 
    select(ID, Reference_AA, Sample_AA, Pos)
}  

# Check if terminal stop present
df$Terminal_stop = substr(df$"Aligned part of query sequence (translated)",
                           nchar(df$"Aligned part of query sequence (translated)"),
                           nchar(df$"Aligned part of query sequence (translated)")) == "*"

# Check for premature stop; loc_stop gives -1 if no stop detected
loc_stop = unlist(gregexpr("\\*", df$"Aligned part of query sequence (translated)"))[1]
df$Premature_stop = (loc_stop != nchar(df$"Aligned part of query sequence (translated)") && loc_stop != -1)

# Check if oprD complete
df$"Aligned query length/Subject length*100" = 
  round(nchar(df$"Aligned part of query sequence (translated)") / df$"Subject length" *100,1)
df$OprD_complete = (df$"Aligned query length/Subject length*100" == 100.0)

# Check if OprD complete excluding subject terminal stop
df$OprD_complete_wo_terminal_stop = (df$Terminal_stop == FALSE) &&
  nchar(df$"Aligned part of query sequence (translated)") == df$"Subject length"-1

# Check if start (Met) present
df$Start_Met = substr(df$"Aligned part of query sequence (translated)",1,1) == "M"

# Check if exact match
df$Exact_match = (df$"Aligned query length/Subject length*100" == 100 && df$"Percentage of identical matches" == 100)

# PorinPredict classification
if(df$OprD_complete == TRUE && df$Exact_match == TRUE){
  df$result_CDS = "intact - exact match"
  df$Allele = paste(df$"Subject ID","--",df$"Aligned query length/Subject length*100","--",df$"Percentage of identical matches","%",sep="")
  df$Mutation = NA
  
  } else if(df$Terminal_stop == TRUE && df$OprD_complete == TRUE && df$Premature_stop == FALSE && df$Start_Met == TRUE && df$Exact_match == FALSE) {
  df$result_CDS = "intact - missense mutation"
  df$Allele = paste(df$"Subject ID","--",df$"Aligned query length/Subject length*100","--",df$"Percentage of identical matches","%",sep="")
  df$Mutation = NA
  
  # Identify substitution
  seqs = c(AAStringSet(df$"Subject sequence"),AAStringSet(df$"Aligned part of query sequence (translated)"))
  names(seqs) = c("ref","query")
  table = bind_rows(mismatches(seqs[2], seqs[1]))

  table$Mutation = paste0(table$Reference_AA,table$Pos,table$Sample_AA,sep="")
  df$Mutation = paste(table$Mutation, collapse =", ")
  
} else if(df$OprD_complete == FALSE && df$OprD_complete_wo_terminal_stop == FALSE) {
  df$result_CDS = "incomplete"
  df$Allele = paste(df$"Subject ID","--",df$"Aligned query length/Subject length*100","--",df$"Percentage of identical matches","%",sep="")
  df$Mutation = NA
  
} else if(df$Premature_stop == TRUE) {
  df$result_CDS = "premature stop"
  df$Allele = paste(df$"Subject ID","--",df$"Aligned query length/Subject length*100","--",df$"Percentage of identical matches","%",sep="")
  
  df$Mutation = paste0(substr(df$"Subject sequence",loc_stop,loc_stop),loc_stop,"X",sep="")
  
} else if(df$OprD_complete == TRUE && df$Start_Met == FALSE) {
  df$result_CDS = "start loss"
  df$Allele = paste(df$"Subject ID","--",df$"Aligned query length/Subject length*100","--",df$"Percentage of identical matches","%",sep="")

  df$Mutation = paste0("M1",substr(df$"Aligned part of query sequence (translated)",1,1),sep="")
  
} else if (df$Terminal_stop == FALSE && df$OprD_complete_wo_terminal_stop == TRUE) {
  df$result_CDS = "stop loss"
  df$Allele = paste(df$"Subject ID","--",df$"Aligned query length/Subject length*100","--",df$"Percentage of identical matches","%",sep="")
  df$Mutation = NA
  
} else {
  df$result_CDS = "Unrecognized mutation"
  df$Mutation = NA
}}

### Classification OprD promoter
if(dfp[1,1] == "missing") {
  names(dfp) = c("result_promoter","Assembly")
  df$promoter = "missing"

  } else {
  names(dfp) = c('qseqid','sseqid','slen','pident','qcovs','length',
                 'mismatch','gaps','gapopen','qstart','qend','sstart',
                 'send','evalue','bitscore','Assembly')
  
  coverage = dfp$length / dfp$slen * 100
  
  if(coverage >=98) {
    df$promoter = paste("complete_",round(coverage,2),"%",sep="")
    df$promoter2 = "complete"

  } else {
    df$promoter = paste("incomplete_",round(coverage,2),"%",sep="")
    df$promoter2 = "incomplete"
}}

### Summarize
if(df$result_CDS == "intact - exact match" && df$promoter2 == "complete") {
  df$"PorinPredict" = "OprD intact"
} else if(df$result_CDS == "intact - missense mutation" && df$promoter2 == "complete") {
  df$"PorinPredict" = "OprD intact but rare AA substitution"
} else {
  df$"PorinPredict" = "OprD inactivated"
}

out = df[,c("Assembly","promoter","Subject ID","result_CDS","Allele","Mutation","PorinPredict")]
names(out) = c("Assembly ID","oprD promoter completeness","OprD variant hit","OprD coding region","variant -- aligned query length/subject length *100 -- subject identity","Missense/nonsense mutation","PorinPredict classification")

write.table(out,outfile,sep="\t",row.names = FALSE)

rm(list=ls())
#!/usr/bin/env Rscript --vanilla

suppressMessages(library(dplyr))
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)

# Import Diamond and Blastn results
df = read.table(args[3],sep="\t",header=FALSE)
dfn = read.table(args[2],sep="\t",header=FALSE)
dfp = read.table(args[1],sep="\t",header=FALSE)

# Define output
outfile = paste(args[5],args[4],"_PorinPredict.tsv",sep="")

# Rename headers
if(dfn[1,1] == "no hit") {
  names(dfn) = c("result_oprD_nuc","Assembly")
} else {
  names(dfn) = c('qseqid','sseqid','slen','pident','qcovs','length',
                 'mismatch','gaps','gapopen','qstart','qend','sstart',
                 'send','evalue','bitscore','Assembly')
}

if(df[1,1] == "no hit") {
  names(df) = c("Result","Assembly")
} else {
  names(df) = c("Query ID","Subject ID","Subject length","Percentage of identical matches",
                "Subject coverage per HSP","Alignment length","Number of mismatches",
                "Number of gap openings","Start of alignment in query",
                "End of alignment in query","Start of alignment in subject",
                "End of alignment in subject","Expected value","Bit score",
                "Query strand","Aligned part of query sequence (translated)",
                "Subject sequence","Assembly")
}

# Function: amino acid substitutions
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

### Presence oprD nucleotide sequence
if(dfn[1,1] == "no hit") {
  df$oprD_nuc = "no hit"
  df$result_CDS = "no hit"
  } else {
  if(dfn$length > dfn$slen) {
    coverage = dfn$length / dfn$slen * 100
    df$oprD_nuc = paste(dfn$length,"/",dfn$slen," (",round(coverage,2),"%)",sep="")       
  } else {
    coverage = (dfn$length - dfn$gaps) / dfn$slen * 100
    df$oprD_nuc = paste((dfn$length - dfn$gaps),"/",dfn$slen," (",round(coverage,2),"%)",sep="")       
  }
  
  if(dfn$gaps %in% c(1,2,4,5,7,8,10,11) && (dfn$length <= (dfn$slen + 12)) && ((dfn$length - dfn$gaps) >= (dfn$slen - 12))) {
    df$gap = "frameshift indel"
  } else if(dfn$gaps %in% c(3,6,9,12) && dfn$gapopen == 1 && dfn$length <= (dfn$slen + 12) & (dfn$length - dfn$gaps) >= (dfn$slen - 12)) {
    df$gap = "in-frame indel"	
  } else if(dfn$length > (dfn$slen + 12)) {
    df$gap = "large insertion"
  } else if((dfn$length - dfn$gaps) < (dfn$slen - 12)) {
    df$gap = "truncated"
  } else {
    df$gap = "multiple indels/other rearrangement"
  }}

### Classification OprD CDS
# Diamond: no hit
if(df[1,1] == "no hit") {
  df$"Subject ID" = "no hit"
  df$Mutation = NA
  df$CDS_coverage = NA
  df$"Percentage of identical matches" = NA
  
  if(dfn[1,1] != "no hit"){
      df$result_CDS = df$gap
  }
# Diamond: hit
} else {
  # alignment length
  df$aln_length = nchar(df$"Aligned part of query sequence (translated)")
  # Presence terminal stop
  df$Terminal_stop = substr(df$"Aligned part of query sequence (translated)",
                            df$aln_length,df$aln_length) == "*"
    
  # Presence premature stop (no stop -> loc_stop = -1)
  loc_stop = unlist(gregexpr("\\*", df$"Aligned part of query sequence (translated)"))[1]
  df$Premature_stop = (loc_stop != df$aln_length) && (loc_stop != -1)
    
  # OprD completeness
  df$"Aligned query length/Subject length*100" = 
    round(df$aln_length / df$"Subject length" *100,1)
  df$OprD_complete = (df$"Aligned query length/Subject length*100" == 100.0)
    
  # OprD completeness, terminal stop (subject) excluded
  df$OprD_complete_wo_terminal_stop = (df$Terminal_stop == FALSE) &&
  nchar(df$"Aligned part of query sequence (translated)") == df$"Subject length"-1
    
  # Start (Met) present
  df$Start_Met = substr(df$"Aligned part of query sequence (translated)",1,1) == "M"
    
  # Exact match
  df$Exact_match = (df$"Aligned query length/Subject length*100" == 100 && df$"Percentage of identical matches" == 100)
    
  # PorinPredict classification
  if(df$OprD_complete == TRUE && df$Exact_match == TRUE){
    df$result_CDS = "intact - exact match"   
    df$Mutation = NA
      
  } else if(df$Terminal_stop == TRUE && df$OprD_complete == TRUE && df$Premature_stop == FALSE && df$Start_Met == TRUE && df$Exact_match == FALSE) {
    df$result_CDS = "intact - missense mutation"
    df$Mutation = NA
      
    # Identify substitution
    seqs = c(AAStringSet(df$"Subject sequence"),AAStringSet(df$"Aligned part of query sequence (translated)"))
    names(seqs) = c("ref","query")
    table = bind_rows(mismatches(seqs[2], seqs[1]))
    table$Mutation = paste0(table$Reference_AA,table$Pos,table$Sample_AA,sep="")
    df$Mutation = paste(table$Mutation, collapse =", ")
      
  } else if(df$OprD_complete == FALSE && df$OprD_complete_wo_terminal_stop == FALSE) {
    df$result_CDS = df$gap
    df$Mutation = NA
      
  } else if(df$Premature_stop == TRUE) {
    df$result_CDS = "premature stop" 
    df$Mutation = paste0(substr(df$"Subject sequence",loc_stop,loc_stop),loc_stop,"*",sep="")
      
  } else if(df$OprD_complete == TRUE && df$Start_Met == FALSE) {
    df$result_CDS = "start loss"
    df$Mutation = paste0("M1",substr(df$"Aligned part of query sequence (translated)",1,1),sep="")
      
  } else if (df$Terminal_stop == FALSE && df$OprD_complete_wo_terminal_stop == TRUE) {
    df$result_CDS = "stop loss"
    df$Mutation = NA
      
  } else {
    df$result_CDS = "Unrecognized mutation"
    df$Mutation = NA
  }
    
  if(df$Terminal_stop == TRUE){
    df$CDS_coverage = paste(nchar(df$"Aligned part of query sequence (translated)")-1,"/",df$"Subject length"-1," (",round((nchar(df$"Aligned part of query sequence (translated)")-1)/(df$"Subject length"-1) *100,1),"%)",sep="")
  } else {
    df$CDS_coverage = paste(nchar(df$"Aligned part of query sequence (translated)"),"/",df$"Subject length"-1," (",round(nchar(df$"Aligned part of query sequence (translated)")/(df$"Subject length"-1) *100,1),"%)",sep="")
  }
}
  
### Classification OprD promoter
if(dfp[1,1] == "no hit") {
  names(dfp) = c("result_promoter","Assembly")
  df$promoter = "no hit"
  df$promoter_ev = "no hit"
    
} else {
  names(dfp) = c('qseqid','sseqid','slen','pident','qcovs','length',
                 'mismatch','gaps','gapopen','qstart','qend','sstart',
                 'send','evalue','bitscore','Assembly')
    
  coverage = dfp$length / dfp$slen * 100
    
  if(coverage >=98) {
    df$promoter = paste("complete (",round(coverage,2),"%)",sep="")
    df$promoter_ev = "complete"
      
  } else {
    df$promoter = paste("incomplete (",round(coverage,2),"%)",sep="")
    df$promoter_ev = "incomplete"
  }
}
  
### Summary
if(df$result_CDS == "intact - exact match" && df$promoter_ev == "complete") {
  df$"PorinPredict" = "OprD intact"
} else if(df$result_CDS == "intact - missense mutation" && df$promoter_ev == "complete") {
  df$"PorinPredict" = "OprD intact but rare AA substitution"
} else if(df$oprD_nuc == "no hit" && df$promoter == "no hit") {
  df$"PorinPredict" =  "WARNING: OprD nucleotide coding and promoter sequences not detected. Verify species affiliation and assembly quality."
} else {
  df$"PorinPredict" = "OprD inactivated"
}

out = df[,c("Assembly","promoter","oprD_nuc","Subject ID","CDS_coverage","Percentage of identical matches","result_CDS","Mutation","PorinPredict")]
names(out) = c("Assembly ID","oprD promoter integrity (alignment coverage)","oprD nucleotide sequence alignment coverage ","OprD variant hit (AA level)","OprD alignment coverage (AA level)","OprD sequence identity (%; AA level)","OprD integrity (AA level)","Missense/nonsense mutations","PorinPredict final classification")

write.table(out,outfile,sep="\t",row.names = FALSE)

rm(list=ls())

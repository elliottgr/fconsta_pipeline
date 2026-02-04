args = commandArgs(trailingOnly=TRUE)
## handling missing filename
if (length(args)==0){
    filename<-paste("outputs/pa_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
} else if (length(args)>=1) {
    filename<-paste("outputs/", args[1], "_pa_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
}


source("scripts/04_filter_blast_results.R")
results<-process_blast_hits(input="outputs/blast_out.tsv", min_pident = 70, min_qcovs = 50, annotation_file = NULL)
write.csv(results$df_pa, filename)
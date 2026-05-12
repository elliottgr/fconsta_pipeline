## Simple script. Takes the BLAST outputs from script 3, invokes 
## script 4 to make a PA table, then saves everything with a nice filename 

args = commandArgs(trailingOnly=TRUE)
## handling missing filename
if (length(args)==0){
    filename_pa<-paste("outputs/pa_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
    filename_filtered<-paste("outputs/filtered_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
    filename_pa_raw<-paste("outputs/pa_raw_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
    filename_ranked<-paste("outputs/ranked_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
} else if (length(args)>=1) {
    filename_pa<-paste("outputs/", args[1], "_pa_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
    filename_filtered<-paste("outputs/", args[1], "_filtered_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
    filename_pa_raw<-paste("outputs/", args[1], "_pa_raw_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
    filename_ranked<-paste("outputs/", args[1], "_ranked_table", format(Sys.Date(),"_%Y-%b-%d"), ".csv",sep="")
}


source("scripts/04_filter_blast_results.R")
results<-process_blast_hits(input="outputs/blast_outputs.tsv", min_pident = 85, min_scov = 75, annotation_file = NULL)
print("Saving tables!")
write.csv(results$df_filtered, filename_filtered)
write.csv(results$df_ranked, filename_ranked)
print("Working :)")
write.csv(results$df_pa, filename_pa)
write.csv(results$df_pa_raw, filename_pa_raw)
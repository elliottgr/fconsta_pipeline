## Rscript to calculate phylogenetic D (Fritz and Orme 2010) on the presence absence data
## In brief, this determines how much of the binary character data is distributed by phylogenetic
## signal. Significant phylo signal would suggest gene variants are being distributed vertically


library("caper")
library("phytools") ## the goat, need it to load phylo tree
## need to load filtered presence absence data (output of scripts/07_format_PA_table.py)
df = read.csv("outputs/labeled_pa_table.csv")# %>% gsub("-", "_", Strain)#%>% rename(names = Strain)
df_filtered <- Filter(function(x)(length(unique(x))>1), df)
tree = read.newick(file = "extra/supplemental_data/Croucher.nwk") ## if you're running this yourself, you'll need to supply your own newick format tree.

## insane behavior please delete this language remove R 
good_columns = colnames(df_filtered)
## removing first few dummy columns
good_columns = good_columns[-1]
good_columns = good_columns[-1]
good_columns = good_columns[-1]
good_columns = good_columns[-1]
l<-c()
i<-1
while (i < length(good_columns)){
for (col in good_columns) {
  l[[i]] <- paste("phylo.d(df_filtered, tree, Strain, binvar=", col , ")", sep = "")
  i<- i+1
}
}
str_eval=function(x){return(eval(parse(text=x)))}

gene_names<-c()
p_no_structure<-c()
p_brownian<-c()
est_d<-c()
num_genes<-1
while (num_genes<length(l)){
for (element in l){
  tryCatch({output = str_eval(element)}, error = function(msg){print("error with gene!")})
  print(output)
  gene_names[num_genes]<-output$binvar
  est_d[num_genes]<-output$DEstimate
  p_no_structure[num_genes]<-output$Pval1
  p_brownian[num_genes]<-output$Pval0
  num_genes<-num_genes+1
}
}

output_list<- list(gene_names, p_no_structure, p_brownian, est_d)
output<- as.data.frame(do.call(cbind, output_list)) %>% rename(gene_variant = V1, p_no_structure = V2, p_brownian_structure = V3, estimated_d = V4)
write.csv(output, file="d_value_estimates.csv")

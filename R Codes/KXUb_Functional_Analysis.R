library(ggplot2)

#Hep2_Jurkat
metascape_enrichment <- read.csv('Plot/functional_analysis/Hep2_Jurkat/HeatmapSelectedGO.csv')
metascape_enrichment <- metascape_enrichment[,2:4]
metascape_enrichment <- reshape2::melt(metascape_enrichment)

ggplot(metascape_enrichment, aes(x= -value, y=reorder(Description,-value), fill=variable)) + geom_bar(stat='identity', position='dodge') + 
  theme_bw() + xlab('-Log10(P-value)') + ylab('')

metascape_cell_type <- read.csv('Plot/functional_analysis/Hep2_Jurkat/HeatmapSelectedGO_PaGenBase.csv')
metascape_cell_type <- metascape_cell_type[,2:4]
metascape_cell_type$X_LogP_Hep2 <- - metascape_cell_type$X_LogP_Hep2
metascape_cell_type <- reshape2::melt(metascape_cell_type)

ggplot(metascape_cell_type, aes(x= -value, y=reorder(Description,-value), fill=variable)) + geom_bar(stat='identity', position='stack') + 
  theme_bw() + xlab('-Log10(P-Value)') + ylab('')

metascape_enrichment <- read.csv('Plot/functional_analysis/Hep2_Jurkat/HeatmapSelectedGO.csv')
temp_a <- metascape_enrichment[order(metascape_enrichment$X_LogP_Hep2)[1:15],]
temp_b <- metascape_enrichment[order(metascape_enrichment$X_LogP_Jurkat)[1:15],]
temp <- rbind(temp_a, temp_b)
temp <- temp[!duplicated(temp),]

temp$X_LogP_Jurkat <- - temp$X_LogP_Jurkat
temp <- temp[,2:4]
temp <- reshape2::melt(temp)

ggplot(temp, aes(x= value, y=reorder(Description,abs(value)), fill=variable)) + geom_bar(stat='identity', position='stack') + 
  scale_fill_manual(values = c('#E7211A','#187C3B')) +
  theme_bw() + xlab('-Log10(P-value)') + ylab('')

#Hep2 inhibitors
metascape_enrichment <- read.csv('Dataset/metascape/Inhibitors/Hep2/Enrichment_heatmap/HeatmapSelectedGOTop100.csv')
temp_a <- metascape_enrichment[order(metascape_enrichment$X_LogP_CTR)[1:5],]
temp_b <- metascape_enrichment[order(metascape_enrichment$X_LogP_Inhibitors)[1:5],]
temp_c <- metascape_enrichment[metascape_enrichment$X_LogP_CTR == 0,]
temp_c <- temp_c[order(temp_c$X_LogP_Inhibitors)[1:5],]
temp_d <- metascape_enrichment[metascape_enrichment$X_LogP_Inhibitors == 0,]
temp_d <- temp_d[order(temp_d$X_LogP_CTR)[1:2],]
temp <- rbind(temp_a, temp_b, temp_c, temp_d)
temp <- temp[!duplicated(temp),]

temp$X_LogP_Inhibitors <- - temp$X_LogP_Inhibitors
temp <- temp[,2:4]
temp <- reshape2::melt(temp)

ggplot(temp, aes(x= value, y=reorder(Description,abs(value)), fill=variable)) + geom_bar(stat='identity', position='stack') + 
  scale_fill_manual(values = c('#187C3B','#E7211A')) +
  theme_bw() + xlab('-Log10(P-value)') + ylab('')

#Jurkat inhibitors
metascape_enrichment <- read.csv('Dataset/metascape/Inhibitors/Jurkat/Enrichment_heatmap/HeatmapSelectedGOTop100.csv')
temp_a <- metascape_enrichment[order(metascape_enrichment$X_LogP_CTR)[1:6],]
temp_b <- metascape_enrichment[order(metascape_enrichment$X_LogP_Inhibitors)[1:5],]
temp_c <- metascape_enrichment[metascape_enrichment$X_LogP_CTR == 0,]
temp_c <- temp_c[order(temp_c$X_LogP_Inhibitors)[1:6],]
temp_d <- metascape_enrichment[metascape_enrichment$X_LogP_Inhibitors == 0,]
temp_d <- temp_d[order(temp_d$X_LogP_CTR)[1:3],]
temp <- rbind(temp_a, temp_b, temp_c, temp_d)
temp <- temp[!duplicated(temp),]

temp$X_LogP_Inhibitors <- - temp$X_LogP_Inhibitors
temp <- temp[,2:4]
temp <- reshape2::melt(temp)

ggplot(temp, aes(x= value, y=reorder(Description,abs(value)), fill=variable)) + geom_bar(stat='identity', position='stack') + 
  scale_fill_manual(values = c('#187C3B','#E7211A')) +
  theme_bw() + xlab('-Log10(P-value)') + ylab('')

#kxub and kub compare
enri_inter <- read.csv('Dataset/metascape/kub_kxub/kxub_kub_intersect.csv')
enri_only <- read.csv('Dataset/metascape/kub_kxub/kxub_only.csv')
enri_inter <- enri_inter[1:10,]
enri_only <- enri_only[1:10,]
enri_only$Log10.P. <- - enri_only$Log10.P.
enri_temp <- rbind(enri_only, enri_inter)
enri_temp <- enri_temp[,c(3,6)]

ggplot(enri_temp, aes(x= Log10.P. , y=reorder(Description,(Log10.P.)))) + geom_bar(stat='identity') + 
  theme_bw() + xlab('-Log10(P-value)') + ylab('')

#Hep2_Jurkat inhibitors enriched
metascape_enrichment <- read.csv('Plot/functional_analysis/Hep2_Jurkat_inhibitors/HeatmapSelectedGO.csv')
temp_a <- metascape_enrichment[order(metascape_enrichment$X_LogP_Hep2)[1:15],]
temp_b <- metascape_enrichment[order(metascape_enrichment$X_LogP_Jurkat)[1:15],]
temp <- rbind(temp_a, temp_b)
temp <- temp[!duplicated(temp),]

temp$X_LogP_Jurkat <- - temp$X_LogP_Jurkat
temp <- temp[,2:4]
temp <- reshape2::melt(temp)

ggplot(temp, aes(x= value, y=reorder(Description,abs(value)), fill=variable)) + geom_bar(stat='identity', position='stack') + 
  scale_fill_manual(values = c('#E7211A','#187C3B')) +
  theme_bw() + xlab('-Log10(P-value)') + ylab('')

#STRING results of kxub and kub analysis
string_result <- read.csv('Dataset/metascape/kub_kxub/string/log10_protein_abundance_data.txt', sep='\t')
ggplot(string_result, aes(x=as.character(inputValue), y=propertyValue, fill=as.character(inputValue))) + 
  geom_boxplot() + theme_bw() +
  xlab('') + ylab('Log10(protein abundance)')
t.test(string_result$propertyValue[string_result$inputValue == 1],string_result$propertyValue[string_result$inputValue == -1])

string_result <- read.csv('Dataset/metascape/kub_kxub/string/log10_nr_of_publications_data.txt', sep='\t')
ggplot(string_result, aes(x=as.character(inputValue), y=propertyValue, fill=as.character(inputValue))) + 
  geom_boxplot() + theme_bw() +
  xlab('') + ylab('Log10(num of publications)')
t.test(string_result$propertyValue[string_result$inputValue == 1],string_result$propertyValue[string_result$inputValue == -1])

string_list <- read.csv('Dataset/metascape/kub_kxub/string/enrichment.all.tsv', sep='\t')
string_list <- string_list[string_list$term.description == 'Disulfide bond',10]
string_list <- strsplit(string_list,',')[[1]]
string_result$propertyValue <- apply(string_result,1,function(x){ifelse(x[3] %in% string_list,1,0)})
ggplot(string_result, aes(x=as.character(inputValue), fill=as.character(propertyValue))) + geom_bar(position = 'fill')

#uniprot analysis
uniprot_result <- read.csv('Dataset/metascape/kub_kxub/uniprot/uniprot.tab', sep='\t')
uniprot_df <- read.csv('Dataset/metascape/kub_kxub/kub_kxub_protein_for_uniprot.csv')
uniprot_df <- dplyr::left_join(uniprot_df,uniprot_result,by = c('ACCID' = 'Entry'))
uniprot_df <- uniprot_df[!is.na(uniprot_df[,3]),]
write.csv(uniprot_df,'Dataset/metascape/kub_kxub/uniprot/uniprot_result.csv')

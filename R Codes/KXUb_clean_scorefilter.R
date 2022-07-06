library(ggplot2)

#read combined KXUb dataset
kxub_site <- read.csv('Dataset/All_combined_remove_false_assignment.csv')

#filter site because of Evidence ID missing
kxub_site <- kxub_site[kxub_site$Evidence.IDs != '',]
write.csv(kxub_site,'Dataset/All_combined_clean.csv',row.names = FALSE)
kxub_site <- read.csv('Dataset/All_combined_clean.csv')

#read kub dataset
kub_site <- read.csv('Dataset/GlyGly (K)Sites_op.csv')
kub_site <- kub_site[kub_site$Evidence.IDs != '',]
write.csv(kub_site,'Dataset/All_GG_clean.csv',row.names = FALSE)
kub_site <- read.csv('Dataset/All_GG_clean.csv')

#score filter by 60
kxub_site <- kxub_site[kxub_site$Score > 60,]
write.csv(kxub_site,'Dataset/All_combined_score_filter.csv',row.names = FALSE)
kxub_site <- read.csv('Dataset/All_combined_score_filter.csv')
nrow(kxub_site)

kub_site <- kub_site[kub_site$Score > 60,]
write.csv(kub_site,'Dataset/All_GG_score_filter.csv',row.names = FALSE)
kub_site <- read.csv('Dataset/All_GG_score_filter.csv')
nrow(kub_site)

#kxub ub intersection label
kub_s <- kub_site[,c('ACCID','Position')]
kxub_site_intersect <- kxub_site

kub_s$Site.in.Ub <- 1
kxub_site_intersect <- dplyr::left_join(kxub_site_intersect, kub_s, by=c('ACCID','Position'))

kub_s <- kub_s[,c('ACCID','Position')]
kub_s$Protein.in.Ub <- 1
kub_s <- aggregate(Protein.in.Ub~ACCID, kub_s, max)
kxub_site_intersect <- dplyr::left_join(kxub_site_intersect, kub_s, by='ACCID')

write.csv(kxub_site_intersect, 'Dataset/All_combined_intersetion.csv', row.names = FALSE)
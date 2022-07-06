library(ggplot2)

#read files
#read meta file
meta <- read.csv('meta.csv')
#read kxub and kub dataset after score filtering
kxub_site <- read.csv('Dataset/All_combined_score_filter.csv')
kub_site <- read.csv('Dataset/All_GG_score_filter.csv')

#score, delta score and intensity distribution
ggplot(kxub_site,aes(x=Score)) + geom_histogram(binwidth = 5, color='black',fill='#F8766D')+theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Score')+ylab('Site Count')
ggplot(kxub_site,aes(x=Delta.score)) + geom_histogram(binwidth = 5, color='black',fill='#F8766D')+theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Delta score')+ylab('Site Count')
ggplot(kxub_site,aes(x=log10(Intensity))) + geom_histogram(binwidth = 0.1, color='black',fill='#F8766D')+theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Log10(Intensity)')+ylab('Site Count')

median(kxub_site$Score)
median(kxub_site$Delta.score)
median(log10(kxub_site$Intensity))

#characteristic fragment plots
kxub_site$Characteristic.Fragment[is.na(kxub_site$Characteristic.Fragment)] <- 0
ggplot(kxub_site,aes(x=reorder(Mod,table(Mod)[Mod]),fill=as.character(Characteristic.Fragment),))+geom_bar()+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Modifications')+ylab('Count')

#site and protein identified number summary
kxub_site_iden <- kxub_site[,meta$Count]
colnames(kxub_site_iden) <- meta$Sample
kxub_site_iden[kxub_site_iden != ''] <- 1
kxub_site_iden[kxub_site_iden == ''] <- 0
kxub_site_iden <- data.frame(sapply(kxub_site_iden, function(x) as.numeric(as.character(x))))
kxub_site_iden$Total <- rowSums(kxub_site_iden[,1:18])
kxub_site_iden$Hep2_Total <- rowSums(kxub_site_iden[,1:9])
kxub_site_iden$Hep2_CTR <- rowSums(kxub_site_iden[,1:3])
kxub_site_iden$Hep2_Bort <- rowSums(kxub_site_iden[,4:6])
kxub_site_iden$Hep2_AP15 <- rowSums(kxub_site_iden[,7:9])
kxub_site_iden$Jurkat_Total <- rowSums(kxub_site_iden[,10:18])
kxub_site_iden$Jurkat_CTR <- rowSums(kxub_site_iden[,10:12])
kxub_site_iden$Jurkat_Bort <- rowSums(kxub_site_iden[,13:15])
kxub_site_iden$Jurkat_AP15 <- rowSums(kxub_site_iden[,16:18])

kxub_site_iden[kxub_site_iden > 0] <- 1

kxub_prot_iden <- cbind(kxub_site$ACCID,kxub_site_iden)
colnames(kxub_prot_iden) <- c('ACCID',colnames(kxub_prot_iden)[2:28])
kxub_prot_iden <- aggregate(cbind(Total,Hep2_Total,Hep2_CTR,Hep2_Bort,Hep2_AP15,Jurkat_Total,Jurkat_CTR,Jurkat_Bort,Jurkat_AP15)~ACCID, kxub_prot_iden, sum)
kxub_prot_iden <- kxub_prot_iden[,2:10]
kxub_prot_iden[kxub_prot_iden > 0] <- 1

kxub_site_iden <- kxub_site_iden[,19:27]

kxub_site_iden_cs <- colSums(kxub_site_iden)
kxub_prot_iden_cs <- colSums(kxub_prot_iden)

kxub_s_p_iden <- cbind(kxub_site_iden_cs,kxub_prot_iden_cs)
kxub_s_p_iden <- data.frame(kxub_s_p_iden)
colnames(kxub_s_p_iden) <- c('Sites','Proteins')
kxub_s_p_iden$Condition <- rownames(kxub_s_p_iden)
kxub_s_p_iden <- reshape2::melt(kxub_s_p_iden)

kxub_s_p_iden$Condition <- factor(kxub_s_p_iden$Condition, levels=c('Total','Hep2_Total','Hep2_CTR','Hep2_Bort','Hep2_AP15','Jurkat_Total','Jurkat_CTR','Jurkat_Bort','Jurkat_AP15'))
ggplot(kxub_s_p_iden,aes(x=Condition,y=value,group=variable,fill=variable))+
  geom_bar(stat = 'identity',colour='black',position='dodge')+
  theme_bw()+theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('')+ylab('Count')+geom_text(mapping = aes(label=value),vjust=-1,position=position_dodge(width=1),size=5)+ylim(0,2200)+
  theme(legend.position = c(.8,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_blank(),legend.text = element_text(size=18))

#compare identification in Hep2 and Jurkat
site_cp <- kxub_site_iden[,c('Hep2_Total','Jurkat_Total')]
prot_cp <- kxub_prot_iden[,c('Hep2_Total','Jurkat_Total')]
site_cp$Jurkat_Total[site_cp$Jurkat_Total == 1] <- 10
prot_cp$Jurkat_Total[prot_cp$Jurkat_Total == 1] <- 10
site_cp$Compare <- rowSums(site_cp)
prot_cp$Compare <- rowSums(prot_cp)
table(site_cp$Compare)
table(prot_cp$Compare)

site_cp_op <- cbind(kxub_site[,c('ACCID','PROTEINS','Position')], kxub_site_iden[,c('Hep2_Total','Jurkat_Total')])
prot_cp_op <- cbind(aggregate(Hep2_Total~ACCID+PROTEINS,site_cp_op,sum),aggregate(Jurkat_Total~ACCID+PROTEINS,site_cp_op,sum))
prot_cp_op$Hep2_Total[prot_cp_op$Hep2_Total > 0] <- 1
prot_cp_op$Jurkat_Total[prot_cp_op$Jurkat_Total > 0] <- 1

write.csv(site_cp_op,'Dataset/site_compare_Hep2_Jurkat.csv',row.names = FALSE)
write.csv(prot_cp_op,'Dataset/protein_compare_Hep2_Jurkat.csv',row.names = FALSE)

#modification num on proteins
mod_num <- kxub_site[,c('ACCID','Mod')]
mod_num$Sum <- 1
mod_num <- aggregate(Sum~ACCID,mod_num,sum)
mod_num$Sum[mod_num$Sum > 5] <- 6
mod_num <- data.frame(table(mod_num$Sum))
colnames(mod_num) <- c('Term','Num')
mod_num$Percentage <- (mod_num$Num / sum(mod_num$Num))

ggplot(mod_num,aes(x='content',y=Percentage,fill=Term))+
  geom_bar(stat = 'identity',position = 'stack')+coord_polar(theta = 'y')+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
  theme(legend.title = element_text(size=20),legend.text = element_text(size=18,hjust=1))+scale_fill_brewer(palette='Accent')+labs(fill='Modifications per Protein')

modtype_num <- kxub_site[,c('ACCID','Mod')]
modtype_num$Sum <- 1
modtype_num <- aggregate(Sum~ACCID+Mod,modtype_num,sum)
modtype_num$Sum[modtype_num$Sum > 3] <- 4
modtype_num <- data.frame(table(modtype_num$Sum))
colnames(modtype_num) <- c('Term','Num')
modtype_num$Percentage <- (modtype_num$Num / sum(modtype_num$Num))

ggplot(modtype_num,aes(x='content',y=Percentage,fill=Term))+
  geom_bar(stat = 'identity',position = 'stack')+coord_polar(theta = 'y')+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
  theme(legend.title = element_text(size=20),legend.text = element_text(size=18,hjust=1))+scale_fill_brewer(palette='Accent')+labs(fill='K-XUb types per Protein')

#compare characteristic fragment of kxub and kub
kub_site$Characteristic.Fragment[is.na(kub_site$Characteristic.Fragment)] <- 0
kxub_site_char <- kxub_site[c('Mod','Characteristic.Fragment')]
kxub_site_char$Mod <- 'XGG'
kub_site_char <- kub_site[c('Mod','Characteristic.Fragment')]
site_char <- rbind(kxub_site_char,kub_site_char)

ggplot(site_char,aes(x=reorder(Mod,table(Mod)[Mod]),fill=as.character(Characteristic.Fragment),))+geom_bar(position='fill')+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Modifications')+ylab('Percentage')

#compare score and intensity distribution
kxub <- kxub_site[,c('ACCID','Position','Score','Intensity')]
kxub$Class <- 'kxub'
kub <- kub_site[,c('ACCID','Position','Score','Intensity')]
kub$Class <- 'kub'
compare_temp <- rbind(kxub,kub)

ggplot(compare_temp,aes(x=Score)) + geom_density(data=subset(compare_temp, Class == 'kxub'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(compare_temp, Class == 'kub'), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,350)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Score')+ylab('Density')

ggplot(compare_temp,aes(x=log10(Intensity))) + geom_density(data=subset(compare_temp, Class == 'kxub'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(compare_temp, Class == 'kub'), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+xlim(5,12)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Log10(Intensity)')+ylab('Density')

#venn of kxub and kub
kxub <- kxub_site[,c('ACCID','Position')]
kxub$Class <- 10
kub <- kub_site[,c('ACCID','Position')]
kub$Class <- 1
compare_temp <- rbind(kxub,kub)
compare_temp <- aggregate(Class~ACCID+Position,compare_temp,sum)
table(compare_temp$Class)

kxub <- kxub_site[,c('ACCID','Position')]
kxub$Class <- 10
kxub <- aggregate(Class~ACCID,kxub,sum)
kxub$Class <- 10
kub <- kub_site[,c('ACCID','Position')]
kub$Class <- 1
kub <- aggregate(Class~ACCID,kub,sum)
kub$Class <- 1
compare_temp <- rbind(kxub,kub)
compare_temp <- aggregate(Class~ACCID,compare_temp,sum)
table(compare_temp$Class)

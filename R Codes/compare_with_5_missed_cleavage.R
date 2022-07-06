library(ggplot2)

#read 5 missed cleavage dataset
kxub_5_miss_site <- read.csv('Dataset/All_combined_5missed.csv')

#dataset process
kxub_5_miss_site_pr <- kxub_5_miss_site[,c('Score','Delta.score','Intensity')]
kxub_5_miss_site_pr$Class <- '5 missed'
kxub_site_pr <- kxub_site[,c('Score','Delta.score','Intensity')]
kxub_site_pr$Class <- '3 missed'

df_cp_pr <- rbind(kxub_5_miss_site_pr,kxub_site_pr)

#score, delta score and intensity distribution compare
ggplot(df_cp_pr,aes(x=Score)) + geom_density(data=subset(df_cp_pr, Class == '5 missed'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(df_cp_pr, Class == '3 missed'), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,300)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Score')+ylab('Density')

ggplot(df_cp_pr,aes(x=Delta.score)) + geom_density(data=subset(df_cp_pr, Class == '5 missed'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(df_cp_pr, Class == '3 missed'), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,200)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Delta Score')+ylab('Density')

ggplot(df_cp_pr,aes(x=log10(Intensity))) + geom_density(data=subset(df_cp_pr, Class == '5 missed'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(df_cp_pr, Class == '3 missed'), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Log10(Intensity)')+ylab('Density')

median(subset(df_cp_pr, Class == '3 missed')$Score)
median(subset(df_cp_pr, Class == '5 missed')$Score)
median(subset(df_cp_pr, Class == '3 missed')$Delta.score)
median(subset(df_cp_pr, Class == '5 missed')$Delta.score)

#compare intersection
site_3 <- kxub_site[,c('ACCID','Position','Mod')]
site_3[site_3 == 'CamCGG'] <- 'CGG'
site_3$Class <- 1
site_5 <- kxub_5_miss_site[,c('ACCID','Position','Mod')]
site_5$Class <- 10
cp_site <- rbind(site_3,site_5)
cp_site <- aggregate(Class~ACCID+Position+Mod,cp_site,sum)
table(cp_site$Class)

prot_3 <- kxub_site[,c('ACCID','Mod')]
prot_3$Class <- 1
prot_3 <- aggregate(Class~ACCID,prot_3,sum)
prot_3$Class <- 1
prot_5 <- kxub_5_miss_site[,c('ACCID','Mod')]
prot_5$Class <- 10
prot_5 <- aggregate(Class~ACCID,prot_5,sum)
prot_5$Class <- 10
cp_prot <- rbind(prot_3,prot_5)
cp_prot <- aggregate(Class~ACCID,cp_prot,sum)
table(cp_prot$Class)

#intersected sites
site_3_inter <- kxub_site[,c('ACCID','Position','Mod','Score','Delta.score','Intensity','Characteristic.Fragment')]
site_3_inter[site_3_inter == 'CamCGG'] <- 'CGG'
site_3_inter <- dplyr::left_join(site_3_inter,site_5,by=c('ACCID','Position','Mod'))
site_3_inter$Class[site_3_inter$Class == 10] <- 1
site_3_inter$Class[is.na(site_3_inter$Class)] <- 0

ggplot(site_3_inter,aes(x=Score)) + geom_density(data=subset(site_3_inter, Class == 0), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(site_3_inter, Class == 1), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,300)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Score')+ylab('Density')

ggplot(site_3_inter,aes(x=Delta.score)) + geom_density(data=subset(site_3_inter, Class == 0), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(site_3_inter, Class == 1), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,200)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Delta Score')+ylab('Density')

ggplot(site_3_inter,aes(x=log10(Intensity))) + geom_density(data=subset(site_3_inter, Class == 0), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(site_3_inter, Class == 1), color='#F8766D', size=1.5, alpha=0.5)+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Log10(Intensity)')+ylab('Density')

ggplot(site_3_inter,aes(x=Class,fill=as.character(Characteristic.Fragment),))+geom_bar(position='fill')+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Modifications')+ylab('Percentage')

median(subset(site_3_inter, Class == 0)$Score)
median(subset(site_3_inter, Class == 1)$Score)
median(subset(site_3_inter, Class == 0)$Delta.score)
median(subset(site_3_inter, Class == 1)$Delta.score)

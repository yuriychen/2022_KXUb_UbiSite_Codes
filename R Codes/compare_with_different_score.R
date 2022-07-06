library(ggplot2)

#copy dataset
kxub_site_sc <- kxub_site
kxub_site_sc$Count <- 1

#label different score range
kxub_site_sc$scorelabel <- ''
kxub_site_sc$scorelabel[(kxub_site_sc$Score >= 40 & kxub_site_sc$Score <= 50)] <- '40-50'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 50 & kxub_site_sc$Score <= 60)] <- '50-60'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 60 & kxub_site_sc$Score <= 70)] <- '60-70'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 70 & kxub_site_sc$Score <= 80)] <- '70-80'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 80 & kxub_site_sc$Score <= 90)] <- '80-90'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 90 & kxub_site_sc$Score <= 100)] <- '90-100'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 100)] <- '100-'
kxub_site_sc$scorelabel <- factor(kxub_site_sc$scorelabel, levels=c('40-50','50-60','60-70','70-80','80-90','90-100','100-'))

#plot characteristic fragment
ggplot(kxub_site_sc,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar(position='fill')+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

ggplot(kxub_site_sc,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar()+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

#detailed score range
kxub_site_sc$scorelabel <- ''
kxub_site_sc$scorelabel[(kxub_site_sc$Score >= 40 & kxub_site_sc$Score <= 45)] <- '40-45'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 45 & kxub_site_sc$Score <= 50)] <- '45-50'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 50 & kxub_site_sc$Score <= 55)] <- '50-55'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 55 & kxub_site_sc$Score <= 60)] <- '55-60'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 60 & kxub_site_sc$Score <= 65)] <- '60-65'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 65 & kxub_site_sc$Score <= 70)] <- '65-70'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 70 & kxub_site_sc$Score <= 75)] <- '70-75'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 75 & kxub_site_sc$Score <= 80)] <- '75-80'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 80 & kxub_site_sc$Score <= 85)] <- '80-85'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 85 & kxub_site_sc$Score <= 90)] <- '85-90'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 90 & kxub_site_sc$Score <= 95)] <- '90-95'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 95 & kxub_site_sc$Score <= 100)] <- '95-100'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 100)] <- '100-'
kxub_site_sc$scorelabel <- factor(kxub_site_sc$scorelabel, levels=c('40-45','45-50','50-55','55-60','60-65','65-70','70-75','75-80','80-85','85-90','90-95','95-100','100-'))

ggplot(kxub_site_sc,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar(position='fill')+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

ggplot(kxub_site_sc,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar()+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

#score and intensity distribution and characteristic fragment of score filtered results
score_40 <- kxub_site_sc[kxub_site_sc$Score >= 40, c('Score','Delta.score','Intensity','Characteristic.Fragment')]
score_40$scorelabel <- '40-'
score_60 <- kxub_site_sc[kxub_site_sc$Score > 60, c('Score','Delta.score','Intensity','Characteristic.Fragment')]
score_60$scorelabel <- '60-'
score_80 <- kxub_site_sc[kxub_site_sc$Score > 80, c('Score','Delta.score','Intensity','Characteristic.Fragment')]
score_80$scorelabel <- '80-'
score_100 <- kxub_site_sc[kxub_site_sc$Score > 100, c('Score','Delta.score','Intensity','Characteristic.Fragment')]
score_100$scorelabel <- '100-'

score_df <- rbind(score_40,score_60,score_80,score_100)
score_df$scorelabel <- factor(score_df$scorelabel, levels=c('40-','60-','80-','100-'))

ggplot(score_df,aes(x=Score)) + geom_density(data=subset(score_df, scorelabel == '40-'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(score_df, scorelabel == '60-'), color='#F8766D', size=1.5, alpha=0.5)+
  geom_density(data=subset(score_df, scorelabel == '80-'), color='red', size=1.5, alpha=0.5)+
  geom_density(data=subset(score_df, scorelabel == '100-'), color='darkred', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,300)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Score')+ylab('Density')

ggplot(score_df,aes(x=Delta.score)) + geom_density(data=subset(score_df, scorelabel == '40-'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(score_df, scorelabel == '60-'), color='#F8766D', size=1.5, alpha=0.5)+
  geom_density(data=subset(score_df, scorelabel == '80-'), color='red', size=1.5, alpha=0.5)+
  geom_density(data=subset(score_df, scorelabel == '100-'), color='darkred', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,200)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Delta Score')+ylab('Density')

ggplot(score_df,aes(x=log10(Intensity))) + geom_density(data=subset(score_df, scorelabel == '40-'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(score_df, scorelabel == '60-'), color='#F8766D', size=1.5, alpha=0.5)+
  geom_density(data=subset(score_df, scorelabel == '80-'), color='red', size=1.5, alpha=0.5)+
  geom_density(data=subset(score_df, scorelabel == '100-'), color='darkred', size=1.5, alpha=0.5)+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Log10(Intensity)')+ylab('Density')

ggplot(score_df,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar(position='fill')+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

ggplot(score_df,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar()+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

#score seperate
kxub_site_sc$scorelabel <- ''
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 40 & kxub_site_sc$PEP <= 50)] <- '40-50'
kxub_site_sc$scorelabel[(kxub_site_sc$Score > 50 & kxub_site_sc$PEP <= 60)] <- '50-60'
kxub_site_sc$scorelabel[( kxub_site_sc$Score > 60)] <- '60-'
kxub_site_sc$scorelabel <- factor(kxub_site_sc$scorelabel, levels=c('40-50','50-60','60-'))

ggplot(kxub_site_sc,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar(position='fill')+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

ggplot(kxub_site_sc,aes(x=scorelabel,fill=as.character(Characteristic.Fragment)))+geom_bar()+theme_bw()+
  theme(legend.position = c(.2,.8),legend.background = element_rect(size = 1,color='black'),legend.title = element_text(size=18),legend.text = element_text(size=18))+
  scale_fill_manual('Characteristic Fragment',labels=c('Without','With'),values=c('darkgray','#F8766D'))+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18,angle = 60,hjust = 1))+
  xlab('Score range')+ylab('Percentage')

ggplot(kxub_site_sc,aes(x=Score)) + geom_density(data=subset(kxub_site_sc, scorelabel == '40-'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(kxub_site_sc, scorelabel == '-60'), color='#F8766D', size=1.5, alpha=0.5)+
  geom_density(data=subset(kxub_site_sc, scorelabel == '60-'), color='darkred', size=1.5, alpha=0.5)+
  theme_bw()+xlim(0,300)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Score')+ylab('Density')

ggplot(kxub_site_sc,aes(x=log10(Intensity))) + geom_density(data=subset(kxub_site_sc, scorelabel == '40-'), color='gray', size=1.5, alpha=0.5) + 
  geom_density(data=subset(kxub_site_sc, scorelabel == '-60'), color='#F8766D', size=1.5, alpha=0.5)+
  geom_density(data=subset(kxub_site_sc, scorelabel == '60-'), color='darkred', size=1.5, alpha=0.5)+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('Log10(Intensity)')+ylab('Density')

ggplot(kxub_site_sc,aes(x=scorelabel,y=Score, fill=scorelabel)) + 
  geom_boxplot()+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('')+ylab('Score')

ggplot(kxub_site_sc,aes(x=scorelabel,y=Delta.score, fill=scorelabel)) + 
  geom_boxplot()+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('')+ylab('Delta score')

ggplot(kxub_site_sc,aes(x=scorelabel,y=log10(Intensity), fill=scorelabel)) + 
  geom_boxplot()+
  theme_bw()+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('')+ylab('Log10(Intensity)')

ggplot(kxub_site_sc,aes(x=scorelabel,y=-log10(PEP), fill=scorelabel)) + 
  geom_boxplot()+
  theme_bw()+ylim(0,10)+
  theme(axis.title = element_text(size=20),axis.text.y = element_text(size=18),axis.text.x = element_text(size=18))+
  xlab('')+ylab('-Log10(PEP)')

nrow(subset(kxub_site_sc, scorelabel == '60-'))
     

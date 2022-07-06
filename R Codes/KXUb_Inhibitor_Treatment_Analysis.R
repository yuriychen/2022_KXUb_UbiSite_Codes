library(ggplot2)
library(pheatmap)
library(limma)

#inhibitor treatment unique
kxub_site_inhibitor <- cbind(kxub_site[,c('ACCID','Position','PROTEINS','Mod')],kxub_site_iden[,2:9])
kxub_site_inhibitor$Hep2_inhibitor <- kxub_site_inhibitor$Hep2_Bort + kxub_site_inhibitor$Hep2_AP15
kxub_site_inhibitor$Hep2_inhibitor[kxub_site_inhibitor$Hep2_inhibitor > 0] <- 10
kxub_site_inhibitor$Jurkat_inhibitor <- kxub_site_inhibitor$Jurkat_Bort + kxub_site_inhibitor$Jurkat_AP15
kxub_site_inhibitor$Jurkat_inhibitor[kxub_site_inhibitor$Jurkat_inhibitor > 0] <- 10

kxub_site_inhibitor$Hep2_inhibitor_ctr <- kxub_site_inhibitor$Hep2_CTR + kxub_site_inhibitor$Hep2_inhibitor
kxub_site_inhibitor$Jurkat_inhibitor_ctr <- kxub_site_inhibitor$Jurkat_CTR + kxub_site_inhibitor$Jurkat_inhibitor

table(kxub_site_inhibitor$Hep2_inhibitor_ctr)
table(kxub_site_inhibitor$Jurkat_inhibitor_ctr)

write.csv(kxub_site_inhibitor,'Dataset/Quantification/hep2_jurkat_inhibitor_ctr_unique.csv', row.names = FALSE)

#function for sam analysis
data_calsmoothcurve = function(x,ta,s0,df){
  if (x > 0) {
    y = ta * (1 + (s0 / ((x/ta)-s0)))
    y = -log10(2*(1-pt(y,df=df)))
    return(y)
  }
  else if (x < 0) {
    y = ta * (1 + (s0 / ((x/-ta)-s0)))
    y = -log10(2*(1-pt(y,df=df)))
    return(y)
  }
}

data_calsmoothcurve_func = function(x,ta,s0,df){
  y <- ifelse(x > 0, (ta * (1 + (s0 / ((x/ta)-s0)))), (ta * (1 + (s0 / ((x/-ta)-s0)))))
  y <- -log10(2*(1-pt(y,df=df)))
  return(y)
}

#volcano plot of hep2 ctr bort ctr ap15 with not unique protein
kxub_quant_hep2 <- kxub_site[kxub_site_inhibitor$Hep2_inhibitor_ctr == 11, c('ACCID','Position','PROTEINS','Mod', meta$Intensity[meta$CellLine == 'Hep2'])]
hep2_bort <- kxub_quant_hep2[,5:10]
hep2_bort[hep2_bort > 0] <- 1
hep2_bort$count <- rowSums(hep2_bort)
hep2_bort <- kxub_quant_hep2[hep2_bort$count > 1,]
df_temp <- hep2_bort[,c(8,9,10,5,6,7)]
df_temp[df_temp == 0] <- (0.5 * min(df_temp[df_temp != min(df_temp)]))
hep2_bort <- cbind(hep2_bort[,1:4],df_temp)

hep2_ap15 <- kxub_quant_hep2[,c(5,6,7,11,12,13)]
hep2_ap15[hep2_ap15 > 0] <- 1
hep2_ap15$count <- rowSums(hep2_ap15)
hep2_ap15 <- kxub_quant_hep2[hep2_ap15$count > 1,]
df_temp <- hep2_ap15[,c(11,12,13,5,6,7)]
df_temp[df_temp == 0] <- (0.5 * min(df_temp[df_temp != min(df_temp)]))
hep2_ap15 <- cbind(hep2_ap15[,1:4],df_temp)

df <- 5
ta <- qt(0.925,df)
s0 <- 0.1

hep2_bort$fc <- rowMeans(hep2_bort[,5:7]) / rowMeans(hep2_bort[,8:10])
hep2_bort$pvalue <- apply(hep2_bort[5:10], 1, function(x){t.test(x[1:3], x[4:5])[["p.value"]]})
hep2_bort$fc <- log2(hep2_bort$fc)
hep2_bort$fc[hep2_bort$fc == 0] <- 0.000001
hep2_bort$pvalue <- -log10(hep2_bort$pvalue)
hep2_bort$fudge <- apply(hep2_bort, 1, function(x){data_calsmoothcurve(as.numeric(x[11]),ta,s0,df)})
hep2_bort$sig <- '0'
hep2_bort$sig <- apply(hep2_bort,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) > ta*s0)),'1',x[14])})
hep2_bort$sig <- apply(hep2_bort,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) < -ta*s0)),'-1',x[14])})
hep2_bort$inhibitor <- 'Bort'

write.csv(hep2_bort,'Dataset/Quantification/hep2_bort_volcano.csv', row.names = FALSE)

hep2_ap15$fc <- rowMeans(hep2_ap15[,5:7]) / rowMeans(hep2_ap15[,8:10])
hep2_ap15$pvalue <- apply(hep2_ap15[5:10], 1, function(x){t.test(x[1:3], x[4:5])[["p.value"]]})
hep2_ap15$fc <- log2(hep2_ap15$fc)
hep2_ap15$fc[hep2_ap15$fc == 0] <- 0.000001
hep2_ap15$pvalue <- -log10(hep2_ap15$pvalue)
hep2_ap15$fudge <- apply(hep2_ap15, 1, function(x){data_calsmoothcurve(as.numeric(x[11]),ta,s0,df)})
hep2_ap15$sig <- '0'
hep2_ap15$sig <- apply(hep2_ap15,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) > ta*s0)),'1',x[14])})
hep2_ap15$sig <- apply(hep2_ap15,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) < -ta*s0)),'-1',x[14])})
hep2_ap15$inhibitor <- 'AP15'

write.csv(hep2_ap15,'Dataset/Quantification/hep2_ap15_volcano.csv', row.names = FALSE)

hep2_q <- rbind(hep2_bort[,c(1,2,3,4,11,12,13,14,15)],hep2_ap15[,c(1,2,3,4,11,12,13,14,15)])

ggplot(hep2_q,aes(x=fc,y=pvalue)) + geom_point(aes(color=sig, shape=inhibitor), size=2) + ylim(0,4.55) + xlim(-15,15) +
  scale_color_manual(values = c('green','grey60','red')) + scale_shape_manual(values= c(16,17)) +
  xlab('Log2(Inhibitor/Control)') + ylab('-Log10(P-Value)') +
  geom_function(fun=function(x){data_calsmoothcurve_func(x,ta,s0,df)},xlim=c(ta*s0,17),size=.75,alpha=.3,color='darkred',n=10000,linetype='dashed') +
  geom_function(fun=function(x){data_calsmoothcurve_func(x,ta,s0,df)},xlim=c(-17,-ta*s0),size=.75,alpha=.3,color='darkgreen',n=10000,linetype='dashed') +
  theme_bw()

#volcano plot of jurkat ctr bort ctr ap15
kxub_quant_jurkat <- kxub_site[kxub_site_inhibitor$Jurkat_inhibitor_ctr == 11, c('ACCID','Position','PROTEINS','Mod', meta$Intensity[meta$CellLine == 'Jurkat'])]
jurkat_bort <- kxub_quant_jurkat[,5:10]
jurkat_bort[jurkat_bort > 0] <- 1
jurkat_bort$count <- rowSums(jurkat_bort)
jurkat_bort <- kxub_quant_jurkat[jurkat_bort$count > 1,]
df_temp <- jurkat_bort[,c(8,9,10,5,6,7)]
df_temp[df_temp == 0] <- (0.5 * min(df_temp[df_temp != min(df_temp)]))
jurkat_bort <- cbind(jurkat_bort[,1:4],df_temp)

jurkat_ap15 <- kxub_quant_jurkat[,c(5,6,7,11,12,13)]
jurkat_ap15[jurkat_ap15 > 0] <- 1
jurkat_ap15$count <- rowSums(jurkat_ap15)
jurkat_ap15 <- kxub_quant_jurkat[jurkat_ap15$count > 1,]
df_temp <- jurkat_ap15[,c(11,12,13,5,6,7)]
df_temp[df_temp == 0] <- (0.5 * min(df_temp[df_temp != min(df_temp)]))
jurkat_ap15 <- cbind(jurkat_ap15[,1:4],df_temp)

df <- 5
ta <- qt(0.925,df)
s0 <- 0.1

jurkat_bort$fc <- rowMeans(jurkat_bort[,5:7]) / rowMeans(jurkat_bort[,8:10])
jurkat_bort$pvalue <- apply(jurkat_bort[5:10], 1, function(x){t.test(x[1:3], x[4:5])[["p.value"]]})
jurkat_bort$fc <- log2(jurkat_bort$fc)
jurkat_bort$fc[jurkat_bort$fc == 0] <- 0.000001
jurkat_bort$pvalue <- -log10(jurkat_bort$pvalue)
jurkat_bort$fudge <- apply(jurkat_bort, 1, function(x){data_calsmoothcurve(as.numeric(x[11]),ta,s0,df)})
jurkat_bort$sig <- '0'
jurkat_bort$sig <- apply(jurkat_bort,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) > ta*s0)),'1',x[14])})
jurkat_bort$sig <- apply(jurkat_bort,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) < -ta*s0)),'-1',x[14])})
jurkat_bort$inhibitor <- 'Bort'

write.csv(jurkat_bort,'Dataset/Quantification/jurkat_bort_volcano.csv', row.names = FALSE)

jurkat_ap15$fc <- rowMeans(jurkat_ap15[,5:7]) / rowMeans(jurkat_ap15[,8:10])
jurkat_ap15$pvalue <- apply(jurkat_ap15[5:10], 1, function(x){t.test(x[1:3], x[4:5])[["p.value"]]})
jurkat_ap15$fc <- log2(jurkat_ap15$fc)
jurkat_ap15$fc[jurkat_ap15$fc == 0] <- 0.000001
jurkat_ap15$pvalue <- -log10(jurkat_ap15$pvalue)
jurkat_ap15$fudge <- apply(jurkat_ap15, 1, function(x){data_calsmoothcurve(as.numeric(x[11]),ta,s0,df)})
jurkat_ap15$sig <- '0'
jurkat_ap15$sig <- apply(jurkat_ap15,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) > ta*s0)),'1',x[14])})
jurkat_ap15$sig <- apply(jurkat_ap15,1,function(x){ifelse(((as.numeric(x[12]) > as.numeric(x[13])) & ((as.numeric(x[11])) < -ta*s0)),'-1',x[14])})
jurkat_ap15$inhibitor <- 'AP15'

write.csv(jurkat_ap15,'Dataset/Quantification/jurkat_ap15_volcano.csv', row.names = FALSE)

jurkat_q <- rbind(jurkat_bort[,c(1,2,3,4,11,12,13,14,15)],jurkat_ap15[,c(1,2,3,4,11,12,13,14,15)])

ggplot(jurkat_q,aes(x=fc,y=pvalue)) + geom_point(aes(color=sig, shape=inhibitor), size=2) + ylim(0,3.5) + xlim(-13,13) +
  scale_color_manual(values = c('green','grey60','red')) + scale_shape_manual(values= c(16,17)) +
  xlab('Log2(Inhibitor/Control)') + ylab('-Log10(P-Value)') +
  geom_function(fun=function(x){data_calsmoothcurve_func(x,ta,s0,df)},xlim=c(ta*s0,15),size=.75,alpha=.3,color='darkred',n=10000,linetype='dashed') +
  geom_function(fun=function(x){data_calsmoothcurve_func(x,ta,s0,df)},xlim=c(-15,-ta*s0),size=.75,alpha=.3,color='darkgreen',n=10000,linetype='dashed') +
  theme_bw()

#combine inhibitor unique and enriched sites
hep2_inhibitor_up <- rbind(kxub_site_inhibitor[kxub_site_inhibitor$Hep2_inhibitor_ctr == 10, c('ACCID','Position','PROTEINS','Mod')], hep2_q[hep2_q$sig == '1', c('ACCID','Position','PROTEINS','Mod')])
hep2_inhibitor_down <- rbind(kxub_site_inhibitor[kxub_site_inhibitor$Hep2_inhibitor_ctr == 1, c('ACCID','Position','PROTEINS','Mod')], hep2_q[hep2_q$sig == '-1', c('ACCID','Position','PROTEINS','Mod')])

jurkat_inhibitor_up <- rbind(kxub_site_inhibitor[kxub_site_inhibitor$Jurkat_inhibitor_ctr == 10, c('ACCID','Position','PROTEINS','Mod')], jurkat_q[jurkat_q$sig == '1', c('ACCID','Position','PROTEINS','Mod')])
jurkat_inhibitor_down <- rbind(kxub_site_inhibitor[kxub_site_inhibitor$Jurkat_inhibitor_ctr == 1, c('ACCID','Position','PROTEINS','Mod')], jurkat_q[jurkat_q$sig == '-1', c('ACCID','Position','PROTEINS','Mod')])

hep2_inhibitor_up$Enrichment <- 'Inhibitor'
hep2_inhibitor_down$Enrichment <- 'CTR'
jurkat_inhibitor_up$Enrichment <- 'Inhibitor'
jurkat_inhibitor_down$Enrichment <- 'CTR'

hep2_inhibitor <- rbind(hep2_inhibitor_up, hep2_inhibitor_down)
jurkat_inhibitor <- rbind(jurkat_inhibitor_up, jurkat_inhibitor_down)

write.csv(hep2_inhibitor, 'Dataset/Quantification/hep2_inhibitor.csv', row.names = FALSE)
write.csv(jurkat_inhibitor, 'Dataset/Quantification/jurkat_inhibitor.csv', row.names = FALSE)

hep2_inhibitor$inhibitors[hep2_inhibitor$Enrichment == 'Inhibitor'] <- 1
hep2_inhibitor$inhibitors[hep2_inhibitor$Enrichment == 'CTR'] <- -1
hep2_inhibitor <- aggregate(hep2_inhibitor,inhibitors~ACCID+Enrichment,max)
hep2_inhibitor <- aggregate(hep2_inhibitor,inhibitors~ACCID,sum)
write.csv(hep2_inhibitor,'Dataset/Quantification/hep2_inhibitor_for_pathview.csv', row.names = FALSE)

#enriched sites after inhibitor treatment overlap in hep2 and jurkat
hep2_q1 <- hep2_inhibitor_up
hep2_q1$cell <- 1
hep2_q1 <- aggregate(cell~ACCID+PROTEINS+Position+Mod, hep2_q1, sum)
hep2_q1$cell <- 1

jurkat_q1 <- jurkat_inhibitor_up
jurkat_q1$cell <- 10
jurkat_q1 <- aggregate(cell~ACCID+PROTEINS+Position+Mod, jurkat_q1, sum)
jurkat_q1$cell <- 10

q_compare <- rbind(hep2_q1,jurkat_q1)
q_compare <- aggregate(cell~ACCID+PROTEINS+Position+Mod, q_compare, sum)
table(q_compare$cell)
write.csv(q_compare, 'Dataset/Quantification/hep2_jurkat_site_inhibitor_enriched.csv', row.names = FALSE)

hep2_q1 <- aggregate(cell~ACCID+PROTEINS, hep2_q1, sum)
hep2_q1$cell <- 1
jurkat_q1 <- aggregate(cell~ACCID+PROTEINS, jurkat_q1, sum)
jurkat_q1$cell <- 10

q_compare <- rbind(hep2_q1,jurkat_q1)
q_compare <- aggregate(cell~ACCID+PROTEINS, q_compare, sum)
table(q_compare$cell)
write.csv(q_compare, 'Dataset/Quantification/hep2_jurkat_protein_inhibitor_enriched.csv', row.names = FALSE)

#heatmap of overall sites with quantification intensity in at least 3 samples
kxub_heatmap_hep2 <- kxub_site[kxub_site_iden$Hep2_Total == 1, c('ACCID','Position','PROTEINS','Mod', meta$Intensity[meta$CellLine == 'Hep2'])]
df_temp <- kxub_heatmap_hep2[,5:13]
df_temp[df_temp > 0] <- 1
df_temp$total <- rowSums(df_temp)
kxub_heatmap_hep2 <- kxub_heatmap_hep2[df_temp$total > 2,]
kxub_heatmap_hep2[,5:13] <- log2(kxub_heatmap_hep2[,5:13] + 1)

ph <- pheatmap(kxub_heatmap_hep2[,5:13], clustering_method = 'ward.D2', scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
kxub_heatmap_hep2 <- kxub_heatmap_hep2[ph$tree_row$order,]
write.csv(kxub_heatmap_hep2, 'Dataset/Quantification/hep2_heatmap_order.csv', row.names = FALSE)

kxub_heatmap_jurkat <- kxub_site[kxub_site_iden$Jurkat_Total == 1, c('ACCID','Position','PROTEINS','Mod', meta$Intensity[meta$CellLine == 'Jurkat'])]
df_temp <- kxub_heatmap_jurkat[,5:13]
df_temp[df_temp > 0] <- 1
df_temp$total <- rowSums(df_temp)
kxub_heatmap_jurkat <- kxub_heatmap_jurkat[df_temp$total > 2,]
kxub_heatmap_jurkat[,5:13] <- log2(kxub_heatmap_jurkat[,5:13] + 1)

pj <- pheatmap(kxub_heatmap_jurkat[,5:13], clustering_method = 'ward.D2', scale = 'row',color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
kxub_heatmap_jurkat <- kxub_heatmap_jurkat[pj$tree_row$order,]
write.csv(kxub_heatmap_jurkat, 'Dataset/Quantification/jurkat_heatmap_order.csv', row.names = FALSE)

#MDS plot of samples
lfq <- kxub_site[,meta$Intensity]
colnames(lfq) <- meta$Sample
plotMDS(log2(lfq+1))


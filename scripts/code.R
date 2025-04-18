rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(devtools))
suppressMessages(library(customLayout))
suppressMessages(library(stringr))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(tidydr))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(corrplot))
suppressMessages(library(colorspace))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(maftools))
suppressMessages(library(vegan))
suppressMessages(library(forcats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(rstatix))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggcor))
suppressMessages(library(ggstance))
suppressMessages(library(tidyverse))
suppressMessages(library(GOplot))
suppressMessages(library(caret))
suppressMessages(library(writexl))
suppressMessages(library(rcartocolor))
suppressMessages(library(ggcorrplot))
suppressMessages(library(psych))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
suppressMessages(library(cols4all))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(scales))
suppressMessages(library(oncoPredict))
suppressMessages(library(gghalves))
suppressMessages(library(cowplot))
suppressMessages(library(IOBR))
suppressMessages(library(estimate))
suppressMessages(library(UpSetR))
suppressMessages(library(ggbiplot))
suppressMessages(library(ggsci))
suppressMessages(library(WGCNA))
suppressMessages(library(circlize))
suppressMessages(library(rJava))
suppressMessages(library(xlsxjars))
suppressMessages(library(xlsx))
suppressMessages(library(glmnet))
suppressMessages(library(tidyr))
suppressMessages(library(pROC))
suppressMessages(library(ROCR))
suppressMessages(library(patchwork))

#01Cluster###########################
tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp.txt')
tcga.n.exp=readMatrix('00_origin_datas/Preprocessed/tcga.n.exp.txt')
tcga.t.cli=readMatrix('00_origin_datas/Preprocessed/tcga.t.cli.txt')
sample_del=setdiff(colnames(tcga.t.exp), tcga.t.cli$SampleID)
tcga.t.exp=tcga.t.exp[,-which(colnames(tcga.t.exp)==sample_del)]
tcga.t.cli=tcga.t.cli[match(colnames(tcga.t.exp),tcga.t.cli$SampleID),]
identical(tcga.t.cli$SampleID,colnames(tcga.t.exp))
m6A <- read.table("00_origin_datas/Data Sheet 1.txt",header = T,check.names = F,fill=T,sep = "\t")
m6A_list <- list()
m6A_list[['m6A']] <- as.vector(unique(m6A$`Polyamines metabolism-associated gene`))
ssgsea_t_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga.t.exp,
                                                  genelist = m6A_list)
ssgsea_n_score <- ssGSEAScore_by_muti_group_genes(gene.exp = tcga.n.exp,
                                                  genelist = m6A_list)
ssgsea_t_score1 <- as.data.frame(t(ssgsea_t_score))
ssgsea_n_score1 <- as.data.frame(t(ssgsea_n_score))
ssgsea_t_score1$"Group"="Tumor"
ssgsea_n_score1$"Group"="Normol"
ssgsea_cli <- rbind(ssgsea_t_score1,ssgsea_n_score1)
ssgsea_cli$Group=factor(ssgsea_cli$Group,levels = c("Tumor","Normol"))

fig1a <- mg_violin(ssgsea_cli[, c("Group", "m6A")]
                   ,melt = T
                   ,xlab = ''
                   ,legend.pos = 'tl'
                   ,ylab = 'Polyamine metabolism')
fig1a
ggsave("01_ConsensusClusterPlus/fig1a.pdf", width = 6, height = 6)



#############
cox.pval=0.05
exp_m6A=as.matrix(tcga.t.exp[which(rownames(tcga.t.exp)%in%as.vector(m6A$`Polyamines metabolism-associated gene`)),])
tcga.t.cli$OS.time=as.numeric(tcga.t.cli$OS.time)
tcga.t.cli$OS=as.numeric(tcga.t.cli$OS)
identical(as.vector(tcga.t.cli$SampleID),colnames(exp_m6A))
tcga.pcd.cox=cox_batch(t(scale(t(as.matrix(exp_m6A))))
                       ,time = tcga.t.cli$OS.time/365
                       ,event = tcga.t.cli$OS)

table(tcga.pcd.cox$p.value<0.05)
table(tcga.pcd.cox$p.value<0.01)
table(tcga.pcd.cox$p.value<0.001)

tcga.pcd.cox=tcga.pcd.cox[order(tcga.pcd.cox$HR,decreasing = T),]
dev.off()

tcga.pcd.cox.sig=tcga.pcd.cox[which(tcga.pcd.cox$p.value<cox.pval),]
nrow(tcga.pcd.cox.sig)
#pdf('PDFs/bioForest.pdf',height = 5,width = 6,onefile = F)
#bioForest(rt = tcga.pcd.cox.sig,col=c('#5C8980','#A5604A'))
#dev.off()

write.table(tcga.pcd.cox.sig,'01_ConsensusClusterPlus/tcga.pcd.cox.sig.txt',sep = "\t",quote = F,row.names = T,col.names = T)

###########
clusterAlg_name=c('hc','pam','km','kmdist')[3]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[1]
tcga_consen_data=as.matrix(tcga.t.exp[which(rownames(tcga.t.exp)%in%(rownames(tcga.pcd.cox.sig))),])
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   

#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))
#tcga_consen_data=as.dist(1-cor(tcga_consen_data,method = 'spearman'))

#tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data, maxK = 10, reps = 500, pItem = 0.8, pFeature = 1, title = "01_ConsensusClusterPlus", clusterAlg = clusterAlg_name, distance = distance_name, plot = "pdf", writeTable = F, seed = 123456)#########
#save(tcga_clust_subtype,file='01_ConsensusClusterPlus/tcga.subtype.RData')
load('01_ConsensusClusterPlus/tcga.subtype.RData')#########
k=2
subtype.cols=c("#003f69","#97824b")
tcga.subtype <- data.frame( Samples=names(tcga_clust_subtype[[k]]$consensusClass),Subtype=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Subtype=paste0('C',tcga.subtype$Subtype)
write.table(tcga.subtype,file = '01_ConsensusClusterPlus/Subtype.txt',sep = '\t',quote = F,row.names = F,col.names = T)


table(tcga.subtype$Subtype)
colnames(tcga.t.cli)[1]='Samples'
tcga.subtype.cli=merge(tcga.subtype,tcga.t.cli,by='Samples')




fig1f=ggplotKMCox(data.frame(time = tcga.subtype.cli$OS.time/365
                             , event = tcga.subtype.cli$OS
                             , tcga.subtype.cli$Subtype) 
                  ,add_text = '',show_confint = F,palette = subtype.cols)


fig1f
ggsave("01_ConsensusClusterPlus/fig1f.pdf", width = 8, height = 8)

#################PCA
tcga.subtype$Samples=as.vector(tcga.subtype$Samples)
tcga_exp_var=t(tcga.t.exp[,tcga.subtype$Samples])

tcga_exp_var=tcga_exp_var[ , which(apply(tcga_exp_var, 2, var) != 0)]##
dim(tcga_exp_var)
cluster.pca <- prcomp(tcga_exp_var, scale=T)
cluster.pca.plot <- ggbiplot(cluster.pca, scale=1, groups = tcga.subtype$Subtype,
                             ellipse = TRUE,ellipse.prob=0.3, circle = F,var.axes=F) +
  scale_color_manual(values = subtype.cols) + 
  theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top',
        panel.grid = element_blank(),text = element_text(family = 'Times')) +
  xlab('PCA1') + ylab('PCA2')+xlim(-3,3)+ylim(-3,3)
cluster.pca.plot

####################
colnames(tcga.subtype.cli)
tcga.subtype.cli$Age=ifelse(tcga.subtype.cli$age>60,'>60','<=60')

tcga.subtype.cli.cmp=list()
for(i in c(14,8:9,12:13)){
  #group.color=subtype.cols[1:4]
  
  p=plotMutiBar_tmp(table(tcga.subtype.cli[,i],tcga.subtype.cli$Subtype)
                    #,fill.color = subtype.cols
                    ,isAuto = F,showValue = F
                    ,legTitle=colnames(tcga.subtype.cli)[i])
  tcga.subtype.cli.cmp=c(tcga.subtype.cli.cmp,list(p$Bar))
}
length(tcga.subtype.cli.cmp)

fig1e=mg_merge_plot(tcga.subtype.cli.cmp,nrow = 1,ncol = length(c(14,8:9,12:13))
                    ,labels='F')
fig1e
ggsave("01_ConsensusClusterPlus/fig1e2.pdf", width = 18, height = 6)




##############DEG##########
C1_sample=as.vector(tcga.subtype.cli$Samples[which(tcga.subtype.cli$Subtype=="C1")]) 
C2_sample=as.vector(tcga.subtype.cli$Samples[which(tcga.subtype.cli$Subtype=="C2")]) 

geo.limma_C1vsC2=mg_limma_DEG(exp = tcga.t.exp[,c(C1_sample,C2_sample)],group = c(rep("C1",length(C1_sample)),rep("C2",length(C2_sample))), ulab = 'C1',dlab = 'C2')
geo.limma_C1vsC2$Summary
df.deg.sig=geo.limma_C1vsC2$DEG[which(geo.limma_C1vsC2$DEG$adj.P.Val<0.05 & abs(geo.limma_C1vsC2$DEG$logFC)>log2(2)),]
write.table(df.deg.sig,'02_DEGs/limma_C1vsC2.txt',sep = "\t",quote = F,row.names = T,col.names = T)
dim(df.deg.sig)
############
fig3a=my_volcano_FDR(geo.limma_C1vsC2,p_cutoff = 0.05,fc_cutoff = log2(2),col = c('#BC3C29FF','#20854EFF','grey'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text = element_text(color = "black",family = 'Times',size = 14),
        axis.text = element_text(color = "black",family = 'Times',size = 14),
        legend.position = 'top')+xlim(-2,2)+ylim(0,20)
fig3a
ggsave("02_DEGs/fig3a.pdf", width = 6, height = 6)

#
########
enrichment=mg_clusterProfiler(as.vector(rownames(df.deg.sig)))
write.table(enrichment$Enrich_tab,file = '02_DEGs/enrichment.txt',sep = '\t',quote = F,row.names = T,col.names = T)
foldchange=df.deg.sig$logFC
names(foldchange) <- rownames(df.deg.sig)
enrichment_kegg=cnetplot(enrichment$KEGG, foldChange=foldchange,showCategory = 5,circular = TRUE, colorEdge = T)
enrichment_GO_BP=cnetplot(enrichment$GO_BP, foldChange=foldchange,showCategory = 5,circular = TRUE, colorEdge = T)
enrichment_GO_CC=cnetplot(enrichment$GO_CC, foldChange=foldchange,showCategory = 5,circular = TRUE, colorEdge = T)
enrichment_GO_MF=cnetplot(enrichment$GO_MF, foldChange=foldchange,showCategory = 5,circular = TRUE, colorEdge = T)

#enrichment_kegg=enrichplot::dotplot(enrichment$KEGG,showCategory=10,)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
#enrichment_GO_BP=enrichplot::dotplot(enrichment$GO_BP,showCategory=10)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
#enrichment_GO_CC=enrichplot::dotplot(enrichment$GO_CC,showCategory=10)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
#enrichment_GO_MF=enrichplot::dotplot(enrichment$GO_MF,showCategory=10)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))
#(enrichment_kegg + enrichment_GO_BP) / (enrichment_GO_CC + enrichment_GO_MF)

#####################
GEO_expression=tcga.t.exp[rownames(df.deg.sig),]
GEO_expression1 <- apply(GEO_expression, 1, scale)
rownames(GEO_expression1) <- colnames(tcga.t.exp)
GEO_expression2 <- as.data.frame(t(GEO_expression1))
gvhd.dist <- dist(GEO_expression1)
gvhd.hclust <- hclust(gvhd.dist, method = "ward.D2")
head(tcga.subtype.cli)
colnames(tcga.subtype.cli)
geo.sample1=tcga.subtype.cli[,c(2,8:9,12:14)]
rownames(geo.sample1)=tcga.subtype.cli$Samples
geo.sample1<-geo.sample1[gvhd.hclust$order,]
geo.sample2<-geo.sample1[str_order(geo.sample1$Subtype,decreasing = T),]

GEO_expression2 <-GEO_expression2[,rownames(geo.sample2)]
identical(rownames(geo.sample2), colnames(GEO_expression2))####
#
colnames(tcga.subtype.cli)[c(2,8:9,12:14)]
subtype.cols=c("#BC3C29","#0072B5",'#E18727','#20854E')
#subtype.cols=c("#008EA07F","#8A41987F", "#FF6F007F" ,"#FF63487F")

Subtype.color=subtype.cols[c(1:2)]
Age.color=subtype.cols[c(1,2)]
Gender.color=subtype.cols[c(1,2)]
Stage.color=subtype.cols[c(1:4)]
pathologic_T.color=subtype.cols[c(1:4)]
Grade.color=subtype.cols[c(1:4)]
#pathologic_N.color=subtype.cols[c(1,2)]
#pathologic_M.color=subtype.cols[c(1,2)]



names(Subtype.color)=c('C1','C2')
names(Age.color)=c('<=60','>60')
names(Gender.color)=c('Male','Female')
names(Stage.color)=c('Stage I','Stage II',"Stage III","Stage IV")
names(pathologic_T.color)=c('T1','T2',"T3","T4")
names(Grade.color)=c('G1','G2',"G3","G4")

#names(pathologic_N.color)=c('N0','N1')
#names(pathologic_M.color)=c('M0','M1')



head(geo.sample2)
column_ha=HeatmapAnnotation(df = geo.sample2
                            , na_col = "grey"
                            , annotation_height = unit(0.01, "mm")
                            , gap = unit(1, 'mm')
                            ,col = list(Subtype=Subtype.color
                                        ,gender=Gender.color
                                        ,Age=Age.color
                                        ,stage=Stage.color
                                        ,"T stage"=pathologic_T.color
                                        ,grade=Grade.color
                                        ))


identical(rownames(geo.sample2), colnames(GEO_expression2))####

heatmap_plot=Heatmap(as.matrix(GEO_expression2),
                     col = colorRamp2(c(-3, -1.5, 0, 1.5, 3),c("#226ED1", "#7AA8E2", "white", "#F09090","#E01010")), 
                     #col = circlize::colorRamp2(c(-3, 0, 3), c('navy', 'white', 'red')),
                     name = "Log2(TPM+1)",top_annotation = column_ha,
                     show_row_names = F,show_column_names = F,
                     clustering_method_rows = "complete",
                     row_names_gp = gpar(fontsize = 10),
                     column_names_gp = gpar(fontsize = 10),
                     column_names_rot = 45,
                     cluster_columns = F,
                     cluster_rows = T,
                     show_row_dend = F,
                     #  width = ncol(exp)*unit(cell_size, "mm"),
                     use_raster = F,
                     ##`use_raster` is automatically set to TRUE for a matrix with more than 2000 rows.
                     #  row_names_max_width = row_name_width,
                     heatmap_legend_param = list(direction = "vertical",
                                                 legend_width = unit(3.05,"cm"),
                                                 legend_height = unit(2.8, "cm"),
                                                 title_position = "lefttop-rot"))
pdf(file="02_DEGs/heatmap_plot.pdf",width=7,height=6)
heatmap_plot
dev.off()
ggsave('02_DEGs/enrichment_kegg.pdf',enrichment_kegg,height = 8,width = 9)
ggsave('02_DEGs/enrichment_GO_BP.pdf',enrichment_GO_BP,height = 8,width = 9)


#############lasso##############
###############


tcga.deg.sig=rownames(df.deg.sig)

length(tcga.deg.sig)

tcga.t.exp_use=as.data.frame(tcga.t.exp)
tcga.subtype.cli=tcga.t.cli
identical(colnames(tcga.t.exp_use),as.vector(tcga.subtype.cli$Samples))
colnames(tcga.subtype.cli)[1]="Samples"
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples


tcga.cox=cox_batch(t(scale(t(tcga.t.exp_use[tcga.deg.sig,])))
                   ,time =  tcga.subtype.cli$OS.time/365
                   ,event =tcga.subtype.cli$OS)
dim(tcga.cox)

table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)
writeMatrix(tcga.cox,outpath = '04_Lasso/tcga.cox.txt')



p.cutoff=0.05
tcga.cox_use=tcga.cox
tcga.cox_use$coef=log(tcga.cox_use$HR)
tcga.cox_use$Gene=rownames(tcga.cox_use)
tcga.cox_use$type=rep('None',nrow(tcga.cox_use))
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef>0)]='Risk'
tcga.cox_use$type[which(tcga.cox_use$p.value<p.cutoff & tcga.cox_use$coef<0)]='Protective'
table(tcga.cox_use$type)

######### lasso
tcga.gene.sig=rownames(tcga.cox)[which(tcga.cox$p.value<p.cutoff)]
length(tcga.gene.sig)

table(tcga.cox_use$type)


tcga.cox_forVis=tcga.cox_use
tcga.cox_forVis=tcga.cox_forVis[which(tcga.cox_forVis$type %in% c('Risk','Protective')),]
tcga.cox_forVis$p.value=-log10(tcga.cox_forVis$p.value)
range(tcga.cox_forVis$p.value)


#################### LASSO
table(tcga.cox$p.value<0.05)
table(tcga.cox$p.value<0.01)
table(tcga.cox$p.value<0.001)

tcga.exp.sig=tcga.t.exp_use[tcga.gene.sig,]
tcga.exp.sig=t(tcga.exp.sig)
dim(tcga.exp.sig)


dim(tcga.exp.sig)
options(ggrepel.max.hnscerlaps = Inf)
tcga.subtype.cli$Samples=as.vector(tcga.subtype.cli$Samples)
identical(rownames(tcga.exp.sig),tcga.subtype.cli$Samples)
tcga.lasso.res=mg_lasso_cox_use(tcga.exp.sig
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , nfolds = 10
                                , lambda.min = T
                                , figLabels=c('B','C'))
tcga.lasso.res$Genes

tcga.lasso.res$lambda

tcga.lasso.res$plot

tcga.exp.for.cox=tcga.t.exp_use[match(tcga.lasso.res$Genes,row.names(tcga.t.exp_use)),]
dim(tcga.exp.for.cox)
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)

lst.modl=createCoxModel_use((t(tcga.exp.for.cox))
                            , time = tcga.subtype.cli$OS.time/365
                            , event = tcga.subtype.cli$OS
                            , isStep =T)
lst.modl$Cox
lst.modl$Genes
lst.modl$fmla

lst.modl.Coef=lst.modl$Coef
names(lst.modl.Coef)=lst.modl$Genes
lst.modl.Coef

tcga.risk.score=lst.modl$Score
#tcga.risk.score=scale(tcga.risk.score)[,1]
tcga.risk.score=mosaic::zscore(tcga.risk.score)

range(tcga.risk.score)

lst.modl$Coef

gene.coef=data.frame(Gene=lst.modl$Genes,Coef=lst.modl$Coef)
gene.coef$Type=ifelse(lst.modl$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
table(gene.coef$Type)

fig3c=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  geom_text(aes(label=round(Coef,digits = 3)),color="black",hjust = "left")+
  
  coord_flip() +
  scale_fill_manual(values=pal_nejm(alpha = 0.9)(8)[c(1,2)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Lasso Cox coefficient") +
  theme_classic()+theme(legend.position = c(0,1))
# theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")


fig3AB=mg_plot_lasso_use(fit = tcga.lasso.res$Mode1
                         , cv_fit = tcga.lasso.res$Model2
                         , show_text = F
                         , figLabels = c('A', 'B'))
fig3AB
fig3abc=mg_merge_plot(fig3AB,fig3c,nrow = 1,ncol = 3,widths = c(1,2,1))
#savePDF('PDFs/Fig7AB.pdf',fig7A,height = 4,width = 9)
#savePDF('PDFs/fig3abc.pdf',fig3abc,height = 5,width = 15)


tcga.exp.forCox<- cbind(time=tcga.subtype.cli$OS.time/365,
                        status=tcga.subtype.cli$OS,
                        t(tcga.t.exp_use)[rownames(tcga.subtype.cli), lst.modl$Genes])



dim(tcga.exp.forCox)

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(lst.modl$Genes,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga.exp.forCox))
fig3d=survminer::ggforest(cox,data=tcga.exp.forCox,noDigits = 3)
fig3abd=mg_merge_plot(fig3AB,fig3d,nrow = 1,ncol =2,widths = c(2,1))

#savePDF('PDFs/fig3abc.pdf',fig3abcd,height = 5,width = 16)


############### TCGA
#tcga.cutoff <- survminer::surv_cutpoint(data.frame(time=tcga.subtype.cli$OS.time/365, event=tcga.subtype.cli$OS, risk=tcga.risk.score),time = "time", event = "event",variables = c("risk"))
#tcga.cutoff=tcga.cutoff$cutpoint$cutpoint
#tcga.cutoff=median(tcga.risk.score)
tcga.cutoff=0
identical(colnames(tcga.exp.for.cox),tcga.subtype.cli$Samples)
risk.group.color=c("#EE0000FF","#3B4992FF")
names(risk.group.color)=c('High','Low')
tcga.roc=plotCoxModel_Batch_use(riskScore = tcga.risk.score
                                ,dat = t(tcga.exp.for.cox[match(lst.modl$Genes, row.names(tcga.exp.for.cox)), ])
                                , time = tcga.subtype.cli$OS.time/365
                                , event = tcga.subtype.cli$OS
                                , cutoff = tcga.cutoff
                                , labs = c('High','Low')
                                , title = 'RiskType'
                                , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                                
                                , pal = risk.group.color
                                , mks = c(1,2,3))
tcga.roc1=tcga.roc[[1]]
#pdf('PDFs/fig3e2.pdf',height = 6,width = 6)
tcga.roc1

dev.off()

tcga.group=ifelse(tcga.risk.score>tcga.cutoff,'High','Low')
tcga.group=data.frame(tcga.group)
colnames(tcga.group)='group'
table(tcga.group$group)
write.table(cbind(tcga.risk.score,tcga.group),file = '04_Lasso/tcga.group.txt',sep='\t',quote = F)



###########################

icgc.t.exp=readMatrix('00_origin_datas/Preprocessed/icgc.t.exp.txt')
icgc.t.cli=readMatrix('00_origin_datas/Preprocessed/icgc.t.cli.txt')
identical(colnames(icgc.t.exp),icgc.t.cli$SampleID)
icgc.HCCDB18.t.exp=icgc.t.exp
match(lst.modl$Genes,row.names(icgc.HCCDB18.t.exp))
length(lst.modl$Genes)
icgc.HCCDB18.t.cli.os=icgc.t.cli
identical(icgc.HCCDB18.t.cli.os$SampleID , colnames(icgc.HCCDB18.t.exp)) ####
icgc.HCCDB18.model.dat=icgc.HCCDB18.t.exp[match(lst.modl$Genes,row.names(icgc.HCCDB18.t.exp)),]

icgc.HCCDB18.risk.score=predictScoreByCoxModel(coxModel = lst.modl
                                               ,(t(icgc.HCCDB18.model.dat)))
icgc.HCCDB18.risk.score=mosaic::zscore(icgc.HCCDB18.risk.score)

lst.modl$fmla
#icgc.HCCDB18.cutoff <- survminer::surv_cutpoint(data.frame(time=icgc.HCCDB18.t.cli.os$OS.time/365, event=icgc.HCCDB18.t.cli.os$OS, risk=icgc.HCCDB18.risk.score),time = "time", event = "event",variables = c("risk"))
#icgc.HCCDB18.cutoff=icgc.HCCDB18.cutoff$cutpoint$cutpoint

icgc.HCCDB18.cutoff=0

icgc.ROC=plotCoxModel_Batch_use(riskScore = icgc.HCCDB18.risk.score
                                , dat = t(icgc.HCCDB18.t.exp[intersect(lst.modl$Genes, row.names(icgc.HCCDB18.t.exp)),])
                                , time = as.numeric(icgc.HCCDB18.t.cli.os$OS.time/365) 
                                , event = as.numeric(icgc.HCCDB18.t.cli.os$OS)
                                , cutoff = icgc.HCCDB18.cutoff
                                , labs = c('High','Low')
                                , title = 'RiskType'
                                , hetColor = c('#3B4992FF', 'white', '#EE0000FF')
                                , pal = risk.group.color
                                , mks = c(1,2,3))
icgc.ROC1=icgc.ROC[[1]]

icgc.ROC1
dev.off()
icgc.HCCDB18.group=ifelse(icgc.HCCDB18.risk.score>icgc.HCCDB18.cutoff,'High','Low')
icgc.HCCDB18.group=data.frame(icgc.HCCDB18.group)
colnames(icgc.HCCDB18.group)='group'
table(icgc.HCCDB18.group)

Fig4_ROC=mg_merge_plot(tcga.roc1,icgc.ROC1,ncol=2,nrow=1)
Fig4=ggpubr::ggarrange(fig3abd,Fig4_ROC, ncol = 1, nrow = 2,heights = c(1,2))

ggsave('PDFs/Fig3.pdf',Fig4,height = 12,width = 15)
ggsave('PDFs/Fig3.jpg',Fig4,height = 12,width = 15)




identical(rownames(as.data.frame(tcga.risk.score)),rownames(tcga.subtype.cli))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=tcga.risk.score)
tcga.risktype.cli$Age=ifelse(tcga.risktype.cli$age>60,'>60','<=60')
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff,'High','Low')
write.table(tcga.risktype.cli,file = '04_Lasso/tcga.risktype.cli.txt',sep='\t',quote = F)
#########################
##########
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)[c(7:13,15)]

table(tcga_cox_datas$M.stage)

table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'


table(tcga_cox_datas$stage)
tcga_cox_datas$stage[tcga_cox_datas$stage=='Stage I'|tcga_cox_datas$stage=='Stage II']<-'Stage I+II'
tcga_cox_datas$stage[tcga_cox_datas$stage=='Stage III'|tcga_cox_datas$stage=='Stage IV']<-'Stage III+IV'


table(tcga_cox_datas$grade)
tcga_cox_datas$grade[tcga_cox_datas$grade=='G1'|tcga_cox_datas$grade=='G2']<-'G1+G2'
tcga_cox_datas$grade[tcga_cox_datas$grade=='G3'|tcga_cox_datas$grade=='G4']<-'G3+G4'

##
univar_res<-unicox(vars=colnames(tcga_cox_datas)[c(7:13,15)],time = tcga_cox_datas$OS.time,event = tcga_cox_datas$OS,data=tcga_cox_datas)
univar_res
univar_res[which(univar_res$pvalue<0.05),]


##
mutivar_res<-multicox(vars=rownames(univar_res[which(univar_res$pvalue<0.05),]),time = tcga_cox_datas$OS.time,event = tcga_cox_datas$OS,data=tcga_cox_datas,forest = F)
mutivar_res
mutivar_res[which(mutivar_res$pvalue<0.05),]


##
#dev.off()
mg_Forestplot(df_m = univar_res,outFile = 'PDFs/tcga.univar.forestplot.pdf',height = 4,width = 6)
mg_Forestplot(df_m = mutivar_res,outFile = 'PDFs/tcga.mutivar.forestplot.pdf',height = 4,width = 6)


###################nomo

dt=data.frame(RiskScore=tcga_cox_datas$Riskscore,
              T.Stage=tcga_cox_datas$T.stage,
              stage =tcga_cox_datas$stage, 
              M.Stage=tcga_cox_datas$M.stage
)

pdf('PDFs/nomogram.pdf', width = 12, height = 10)

nom.plot=mg_nomogram(clinical_riskscore=dt,
                     os = as.numeric(tcga_cox_datas$OS.time/365),
                     status = as.numeric(tcga_cox_datas$OS))


dev.off()

#pdf('PDFs/fig4b.pdf', width = 9, height = 8)
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1,3,5))
dev.off()



########################################
##############
tcga.risktype.cli=read.table('04_Lasso/tcga.risktype.cli.txt',header = T,check.names = F,fill=T,sep = "\t")
tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp.txt')
tcga.t.cli=readMatrix('00_origin_datas/Preprocessed/tcga.t.cli.txt')
sample_del=setdiff(colnames(tcga.t.exp), tcga.t.cli$SampleID)
tcga.t.exp=tcga.t.exp[,-which(colnames(tcga.t.exp)==sample_del)]
identical(as.vector(tcga.risktype.cli$Samples),colnames(tcga.t.exp))

risk.group.color1=c("#EE0000FF","#3B4992FF")
names(risk.group.color1)=c('High','Low')

library(estimate)
#### ESTIMATE
#tcga.exp.estimate<-deconvo_estimate(eset=tcga.t.exp)
#save(tcga.exp.estimate,file='06_imm/tcga.exp.estimate.RData')
load('06_imm/tcga.exp.estimate.RData')
tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)


############ TIMER 
#tcga.exp.timer<-deconvo_timer(eset=as.matrix(tcga.t.exp),indications=rep('LIHC',ncol(tcga.t.exp)))
#save(tcga.exp.timer,file='06_imm/tcga.exp.timer.RData')
load('06_imm/tcga.exp.timer.RData')
tcga.exp.timer=get.IOBR.immu.format(tcga.exp.timer)

library(EPIC)
############ EPIC 
#tcga.exp.epic<-deconvo_epic(eset=as.matrix(tcga.t.exp),tumor = TRUE)
#save(tcga.exp.epic,file='06_imm/tcga.exp.epic.RData')
load('06_imm/tcga.exp.epic.RData')
tcga.exp.epic=get.IOBR.immu.format(tcga.exp.epic)

############ MCP-counter 
#tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(tcga.t.exp))
#save(tcga.exp.mcp,file='06_imm/tcga.exp.mcp.RData')
load('06_imm/tcga.exp.mcp.RData')
tcga.exp.mcp=get.IOBR.immu.format(tcga.exp.mcp)

### CIBERSORT
#tcga.exp.cibersort<-deconvo_cibersort(eset=tcga.t.exp,arrays=F)
#save(tcga.exp.cibersort,file='06_imm/tcga.exp.cibersort.RData')
load('06_imm/tcga.exp.cibersort.RData')
tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)

#######sssGSEA#######
#geo.immu.ssgsea=immu_ssgsea(exp = tcga.t.exp)
#save(geo.immu.ssgsea,file='06_imm/geo.immu.ssgsea.RData')
#load('06_imm/geo.immu.ssgsea.RData')



#
tcga.t.estimate=tcga.exp.estimate[rownames(tcga.risktype.cli),1:3]
tcga.t.timer=tcga.exp.timer[rownames(tcga.risktype.cli),]
tcga.t.epic=tcga.exp.epic[rownames(tcga.risktype.cli),]
tcga.t.mcp=tcga.exp.mcp[rownames(tcga.risktype.cli),]
tcga.t.cibersort=tcga.exp.cibersort[rownames(tcga.risktype.cli),1:22]
#tcga.t.ssGSEA28=as.data.frame(geo.immu.ssgsea[tcga.risktype.cli$Samples,])




fig5a=get_PlotMutiBoxplot(tcga.t.estimate,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')

fig5a

fig5b=get_PlotMutiBoxplot(tcga.t.timer,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5b

fig5c=get_PlotMutiBoxplot(tcga.t.epic,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5c

fig5d=get_PlotMutiBoxplot(tcga.exp.mcp,tcga.risktype.cli
                          ,group_cols = risk.group.color1
                          #,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5d

groupViolin(tcga.t.estimate,
            tcga.risktype.cli$Risktype,
            ylab = 'ssgsea Immune Score',
            group_col=risk.group.color1)




fig5e=get_PlotMutiBoxplot(tcga.t.cibersort,tcga.risktype.cli
                          ,group_cols = risk.group.color
                          ,legend.pos = NULL
                          ,ylab = 'Score'
                          ,group.val = 'Risktype',xangle=45)+labs(color='Risktype')
fig5e


################imm############

#######
#library('oncoPredict')

#tcga.t.exp=readMatrix('00_origin_datas/Preprocessed/tcga.t.exp.txt')
#tcga.t.cli=readMatrix('00_origin_datas/Preprocessed/tcga.t.cli.txt')
#sample_del=setdiff(colnames(tcga.t.exp), tcga.t.cli$SampleID)
#tcga.t.exp=tcga.t.exp[,-which(colnames(tcga.t.exp)==sample_del)]
#tcga.t.cli=tcga.t.cli[match(colnames(tcga.t.exp),tcga.t.cli$SampleID),]
#identical(tcga.t.cli$SampleID,colnames(tcga.t.exp))
#drug_exp=as.matrix(tcga.t.exp)

#GDSC2_Expr = readRDS(file=file.path(dir,'Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
#GDSC2_Res = readRDS(file = file.path(dir,"Training Data/GDSC2_Res.rds"))
#GDSC2_Res <- exp(GDSC2_Res)
calcPhenotype(trainingExprData = as.matrix(GDSC2_Expr),
              trainingPtype = as.matrix(GDSC2_Res),
              testExprData = as.matrix(drug_exp),
             batchCorrect = 'eb',  #   "eb" for ComBat
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = TRUE,
             removeLowVaringGenesFrom = 'rawData' )

tcga_durg_ic50_res=read.csv('06_imm/calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
tcga.risktype.cli=read.table('04_Lasso/tcga.risktype.cli.txt',header = T,check.names = F,fill=T,sep = "\t")

tcga.risktype.cli=tcga.risktype.cli[rownames(tcga_durg_ic50_res),]
dim(tcga_durg_ic50_res)


IC50.mat=data.frame(Riskscore=tcga.risktype.cli$Riskscore,tcga_durg_ic50_res[as.vector(tcga.risktype.cli$Samples),])

IC50_RS_cor <- corr.test(x =IC50.mat$Riskscore,
                         y = IC50.mat[,-1])


IC50_RS_cor_res=data.frame(drugs=colnames( IC50.mat[,-1]))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.05,abs(IC50_RS_cor_res$cor)>0.4)
IC50_RS_cor_res=IC50_RS_cor_res[IC50_RS_cor_res$p.adj<0.05 & abs(IC50_RS_cor_res$cor)>0.4,]
IC50_RS_cor_res=IC50_RS_cor_res[order(IC50_RS_cor_res$cor),]
head(IC50_RS_cor_res)
IC50_plot <- ggplot(data=IC50_RS_cor_res,aes(x=cor,y=reorder(drugs,cor),
                                             color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_colour_gradient(low ='#ffc7c7' ,high = "#8785a2")+
  geom_segment(aes(yend=drugs,xend=0),size=1) +
  labs(x='spearman Correlation',y='')+theme_bw()+
  theme(text = element_text(family = 'Times'),legend.position = "bottom")

IC50_RS_cor_res$drugs
#[1] MK.1775_1179      BPD.00008900_1998 ML323_1629        Lapatinib_1558    WIKI4_1940       
#[6] GDC0810_1925      AZD6738_1917      BDP.00009066_1866 VE821_2111        AZD7762_1022     
#[11] PF.4708671_1129   Doramapimod_1042  SB505124_1194    
fig6a=mg_merge_plot(fig5b,fig5d,nrow = 1,ncol = 2,labels = LETTERS[1:2])
fig6b=mg_merge_plot(fig5e,IC50_plot,nrow = 1,ncol =2,labels = LETTERS[3:4],widths = c(2,1))
fig6=mg_merge_plot(fig6a,fig6b,nrow = 2,ncol =1)

ggsave("PDFs/Fig5.pdf", fig6, width = 14, height = 10)

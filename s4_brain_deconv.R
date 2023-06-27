library("DESeq2")
library("limma")
library("edgeR")
library(ggpubr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(immunedeconv)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggsci)
library(ggpubr)
options(stringsAsFactors = F)


###############################################
#################### data # ###################

# rosmap
ROSMAP_norm <- read.table('/Users/songlt/Documents/AD/cTI/rosmap/ROSMAP_arrayExpression_normalized.tsv',header = T)[,c(-1,-2)]
ros_meta <- read.csv('/Users/songlt/Documents/AD/cTI/rosmap/rosmap_clinical_2019-05_v3.csv')
rownames(ros_meta) <- paste(ros_meta$Study,ros_meta$projid, sep='_')

# cogdx : Final consensus cognitive diagnosis: Clinical consensus diagnosis of cognitive status at time of death
# 1: NC, 2-3 mci; 4-5: ad; 6 other dementia

ros_meta$Diagnosis <- car::recode(ros_meta$cogdx,"1='NC';2:3='MCI';4:5='AD';6='Other dementia'")
ros_meta <- subset(ros_meta, cogdx%in%c(1:5))

# ceradsc: CERAD score ; Semiquantitative measure of neuritic plaques(1-4; 4-1, 越小越严重)
# braaksc: Semiquantitative measure of neurofibrillary tangles (0-6; 0-VII)
# average replicates
ROSMAP_norm <- avereps(ROSMAP_norm[,-1],ID = ROSMAP_norm$Symbol)
ROSMAP_norm <- ROSMAP_norm[,intersect(colnames(ROSMAP_norm),rownames(ros_meta))]
ros_meta <- ros_meta[colnames(ROSMAP_norm),]
ros_meta$Diagnosis <- factor(ros_meta$Diagnosis, levels = c('NC','MCI','AD'))

# mayo
ensg_sym <- read.table('~/Desktop/AD2/data/ensg_sym_autosome.txt',header = T, stringsAsFactors = F)
ensg_sym <- ensg_sym[grepl('ENSG',rownames(ensg_sym)),]
rownames(ensg_sym) <- str_split_fixed(rownames(ensg_sym), '\\.',2)[,1]

TCX_exp <- read.table('/Users/songlt/Documents/AD/cTI/mayo/MayoRNAseq_RNAseq_TCX_geneCounts_normalized.tsv',
                      header = T, row.names = 1, check.names = F)
CBE_exp <- read.table('/Users/songlt/Documents/AD/cTI/mayo/MayoRNAseq_RNAseq_CBE_geneCounts_normalized.tsv',
                      header = T, row.names = 1, check.names = F)

TCX_exp <- TCX_exp[intersect(rownames(TCX_exp), rownames(ensg_sym)),]
CBE_exp <- CBE_exp[intersect(rownames(CBE_exp), rownames(ensg_sym)),]

TCX_exp$gene_symbol <- ensg_sym[rownames(TCX_exp),'gene_symbol']
TCX_exp <- avereps(TCX_exp[,-ncol(TCX_exp)],ID = TCX_exp$gene_symbol)

CBE_exp$gene_symbol <- ensg_sym[rownames(CBE_exp),'gene_symbol']
CBE_exp <- avereps(CBE_exp[,-ncol(CBE_exp)],ID = CBE_exp$gene_symbol)

mayo_meta <- read.csv('/Users/songlt/Documents/AD/cTI/mayo/MayoRNAseq_individual_metadata.csv')
mayo_meta$Diagnosis <- car::recode(mayo_meta$diagnosis,"'control'='NC';'Alzheimer Disease'='AD';else='others'")

mayo_meta_tcx <- mayo_meta
rownames(mayo_meta_tcx) <- paste(mayo_meta_tcx$individualID, 'TCX',sep='_')

mayo_meta_cbe <- mayo_meta
rownames(mayo_meta_cbe) <- paste(mayo_meta_cbe$individualID, 'CER',sep='_')

mayo_meta_tcx <- mayo_meta_tcx[colnames(TCX_exp),]
mayo_meta_cbe <- mayo_meta_cbe[colnames(CBE_exp),]

# MSBB
#MSBB_cli <- read.table('../Metadata/AMP-AD_MSBB_MSSM_covariates_mRNA_AffymetrixU133AB.tsv',header = T,row.names = 1)
MSBB_cli <- read.csv('/Users/songlt/Documents/AD/cTI/MSBB/Metadata/MSBB_individual_metadata.csv',header = T,row.names = 1)
MSBB_cli$Diagnosis <- car::recode(MSBB_cli$CDR,"0='NC';0.5='MCI';else='AD'")

###############################################
#################### function # ###############

color_code <- c( "#979797","#c4c4c4", "#af1319")

names(color_code) <- c('NC','MCI','AD')

Deconv_immune <- function(GEM, meta_data, array , deconvol_method,dataset){
  # GEM: expression matrix 
  # array T/F
  # meta_data: meta_data
  
  meta_data <- subset(meta_data,Diagnosis%in%c('NC','SCD','MCI','AD') )
  GEM <- GEM[,rownames(meta_data)]
  
  res = deconvolute(GEM, method=deconvol_method, tumor = F,arrays=array)
  
  res <- as.data.frame(res)## 检验不同诊断人群免疫细胞组成是否有差别
  rownames(res) <- res$cell_type
  res <- as.data.frame(t(res[,-1]))
  
  if (deconvol_method=='quantiseq'){res$Macrophage <- res$`Macrophage M1`+res$`Macrophage M2`}
  
  res2 <- melt(as.matrix(res))
  res2$diagnosis  <- meta_data[as.character(res2$Var1),'Diagnosis']
  colnames(res2) <- c('sample','cell_type','value','diagnosis')
  res2$diagnosis <- factor(res2$diagnosis, levels=c('NC','SCD','MCI','AD'))
  res2$dataset <- dataset
  res2$method <- deconvol_method
  
  return(list(res_m=res,res_p=res2))
  
}

plot_box <-function(cell_perc,nr){
  
  #cell_perc <- subset(cell_perc, diagnosis!='SCD' & !cell_type%in%c("uncharacterized cell","Macrophage M1",'Macrophage M2'))
  cell_perc <- subset(cell_perc, diagnosis!='SCD' & !cell_type%in%c("Macrophage M1",'Macrophage M2'))
  
  #p <- ggplot(subset(cell_perc,  cell_type%in%unique(cell_perc$cell_type[cell_perc$value>0])), aes(x=diagnosis, y= value )) +
    p <- ggplot(subset(cell_perc,  cell_type%in%unique(cell_perc$cell_type[cell_perc$value>0])), aes(x=diagnosis, y= value*100 )) +
      geom_violin(width=1,size=0.3)+
    geom_boxplot(width=0.2,aes(fill=diagnosis),outlier.shape ='',size=0.3)+
    scale_fill_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    scale_color_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    #scale_fill_jco()+scale_color_jco()+
    facet_wrap(~cell_type,nrow = nr,scale='free_y')+ylab('Percentage (%)')+theme_bw()+xlab('')+
    stat_compare_means(#comparisons = list(c('AD','NC'),c('AD','MCI'),c('NC','MCI')), method = "wilcox.test",label = "p.signif",
                       #method.args = list(alternative = "greater"),
      comparisons = list(c('AD','NC')), method = "wilcox.test",label = "p.signif",
                       tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                       map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
    theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =8))
  
  return(p)
}

plot_box_dataset <-function(cell_perc,nr,celltp){
  
  cell_perc <- subset(cell_perc,  cell_type==celltp)
  #ggplot(subset(cell_perc,  cell_type%in%unique(cell_perc$cell_type[cell_perc$value>0])), aes(x=diagnosis, y= value )) +
  p <- ggplot(subset(cell_perc), aes(x=diagnosis, y= value )) +
    
    geom_violin(alpha=0.3,width=0.5,size=0.3)+
    geom_boxplot(width=0.2,aes(fill=diagnosis),outlier.shape ='',size=0.3)+#scale_fill_jco()+scale_color_jco()+
    scale_fill_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    scale_color_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    facet_wrap(~dataset,nrow = nr,scale='free')+ylab('MCP-counter score')+theme_bw()+xlab('')+
    stat_compare_means(comparisons = list(c('AD','NC')),method = "wilcox.test",
                       #method.args = list(alternative = "less"),
                       tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                       map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
    theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =9))
  
  return(p)
}

#ggscater
ggscater_plot <- function(cell_p,measure,dataset){
  
  ggscater_p <- ggscatter(subset(cell_p, cell_type=='Neutrophil' & value >0 & !is.na(measure)), 
                          x='value',y=measure, size = 0.5, color="diagnosis" ,# palette = 'jco',# Points color, shape and size
                          palette=color_code,cor.coef.size = 3,
                          add = "reg.line",  # Add regressin line
                          add.params = list(color = "#3D59A5", fill = "lightgray"), # Customize reg. line
                          conf.int = TRUE, # Add confidence interval
                          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                          cor.coeff.args = list(method = "spearman", label.sep = "\n"))+
    theme_bw()+facet_grid(.~dataset)+
    theme(legend.position = '', text = element_text(size = 9))+xlab('Neutrophil Score')#+ggtitle(paste(dataset,measure,sep=':'))
  
  return(ggscater_p)
}


###############################################
#################### main # ###################

# rosmap 比例计算

###############################################
#################### quantiseq # ##############

ROSMAP_Deconv <- Deconv_immune(GEM=2^ROSMAP_norm, meta_data=ros_meta, array =T, deconvol_method='quantiseq',dataset='ROSMAP')
Mayo_cbe_Deconv <- Deconv_immune(GEM=CBE_exp, meta_data=mayo_meta_cbe, array =F, deconvol_method='quantiseq',dataset='Mayo_cbe')
Mayo_tcx_Deconv <- Deconv_immune(GEM=TCX_exp, meta_data=mayo_meta_tcx, array =F, deconvol_method='quantiseq',dataset='Mayo_tcx')

###############################################
#################### epic # ###################
###############################################

ROSMAP_Deconv_epic <- Deconv_immune(GEM=2^ROSMAP_norm, meta_data=ros_meta, array =T, deconvol_method='epic',dataset='ROSMAP')
Mayo_cbe_Deconv_epic <- Deconv_immune(GEM=CBE_exp, meta_data=mayo_meta_cbe, array =F, deconvol_method='epic',dataset='Mayo_cbe')
Mayo_tcx_Deconv_epic <- Deconv_immune(GEM=TCX_exp, meta_data=mayo_meta_tcx, array =F, deconvol_method='epic',dataset='Mayo_tcx')

plot_box(cell_perc=ROSMAP_Deconv_epic$res_p, nr=1)
Mayo_cbe_epic <- plot_box(cell_perc=Mayo_cbe_Deconv_epic$res_p, nr=1)
Mayo_tcx_epic <- plot_box(cell_perc=Mayo_tcx_Deconv_epic$res_p, nr=1)
plot_grid(Mayo_cbe_epic,Mayo_tcx_epic, nrow=2,labels = 'AUTO',label_size = 9)
dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/Fig4s_mayo_brain_epic.pdf', width = 7.3,height = 4)

###############################################
#################### mcp_counter # ############
###############################################

ROSMAP_Deconv_mcpcounter <- Deconv_immune(GEM=ROSMAP_norm, meta_data=ros_meta, array =T, deconvol_method='mcp_counter',dataset='ROSMAP')
Mayo_cbe_Deconv_mcpcounter <- Deconv_immune(GEM=CBE_exp, meta_data=mayo_meta_cbe, array =F, deconvol_method='mcp_counter',dataset='Mayo_cbe')
Mayo_tcx_Deconv_mcpcounter_NC <- Deconv_immune(GEM=TCX_exp, meta_data=subset(mayo_meta_tcx,Diagnosis=='NC'), array =F, deconvol_method='mcp_counter',dataset='Mayo_tcx')
Mayo_tcx_Deconv_mcpcounter_AD <- Deconv_immune(GEM=TCX_exp, meta_data=subset(mayo_meta_tcx,Diagnosis=='AD'), array =F, deconvol_method='mcp_counter',dataset='Mayo_tcx')

plot_box(cell_perc=ROSMAP_Deconv_mcpcounter$res_p, nr=1)
plot_box(cell_perc=Mayo_cbe_Deconv_mcpcounter$res_p, nr=2)
plot_box(cell_perc=rbind(Mayo_tcx_Deconv_mcpcounter_NC$res_p,Mayo_tcx_Deconv_mcpcounter_AD$res_p), nr=1)

###############################################
#################### Cibersort # ##############
###############################################


source('/Users/songlt/Desktop/AD2/code2/deconv_code/Cibersort.R')
lm6 <- read.table('~/Desktop/AD2/data/signature_rnaseq_geo60424_LM6.txt',sep='\t',header = T, row.names = 1)

ROSMAP_ciber <- CIBERSORT(sig_matrix=lm6,mixture_file =2^ROSMAP_norm , QN = T)
CBE_ciber <- CIBERSORT(sig_matrix=lm6,mixture_file =CBE_exp , QN = F)
TCX_ciber <- CIBERSORT(sig_matrix=lm6,mixture_file =TCX_exp , QN = F)

cibersort_deconv <- function(res_cibersort,dataset,LM, meta_data){
  res_cibersort <- as.data.frame(res_cibersort)
  #res_cibersort <- read.table(paste0('~/Desktop/AD2/data/CIBERSORT.',dataset,'_', LM, '.txt' ), sep='\t', row.names = 1, header = T, check.names = F )
  res_cibersort <- res_cibersort[,1:(ncol(res_cibersort)-3)]
  
  res_cibersort <- res_cibersort[,colSums(res_cibersort)!=0]
  res_cibersort$sample <- rownames(res_cibersort)
  res_cibersort <- melt(res_cibersort)
  res_cibersort$diagnosis <- meta_data[res_cibersort$sample,'Diagnosis']
  
  colnames(res_cibersort) <- c('sample','cell_type','value','diagnosis')
  res_cibersort <- subset(res_cibersort,diagnosis%in%c('NC','SCD','MCI','AD') )
  res_cibersort$diagnosis <- factor(res_cibersort$diagnosis, levels=c('NC','SCD','MCI','AD'))
  
  res_cibersort$dataset <- dataset
  res_cibersort$method <- 'cibersort'
  
  return(res_cibersort)
}

ROSMAP_cibersort_LM6 <- cibersort_deconv(ROSMAP_ciber,dataset='ROSMAP', LM='LM6', meta_data=ros_meta)
CBE_cibersort_LM6 <- cibersort_deconv(CBE_ciber,dataset='CBE', LM='LM6', meta_data=mayo_meta_cbe)
TCX_cibersort_LM6 <- cibersort_deconv(TCX_ciber,dataset='TCX', LM='LM6', meta_data=mayo_meta_tcx)

ROSMAP_cibersort_LM6_p <- plot_box(cell_perc=ROSMAP_cibersort_LM6, nr=1)
CBE_cibersort_LM6_p <- plot_box(cell_perc=CBE_cibersort_LM6, nr=1)
TCX_cibersort_LM6_p <- plot_box(cell_perc=TCX_cibersort_LM6, nr=1)


#ROSMAP cibersort web page
cibersort <- function(dataset,LM, meta_data){
  
  res_cibersort <- read.table(('~/Downloads/CIBERSORT.Output_Job17.txt' ), sep='\t', row.names = 1, header = T, check.names = F )
  res_cibersort <- res_cibersort[,1:(ncol(res_cibersort)-3)]
  
  res_cibersort <- res_cibersort[,colSums(res_cibersort)!=0]
  res_cibersort$sample <- rownames(res_cibersort)
  res_cibersort <- melt(res_cibersort)
  res_cibersort$diagnosis <- meta_data[res_cibersort$sample,'Diagnosis']
  
  colnames(res_cibersort) <- c('sample','cell_type','value','diagnosis')
  res_cibersort <- subset(res_cibersort,diagnosis%in%c('NC','SCD','MCI','AD') )
  res_cibersort$diagnosis <- factor(res_cibersort$diagnosis, levels=c('NC','SCD','MCI','AD'))
  
  res_cibersort$dataset <- dataset
  res_cibersort$method <- 'cibersort'
  
  return(res_cibersort)
}
ros_cibersort_LM6 <- cibersort(dataset='rosmap', LM='LM6', meta_data=ros_meta)

#save(MSBB_cell_perc, ROSMAP_Deconv, Mayo_cbe_Deconv , Mayo_tcx_Deconv, file='~/Desktop/AD2/data/brain_deconv2.RData')

load('~/Desktop/AD2/data/brain_deconv2.RData')

# #################### 临床指标 # ###################

# rosmap 
ROSMAP_p <- ROSMAP_Deconv$res_p
ROSMAP_p$value <- ROSMAP_p$value*100
ROSMAP_p$CERAD <- ros_meta[(ROSMAP_p$sample),'ceradsc']
ROSMAP_p$Braak <- ros_meta[(ROSMAP_p$sample),'braaksc']
ROSMAP_p$MMSE <- ros_meta[(ROSMAP_p$sample),'cts_mmse30_lv']

# Mayo_cbe

Mayo_cbe_p <- Mayo_cbe_Deconv_mcpcounter$res_p
#Mayo_cbe_p$value <- Mayo_cbe_p$value*100
Mayo_cbe_p$Braak <- mayo_meta_cbe[as.character(Mayo_cbe_p$sample),'Braak']
Mayo_cbe_p$thal <- mayo_meta_cbe[as.character(Mayo_cbe_p$sample),'thal']
#Mayo_cbe_p$thal[Mayo_cbe_p$diagnosis=='NC' & is.na(Mayo_cbe_p$thal)  ] <- 0

# Mayo_tcx
Mayo_tcx_p <-rbind(Mayo_tcx_Deconv_mcpcounter_NC$res_p,Mayo_tcx_Deconv_mcpcounter_AD$res_p)
#Mayo_tcx_p$value <- Mayo_tcx_p$value*100
Mayo_tcx_p$Braak <- mayo_meta_tcx[as.character(Mayo_tcx_p$sample),'Braak']
Mayo_tcx_p$thal <- mayo_meta_tcx[as.character(Mayo_tcx_p$sample),'thal']
#Mayo_tcx_p$thal[Mayo_tcx_p$diagnosis=='NC' & is.na(Mayo_tcx_p$thal)  ] <- 0


## plot
## neutrophils abundance


Mayo_mcp <- plot_box_dataset(cell_perc=rbind(Mayo_cbe_Deconv_mcpcounter$res_p,Mayo_tcx_Deconv_mcpcounter_NC$res_p,Mayo_tcx_Deconv_mcpcounter_AD$res_p),
                                     nr=1,celltp='Neutrophil')

#dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/Mayo_mcp.pdf', width=3, height=6)

  
# Correlation
library(cowplot)

Mayo_tcx_Braak <- ggscater_plot(cell_p=Mayo_tcx_p, measure='Braak',dataset='Mayo_tcx')
Mayo_tcx_thal <- ggscater_plot(cell_p=Mayo_tcx_p, measure='thal',dataset='Mayo_tcx')

Mayo_cbe_Braak <- ggscater_plot(cell_p=Mayo_cbe_p, measure='Braak',dataset='Mayo_cbe')
Mayo_cbe_thal <- ggscater_plot(cell_p=Mayo_cbe_p, measure='thal',dataset='Mayo_cbe')

plot_grid( Mayo_cbe_Braak,Mayo_cbe_thal,Mayo_tcx_Braak,Mayo_tcx_thal,nrow = 2)
#dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/Mayo_score_scater.pdf', width=6, height=6)

plot_grid( Mayo_cbe_Braak,Mayo_cbe_thal,nrow = 1,rel_widths = c(1,1), labels = 'AUTO',label_size = 9)
dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/Fig4s_mayo_brain_cbe.pdf', width = 6,height = 3)

plot_grid( Mayo_mcp,Mayo_tcx_Braak,Mayo_tcx_thal,nrow = 1,rel_widths = c(2,1,1), labels = 'AUTO',label_size = 9)
dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/Fig4_mayo_brain.pdf', width = 7.3,height = 2.5)


ROSMAP_CERAD <- ggscater_plot(cell_p=ROSMAP_p, measure='CERAD',dataset='ROSMAP')
ROSMAP_Braak <- ggscater_plot(cell_p=ROSMAP_p, measure='Braak',dataset='ROSMAP')
ROSMAP_MMSE <- ggscater_plot(cell_p=ROSMAP_p, measure='MMSE',dataset='ROSMAP')


plot_grid(ROS_Neutrophil, ROSMAP_CERAD,ROSMAP_Braak, ROSMAP_MMSE,nrow = 1)

dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/ROSMAP_scater.pdf', width=8, height=3)


#MSBB
setwd('/Users/songlt/Documents/AD/cTI/MSBB/Processed/')
MSBB_cell_perc <- c()
for (region in dir('./')) {

  MSBB <- read.table(region, sep='\t',header = T) [,c(-1,-2,-4)]
  MSBB <- avereps(MSBB[,-1],ID = MSBB$Gene.Symbol)

  ROSMAP_Deconv <- Deconv_immune(GEM=2^MSBB, meta_data=MSBB_cli[colnames(MSBB),], array =T, deconvol_method='quantiseq',dataset=region)
  MSBB_cell_perc <- rbind(ROSMAP_Deconv$res_p,MSBB_cell_perc)
}
cell_perc <- MSBB_cell_perc
cell_perc <- subset(cell_type=='B cell')
cell_perc$region <- str_split_fixed(cell_perc$dataset,'_',5)[,4]
ggplot(cell_perc, aes(x=diagnosis, y= value*100 )) +
  geom_violin(alpha=0.3,width=0.5,size=0.3)+
  geom_boxplot(alpha=0.6,width=0.2,aes(fill=diagnosis),outlier.shape ='',size=0.3)+scale_fill_jco()+scale_color_jco()+
  facet_wrap(~region,nrow = 2,scale='free_y')+ylab('Percentage (%)')+theme_bw()+
  stat_compare_means(comparisons = list(c('AD','NC')),label = "p.signif", method = "wilcox.test",
                     #method.args = list(alternative = "greater"),
                     tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                     map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
  theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =8))

dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/MSBB.pdf', width=11, height=5)



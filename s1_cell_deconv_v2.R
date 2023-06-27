library(immunedeconv)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(stringr)

load('~/Desktop/AD2/data/exp_meta.RData')

coldata$Diagnosis[coldata$Diagnosis=='SCD'] <- 'NC'
label_DXSUM_ec$Diagnosis[label_DXSUM_ec$Diagnosis=='SCD'] <- 'NC'

## immune cell deconvolution

Deconv_immune <- function(GEM, meta_data, array , deconvol_method,dataset){
  # GEM: expression matrix 
  # array T/F
  # meta_data: meta_data
  
  meta_data <- subset(meta_data,Diagnosis%in%c('NC','SCD','MCI','AD') )
  GEM <- GEM[,rownames(meta_data)]
  
  res = deconvolute(GEM, method=deconvol_method, tumor = F,arrays=array)
  
  res <- as.data.frame(res)## 
  rownames(res) <- res$cell_type
  res <- as.data.frame(t(res[,-1]))
  
  if (deconvol_method=='quantiseq'){res$Macrophage <- res$`Macrophage M1`+res$`Macrophage M2`}
  
  res2 <- melt(as.matrix(res))
  colnames(res2) <- c('sample','cell_type','value')
  res2$sample <- as.character(res2$sample)
  res2$diagnosis  <- meta_data[res2$sample,'Diagnosis']
  res2$diagnosis <- factor(res2$diagnosis, levels=c('NC','SCD','MCI','AD'))
  res2$Gender <- meta_data[res2$sample,'Gender']
  res2$Age <- meta_data[res2$sample,'Age']
  res2$dataset <- dataset
  res2$method <- deconvol_method

  
  res$NLR <- res[,grepl('Neutro',colnames(res))]/rowSums(res[,grepl('T|B',colnames(res))])
  res$diagnosis <- meta_data[rownames(res),'Diagnosis']
  res$dataset <- dataset
  res$method <- deconvol_method
  res$Gender <- meta_data[rownames(res),'Gender']
  res$Age <- meta_data[rownames(res),'Age']
  
  return(list(res_m=res,res_p=res2))
  
}

plot_box <-function(cell_perc,nr){
  #cell_perc$diagnosis[cell_perc$diagnosis=='SCD'] <- 'NC'
  cell_perc <- subset(cell_perc, !diagnosis%in%c('SCD') & !cell_type%in%c("uncharacterized cell","Macrophage M1",'Macrophage M2','NLR'))
  p <-  ggplot(subset(cell_perc,  cell_type%in%unique(cell_perc$cell_type[cell_perc$value>0])), 
               aes(x=diagnosis, y= value*100 )) +
    geom_violin(alpha=0.3,width=0.5,size=0.3)+
    geom_boxplot(width=0.2,aes(fill=diagnosis),outlier.shape ='',size=0.3)+
    scale_fill_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    scale_color_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    #scale_fill_jco()+scale_color_jco()+
    facet_wrap(~cell_type,nrow = nr)+ylab('Percentage (%)')+theme_bw()+
    stat_compare_means(comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI')),label = "p.value", method = "wilcox.test",
                       tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                       map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
    theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =8))
  
  # geom_signif(test='wilcox.test', comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI'))
  #             ,map_signif_level=F, size = 0.2, textsize =2,step_increase = 0.05 ,tip_length = 0.01)
  
  
  return(p)
}

plot_box_nlr <-function(cell_perc){
  #cell_perc$diagnosis[cell_perc$diagnosis=='SCD'] <- 'NC'
  cell_perc <- subset(cell_perc, !diagnosis%in%c('SCD') )
  p <-  ggplot(cell_perc, 
               aes(x=diagnosis, y= NLR)) +
    geom_violin(alpha=0.3,width=0.5,size=0.3)+
    geom_boxplot(width=0.2,aes(fill=diagnosis),outlier.shape ='',size=0.3)+
    scale_fill_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    scale_color_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    #scale_fill_jco()+scale_color_jco()+
    ylab('NLR')+theme_bw()+
    stat_compare_means(comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI')),label = "p.value", method = "wilcox.test",
                       tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                       map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
    theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =8))
  
  # geom_signif(test='wilcox.test', comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI'))
  #             ,map_signif_level=F, size = 0.2, textsize =2,step_increase = 0.05 ,tip_length = 0.01)
  
  
  return(p)
}

plot_box_mcpcounter <-function(cell_perc,nr){
  #cell_perc$diagnosis[cell_perc$diagnosis=='SCD'] <- 'NC'
  cell_perc <- subset(cell_perc, !diagnosis%in%c('SCD') & !cell_type%in%c("uncharacterized cell","Macrophage M1",'Macrophage M2'))
  p <-  ggplot(subset(cell_perc,  cell_type%in%unique(cell_perc$cell_type[cell_perc$value>0])), 
               aes(x=diagnosis, y= value )) +
    geom_violin(alpha=0.3,width=0.5,size=0.3)+
    geom_boxplot(width=0.2,aes(fill=diagnosis),outlier.shape ='',size=0.3)+
    scale_fill_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    scale_color_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    #scale_fill_jco()+scale_color_jco()+
    facet_wrap(~cell_type,nrow = nr)+ylab('Score')+theme_bw()+
    stat_compare_means(comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI')),label = "p.signif", method = "wilcox.test",
                       tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                       map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
    theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =8))
  
  # geom_signif(test='wilcox.test', comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI'))
  #             ,map_signif_level=F, size = 0.2, textsize =2,step_increase = 0.05 ,tip_length = 0.01)
  
  
  return(p)
}

plot_box_wrap <-function(cell_perc,nr){
  cell_perc <- subset(cell_perc, !diagnosis%in%c('SCD') & !cell_type%in%c("uncharacterized cell","Macrophage M1",'Macrophage M2'))
  
  percentage_plot <- ggplot(subset(cell_perc,  cell_type%in%unique(cell_perc$cell_type[cell_perc$value>0])), 
                            aes(x=diagnosis, y= value*100 )) +
    geom_boxplot(alpha=1,width=0.2,aes(fill=diagnosis),outlier.shape = '',size=0.3,position = position_dodge(0.5))+
    geom_violin(aes(fill=diagnosis),alpha=0,position = position_dodge(0.5),size=0.3 )+
    scale_fill_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    scale_color_manual(values=color_code[ as.character(sort(unique(cell_perc$diagnosis)))])+
    theme_bw()+
    facet_wrap(.~cell_type,scale='free_y',nrow = 1)+
    theme( legend.position = '') +ylab('Percentage (%)')+xlab('')+theme_bw()+
    stat_compare_means(comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI')),label = "p.signif", method = "wilcox.test",
                       tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,0.1,1), symbols = c("****", "***", "**", "*",'.',"ns")),
                       size=2.8,family='sans',face='plain')+
    theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =8))
  
  return(percentage_plot)
  
}

#color_code <- c("#c4c4c4", "#979797", "#af1319")
color_code <- c( "#979797","#c4c4c4", "#af1319")

names(color_code) <- c('NC','MCI','AD')

ADNI_exp <- 2^adni
ANM1_exp <- 2^ANM1_norm
ANM2_exp <- 2^ANM2_norm

######################################################     epic    ###############################################################################

ZIB_epic <- Deconv_immune(GEM=ZIB, dataset='ZIB', meta_data=coldata, array=F, deconvol_method='epic')
ADNI_epic <- Deconv_immune(GEM=ADNI_exp, dataset='ADNI',meta_data=label_DXSUM_ec, array=T, deconvol_method='epic')
ANM1_epic <- Deconv_immune(GEM=ANM1_exp, dataset='ANM1',meta_data=anm1_meta, array=T, deconvol_method='epic')
ANM2_epic <- Deconv_immune(GEM=ANM2_exp, dataset='ANM2',meta_data=anm2_meta, array=T, deconvol_method='epic')

ZIB_epic_p <- plot_box(cell_perc=ZIB_epic$res_p, nr=1)
ADNI_epic_p <- plot_box(cell_perc=ADNI_epic$res_p, nr=1)
ANM1_epic_p <- plot_box(cell_perc=ANM1_epic$res_p, nr=1)
ANM2_epic_p <- plot_box(cell_perc=ANM2_epic$res_p, nr=1)

plot_box_nlr(ZIB_epic$res_m)
plot_box_nlr(ADNI_epic$res_m)
plot_box_nlr(ANM1_epic$res_m)
plot_box_nlr(ANM2_epic$res_m)

# ggsave(ZIB_epic_p, file='~/Desktop/AD2/figure/CELL_Deconv/ZIB_epic_p.pdf', width = 8,height = 4)
# ggsave(ADNI_epic_p, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_epic_p.pdf', width = 8,height = 3)
# ggsave(ANM1_epic_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM1_epic_p.pdf', width = 8,height = 4)
# ggsave(ANM2_epic_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM2_epic_p.pdf', width = 8,height = 4)


######################################################     quantiseq   ###############################################################################
ZIB_quantiseq <- Deconv_immune(GEM=ZIB, dataset='ZIB', meta_data=coldata, array=F, deconvol_method='quantiseq')
ADNI_quantiseq <- Deconv_immune(GEM=ADNI_exp, dataset='ADNI', meta_data=label_DXSUM_ec, array=T, deconvol_method='quantiseq')
ANM1_quantiseq <- Deconv_immune(GEM=ANM1_exp, dataset='ANM1', meta_data=anm1_meta, array=T, deconvol_method='quantiseq')
ANM2_quantiseq <- Deconv_immune(GEM=ANM2_exp, dataset='ANM2', meta_data=anm2_meta, array=T, deconvol_method='quantiseq')

ZIB_quantiseq_p <- plot_box(cell_perc=(ZIB_quantiseq$res_p), nr=1)
ADNI_quantiseq_p <- plot_box(cell_perc=ADNI_quantiseq$res_p, nr=1)
ANM1_quantiseq_p <- plot_box(cell_perc=ANM1_quantiseq$res_p, nr=1)
ANM2_quantiseq_p <- plot_box(cell_perc=ANM2_quantiseq$res_p, nr=1)

plot_box_nlr(ZIB_quantiseq$res_m)
plot_box_nlr(ADNI_quantiseq$res_m)
plot_box_nlr(ANM1_quantiseq$res_m)
plot_box_nlr(ANM2_quantiseq$res_m)

# ggsave(ZIB_quantiseq_p, file='~/Desktop/AD2/figure/CELL_Deconv/ZIB_quantiseq_p.pdf', width = 8,height = 4)
# ggsave(ADNI_quantiseq_p, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_quantiseq_p.pdf', width = 8,height = 4)
# ggsave(ANM1_quantiseq_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM1_quantiseq_p.pdf', width = 8,height = 4)
# ggsave(ANM2_quantiseq_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM2_quantiseq_p.pdf', width = 8,height = 4)
 
######################################################     mcp_counter   ###############################################################################
ZIB_mcp_counter <- Deconv_immune(GEM=ZIB, dataset='ZIB', meta_data=coldata, array=F, deconvol_method='mcp_counter')
ADNI_mcp_counter <- Deconv_immune(GEM=ADNI_exp, dataset='ADNI', meta_data=label_DXSUM_ec, array=T, deconvol_method='mcp_counter')
ANM1_mcp_counter <- Deconv_immune(GEM=ANM1_exp, dataset='ANM1', meta_data=anm1_meta, array=T, deconvol_method='mcp_counter')
ANM2_mcp_counter <- Deconv_immune(GEM=ANM2_exp, dataset='ANM2', meta_data=anm2_meta, array=T, deconvol_method='mcp_counter')

ZIB_mcp_counter_p <- plot_box_mcpcounter(cell_perc=ZIB_mcp_counter$res_p, nr=1)
ADNI_mcp_counter_p <- plot_box_mcpcounter(cell_perc=ADNI_mcp_counter$res_p, nr=1)
ANM1_mcp_counter_p <- plot_box_mcpcounter(cell_perc=ANM1_mcp_counter$res_p, nr=1)
ANM2_mcp_counter_p <- plot_box_mcpcounter(cell_perc=ANM2_mcp_counter$res_p, nr=1)

# ggsave(ZIB_mcp_counter_p, file='~/Desktop/AD2/figure/CELL_Deconv/ZIB_mcp_counter_p.pdf', width = 8,height = 4)
# ggsave(ADNI_mcp_counter_p, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_mcp_counter_p.pdf', width = 8,height = 4)
# ggsave(ANM1_mcp_counter_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM1_mcp_counter_p.pdf', width = 8,height = 4)
# ggsave(ANM2_mcp_counter_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM2_mcp_counter_p.pdf', width = 8,height = 4)


######################################################     cibersort deconv       ###############################################################################

source('/Users/songlt/Desktop/AD2/code2/deconv_code/Cibersort.R')
lm6 <- read.table('~/Desktop/AD2/data/signature_rnaseq_geo60424_LM6.txt',sep='\t',header = T, row.names = 1)
lm22 <- read.table('~/Downloads/LM22.txt',sep='\t',header = T, row.names = 1)

adni_ciber_LM6 <- CIBERSORT(sig_matrix=lm6,mixture_file =ADNI_exp , QN = T)
ZIB_ciber_LM6 <- CIBERSORT(sig_matrix=lm6,mixture_file =ZIB , QN = F)
ANM1_ciber_LM6 <- CIBERSORT(sig_matrix=lm6,mixture_file =ANM1_exp , QN = T)
ANM2_ciber_LM6 <- CIBERSORT(sig_matrix=lm6,mixture_file =ANM2_exp , QN = T)

adni_ciber_LM22 <- CIBERSORT(sig_matrix=lm22,mixture_file =ADNI_exp , QN = T)
ZIB_ciber_LM22 <- CIBERSORT(sig_matrix=lm22,mixture_file =ZIB , QN = F)
ANM1_ciber_LM22 <- CIBERSORT(sig_matrix=lm22,mixture_file =ANM1_exp , QN = T)
ANM2_ciber_LM22 <- CIBERSORT(sig_matrix=lm22,mixture_file =ANM2_exp , QN = T)

cibersort_deconv <- function(res_cibersort,dataset,LM, meta_data){
  res_cibersort <- as.data.frame(res_cibersort)
  #res_cibersort <- read.table(paste0('~/Desktop/AD2/data/CIBERSORT.',dataset,'_', LM, '.txt' ), sep='\t', row.names = 1, header = T, check.names = F )
  res_cibersort <- res_cibersort[,1:(ncol(res_cibersort)-3)]
  
  res_cibersort <- res_cibersort[,colSums(res_cibersort)!=0]
  res_cibersort$NLR <- res_cibersort[,grepl('Neutro',colnames(res_cibersort))]/rowSums(res_cibersort[,grepl('T|B',colnames(res_cibersort))])

  res_cibersort$sample <- rownames(res_cibersort)
  res_cibersort <- melt(res_cibersort)
  res_cibersort$diagnosis <- meta_data[res_cibersort$sample,'Diagnosis']
  res_cibersort$Gender <- meta_data[res_cibersort$sample,'Gender']
  res_cibersort$Age <- meta_data[res_cibersort$sample,'Age']
  
  
  colnames(res_cibersort) <- c('sample','cell_type','value','diagnosis','Gender','Age')
  res_cibersort <- subset(res_cibersort,diagnosis%in%c('NC','SCD','MCI','AD') )
  res_cibersort$diagnosis <- factor(res_cibersort$diagnosis, levels=c('NC','SCD','MCI','AD'))
  
  res_cibersort$dataset <- dataset
  res_cibersort$method <- 'cibersort'
  
  
  return(res_cibersort)
}

ZIB_cibersort_LM6 <- cibersort_deconv(ZIB_ciber_LM6, dataset='ZIB', LM='LM6', meta_data=coldata)
ADNI_cibersort_LM6 <- cibersort_deconv(adni_ciber_LM6, dataset='ADNI', LM='LM6', meta_data=label_DXSUM_ec)
ANM1_cibersort_LM6 <- cibersort_deconv(ANM1_ciber_LM6, dataset='ANM1', LM='LM6', meta_data=anm1_meta)
ANM2_cibersort_LM6 <- cibersort_deconv(ANM2_ciber_LM6, dataset='ANM2', LM='LM6', meta_data=anm2_meta)

ZIB_cibersort_LM22 <- cibersort_deconv(ZIB_ciber_LM22, dataset='ZIB', LM='LM22', meta_data=coldata)
ADNI_cibersort_LM22 <- cibersort_deconv(adni_ciber_LM22, dataset='ADNI', LM='LM22', meta_data=label_DXSUM_ec)
ANM1_cibersort_LM22 <- cibersort_deconv(ANM1_ciber_LM22, dataset='ANM1', LM='LM22', meta_data=anm1_meta)
ANM2_cibersort_LM22 <- cibersort_deconv(ANM2_ciber_LM22, dataset='ANM2', LM='LM22', meta_data=anm2_meta)

ZIB_cibersort_LM6_p <- plot_box(cell_perc=ZIB_cibersort_LM6, nr=1)
ADNI_cibersort_LM6_p <- plot_box(cell_perc=ADNI_cibersort_LM6, nr=1)
ANM1_cibersort_LM6_p <- plot_box(cell_perc=ANM1_cibersort_LM6, nr=1)
ANM2_cibersort_LM6_p <- plot_box(cell_perc=ANM2_cibersort_LM6, nr=1)

ZIB_cibersort_LM22_p <- plot_box(cell_perc=ZIB_cibersort_LM22, nr=1)
ADNI_cibersort_LM22_p <- plot_box(cell_perc=ADNI_cibersort_LM22, nr=1)
ANM1_cibersort_LM22_p <- plot_box(cell_perc=ANM1_cibersort_LM22, nr=1)
ANM2_cibersort_LM22_p <- plot_box(cell_perc=ANM2_cibersort_LM22, nr=1)


# ggsave(ZIB_cibersort_LM6_p, file='~/Desktop/AD2/figure/CELL_Deconv/ZIB_cibersort_LM6_p.pdf', width = 8,height = 4)
# ggsave(ADNI_cibersort_LM6_p, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_cibersort_LM6_p.pdf', width = 8,height = 4)
# ggsave(ANM1_cibersort_LM6_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM1_cibersort_LM6_p.pdf', width = 8,height = 4)
# ggsave(ANM2_cibersort_LM6_p, file='~/Desktop/AD2/figure/CELL_Deconv/ANM2_cibersort_LM6_p.pdf', width = 8,height = 4)
# 

save(ZIB_epic, ADNI_epic, ANM1_epic, ANM2_epic,
     ZIB_quantiseq,ADNI_quantiseq,ANM1_quantiseq,ANM2_quantiseq,
     ZIB_cibersort_LM6,ADNI_cibersort_LM6,ANM1_cibersort_LM6,ANM2_cibersort_LM6,
     ZIB_mcp_counter,ADNI_mcp_counter, ANM1_mcp_counter, ANM2_mcp_counter,
     file='~/Desktop/AD2/data/blood_cell_deconv.RData')

load( '~/Desktop/AD2/data/blood_cell_deconv.RData' )




## PCA 
library(FactoMineR)
library(factoextra)
library(cowplot)

make_pca_cell <- function(res_m, meta_data){
  
  #cell.pca <- PCA(res_m[rownames(subset(meta_data,Diagnosis%in%c('NC','AD'))),]%>%mutate(NL=Neutrophil/`B cell`), graph = FALSE)
  cell.pca <- PCA(res_m[rownames(subset(meta_data,Diagnosis%in%c('NC','AD'))),], graph = FALSE)
  
  pca_cell <- fviz_pca_biplot(cell.pca, 
                              col.ind = subset(meta_data,Diagnosis%in%c('NC','AD'))[,'Diagnosis'] ,
                              palette =as.character(color_code[ c('AD','NC')]),labelsize = 3,
                              addEllipses = TRUE, label='var',ellipse.level=0.8,pointsize=0.3,
                              repel = T,title=c("")) +
    xlab('PC1 (30.0%)')+ylab('PC2 (21.3%)')+
    theme_bw()+
    theme(legend.position = c(0.1,0.8),legend.margin=margin(c(0,0,0,0)),
          panel.grid.major = element_blank(),legend.key.size = unit(0.1,'in'),
          panel.grid.minor = element_blank(),
          legend.title  = element_blank(), text = element_text(size = 8))
  
  
  cell_contrib <- fviz_contrib(cell.pca, choice = "var", axes = 1, top = 10,,title=c(""))+
    theme_bw()+xlab('')+ylab('Contribution to PC1')+
    theme(legend.position = c(0.8,0.8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1 ,vjust = 1),
          legend.title  = element_blank(), text = element_text(size = 8))
  
  
  pca_p <- plot_grid(pca_cell,cell_contrib,nrow = 2  )
  
  #pca_p <- plot_grid(pca_cell)
  return(pca_p)
}

ADNI_epic_pca <- make_pca_cell(res_m=ADNI_epic$res_m[,c(1:6)], meta_data=label_DXSUM_ec)

ZIB_quantiseq_pca <- make_pca_cell(res_m=ZIB_quantiseq$res_m[,c(-2,-3,-11)], meta_data=coldata)
ADNI_quantiseq_pca <- make_pca_cell(res_m=ADNI_quantiseq$res_m[,c(-2,-3,-11)], meta_data=label_DXSUM_ec)
ANM1_quantiseq_pca <- make_pca_cell(res_m=ANM1_quantiseq$res_m[,c(-2,-3,-11)], meta_data=anm1_meta)
ANM2_quantiseq_pca <- make_pca_cell(res_m=ANM2_quantiseq$res_m[,c(-2,-3,-11)], meta_data=anm2_meta)

ggsave(ADNI_epic_pca, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_epic_pca.pdf', width = 4,height = 4)
ggsave(ADNI_quantiseq_pca, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_quantiseq_pca.pdf', width = 4,height = 4)

ggsave(ZIB_quantiseq_pca, file='~/Desktop/AD2/figure/CELL_Deconv/ZIB_quantiseq_pca.pdf', width = 8,height = 4)
ggsave(ADNI_quantiseq_pca, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_quantiseq_pca.pdf', width = 8,height = 4)
ggsave(ANM1_quantiseq_pca, file='~/Desktop/AD2/figure/CELL_Deconv/ANM1_quantiseq_pca.pdf', width = 8,height = 4)
ggsave(ANM2_quantiseq_pca, file='~/Desktop/AD2/figure/CELL_Deconv/ANM2_quantiseq_pca.pdf', width = 8,height = 4)


############################################################################################################################################
# Correlations between immune cell proportion, cognitive measurements, brain‐area volumes, and demographic variables 


# Age and gender
adni_meta <- label_DXSUM_ec
#adni_meta$DX[adni_meta$DX=='SCD'] <- 'NC'
demo_sex <- adni_meta %>%
  group_by(DX,PTGENDER) %>%
  dplyr::summarise(N=n()) %>%
  group_by(DX) %>%
  dplyr::mutate(freq = (N/sum(N)*100))

chisq.test(matrix(unlist(demo_sex[c(1:4),3]),nrow=2))$p.value
chisq.test(matrix(unlist(demo_sex[-3:-4,3]),nrow=2))$p.value

adni_meta %>% group_by(DX) %>% dplyr::summarise(mean(Age),na.rm=T)
adni_meta %>% group_by(DX) %>% dplyr::summarise(sd(age))

t.test(adni_meta$age[adni_meta$DX=='NC'],adni_meta$age[adni_meta$DX=='AD'])$p.value


#ADNI_moca
ADNI_moca <- read.csv('/Users/songlt/Documents/AD/ADNI/Neuropsychological/MOCA.csv')
ADNI_moca$SERIAL_sum <-  rowSums(ADNI_moca[,grepl('SERIAL',colnames(ADNI_moca))])
ADNI_moca$SERIAL_score <- car::recode(ADNI_moca$SERIAL_sum, "0=0;1=1;2:3=2;4:5=3")
ADNI_moca$FFLUENCY <- ifelse(ADNI_moca$FFLUENCY >=11,1,0)
ADNI_moca$DELW_score <- rowSums(ADNI_moca[,grepl('DELW',colnames(ADNI_moca))]==1)
ADNI_moca$LETTERS <- ifelse(ADNI_moca$LETTERS <=1,1,0)
ADNI_moca$MOCA <- rowSums(ADNI_moca[,c(-1:-8,-17:-26,-30:-34,-40:-44,-51:-53)])
ADNI_moca$id_v <-  paste(ADNI_moca$RID,ADNI_moca$VISCODE,sep='--')

# ADNI mmse
ADNI_mmse <- read.csv('/Users/songlt/Documents/AD/ADNI/Neuropsychological/MMSE.csv')
ADNI_mmse$id_v <- paste(ADNI_mmse$RID,ADNI_mmse$VISCODE,sep='--')

# ADNI TAU
# ADNI_tau <- read.csv('/Users/songlt/Documents/AD/ADNI/Neuropsychological/BLENNOWPLASMATAU.csv')
# ADNI_tau$id_v <- paste(ADNI_tau$RID,ADNI_tau$VISCODE,sep='--')
label_DXSUM_ec$sample_vis <- rownames(label_DXSUM_ec)
label_DXSUM_ec1 <- merge(ADNI_mmse[,c('id_v','MMSCORE')], label_DXSUM_ec, by='id_v', all.y=T)
label_DXSUM_ec2 <- merge(ADNI_moca[,c('id_v','MOCA')], label_DXSUM_ec1, by='id_v', all.y=T)
#label_DXSUM_ec2 <- merge(ADNI_tau[,c('id_v','PLASMATAU')], label_DXSUM_ec2, by='id_v', all.y=T)

rownames(label_DXSUM_ec2) <- label_DXSUM_ec2$sample_vis

# logistic regression to adjust age and sex
ADNI_P <- ADNI_epic$res_m
#ADNI_P$diagnosis[ADNI_P$diagnosis=='SCD'] <- 'NC'
#ADNI_P <- ADNI_quantiseq$res_m
#ADNI_P$age <- label_DXSUM_ec2[rownames(ADNI_P),'age']
ADNI_P$Sex <- ifelse(ADNI_P$Gender=='Female',0,1) 
ADNI_P$MoCA_B <- label_DXSUM_ec2[rownames(ADNI_P),'MOCA']
ADNI_P$MMSE <- label_DXSUM_ec2[rownames(ADNI_P),'MMSCORE']

ADNI_P <- subset(ADNI_P, diagnosis%in%c('NC','AD'))
ADNI_P$Diagnosis <- ifelse(ADNI_P$diagnosis=='AD',1,0)

fit <- glm( Diagnosis~Neutrophil+Age+Sex,family = gaussian,data=ADNI_P)
summary(fit)

fit <- glm( Diagnosis~`B cell`+Age+Sex,family = gaussian,data=ADNI_P)
summary(fit)

fit <- glm( Diagnosis~NLR+Age+Sex,family = gaussian,data=ADNI_P)
summary(fit)

# cognitive function
library(corrplot)
library(Hmisc)

#quantiseq 
#cor_data <- ADNI_P[,!grepl('diagnosis|unchara|CD4|NL|M1|M2',colnames(ADNI_P)) ]

#epic 
cor_data <- ADNI_P[,!grepl('diagnosis|unchara|M1|M2|dataset|method|Gender',colnames(ADNI_P)) ]

corp <- rcorr(as.matrix(cor_data),type='spearman')

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

ADNI_corrplot <- corrplot(cor(cor_data,use = "pairwise.complete.obs",method ='spearman' ),        # Correlation matrix
                          method = "square", # Correlation plot method
                          type = "full",    # Correlation plot style (also "upper" and "lower")
                          diag = TRUE,      # If TRUE (default), adds the diagonal
                          bg = "white",     # Background color
                          title = "" ,
                          order='original',
                          tl.srt=45, tl.cex=0.8, tl.col = "black",
                          p.mat= corp$P,
                          insig = "blank",sig.level = 0.1,
                          col = col2(20), cl.length = 21,cl.align = "c",cl.ratio = 0.2,cl.cex = 0.6)

#dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_epic_corrplot.pdf', width = 8,height = 8)
dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/fig1c_epic_corrplot.pdf', width = 7.3,height = 7.3)

####### figure 1 ####### ####### ####### 
ADNI_epic_box <- plot_box_wrap(cell_perc=ADNI_epic$res_p, nr=1)

# nlr
nlr_cols <- c('NLR', 'diagnosis','dataset', 'method', 'Gender', 'Age')
NLR_epic_quati <- rbind(ZIB_epic$res_m[,nlr_cols], ADNI_epic$res_m[,nlr_cols], ANM1_epic$res_m[,nlr_cols], ANM2_epic$res_m[,nlr_cols],
                  ZIB_quantiseq$res_m[,nlr_cols],ADNI_quantiseq$res_m[,nlr_cols],ANM1_quantiseq$res_m[,nlr_cols],ANM2_quantiseq$res_m[,nlr_cols],
                  ZIB_mcp_counter$res_m[,nlr_cols],ADNI_mcp_counter$res_m[,nlr_cols],ANM1_mcp_counter$res_m[,nlr_cols],ANM2_mcp_counter$res_m[,nlr_cols])

NLR_ciber <- rbind(ZIB_cibersort_LM6, ADNI_cibersort_LM6, ANM1_cibersort_LM6, ANM2_cibersort_LM6)%>%subset(cell_type=='NLR')%>%rename(NLR=value)

NLR_multi <- rbind( NLR_epic_quati, NLR_ciber[,nlr_cols] )
NLR_multi$diagnosis <- factor(NLR_multi$diagnosis, levels = c('NC','MCI','AD'))
NLR_multi$method <- factor(NLR_multi$method, levels = c("epic", "quantiseq", "cibersort",'mcp_counter'))
NLR_box <- ggplot(subset(NLR_multi, diagnosis%in%c('NC','AD') & NLR<15 ), aes(x=dataset, y= NLR,fill=diagnosis )) +
  geom_violin(alpha=0,width=0.8,size=0.2,position = position_dodge(0.8))+
  geom_boxplot(width=0.15,outlier.shape ='',size=0.2,position = position_dodge(0.8))+
  facet_wrap(.~method,scale='free_y',nrow = 1)+
  #scale_fill_jco()+scale_color_jco() +
  scale_fill_manual(values=color_code[ c('NC','AD')])+
  scale_color_manual(values=color_code[ c('NC','AD')])+
  ylab('NLR')+xlab('')+theme_bw() +#ylim(c(0,20))+
  stat_compare_means(label = "p.signif", method = "wilcox.test",vjust=0.5,
                     method.args = list(alternative = "greater"),
                     tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,0.1,1), symbols = c("****", "***", "**", "*",'.',"ns"))
                     ,size=2.8,family='sans')+
  theme(panel.spacing = unit(0.2, "lines"),text = element_text(size =9),legend.title = element_blank(), legend.key.size = unit(0.2,'in'))

#plot_grid(ADNI_epic_box, NLR_box,nrow = 2, labels = c('A','B'),label_size = 9)

#dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/Fig1ab.pdf', width = 7.3,height = 4)

plot_grid(plot_grid(ADNI_epic_box, NLR_box,nrow = 2, labels = c('A','B'),label_size = 9),
          plot_grid(ADNI_epic_pca,'',labels = c('C','D'),nrow = 1,label_size = 9,rel_widths =c(1,2) ), nrow = 2)

dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/Fig1abc2.pdf', width = 7.3,height = 8)

NLR_box_mci <- ggplot(subset(NLR_multi, diagnosis%in%c('NC','MCI') & NLR<15 ), aes(x=dataset, y= NLR,fill=diagnosis )) +
  geom_violin(alpha=0,width=0.8,size=0.2,position = position_dodge(0.8))+
  geom_boxplot(width=0.15,outlier.shape ='',size=0.2,position = position_dodge(0.8))+
  facet_wrap(.~method,scale='free_y',nrow = 1)+
  #scale_fill_jco()+scale_color_jco() +
  scale_fill_manual(values=color_code[ c('NC','MCI')])+
  scale_color_manual(values=color_code[ c('NC','MCI')])+
  ylab('NLR')+xlab('')+theme_bw() +#ylim(c(0,20))+
  stat_compare_means(label = "p.signif", method = "wilcox.test",vjust=0.5,
                     tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,0.1,1), symbols = c("****", "***", "**", "*",'.',"ns"))
                     ,size=2.8,family='sans')+
  theme(panel.spacing = unit(0.2, "lines"),text = element_text(size =9),legend.title = element_blank(), legend.key.size = unit(0.2,'in'))

ggsave(NLR_box_mci, file='~/Desktop/AD2/figure/CELL_Deconv/NLR_box_mci.pdf', width = 7.3,height = 2.5)

# figure3


### 不同方法得到的中性粒细胞的比例
cell_p <- rbind(ZIB_epic$res_p,ADNI_epic$res_p, ANM1_epic$res_p,ANM2_epic$res_p,
                ZIB_quantiseq$res_p,ADNI_quantiseq$res_p, ANM1_quantiseq$res_p,ANM2_quantiseq$res_p,
                ZIB_cibersort_LM6,ADNI_cibersort_LM6, ANM1_cibersort_LM6,ANM2_cibersort_LM6)#,
cell_p$method <- factor(cell_p$method, levels = c('epic','quantiseq','cibersort'))

Neutrophil_p <- ggplot(subset(cell_p, diagnosis%in%c('NC','AD') & grepl('Neutrophil',cell_type)), aes(x=dataset, y= value*100,fill=diagnosis )) +
  geom_violin(alpha=0,width=0.8,size=0.3,position = position_dodge(0.8))+
  geom_boxplot(width=0.2,outlier.shape ='',size=0.3,position = position_dodge(0.8))+
  facet_wrap(.~method)+
  #scale_fill_jco()+scale_color_jco() +
  scale_fill_manual(values=color_code[ c('NC','AD')])+
  scale_color_manual(values=color_code[ c('NC','AD')])+
  ylab('Proportion (%)')+theme_bw() +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     method.args = list(alternative = "greater"),vjust = 1,
                     tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,0.1,1), symbols = c("****", "***", "**", "*",'.',"ns"))
                     ,size=2.5,family='sans',face='plain')+
  theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =9))


ggsave(Neutrophil_p, file='~/Desktop/AD2/figure/CELL_Deconv/Neutrophil_datasets.pdf', width = 8,height = 4)

# Bcell
Bcell_p <- ggplot(subset(cell_p, diagnosis%in%c('NC','AD') & grepl('B',cell_type)), aes(x=dataset, y= value*100,fill=diagnosis )) +
  geom_violin(alpha=0,width=0.8,size=0.3,position = position_dodge(0.8))+
  geom_boxplot(width=0.2,outlier.shape ='',size=0.3,position = position_dodge(0.8))+
  facet_wrap(.~method)+
  #scale_fill_jco()+scale_color_jco() +
  scale_fill_manual(values=color_code[ c('NC','AD')])+
  scale_color_manual(values=color_code[ c('NC','AD')])+
  ylab('Proportion (%)')+theme_bw() +
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     method.args = list(alternative = "less"),vjust = 1,
                     tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,0.1,1), symbols = c("****", "***", "**", "*",'.',"ns"))
                     ,size=2.5,family='sans',face='plain')+
  theme(legend.key.size = unit(0.15,'in'),legend.title = element_blank(),panel.spacing = unit(0.1, "lines"),text = element_text(size =9))+
  ylim(c(0,25))

ggsave(Bcell_p, file='~/Desktop/AD2/figure/CELL_Deconv/Bcell_datasets.pdf', width = 8,height = 4)

plot_grid(plot_grid(Neutrophil_p,'',rel_widths =  c(1,0.1),labels = 'A',label_size = 9),
          Bcell_p, nrow=2,labels = c('','B'),label_size = 9)

dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_NB.pdf', width = 7.3,height =4.3 )

# mcp counter
mcp_mult <-  rbind(ZIB_mcp_counter$res_p,ADNI_mcp_counter$res_p, ANM1_mcp_counter$res_p,ANM2_mcp_counter$res_p)%>%
  subset(cell_type%in%c('B cell','Neutrophil'))


mcp_mult$diagnosis <- factor(mcp_mult$diagnosis, levels = c('NC','MCI','AD'))
mcp_box <- ggplot(subset(mcp_mult, diagnosis%in%c('NC','AD')), aes(x=dataset, y= value,fill=diagnosis )) +
  geom_violin(alpha=0,width=0.8,size=0.2,position = position_dodge(0.8))+
  geom_boxplot(width=0.15,outlier.shape ='',size=0.2,position = position_dodge(0.8))+
  facet_wrap(.~cell_type,scale='free_y',nrow = 1)+
  #scale_fill_jco()+scale_color_jco() +
  scale_fill_manual(values=color_code[ c('NC','AD')])+
  scale_color_manual(values=color_code[ c('NC','AD')])+
  ylab('MCP-counter score')+xlab('')+theme_bw() +#ylim(c(0,20))+
  stat_compare_means(label = "p.signif", method = "wilcox.test",vjust=0.5,
                     tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,0.1,1), symbols = c("****", "***", "**", "*",'.',"ns"))
                     ,size=2.8,family='sans')+
  theme(panel.spacing = unit(0.2, "lines"),text = element_text(size =9),legend.title = element_blank(), legend.key.size = unit(0.2,'in'))

ggsave(mcp_box, file='~/Desktop/AD2/figure/CELL_Deconv/mco_counter_nb.pdf', width = 7.3,height =2.5 )

# correlation between cell proportion and cognition

ADNI_ggscatter_plot <- function(ADNI_P,X_LAB,y_LAB){
  
  ADNI_P$MoCA_B <- label_DXSUM_ec2[rownames(ADNI_P),'MOCA']
  ADNI_P$MMSE <- label_DXSUM_ec2[rownames(ADNI_P),'MMSCORE']
  ADNI_P$NL <- ADNI_P$Neutrophil/ADNI_P$`B cell`
  ADNI_P$diagnosis <- label_DXSUM_ec2[rownames(ADNI_P),'Diagnosis']
  
  P1 <- ggscatter(subset(ADNI_P,diagnosis%in%c('NC','AD','MCI') ), 
                  x=X_LAB,y=y_LAB, size = 0.2, color="diagnosis" , palette = as.character(color_code)[c(3,2,1)],# Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.sep = "\n"))+theme(legend.title = element_blank(),legend.position = '', text = element_text(size=8))
    
  return(P1)
}
ADNI_ggscatter_plot <- function(ADNI_P,X_LAB,y_LAB){
  
  ADNI_P$MoCA_B <- label_DXSUM_ec2[rownames(ADNI_P),'MOCA']
  ADNI_P$MMSE <- label_DXSUM_ec2[rownames(ADNI_P),'MMSCORE']
  ADNI_P$diagnosis <- label_DXSUM_ec2[rownames(ADNI_P),'Diagnosis']
  
  P1 <- ggscatter(subset(ADNI_P,diagnosis%in%c('NC','AD') ), color='grey',
                  x=X_LAB,y=y_LAB, size = 0.2,# Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.sep = "\n"))+theme(legend.title = element_blank(),legend.position = '', text = element_text(size=8))
  
  return(P1)
}

## ADNI_epic
P1 <- ADNI_ggscatter_plot(ADNI_P=ADNI_epic$res_m,X_LAB='Neutrophil',y_LAB='MMSE') 
P2 <- ADNI_ggscatter_plot(ADNI_P=ADNI_epic$res_m,X_LAB='B cell',y_LAB='MMSE') 
P3 <- ADNI_ggscatter_plot(ADNI_P=ADNI_epic$res_m,X_LAB='NLR',y_LAB='MMSE') 

P4 <- ADNI_ggscatter_plot(ADNI_P=ADNI_epic$res_m,X_LAB='Neutrophil',y_LAB='MoCA_B') 
P5 <- ADNI_ggscatter_plot(ADNI_P=ADNI_epic$res_m,X_LAB='B cell',y_LAB='MoCA_B') 
P6 <- ADNI_ggscatter_plot(ADNI_P=ADNI_epic$res_m,X_LAB='NLR',y_LAB='MoCA_B') 
ADNI_epic_corr <- plot_grid(P1,P2,P3,P4 ,P5,P6, nrow = 2,labels = 'AUTO',label_size = 9)

ggsave(ADNI_epic_corr, file='~/Desktop/AD2/figure/CELL_Deconv/ADNI_corr_ad.pdf', width = 7.3,height =4.3 )



# ZIB
ZIB_ggscatter_plot <- function(ADNI_P,X_LAB,y_LAB){
  
  ZIB_P <-  ZIB_quantiseq$res_m
  ZIB_P$diagnosis <- coldata[rownames(ZIB_P),'Diagnosis']
  ZIB_P$MMSE  <- coldata[rownames(ZIB_P),'MMSE']
  ZIB_P$MoCA_B <- coldata[rownames(ZIB_P),'MoCA_B']
  ZIB_P$ACEIII_score <- coldata[rownames(ZIB_P),'ACEIII_score']
  
  ZIB_P$NL <- ZIB_P$Neutrophil/ZIB_P$`B cell`
  
  P1 <- ggscatter(subset(ZIB_P,diagnosis%in%c('NC','AD') ), 
                  x=X_LAB,y=y_LAB, size = 1, color="diagnosis" , palette = 'jco',# Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "spearman", label.sep = "\n"))
  return(P1)
}

P1 <- ZIB_ggscatter_plot(ADNI_P=ZIB_epic$res_m,X_LAB='Neutrophil',y_LAB='MoCA_B') 
P2 <- ZIB_ggscatter_plot(ADNI_P=ZIB_epic$res_m,X_LAB='Neutrophil',y_LAB='ACEIII_score') 
P3 <- ZIB_ggscatter_plot(ADNI_P=ZIB_epic$res_m,X_LAB='Neutrophil',y_LAB='MMSE') 
P4 <- ZIB_ggscatter_plot(ADNI_P=ZIB_epic$res_m,X_LAB='B cell',y_LAB='MMSE') 
P5 <- ZIB_ggscatter_plot(ADNI_P=ZIB_epic$res_m,X_LAB='B cell',y_LAB='MoCA_B') 
P6 <- ZIB_ggscatter_plot(ADNI_P=ZIB_epic$res_m,X_LAB='B cell',y_LAB='ACEIII_score') 

plot_grid(P1,P2,P3,P4 ,P5,P6, nrow = 2)

P1 <- ZIB_ggscatter_plot(ADNI_P=ZIB_quantiseq$res_m,X_LAB='Neutrophil',y_LAB='MoCA_B') 
P2 <- ZIB_ggscatter_plot(ADNI_P=ZIB_quantiseq$res_m,X_LAB='Neutrophil',y_LAB='ACEIII_score') 
P3 <- ZIB_ggscatter_plot(ADNI_P=ZIB_quantiseq$res_m,X_LAB='Neutrophil',y_LAB='MMSE') 
P4 <- ZIB_ggscatter_plot(ADNI_P=ZIB_quantiseq$res_m,X_LAB='B cell',y_LAB='MMSE') 
P5 <- ZIB_ggscatter_plot(ADNI_P=ZIB_quantiseq$res_m,X_LAB='B cell',y_LAB='MoCA_B') 
P6 <- ZIB_ggscatter_plot(ADNI_P=ZIB_quantiseq$res_m,X_LAB='B cell',y_LAB='ACEIII_score') 
plot_grid(P1,P2,P3,P4 ,P5,P6, nrow = 2)





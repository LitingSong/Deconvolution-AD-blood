# Deconvolution using Alzmap (spatial sequencing data of AD model mouse)

library(homologene)
library(limma)

# Expression matrix

ALZMAP <- readr::read_csv('~/Downloads/GSE152506_logCPM_counts.txt')
ALZMAP <- as.data.frame(ALZMAP)
rownames(ALZMAP) <- ALZMAP$X1
ALZMAP <- ALZMAP[,-1]

ALZ_homo_gene <- homologene(colnames(ALZMAP), inTax = 10090, outTax = 9606)
ALZMAP <- ALZMAP[,ALZ_homo_gene$`10090`]
ALZMAP <- t(ALZMAP)
ALZMAP <- avereps( ALZMAP, ALZ_homo_gene$`9606` )

atp11b <- read.table('~/Desktop/atp11b/GSE152506_atp11b.txt', sep=',', row.names = 1, header = T)
plaque <- read.table('~/Desktop/atp11b/plaqueScore.txt', sep=',', row.names = 2,header = T)

spot_meta <- read.table('~/Desktop/atp11b/spot_metadata(1).tsv',sep='\t', header = T)
rownames(spot_meta) <- spot_meta$Spot
spot_meta$Atp11b <- atp11b[rownames(spot_meta),1]
spot_meta$plaque <- plaque[rownames(spot_meta),'plaque_index']
spot_meta$Diagnosis <- ifelse(spot_meta$Genotype=='AD','AD','NC')
spot_meta$Diagnosis <- factor(spot_meta$Diagnosis, levels = c('NC','AD'))
spot_meta <- spot_meta[colnames(ALZMAP),]

save(ALZMAP, spot_meta, file='~/Desktop/AD2/data/ALZMAP.RData')


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


# # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # #
# # ## # ##                       EPIC                     ## # ## # ## # ## # ## # ## # #
# # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # #


ALZMAP_epic <- Deconv_immune(GEM=2^ALZMAP, dataset='ALZMAP',meta_data=spot_meta, array=F, deconvol_method='epic')

ALZMAP_epic_res <- ALZMAP_epic$res_p

# # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # #
# # ## # ##                       quantiseq                     ## # ## # ## # ## # ## # ## # #
# # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # #

ALZMAP_quantiseq_3 <- Deconv_immune(GEM=2^ALZMAP, dataset='ALZMAP',meta_data=subset(spot_meta,Age=='3'), array=F, deconvol_method='quantiseq')
ALZMAP_quantiseq_6 <- Deconv_immune(GEM=2^ALZMAP, dataset='ALZMAP',meta_data=subset(spot_meta,Age=='6'), array=F, deconvol_method='quantiseq')
ALZMAP_quantiseq_12 <- Deconv_immune(GEM=2^ALZMAP, dataset='ALZMAP',meta_data=subset(spot_meta,Age=='12'), array=F, deconvol_method='quantiseq')
ALZMAP_quantiseq_18 <- Deconv_immune(GEM=2^ALZMAP, dataset='ALZMAP',meta_data=subset(spot_meta,Age=='18'), array=F, deconvol_method='quantiseq')


ALZMAP_quantiseq_res <- rbind(ALZMAP_quantiseq_3$res_p,ALZMAP_quantiseq_6$res_p,
                              ALZMAP_quantiseq_12$res_p, ALZMAP_quantiseq_18$res_p)


# # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # #
# # ## # ##                       mcp_counter                     ## # ## # ## # ## # ## # ## #
# # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # ## # #

ALZMAP_mcp_counter <- Deconv_immune(GEM=ALZMAP, dataset='ALZMAP', meta_data=spot_meta, array=F, deconvol_method='mcp_counter')
ALZMAP_mcp_res <- ALZMAP_mcp_counter$res_p

#xx <- deconvolute_timer(2^ALZMAP[,rownames(subset(spot_meta,Age=='6'))], indications=rep('gbm',nrow((subset(spot_meta,Age=='6'))) ))
xx <- deconvolute_timer(2^ALZMAP[,rownames(subset(spot_meta,Age=='3'))], indications=rep('gbm',nrow((subset(spot_meta,Age=='6'))) ))
xx <- deconvolute_timer(2^ALZMAP[,rownames(subset(spot_meta,Age=='6'))], indications=rep('gbm',nrow((subset(spot_meta,Age=='6'))) ))
xx <- deconvolute_timer(2^ALZMAP[,rownames(subset(spot_meta,Age=='12'))], indications=rep('gbm',nrow((subset(spot_meta,Age=='6'))) ))
xx <- deconvolute_timer(2^ALZMAP[,rownames(subset(spot_meta,Age=='18'))], indications=rep('gbm',nrow((subset(spot_meta,Age=='6'))) ))

# Comparison of the absolute abundance of neutrophils between AD and NC model mice in the hippocampus and brain cortex of different ages
plot_pox <- function(ALZMAP_res){
  
  ALZMAP_res <- merge(ALZMAP_res, spot_meta, by.x='sample',by.y='Spot')
  
  ALZMAP_res$Age <- as.character(ALZMAP_res$Age )
  ALZMAP_res$Age <- factor(ALZMAP_res$Age, levels = c('3','6','12','18'))
  
  
  yl <- 'Score'
  if(max(ALZMAP_res$value) <= 1){ALZMAP_res$value <- ALZMAP_res$value*100; yl='Percentage (%)'}
  
  ggplot(subset(ALZMAP_res,  !is.na(Level_01) & cell_type=='Neutrophil'), 
         aes(x=diagnosis, y= value)) +
    geom_violin(alpha=0.3,width=0.5,size=0.3)+
    geom_boxplot(width=0.2,aes(fill=diagnosis),outlier.shape ='',size=0.3)+
    scale_fill_manual(values=color_code[ c('NC','AD')])+
    scale_color_manual(values=color_code[  c('NC','AD')])+
    #scale_fill_jco()+scale_color_jco()+
    facet_wrap(~Age,nrow = 1)+ylab(yl)+theme_bw()+
    stat_compare_means(label.y = 5,comparisons = list(c('NC','AD')),label = "p.label", method = "wilcox.test",
                       tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                       map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
    theme(legend.position = '',panel.spacing = unit(0.1, "lines"),text = element_text(size =8))  #+ylim(c(0,5))
  
}

plot_pox(ALZMAP_quantiseq_res)
plot_pox(ALZMAP_epic_res)
plot_pox(ALZMAP_mcp_counter$res_p)


dev.print(pdf, file='~/Desktop/AD2/figure/CELL_Deconv/ALZMAP_age.pdf', width=6, height=3)





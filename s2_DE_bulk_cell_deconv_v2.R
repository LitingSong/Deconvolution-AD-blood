# differential expression analysis at bulk level
# NC vs AD


library(limma)
library(DESeq2)
require(DOSE)
require(clusterProfiler)
require(org.Hs.eg.db)
library(Vennerable)


# limma for array
# ANM1_norm, ANM2_norm, ZIB, adni
# coldata, label_DXSUM_ec, anm1_meta, anm2_meta
load('~/Desktop/AD2/data/DEs.RData')

load('~/Desktop/AD2/data/exp_meta.RData')

ensg_sym <- read.table('~/Desktop/AD2/data/ensg_sym_autosome.txt',header = T, stringsAsFactors = F)

ZIB <- as.data.frame(gene_exp$vsd)

ZIB$gene_symbol <- ensg_sym[rownames(ZIB),'gene_symbol']
ZIB <- aggregate( ZIB[,-ncol(ZIB)], by = list(ZIB$gene_symbol),FUN = mean)
rownames(ZIB) <- ZIB$Group.1
ZIB <- ZIB[,-1]
colnames(ZIB)==rownames(coldata)

#
coldata$Diagnosis[coldata$Diagnosis=='SCD'] <- 'NC'
label_DXSUM_ec$Diagnosis[label_DXSUM_ec$Diagnosis=='SCD'] <- 'NC'
label_DXSUM_ec$Age[is.na(label_DXSUM_ec$Age)] <- label_DXSUM_ec$AGE[is.na(label_DXSUM_ec$Age)]

# limma
de_limma <- function(GEM,meta_data, dataset){
  
  meta_data <- subset(meta_data,Diagnosis%in%c('NC','MCI','AD'))
  meta_data$Gender <- ifelse(meta_data$Gender%in%c('Female','F'),0,1)
  GEM <- GEM[,rownames(meta_data)]

  design <- model.matrix(~0+Diagnosis+Age+Gender,meta_data )
  #design <- model.matrix(~0+Diagnosis,meta_data )
  
  fit <- lmFit(GEM,design)
  cont.matrix <- makeContrasts(P1="DiagnosisAD-DiagnosisNC",P2="DiagnosisMCI-DiagnosisNC", levels=design)#P2="DXMCI-DXNC",#P3="DXAD-DXNC"
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit3  <- eBayes(fit2)
  
  tT_ad <- topTable(fit3, coef=1, n=Inf, adjust="BH")
  tT_ad = subset(tT_ad, select=c("adj.P.Val","P.Value","logFC"))
  colnames(tT_ad)=c("FDR","P.Value","logFC")
  tT_ad$gene <- rownames(tT_ad)
  tT_ad$dataset <- dataset
  
  tT_ad$dirc <- 'up'
  tT_ad$dirc[tT_ad$logFC<0] <- 'down'
  
  tT_mci <- topTable(fit3, coef=2, n=Inf, adjust="BH")
  tT_mci = subset(tT_mci, select=c("adj.P.Val","P.Value","logFC"))
  colnames(tT_mci)=c("FDR","P.Value","logFC")
  tT_mci$gene <- rownames(tT_mci)
  tT_mci$dataset <- dataset
  
  tT_mci$dirc <- 'up'
  tT_mci$dirc[tT_mci$logFC<0] <- 'down'
  
  tT <- list(DE_AD=tT_ad,DE_MCI=tT_mci)
  
  return(tT)
}

# DE analysis
DE_ADNI <- de_limma(GEM=adni, dataset='ADNI', meta_data=label_DXSUM_ec)
DE_ANM1 <- de_limma(GEM=ANM1_norm, dataset='ANM1', meta_data=anm1_meta)
DE_ANM2 <- de_limma(GEM=ANM2_norm, dataset='ANM2', meta_data=anm2_meta)
DE_ZIB <- de_limma(GEM=ZIB, dataset='ZIB', meta_data=coldata)


fdr <- 0.01
logfc <- log2(1)
DeSig_ADNI <- subset(DE_ADNI$DE_AD, P.Value < fdr & abs(logFC) >  logfc)
DeSig_ANM1 <- subset(DE_ANM1$DE_AD, P.Value < fdr & abs(logFC) >  logfc)
DeSig_ANM2 <- subset(DE_ANM2$DE_AD, P.Value < fdr & abs(logFC) >  logfc)
DeSig_ZIB <- subset(DE_ZIB$DE_AD, P.Value < fdr & abs(logFC) >  logfc)
#DeSig_ANM2$gene[DeSig_ANM2$gene=='HLA-A29.1'] <- 'HLA-A'

#  MCI vs NC
# DeSig_ADNI <- subset(DE_ADNI$DE_MCI, P.Value < fdr & abs(logFC) >  logfc)
# DeSig_ANM1 <- subset(DE_ANM1$DE_MCI, P.Value < fdr & abs(logFC) >  logfc)
# DeSig_ANM2 <- subset(DE_ANM2$DE_MCI, P.Value < fdr & abs(logFC) >  logfc)
# DeSig_ZIB <- subset(DE_ZIB$DE_MCI, P.Value < fdr & abs(logFC) >  logfc)

#write.table(DeSig_ADNI,file='~/Desktop/DeSig_ADNI.txt',sep='\t',quote = F)


library(VennDiagram)
de_gene <- venn.diagram(list(ADNI=DeSig_ADNI$gene,
                             ANM1=DeSig_ANM1$gene,
                             ANM2=DeSig_ANM2$gene,
                             ZIB=DeSig_ZIB$gene),fill =  pal_jco("default")(4), cex=0.5,cat.cex = 0.5,
                        alpha = 0.8, filename = NULL,col='white')
plot_grid(de_gene)
#dev.print(pdf, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/de_gene_venn.pdf',width = 3, height = 3)

bulk_DEs <- rbind(DeSig_ADNI, DeSig_ANM1,DeSig_ANM2,DeSig_ZIB)

bulk_DEs <- bulk_DEs%>%arrange(dataset,dirc)

write.csv(bulk_DEs, file='~/Desktop/de/FC12P5/bulk_DEs.csv')

library(gprofiler2)

pat_Up <- (gprofiler2::gost(list(ADNI=DeSig_ADNI$gene[DeSig_ADNI$logFC>0],
                                 ANM1=DeSig_ANM1$gene[DeSig_ANM1$logFC>0],
                                 ANM2=DeSig_ANM2$gene[DeSig_ANM2$logFC>0],
                                 ZIB=DeSig_ZIB$gene[DeSig_ZIB$logFC>0]), 
                            sources=c('GO:BP'), correction_method='fdr'))$result %>% subset(term_size<500 & term_size> 10)

pat_Down <- (gprofiler2::gost(list(ADNI=DeSig_ADNI$gene[DeSig_ADNI$logFC<0],
                                   ANM1=DeSig_ANM1$gene[DeSig_ANM1$logFC<0],
                                   ANM2=DeSig_ANM2$gene[DeSig_ANM2$logFC<0],
                                   ZIB=DeSig_ZIB$gene[DeSig_ZIB$logFC<0]), 
                              sources=c('GO:BP'), correction_method='fdr'))$result %>% subset(term_size<500 & term_size> 10)


top_path_up <- pat_Up%>%group_by(query)%>%slice_head(n = 3, order_by = p_value)
top_path_up <- subset(pat_Up, term_name%in%top_path_up$term_name)
top_path_up$term_name <- factor(top_path_up$term_name,levels = unique(top_path_up$term_name))
top_path_up$gene_ratio <- top_path_up$intersection_size/top_path_up$query_size
colnames(top_path_up)[6] <- 'GeneCount'

top_path_Down <- pat_Down%>%group_by(query)%>%slice_head(n = 3, order_by = p_value)
top_path_Down <- subset(pat_Down, term_name%in%top_path_Down$term_name)
top_path_Down$term_name <- factor(top_path_Down$term_name,levels = unique(top_path_Down$term_name))
top_path_Down$gene_ratio <- top_path_Down$intersection_size/top_path_up$query_size
colnames(top_path_Down)[6] <- 'GeneCount'

path_up <- ggplot(top_path_up, aes(x=query,y=term_name,fill=-log10(p_value)))+
  geom_point(aes(size=GeneCount),shape=21)+theme_bw()+
  scale_fill_gradient2(low='#0073C2FF',mid='grey',high="#af1319",name='-log10(FDR)') +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 32) )+
  theme(text = element_text(size = 7),legend.key.size  = unit(0.3,'cm'),
        legend.margin = margin(c(0,0,0,0)),axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6)) +  
  xlab('')+ylab('')


#dev.print(pdf, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/bulk_DE_path_up_0.05.pdf',width = 4.2,height = 3 )

path_down <-ggplot(top_path_Down, aes(x=query,y=term_name,fill=-log10(p_value)))+
  geom_point(aes(size=GeneCount),shape=21)+theme_bw()+
  scale_fill_gradient2(low='#EFC000FF',mid='grey',high="#0073C2FF",name='-log10(FDR)') +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 32) )+
  theme(text = element_text(size = 7),legend.key.size  = unit(0.3,'cm'), 
        legend.margin = margin(c(0,0,0,0)),axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6)) +  
  xlab('')+ylab('')

#dev.print(pdf, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/bulk_DE_path_down_0.05.pdf',width = 4.2,height = 3 )


plot_grid(de_gene,nrow = 1,
          labels = 'A',label_size = 9)
dev.print(pdf, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/de_bulk_VENN.pdf',
          width = 1.7,height = 1.7 )


plot_grid(path_up,path_down, nrow = 1,rel_widths = c(2,2),scale = c(1,1),
          labels = c('B','C'),label_size = 9)

dev.print(pdf, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/de_bulk_ad.pdf',
          width = 6,height = 2.3 )


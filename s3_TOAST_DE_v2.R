# Differentially expressed genes (DEG) in AD at cellular level. 

library(TOAST)
load('~/Desktop/AD2/data/exp_meta.RData')
load('~/Desktop/AD2/data/blood_cell_deconv.RData')


coldata$Diagnosis[coldata$Diagnosis=='SCD'] <- 'NC'
label_DXSUM_ec$Diagnosis[label_DXSUM_ec$Diagnosis=='SCD'] <- 'NC'
label_DXSUM_ec$Age[is.na(label_DXSUM_ec$Age)] <- label_DXSUM_ec$AGE[is.na(label_DXSUM_ec$Age)]


ADNI_exp <- 2^adni
ANM1_exp <- 2^ANM1_norm
ANM2_exp <- 2^ANM2_norm


#lm6 in ciborsort
Ref <- read.table('~/Desktop/AD2/data/signature_rnaseq_geo60424_LM6.txt',sep='\t',header = T, row.names = 1)
TOAST_DE <- function(GEM,dataset,meta_data, Prop, deconvol_method){
  
  meta_data <- subset(meta_data,Diagnosis%in%c('NC','AD') )
  GEM <- GEM[,rownames(meta_data)]
  
  meta_data$Diag <- ifelse(meta_data$Diagnosis=='NC',0,1)
  meta_data$Gender <- ifelse(meta_data$Gender%in%c('Female','F'),0,1)
  
  design <- data.frame(disease = as.factor(meta_data$Diag), Gender=meta_data$Gender, Age=meta_data$Age)
  
  if(deconvol_method=='epic'){Prop <-  Prop[rownames(meta_data),-7][,c(1,3,2,5,4,6)] }
  if(deconvol_method=='quantiseq'){Prop <-  Prop[rownames(meta_data),c(-2,-3,-9:-12)][,c(1,6,5,4,2,3)] }
  #if(deconvol_method=='quantiseq'){Prop <-  Prop[rownames(meta_data),c(-11:-12)] }
  Prop <-  Prop[rownames(meta_data),]
  
  colnames(Prop) <- colnames(Ref) 
  
  Design_out <- makeDesign(design, Prop)
  #Design_out$design$Condition_num <- factor(ifelse(Design_out$design$disease == "AD", 1, 0))
  
  fitted_model <- fitModel(Design_out, as.matrix(GEM))
  fitted_model$all_coefs # list all phenotype names
  fitted_model$all_cell_types # list all cell type names
  
  Ref <- Ref[intersect(rownames(Ref),rownames(GEM) ),]
  
  res_tables <- c()
  for(cell_ty in colnames(Ref) ){
    
    res_table <- csTest(fitted_model, coef = "disease", 
                        cell_type = cell_ty , contrast_matrix = NULL)
    
    res_table$gene <- rownames(res_table)
    res_table$cell_type <- cell_ty
    res_table$dataset <- dataset
    res_tables <- rbind(res_tables, res_table)
    
  }
  
  #sig_res_tables <- subset(res_tables,p_value < 0.01)
  
  return(res_tables)
  
}


# DEG of each cell type in AD
ADNI_epic_TOAST_DE <- TOAST_DE(GEM=ADNI_exp, dataset='ADNI',meta_data=label_DXSUM_ec, deconvol_method='epic', Prop = ADNI_epic$res_m)
ZIB_epic_TOAST_DE <- TOAST_DE(GEM=ZIB, dataset='ZIB', meta_data=coldata, deconvol_method='epic', Prop = ZIB_epic$res_m)
ANM1_epic_TOAST_DE <- TOAST_DE(GEM=ANM1_exp, dataset='ANM1', meta_data=anm1_meta, deconvol_method='epic', Prop = ANM1_epic$res_m)
ANM2_epic_TOAST_DE <- TOAST_DE(GEM=ANM2_exp, dataset='ANM2', meta_data=anm2_meta, deconvol_method='epic', Prop = ANM2_epic$res_m)

# Overlap between cell-type aware DEGs and buk DEGs 
sig_ADNI_epic_TOAST <- subset(ADNI_epic_TOAST_DE,  p_value < 0.01 & abs(effect_size) > log2(1))
sig_ADNI_epic_TOAST$dir <- ifelse(sig_ADNI_epic_TOAST$effect_size<0,'down','up')
write.table(sig_ADNI_epic_TOAST, file='~/Desktop/de/FC12P5/sig_ADNI_epic_TOAST.txt',sep='\t',quote = F)


library(VennDiagram)
cell_types <- unique(sig_ADNI_epic_TOAST$cell_type)

down_venn <- list()
for(ct in cell_types){
  
  down_gene <- venn.diagram(list(bulk =  DeSig_ADNI$gene[DeSig_ADNI$dirc=='down'],
                                 cell_type = sig_ADNI_epic_TOAST$gene[sig_ADNI_epic_TOAST$cell_type==ct & sig_ADNI_epic_TOAST$dir=='down' ]),
                            fill =  pal_jco("default")(4)[c(1,4)], cex=0.5,cat.cex = 0,
                            alpha = 0.8, filename = NULL,col='white')
  down_venn[[ct]]  <- down_gene
}
plot_grid(down_venn$B.cells,down_venn$CD8.T.cells,down_venn$CD4.T.cells,down_venn$NK.cells,down_venn$Monocytes,down_venn$Neutrophils, labels =  cell_types,label_size = 8)

up_venn <- list()
for(ct in cell_types){
  
  up_gene <- venn.diagram(list(bulk =  DeSig_ADNI$gene[DeSig_ADNI$dirc=='up'],
                                 cell_type=sig_ADNI_epic_TOAST$gene[sig_ADNI_epic_TOAST$cell_type==ct & sig_ADNI_epic_TOAST$dir=='up' ]),
                            fill =  pal_jco("default")(4)[c(1,4)], cex=0.5,cat.cex = 0,
                            alpha = 0.8, filename = NULL,col='white')
  up_venn[[ct]]  <- up_gene
}
plot_grid(up_venn$B.cells,up_venn$CD8.T.cells,up_venn$CD4.T.cells,up_venn$NK.cells,up_venn$Monocytes,up_venn$Neutrophils, labels =  cell_types,label_size = 8)

all_venn <- list()
p_hyper <- c()

for(ct in cell_types){
  
  up_gene <- venn.diagram(list(bulk =  DeSig_ADNI$gene,
                               cell_type=sig_ADNI_epic_TOAST$gene[sig_ADNI_epic_TOAST$cell_type==ct  ]),
                          fill =  pal_jco("default")(4)[c(1,4)], cex=0.5,cat.cex = 0,
                          alpha = 0.8, filename = NULL,col='white')
  up_venn[[ct]]  <- up_gene
  
  a <- length(DeSig_ADNI$gene)
  b <- length(sig_ADNI_epic_TOAST$gene[sig_ADNI_epic_TOAST$cell_type==ct  ])
  inter <- length(intersect(DeSig_ADNI$gene, sig_ADNI_epic_TOAST$gene[sig_ADNI_epic_TOAST$cell_type==ct  ]))

  p_hyper <-  signif(c(p_hyper,phyper(inter-1, a, 20092-a, b, lower.tail=F)),2)
}

plot_grid(up_venn$B.cells,up_venn$CD8.T.cells,up_venn$CD4.T.cells,up_venn$NK.cells,up_venn$Monocytes,up_venn$Neutrophils, 
          labels =  paste(cell_types,' (',p_hyper,')',sep=''),label_size = 8)

dev.print(pdf, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/venn_bulk_celltp_de.pdf',
          width = 7.3,height = 4 )


#AD_GENE <- read.table('~/Desktop/AD2/data/AD.combined2.txt',header = F, stringsAsFactors = F)[,1]
AD_GENE <- read.table('/Users/songlt/Desktop/AD2/data/AD_LIST/Combined.txt',header = F, stringsAsFactors = F)[,1]

sig_ADNI_epic_TOAST$AD <- 'N'
sig_ADNI_epic_TOAST$AD[sig_ADNI_epic_TOAST$gene%in%AD_GENE] <- 'Y'
xtabs(~AD+cell_type,sig_ADNI_epic_TOAST)


# Neutrophils DEGs

sig_ADNI_Neutrophils <- subset(sig_ADNI_epic_TOAST,cell_type== 'Neutrophils')
sig_ADNI_Neutrophils$dir <- ifelse(sig_ADNI_Neutrophils$effect_size<0,'down','up')
#write.table(sig_ADNI_Neutrophils[,c('gene','AD','dir')], file='~/Desktop/sig_ADNI_Neutrophils.txt',sep='\t', row.names = F,quote = F)

# PPI module (0.5)
string_ppi <- read.csv('~/Desktop/AD2/data/de/FC1P1/ADNI_neutro_FIG3/Enrichment_PPI/ADNI_MCODE_ALL_PPIColorByCluster default node.csv')#[,c('query.term','MCODE_Cluster','MCODE_Node_Status')]

# pathway (metascape)

cluster_patw <- c()
for(cl in 1:5){
  #ptw <- gProfileR::gprofiler(string_ppi$query.term[string_ppi$MCODE_Cluster==paste('Cluster',cl)], hier_filtering = 'moderate',max_set_size = 500, src_filter = c('GO:BP'),correction_method = 'fdr')
  ptw <- gProfileR::gprofiler(string_ppi$query.term[grep(paste('Cluster',cl), string_ppi$MCODE_Cluster)], hier_filtering = 'moderate',max_set_size = 500, src_filter = c('GO:BP'),correction_method = 'fdr')
  
  ptw$cluster <- paste('Cluster',cl)
  cluster_patw <- rbind(ptw, cluster_patw)
}


# the enrichment of modular clusters for AD risk gene sets

fisher_p<- c()
for(cl in (1:6)){
  #for (AD_gene_list in dir('/Users/songlt/Desktop/AD2/data/disease_gene_sets20190906',pattern = 'AD_',full.names = F)){
    for (AD_gene_list in dir('/Users/songlt/Desktop/AD2/data/AD_LIST/',pattern = '.txt',full.names = F)){
      
    AD_GENE <- read.table(paste0('/Users/songlt/Desktop/AD2/data/AD_LIST/',AD_gene_list),header = F, stringsAsFactors = F)[,1]
    
    AD_adni_gene <- intersect(rownames(adni),AD_GENE)
    
    aa <- length(intersect(AD_adni_gene,string_ppi$Symbol[string_ppi$MCODE_CLUSTER_ID==cl] ))
    ab <- length(setdiff(AD_adni_gene,string_ppi$Symbol[string_ppi$MCODE_CLUSTER_ID==cl] ))
    ba <- length(setdiff(string_ppi$Symbol[string_ppi$MCODE_CLUSTER_ID==cl],AD_adni_gene ))
    bb <- nrow(adni)- aa-ab-ba
    fis <- fisher.test(matrix(c(aa,ab,ba,bb),nrow=2))
    
    fisher_p <- as.data.frame(rbind(fisher_p,c(fis$p.value,fis$estimate, str_split_fixed(AD_gene_list,'\\.',3)[1],cl)))
  }
}

colnames(fisher_p) <- c('p_value','Odds_Ratio','dataset','cluster')
fisher_p$p_value <- as.numeric(fisher_p$p_value)
fisher_p$Odds_Ratio <- as.numeric(fisher_p$Odds_Ratio)
fisher_p$signif <- factor(ifelse(fisher_p$p_value<0.05,'Y','N'),levels = c('Y','N'))
#fisher_p$cluster <- factor(fisher_p$cluster, levels = paste('Cluster',(1:10)))

fisher_point <- ggplot(fisher_p, aes(y=dataset,x=cluster))+
  geom_point(aes(size=Odds_Ratio,fill=-log10(p_value), color=signif),shape=22)+theme_bw()+
  scale_fill_gradientn(colours= c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404'),name="-log10(P-value)" )+
  scale_color_manual(values = c('black','grey'))+
  theme(axis.text.x = element_text(angle=0,size=9),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'lines'), 
        legend.key.height = unit (0.3, 'cm'),legend.key.width = unit (0.3, 'cm'),
        legend.title = element_text(size=7), 
        axis.text.y = element_text(size=8),
        text = element_text(size = 9)) +   xlab('Cluster')+ylab('')

ggsave(fisher_point, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/ADNI_cluster_fisher.pdf',width = 4.3, height = 4)



# the DESs of Neutrophils in AD across diverse datasets

Neu_venn <- venn.diagram(list(ADNI=subset(ADNI_epic_TOAST_DE, cell_type== 'Neutrophils' & p_value < 0.01 )[,'gene'],
                      ANM1=subset(ANM1_epic_TOAST_DE, cell_type== 'Neutrophils' & p_value < 0.01)[,'gene'],
                      ANM2=subset(ANM2_epic_TOAST_DE, cell_type== 'Neutrophils' & p_value < 0.01)[,'gene'],
                      ZIB=subset(ZIB_epic_TOAST_DE, cell_type== 'Neutrophils' & p_value < 0.01)[,'gene']),
                      fill =  pal_jco("default")(4), cex=0.5,cat.cex = 0.5,
                      alpha = 0.8, filename = NULL,col='white')

plot_grid(Neu_venn)

dev.print(pdf, file='/Users/songlt/Desktop/AD2/figure/CELL_Deconv/Neu_venn.pdf',width = 2.2,height = 2.2 )



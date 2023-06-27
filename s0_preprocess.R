# data preprocess 
# demographic data 

library(stringr)
library(ggplot2)
library(ggpubr)
library(cowplot)

##################################################################################################################################
#############                        part1        data preprocess                                                         ########
##################################################################################################################################

# meta data and gene expression matrix

## ZIB
load('~/Desktop/AD2/data/DEs.RData')
ensg_sym <- read.table('~/Desktop/AD2/data/ensg_sym_autosome.txt',header = T, stringsAsFactors = F)

#ZIB <- as.data.frame(gene_exp$vsd)
ZIB <- as.data.frame(gene_exp$tpm)

ZIB$gene_symbol <- ensg_sym[rownames(ZIB),'gene_symbol']
ZIB <- aggregate( ZIB[,-ncol(ZIB)], by = list(ZIB$gene_symbol),FUN = mean)
rownames(ZIB) <- ZIB$Group.1
ZIB <- ZIB[,-1]
colnames(ZIB)==rownames(coldata)

#write.table(ZIB, file='~/Desktop/AD2/data/ZIB.txt',sep='\t',row.names = T, quote = F)

##  adni 


ADNI_expre <- read.csv('/Users/songlt/Desktop/ADNI/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv',header = F)

# average for same gene
ADNI_exp <- ADNI_expre
ADNI_exp <- subset(ADNI_exp,)
colnames(ADNI_exp) <- ADNI_exp[3,]
ADNI_exp <- ADNI_exp[-1:-9,c(-1,-2,-ncol(ADNI_exp))]
colnames(ADNI_exp)[1] <- 'gene'
ADNI_exp <- aggregate(.~gene,ADNI_exp,function(x){mean(as.numeric(x))})
rownames(ADNI_exp) <- ADNI_exp[,1]
ADNI_exp <- ADNI_exp[,-1]

write.table(ADNI_exp, file='/Users/songlt/Documents/AD/ADNI/ADNI_recode1.txt',sep = '\t',quote = F)
adni <- read.table('/Users/songlt/Documents/AD/ADNI/ADNI_recode1.txt', row.names = 1,check.names = F,sep='\t',header = T)

load("/Users/songlt/Desktop/ADNI/ADNIMERGE/data/ptdemog.rdata")
ptdemog <- unique(ptdemog[,c('RID','PTEDUCAT','PTETHCAT','PTGENDER','AGE')])
ptdemog <- ptdemog[!duplicated(ptdemog$RID),]
row.names(ptdemog) <- as.character(ptdemog$RID)

label_DXSUM_ec <- read.table('/Users/songlt/Documents/AD/ADNI/label_DXSUM_ec.txt',sep='\t',header = T,comment.char = '', row.names = 3)
label_DXSUM_ec$Diagnosis <- label_DXSUM_ec$DX
label_DXSUM_ec$age <- as.numeric(str_split_fixed(label_DXSUM_ec$USERDATE,'-',2)[,1])-label_DXSUM_ec$PTDOBYY
label_DXSUM_ec$RID.x <- as.character(as.numeric(str_split_fixed(label_DXSUM_ec$SubjectID,'_',3)[,3]))
label_DXSUM_ec <- cbind(label_DXSUM_ec[,-which(colnames(label_DXSUM_ec)=='PTGENDER')],ptdemog[label_DXSUM_ec$RID.x,])
label_DXSUM_ec <- label_DXSUM_ec[colnames(adni),]


colnames(adni)==rownames(label_DXSUM_ec)


# ANM DATA

##  ANM1

anm1_meta <- as.data.frame(t(read.table('~/Desktop/AD2/data/ANM/GSE63060_series_matrix.txt',sep = '\t',comment.char = '!',header = T,row.names = 1)[1:4,]))
anm1_meta$Diagnosis <- gsub('status: ','',anm1_meta$Diagnosis)
anm1_meta$Age <- as.numeric(gsub('age: ','',anm1_meta$Age))
anm1_meta$Gender <- gsub('gender: ','',anm1_meta$Gender)
anm1_meta$Diagnosis[anm1_meta$Diagnosis=='CTL'] <- 'NC'

#ANM1_norm <- read.table('~/Desktop/AD2/data/ANM/GSE63060_normalized.txt',sep='\t',row.names = 1)
ANM1_norm <- read.table('~/Desktop/AD2/data/ANM/GSE63060_series_matrix.txt',sep = '\t',comment.char = '!',header = T,row.names = 1)[-1:-4,]
samples_ANM1 <- colnames(ANM1_norm)
ANM1_norm <- t(apply(ANM1_norm,1, as.numeric))
colnames(ANM1_norm) <- samples_ANM1
ANM1_norm <- as.data.frame(ANM1_norm)

anm1_gene_sym <- read.table('~/Desktop/AD2/data/ANM/GPL6947-13512.txt',sep='\t',quote = '"',comment.char = '',header = T, row.names = 1)

ANM1_norm$gene_symbol <- anm1_gene_sym[rownames(ANM1_norm),'ILMN_Gene']
ANM1_norm <- aggregate( ANM1_norm[,-ncol(ANM1_norm)], by = list(ANM1_norm$gene_symbol),FUN = mean)
rownames(ANM1_norm) <- ANM1_norm$Group.1
ANM1_norm <- ANM1_norm[,-1]


## ANM2
anm2_meta <- as.data.frame(t(read.table('~/Desktop/AD2/data/ANM/GSE63061_series_matrix.txt',sep = '\t',comment.char = '!',header = T,row.names = 1)[1:4,]))
anm2_meta$Diagnosis <- gsub('status: ','',anm2_meta$Diagnosis)
anm2_meta$Age <- as.numeric(gsub('age: ','',anm2_meta$Age))
anm2_meta$Gender <- gsub('gender: ','',anm2_meta$Gender)
anm2_meta$Diagnosis[anm2_meta$Diagnosis=='CTL'] <- 'NC'

#ANM2_norm <- read.table('~/Desktop/AD2/data/ANM/GSE63061_normalized.txt',sep='\t',row.names = 1)
ANM2_norm <- read.table('~/Desktop/AD2/data/ANM/GSE63061_series_matrix.txt',sep = '\t',comment.char = '!',header = T,row.names = 1)[-1:-4,]
samples_ANM2 <- colnames(ANM2_norm)
ANM2_norm <- t(apply(ANM2_norm,1, as.numeric))
colnames(ANM2_norm) <- samples_ANM2
ANM2_norm <- as.data.frame(ANM2_norm)
anm2_gene_sym <- read.table('~/Desktop/AD2/data/ANM/GPL10558-50081.txt',sep='\t',quote = '"',comment.char = '',header = T, row.names = 1)

ANM2_norm$gene_symbol <- anm2_gene_sym[rownames(ANM2_norm),'ILMN_Gene']
ANM2_norm <- aggregate( ANM2_norm[,-ncol(ANM2_norm)], by = list(ANM2_norm$gene_symbol),FUN = mean)
rownames(ANM2_norm) <- ANM2_norm$Group.1
ANM2_norm <- ANM2_norm[,-1]

save(ANM1_norm,ANM2_norm,anm2_meta,anm1_meta, file='~/Desktop/AD2/data/ANM/ANM.RData')

load('~/Desktop/AD2/data/ANM/ANM.RData')
#write.table(ANM1_norm, file='~/Desktop/AD2/data/ANM/ANM1_norm.txt',sep='\t',quote = F)
#write.table(ANM2_norm, file='~/Desktop/AD2/data/ANM/ANM2_norm.txt',sep='\t',quote = F)
#ANM1_norm,ANM2_norm,anm2_meta,anm1_meta

rownames(anm1_meta)==colnames( ANM1_norm)

rownames(anm2_meta)==colnames( ANM2_norm)

coldata$Gender <- ifelse(coldata$Gender=='F','Female','Male')

label_DXSUM_ec <- reshape::rename(label_DXSUM_ec,c(age='Age',PTGENDER='Gender'))

save(ANM1_norm, ANM2_norm, ZIB, adni, coldata, label_DXSUM_ec, anm1_meta, anm2_meta, file='~/Desktop/AD2/data/exp_meta.RData')

##################################################################################################################################
#############                        Part2.         demographic data analysis                                             ########
##################################################################################################################################



color_code <- c( "#979797","#c4c4c4", "#af1319")
names(color_code) <- c('NC','MCI','AD')


load('~/Desktop/AD2/data/exp_meta.RData')

coldata <- coldata[,c('Diagnosis','Gender','Age')]
coldata$dataset <- 'ZIB'

anm1_meta <- anm1_meta[,c('Diagnosis' ,'Gender','Age' )]
anm1_meta$dataset <- 'ANM1'

anm2_meta <- anm2_meta[,c('Diagnosis' ,'Gender','Age' )]
anm2_meta$dataset <- 'ANM2'

label_DXSUM_ec <- label_DXSUM_ec[,c('Diagnosis' ,'Gender','Age' )]
label_DXSUM_ec$dataset <- 'ADNI'

datasets_demog <-  rbind(coldata, anm1_meta, anm2_meta, label_DXSUM_ec)
datasets_demog <- subset(datasets_demog,Diagnosis%in%c('NC','MCI','AD'))
datasets_demog$Diagnosis <- factor(datasets_demog$Diagnosis, levels = c('NC','MCI','AD'))
datasets_demog$Gender <- ifelse(datasets_demog$Gender%in%c('M','Male'),'M','F')

# number of individuals
diag_bar1 <- ggplot(datasets_demog,aes(x=Diagnosis,fill=Diagnosis))+
  geom_bar(width = 0.7,position = 'dodge')+
  scale_fill_manual(values=color_code)+
  facet_grid(.~dataset) + 
  theme_bw()+
  geom_text(aes(label=after_stat(count)), vjust=0, stat = "count", size=2.5)+
  theme(text = element_text(size =8), legend.position = '')+
  ylab('Number of individuals')

diag_bar2 <- ggplot(datasets_demog,aes(x=dataset,fill=Diagnosis))+
  geom_bar(width = 0.7,position = 'stack')+
  scale_fill_manual(values=color_code)+
  theme_bw()+
  geom_text(aes(label=after_stat(count)),
            stat = "count", size=2.5, position='stack')+
  theme(text = element_text(size =8),legend.margin = margin(l = -0.2, unit='cm'))+
  ylab('Number of individuals')


# age distribution
age_distribution <- ggplot(datasets_demog,aes(x=Diagnosis,y=Age,fill=Diagnosis))+
  facet_wrap(~dataset,nrow = 1)+
  scale_fill_manual(values=color_code)+
  geom_boxplot(width=0.5,size=0.2,outlier.size = 0.2)+ylab('Age')+theme_bw()+
  stat_compare_means(comparisons = list(c('NC','AD'),c('MCI','AD') ,c('NC','MCI')),label = "p.signif", method = "wilcox.test",
                     tip.length = 0.01,bracket.size = 0.2,step.increase = 0.05,hide.ns = F,
                     map_signif_level	= c("***"=0.001, "**"=0.01, "*"=0.05,"."=0.1),size=2.5,family='sans',face='plain')+
  theme(panel.spacing = unit(0.1, "lines"),text = element_text(size =8), legend.position = '')


# sex proportion
demo_sex <- datasets_demog %>%
  group_by(dataset,Diagnosis,Gender) %>%
  dplyr::summarise(N=n()) %>%
  group_by(Diagnosis,dataset) %>%
  dplyr::mutate(freq = (N/sum(N)*100))


gender_stat <- data.frame(p=paste('P =',round(c(chisq.test(matrix(unlist(demo_sex[c(1:6),4]),nrow=2))$p.value
                                                ,chisq.test(matrix(unlist(demo_sex[c(7:12),4]),nrow=2))$p.value
                                                ,chisq.test(matrix(unlist(demo_sex[c(13:18),4]),nrow=2))$p.value
                                                ,chisq.test(matrix(unlist(demo_sex[c(19:24),4]),nrow=2))$p.value ),2)),
                          sig=c('*',"*",'ns','ns'),
                          dataset=c('ADNI','ANM1','ANM2','ZIB'))

sex_proportion <- ggplot(datasets_demog,aes(x=Diagnosis))+
  geom_bar(width = 0.7,position = position_fill(),aes(fill=Gender))+
  scale_fill_manual(values=as.character(color_code)[c(1,3)])+
  facet_grid(.~dataset) + 
  theme_bw()+
  ylab('Proportion')+
  theme(panel.spacing = unit(0.1, "lines"),text = element_text(size =8),
        legend.margin = margin(l = -0.2, unit='cm')) +
  geom_text(data = gender_stat,
    mapping = aes(x = 2, y = 1.03, label = sig),size=2.5)+ylim(c(0,1.05))


plot_grid(diag_bar2,age_distribution,sex_proportion,labels = 'AUTO',label_size = 9, rel_widths = c(1,1.25,1.45), nrow = 1)
dev.print(pdf, file= '~/Desktop/AD2/figure/CELL_Deconv/demog_stat.pdf', width = 8, height = 2.8)

plot_grid(plot_grid(diag_bar1,age_distribution, labels=c('A','B'),label_size = 9),
          plot_grid(sex_proportion,'',labels = 'C',label_size = 9, rel_widths = c(1.25,1), nrow = 1),nrow = 2)
dev.print(pdf, file= '~/Desktop/AD2/figure/CELL_Deconv/demog_stat2.pdf', width = 7.3, height = 5)





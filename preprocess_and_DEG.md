## reference 1 https://blog.csdn.net/weixin_48275332/article/details/124151180
## reference 2 https://mp.weixin.qq.com/s?__biz=MzAxMDkxODM1Ng==&mid=2247512042&idx=1&sn=acf55e1a4090661b75dca179213ecfd5&chksm=9b4bef51ac3c6647e387e8c9ab057430ea2c50487babf1f0eef7a62cdb7194d4f6a5077699db&mpshare=1&scene=24&srcid=0331JBdCxi66Ja7KQMqpBruA&sharer_sharetime=1648731156544&sharer_shareid=0cdcf284bc426c9af2f9246b33b7ac7f&ascene=14&devicetype=android-34&version=28003253&nettype=WIFI&abtest_cookie=AAACAA%3D%3D&lang=zh_CN&countrycode=CN&exportkey=n_ChQIAhIQO%2FfwQ4%2BRbdM4O%2BA%2BkjpkxRL4AQIE97dBBAEAAAAAAObtGKHduI0AAAAOpnltbLcz9gKNyK89dVj0k0z17HPsG86D0ggeZzYqReJVxvv9fISgiN%2BzqGPYZJBemSVvyV2qzjKxHSGJQTmwTWJPgUryGoewhtIvZJ0ZuQYTZqrNEV7HOTqmB03sqDr3KOgZ8uP1m6R9G0gYCXlgshbEU4B%2F0GJAHeZL4amUPn%2BgeVVOYLnBsuUqiD5gAB%2FG6Qt9YqsKjQLOoGwCcaIXp%2BW6J%2BhaOVUwOp51udq9S30obNdxG64LHitkU40PxS5Dn394QNJp4LG42ecOH%2FL3EG0xctmwkBleHIRfjf5v9A2y&pass_ticket=gmf50FYXyrIdG00kmhgKHDgB8cYmP4788s7gY%2FtZSNeP04zt6hq%2B21oSl5vX1QGz&wx_header=3

# Download data from UCSC Xena 
### htseq_counts, GDC_phenotype, survival


# 1. Load data
clinical = read.delim(paste0(proj,".GDC_phenotype.tsv.gz"),fill = T,header = T,sep = "\t")
surv = read.delim(paste0(proj,".survival.tsv"),header = T)
dat = read.table(paste0(proj,".htseq_counts.tsv.gz"),check.names = F,row.names = 1,header = T)


# 2. preprocess data
#逆转log
dat = as.matrix(2^dat - 1)
#  用apply转换为整数矩阵
exp = apply(dat, 2, as.integer)
rownames(exp) = rownames(dat)


# 行名ID转换：方法1(推荐)
library(stringr)
library(AnnoProbe)
rownames(exp) = str_split(rownames(exp),"\\.",simplify = T)[,1];head(rownames(exp))
re = annoGene(rownames(exp),ID_type = "ENSEMBL");head(re)
library(tinyarray)
exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")


# Filter out gene with less than 50% blank
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ]
## deprecate!!! Filter out gene with 0
# exp1 = exp[rowSums(exp)>0,]

# 3.Group
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
Group = factor(Group,levels = c("normal","tumor"))

# 4.Saving data
save(exp,Group,proj,clinical,surv,file = paste0(proj,".Rdata"))


# 5.DEG analysis
y <- DGEList(counts=exp,group=Group)
##Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
##Perform TMM normalization and transfer to CPM (Counts Per Million)
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)
count_norm<-as.data.frame(count_norm)


pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),Group)
  p=wilcox.test(gene~Group, data)$p.value
  return(p)
})
fdr=p.adjust(pvalues,method = "fdr")


conditionsLevel<-levels(Group)
dataCon1=count_norm[,c(which(Group==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(Group==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))


outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
rownames(outRst)=rownames(count_norm)
outRst=na.omit(outRst)
fdrThres=0.05


# 6.PLot
# Load the necessary library
library(ggplot2)

# Assign groups based on log2foldChange
outRst$Group <- ifelse(outRst$log2foldChange > 2, "Tumor",
                       ifelse(outRst$log2foldChange < -2, "Normal", "Intermediate"))

# Determine significance based on FDR
outRst$Significant <- ifelse(outRst$FDR < 0.05, "Yes", "No")

# Assign colors for plotting
outRst$Color <- ifelse(outRst$Significant == "Yes" & outRst$Group == "Tumor", "Tumor_Sig",
                       ifelse(outRst$Significant == "Yes" & outRst$Group == "Normal", "Normal_Sig", "Not_Sig"))

# Find the TIMP3 gene
timp3_data <- outRst[rownames(outRst) == "TIMP3", ]

# Define colors for the plot
colors <- c("Tumor_Sig" = "red", "Normal_Sig" = "blue", "Not_Sig" = "grey")

# Create the volcano plot
ggplot(outRst, aes(x = log2foldChange, y = -log10(FDR), color = Color)) +
  geom_point(alpha = 0.8, size = 1.5) +  # Points for all data
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 FDR", color = "Group") +
  theme(legend.position = "right") +
  geom_point(data = timp3_data, aes(x = log2foldChange, y = -log10(FDR)), color = "black", size = 3) +  # Black dot for TIMP3
  annotate("text", x = timp3_data$log2foldChange + 0.2, y = -log10(timp3_data$FDR) + 0.2,
           label = sprintf("TIMP3\nFDR: %.3f\nFC: %.3f", timp3_data$FDR, timp3_data$log2foldChange), 
           color = "black", size = 3, hjust = 0)

# test GSVA
# source("https://bioconductor.org/biocLite.R")
# biocLite("GSVA")
# biocLite("GSEABase")
# biocLite("GSVAdata")

# load the required libraries
library(gplots)
library(tidyverse)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)

# this is a set of signatures
data("c2BroadSets")
str(c2BroadSets)

# these are the data already normalized (ExpressionSet object)
data(leukemia)
leukemia_eset

#We carry out a non-specific filtering step by discarding the 50% of the probesets with smaller variability, probesets without Entrez ID annotation, probesets whose associated Entrez ID is duplicated in the annotation, and Affymetrix quality control probes:

filtered_eset <- nsFilter(leukemia_eset,
                          require.entrez=TRUE,
                          remove.dupEntrez=TRUE,
                          var.func=IQR, var.filter=TRUE,
                          var.cutoff=0.5,
                          filterByQuantile=TRUE,
                          feature.exclude="^AFFX")

filtered_eset

# still this is a expression set
leukemia_filtered_eset <- filtered_eset$eset

#################################################################################
# according to the vignette is possible to work with both the expression set and with a matrix of normalized expression.
# extract the data for the expression
exp <- exprs(leukemia_filtered_eset)
# save the pheno organization
pheno <- leukemia_filtered_eset@phenoData@data
####################################################################################

#The calculation of GSVA enrichment scores is performed in one single call to the gsva() function. However, one should take into account that this function performs further non-specific filtering steps prior to the actual calculations. On the one hand, it matches gene identifiers between gene sets and gene expression values. On the other hand, it discards gene sets that do not meet minimum and maximum gene set size requirements specified with the arguments min.sz and max.sz, respectively, which, in the call below, are set to 10 and 500 genes. Because we want to use limma on the resulting GSVA enrichment scores, we leave deliberately unchanged the default argument mx.diff=TRUE to obtain approximately normally distributed ES.

#cache(leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,min.sz=10, max.sz=500, verbose=TRUE),dir=cacheDir, prefix=cachePrefix)


# does it works with more than two classes?
leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,min.sz=10, max.sz=500, verbose=TRUE)

# the object is another expression set. but now instead of the genes there are the signatures as rownames
head(leukemia_es@assayData$exprs)

# We test whether there is a difference between the GSVA enrichment scores from each pair of phenotypes using a simple linear model and moderated t-statistics computed by the limma package using an empirical Bayes shrinkage method (see Smyth, 2004). We are going to examine both, changes at gene level and changes at pathway level and since, as we shall see below, there are plenty of them, we are going to employ the following stringent cut-offs to attain a high level of statistical and biological significance:
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)

# where we will use the latter only for the gene-level differential expression analysis.

design <- model.matrix(~ factor(leukemia_es$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgeneSets <- topTable(fit, coef="MLLvsALL", number=Inf,p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)
summary(res)

# vulcono plot for enriched pathways
allGeneSets%>%
  mutate(term=rownames(.))%>%
  # label the one that are below 0.05
  mutate(col=ifelse(adj.P.Val<0.001 ,yes = T,no = F))%>%
  ggplot(aes(x = logFC,y = -log10(adj.P.Val),col=col))+geom_point()+scale_color_manual("factor",values = c("grey","red"))

# make an heatmap of the most common disregulated pathways
# from the differential geneset collect the name of the differential sets

df_expression <- leukemia_es@assayData$exprs%>%
  data.frame()%>%
  mutate(term=rownames(.))%>%
  filter(term%in%rownames(DEgeneSets))%>%
  select(term,CL2001011101AA.CEL:CL2001011152AA.CEL )

# save the object ot recall the ordering of the sample
heatmap <- try(heatmap.2(mat),silent = T)

# from map retrieve the clustering of the sample
col_factor <- pData(leukemia_es)$subtype[heatmap$colInd]
# make it as color
unique(col_factor)
# produce a function for colors
color <- colorRampPalette(c("green","black","yellow"),space="rgb")
# priduce the colors
palette <- color(length(unique(col_factor)))
#palette <- color(10)
# plot the colors
#ggplot(data.frame(col=palette),aes(x=factor(1),fill=col))+geom_bar()+scale_fill_manual(values = palette)+xlab("")
# make the vector of factors as vector of colors
col_columns <- as.character(factor(col_factor,levels = unique(col_factor),labels = palette))

  
# set the color scale for the tiles in the heatmap
colfunc <- colorRampPalette(c("blue","white","red"))
mat <- as.matrix(df_expression[,-1])
rownames(mat) <- df_expression$term  

# make the heatmap with the colgroup
heatmap.2(
  # matrix of data
  mat,
  # color scheme
  col=colfunc(20),
  # other parameters
  density.info="none", trace="none", symm=F,symkey=T,symbreaks=T,
  # scale the data by row (gene expression)
  scale="row",
  # size of the font
  cexRow = 0.7, cexCol = 0.8,
  # rotation and position of the labels
  # srtCol=45, 
  #adjCol = c(1,1),
  # change the size of the margins
  margins = c(10,25),
  # add the marginal grouping
  ColSideColors = col_columns)


# Thus, there are 38 MSigDB C2 curated pathways that are differentially activated between MLL and ALL at 0.1% FDR. When we carry out the corresponding differential expression analysis at gene level:

logFCcutoff <- log2(2)
design <- model.matrix(~ factor(leukemia_eset$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_filtered_eset, design)
fit <- eBayes(fit)
allGenes <- topTable(fit, coef="MLLvsALL", number=Inf)
DEgenes <- topTable(fit, coef="MLLvsALL", number=Inf,p.value=adjPvalueCutoff, adjust="BH", lfc=logFCcutoff)
res <- decideTests(fit, p.value=adjPvalueCutoff, lfc=logFCcutoff)
summary(res)

# vulcano plot for enriched genes
# vulcono plot for enriched pathways
allGenes%>%
  mutate(term=rownames(.))%>%
  # label the one that are below 0.05
  mutate(col=ifelse(adj.P.Val<0.001 & abs(logFC)>1,yes = T,no = F))%>%
  ggplot(aes(x = logFC,y = -log10(adj.P.Val),col=col))+geom_point()+scale_color_manual("factor",values = c("grey","red"))+
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")

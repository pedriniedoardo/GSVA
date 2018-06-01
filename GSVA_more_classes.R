# GSVA using more than two classes

# load the required libraries
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

# We carry out a non-specific filtering step by discarding the 50% of the probesets with smaller variability, probesets without Entrez ID annotation, probesets whose associated Entrez ID is duplicated in the annotation, and Affymetrix quality control probes:

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

# 1
# extract the data for the expression
exp <- exprs(leukemia_filtered_eset)

# 2
# modify the pheno df in order to have more classes
df_pheno <- leukemia_filtered_eset@phenoData@data
# create a new col with more classes and overwrite it
newCol <- rep(c("A","B","C"),13)[1:37]
df_pheno$subtype <- newCol
df_pheno
# the object is a data frame
class(matrix_pheno)
# here is the code on how to built it from scratch. notice that the rownames whould be the same of the colnames of the expression set
# if wanted is possible to add metadata to the object
# phenoData <- new("AnnotatedDataFrame",data=df_pheno, varMetadata=metadata)
pheno <- new("AnnotatedDataFrame",data=df_pheno)
head(pheno)

# 3
# build experiment data 
expData <- new("MIAME",
               name="test",
               lab="test",
               contact="test",
               title="test",
               abstract="test",
               url="test",
               other=list(
                 notes="test"
               ))

class(expData)

# 4
# the annotation holds the info about the probe set, is just a character string
leukemia_filtered_eset@annotation

# 5
# build a new expression set from its component
exampleSet <- ExpressionSet(assayData = exp,
                            phenoData = pheno,
                            experimentData = expData,
                            annotation = "hgu95a")

exampleSet
####################################################################################

# does it works with more than two classes?
leukemia_es_test <- gsva(exampleSet, c2BroadSets,min.sz=10, max.sz=500, verbose=TRUE)

# it seems it worked anyway. interestingly the presence of the new phenotype do not change the estimation of the score of GSVA in the data.
identical(leukemia_es_test@assayData$exprs,leukemia_es@assayData$exprs)

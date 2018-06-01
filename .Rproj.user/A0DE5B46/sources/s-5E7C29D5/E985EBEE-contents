# GSVA using expression set built from the component

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
# exp is a matrix of normalized expression after filtering of constant genes
head(exp)

# 2
# pheno is a complex object build from a df with the organization of the samples (the classes of the sample)
df_pheno <- leukemia_filtered_eset@phenoData@data
# the object is a data frame
class(matrix_pheno)
# here is the code on how to built it from scratch. notice that the rownames whould be the same of the colnames of the expression set
# if wanted is possible to add metadata to the object
# phenoData <- new("AnnotatedDataFrame",data=df_pheno, varMetadata=metadata)
pheno <- new("AnnotatedDataFrame",data=df_pheno)
head(pheno)

# 3
# experiment data are just annotations
leukemia_filtered_eset@experimentData
# here is how to build the object from scratch
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

# now technically we should obtain the same reuslt by running the GSVA on both the new and old object

# does it works with more than two classes?
leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets,min.sz=10, max.sz=500, verbose=TRUE)
leukemia_es_test <- gsva(exampleSet, c2BroadSets,min.sz=10, max.sz=500, verbose=TRUE)

# quickly compare the two matricies of the set
identical(leukemia_es@assayData$exprs,leukemia_es_test@assayData$exprs)

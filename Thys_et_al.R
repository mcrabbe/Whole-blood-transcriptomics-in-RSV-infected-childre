
# Start from Bioconductor ExpressionSet "eset"

#===============================================================================
# Batch effects correction with ComBat
#===============================================================================

library(sva)

com <- ComBat(exprs(eset), batch=pData(eset)$BATCH, mod=NULL, par.prior = TRUE)

#===============================================================================
# Limma analysis
#===============================================================================

library(limma)

esetLimma <- eset[,!is.na(eset$GENDER) & !is.na(eset$AGEMONTHS)]

gps <- paste(pData(esetLimma)[,"GROUP"],pData(esetLimma)[,"SEVERITY"],pData(esetLimma)[,"AGECAT2"],sep=".")

myDesign <- model.matrix(~ 0 + factor(gps) + GENDER + AGEMONTHS, data = pData(esetLimma))
colnames(myDesign) <- sub("factor\\(gps\\)", "", colnames(myDesign))

corfit <- duplicateCorrelation(esetLimma, myDesign, block = pData(esetLimma)$SUBJECTID)
fit    <- lmFit(esetLimma, myDesign, block = pData(esetLimma)$SUBJECTID, correlation = corfit$consensus)

contMatrix <- makeContrasts(
  
  D1.vs.HC    = (DISDAY1.MILD.HIGH + DISDAY1.MILD.LOW + DISDAY1.MOD.HIGH + DISDAY1.MOD.LOW + DISDAY1.SEVERE.HIGH + DISDAY1.SEVERE.LOW)/6 - 
    (HC.HC.HIGH + HC.HC.LOW)/2,
  D30.vs.HC   = (DISDAY30.MILD.HIGH + DISDAY30.MILD.LOW + DISDAY30.MOD.HIGH + DISDAY30.MOD.LOW + DISDAY30.SEVERE.HIGH + DISDAY30.SEVERE.LOW)/6 - 
    (HC.HC.HIGH + HC.HC.LOW)/2,
  D180.vs.HC  = (DISDAY180.MILD.HIGH + DISDAY180.MILD.LOW + DISDAY180.MOD.HIGH + DISDAY180.MOD.LOW + DISDAY180.SEVERE.HIGH + DISDAY180.SEVERE.LOW)/6 - 
    (HC.HC.HIGH + HC.HC.LOW)/2,
  D30.vs.D1   = (DISDAY30.MILD.HIGH + DISDAY30.MILD.LOW + DISDAY30.MOD.HIGH + DISDAY30.MOD.LOW + DISDAY30.SEVERE.HIGH + DISDAY30.SEVERE.LOW)/6 - 
    (DISDAY1.MILD.HIGH + DISDAY1.MILD.LOW + DISDAY1.MOD.HIGH + DISDAY1.MOD.LOW + DISDAY1.SEVERE.HIGH + DISDAY1.SEVERE.LOW)/6,
  D180.vs.D1  = (DISDAY180.MILD.HIGH + DISDAY180.MILD.LOW + DISDAY180.MOD.HIGH + DISDAY180.MOD.LOW + DISDAY180.SEVERE.HIGH + DISDAY180.SEVERE.LOW)/6 - 
    (DISDAY1.MILD.HIGH + DISDAY1.MILD.LOW + DISDAY1.MOD.HIGH + DISDAY1.MOD.LOW + DISDAY1.SEVERE.HIGH + DISDAY1.SEVERE.LOW)/6,
  
  levels = myDesign)

lmFit1 <- eBayes(contrasts.fit(fit, contMatrix))

#===============================================================================
# WGCNA
#===============================================================================

library(WGCNA)

datExpr <- t(exprs(eset))

softPower <- 6
minModuleSize <- 20

adjacency <- adjacency(datExpr, power = softPower, type = "unsigned")

TOM <- TOMsimilarity(adjacency, TOMType = "unsigned", verbose=5)
dissTOM <- 1-TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

net_unsigned <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              method = "hybrid", pamStage = TRUE,
                              deepSplit = 2, maxDistToLabel = 0,
                              minClusterSize = minModuleSize)

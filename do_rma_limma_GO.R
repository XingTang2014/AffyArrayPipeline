library(affy)
library(limma)
library("plyr")
library("GO.db")
library(GOstats)
options(stringsAsFactors=FALSE)

parameters <- read.table("parameters.txt", blank.lines.skip=TRUE, header=FALSE, row.names =1, comment.char="#")
arrayType <- parameters["array",1]
specie <- parameters["specie",1]
DEG_cutoff <- parameters["DEG_cutoff",1]
GO_cutoff <- parameters["GO_cutoff",1]
library(arrayType, character.only = TRUE)
if(specie=="mouse"){ 
  library(org.Mm.eg.db) 
  goDB <- "org.Mm.eg"}
if(specie=="human"){ 
  library(org.Hs.eg.db) 
  goDB <- "org.Hs.eg"}

# prepare pheno data
sample_desc <- read.table("sample_desc.txt", sep="\t", header = TRUE)
phenoData <- AnnotatedDataFrame(sample_desc)

# rma call expression
eset <- rma(read.affybatch(paste("CELs/", sample_desc[,"cel_file"], sep="")))
exps <- exprs(eset)  
colnames(exps) <- as.character(sample_desc[, "sample_name", drop=TRUE])
Symbols <-  unlist(mget(rownames(exps), eval(parse(text = sub(".db$", "SYMBOL",arrayType))), ifnotfound=NA))
exp_table <- data.frame(Probes=rownames(exps), Symbols=Symbols, exps)
write.table(exp_table, file="results/AllExp.xls", sep="\t", quote=FALSE, row.names=FALSE)

# limma call DEG
combn <- factor(paste(pData(phenoData)[,"condition"],
                      pData(phenoData)[,"tissue"], sep = "_"))
design <- model.matrix(~combn) # describe model to be fit
fit <- lmFit(eset, design)  # fit each probeset to model
efit <- eBayes(fit)        # empirical Bayes adjustment
DEG_table <- topTable(efit, coef=2, number=nrow(exp_table))      # table of differentially expressed probesets
Symbols <-  unlist(mget(rownames(exps), eval(parse(text = sub(".db$", "SYMBOL",arrayType))), ifnotfound=NA))
DEG_table <- data.frame(Probes=rownames(DEG_table), Symbols=Symbols, DEG_table)
write.table(DEG_table, file="results/DEG.xls", sep="\t", quote=FALSE, row.names=FALSE)

# GO analysis
annoDB <- paste(goDB, ".db", sep="")
SYMBOL2EG <- eval(parse(text = paste(goDB, "SYMBOL2EG", sep="")))
GO2ALLEGS <- eval(parse(text = paste(goDB, "GO2ALLEGS", sep="")))
SYMBOL <- eval(parse(text = paste(goDB, "SYMBOL", sep="")))
all.genes <- intersect(unique(exp_table[,"Symbols"]), mappedkeys( SYMBOL2EG ))
univ <- unique(unlist(mget(all.genes, SYMBOL2EG  )))
DEGs <- unique(DEG_table[ DEG_table[,"adj.P.Val"] < 0.01, "Symbols"])
gene_list <- intersect(DEGs, all.genes)
entrezList <- unique(unlist(mget(gene_list, SYMBOL2EG)))

ParamObjs <- new("GOHyperGParams", geneIds = entrezList, universeGeneIds = univ, annotation = annoDB, ontology = "BP")
hyper_res <- summary( hyperGTest(ParamObjs), pvalue=1)
adjpValueList <- p.adjust( hyper_res[,2], method="BH" )
output <- cbind( hyper_res, adjpValueList )[ adjpValueList < GO_cutoff, ]
if( nrow(output) > 0 ) {
  colnames( output ) <- c( colnames( hyper_res), "BH-adj.Pvalue" )
  goMaps <- lapply(output[["GOBPID"]], function(x) unlist(mget(x, GO2ALLEGS)))
  goSelected <- lapply(goMaps, function(x) { temp <- entrezList[entrezList %in% x]; unlist(mget(temp, SYMBOL))} )
  output <- data.frame( output, Genes=unlist(lapply(goSelected, function(x) paste(x, collapse=";"))) )
  write.table( output, file="results/GO.xls", sep="\t", quote=FALSE , row.names=FALSE)
}

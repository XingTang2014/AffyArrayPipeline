###---------------------------------------------------------------------------------------------------------------------
### Affy + limma + GOstat	    pipeline for analysing affymetrix expression array data.
### 06/25/2015            	    Xing Tang, The Ohio State University Comprehensive Cancer Center
###				    tangx1986@gmail.com
###---------------------------------------------------------------------------------------------------------------------
How to use this pipeline?
1. customize you parameters by modifying sample_desc.txt file
2. type "Rscript do_rma_limma_GO.R" to get the results
3. check results/ folder. Four tables are generated.
   a. AllExp.xls: log2 signal generated by rma function from affy
   b. DEG.xls: differential expression analysis results from limma
   c. GO.xls: enriched GO terms from GOstat package
 
*Note
1. you have to modify the do_rma_limma_GO.R scripts to fit the pipeline to multiple groups comparison

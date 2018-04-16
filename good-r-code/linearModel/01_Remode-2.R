#!usr/bin/Rscript
########### load in the library
library(doParallel)
library(doMC)
registerDoMC()
library(foreach)
library(MASS)
library(dplyr)
library(methods)
###########  load in the data for each of cancer types to list object
expressionPath = '/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/'
expressionFileAll <- list.files(path = expressionPath, pattern = "_mRNA_each_exp_20160513")
APAFilepath = '/home/yxiang/data/Project/refine_APA_coverage_30off/01_refined_data_using_30_coverage//TCGA/refined_data/'
APAFileAll <- list.files(path =APAFilepath,   pattern = "_result.txt")
APAfactorID = read.table('/tmp/APA_geneID.txt')
APAfactorID = APAfactorID$V1 %>% as.character()
getExpressionID = function(x){
    strsplit(x, split = '|', fixed = T)[[1]][2]
}

for (curr_mRNAFile in expressionFileAll) {
    Cond <- gsub("_mRNA_each_exp_20160513", "", curr_mRNAFile)
    APA_file <- APAFileAll[ grepl(paste0('_',Cond,"_30cutoffRefined_result.txt"), APAFileAll )]
    if(length(APA_file)==0) next
    Matrix_tumor_t <- read.table(file.path(APAFilepath,APA_file), row.names = 1, sep = "\t", comment.char = "", quote = "", fill= TRUE, header = T)
    Seed_tumor_t <- read.table(file.path(expressionPath,curr_mRNAFile), row.names = 1, sep = "\t", header = T )
    rowID = lapply(rownames(Seed_tumor_t), getExpressionID) %>% unlist()    ## filter our the 22 apa factors
    Seed_tumor_t = Seed_tumor_t[ rowID %in%APAfactorID ,]    ## filter our the 22 apa factors
    Matrix_tumor_t <- Matrix_tumor_t[,grepl("PDUI", colnames(Matrix_tumor_t))]
    pathToAnnotation <- system(paste0("echo /extraspace/TCGA/Depth/TCGA/",Cond, "/annotations/*.txt"), intern = T)
    Annotation <- read.table(pathToAnnotation, stringsAsFactors = F)
    colnames(Matrix_tumor_t) <- as.vector(Annotation[,1])
    colnames(Matrix_tumor_t) <- gsub("\\-", ".", colnames(Matrix_tumor_t))
    print (paste(curr_mRNAFile, APA_file, "analysis \n"))
    ############  pull out the matched tumor data
    ############  for the APA data
    keeplink <- grep("Tumor", colnames(Matrix_tumor_t))
    Matrix_tumor_t <- Matrix_tumor_t[, keeplink]
    colnames(Matrix_tumor_t) <- gsub(paste0(Cond,".Tumor."), "", colnames(Matrix_tumor_t))
    ############  for the mRNA data
    temp <- substr(colnames(Seed_tumor_t), 14, 15)
    table(temp)
    if( Cond == 'SKCM'){ # 
        Seed_tumor_t <- Seed_tumor_t[, temp == "06"]
    }else{
        Seed_tumor_t <- Seed_tumor_t[, temp == "01"]
    }
    temp <- substr(colnames(Seed_tumor_t), 1, 12)
    Seed_tumor_t <- Seed_tumor_t[, !duplicated(temp)]
    #if(!(length(temp) == length(unique(temp))))  stop("Duplicated mRNA ID", call. = F)
    colnames(Seed_tumor_t) <-  substr(colnames(Seed_tumor_t), 1, 12)
    ############  As both data set now has unique patient for each ID, we do not remove the duplicates
    ############  otherwise we will need to remove the duplicates for each patient
    
    ############  we find common samples
    commonsamples <- intersect(colnames(Seed_tumor_t), colnames(Matrix_tumor_t))
    Matrix_tumor_t <- Matrix_tumor_t[, commonsamples]
    Seed_tumor_t <- Seed_tumor_t[,commonsamples]
    if (!(identical(colnames(Seed_tumor_t), colnames(Matrix_tumor_t)))) {
        stop("mRNA and APA not identical", call. = F)
    }
    ############  transform
    Matrix_tumor_t <- t(Matrix_tumor_t)
    Seed_tumor_t <- t(Seed_tumor_t)
    ######################## clean mRNA name
    genenameused <- unlist(lapply(strsplit(colnames(Seed_tumor_t), "\\|"), function(x){x[[1]]}))
    colnames(Seed_tumor_t) <- genenameused
    ########################  scale the data, here we suppose nothing need log2 transformation
    Matrix_tumor_t_scaled <- scale(Matrix_tumor_t)
    Seed_tumor_t_scaled <- scale(Seed_tumor_t)
    set.seed(123)
    
    ########################  model selection
    output_tumor <- foreach(i = 1:ncol(Matrix_tumor_t_scaled),.errorhandling='pass',.verbose=F) %dopar% {
        effectN <- sum(!is.na(Matrix_tumor_t_scaled[,i]))	# sample size effective
        keeplink <- which(!is.na(Matrix_tumor_t_scaled[,i]))
        lm_df <- data.frame(editting = scale(Matrix_tumor_t[keeplink,i]), scale(Seed_tumor_t[keeplink,]), check.names= F)
        result <- lm( editting ~ ., data = lm_df)
        #  model selection
        result_update <- stepAIC(result,~ ., trace=F)
        result_smmy <- summary(result_update)
        if (is.null(result_smmy$fstatistic)) {
            model_pvalue <- NA
        } else model_pvalue <- with(as.data.frame(t(result_smmy$fstatistic)),pf(value,numdf,dendf,lower.tail=F))
        # calculate the partial correlation
        PartialCorr <- coef(result_smmy)[-1,'Estimate']
        Pr_coefs <- coef(result_smmy)[-1,'Pr(>|t|)']	# remove intercept estimate
        # extract a neat name essence for usage later
        nameEssence <- rownames(coef(result_smmy))[-1]
        # neat the name
        names(Pr_coefs) <- sprintf("Pvalue_%s",nameEssence)
        names(PartialCorr) <- sprintf("PartialCorr_%s",nameEssence)
        adjRsq <- result_smmy$adj.r.squared
        # then the spearman corr with no adjustments on other covariates
        spearCors <- apply(Seed_tumor_t[, Seed_tumor_t_scaled %>% colnames],2,cor.test,y=Matrix_tumor_t[,i],method="spearman",na.action=na.omit)	
        # only the subset as the same as in Seed_tumor_t_scaled
        spearCors_vec <- lapply(spearCors,function(x)  with(x,c(spearCorEstimate = estimate,spearCorP = p.value)))
        spearCors_vec <- unlist(spearCors_vec)
        return(c(effectN = effectN,model_pvalue = model_pvalue,PartialCorr,Pr_coefs,adjRsq = adjRsq,spearCors_vec ) ) 
    } # modelSelection end
    names(output_tumor) <- colnames(Matrix_tumor_t)
    my.fig <- "ModelSelection"
    if (!file.exists(my.fig))   dir.create(my.fig) 
    save(output_tumor,file = file.path(my.fig,sprintf("%s_regressOn_%s.Rdata",Cond,"AllSeedGene")))	# that is for the original
    
    ############## remove those with models not run or with less than 50 samples have APA 
    output_tumor <- output_tumor[unlist(sapply(output_tumor,function(x) !any( grepl("simpleError",class(x)) ) ) )]	# delete those are simpleError returned by foreach.
    output_tumor <- output_tumor[ unlist(sapply(output_tumor,function(x) x['effectN'] >= 100 ) ) ] # effective N!!
    if (length(output_tumor)==0) next
    pvalues <- sapply(output_tumor,function(x) x['model_pvalue'] )
    qvalues <- p.adjust(pvalues,method='fdr')
    tempused <- c('effectN','model_pvalue', paste('PartialCorr_', genenameused,  sep = ""), 
                  paste('Pvalue_', genenameused,  sep = ""), 'adjRsq',
                  paste(genenameused, '.spearCorEstimate.rho', sep = ""),  
                  paste(genenameused, '.spearCorP', sep = ""))
    
    ############  in some case, the name as PartialCorr_'XXX', so need the code below
    ##	tempused <- c('effectN','model_pvalue', paste('PartialCorr_`', genenameused, '`', sep = ""), 
    #						paste('Pvalue_`', genenameused, '`', sep = ""), 'adjRsq',
    #						paste(genenameused, '.spearCorEstimate.rho', sep = ""),  
    #						paste(genenameused, '.spearCorP', sep = ""))
    
    ############ this will also affect low other code to define the name
    ############ each time, will check the name first before the run
    
    output_unify <- function(x){ x[tempused] }
    output_tumor_refined <- sapply(output_tumor,output_unify)
    output_tumor_refined_df <- as.data.frame(t( output_tumor_refined ) )
    colnames(output_tumor_refined_df) <- tempused
    output_tumor_refined_df$model_qvalue <- qvalues
    #
    Pvalues_enzymes <- as.relistable( output_tumor_refined_df[paste('Pvalue_', genenameused,  sep = "")] )
    ############### multiple test adjustment for the partial correlation p value 
    qvalues_enzymes <- p.adjust( unlist(Pvalues_enzymes),method='fdr' )
    qvalues_enzymes_df <- matrix(qvalues_enzymes, ncol = length(genenameused), byrow=F )
    colnames(qvalues_enzymes_df) <- sub("Pvalue","qvalue",paste('Pvalue_', genenameused,  sep = ""),perl=T)
    output_tumor_refined_df2 <- data.frame(output_tumor_refined_df,qvalues_enzymes_df ,check.names=F)
    my.fig2 <- "OutputClean"
    if (!file.exists(my.fig2))   dir.create(my.fig2) 
    write.table(output_tumor_refined_df2,file=file.path(my.fig2, paste0(Cond,"_Tumor_regressOn_AllSeedGene.txt")),sep='\t',row.names=T,col.names=NA,quote=F)
    Result <- output_tumor_refined_df2
    save(Result, file = file.path(my.fig2, paste0(Cond, "_Result.RData")))
    
}

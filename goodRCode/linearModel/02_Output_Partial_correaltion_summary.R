# library(doParallel)
# library(doMC)
# registerDoMC()

dataFiles <- list.files(path = "OutputClean",
                        pattern = "AllSeedGene.txt")
temp <- read.csv("APA_factor_geneID.csv")
genenameused <- na.omit(as.vector(temp[,"symbol"]))
genenameused <- genenameused[which(genenameused != "")]
genes <- toupper(gsub(" ", "", genenameused))
output <- data.frame(factor = genes)
outputPercentage <- data.frame(factor = genes)
# foreach ( currFile = dataFiles) %dopar%{

if(!file.exists("reg_summary")) dir.create("reg_summary")

for (currFile in dataFiles){
    cancerType <- gsub("_Tumor_regressOn_AllSeedGene.txt", "", currFile)
    rawData <- read.table(file.path('OutputClean',currFile), header = T, quote = "")
    data <- rawData[ rawData$model_qvalue < 0.05 &  !is.na(rawData$model_qvalue),]
    geneNum <- vector()
    geneTotal <- vector()
    for (gene in genes){
        Corr <- (data[ ,paste0("PartialCorr_", gene)] > 0.3 )
        qvalue <- (data[,paste0("qvalue_", gene)] < 0.05)
        signGene <- sum( Corr & qvalue, na.rm = T )
        geneSum <- sum(!is.na(data[ ,paste0("PartialCorr_", gene)]) & !is.na(data[,paste0("qvalue_", gene)]))
        geneNum <- c(geneNum, signGene)
        geneTotal <- c(geneTotal,geneSum )
    }
    genePer <- geneNum/geneTotal *100
    output <- cbind(output, geneNum)
    outputPercentage <- cbind(outputPercentage,genePer )
    colnames(output)[ncol(output)] <- cancerType
    colnames(outputPercentage)[ncol(outputPercentage)] <- cancerType
}
write.table(output, "reg_summary/APAfactor_regresson_Pos.xls", sep = "\t", row.names = F, quote = F)
write.table(outputPercentage, "reg_summary/APAfactor_regresson_percentage_Pos.xls", sep = "\t", row.names = F, quote = F)

output <- data.frame(factor = genes)
outputPercentage <- data.frame(factor = genes)
for (currFile in dataFiles){
    cancerType <- gsub("_Tumor_regressOn_AllSeedGene.txt", "", currFile)
    rawData <- read.table(file.path('OutputClean',currFile), header = T, quote = "")
    data <- rawData[ rawData$model_qvalue < 0.05 &  !is.na(rawData$model_qvalue),]
    geneNum <- vector()
    geneTotal <- vector()
    for (gene in genes){
        Corr <- (data[ ,paste0("PartialCorr_", gene)] < ( -0.3) )
        qvalue <- (data[,paste0("qvalue_", gene)] < 0.05)
        signGene <- sum( Corr & qvalue, na.rm = T )
        geneSum <- sum(!is.na(data[ ,paste0("PartialCorr_", gene)]) & !is.na(data[,paste0("qvalue_", gene)]))
        geneNum <- c(geneNum, signGene)
        geneTotal <- c(geneTotal,geneSum )
    }
    genePer <- geneNum/geneTotal *100
    output <- cbind(output, geneNum)
    outputPercentage <- cbind(outputPercentage,genePer )
    colnames(output)[ncol(output)] <- cancerType
    colnames(outputPercentage)[ncol(outputPercentage)] <- cancerType
}

write.table(output, "reg_summary/APAfactor_regresson_Neg.xls", sep = "\t", row.names = F, quote = F)
write.table(outputPercentage, "reg_summary/APAfactor_regresson_percentage_Neg.xls", sep = "\t", row.names = F, quote = F)

output <- data.frame(factor = genes)
outputPercentage <- data.frame(factor = genes)
for (currFile in dataFiles){
    cancerType <- gsub("_Tumor_regressOn_AllSeedGene.txt", "", currFile)
    rawData <- read.table(file.path('OutputClean',currFile), header = T, quote = "")
    data <- rawData[ rawData$model_qvalue < 0.05 &  !is.na(rawData$model_qvalue),]
    geneNum <- vector()
    geneTotal <- vector()
    for (gene in genes){
        Corr <- (abs(data[ ,paste0("PartialCorr_", gene)]) >  0.3 )
        qvalue <- (data[,paste0("qvalue_", gene)] < 0.05)
        signGene <- sum( Corr & qvalue, na.rm = T )
        geneSum <- sum(!is.na(data[ ,paste0("PartialCorr_", gene)]) & !is.na(data[,paste0("qvalue_", gene)]))
        geneNum <- c(geneNum, signGene)
        geneTotal <- c(geneTotal,geneSum )
    }
    genePer <- geneNum/geneTotal *100
    output <- cbind(output, geneNum)
    outputPercentage <- cbind(outputPercentage,genePer )
    colnames(output)[ncol(output)] <- cancerType
    colnames(outputPercentage)[ncol(outputPercentage)] <- cancerType
}

write.table(output, "reg_summary/APAfactor_regresson_Total.xls", sep = "\t", row.names = F, quote = F)
write.table(outputPercentage, "reg_summary/APAfactor_regresson_percentage_Total.xls", sep = "\t", row.names = F, quote = F)



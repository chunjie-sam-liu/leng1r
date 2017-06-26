###Quantiling Circadian genes expression, and keep upper(2) and bottom(1) quantile patient samples
setwd("/extraspace/yye1/analysis/Circadian/RPPA")
filenameall <- list.files("/extraspace/TCGA/TCGA_protein/", pattern="*re_20160627")
tumornameall <-  gsub("_protein.re_20160627","",filenameall)
circadian.genes <- read.table("~/Circadian/circadian.genes.txt",header=F) 
expfileall <- list.files("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="*_mRNA_each_exp_20160513")
expfileallname <- gsub("_mRNA_each_exp_20160513","",expfileall)
tumornameall <- intersect(expfileallname,tumornameall)
#previous <- list.files("/extraspace/TCGA/TCGA_RPPA/",pattern="*csv")
#previousname <- toupper(gsub("rppa_|\\.csv","",previous))
#tumornameall <- setdiff(tumornameall,previousname)
if (!file.exists("Circadian.gene.cluster_median")) {
  dir.create("Circadian.gene.cluster_median")
}
for(m in 1:length(tumornameall)){
  if (!file.exists(file.path("Circadian.gene.cluster_median", tumornameall[m]))) {
    dir.create(file.path("Circadian.gene.cluster_median", tumornameall[m]))
  }
  sub <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",toupper(tumornameall[m]),"_mRNA_each_exp_20160513",sep=""),header=T)
  sub$gene <- data.frame(do.call(rbind,strsplit(as.character(sub$gene),"\\|")))$X1
  sub <- sub[which(sub$gene %in% circadian.genes$V1),] ##get circadian genes
  ###consider tumor type contain Tumor and normal genes
  sub <- sub[,c("gene",colnames(sub)[as.numeric(substr(colnames(sub),14,15)) %in% c(1,6)])]
  colnames(sub) <- gsub("\\.","\\-",substr(colnames(sub),1,12))
  for(n in sub$gene){
    sub1 <- sub[which(sub$gene==n),2:ncol(sub)]
    sub1 <- data.frame(t(sub1))
    colnames(sub1) <- n
    sub1$Sample <-  rownames(sub1)
    #sub1 <- sub1[which(sub1[n] <= quantile(sub1[,n])[2] | sub1[n] >= quantile(sub1[,n])[4]),]
    sub1 <- sub1[order(sub1[n]),]
    sub1$cluster <- rep(1,times=nrow(sub1))
    sub1["cluster"][sub1[n] >= quantile(sub1[,n])[3]] <- 2
    write.csv(sub1,file.path("Circadian.gene.cluster_median/", tumornameall[m], paste(n, "-bottom_upper_quantile_class.csv", sep = "")),row.names = F)
    
  }
}

library(gdata)
library(survival)
library(plotrix)  
filenameall <- list.files("/extraspace/TCGA/TCGA_protein/", pattern="*_protein.re_20160627")
filenameall <- filenameall[gsub("_protein.re_20160627","",filenameall) %in% tumornameall]
Output <- list()
usedid <- 1
for(kkk1 in 1:length(filenameall))
{
  
  
  filename <- filenameall[kkk1]
  
  Cond <- gsub("_protein.re_20160627","",filename)
  Fulldata <- read.delim(file.path("/extraspace/TCGA/TCGA_protein/", filename))
  if(Cond == "SKCM")
  {
    keeplink <- which(as.numeric(substr(as.character(colnames(Fulldata)),14,15)) == 6)
    Fulldata <- Fulldata[,c(1,keeplink)]
  }else
  {
    keeplink <- which(as.numeric(substr(as.character(colnames(Fulldata)),14,15)) == 1)
    Fulldata <- Fulldata[,c(1,keeplink)]
    
    
  }
  rownames(Fulldata) <- Fulldata$protein
  
  FullRPPAData <- Fulldata[,2:ncol(Fulldata)]
  colnames(FullRPPAData) <- gsub("\\.","\\-",substr(as.character(colnames(FullRPPAData)),1,12))
  ##############   calculate signature
  
  ScoreUsedData <- sweep(FullRPPAData, 1, apply(FullRPPAData, 1 , median), FUN="-")
  ScoreUsedData <- sweep(ScoreUsedData, 1, apply(FullRPPAData, 1 , sd), FUN="/")
  
  
  usedlink <- grep("Akt_", rownames(ScoreUsedData))
  AKT <- cor(t(ScoreUsedData[usedlink,]))[2,1]
  usedlink <- c(grep("GSK3.alpha.beta_", rownames(ScoreUsedData)),
                grep("GSK3_pS9", rownames(ScoreUsedData)))
  GSK <- cor(t(ScoreUsedData[usedlink,]))[2,1]
  
  usedlink <- grep("p27_", rownames(ScoreUsedData))
  p27 <- cor(t(ScoreUsedData[usedlink,]))[2,1]
  
  
  usedlink <- grep("EGFR_", rownames(ScoreUsedData))
  pEGFR <- cor(t(ScoreUsedData[usedlink,]))[2,1]
  
  usedlink <- grep("S6_", rownames(ScoreUsedData))
  pS6 <- cor(t(ScoreUsedData[usedlink,]))[2,1]
  
  usedlink <- grep("Src_", rownames(ScoreUsedData))
  pSrc <- cor(t(ScoreUsedData[usedlink,]))[2,1]
  
  usedlink <- grep("ER.alpha", rownames(ScoreUsedData))
  ERalpha <- cor(t(ScoreUsedData[usedlink,]))[2,1]
  
  
  usedlink <- grep("4E.BP1_p", rownames(ScoreUsedData))
  tt <- cor(t(ScoreUsedData[usedlink,]))
  EBP1.65.37 <- tt[grep("65", rownames(tt)), grep("37", colnames(tt))]
  EBP1.65.70 <- tt[grep("65", rownames(tt)), grep("70", colnames(tt))]
  EBP1.37.70 <- tt[grep("37", rownames(tt)), grep("70", colnames(tt))]
  
  score <- c(AKT, GSK, p27, pEGFR, pS6, EBP1.65.37, EBP1.65.70, EBP1.37.70, pSrc, ERalpha)
  names(score) <- c("AKT", "GSK", "p27", "pEGFR", "pS6", "EBP1.65.37", "EBP1.65.70", "EBP1.37.70", "pSrc", "ERalpha")
  Output[[usedid]] <- score
  names(Output)[usedid] <- Cond
  usedid <- usedid +1
  
}

pdf(file.path("/extraspace/yye1/analysis/Circadian/RPPA/", paste("Output-RPPA-Cor-",  ".pdf", sep = "")))
library(plyr)
for(i1 in 1:length(Output[[1]]))
{
  Cond <- names(Output[[1]])[i1]
  
  useddata <- ldply(Output, function(x){x[i1]})
  plot(c(1:dim(useddata)[1]), useddata[,2], xaxt = "n", main = colnames(useddata)[2], pch =16, xlab= "", ylab = "")
  axis(1, at = c(1:dim(useddata)[1]), names(Output), las = 2)
  abline(h = 0.85, lty =2)
  
}
dev.off()


write.table(c(PI3KAKTSname,RASMAPKSname,RTKname,TSCmTORname, Apoptosisnames,CellCyclename,DNADamageResponsename,EMTname,Hormoneaname,Hormonebname),
            file="/extraspace/yye1/analysis/Circadian/RPPA/pathway_protein.txt",quote = F,row.names = F,col.names = F)
pathway_protein <- c(PI3KAKTSname,RASMAPKSname,RTKname,TSCmTORname, Apoptosisnames,CellCyclename,DNADamageResponsename,EMTname,Hormoneaname,Hormonebname)
####################################
##################  read in RPPA data
library(gdata)
library(survival)
library(plotrix)  
filenameall <- list.files("/extraspace/TCGA/TCGA_protein//", pattern="*re_20160627")
filenameall <- filenameall[gsub("_protein.re_20160627","",filenameall) %in% tumornameall]
#library(foreach)

for(kkk1 in 1:length(filenameall))
{
  
  filename <- filenameall[kkk1]
  
  Cond <- gsub("_protein.re_20160627","",filename)
  Fulldata <- read.delim(file.path("/extraspace/TCGA/TCGA_protein/", filename))
  if(Cond == "SKCM")
  {
    keeplink <- which(as.numeric(substr(as.character(colnames(Fulldata)),14,15)) == 6)
    Fulldata <- Fulldata[,c(1,keeplink)]
  }else
  {
    keeplink <- which(as.numeric(substr(as.character(colnames(Fulldata)),14,15)) == 1)
    Fulldata <- Fulldata[,c(1,keeplink)]
    
    
  }
  rownames(Fulldata) <- Fulldata$protein
  
  FullRPPAData <- Fulldata[,2:ncol(Fulldata)]
  colnames(FullRPPAData) <- gsub("\\.","\\-",substr(as.character(colnames(FullRPPAData)),1,12))
  FullRPPAData <- FullRPPAData[,!duplicated(colnames(FullRPPAData))]
  outputp <- NA
  countid <- 2
  
  
  ##############   calculate signature
  FullRPPAData$protein <- gsub("-R-C|-R-V|-M-C|-M-V|-R-E|-G-C|-R-E|-M-E","",rownames(FullRPPAData))
  FullRPPAData[is.na(FullRPPAData) == T] <- 0
  FullRPPAData.m <- sapply(split(FullRPPAData[,1:(ncol(FullRPPAData)-1)],FullRPPAData$protein),colMeans)
  FullRPPAData.m <- t(FullRPPAData.m)
  # ScoreUsedData <- sweep(FullRPPAData, 1, apply(FullRPPAData, 1 , median), FUN="-")
  #ScoreUsedData <- sweep(ScoreUsedData, 1, apply(FullRPPAData, 1 , sd), FUN="/")
 # ScoreUsedData <- FullRPPAData[which(rownames(FullRPPAData) %in% pathway_protein),]
  ScoreUsedData <- FullRPPAData.m[apply(FullRPPAData.m, 1, function(y) !all(is.na(y))),apply(FullRPPAData.m, 2, function(y) !all(is.na(y)))]
  
  
  ##############  check whether we need to average between replicates
  meanout <- function(usedlink)
  {
    
    tt <- cor(t(ScoreUsedData[usedlink,]))[2,1]
    
    if(tt > 0.85)
    {
      tt <- apply(ScoreUsedData[usedlink,],2, mean)
    }else
    {
      tt <- apply(ScoreUsedData[usedlink,],2, sum)
    }
    
    return (tt)
  }
  
  tmp <- grep("Akt_", rownames(ScoreUsedData))
  AKT <- meanout(tmp)
  tmp <- c(grep("GSK3.alpha.beta_", rownames(ScoreUsedData)),
           grep("GSK3_pS9", rownames(ScoreUsedData)))
  GSK <- meanout(tmp)
  tmp <- grep("EGFR_", rownames(ScoreUsedData))
  pEGFR <- meanout(tmp)
  tmp <- grep("S6_", rownames(ScoreUsedData))
  pS6 <- meanout(tmp)
  
  
  
  ##########  PI3KAKTScore
  linkused <- c(grep("p27_", rownames(ScoreUsedData)),
                grep("PRAS40_", rownames(ScoreUsedData)),
                grep("Tuberin_pT", rownames(ScoreUsedData)))
  
  
  PI3KAKTSname <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  
  linkused <- c(grep("INPP4B", rownames(ScoreUsedData)),
                grep("PTEN", rownames(ScoreUsedData)))
  PI3KAKTSname <- c(PI3KAKTSname,rownames(ScoreUsedData)[linkused])
  NegativeSingle <- ScoreUsedData[linkused , ]
  NegativeSingle[is.na(NegativeSingle)] <- 0
  
  PI3KAKTScore <- apply(PositiveSingle, 2, sum) - apply(NegativeSingle, 2, sum) +AKT + GSK 
  
  
  
  
  ##########  RASMAPKScore
  linkused <- c(grep("A.Raf_", rownames(ScoreUsedData)),
                grep("c\\.Jun_pS73", rownames(ScoreUsedData)),
                grep("C.Raf_", rownames(ScoreUsedData)),
                grep("JNK_", rownames(ScoreUsedData)),
                grep("MAPK_", rownames(ScoreUsedData)),
                grep("MEK1_", rownames(ScoreUsedData)),
                grep("p38_p", rownames(ScoreUsedData)), 
                grep("p90RSK_p", rownames(ScoreUsedData)), 
                grep("Shc_p", rownames(ScoreUsedData)),
                grep("YB-1_p", rownames(ScoreUsedData)))
  
  RASMAPKSname <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  
  RASMAPKScore <- apply(PositiveSingle, 2, sum)
  
  
  
  ##########  RTKScore
  linkused <- c(grep("HER2_", rownames(ScoreUsedData)),
                grep("HER3_", rownames(ScoreUsedData)),
                grep("Ret_", rownames(ScoreUsedData)),
                grep("Shc_p", rownames(ScoreUsedData)), 
                grep("Src_p", rownames(ScoreUsedData)))
  
  RTKname <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  
  RTKScore <- apply(PositiveSingle, 2, sum) + pEGFR
  
  
  
  ##########  TSCmTORScore
  linkused <- c(grep("4E.BP1_p", rownames(ScoreUsedData)),
                grep("mTOR_", rownames(ScoreUsedData)),
                grep("p70S6K_", rownames(ScoreUsedData)),
                grep("Rictor_", rownames(ScoreUsedData)))
  
  
  TSCmTORname <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  
  TSCmTORScore <- apply(PositiveSingle, 2, sum) + pS6
  
  
  
  
  
  ##########  ApoptosisScore
  linkused <- c(grep("Bak", rownames(ScoreUsedData)),
                grep("Bid", rownames(ScoreUsedData)),
                grep("Bim", rownames(ScoreUsedData)),
                grep("Caspase.", rownames(ScoreUsedData)),
                grep("Bax", rownames(ScoreUsedData)))
  
  Apoptosisnames <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  linkused <- c(grep("Bcl.2", rownames(ScoreUsedData)),
                grep("Bcl.xL", rownames(ScoreUsedData)),
                grep("Bad_", rownames(ScoreUsedData)),
                grep("cIAP", rownames(ScoreUsedData)))
  
  Apoptosisnames <- c( Apoptosisnames,rownames(ScoreUsedData)[linkused])
  NegativeSingle <- ScoreUsedData[linkused , ]
  NegativeSingle[is.na(NegativeSingle)] <- 0
  
  ApoptosisScore <- apply(PositiveSingle, 2, sum) - apply(NegativeSingle, 2, sum)
  
  
  
  ##########  CellCycleScore
  linkused <- c(grep("CDK1", rownames(ScoreUsedData)),
                grep("Cyclin", rownames(ScoreUsedData)),
                grep("p27_", rownames(ScoreUsedData)), 
                grep("PCNA", rownames(ScoreUsedData)))
  
  CellCyclename <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  
  
  CellCycleScore <- apply(PositiveSingle, 2, sum) 
  
  
  
  
  ##########  DNADamageResponseScore
  linkused <- c(grep("53BP1", rownames(ScoreUsedData)),
                grep("ATM", rownames(ScoreUsedData)),
                grep("BRCA2", rownames(ScoreUsedData)),	
                grep("Chk1_", rownames(ScoreUsedData)),
                grep("Chk2_", rownames(ScoreUsedData)),
                grep("Ku80", rownames(ScoreUsedData)),
                grep("Mre11", rownames(ScoreUsedData)),
                grep("PARP", rownames(ScoreUsedData)),
                grep("Rad50", rownames(ScoreUsedData)),
                grep("Rad51", rownames(ScoreUsedData)),
                grep("XRCC1", rownames(ScoreUsedData)),
                grep("p53", rownames(ScoreUsedData)))
  
  DNADamageResponsename <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  DNADamageResponseScore <- apply(PositiveSingle, 2, sum) 
  
  ##########  EMTScore
  linkused <- c(grep("Collagen", rownames(ScoreUsedData)),
                grep("Fibronectin", rownames(ScoreUsedData)),
                grep("N.Cadherin", rownames(ScoreUsedData)))
  
  EMTname <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  linkused <- c(grep("Claudin.7", rownames(ScoreUsedData)),
                grep("E.Cadherin", rownames(ScoreUsedData)))
  
  EMTname <- c(EMTname,rownames(ScoreUsedData)[linkused])
  NegativeSingle <- ScoreUsedData[linkused , ]
  NegativeSingle[is.na(NegativeSingle)] <- 0
  
  EMTScore <- apply(PositiveSingle, 2, sum) - apply(NegativeSingle, 2, sum)
  
  
  
  ##########  HormoneaScore
  linkused <- c(grep("ER-", rownames(ScoreUsedData)),
                match("PR", rownames(ScoreUsedData)))
  
  Hormoneaname <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  HormoneaScore <- apply(PositiveSingle, 2, sum) 
  
  
  
  
  ##########  HormonebScore
  linkused <- c(grep("AR", rownames(ScoreUsedData)),
                grep("INPP4B", rownames(ScoreUsedData)),
                grep("GATA3", rownames(ScoreUsedData)),
                grep("Bcl.2", rownames(ScoreUsedData)))
  
  Hormonebname <- rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  PositiveSingle[is.na(PositiveSingle)] <- 0
  HormonebScore <- apply(PositiveSingle, 2, sum) 
  
  
  AllScore <- data.frame(PI3KAKTScore,RASMAPKScore, RTKScore, TSCmTORScore,
                         ApoptosisScore, CellCycleScore, DNADamageResponseScore,
                         EMTScore, HormoneaScore, HormonebScore)
  
  AllScore$Sample <- rownames(AllScore)
  write.csv(AllScore, file=file.path("RPPARelated", paste(Cond, "-Pathway-Scores.csv", sep = "")))
 # Cond="BRCA"
# AllScore <- read.csv(file=file.path("RPPARelated", paste(Cond, "-Pathway-Scores.csv", sep = "")),header=T)
# rownames(AllScore ) <- AllScore[,1]
# AllScore <- AllScore[,2:ncol(AllScore)]
  ##################  load in cluster information
  clusterfilenameall <- list.files(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/",Cond,"/",sep=""), pattern="-bottom_upper_quantile_class.csv")
  consider.genes <- gsub("-bottom_upper_quantile_class.csv","", clusterfilenameall)
  for(ttt1 in consider.genes){
    Fulldata <- read.csv(file.path("Circadian.gene.cluster_median/",Cond, paste(ttt1,"-bottom_upper_quantile_class.csv",sep="")), header = T)
    Fulldata_Score <- merge(Fulldata,AllScore,by="Sample") 
    Result <- matrix(NA, ncol(AllScore)-1, 2, byrow=TRUE)
    rownames(Result) <- colnames(AllScore)[1:(ncol(AllScore)-1)]
    colnames(Result) <- c("p.value","Diff")
    for(nnn1 in 1:length(rownames(Result))){
      
      tempresult<-try(fit1<-t.test(Fulldata_Score[which(Fulldata_Score$cluster==1),rownames(Result)[nnn1]], 
                                   Fulldata_Score[which(Fulldata_Score$cluster==2),rownames(Result)[nnn1]],silent=TRUE))
      
      if(length(grep("error",tempresult))==0 &  length(grep("Error",tempresult))==0){  
        
        Result[nnn1,"p.value"] <- tempresult$p.value
      }
      
      Result[nnn1,"Diff"] <- mean(Fulldata_Score[which(Fulldata_Score$cluster==2),rownames(Result)[nnn1]])-mean(Fulldata_Score[which(Fulldata_Score$cluster==1),rownames(Result)[nnn1]])
      
    }
    FDRAdjustPvalue <- signif(p.adjust(as.numeric(as.vector(Result[, "p.value"])), method="fdr"),digits = 2)
    OtherAdjustPvalue <- signif(p.adjust(as.numeric(as.vector(Result[, "p.value"])), method="bonferroni"),digits = 2)
    Result <- cbind(Pathway = rownames(Result), Result,  FDR = FDRAdjustPvalue, Bonferroni =OtherAdjustPvalue,Gene = rep(ttt1,times=nrow(Result)))
    Fulldata_Score_eachGene <- Fulldata_Score[,3:ncol(Fulldata_Score)]
    Fulldata_Score_eachGene$Gene <- rep(ttt1,times=nrow(Fulldata_Score_eachGene))
    Fulldata_Score_eachGene$Expression <- Fulldata_Score[,2]
    if(ttt1 == consider.genes[1]){
      ResultAll <- Result
      Fulldata_Score_AllGene <- Fulldata_Score_eachGene
    }else{
      ResultAll <- rbind(ResultAll,Result)
      Fulldata_Score_AllGene <- rbind(Fulldata_Score_AllGene,Fulldata_Score_eachGene)
    }
    
  } 
  ResulteachCancer <- data.frame(ResultAll)
  ResulteachCancer$type <- rep(Cond,times=nrow(ResulteachCancer))
  Fulldata_Score_eachCancer <- Fulldata_Score_AllGene
  Fulldata_Score_eachCancer$type <- rep(Cond,times=nrow(Fulldata_Score_eachCancer))
  if(Cond == gsub("_protein.re_20160627","",filenameall)[1]){
    Fulldata_Score_AllCancer <- Fulldata_Score_eachCancer
    ResultAllCancer <- ResulteachCancer 
  }else{
    Fulldata_Score_AllCancer <- rbind(Fulldata_Score_AllCancer,Fulldata_Score_eachCancer)
    ResultAllCancer <- rbind(ResultAllCancer,ResulteachCancer)
  }
} ########################################  end of RPPA file loop

write.csv(Fulldata_Score_AllCancer,file="AllCancer_AllCir_geneExp_pathwayScore.csv",row.names = F)
write.csv( ResultAllCancer,file="AllCancer_AllCir_pathwayScore_pvalue.csv",row.names = F)
###upper and bottom quantile
#ResultAllCancer_quantile <- read.csv("AllCancer_AllCir_pathwayScore_pvalue_upperbottom_quantile.csv",header=T)
###median
ResultAllCancer <- read.csv("/extraspace/yye1/analysis/Circadian/RPPA/AllCancer_AllCir_pathwayScore_pvalue.csv",header = T)
ResultAllCancer <- ResultAllCancer[(!is.na(ResultAllCancer$p.value)),]
ResultAllCancer$type <- toupper(ResultAllCancer$type)
ResultAllCancer$Pathway <- gsub("Score","",ResultAllCancer$Pathway)
ResultAllCancer$class <- rep("None",times=nrow(ResultAllCancer))
ResultAllCancer["Pathway"][ResultAllCancer["Pathway"]=="PI3KAKT"] <- "PI3K/AKT"
ResultAllCancer["Pathway"][ResultAllCancer["Pathway"]=="RASMAPK"] <- "RAS/MAPK"
ResultAllCancer["Pathway"][ResultAllCancer["Pathway"]=="TSCmTOR"] <- "TSC/mTOR"
ResultAllCancer["Pathway"][ResultAllCancer["Pathway"]=="CellCycle"] <- "Cell Cycle"
ResultAllCancer["Pathway"][ResultAllCancer["Pathway"]=="DNADamageResponse"] <- "DNA Damage Response"
ResultAllCancer["Pathway"][ResultAllCancer["Pathway"]=="Hormonea"] <- "Hormone a"
ResultAllCancer["Pathway"][ResultAllCancer["Pathway"]=="Hormoneb"] <- "Hormone b"



ResultAllCancer["class"][ResultAllCancer["Diff"] < 0 & ResultAllCancer["FDR"] < 0.05] <- "Inhibition"
ResultAllCancer["class"][ResultAllCancer["Diff"] > 0 & ResultAllCancer["FDR"] < 0.05] <- "Activation"

write.csv(ResultAllCancer,file="/extraspace/yye1/analysis/Circadian/RPPA/AllCancer_AllCir_pathwayScore_pvalue.final.csv",quote = F,row.names = F)
pathway_name <- sort(c("PI3K/AKT","RAS/MAPK","RTK", "TSC/mTOR","Apoptosis","Cell Cycle","DNA Damage Response", "EMT","Hormone a","Hormone b" ))

for( i in unique(ResultAllCancer$Gene)){
  cir_pathway <- ResultAllCancer[which(ResultAllCancer$Gene == i),c("Pathway","type","class")]
  if(length(unique(cir_pathway$type)) !=length(tumornameall)){
    cir_pathway <- rbind(data.frame(Pathway=rep(unique(ResultAllCancer$Pathway),length(setdiff(tumornameall,unique(cir_pathway$type)))),
                                    type=rep(setdiff(tumornameall,unique(cir_pathway$type)),each=10),class=rep("None",times=length(setdiff(tumornameall,unique(cir_pathway$type)))*10)),cir_pathway)
  }
  #cir_pathway$FDR <- -log10(cir_pathway$FDR)
  #cir_pathway["FDR"][cir_pathway["FDR"] >= 4] <- 4
  pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/",i,"_pathway FDR sample morethan 50.pdf",sep=""),width = 9,height = 4.5)
  p <- ggplot(cir_pathway,aes(x=type,y=Pathway))+
    geom_tile(aes(fill=factor(class)),col="white")+
    scale_fill_manual(limits=c("Activation","Inhibition","None"),values = c("red","blue","lightgray"),na.value="white",labels=c("Activation","Inhibition","None"),name=i)+
    scale_x_discrete(limits=tumornameall[!(tumornameall %in% c("ACC","CHOL","DLBC","UCS"))])+
    scale_y_discrete(limits=pathway_name,label=pathway_name)+
    theme(panel.background=element_rect(colour=NA,fill="white",size=2),
          panel.grid=element_line(colour="white",linetype="dashed"),
          panel.grid.major=element_line(colour="white",linetype="dashed"),
          axis.title=element_blank(),
          axis.text.y=element_text(size=16,colour = "black"),
          axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
          axis.ticks=element_line(color="black"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=14),legend.position="bottom",legend.direction="horizontal")

  print(p)
  dev.off()
 # print(c(i,table(cir_pathway$class)))
}
###Only for CRY1/2, PER1/2/3, RORA/B/C
Select_core <- c("CRY1","CRY2","PER1","PER2","PER3","RORA","RORB","RORC")
Select_cir_pathway <- ResultAllCancer[which(ResultAllCancer$Gene %in% Select_core),]

pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/","CRY_PER_ROR_pathway FDR sample morethan 50.pdf",sep=""),width = 7,height = 10)
ggplot(Select_cir_pathway,aes(x=type,y=Pathway))+
  geom_tile(aes(fill=factor(class)),col="white")+
  scale_fill_manual(limits=c("Activation","Inhibition","None"),values = c("red","blue","lightgray"),na.value="white",labels=c("Activation","Inhibition","None"),name="")+
  scale_x_discrete(limits=tumornameall[!(tumornameall %in% c("ACC","CHOL","DLBC","UCS"))])+
  scale_y_discrete(limits=pathway_name,label=pathway_name)+
  facet_grid(Gene~.)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=9,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),strip.background = element_blank(),strip.text = element_text(size=14),
        legend.title=element_text(size=10),legend.position="bottom",legend.direction="horizontal",
        legend.key.size = unit(0.3, "cm"),
        panel.spacing = unit(0.1, "lines"))
dev.off()
###Only for CRY1/2
Select_core <- c("CRY1","CRY2")
Select_cir_pathway <- ResultAllCancer[which(ResultAllCancer$Gene %in% Select_core),]

pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/","CRY_pathway FDR sample morethan 50.pdf",sep=""),width = 7,height = 4)
ggplot(Select_cir_pathway,aes(x=type,y=Pathway))+
  geom_tile(aes(fill=factor(class)),col="white")+
  scale_fill_manual(limits=c("Activation","Inhibition","None"),values = c("red","blue","lightgray"),na.value="white",labels=c("Activation","Inhibition","None"),name="")+
  scale_x_discrete(limits=tumornameall[!(tumornameall %in% c("ACC","CHOL","DLBC","UCS"))])+
  scale_y_discrete(limits=pathway_name,label=pathway_name)+
  facet_grid(Gene~.)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=9,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),strip.background = element_blank(),strip.text = element_text(size=14),
        legend.title=element_text(size=10),legend.position="bottom",legend.direction="horizontal",
        legend.key.size = unit(0.3, "cm"),
        panel.spacing = unit(0.1, "lines"))
dev.off()
###Only for RORA/B/C
Select_core <- c("RORA","RORB","RORC")
Select_cir_pathway <- ResultAllCancer[which(ResultAllCancer$Gene %in% Select_core),]
Select_cir_pathway_EMT <- Select_cir_pathway[which(Select_cir_pathway$Pathway=="EMT"),]
pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/","RORs_EMT_pathway FDR sample morethan 50.pdf",sep=""),width = 7,height = 2)
ggplot(Select_cir_pathway_EMT,aes(x=type,y=Gene))+
  geom_tile(aes(fill=factor(class)),col="white")+
  scale_fill_manual(limits=c("Activation","Inhibition","None"),values = c("red","blue","lightgray"),na.value="white",labels=c("Activation","Inhibition","None"),name="")+
  scale_x_discrete(limits=tumornameall[!(tumornameall %in% c("ACC","CHOL","DLBC","UCS"))])+
#  scale_y_discrete(limits=pathway_name,label=pathway_name)+
 # facet_grid(Gene~.)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=9,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=10),legend.position="bottom",legend.direction="horizontal",
        legend.key.size = unit(0.3, "cm"))
dev.off()


###Each cancer types
for( i in unique(ResultAllCancer$type)){
  cir_pathway <- ResultAllCancer[which(ResultAllCancer$type == i),c("Pathway","Gene","class")]
  if(length(unique(cir_pathway$Gene)) !=length(circadian.genes)){
    cir_pathway <- rbind(data.frame(Pathway=rep(unique(ResultAllCancer$Pathway),length(setdiff(circadian.genes,unique(cir_pathway$Gene)))),
                                    Gene=rep(setdiff(circadian.genes,unique(cir_pathway$Gene)),each=10),class=rep("None",times=length(setdiff(circadian.genes,unique(cir_pathway$Gene)))*10)),cir_pathway)
  }
  #cir_pathway$FDR <- -log10(cir_pathway$FDR)
  #cir_pathway["FDR"][cir_pathway["FDR"] >= 4] <- 4
  pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/",i,"_circadianGenes_pathway FDR sample morethan 50.pdf",sep=""),width = 9,height = 4.5)
  p <- ggplot(cir_pathway,aes(x=Gene,y=Pathway))+
    geom_tile(aes(fill=factor(class)),col="white")+
    scale_fill_manual(limits=c("Activation","Inhibition","None"),values = c("red","blue","lightgray"),na.value="white",labels=c("Activation","Inhibition","None"),name=i)+
    scale_x_discrete(limits=circadian.genes)+
    scale_y_discrete(limits=pathway_name,label=pathway_name)+
    theme(panel.background=element_rect(colour=NA,fill="white",size=2),
          panel.grid=element_line(colour="white",linetype="dashed"),
          panel.grid.major=element_line(colour="white",linetype="dashed"),
          axis.title=element_blank(),
          axis.text.y=element_text(size=16,colour = "black"),
          axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
          axis.ticks=element_line(color="black"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=14),legend.position="bottom",legend.direction="horizontal")
  
  print(p)
  dev.off()
  # print(c(i,table(cir_pathway$class)))
}

ResultAllCancer.sig <- ResultAllCancer[which(ResultAllCancer$class !="None"),]

ggplot(ResultAllCancer.sig ,aes(x=CancerType,fill=factor(Class)))+
  geom_bar(color=NA,width = 0.5)+
  scale_fill_manual(limit=c("neg","pos"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(expand = c(0.02,0),breaks = seq(0,2500,length.out = 6))+
  scale_x_discrete(limit= names(table(Ac_cirAll$CancerType))[order(table(Ac_cirAll$CancerType))],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_blank(),#line(linetype="dashed",color="lightgray"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))



###Example of CRY2
i= "CRY2"#,"CRY2"
change_cancer <- cir.exp.pval[which(cir.exp.pval$geneSymbol==i & abs(cir.exp.pval$Fold) >= log2(1.5) & cir.exp.pval$Pvalue <= 0.05 ),"Type"]
change_cancer <- as.character(change_cancer)
change_cancer_FDR <- paste(change_cancer,"_FDR",sep="")
change_cancer_pval <- paste(change_cancer,"_pval",sep="")
cor_data <- read.delim(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/sign_cor/",i,"correlation_across_cancer.txt",sep=""),header = T)
change_cor_data <- cor_data#[,c("gene",change_cancer,change_cancer_FDR)]
DNADamageGenes <- c("TP53BP1","ATM","BRCA2","CHEK1","CHEK2","XRCC5","MRE11A","TP53","RAD50","RAD51","XRCC1")
DNADamage.cor <- change_cor_data[which(change_cor_data$gene %in% DNADamageGenes),]
DNADamage.exp <- cancer_ty.genes.exp.All[which(cancer_ty.genes.exp.All$gene %in% DNADamageGenes),]
CellCycleGenes <- c("CDK1","CCNB1","CCND1","CCNE1","CCNE2","CDKN1B","PCNA")
CellCycle.cor <- change_cor_data[which(change_cor_data$gene %in% CellCycleGenes),]
CellCycle.exp <- cancer_ty.genes.exp.All[which(cancer_ty.genes.exp.All$gene %in% CellCycleGenes),]


###
LUAD<- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/","LUAD_mRNA_each_exp_20160513",sep=""),header=T)
LUAD$gene <- data.frame(do.call(rbind,strsplit(as.character(LUAD$gene),"\\|")))$X1
filter <- colnames(LUAD)[as.numeric(substr(colnames(LUAD),14,15)) ==1]
filter <- filter[!is.na(filter)]
LUAD_DNADamage <- LUAD[which(LUAD$gene %in% c("CRY2","BRCA2","CHEK1","CHEK2","RAD51")),c("gene",filter)]
LUAD_DNADamage <- data.frame(t(as.matrix(LUAD_DNADamage)))
colnames(LUAD_DNADamage) <- LUAD_DNADamage[1,]
LUAD_DNADamage <- LUAD_DNADamage[2:nrow(LUAD_DNADamage),]
LUAD_DNADamage.m <- melt(LUAD_DNADamage,id.vars = "CRY2",measure.vars = c("BRCA2","CHEK1","CHEK2","RAD51"))
levels(LUAD_DNADamage.m$variable) <- paste(levels(LUAD_DNADamage.m$variable), " ( R=",signif(DNADamage.cor[which(DNADamage.cor$gene %in% c("BRCA2","CHEK1","CHEK2","RAD51")),"LUAD"],digits=2),")",sep="")
pdf("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/LUAD_CRY2_DNADamage/LUAD_CRY2_DNADamgeGenes.pdf",width = 8,height = 8)
ggplot(data = LUAD_DNADamage.m , aes(x = log2(as.numeric(CRY2)), 
                                     y = log2(as.numeric(value))))+
  geom_point(size=3)+
  stat_smooth(method = "lm",color="red",size=1)+
  facet_wrap(~variable,scales = "free_y")+
  ylab('Gene Expression')+xlab("Gene Expression (CRY2)")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=16),
        strip.background = element_rect(fill=NA))
dev.off()

pair_sampleAll <- read.delim("~/Circadian/expression/pair_sampleAll.txt",header=T)
LUAD_DNADamage.exp <- LUAD[which(LUAD$gene %in% c("CRY2","BRCA2","CHEK1","CHEK2","RAD51")),]
LUAD_DNADamage.exp <- data.frame(t(as.matrix(LUAD_DNADamage.exp)))
colnames(LUAD_DNADamage.exp) <- LUAD_DNADamage.exp[1,]
LUAD_DNADamage.exp <- LUAD_DNADamage.exp[2:nrow(LUAD_DNADamage.exp),]
LUAD_DNADamage.exp$barcode <- gsub("\\.","\\-",row.names(LUAD_DNADamage.exp))
LUAD_DNADamage.exp <- merge(LUAD_DNADamage.exp,pair_sampleAll[,1:2],by="barcode")
LUAD_DNADamage.exp.m <- melt(LUAD_DNADamage.exp,id.vars = "class",measure.vars =c("CRY2","BRCA2","CHEK1","CHEK2","RAD51") )
pdf("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/LUAD_CRY2_DNADamage/boxplot_LUAD_CRY2_DNADamgeGenes.pdf",width = 6,height = 2.5)
ggplot(LUAD_DNADamage.exp.m, aes( x = variable, y = log2(as.numeric(value)), fill = class))+  #stat_boxplot(geom = "errorbar",group = class, width = 0.2 )+
  geom_boxplot(width = 0.6,outlier.shape = NA)+
  scale_fill_manual( limits= c("tumor","normal"),values = c("red","blue"),labels = c("T","N"),name="Class")+
  # geom_text(data = LGG_sig_Pval,aes(y = 11.5,x=1, label = paste("p = ",t.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
  # scale_y_continuous(breaks=seq(0,12,by=4))+
  ylab('Gene Expression')+xlab("")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour="lightgray",size=0.1,linetype = "dashed"),panel.grid.minor=element_line(colour=NA),
        axis.text=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.ticks.length = unit(.25, "cm"))
dev.off()

pathwayGenes <- c(CellCycleGenes = c("CDK1","CCNB1","CCND1","CCNE1","CCNE2","CDKN1B","PCNA"),
                  DNADamageGenes = c("TP53BP1","ATM","BRCA2","CHEK1","CHEK2","XRCC5","MRE11A","TP53","RAD50","RAD51","XRCC1"),
                  ApoptosisGenes = c("BAK1","BAX","BID","BCL2L11","CASP7","BAD","BCL2","BCL2L1","BIRC2"),
                  EMTGenes = c("CDH1","CLDN7","COL6A1","CDH2","FN1"),
                  HormoneaGenes = c("ESR1","PGR"),
                  HormonebGenes = c("AR","BCL2","INPP4B","GATA3"),
                  RASMAPKGenes = c("ARAF","JUN","RAF1","MAPK8","MAPK1","MAP2K1","MAPK14","RPS6KA1","SHC1","YBX1"),
                  TSCmTORGenes = c("EIF4EBP1","MTOR","RPS6KB1","RICTOR","RPS6"))
pathwayGenes <- data.frame(pathwayGenes)
pathwayGenes$Pathway <- gsub("[1-9]","",rownames(pathwayGenes))
pathwayGenescor <- change_cor_data[which(change_cor_data$gene %in%pathwayGenes$pathwayGenes),c("gene","BRCA")]
pathwayGenescor <- merge(pathwayGenescor,pathwayGenes,by.x="gene",by.y="pathwayGenes")
CellCycle.cor <- change_cor_data[which(change_cor_data$gene %in% CellCycleGenes),]
CellCycle.exp <- cancer_ty.genes.exp.All[which(cancer_ty.genes.exp.All$gene %in% CellCycleGenes),]


###
BRCA<- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/","BRCA_mRNA_each_exp_20160513",sep=""),header=T)
BRCA$gene <- data.frame(do.call(rbind,strsplit(as.character(BRCA$gene),"\\|")))$X1
filter <- colnames(BRCA)[as.numeric(substr(colnames(BRCA),14,15)) ==1]
filter <- filter[!is.na(filter)]
BRCA_DNADamage <- BRCA[which(BRCA$gene %in% c("CRY2","BRCA2","CHEK1","CHEK2","RAD51")),c("gene",filter)]
BRCA_DNADamage <- data.frame(t(as.matrix(BRCA_DNADamage)))
colnames(BRCA_DNADamage) <- BRCA_DNADamage[1,]
BRCA_DNADamage <- BRCA_DNADamage[2:nrow(BRCA_DNADamage),]
BRCA_DNADamage.m <- melt(BRCA_DNADamage,id.vars = "CRY2",measure.vars = c("BRCA2","CHEK1","CHEK2","RAD51"))
levels(BRCA_DNADamage.m$variable) <- paste(levels(BRCA_DNADamage.m$variable), " ( R=",signif(DNADamage.cor[which(DNADamage.cor$gene %in% c("BRCA2","CHEK1","CHEK2","RAD51")),"BRCA"],digits=2),")",sep="")
pdf("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/BRCA/BRCA_CRY2_DNADamgeGenes.pdf",width = 8,height = 8)
ggplot(data = BRCA_DNADamage.m , aes(x = log2(as.numeric(CRY2)), 
                                     y = log2(as.numeric(value))))+
  geom_point(size=3)+
  stat_smooth(method = "lm",color="red",size=1)+
  facet_wrap(~variable,scales = "free_y")+
  ylab('Gene Expression')+xlab("Gene Expression (CRY2)")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=16),
        strip.background = element_rect(fill=NA))
dev.off()

pair_sampleAll <- read.delim("~/Circadian/expression/pair_sampleAll.txt",header=T)
BRCA_DNADamage.exp <- BRCA[which(BRCA$gene %in% c("CRY2","BRCA2","CHEK1","CHEK2","RAD51")),]
BRCA_DNADamage.exp <- data.frame(t(as.matrix(BRCA_DNADamage.exp)))
colnames(BRCA_DNADamage.exp) <- BRCA_DNADamage.exp[1,]
BRCA_DNADamage.exp <- BRCA_DNADamage.exp[2:nrow(BRCA_DNADamage.exp),]
BRCA_DNADamage.exp$barcode <- gsub("\\.","\\-",row.names(BRCA_DNADamage.exp))
BRCA_DNADamage.exp <- merge(BRCA_DNADamage.exp,pair_sampleAll[,1:2],by="barcode")
BRCA_DNADamage.exp.m <- melt(BRCA_DNADamage.exp,id.vars = "class",measure.vars =c("CRY2","BRCA2","CHEK1","CHEK2","RAD51") )
pdf("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/BRCA/boxplot_BRCA_CRY2_DNADamgeGenes.pdf",width = 6,height = 2.5)
ggplot(BRCA_DNADamage.exp.m, aes( x = variable, y = log2(as.numeric(value)), fill = class))+  #stat_boxplot(geom = "errorbar",group = class, width = 0.2 )+
  geom_boxplot(width = 0.6,outlier.shape = NA)+
  scale_fill_manual( limits= c("tumor","normal"),values = c("red","blue"),labels = c("T","N"),name="Class")+
  # geom_text(data = LGG_sig_Pval,aes(y = 11.5,x=1, label = paste("p = ",t.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
  # scale_y_continuous(breaks=seq(0,12,by=4))+
  ylab('Gene Expression')+xlab("")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour="lightgray",size=0.1,linetype = "dashed"),panel.grid.minor=element_line(colour=NA),
        axis.text=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.ticks.length = unit(.25, "cm"))
dev.off()







###
#circadian pathway network
for( i in unique(ResultAllCancer.sig$type)){
  BRCA.nodes <- ResultAllCancer.sig[which(ResultAllCancer.sig$type==i),]
  if(nrow(BRCA.nodes) >= 10){
    BRCA.pathway <- as.vector(t(BRCA.nodes[,1]))
    usedpathway <- names(table(BRCA.pathway))[table(BRCA.pathway) >=1]
    BRCA.Cir <- as.vector(t(BRCA.nodes[,"Gene"]))
    usedcir <- names(table(BRCA.Cir))[table(BRCA.Cir) >=1]
    #link <- which(BRCA.nodes[,"Pathway"] %in% )
    temp <- BRCA.nodes[,c("Gene","Pathway")]
   
    df <- graph.data.frame(BRCA.nodes[,c("Gene","Pathway","class","FDR")],directed=T)
    df1 <- graph.data.frame(data.frame(First = unique(temp[,1]), second = unique(temp[,1])))
    tt1 <- layout.circle(df1)
    df2 <- graph.data.frame(data.frame(First = unique(temp[,2]), second = unique(temp[,2])))
    tt2 <- layout.circle(df2)
    tt <- layout.circle(df)
    #########  rescale two layout
    n1=length(unique(temp[,1]))
    n2=length(unique(temp[,2]))
    tt[1:n1,1] <- tt1[,1]*1.7
    tt[1:n1,2] <- tt1[,2]*1.7
    tt[(n1+1):(n1+n2),1] <- tt2[,1]
    tt[(n1+1):(n1+n2),2] <- tt2[,2]
    sig.cir <- as.vector(unique(BRCA.nodes$Gene)[unique(BRCA.nodes$Gene) %in% core.circadian])
    ttname <- V(df)$name
    V(df)$color <- "gold"
    V(df)[1:n1]$color <- "yellowgreen"
    V(df)[ttname %in% sig.cir]$color <- "tomato"
    pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/",i,"_circadian_pathway_network.pdf",sep=""),width=10,height=10)
    P <- plot(df,edge.color=c("red","blue")[(E(df)$class=="Inhibition")+1],
              layout=tt,vertex.label.color="black",edge.arrow.size=.4,vertex.size=15,vertex.label.cex=0.8)
    dev.off()
  }
}


for( i in unique(ResultAllCancer.sig$type)){
  BRCA.nodes <- ResultAllCancer.sig[which(ResultAllCancer.sig$type==i & ResultAllCancer.sig$Gene %in% core.circadian),]
  if(nrow(BRCA.nodes) >= 10){
    BRCA.pathway <- as.vector(t(BRCA.nodes[,1]))
    usedpathway <- names(table(BRCA.pathway))[table(BRCA.pathway) >=1]
    BRCA.Cir <- as.vector(t(BRCA.nodes[,"Gene"]))
    usedcir <- names(table(BRCA.Cir))[table(BRCA.Cir) >=1]
    #link <- which(BRCA.nodes[,"Pathway"] %in% )
    temp <- BRCA.nodes[,c("Gene","Pathway")]
    
    df <- graph.data.frame(BRCA.nodes[,c("Gene","Pathway","class","FDR")],directed=T)
    df1 <- graph.data.frame(data.frame(First = unique(temp[,1]), second = unique(temp[,1])))
    tt1 <- layout.circle(df1)
    df2 <- graph.data.frame(data.frame(First = unique(temp[,2]), second = unique(temp[,2])))
    tt2 <- layout.circle(df2)
    tt <- layout.circle(df)
    #########  rescale two layout
    n1=length(unique(temp[,1]))
    n2=length(unique(temp[,2]))
    tt[1:n1,1] <- tt1[,1]*1.7
    tt[1:n1,2] <- tt1[,2]*1.7
    tt[(n1+1):(n1+n2),1] <- tt2[,1]
    tt[(n1+1):(n1+n2),2] <- tt2[,2]
    sig.cir <- as.vector(unique(BRCA.nodes$Gene)[unique(BRCA.nodes$Gene) %in% core.circadian])
    ttname <- V(df)$name
    V(df)$color <- "gold"
    V(df)[1:n1]$color <- "tomato"
  #  V(df)[ttname %in% sig.cir]$color <- "tomato"
    pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/",i,"_core_circadian_pathway_network.pdf",sep=""),width=10,height=10)
    P <- plot(df,edge.color=c("red","blue")[(E(df)$class=="Inhibition")+1],
              layout=tt,vertex.label.color="black",edge.arrow.size=.4,vertex.size=15,vertex.label.cex=0.8)
    dev.off()
  }
}

df <- data.frame(time = c("201501", "201501", "201502", "201502", "201503", "201503"), 
factors = c("x", "y", "x", "y", "x", "y"),
values = c(0.33, 0.15, 0.36, 0.15, 0.39, 0.18))

seq_adj <- seq(from = 1, to = nrow(df)-1, by = 2)

df$val_mod <- ""
df$val_mod[seq_adj] <- paste(as.character(df$values[seq_adj] * 100), "%(",
                             as.character(df$values[seq_adj + 1] * 100), "%)", sep = "")
df$val_mod2 <- df$values
df$val_mod2[seq_adj] <- df$values[seq_adj] - df$values[seq_adj + 1]

ggplot(data = df, aes(x = factor(1), y = val_mod2, fill = factor(factors))) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_y_continuous(limits = c(0,1)) +
  coord_polar("y") +
  facet_grid(.~time) +
  geom_text(aes(x = factor(1), y= .5, label = val_mod, vjust = 4.5)) +
  theme(axis.text.x=element_blank()) +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(legend.title=element_blank()) +
  scale_fill_discrete(labels=c("x","y"))

Cir_Pathway_Gene_sig <- ResultAllCancer.sig[(ResultAllCancer.sig$class != "None" & ResultAllCancer.sig$Gene %in% core.circadian),]
ResultAllCancer.m <- ResultAllCancer#[!(ResultAllCancer$type %in% c("ACC","CHOL","DLBC","UCS")),]
ResultAllCancer$
Cir_pathway_sig <- t(sapply(split(ResultAllCancer.m[,"class"],list(ResultAllCancer.m$Pathway,ResultAllCancer.m$Gene)),function(x){c(length(x[x=="Activation"]),length(x[x=="Inhibition"]),length(x[x=="None"]))}))
Cir_pathway_sig <- data.frame(Cir_pathway_sig)
colnames(Cir_pathway_sig) <- c("Activation","Inhibition","None")
Cir_pathway_sig$Pathway <- data.frame(do.call(rbind,strsplit(as.character(rownames(Cir_pathway_sig)),"\\.")))$X1
Cir_pathway_sig$Gene <- data.frame(do.call(rbind,strsplit(as.character(rownames(Cir_pathway_sig)),"\\.")))$X2
Cir_pathway_sig.m <- melt(Cir_pathway_sig,id.vars = c("Pathway","Gene"),measure.vars = c("Activation","Inhibition","None"))
Cir_pathway_sig.m$per <- Cir_pathway_sig.m$value/31
#For only circadian genes, width=10
Cir_pathway_sig.mm <- Cir_pathway_sig.m[which(Cir_pathway_sig.m$Gene %in% core.circadian),]
strip.x.color <- rep("black",times=length(unique(Cir_pathway_sig.m$Gene)))
strip.x.color[sort(unique(Cir_pathway_sig.m$Gene)) %in% core.circadian] <- "red"
pdf("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/Percentage of cancer types in core circadian genes with pathway.pdf",width = 9,height=5)
ggplot(data = Cir_pathway_sig.mm, aes(x = factor(1), y = per, fill = factor(variable))) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_y_continuous(limits = c(0,1)) +
  coord_polar("y") +
  facet_grid(Pathway~Gene) +
 # geom_text(aes(x = factor(1), y= .5, label = val_mod, vjust = 4.5)) +
  theme(axis.text=element_blank(),axis.title = element_blank(), 
        panel.background = element_blank(),
        legend.title=element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle =0,hjust=0,color="black",size=11),strip.background = element_blank(),
        strip.text.x = element_text(color=strip.x.color,size=11,angle = 90,vjust = 0),legend.text = element_text(size=14),
        legend.position = "bottom",panel.spacing  = unit(0.02, "lines")) +
  scale_fill_manual(limits=c("Activation","Inhibition","None"),values=c("red","blue","lightgray"))
dev.off()
Cir_pathway_sig.mmm <- Cir_pathway_sig.m[!(Cir_pathway_sig.m$Gene %in% core.circadian),]
pdf("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/Percentage of cancer types in 37 related circadian genes with pathway.pdf",width = 12,height=5)
ggplot(data = Cir_pathway_sig.mmm, aes(x = factor(1), y = per, fill = factor(variable))) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_y_continuous(limits = c(0,1)) +
  coord_polar("y") +
  facet_grid(Pathway~Gene) +
  # geom_text(aes(x = factor(1), y= .5, label = val_mod, vjust = 4.5)) +
  theme(axis.text=element_blank(),axis.title = element_blank(), 
        panel.background = element_blank(),
        legend.title=element_blank(), axis.ticks = element_blank(),
        strip.text.y = element_text(angle =0,hjust=0,color="black",size=11),strip.background = element_blank(),
        strip.text.x = element_text(color=strip.x.color,size=11,angle = 90,vjust = 0),legend.text = element_text(size=14),
        legend.position = "bottom",panel.spacing  = unit(0.02, "lines")) +
  scale_fill_manual(limits=c("Activation","Inhibition","None"),values=c("red","blue","lightgray"))
dev.off()
pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/","Network of circadian_pathway nodes than 40 beside BRCA.pdf",sep=""),width=12,height=15)
par(mfrow = c(5,4),mai = c(1, 0.1, 0.1, 0.1),oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( i in unique(ResultAllCancer.sig$type)[!(unique(ResultAllCancer.sig$type) %in% "BRCA")]){
  BRCA.nodes <- ResultAllCancer.sig[which(ResultAllCancer.sig$type==i),]
  if(nrow(BRCA.nodes) >= 40){  #print(c(nrow(BRCA.nodes),i))}}
    BRCA.pathway <- as.vector(t(BRCA.nodes[,1]))
    usedpathway <- names(table(BRCA.pathway))[table(BRCA.pathway) >=1]
    BRCA.Cir <- as.vector(t(BRCA.nodes[,"Gene"]))
    usedcir <- names(table(BRCA.Cir))[table(BRCA.Cir) >=1]
    #link <- which(BRCA.nodes[,"Pathway"] %in% )
    temp <- BRCA.nodes[,c("Gene","Pathway")]
    
    df <- graph.data.frame(BRCA.nodes[,c("Gene","Pathway","class","FDR")],directed=T)
    df1 <- graph.data.frame(data.frame(First = unique(temp[,1]), second = unique(temp[,1])))
    tt1 <- layout.circle(df1)
    df2 <- graph.data.frame(data.frame(First = unique(temp[,2]), second = unique(temp[,2])))
    tt2 <- layout.circle(df2)
    tt <- layout.circle(df)
    #########  rescale two layout
    n1=length(unique(temp[,1]))
    n2=length(unique(temp[,2]))
    tt[1:n1,1] <- tt1[,1]*1.7
    tt[1:n1,2] <- tt1[,2]*1.7
    tt[(n1+1):(n1+n2),1] <- tt2[,1]
    tt[(n1+1):(n1+n2),2] <- tt2[,2]
    sig.cir <- as.vector(unique(BRCA.nodes$Gene)[unique(BRCA.nodes$Gene) %in% core.circadian])
    ttname <- V(df)$name
    V(df)$color <- "gold"
    V(df)[1:n1]$color <- "yellowgreen"
    V(df)[ttname %in% sig.cir]$color <- "tomato"
    
    plot(df,edge.color=c("red","blue")[(E(df)$class=="Inhibition")+1],
              layout=tt,vertex.label.color="black",edge.arrow.size=.4,vertex.size=15,vertex.label.cex=0.8,main=i)
   # dev.off()
  }
}
dev.off()
data.frame(sapply(split(ResultAllCancer.sig[,"type"],ResultAllCancer.sig$class),table))
pdf("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/count the number of circadian genes and pathway interaction in each cancer type.pdf",width=2,height = 5)#,width=4000,height=1500,res=400)
ggplot(ResultAllCancer.sig,aes(x=type,fill=factor(class)))+
  geom_bar(color=NA,width = 0.5)+
  coord_flip()+
  scale_fill_manual(limit=c("Inhibition","Activation"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(limit= c(-1,290),expand = c(0,0),breaks = seq(0,250,length.out = 6))+
  scale_x_discrete(limit= names(table(ResultAllCancer.sig$type))[order(table(ResultAllCancer.sig$type))],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=0.5),
        panel.grid.major=element_line(linetype="dashed",color="lightgray",size=0.1),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black"),
        axis.ticks.x=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_hline(yintercept = mean(table(ResultAllCancer.sig$type)), colour = "black",linetype="dashed")
dev.off()
#######






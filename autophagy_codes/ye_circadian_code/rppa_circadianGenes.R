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
    write.csv(sub1,file.path("Circadian.gene.cluster_median", tumornameall[m], paste(n, "-bottom_upper_quantile_class.csv", sep = "")),row.names = F)
    
  }
}

library(gdata)
library(survival)
library(plotrix)  
filenameall <- list.files("/extraspace/TCGA/TCGA_RPPA/", pattern="*csv")

Output <- list()
usedid <- 1
for(kkk1 in 1:length(filenameall))
{
  
  
  filename <- filenameall[kkk1]
  
  Cond <- gsub("rppa_|\\.csv","",filename)
  Fulldata <- read.csv(file.path("/extraspace/TCGA/TCGA_RPPA/", filename))
  if(Cond == "SKCM")
  {
    keeplink <- which(toupper(Fulldata[, "Sample_type"]) == "MET")
    Fulldata <- Fulldata[keeplink,]
  }else
  {
    keeplink <- which(toupper(Fulldata[, "Sample_type"]) == "PRIMARY")
    Fulldata <- Fulldata[keeplink,]
    
    
  }
  rownames(Fulldata) <- as.vector(Fulldata)[,1]
  
  FullRPPAData <- Fulldata[, c(grep("14\\.3\\.3_beta", colnames(Fulldata)):dim(Fulldata)[2])]
  FullRPPAData <- t(FullRPPAData)
  
  
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



####################################
##################  read in RPPA data
library(gdata)
library(survival)
library(plotrix)  
filenameall <- list.files("/extraspace/TCGA/TCGA_RPPA/", pattern="*csv")
#library(foreach)

for(kkk1 in 1:length(filenameall))
{
  
  
  filename <- filenameall[kkk1]
  
  #filename <- "TCGA-BLCA-rnaexpr.tsv"
  Cond <- gsub("rppa_|\\.csv","",filename)
  
  #if (!file.exists(file.path("RPPARelated", Cond))) {
  #  dir.create(file.path("RPPARelated", Cond))
  #}
  
  
  Fulldata <- read.csv(file.path("/extraspace/TCGA/TCGA_RPPA/", filename))
  if(Cond == "skcm")
  {
    keeplink <- which(toupper(Fulldata[, "Sample_type"]) == "MET")
    Fulldata <- Fulldata[keeplink,]
  }else
  {
    keeplink <- which(toupper(Fulldata[, "Sample_type"]) == "PRIMARY")
    Fulldata <- Fulldata[keeplink,]
    
    
  }
  rownames(Fulldata) <- as.vector(Fulldata[,1])
  FullRPPAData <- Fulldata[, c(grep("14\\.3\\.3_beta", colnames(Fulldata)):dim(Fulldata)[2])]
  FullRPPAData <- t(FullRPPAData)
  
  outputp <- NA
  countid <- 2
  
  
  ##############   calculate signature
  
  ScoreUsedData <- sweep(FullRPPAData, 1, apply(FullRPPAData, 1 , median), FUN="-")
  ScoreUsedData <- sweep(ScoreUsedData, 1, apply(FullRPPAData, 1 , sd), FUN="/")
  
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
  
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  linkused <- c(grep("INPP4B", rownames(ScoreUsedData)),
                grep("PTEN", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  NegativeSingle <- ScoreUsedData[linkused , ]
  
  
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
                grep("YB.1_p", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  
  RASMAPKScore <- apply(PositiveSingle, 2, sum)
  
  
  
  ##########  RTKScore
  linkused <- c(grep("HER2_", rownames(ScoreUsedData)),
                grep("HER3_", rownames(ScoreUsedData)),
                grep("Ret_", rownames(ScoreUsedData)),
                grep("Shc_p", rownames(ScoreUsedData)), 
                grep("Src_p", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  
  RTKScore <- apply(PositiveSingle, 2, sum) + pEGFR
  
  
  
  ##########  TSCmTORScore
  linkused <- c(grep("4E.BP1_p", rownames(ScoreUsedData)),
                grep("mTOR_", rownames(ScoreUsedData)),
                grep("p70S6K_", rownames(ScoreUsedData)),
                grep("Rictor_", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  
  TSCmTORScore <- apply(PositiveSingle, 2, sum) + pS6
  
  
  
  
  
  ##########  ApoptosisScore
  linkused <- c(grep("Bak", rownames(ScoreUsedData)),
                grep("Bid", rownames(ScoreUsedData)),
                grep("Bim", rownames(ScoreUsedData)),
                grep("Caspase.", rownames(ScoreUsedData)),
                grep("Bax", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  linkused <- c(grep("Bcl.2", rownames(ScoreUsedData)),
                grep("Bcl.xL", rownames(ScoreUsedData)),
                grep("Bad_", rownames(ScoreUsedData)),
                grep("cIAP", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  NegativeSingle <- ScoreUsedData[linkused , ]
  
  
  ApoptosisScore <- apply(PositiveSingle, 2, sum) - apply(NegativeSingle, 2, sum)
  
  
  
  ##########  CellCycleScore
  linkused <- c(grep("CDK1", rownames(ScoreUsedData)),
                grep("Cyclin", rownames(ScoreUsedData)),
                grep("p27_", rownames(ScoreUsedData)), 
                grep("PCNA", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  
  
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
                match("p53", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  DNADamageResponseScore <- apply(PositiveSingle, 2, sum) 
  
  ##########  EMTScore
  linkused <- c(grep("Collagen", rownames(ScoreUsedData)),
                grep("Fibronectin", rownames(ScoreUsedData)),
                grep("N.Cadherin", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  linkused <- c(grep("Claudin.7", rownames(ScoreUsedData)),
                grep("E.Cadherin", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  NegativeSingle <- ScoreUsedData[linkused , ]
  
  
  EMTScore <- apply(PositiveSingle, 2, sum) - apply(NegativeSingle, 2, sum)
  
  
  
  ##########  HormoneaScore
  linkused <- c(grep("ER\\.", rownames(ScoreUsedData)),
                match("PR", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  HormoneaScore <- apply(PositiveSingle, 2, sum) 
  
  
  
  
  ##########  HormonebScore
  linkused <- c(match("AR", rownames(ScoreUsedData)),
                grep("INPP4B", rownames(ScoreUsedData)),
                grep("GATA3", rownames(ScoreUsedData)),
                grep("Bcl.2", rownames(ScoreUsedData)))
  
  rownames(ScoreUsedData)[linkused]
  PositiveSingle <- ScoreUsedData[linkused , ]
  
  HormonebScore <- apply(PositiveSingle, 2, sum) 
  
  
  AllScore <- data.frame(PI3KAKTScore,RASMAPKScore, RTKScore, TSCmTORScore,
                         ApoptosisScore, CellCycleScore, DNADamageResponseScore,
                         EMTScore, HormoneaScore, HormonebScore)
  
  AllScore$Sample <- rownames(AllScore)
  write.csv(AllScore, file=file.path("RPPARelated", paste(Cond, "-Pathway-Scores.csv", sep = "")))
  
  
  
  ##################  load in cluster information
  clusterfilenameall <- list.files(paste("Circadian.gene.cluster/",Cond,"/",sep=""), pattern="-bottom_upper_quantile_class.csv")
  consider.genes <- gsub("-bottom_upper_quantile_class.csv","", clusterfilenameall)
  for(ttt1 in consider.genes){
    Fulldata <- read.csv(file.path("Circadian.gene.cluster/",Cond, paste(ttt1,"-bottom_upper_quantile_class.csv",sep="")), header = T)
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
  if(Cond == gsub("rppa_|\\.csv","",filenameall)[1]){
    Fulldata_Score_AllCancer <- Fulldata_Score_eachCancer
    ResultAllCancer <- ResulteachCancer 
  }else{
    Fulldata_Score_AllCancer <- rbind(Fulldata_Score_AllCancer,Fulldata_Score_eachCancer)
    ResultAllCancer <- rbind(ResultAllCancer,ResulteachCancer)
  }
} #################################  end of RPPA file loop

write.csv(Fulldata_Score_AllCancer,file="AllCancer_AllCir_geneExp_pathwayScore.csv",row.names = F)
write.csv( ResultAllCancer,file="AllCancer_AllCir_pathwayScore_pvalue.csv",row.names = F)
###upper and bottom quantile
Fulldata_Score_AllCance_quantile <- read.csv("AllCancer_AllCir_pathwayScore_pvalue_upperbottom_quantile.csv",header=T)
###median
Fulldata_Score_AllCancer <- read.csv("AllCancer_AllCir_pathwayScore_pvalue_median.csv",header = T)
Fulldata_Score_AllCancer <- Fulldata_Score_AllCancer[(!is.na(Fulldata_Score_AllCancer$p.value)),]
Fulldata_Score_AllCancer$type <- toupper(Fulldata_Score_AllCancer$type)
Fulldata_Score_AllCancer$Pathway <- gsub("Score","",Fulldata_Score_AllCancer$Pathway)
for( i in unique(Fulldata_Score_AllCancer$Gene)){
  cir_pathway <- Fulldata_Score_AllCancer[which(Fulldata_Score_AllCancer$Gene == i),]
  cir_pathway$FDR <- -log10(cir_pathway$FDR)
  cir_pathway["FDR"][cir_pathway["FDR"] >= 4] <- 4
  pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/",i,"_pathway FDR.pdf",sep=""),width = 6,height = 2.2)
  p <- ggplot(cir_pathway,aes(x=type,y=Pathway,fill=FDR))+
    geom_tile()+ #values=c(0,0.1,seq(0.1001,1,length.out=7))
    scale_fill_gradientn(colours=c(colorRampPalette(c("green","red"),space="rgb")(100)[c(seq(0,60,by=3),61:100)],"red","red"),
                         values=c(seq(0,-log10(0.05),length.out = 60),1.35,4),limits=c(0,4),
                         breaks=c(0,-log10(0.05),4),labels=c(1,0.05,1E-4),name="")+
    theme(panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),#line(colour="grey",linetype="dashed"),
          axis.title=element_blank(),
          axis.text.y=element_text(size=10,color="black"),
          axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
          axis.ticks=element_line(color="black"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10),
          legend.key.width = unit(0.3,"cm"),
          legend.key.heigh = unit(0.9,"cm"))#,
         # legend.key = element_rect(fill="white",colour = "black"),
        #  legend.position = "bottom",legend.direction = "horizontal")
    print(p)
    dev.off()
    print(c(i,nrow(cir_pathway[which(cir_pathway$FDR > 1.3 ),])))
}

###Example of CRY2
i="CRY2"
change_cancer <- cir.exp.pval[which(cir.exp.pval$geneSymbol==i & abs(cir.exp.pval$Fold) >= log2(1.5) & cir.exp.pval$Pvalue <= 0.05 & cir.exp.pval$Type != "ESCA"),"Type"]
change_cancer <- as.character(change_cancer)
change_cancer_FDR <- paste(change_cancer,"_FDR",sep="")
change_cancer_pval <- paste(change_cancer,"_pval",sep="")
cor_data <- read.delim(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/sign_cor/",i,"correlation_across_cancer.txt",sep=""),header = T)
change_cor_data <- cor_data[,c("gene",change_cancer,change_cancer_FDR)]
DNADamageGenes <- c("TP53BP1","ATM","BRCA2","CHEK1","CHEK2","XRCC5","MRE11A","TP53","RAD50","RAD51","XRCC1")
DNADamage.cor <- change_cor_data[which(change_cor_data$gene %in% DNADamageGenes),]
DNADamage.exp <- change_cor_data[which(cancer_ty.genes.exp.All$gene %in% DNADamageGenes),]

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



for( i in c("BRCA2","CHEK1","CHEK2","RAD51")){
  LUAD_sub <- LUAD[which(LUAD$gene %in% c("CRY2",i)),]
  
  LUAD_sub <- LUAD_sub[,filter]
  LUAD_sub.m <- data.frame(t(LUAD_sub))
  colnames(LUAD_sub.m) <- c("CRY2",i)
 # LUAD_sub.m <- LUAD_sub.m[which(LUAD_sub.m$CRY2 !=41.6069),]
  PT_cor <- signif(DNADamage.cor[which(DNADamage.cor$gene == i),"LUAD"],digits=2)
  #P_val <-  signif(cor.test(LUAD_sub.m$CRY2,LUAD_sub.m[,i])$p.value,digits = 4)
  pdf(paste("/extraspace/yye1/analysis/Circadian/RPPA/Circadian.gene.cluster_median/pathway_sig_plot/LUAD_CRY2_DNADamage/LUAD_CRY2_",i,".expression.pdf",sep=""),width=5,height=5)
  p <- ggplot(LUAD_sub.m,aes(x=log2(as.numeric(as.vector(LUAD_sub.m[,i]))),y=log2(as.numeric(as.vector(CRY2)))))+
    geom_point(color="black")+
    stat_smooth(method = "lm",color="red",size=1)+
    labs(x=paste("CRY2 (R = ",PT_cor,")",sep=""),
         y=i)+
    #annotate("text",x=-2,y=12,label=paste("n=",NB_num,sep=""),size=8,color="darkorange2")+
    # annotate("text",x=10.1,y=-2,hjust=0,label=paste("n=",n_num,sep=""),size=8,color="steelblue2")+
    theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
          panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
          axis.text=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
          axis.ticks.length = unit(.25, "cm"))
  print(p)
  dev.off()
}



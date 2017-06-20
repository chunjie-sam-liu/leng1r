#### circadian genes expression across different tumor subtype
setwd("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/")
###mRNA expression data from TCGA mRNA expression
#exp.files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
clinic.files.names <- list.files(path="/extraspace/TCGA/TCGA_clinical/",pattern = "_clinical_clean.txt")
cancer_subtype <- read.table("cancer_subtype.txt",header=T)
#commonsamples <- intersect(gsub("_clinical_clean.txt","",clinic.files.names),gsub("_mRNA_each_exp_20160513","",exp.files.names))
#as.character(cancer_subtype$category[cancer])
for(cancer in 1:nrow(cancer_subtype)){
  clinic <- read.delim(paste("/extraspace/TCGA/TCGA_clinical/",as.character(cancer_subtype$cancer[cancer]), "_clinical_clean.txt",sep=""),header=T)
  subtype1 <- clinic[!(clinic[,as.character(cancer_subtype$category[cancer])] %in% names(table(clinic[,as.character(cancer_subtype$category[cancer])]))),c("barcode",as.character(cancer_subtype$category[cancer]))]
  subtype1[,2] <- rep("Normal",times=nrow(subtype1)) 
  colnames(subtype1) <- c("barcode","subtype")
  subtype <- clinic[which(clinic[,as.character(cancer_subtype$category[cancer])] %in% names(table(clinic[,as.character(cancer_subtype$category[cancer])]))),c("barcode",as.character(cancer_subtype$category[cancer]))]
  colnames(subtype) <- c("barcode","subtype")
  subtype <- rbind(subtype,subtype1)
  exp <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",as.character(cancer_subtype$cancer[cancer]),"_mRNA_each_exp_20160513",sep=""),header=T)
  exp <- exp[which(data.frame(do.call(rbind, strsplit(as.character(exp$gene),'\\|')))$X1 %in% circadian.genes$V1),]
  commonbarcode <- intersect(clinic$barcode,gsub("\\.","\\-",substr(colnames(exp),1,12)))
  commonexpnames <- colnames(exp)[gsub("\\.","\\-",substr(colnames(exp),1,12)) %in% commonbarcode]
  commonexpnames <- commonexpnames[!duplicated(gsub("\\.","\\-",substr(commonexpnames,1,12)))]
  exp <- exp[,c("gene",commonexpnames)]
  colnames(exp) <- c("gene",gsub("\\.","\\-",substr(commonexpnames,1,12)))
  exp$GeneSymbol <- data.frame(do.call(rbind,strsplit(as.character(exp$gene),'\\|')))$X1
  exp.m <- melt(exp,id.vars = "GeneSymbol",measure.vars =gsub("\\.","\\-",substr(commonexpnames,1,12)))
  colnames(exp.m) <- c("GeneSymbol","barcode","Exp_val")
  exp.mm <- merge(exp.m,subtype,by="barcode")
  exp.mm$type <- rep(cancer_subtype$cancer[cancer],times=nrow(exp.mm))
  exp.mm$subtype <- factor(exp.mm$subtype,levels = unique(exp.mm$subtype))
  if(cancer==1){
    clock.exp.subtype <- exp.mm
  }else{
    clock.exp.subtype <- rbind(clock.exp.subtype,exp.mm)
  }
}

write.table(clock.exp.subtype,file="/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/clock.exp.in.different.subtype.txt",quote = F,row.names = F,sep="\t")
a <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/clock.exp.in.different.subtype.txt")
cancer_f <- as.character(cancer_subtype$cancer)
cancer_f <- cancer_f[!(cancer_f %in% c("READ","UCEC"))]
trend <- c()
for( i in cancer_f) {
 # if(nrow(clock.exp.subtype[which( clock.exp.subtype$type==i),]) >= 51*20){
    for(j in unique(clock.exp.subtype$GeneSymbol)){
      sub <- clock.exp.subtype[which(clock.exp.subtype$GeneSymbol==as.character(j) & clock.exp.subtype$type==i),]
      subtypeName <-  unique(sub$subtype)[!(unique(sub$subtype)=="Normal")]
      sub <- sub[which(sub$subtype != "Normal"),]
      sub$subtype <- factor(sub$subtype,levels = subtypeName)
      mean_val <- sapply(split(sub[,"Exp_val"],sub$subtype),mean)
      FC <- log2((max(mean_val)+1)/(min(mean_val)+1))
      aovPval <- oneway.test(Exp_val~subtype,data = sub)$p.value
      aovPval[is.nan(aovPval) == T] <- 1
    #  aovPval <- signif(aovPval,digits = 2)
      trend <- rbind(trend,c(i,j,FC,aovPval))
      
    }
  #}
  
}
library(varhandle)
trend <- data.frame( trend)
colnames(trend) <- c("tumor","geneSymbol","FC_log2","p_val")
trend$FDR <- signif(p.adjust(as.numeric(trend$p_val),method = "fdr"),digits = 2)
trend_sign <- trend[which(abs(as.numeric(trend$FC_log2))>= log2(1.5) & trend$FDR <= 0.05),]
write.csv(trend_sign,file="/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/significant.DE.clock.genes.in.different.subtype.csv",quote = F,row.names = F)

tumor.subtype <- unique(clock.exp.subtype[,c("barcode","subtype","type")])
tumor.subtype_sum <- sapply(split(tumor.subtype[,"subtype"],tumor.subtype$type),table)
tumor.subtype_sum <- data.frame(tumor.subtype_sum)
tumor.subtype_sum$subtype <- rownames(tumor.subtype_sum)
tumor.subtype_sum.m <- melt(tumor.subtype_sum,id.vars = "subtype",measure.vars=as.character(unique(trend$tumor)))
tumor.subtype_sum.m <- tumor.subtype_sum.m[which(tumor.subtype_sum.m$value !=0),]
tumor.subtype_sum.m <- merge(cancer_subtype[,c(1,3)],tumor.subtype_sum.m,by.y="variable",by.x="cancer")
colnames(tumor.subtype_sum.m) <- c("Cancer","Category","Subtype","Samples Number")
write.csv(tumor.subtype_sum.m,file="tumor subtype statistics.csv",quote = F,row.names = F)



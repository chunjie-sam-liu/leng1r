setwd("/extraspace/yye1/Circadian/DNA_methy")
probe <-  read.delim("/extraspace/TCGA/TCGA_methylation450k/HumanMethylation450_annotation",header=T)

colnames(probe) <- c("ProbeID","GeneSymbol","Chr","Coordinate")
cir <- read.delim("/extraspace/yye1/Circadian/circadian.genes.txt",header=F)
probe <- probe[which(probe$GeneSymbol != ""),]
probe$GeneSymbol1 <- data.frame(do.call(rbind, strsplit(as.character(probe$GeneSymbol),'\\;')))$X1
probe$GeneSymbol2 <- data.frame(do.call(rbind, strsplit(as.character(probe$GeneSymbol),'\\;')))$X2
probe.f1 <- probe[which(probe$GeneSymbol1 %in% cir$V1),]
probe.f1$GeneSymbol <- probe.f1$GeneSymbol1
probe.f2 <- probe[which((probe$GeneSymbol2 %in%  cir$V1) & (as.character(probe$GeneSymbol2) != as.character(probe$GeneSymbol1))),]
probe.f2$GeneSymbol <- probe.f2$GeneSymbol2
probe.f <- rbind(probe.f1[,1:4],probe.f2[,1:4])

ann <- read.delim("/extraspace/yye1/annotation/hg19_refGene_location.txt",header=T)

GNB2L1 <- ann[which(toupper(ann$name2) == "RACK1"),]
GNB2L1["name2"] <- "GNB2L1"
ann.cir <- ann[which(ann$name2 %in% cir$V1),]
ann.cir <- rbind(ann.cir,GNB2L1)
ann.cir <- unique(ann.cir[!(ann.cir$chrom %in% grep("_alt",ann.cir$chrom,value=T)),c(1:3,5)])

ann.cir1 <- ann.cir[which(ann.cir$strand=="+"),]
ann.cir1$pro.start <- ann.cir1$txStart - 1500
ann.cir1$pro.end <- ann.cir1$txStart + 500
ann.cir2 <- ann.cir[which(ann.cir$strand=="-"),]
ann.cir2$pro.start <- ann.cir2$txStart - 500
ann.cir2$pro.end <- ann.cir2$txStart + 1500
ann.cir3 <- rbind(ann.cir1,ann.cir2)
Merge <- merge(probe.f,ann.cir3,by.x="GeneSymbol",by.y="name2",all=T)
promoter1500.Merge <- Merge[which(as.numeric(Merge$Coordinate) >= as.numeric(Merge$pro.start) & as.numeric(Merge$Coordinate) <= as.numeric(Merge$pro.end)),]
#setdiff(as.character(unique(promoter2000.Merge$GeneSymbol)),as.character(unique(promoter.Merge$GeneSymbol)))
promoter.cir.probe <- unique(promoter1500.Merge[,c(1,2)])
write.table(promoter.cir.probe,file="promoter1000.cir.pro.txt",quote = F,row.names = F,sep="\t")

library(varhandle)
Merge$Coordinate <- unfactor(Merge$Coordinate)
promoter3000.Merge <- Merge[which(as.numeric(Merge$Coordinate) >= as.numeric(Merge$pro.start) & as.numeric(Merge$Coordinate) <= as.numeric(Merge$pro.end)),]
#setdiff(as.character(unique(promoter2000.Merge$GeneSymbol)),as.character(unique(promoter.Merge$GeneSymbol)))
promoter.cir.probe <- unique(promoter3000.Merge[,c(1,2)])
write.table(promoter.cir.probe,file="promoter.cir.pro.txt",quote = F,row.names = F,sep="\t")

###
pair.samples <- read.delim("../pair_sampleAll.txt",header=T)
pair.Methy <- promoter.cir.probe
circadian.genes <- read.delim("../circadian.genes.txt",header=F)
for( i in as.character(unique(pair.samples$type))[6:14]){
  methy <- read.delim(paste("/extraspace/TCGA/TCGA_methylation450k/",i,"_methy_probe_20160801",sep=""),header=T)
  methy_sub <- methy[which(data.frame(do.call(rbind,strsplit(as.character(methy$gene),"\\_")))$X1 %in% promoter.cir.probe$ProbeID),]
  pair <- intersect(substr(gsub("-",".",pair.samples[which(pair.samples$type == i),"barcode"]),1,15),substr(colnames(methy_sub),1,15))
  methy_sub1 <- methy_sub[,c("gene",colnames(methy_sub)[substr(colnames(methy_sub),1,15) %in% pair])]
  methy_sub1$gene <- data.frame(do.call(rbind,strsplit(as.character(methy_sub1$gene),"\\_")))$X1
  methy_sub1 <- merge(promoter.cir.probe,methy_sub1,by.x="ProbeID",by.y="gene")
  methy_sub1[is.na(methy_sub1)==T] <- 0
  t <- colnames(methy_sub1[3:ncol(methy_sub1)])[as.numeric(substr(colnames(methy_sub1[3:ncol(methy_sub1)]),14,15)) ==1]
  n <- colnames(methy_sub1[3:ncol(methy_sub1)])[as.numeric(substr(colnames(methy_sub1[3:ncol(methy_sub1)]),14,15)) >=11]
  if(length(t) >= 3 & length(n) >=3){
    DMR <- apply(methy_sub1[,t],1,mean) - apply(methy_sub1[,n],1,mean)
    DMR <- data.frame(DMR)
    colnames(DMR) <- i
    DMR$ProbeID <- methy_sub1$ProbeID
    DMR[paste(i,"_pval",sep="")] <-  apply(as.matrix(methy_sub1[3:ncol(methy_sub1)]),1,function(x){t.test(x[t],x[n])$p.value})
    pair.Methy <- merge(pair.Methy, DMR,by="ProbeID")
  }
  
}
pair.Methy$GeneSymbol <- factor(pair.Methy$GeneSymbol,levels=unique(pair.Methy$GeneSymbol))
#write.table(pair.Methy,file="pair.Methy.beta.value.FC_pval.txt",quote = F,row.names = F,sep="\t",col.names = T)
#pair.Methy <- read.delim("pair.Methy.beta.value.FC_pval.txt")
##-1000bp - +500bp
#pair.Methy <- pair.Methy[which(pair.Methy$ProbeID %in% promoter.cir.probe$ProbeID),]

filter.pair <- gsub("_pval","",grep("_pval",colnames(pair.Methy),value=T))
pair.Methy.m <- as.matrix(pair.Methy[,filter.pair])
row.names(pair.Methy.m) <- pair.Methy$ProbeID
rgb.palette <- colorRampPalette(c("blue","white","red"), space="rgb")

par(mai=c(1,1.5,1.5,2),mfcol=c(1,2)) #labRow = FALSE,
#heatmap.2(pair.Methy.m[order(pair.Methy$GeneSymbol),],Rowv = FALSE,labRow = FALSE,Colv = FALSE,dendrogram="none",col = rgb.palette(100),zlim=c(-1,1),trace="none",key=F)
colMax <- function(data) sapply(data, max, na.rm = TRUE)
pair.Methy.1 <- pair.Methy[which(as.character(pair.Methy$GeneSymbol) %in% as.character(unique(pair.Methy$GeneSymbol)[table(pair.Methy$GeneSymbol) == 1])),]
pair.Methy.1$GeneSymbol <- factor(pair.Methy.1$GeneSymbol,levels=unique(pair.Methy.1$GeneSymbol))
pair.Methy.1 <- pair.Methy.1[,c("GeneSymbol",filter.pair,paste(filter.pair,"_pval",sep=""))]
pair.Methy.2 <- pair.Methy[which(pair.Methy$GeneSymbol %in% unique(pair.Methy$GeneSymbol)[table(pair.Methy$GeneSymbol) > 1]),]
pair.Methy.2$GeneSymbol <- factor(pair.Methy.2$GeneSymbol,levels=unique(pair.Methy.2$GeneSymbol))
for( i in unique(pair.Methy.2$GeneSymbol)){
  sub <- pair.Methy.2[which(pair.Methy.2$GeneSymbol==i),]
  beta <- apply(sub[,filter.pair],2,function(x) { x[which.max( abs(x))]})
  pval <- apply(sub[,paste(filter.pair,"_pval",sep="")],2,function(x) { x[which.min( abs(x))]})
  beta_pval <- c(i,beta,pval)
  names(beta_pval) <- c("GeneSymbol",filter.pair,paste(filter.pair,"_pval",sep=""))
  beta_pval <- data.frame(t(as.matrix(beta_pval)))
  pair.Methy.1 <- rbind(pair.Methy.1,beta_pval)
  
}
cir.beta.m <- melt(pair.Methy.1,id.vars="GeneSymbol",measure.vars=filter.pair)
cir.pval.m <- melt(pair.Methy.1,id.vars="GeneSymbol",measure.vars=paste(filter.pair,"_pval",sep=""))
cir.exp.pval <-cbind(cir.beta.m,cir.pval.m)[c(1:3,6)]
colnames(cir.exp.pval) <- c("GeneSymbol","Type","BValue.Change","Pvalue")
#cir.exp.pval$Fold_label <- cir.exp.pval$Fold
cir.exp.pval$BValue.Change <- signif(as.numeric(cir.exp.pval$BValue.Change),digits = 2)
cir.exp.pval$Pvalue <- signif(as.numeric(cir.exp.pval$Pvalue),digits = 2)

cir.exp.pval <- cir.exp.pval[which((abs(as.numeric(cir.exp.pval$BValue.Change)) >= 0.15) & (as.numeric(cir.exp.pval$Pvalue) <= 0.05)),]
#cir.exp.pval$Fold_label <- as.character(signif(cir.exp.pval$Fold_label,digits = 2))
cir.exp.pval[3][abs(cir.exp.pval[3]) < 0.15] <- 0


cir.exp.pval$Pvalue <- -log10(cir.exp.pval$Pvalue)
cir.exp.pval[4][cir.exp.pval[4] > 15] <- 15
x_label_color <- read.delim("x_label_color.txt")
y_label <- read.delim("expression_y_label.txt",header=T)
#x_label_color <- x_label_color[which(x_label_color$x_label %in% unique(cir.exp.pval$GeneSymbol)),]
pdf("circadian.genes.methylation.Change.tumor.normal.morethan_10pairs_0.15_1500bp.pdf",width=12,height = 5)#,width=4000,height=1500,res=400)
ggplot(cir.exp.pval,aes(x=GeneSymbol,y=Type))+
  geom_point(aes(color=BValue.Change,size=Pvalue))+
  scale_color_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,na.value="white",breaks=seq(-0.4,0.4,length.out = 5),labels=seq(-0.4,0.4,length.out = 5),name="B-Value Change")+
  scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 6),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit=y_label$x,expand = c(0.04,0.04))+ #names(table(cir.exp.pval$Type))[order(table(cir.exp.pval$Type))]
  scale_x_discrete(limit=x_label_color$x_label,expand = c(0.03,0.2))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = x_label_color$x_color,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()



cir.exp.pval <-cbind(cir.beta.m,cir.pval.m)[c(1:3,6)]
colnames(cir.exp.pval) <- c("GeneSymbol","Type","BValue.Change","Pvalue")
#cir.exp.pval$Fold_label <- cir.exp.pval$Fold
cir.exp.pval$BValue.Change <- signif(as.numeric(cir.exp.pval$BValue.Change),digits = 2)
cir.exp.pval$Pvalue <- signif(as.numeric(cir.exp.pval$Pvalue),digits = 2)

cir.exp.pval <- cir.exp.pval[which((abs(as.numeric(cir.exp.pval$BValue.Change)) >= 0.1) & (as.numeric(cir.exp.pval$Pvalue) <= 0.05)),]
#cir.exp.pval$Fold_label <- as.character(signif(cir.exp.pval$Fold_label,digits = 2))
cir.exp.pval[3][abs(cir.exp.pval[3]) < 0.1] <- 0


cir.exp.pval$Pvalue <- -log10(cir.exp.pval$Pvalue)
cir.exp.pval[4][cir.exp.pval[4] > 15] <- 15
x_label_color <- read.delim("x_label_color.txt")
x_label_color <- x_label_color[which(x_label_color$x_label %in% unique(cir.exp.pval$GeneSymbol)),]
pdf("circadian.genes.methylation.Change.tumor.normal.morethan_10pairs_0.1.pdf",width=10,height = 5)#,width=4000,height=1500,res=400)
ggplot(cir.exp.pval,aes(x=GeneSymbol,y=Type))+
  geom_point(aes(color=BValue.Change,size=Pvalue))+
  scale_color_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,na.value="white",breaks=seq(-0.4,0.4,length.out = 5),labels=seq(-0.4,0.4,length.out = 5),name="B-Value Change")+
  scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 7),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit=names(table(cir.exp.pval$Type))[order(table(cir.exp.pval$Type))],expand = c(0.03,1))+
  scale_x_discrete(limit=x_label_color$x_label,expand = c(0.02,0))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = x_label_color$x_color,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()

####Use probe mean
pair.Methy.1 <- pair.Methy[which(as.character(pair.Methy$GeneSymbol) %in% as.character(unique(pair.Methy$GeneSymbol)[table(pair.Methy$GeneSymbol) == 1])),]
pair.Methy.1$GeneSymbol <- factor(pair.Methy.1$GeneSymbol,levels=unique(pair.Methy.1$GeneSymbol))
pair.Methy.1 <- pair.Methy.1[,c("GeneSymbol",filter.pair,paste(filter.pair,"_pval",sep=""))]
pair.Methy.2 <- pair.Methy[which(pair.Methy$GeneSymbol %in% unique(pair.Methy$GeneSymbol)[table(pair.Methy$GeneSymbol) > 1]),]
pair.Methy.2$GeneSymbol <- factor(pair.Methy.2$GeneSymbol,levels=unique(pair.Methy.2$GeneSymbol))
for( i in unique(pair.Methy.2$GeneSymbol)){
  sub <- pair.Methy.2[which(pair.Methy.2$GeneSymbol==i),]
  beta <- apply(sub[,filter.pair],2,mean)
  pval <- apply(sub[,paste(filter.pair,"_pval",sep="")],2,mean})
  beta_pval <- c(i,beta,pval)
  names(beta_pval) <- c("GeneSymbol",filter.pair,paste(filter.pair,"_pval",sep=""))
  beta_pval <- data.frame(t(as.matrix(beta_pval)))
  pair.Methy.1 <- rbind(pair.Methy.1,beta_pval)
  
}
cir.beta.m <- melt(pair.Methy.1,id.vars="GeneSymbol",measure.vars=filter.pair)
cir.pval.m <- melt(pair.Methy.1,id.vars="GeneSymbol",measure.vars=paste(filter.pair,"_pval",sep=""))
cir.exp.pval <-cbind(cir.beta.m,cir.pval.m)[c(1:3,6)]
colnames(cir.exp.pval) <- c("GeneSymbol","Type","BValue.Change","Pvalue")
#cir.exp.pval$Fold_label <- cir.exp.pval$Fold
cir.exp.pval$BValue.Change <- signif(as.numeric(cir.exp.pval$BValue.Change),digits = 2)
cir.exp.pval$Pvalue <- signif(as.numeric(cir.exp.pval$Pvalue),digits = 2)

cir.exp.pval <- cir.exp.pval[which((abs(as.numeric(cir.exp.pval$BValue.Change)) >= 0.2) & (as.numeric(cir.exp.pval$Pvalue) <= 0.05)),]
#cir.exp.pval$Fold_label <- as.character(signif(cir.exp.pval$Fold_label,digits = 2))
cir.exp.pval[3][abs(cir.exp.pval[3]) < 0.2] <- 0


cir.exp.pval$Pvalue <- -log10(cir.exp.pval$Pvalue)
cir.exp.pval[4][cir.exp.pval[4] > 15] <- 15
x_label_color <- read.delim("x_label_color.txt")
x_label_color <- x_label_color[which(x_label_color$x_label %in% unique(cir.exp.pval$GeneSymbol)),]
pdf("circadian.genes.methylation.Change.Mean_0.2.pdf",width=8,height = 5)#,width=4000,height=1500,res=400)
ggplot(cir.exp.pval,aes(x=GeneSymbol,y=Type))+
  geom_point(aes(color=BValue.Change,size=Pvalue))+
  scale_color_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,na.value="white",breaks=seq(-0.4,0.4,length.out = 5),labels=seq(-0.4,0.4,length.out = 5),name="B-Value Change")+
  scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 8),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit=names(table(cir.exp.pval$Type))[order(table(cir.exp.pval$Type))],expand = c(0.03,1))+
  scale_x_discrete(limit=x_label_color$x_label,expand = c(0.08,0))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = x_label_color$x_color,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()


###Correlation of prmoter DNA methylation and gene expression
MethyFiles <- list.files("/extraspace/TCGA/TCGA_methylation450k/",pattern="_methy_probe_20160801")
ExpFiles <- list.files("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_mRNA_each_exp_20160513")
Methy.ExpCommon <- intersect(gsub("_mRNA_each_exp_20160513","",ExpFiles),gsub("_methy_probe_20160801","",MethyFiles))
pair.samples <- read.delim("../pair_sampleAll.txt",header=T)
methyexp.corAll <- data.frame()
for( i in Methy.ExpCommon){
  methy <- read.delim(paste("/extraspace/TCGA/TCGA_methylation450k/",i,"_methy_probe_20160801",sep=""),header=T)
  methy_sub <- methy[which(data.frame(do.call(rbind,strsplit(as.character(methy$gene),"\\_")))$X1 %in% promoter.cir.probe$ProbeID),]
  methy_sub1 <- methy_sub
  methy_sub1$gene <- data.frame(do.call(rbind,strsplit(as.character(methy_sub1$gene),"\\_")))$X1
  methy_sub1 <- merge(promoter.cir.probe,methy_sub1,by.x="ProbeID",by.y="gene")
  methy_sub1[is.na(methy_sub1)==T] <- 0
  Exp <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",i,"_mRNA_each_exp_20160513",sep=""),header=T)
  Exp$gene <- data.frame(do.call(rbind,strsplit(as.character(Exp$gene),"\\|")))$X1
  Exp <-Exp[which(Exp$gene %in% helicases$V1),] ##get circadian genes
  colnames(Exp) <- substr(colnames(Exp),1,15)
  colnames(methy_sub1) <- substr(colnames(methy_sub1),1,15)
  CommonSamples <- intersect(colnames(Exp),colnames(methy_sub1))
  CommonSamples <- CommonSamples[as.numeric(substr(CommonSamples,14,15))==1]
  CommonSamples <- CommonSamples[!is.na(CommonSamples)]
  if(length(CommonSamples) >= 50){
    methy_sub1$GeneSymbol <- as.character(methy_sub1$GeneSymbol)
    Exp$gene <- as.character(Exp$gene)
    methyexp <- merge(methy_sub1[,c("ProbeID","GeneSymbol",CommonSamples)],Exp[,c("gene", CommonSamples)],by.x="GeneSymbol",by.y="gene")
    methyexp.cor <- apply(methyexp[,3:ncol(methyexp)],1,function(x){cor(x[1:length(CommonSamples)],x[(length(CommonSamples)+1) : (2*length(CommonSamples))])})
    methyexp.cor <- data.frame(methyexp.cor)
    colnames(methyexp.cor) <- i
    methyexp.cor$ProbeID <- methyexp$ProbeID
    methyexp.cor$GeneSymbol <- methyexp$GeneSymbol  
    if(nrow(methyexp.corAll)==0){
      methyexp.corAll <- methyexp.cor
    }else{
      methyexp.corAll <- merge(methyexp.corAll,methyexp.cor,by=c("ProbeID","GeneSymbol"),all=T)
    }
  }
}
methyexp.corAll[is.na(methyexp.corAll)==T] <- 0
min_function <- function(x){
  if(nrow(x) == 1){
    a <- as.numeric(x)
    names(a) <- names(x)
    
    return(a)
  }else{
    a <- as.numeric(apply(x,2,min))
    names(a) <- names(x)
    return(a)
  }
}
methyexp.corAll_select <- t(sapply(split(methyexp.corAll[,3:31],methyexp.corAll$GeneSymbol),min_function))
methyexp.corAll_select <- as.matrix(methyexp.corAll_select)

library(gplots)
library(RColorBrewer)
rgb.palette <- colorRampPalette(c("green","black","red"), space="rgb")
methyexp.corAll_select[methyexp.corAll_select > 0.5] <- 0.5
methyexp.corAll_select[methyexp.corAll_select < -0.5] <- -0.5

pdf("/extraspace/yye1/Circadian/DNA_methy/correlation of DNA methylation and circadian gene expression_0.3.pdf",width = 7,height = 7)
heatmap.2(methyexp.corAll_select,col = rgb.palette(100),breaks = c(seq(-0.5,0.1,length.out = 20),seq(-0.1,0.5,length.out = 61),seq(-0.1,0.5,length.out = 20)),
          trace="none",key=T,key.title = "none",key.ylab = "",key.xlab ="Correlation",density.info="none", keysize=1,cexCol = 1.5)
dev.off()
###1 CTRP Drug and CCLE expression data process
###1.1 CTRP Durg target gene pre-process
###1.2 Get CTRP each drug correlated genes, including correlation cofficient, fisher's Z-transformation (Z-score), FDR
###1.3  Get CTRP each drug correlated genes, including correlation cofficient, fisher's Z-transformation (Z-score), FDR
###2 GDSC Drug and CCLE expression data process
##2.1 Get GDSC each drug correlated genes, including correlation cofficient, fisher's Z-transformation (Z-score), FDR
##2.2 Calculating the threshold for Z-score for GDSC data

circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
###1 CTRP Drug and CCLE expression data process
##1.1 CCLE expression data
setwd("~/Circadian/Drug/")
CCLE_exp <- read.delim("/extraspace/yye1/share_data/CCLE/CCLE_Expression_2012-09-29.res",header=T,skip=2)
CCLE_name <- read.delim("/extraspace/yye1/share_data/CCLE/CCLE_Expression_2012-09-29.title.txt",header = T)
colnames(CCLE_exp) <- colnames(CCLE_name)
CCLE_exp <- CCLE_exp[,c(1,2,c(2:1038)*2-1)]

colMax <- function(data) sapply(data, max, na.rm = TRUE)
CCLE_exp.f <- t(sapply(split(CCLE_exp[,3:ncol(CCLE_exp)],CCLE_exp$Description),colMeans))
CCLE_exp.f<-data.frame(CCLE_exp.f)
CCLE_exp.f$Gene_Symbol <- rownames(CCLE_exp.f)
CCLE_exp.f <- CCLE_exp.f[,c(ncol(CCLE_exp.f),1:(ncol(CCLE_exp.f)-1))]
write.table(CCLE_exp.f,file="/extraspace/yye1/CCLE/CCLE_Expression_09292016.txt",quote=F,row.names = F,sep="\t")
##1.2 CTRP Durg target gene pre-process 
setwd("~/Circadian/Drug/")

CTRP_Drug_target <- read.csv("CTRP/CTRP_Drug.target.genes.csv",header=T)
CTRP_targetAll <- data.frame()
for(i in 1:nrow(CTRP_Drug_target)){
  if(length(CTRP_Drug_target[i,"gene_symbol_of_protein_target"]) > 0){
    target <-  data.frame(do.call(cbind, strsplit(as.character(CTRP_Drug_target[i,"gene_symbol_of_protein_target"]),'\\;')))
    colnames(target) <- "target.genes"
    target$Compound <- rep(CTRP_Drug_target[i,"cpd_name"],times=nrow(target))
    target$Compound_id <- rep(CTRP_Drug_target[i,"master_cpd_id"],times=nrow(target))
  }
  
  if(nrow(CTRP_targetAll)==0){
    CTRP_targetAll <- target
  }else{
    CTRP_targetAll <- rbind(CTRP_targetAll,target)
  }
}

write.table(CTRP_targetAll,file="~/Circadian/Drug/CTRP/CTRP_target.genes_uniq.txt",quote=F,row.names=F,sep="\t")
###1.3 Get CTRP each drug correlated genes, including correlation cofficient, fisher's Z-transformation (Z-score), FDR
#  CCLE name and CCLE master ccl id convertion
setwd("/home/yye1/Circadian/Drug/")
CCLE_CellLine_ann <- read.csv("/extraspace/yye1/share_data/CCLE/CCLE_CellLine_annotation.csv")
CCLE_exp <- read.delim("/extraspace/yye1/share_data//CCLE/CCLE_Expression_Entrez_2012-09-29.gct",header=T,skip=2)
CTRP.Drug <- read.csv("CTRP/CTRP_Drug.AUC.csv",header=T)
CTRP.Drug <- merge(CTRP.Drug,CCLE_CellLine_ann[,1:2],by="master_ccl_id")
CTRP.Drug.ann <- read.csv("CTRP/CTRP_Drug.name..csv",header=T)
CTRP.Drug <- merge(CTRP.Drug,CTRP.Drug.ann,by="master_cpd_id")

ccl_nameAll <- c()
for( i in 3:ncol(CCLE_exp)){
  ccl_name <- as.character(data.frame(do.call(rbind, strsplit(as.character(colnames(CCLE_exp)[i]),'\\_')))$X1)
  ccl_nameAll <- c(ccl_nameAll,ccl_name)
}
colnames(CCLE_exp) <- c("Name","Gene_Symbol",ccl_nameAll)
colnames_CCLEExp <- data.frame(colnames(CCLE_exp)[3:ncol(CCLE_exp)])
colnames(colnames_CCLEExp) <- "ccl_name"
CTRP_CCLE_CellLine <- merge(colnames_CCLEExp,CCLE_CellLine_ann[1:2],by="ccl_name")
CCLE_exp.f <- t(sapply(split(CCLE_exp[,ccl_nameAll],unique(CCLE_exp$Gene_Symbol)),colMeans))
CCLE_exp.f<-data.frame(CCLE_exp.f)
CCLE_exp.f$GeneSymbol <- rownames(CCLE_exp.f)
###
CTRP_Drug_Exp.cor <- data.frame(CCLE_exp$Gene_Symbol)
colnames(CTRP_Drug_Exp.cor) <- "Gene_Symbol"
for(Compound in  unique(CTRP.Drug$cpd_name)){
  EachDrug <- CTRP.Drug[which(CTRP.Drug$cpd_name == Compound),]
  EachDrug <- EachDrug[which(EachDrug$master_ccl_id %in% CTRP_CCLE_CellLine$master_ccl_id),]
  filter <- as.character(EachDrug$ccl_name)
  gene.cor <- apply( as.matrix(2^CCLE_exp[,sort(filter)]) , 1 , cor , y = EachDrug[order(EachDrug$ccl_name),"area_under_curve"])
  zscore <- log((1+gene.cor)/(1-gene.cor))/2*sqrt(nrow(EachDrug)-3)
  gene.cor <- data.frame(gene.cor)
  colnames(gene.cor) <- Compound
  zscore <- data.frame(zscore)
  colnames(zscore) <- paste(Compound,"_zscore",sep="")
  p_value <- apply(as.matrix(2^CCLE_exp[,sort(filter)]),1,function(x){cor.test(x, EachDrug[order(EachDrug$ccl_name),"area_under_curve"])$p.value})
  #FDR_value <- p.adjust(p_value,method="fdr")
  p_value <- data.frame(p_value)
  colnames(p_value) <- paste(Compound,"_p",sep = "")
  CTRP_Drug_Exp.cor <- cbind(CTRP_Drug_Exp.cor,gene.cor,zscore,p_value)
}


write.table(CTRP_Drug_Exp.cor,file="CTRP/CTRP.Drug_Exp.cor.txt",quote=F,row.names=F,sep="\t")
CTRP_Drug_Exp.cor <- read.delim("/home/yye1/Circadian/Drug/CTRP/CTRP.Drug_Exp.cor.txt",header=T)
colnames(CTRP_Drug_Exp.cor) <- c("Gene_Symbol",paste(rep(unique(CTRP.Drug$cpd_name),each=3),c("","_zscore","_p"),sep=""))
CTRP_Drug_Exp.cor.m <- melt(CTRP_Drug_Exp.cor,id.vars="Gene_Symbol",measure.vars=grep("_zscore",colnames(CTRP_Drug_Exp.cor)[2:ncol((CTRP_Drug_Exp.cor))],value=T))
colnames(CTRP_Drug_Exp.cor.m) <- c("Gene_Symbol","Compound","zscore")
CTRP_Drug_Exp.cor.m <- CTRP_Drug_Exp.cor.m[complete.cases(CTRP_Drug_Exp.cor.m),]
CTRP_cpd.target <- read.delim("/home/yye1/Circadian/Drug/CTRP/CTRP_target.genes_uniq.txt",header = T)
CTRP_Drug_Exp.cor.m$Compound <- gsub("_zscore","",CTRP_Drug_Exp.cor.m$Compound)
#CTRP_cpd.target$Compound_name <- gsub("\\-|\\.||\\ ","",CTRP_cpd.target$Compound)
CTRP_cpd.target.cor <- merge(CTRP_Drug_Exp.cor.m,CTRP_cpd.target,by.x=c("Gene_Symbol","Compound"),by.y=c("target.genes","Compound"))
set.seed(1)
random_zscore <- sample(CTRP_Drug_Exp.cor.m$zscore,nrow(CTRP_cpd.target.cor))
CTRP_cpd.target.cor$random <- random_zscore
#threshold <- CTRP_cpd.target.cor[order(abs(CTRP_cpd.target.cor$random)),][as.integer(0.95*nrow(CTRP_cpd.target.cor)),"random"]
x <- c(abs(random_zscore),abs(CTRP_cpd.target.cor$zscore)) #
p <- pnorm(x,lower.tail = F)
bonftest <- p > 0.025/length(x)
summary(bonftest[1:length(random_zscore)])
summary(bonftest[(length(random_zscore)):(2*length(random_zscore))])
threshold_CTRP <- min(c(abs(random_zscore),abs(CTRP_cpd.target.cor$zscore))[!bonftest]) #
pdf("CTRP/CTRP_Drug.target_random_correlation1.pdf",width=8,height=8)
par(mar=c(5,5,1,1))
d <- density(CTRP_cpd.target.cor$zscore) # returns the density data 
random <- density(CTRP_cpd.target.cor$random) 
plot(range(d$x,random$x),range(d$y,random$y),ylab="Probability density estimate",xlab = "z-scored correlation strength",main="",xlim=c(-22,22),type = "n",cex.lab=2,cex.axis=1.5) # plots the results
lines(d,col="darkgreen",lwd=3)
lines(random,col="black",lwd=3)
abline(v=c(-threshold_CTRP,threshold_CTRP),lty="dashed",lwd=2)

dev.off()
Sign.CTRP_Drug_Exp.cor <- CTRP_Drug_Exp.cor.m[which(abs(CTRP_Drug_Exp.cor.m$zscore) >= threshold_CTRP),]
Sign.CTRP_Drug_Exp.cor <- Sign.CTRP_Drug_Exp.cor[,c(1,3,4)]

Circadian_Sign.CTRP_Drug_Exp.cor <- Sign.CTRP_Drug_Exp.cor[which(Sign.CTRP_Drug_Exp.cor$Gene_Symbol %in% circadian.genes$V1),]
###
CTRP_Drug_Exp.cor.p <- melt(CTRP_Drug_Exp.cor,id.vars="Gene_Symbol",measure.vars=grep("_p",colnames(CTRP_Drug_Exp.cor)[2:ncol((CTRP_Drug_Exp.cor))],value=T))
colnames(CTRP_Drug_Exp.cor.p) <- c("Gene_Symbol","Compound","p")
CTRP_Drug_Exp.cor.p$Compound <- gsub("_p","",CTRP_Drug_Exp.cor.p$Compound)
Circadian_Sign.CTRP_Drug_Exp.cor.p <- merge(Circadian_Sign.CTRP_Drug_Exp.cor,CTRP_Drug_Exp.cor.p,by=c("Gene_Symbol","Compound"))
Circadian_Sign.CTRP_Drug_Exp.cor.m <- Circadian_Sign.CTRP_Drug_Exp.cor.p[which(Circadian_Sign.CTRP_Drug_Exp.cor.p$p < 0.05),]
write.table(Circadian_Sign.CTRP_Drug_Exp.cor.m,file="CTRP/51Circadian_Sign.CTRP_Drug_Exp.cor.txt",quote = F,sep="\t",row.names = F)
Circadian_Sign.CTRP_Drug_Exp.cor.m <- read.delim("CTRP/51Circadian_Sign.CTRP_Drug_Exp.cor.txt",header=T)
Circadian.Drug.cor.num.pos <- table(Circadian_Sign.CTRP_Drug_Exp.cor.m[which(Circadian_Sign.CTRP_Drug_Exp.cor.m$zscore > 0),]$Gene_Symbol)[table(Circadian_Sign.CTRP_Drug_Exp.cor.m[which(Circadian_Sign.CTRP_Drug_Exp.cor.m$zscore > 0),]$Gene_Symbol)>=1]
Circadian.Drug.cor.num.pos <- Circadian.Drug.cor.num.pos[rev(order(Circadian.Drug.cor.num.pos))]
Circadian.Drug.cor.num.pos <- data.frame(Circadian.Drug.cor.num.pos)
Circadian.Drug.cor.num.pos$GeneSymbol <- rownames(Circadian.Drug.cor.num.pos)
colnames(Circadian.Drug.cor.num.pos) <- c("Pos","GeneSymbol")
Circadian.Drug.cor.num.neg <- table(Circadian_Sign.CTRP_Drug_Exp.cor.m[which(Circadian_Sign.CTRP_Drug_Exp.cor.m$zscore < 0),]$Gene_Symbol)[table(Circadian_Sign.CTRP_Drug_Exp.cor.m[which(Circadian_Sign.CTRP_Drug_Exp.cor.m$zscore < 0),]$Gene_Symbol)>=1]
Circadian.Drug.cor.num.neg <- Circadian.Drug.cor.num.neg[rev(order(Circadian.Drug.cor.num.neg))]
Circadian.Drug.cor.num.neg <- data.frame(Circadian.Drug.cor.num.neg)
Circadian.Drug.cor.num.neg$GeneSymbol <- rownames(Circadian.Drug.cor.num.neg)
colnames(Circadian.Drug.cor.num.neg) <- c("Neg","GeneSymbol")
Circadian.Drug.cor.num.both <- merge(Circadian.Drug.cor.num.pos,Circadian.Drug.cor.num.neg,by="GeneSymbol",all=T)
Circadian.Drug.cor.num.both[is.na(Circadian.Drug.cor.num.both)] <- 0
Circadian.Drug.cor.num.both.m <- melt(Circadian.Drug.cor.num.both,id.vars = "GeneSymbol",measure.vars = c("Pos","Neg"))
Num_xlimit <- Circadian.Drug.cor.num.both[order(Circadian.Drug.cor.num.both$Pos - Circadian.Drug.cor.num.both$Neg),"GeneSymbol"]
xcolor <- rep("black",times=length(Num_xlimit))
xcolor[Num_xlimit %in% core.circadian] <- "red"

####
pdf("CTRP/51Circadian.Drug.cor.num1.pdf",width=6,height=3.2)
ggplot(Circadian.Drug.cor.num.both.m,aes(x=GeneSymbol,y=value,fill=factor(variable)))+
  geom_bar(color=NA,width = 0.8,stat = "identity")+
  scale_y_continuous(expand=c(0,0),limit=c(-3,(max(Circadian.Drug.cor.num.both[,2:3])+10)))+
  scale_x_discrete(limits=  Num_xlimit)+
  scale_fill_manual(limits=c("Pos","Neg"),values = c("red","blue"),labels=c("Positive-Cor","Negative-Cor"))+
  ylab("Compounds Count")+
  theme(axis.text.x=element_text(size=10,colour =xcolor,angle=90,hjust=1,vjust=0.5),
        panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor=element_line(colour=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10,color="black",angle=90,vjust=0.5,hjust=0.5),
        legend.title=element_blank(),legend.text=element_text(size=14),
        legend.position="bottom",legend.direction="horizontal")#+
dev.off()
###Compound effect at least 10 Circadian genes 
GeneSymbol <- as.character(circadian.genes$V1)
CTRP.cpd.f <- table(Circadian_Sign.CTRP_Drug_Exp.cor.m$Compound)[table(Circadian_Sign.CTRP_Drug_Exp.cor.m$Compound) >=10]
CTRP.cpd.f <- CTRP.cpd.f[rev(order(CTRP.cpd.f))]
Circadian_Sign.CTRP_Drug_Exp.cor.mm <- Circadian_Sign.CTRP_Drug_Exp.cor.m[which(Circadian_Sign.CTRP_Drug_Exp.cor.m$Compound %in% names(CTRP.cpd.f)),]
Pvalue <- -log10(Circadian_Sign.CTRP_Drug_Exp.cor.mm$p)
Pvalue[is.infinite(Pvalue)] <- 15
Pvalue[Pvalue>=15] <- 15
Circadian_Sign.CTRP_Drug_Exp.cor.mm$Pvalue <- Pvalue
CTRP.cir.f <- CTRP_Drug_Exp.cor[which(CTRP_Drug_Exp.cor$Gene_Symbol %in% GeneSymbol),c("Gene_Symbol",paste(names(CTRP.cpd.f),"_zscore",sep=""))]
CTRP.cir.f[2:ncol(CTRP.cir.f)][CTRP.cir.f[2:ncol(CTRP.cir.f)] < threshold & CTRP.cir.f[2:ncol(CTRP.cir.f)] > -threshold] <- 0
CTRP.cir.f[2:ncol(CTRP.cir.f)][CTRP.cir.f[2:ncol(CTRP.cir.f)] >= threshold] <- 1
CTRP.cir.f[2:ncol(CTRP.cir.f)][CTRP.cir.f[2:ncol(CTRP.cir.f)] <= -threshold] <- -1
CTRP_ylimit <- CTRP.cir.f[order(apply(CTRP.cir.f[2:ncol(CTRP.cir.f)],1,sum)),'Gene_Symbol']
CTRP_ycolor <- rep("black",times=length(CTRP_ylimit ))
CTRP_ycolor[CTRP_ylimit %in% core.circadian] <- "red"
Circadian_Sign.CTRP_Drug_Exp.cor.mm["zscore"][Circadian_Sign.CTRP_Drug_Exp.cor.mm["zscore"] > 11] <- 11
Circadian_Sign.CTRP_Drug_Exp.cor.mm["zscore"][Circadian_Sign.CTRP_Drug_Exp.cor.mm["zscore"] < -11] <- -11

pdf("CTRP/51Circadian_Sign.CTRP_Drug_Exp.cor.morethan_10cpds1.pdf",width=14,height = 8)#,width=4000,height=1500,res=400)
ggplot(Circadian_Sign.CTRP_Drug_Exp.cor.mm,aes(x=Compound,y=Gene_Symbol))+
  geom_point(aes(size=Pvalue,col=zscore))+
  scale_color_gradientn(colours=colorRampPalette(c("green2","white","magenta"),space="rgb")(100))+#,
  scale_size_continuous(limit=c(-log10(0.05),15),range = c(0.5, 2.5),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit=CTRP_ylimit,expand = c(0.012,0.012))+
  scale_x_discrete(limit=names(CTRP.cpd.f),expand = c(0.005,0.005))+
  theme(panel.background=element_rect(colour="black",fill="white",size=0.1),
       # panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.05),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour =CTRP_ycolor),
        axis.text.x=element_text(size=7,colour = "black",angle=90,hjust=1,vjust=0.5,face = "bold"),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(fill="white",colour = "black"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()
# 

CTRP.cpd.fm <- data.frame(CTRP.cpd.f)
CTRP.cpd.fm$DRUG.NAME <- rownames(CTRP.cpd.fm)
CTRP_screen_cpd <- read.csv("/home/yye1/Circadian/Drug/CTRP/CTRP_Drug.target.genes.pathway.csv",header=T)
#CTRP_screen_cpd$cpd_name <- gsub(" |-|,",".",CTRP_screen_cpd$cpd_name)
CTRP_screen_cpd.m <- merge(CTRP.cpd.fm,CTRP_screen_cpd,by.x="DRUG.NAME",by.y="cpd_name")
CTRP_screen_cpd.m <- CTRP_screen_cpd.m[rev(order(CTRP_screen_cpd.m$CTRP.cpd.f)),]
#CTRP_screen_cpd.m <- unique(CTRP_screen_cpd.m[which(CTRP_screen_cpd.m$target_or_activity_of_compound != "other"),c("DRUG.NAME","target_or_activity_of_compound","CTRP.cpd.f")])

CTRP_DRUG_circadian_pathway_Count <- table(CTRP_screen_cpd.m$target_or_activity_of_compound)
CTRP_DRUG_circadian_pathway_Count <- data.frame(CTRP_DRUG_circadian_pathway_Count)

#Circadian_Sign.CTRP_Drug_Exp.cor.mmm <- merge(Circadian_Sign.CTRP_Drug_Exp.cor.mm,CTRP_screen_cpd[,c(1,7)],by.x="Compound",by.y="cpd_name")
#write.table(Circadian_Sign.CTRP_Drug_Exp.cor.mmm ,file = "/home/yye1/Circadian/Drug/CTRP/51Circadian_Sign.CTRP_Drug_Exp.cor.mmm.txt",quote = F,row.names = F,sep="\t")
Circadian_Sign.CTRP_Drug_Exp.cor.mmm <- read.delim("/home/yye1/Circadian/Drug/CTRP/51Circadian_Sign.CTRP_Drug_Exp.cor.mmm.txt",header = T)
#write.csv(unique(Circadian_Sign.CTRP_Drug_Exp.cor.mmm[,c(1,6)]),file="/home/yye1/Circadian/Drug/CTRP/51correlated.cpd.target.csv",quote = F,row.names = F)

correlated.cpd.target <- read.csv("/home/yye1/Circadian/Drug/CTRP/51correlated.cpd.target.csv")
Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_name <- gsub(" |-|,",".",Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound)
Circadian_Sign.CTRP_Drug_Exp.cor.mmm <- merge(Circadian_Sign.CTRP_Drug_Exp.cor.mmm[c(1:5,7)],correlated.cpd.target,by="Compound_name")
#pathway_order <- c("metabolism","chromain  histone acetylation","cell cycle","PI3K signaling","mitosis","chromatin  histone methylation","ERK MAPK signaling","DNA replication","apoptosis regulation","TOR signaling","cytoskeleton","chromatin  other","WNT signaling","Genome integrity","RTK signaling","JNK and p38 signaling","EGFR signaling","other")
Circadian_Sign.CTRP_Drug_Exp.cor.mmm$order <- rep(0,times=nrow(Circadian_Sign.CTRP_Drug_Exp.cor.mmm))
pathway_order <- data.frame(table(correlated.cpd.target$target_or_activity_of_compound))
pathway_order <- pathway_order[rev(order(pathway_order$Freq)),]
for( i in 1:length(pathway_order$Var1)){
  Circadian_Sign.CTRP_Drug_Exp.cor.mmm["order"][Circadian_Sign.CTRP_Drug_Exp.cor.mmm["target_or_activity_of_compound"]== as.character(pathway_order$Var1[i])] <- i
}

Circadian_Sign.CTRP_Drug_Exp.cor.mmm$cor <- rep(0,times=nrow(Circadian_Sign.CTRP_Drug_Exp.cor.mmm))
Circadian_Sign.CTRP_Drug_Exp.cor.mmm["cor"][Circadian_Sign.CTRP_Drug_Exp.cor.mmm["zscore"] >= threshold] <- 1
Circadian_Sign.CTRP_Drug_Exp.cor.mmm["cor"][Circadian_Sign.CTRP_Drug_Exp.cor.mmm["zscore"] <= -threshold] <- -1
Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Gene_Symbol <- factor(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Gene_Symbol,levels=unique(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Gene_Symbol))
CTRP_x <- sapply(split(Circadian_Sign.CTRP_Drug_Exp.cor.mmm[,"cor"],Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Gene_Symbol),sum)
CTRP_x <- names(CTRP_x[order(CTRP_x)])
CTRP_y <- unique(Circadian_Sign.CTRP_Drug_Exp.cor.mmm[rev(order(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$order)),1])
CTRP_yy <- unique(Circadian_Sign.CTRP_Drug_Exp.cor.mmm[rev(order(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$order)),c(1,6,7)])
CTRP_xcolor <- rep("black",times=length(CTRP_x))
CTRP_xcolor[CTRP_x %in% core.circadian] <- "red"
Circadian_Sign.CTRP_Drug_Exp.cor.mmm["zscore"][Circadian_Sign.CTRP_Drug_Exp.cor.mmm["zscore"] > 11] <- 11
Circadian_Sign.CTRP_Drug_Exp.cor.mmm["zscore"][Circadian_Sign.CTRP_Drug_Exp.cor.mmm["zscore"] < -11] <- -11

library(RColorBrewer)
CTRP_ycol <- rep(brewer.pal(9, "Set1")[c(1:5,7:9)],times=16)[CTRP_yy$order]
pdf("CTRP/51Circadian_Sign.CTRP_Drug_Exp.cor.morethan_10cpds_pathway_ycolors_2.pdf",width=12,height = 8)#,width=4000,height=1500,res=400)
ggplot(Circadian_Sign.CTRP_Drug_Exp.cor.mmm[order(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$order),],aes(x=Compound_name,y=Gene_Symbol))+
  geom_point(aes(size=Pvalue,col=zscore))+
  scale_color_gradientn(limit=c(-11,11),colours=colorRampPalette(c("green2","white","magenta"),space="rgb")(100))+#,
  # values=c(seq(min(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$zscore),-3,length.out = 10),seq(-3,3,length.out = 81),seq(3,max(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$zscore),length.out = 10)))+
 scale_size_continuous(limit=c(-log10(0.05),15),range = c(0.2, 4),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_x_discrete(limit=CTRP_yy$Compound_name,expand = c(0.01,0.01))+
  scale_y_discrete(limit=CTRP_x,expand = c(0.02,0.02))+
  theme(panel.background=element_rect(colour="black",fill="white",size=0.2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size = 0.1),
        axis.title=element_blank(),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_text(size=10,colour = CTRP_xcolor),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key = element_rect(fill="white",colour = "black"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()

CTRP_yy1 <- data.frame(table(CTRP_yy$target_or_activity_of_compound))
CTRP_yy2 <- merge(CTRP_yy1,unique(Circadian_Sign.CTRP_Drug_Exp.cor.mmm[,c(6,7)]),by.x="Var1",by.y="target_or_activity_of_compound")
CTRP_yy2 <- CTRP_yy2[order(CTRP_yy2$order),]

CTRP_yy2$order <- factor(CTRP_yy2$order,levels=121:1)
plus.vector<-c()
plus<-0
plus1.vector<-c()
for(i in rev(CTRP_yy2$Freq)){
  plus1<-plus+as.integer(i*0.5)+0.5
  plus<-plus+i
  plus1.vector<-c(plus1.vector,plus1)
  plus.vector<-c(plus.vector,plus)
}
pdf("CTRP/51_CTRP_bar.pdf",width=6.3,height=18)
ggplot(CTRP_yy2,aes(x=1,y=rev(Freq),fill=factor(order)))+
  geom_bar(stat="identity",position="stack",width=0.001)+
  scale_fill_manual(limits=seq(1,121),values= rev(rep(brewer.pal(9, "Set1")[c(1:5,7:9)],times=16)[1:121]),guide=F)+
  scale_y_continuous(breaks=plus1.vector,labels=rev(CTRP_yy2$Var1))+
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(colour=NA),
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=10,colour =rev(rep(brewer.pal(9, "Set1")[c(1:5,7:9)],times=16)[1:121])),
        axis.ticks=element_blank(),
        legend.title=element_blank(),legend.text=element_text(size=16))
dev.off()



#CTRP_cpd.target.cor[c("zscore","random")] <- abs(CTRP_cpd.target.cor[,c("zscore","random")])
#CTRP_cpd.target.cor.m <- melt(CTRP_cpd.target.cor,id.vars = "Gene_Symbol",measure.vars = c("zscore","random"))


###2 GDSC Drug and CCLE expression data process
##2.1 Get GDSC each drug correlated genes, including correlation cofficient, fisher's Z-transformation, p
GDSC.exp <- read.delim("/home/yye1/Circadian/Drug/GDSC/Cell_line_RMA_proc_basalExp.txt",header=T)
GDSC.Drug <- read.csv("/home/yye1/Circadian/Drug/GDSC/v17_fitted_dose_response.csv",header=T)
DB_conversion <- read.csv("/home/yye1/Circadian/Drug/GDSC/GDSC_CCLE_CTRP_conversion.csv",header=T)
GDSC.compound <- read.csv("/home/yye1/Circadian/Drug/GDSC/GDSC_screen_compound.csv",header=T)
GDSC.Drug <- merge(GDSC.Drug,GDSC.compound[,1:2],by.y="DRUG.ID",by.x="DRUG_ID")
GDSC.Drug$COSMIC_ID <- paste("DATA.",GDSC.Drug$COSMIC_ID,sep="")

GDSC_Drug_Exp.cor <- data.frame(GDSC.exp$GENE_SYMBOLS)
colnames(GDSC_Drug_Exp.cor) <- "Gene_Symbol"
for( Compound in unique(GDSC.Drug$DRUG.NAME)){
  EachDrug <-GDSC.Drug[which(GDSC.Drug$DRUG.NAME == Compound & (GDSC.Drug$COSMIC_ID %in% colnames(GDSC.exp))),]
  if(length(EachDrug$COSMIC_ID) > length(unique(EachDrug$COSMIC_ID))){
    EachDrug<- t(sapply(split(EachDrug[,c("LN_IC50","AUC")],EachDrug$COSMIC_ID),colMeans))
    EachDrug<-data.frame(EachDrug)
    EachDrug$COSMIC_ID <- rownames(EachDrug)
  }
  gene.cor <- apply( as.matrix(GDSC.exp[, EachDrug$COSMIC_ID]) , 1 , cor , y = EachDrug$AUC) 
  zscore <- log((1+gene.cor)/(1-gene.cor))/2*sqrt(nrow(EachDrug)-3)
  gene.cor <- data.frame(gene.cor)
  colnames(gene.cor) <- Compound
  zscore <- data.frame(zscore)
  colnames(zscore) <- paste(Compound,"_zscore",sep="")
  p_value <- apply(as.matrix(GDSC.exp[, EachDrug$COSMIC_ID]),1,function(x){cor.test(x, EachDrug$AUC)$p.value})
  #FDR_value <- p.adjust(p_value,method="fdr")
  p_value <- data.frame(p_value)
  colnames(p_value) <- paste(Compound,"_p",sep = "")
  GDSC_Drug_Exp.cor <- cbind(GDSC_Drug_Exp.cor,gene.cor,zscore,p_value)
}


write.table(GDSC_Drug_Exp.cor,file="/home/yye1/Circadian/Drug/GDSC/GDSC_Drug_Exp_correlation.txt",sep="\t",row.names=F,col.names=T)
##2.2 Calculating the threshold for Z-score 
GDSC_Drug_Exp.cor <- read.delim("/home/yye1/Circadian/Drug/GDSC/GDSC_Drug_Exp_correlation.txt",header=T)
GDSC_Drug_Exp.cor.m <- melt(GDSC_Drug_Exp.cor,id.vars="Gene_Symbol",measure.vars=grep("_zscore",colnames(GDSC_Drug_Exp.cor)[2:ncol((GDSC_Drug_Exp.cor))],value=T))
colnames(GDSC_Drug_Exp.cor.m) <- c("Gene_Symbol","Compound","zscore")
GDSC_Drug_Exp.cor.m <- GDSC_Drug_Exp.cor.m[complete.cases(GDSC_Drug_Exp.cor.m),]
GDSC_cpd.target <- read.csv("GDSC/GDSC.Drug.target.genes.csv",header = T)
GDSC_cpd.target$TARGET <- gsub("\\ ","",GDSC_cpd.target$TARGET)
GDSC_cpd.target$Compound_name <- gsub("\\-|\\ |\\(|\\)|\\/","\\.",GDSC_cpd.target$DRUG.NAME)
GDSC_cpd.target$Compound_name <- gsub("17\\.AAG","X17\\.AAG",GDSC_cpd.target$Compound_name)
GDSC_cpd.target$Compound_name <- gsub("681640","X681640",GDSC_cpd.target$Compound_name)
GDSC_cpd.target$Compound_name <- gsub("\\.5Z\\.\\.7\\.Oxozeaenol","X\\.5Z\\.\\.7\\.Oxozeaenol",GDSC_cpd.target$Compound_name)
GDSC_cpd.target$Compound_name <- gsub("5\\.Fluorouracil","X5\\.Fluorouracil",GDSC_cpd.target$Compound_name)

GDSC_Drug_Exp.cor.m$Compound_name <- gsub("_zscore","",GDSC_Drug_Exp.cor.m$Compound)

GDSC_cpd.target.cor <- merge(GDSC_Drug_Exp.cor.m,GDSC_cpd.target,by.x=c("Gene_Symbol","Compound_name"),by.y=c("TARGET","Compound_name"))
set.seed(1001)
random_zscore <- sample(GDSC_Drug_Exp.cor.m$zscore,nrow(GDSC_cpd.target.cor))
GDSC_cpd.target.cor$random <- random_zscore
x <- c(abs(random_zscore),abs(GDSC_cpd.target.cor$zscore))
p <- pnorm(x,lower.tail = F)
bonftest <- p > 0.025/length(x)
summary(bonftest[1:length(random_zscore)])
summary(bonftest[(length(random_zscore)):(2*length(random_zscore))])
#threshold <- min(c(abs(random_zscore),abs(GDSC_cpd.target.cor$zscore))[!bonftest])
a <- c(random_zscore,GDSC_cpd.target.cor$zscore)[!bonftest]
threshold_GDSC <- signif(max(abs(max(a[a<0])),min(a[a<0])),digits = 3)
pdf("GDSC/GDSC_Drug.target_random_correlation.pdf",width=8,height=8)
par(mar=c(5,5,1,1))
d <- density(GDSC_cpd.target.cor$zscore) # returns the density data 
random <- density(GDSC_cpd.target.cor$random) 
plot(range(d$x,random$x),range(d$y,random$y),ylab="Probalility density estimate",xlab = "z-scored correlation strength",main="",xlim=c(-22,22),type = "n",cex.lab=2,cex.axis=1.5) # plots the results
lines(d,col="darkgreen",lwd=3)
lines(random,col="black",lwd=3)
abline(v=c(-threshold_GDSC,threshold_GDSC),lty="dashed",lwd=2)
dev.off()
Sign.GDSC_Drug_Exp.cor <- GDSC_Drug_Exp.cor.m[which(abs(GDSC_Drug_Exp.cor.m$zscore) > threshold_GDSC),]
Sign.GDSC_Drug_Exp.cor <- Sign.GDSC_Drug_Exp.cor[,c(1,3,4)]
Circadian_Sign.GDSC_Drug_Exp.cor <- Sign.GDSC_Drug_Exp.cor[which(Sign.GDSC_Drug_Exp.cor$Gene_Symbol %in% circadian.genes),]
write.table(Circadian_Sign.GDSC_Drug_Exp.cor,file="GDSC/Circadian_Sign.GDSC_Drug_Exp.cor.txt",quote=F,row.names = F,col.names = T,sep="\t")
a <- read.delim("GDSC/Circadian_Sign.GDSC_Drug_Exp.cor.txt",header=T)
###
GDSC_Drug_Exp.cor.p <- melt(GDSC_Drug_Exp.cor,id.vars="Gene_Symbol",measure.vars=grep("_p",colnames(GDSC_Drug_Exp.cor)[2:ncol((GDSC_Drug_Exp.cor))],value=T))
colnames(GDSC_Drug_Exp.cor.p) <- c("Gene_Symbol","Compound_name","p")
GDSC_Drug_Exp.cor.p$Compound_name <- gsub("_p","",GDSC_Drug_Exp.cor.p$Compound_name)
Circadian_Sign.GDSC_Drug_Exp.cor.p <- merge(Circadian_Sign.GDSC_Drug_Exp.cor[,c(1,3)],GDSC_Drug_Exp.cor.p,by=c("Gene_Symbol","Compound_name"))
Pvalue <- -log10(Circadian_Sign.GDSC_Drug_Exp.cor.p$p)
Pvalue[is.infinite(Pvalue)] <- 30
Pvalue[Pvalue>=30] <- 30
Circadian_Sign.GDSC_Drug_Exp.cor.p$pvalue <- Pvalue
Circadian_Sign.GDSC_Drug_Exp.cor.m <-merge(Circadian_Sign.GDSC_Drug_Exp.cor,Circadian_Sign.GDSC_Drug_Exp.cor.p,by=c("Gene_Symbol","Compound_name"))
write.table(Circadian_Sign.GDSC_Drug_Exp.cor.m,file="GDSC/51Circadian_Sign.GDSC_Drug_Exp.cor.txt",quote = F,sep="\t",row.names = F)
Circadian_Sign.GDSC_Drug_Exp.cor.m <- read.delim("GDSC/51Circadian_Sign.GDSC_Drug_Exp.cor.txt",header=T)
GDSC.Circadian.Drug.cor.num.pos <- table(Circadian_Sign.GDSC_Drug_Exp.cor.m[which(Circadian_Sign.GDSC_Drug_Exp.cor.m$zscore > 0),]$Gene_Symbol)[table(Circadian_Sign.GDSC_Drug_Exp.cor.m[which(Circadian_Sign.GDSC_Drug_Exp.cor.m$zscore > 0),]$Gene_Symbol)>=1]
GDSC.Circadian.Drug.cor.num.pos <- GDSC.Circadian.Drug.cor.num.pos[rev(order(GDSC.Circadian.Drug.cor.num.pos))]
GDSC.Circadian.Drug.cor.num.pos <- data.frame(GDSC.Circadian.Drug.cor.num.pos)
GDSC.Circadian.Drug.cor.num.pos$GeneSymbol <- rownames(GDSC.Circadian.Drug.cor.num.pos)
colnames(GDSC.Circadian.Drug.cor.num.pos) <- c("Pos","GeneSymbol")
GDSC.Circadian.Drug.cor.num.neg <- table(Circadian_Sign.GDSC_Drug_Exp.cor.m[which(Circadian_Sign.GDSC_Drug_Exp.cor.m$zscore < 0),]$Gene_Symbol)[table(Circadian_Sign.GDSC_Drug_Exp.cor.m[which(Circadian_Sign.GDSC_Drug_Exp.cor.m$zscore < 0),]$Gene_Symbol)>=1]
GDSC.Circadian.Drug.cor.num.neg <- GDSC.Circadian.Drug.cor.num.neg[rev(order(GDSC.Circadian.Drug.cor.num.neg))]
GDSC.Circadian.Drug.cor.num.neg <- data.frame(GDSC.Circadian.Drug.cor.num.neg)
GDSC.Circadian.Drug.cor.num.neg$GeneSymbol <- rownames(GDSC.Circadian.Drug.cor.num.neg)
colnames(GDSC.Circadian.Drug.cor.num.neg) <- c("Neg","GeneSymbol")
GDSC.Circadian.Drug.cor.num.both <- merge(GDSC.Circadian.Drug.cor.num.pos,GDSC.Circadian.Drug.cor.num.neg,by="GeneSymbol",all=T)
GDSC.Circadian.Drug.cor.num.both[is.na(GDSC.Circadian.Drug.cor.num.both)] <- 0
GDSC.Circadian.Drug.cor.num.both.m <- melt(GDSC.Circadian.Drug.cor.num.both,id.vars = "GeneSymbol",measure.vars = c("Pos","Neg"))
#GDSC.Num_xlimit <- GDSC.Circadian.Drug.cor.num.both[order(GDSC.Circadian.Drug.cor.num.both$Pos - GDSC.Circadian.Drug.cor.num.both$Neg),"GeneSymbol"]
N_than_P <- GDSC.Circadian.Drug.cor.num.both[which(GDSC.Circadian.Drug.cor.num.both$Neg>GDSC.Circadian.Drug.cor.num.both$Pos),]
P_than_N <- GDSC.Circadian.Drug.cor.num.both[which(GDSC.Circadian.Drug.cor.num.both$Pos>GDSC.Circadian.Drug.cor.num.both$Neg),]
GDSC.Num_xlimit <- rbind(N_than_P[rev(order(apply(N_than_P[,2:3],1,sum))),],P_than_N[order(apply(P_than_N[,2:3],1,sum)),])$GeneSymbol
GDSC.Num_xcol <- rep("black",times=length(GDSC.Num_xlimit))
GDSC.Num_xcol[GDSC.Num_xlimit %in% core.circadian] <- "red"
pdf("GDSC/51_GDSC.Circadian.Drug.cor.num1.pdf",width=6,height=3.2)
ggplot(GDSC.Circadian.Drug.cor.num.both.m,aes(x=GeneSymbol,y=value,fill=factor(variable)))+
  geom_bar(color=NA,width = 0.8,stat = "identity")+
  scale_y_continuous(expand=c(0,0),limit=c(-1,(max(apply(GDSC.Circadian.Drug.cor.num.both[,2:3],1,sum))+10)))+
  scale_x_discrete(limits=GDSC.Num_xlimit)+
  scale_fill_manual(limits=c("Pos","Neg"),values = c("red","blue"),labels=c("Positive Cor","Negative Cor"))+
  ylab("Compounds Count")+
  theme(axis.text.x=element_text(size=10,colour = GDSC.Num_xcol,angle=90,hjust=1,vjust=0.5),
        panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor=element_line(colour=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=10,color="black",angle=90,vjust=0.5,hjust=0.5),
        legend.title=element_blank(),legend.text=element_text(size=14),
        legend.position="bottom",legend.direction="horizontal")#+
dev.off()
###Compound effect at least 10 Circadian genes 
GeneSymbol <- as.character(circadian.genes)
#GeneSymbol <- GeneSymbol[!(GeneSymbol %in% c("FOXO1","FOXM1","FOXO3","FOXA2"))]
cpd.f <- table(Circadian_Sign.GDSC_Drug_Exp.cor.m$Compound_name)[table(Circadian_Sign.GDSC_Drug_Exp.cor.m$Compound_name) >=10]
cpd.f <- cpd.f[rev(order(cpd.f))]
Circadian_Sign.GDSC_Drug_Exp.cor.mm <- Circadian_Sign.GDSC_Drug_Exp.cor.m[which(Circadian_Sign.GDSC_Drug_Exp.cor.m$Compound_name %in% names(cpd.f)),]
#Circadian_Sign.GDSC_Drug_Exp.cor.mm <- Circadian_Sign.GDSC_Drug_Exp.cor.mm[which(Circadian_Sign.GDSC_Drug_Exp.cor.mm$Gene_Symbol),]
Circadian_Sign.GDSC_Drug_Exp.cor.mm$Gene_Symbol <- factor(Circadian_Sign.GDSC_Drug_Exp.cor.mm$Gene_Symbol,levels=GeneSymbol)
cir.f <- GDSC_Drug_Exp.cor[which(GDSC_Drug_Exp.cor$Gene_Symbol %in% GeneSymbol),c("Gene_Symbol",paste(names(cpd.f),"_zscore",sep=""))]
cir.f[2:ncol(cir.f)][cir.f[2:ncol(cir.f)] < threshold & cir.f[2:ncol(cir.f)] > -threshold] <- 0
cir.f[2:ncol(cir.f)][cir.f[2:ncol(cir.f)] >= threshold] <- 1
cir.f[2:ncol(cir.f)][cir.f[2:ncol(cir.f)] <= -threshold] <- -1
GDSC_ylimit <- cir.f[order(apply(cir.f[2:ncol(cir.f)],1,sum)),'Gene_Symbol']
GDSC_ycolor <- rep("black",times=length(GDSC_ylimit))
GDSC_ycolor[GDSC_ylimit %in% core.circadian] <- "red"
pdf("GDSC/51Circadian_Sign.GDSC_Drug_Exp.cor.morethan_10cpds2.pdf",width=16,height = 9)#,width=4000,height=1500,res=400)
ggplot(Circadian_Sign.GDSC_Drug_Exp.cor.mm,aes(x=Compound_name,y=Gene_Symbol))+
  geom_point(aes(size=pvalue,col=zscore))+
  scale_color_gradient2(low = "green2",mid="white", high = "magenta",midpoint = 0,na.value="white",name="zscore")+ #breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),
  scale_size_continuous(limit=c(-log10(0.05),30),range = c(1, 5),breaks=c(-log10(0.05),10,20,30),labels=c("0.05","1e-10","1e-20","<1e-30"))+
  scale_y_discrete(limit=GDSC_ylimit,expand = c(0.03,0.03))+
  scale_x_discrete(limit=names(cpd.f))+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=12,colour = GDSC_ycolor),
        axis.text.x=element_text(size=12,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()
cpd.fm <- data.frame(cpd.f)
cpd.fm$DRUG.NAME <- rownames(cpd.fm)

####target pathway
GDSC_screen_cpd <- read.csv("GDSC/GDSC_screen_compound.csv",header=T)
GDSC_screen_cpd["TARGET.PATHWAY"][GDSC_screen_cpd["DRUG.NAME"]=="WZ3105"] <- "RTK signaling"
GDSC_screen_cpd["TARGET.PATHWAY"][GDSC_screen_cpd["DRUG.NAME"]=="GL-X-138"] <- "TOR signaling"
GDSC_screen_cpd["TARGET.PATHWAY"][GDSC_screen_cpd["DRUG.NAME"]=="TG101348"] <- "RTK signaling"


GDSC_screen_cpd$DRUG.NAME <- gsub(" |-|,|\\/",".",GDSC_screen_cpd$DRUG.NAME)
GDSC_screen_cpd["DRUG.NAME"][GDSC_screen_cpd["DRUG.NAME"] == "5.Fluorouracil"] <- "X5.Fluorouracil"
GDSC_screen_cpd["DRUG.NAME"][GDSC_screen_cpd["DRUG.NAME"] == "(5Z).7.Oxozeaenol"] <- "X.5Z..7.Oxozeaenol" 
GDSC_screen_cpd["DRUG.NAME"][GDSC_screen_cpd["DRUG.NAME"] == "17.AAG"] <- "X17.AAG"
GDSC_screen_cpd["DRUG.NAME"][GDSC_screen_cpd["DRUG.NAME"] == "Bleomycin.(50.uM)"] <- "Bleomycin..50.uM."
GDSC_screen_cpd.m <- merge(GDSC_screen_cpd,cpd.fm,by="DRUG.NAME")
GDSC_screen_cpd.m <- GDSC_screen_cpd.m[rev(order(GDSC_screen_cpd.m$cpd.f)),]
#GDSC_screen_cpd.m <- unique(GDSC_screen_cpd.m[which(GDSC_screen_cpd.m$TARGET.PATHWAY != "other"),c("DRUG.NAME","TARGET.PATHWAY","cpd.f")])
GDSC_DRUG_circadian_pathway_Count <- table(GDSC_screen_cpd.m$TARGET.PATHWAY)
GDSC_DRUG_circadian_pathway_Count <- data.frame(GDSC_DRUG_circadian_pathway_Count)

write.csv(GDSC_DRUG_circadian_pathway_Count,file="GDSC/51_GDSC_DRUG_circadian_pathway_Count.csv",quote = F,row.names = F)


Circadian_Sign.GDSC_Drug_Exp.cor.mmm <- merge(Circadian_Sign.GDSC_Drug_Exp.cor.mm,GDSC_screen_cpd,by.x="Compound_name",by.y="DRUG.NAME")
unique_cir_cpd <- unique(Circadian_Sign.GDSC_Drug_Exp.cor.mmm[,c(1,9)])
GDSC_screen_cpd.f <- GDSC_screen_cpd[which(GDSC_screen_cpd$TARGET.PATHWAY %in% unique_cir_cpd$TARGET.PATHWAY),]
pathway_order <- names(table(unique_cir_cpd$TARGET.PATHWAY))[rev(order((table(unique_cir_cpd$TARGET.PATHWAY)/table(GDSC_screen_cpd.f$TARGET.PATHWAY))))]
pathway_order <- pathway_order[pathway_order %in% names(table(unique_cir_cpd$TARGET.PATHWAY)[table(unique_cir_cpd$TARGET.PATHWAY) >0])]
Circadian_Sign.GDSC_Drug_Exp.cor.mmm$order <- rep(0,times=nrow(Circadian_Sign.GDSC_Drug_Exp.cor.mmm))
for( i in 1:length(pathway_order)){
  Circadian_Sign.GDSC_Drug_Exp.cor.mmm["order"][Circadian_Sign.GDSC_Drug_Exp.cor.mmm["TARGET.PATHWAY"]==pathway_order[i]] <- i
}

write.table(Circadian_Sign.GDSC_Drug_Exp.cor.mmm,file="Circadian_Sign.GDSC_Drug_Exp.cor.mmm.txt",quote = F,row.names = F,col.names = T,sep="\t")
write.table(cir.f,file="cir.f.txt",quote = F,row.names = F,col.names = T,sep="\t")
GDSC_x <- cir.f[order(apply(cir.f[2:ncol(cir.f)],1,sum)),'Gene_Symbol']
GDSC_y <- unique(Circadian_Sign.GDSC_Drug_Exp.cor.mmm[rev(order(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$order)),1])
GDSC_yyyy <- unique(GDSC_cpd.target[which(GDSC_cpd.target$Compound_name %in% GDSC_y),c(1,3)])
GDSC_yyyy <- GDSC_yyyy[match(GDSC_y,GDSC_yyyy$Compound_name),]
GDSC_yy <- unique(Circadian_Sign.GDSC_Drug_Exp.cor.mmm[rev(order(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$order)),c(1,10)])

GDSC_xcolor <- rep("black",times=length(GDSC_x))
GDSC_xcolor[GDSC_x %in% core.circadian] <- "red"
library(RColorBrewer)
#GDSC_ycol <-  rep(brewer.pal(9, "Set1")[c(1:5,7:9)],times=3)[GDSC_yy$order]
col_select <-  colors()[c(552,555,645,451,640,491,26,562,121,551,31,258,578,11,428,150,142,175,200)]
GDSC_ycol <- col_select[GDSC_yy$order]
Circadian_Sign.GDSC_Drug_Exp.cor.mmm["zscore"][Circadian_Sign.GDSC_Drug_Exp.cor.mmm["zscore"] > 11] <- 11
Circadian_Sign.GDSC_Drug_Exp.cor.mmm["zscore"][Circadian_Sign.GDSC_Drug_Exp.cor.mmm["zscore"] < -11] <- -11
write.csv(Circadian_Sign.GDSC_Drug_Exp.cor.mmm,file="/home/yye1/Circadian/Drug/GDSC/Circadian_Sign.GDSC_Drug_Exp.cor.mmm.csv",quote = F,row.names = F)
pdf("~/Circadian/Drug/GDSC/51Circadian_Sign.GDSC_Drug_Exp.cor.morethan_10cpds_pathway_ycolors.pdf",width=9,height = 12)#,width=4000,height=1500,res=400)
#tiff("GDSC/51Circadian_Sign.GDSC_Drug_Exp.cor.morethan_10cpds_pathway_ycolors2.tiff",width=3000,height = 4300)#,width=4000,height=1500,res=400)
ggplot(Circadian_Sign.GDSC_Drug_Exp.cor.mmm[order(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$order),],aes(y=Compound_name,x=Gene_Symbol))+
  geom_point(aes(size=pvalue,col=zscore))+
# scale_color_gradient2(low = "green2",mid="white", high = "magenta",midpoint = 0,na.value="white",name="z-score")+ #breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),
  scale_color_gradientn(colours=colorRampPalette(c("green2","white","magenta"),space="rgb")(100))+#,
                      # values=c(seq(min(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$zscore),-3,length.out = 10),seq(-3,3,length.out = 81),seq(3,max(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$zscore),length.out = 10)))+
  
   scale_size_continuous(limit=c(-log10(0.05),30),range = c(1, 4.5),breaks=c(-log10(0.05),10,20,30),labels=c("0.05","1e-10","1e-20","<1e-30"),name="p-value")+
  scale_x_discrete(limit=GDSC_x,expand = c(0.02,0.02))+
  scale_y_discrete(limit=GDSC_yyyy$Compound_name,labels=GDSC_yyyy$DRUG.NAME,expand = c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = GDSC_ycol),
        axis.text.x=element_text(size=12,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"),
        #legend.key.width = unit(1,"cm"),
      #  legend.key.heigh = unit(0.2,"cm"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()



per <-  table(unique_cir_cpd$TARGET.PATHWAY)/table(GDSC_screen_cpd.f$TARGET.PATHWAY)
per <- data.frame(per)
per <- per[which(per$Freq > 0),]
colnames(per) <- c("TARGET.PATHWAY","Per")
per$Per <- signif(per$Per*100,digit=3)


GDSC_yy1 <- data.frame(table(GDSC_yy$order))
GDSC_yy2 <- merge(GDSC_yy1,unique(Circadian_Sign.GDSC_Drug_Exp.cor.mmm[,c(9,10)]),by.x="Var1",by.y="order")
GDSC_yy2 <- merge(GDSC_yy2,per,by="TARGET.PATHWAY")
GDSC_yy2$Var1 <- factor(GDSC_yy2$Var1,levels=19:1)
GDSC_yy2$TARGET.PATHWAY_per <- paste(GDSC_yy2$TARGET.PATHWAY," (",GDSC_yy2$Per,"%)",sep="")
GDSC_yy2 <- GDSC_yy2[order(GDSC_yy2$Var1),]
GDSC_yy2$Var1 <- factor(GDSC_yy2$Var1,levels=19:1)

write.table(GDSC_yy2,file="GDSC_yy2.txt",quote = F,row.names = F,col.names = T,sep="\t")
plus.vector<-c()
plus<-0
plus1.vector<-c()
for(i in GDSC_yy2$Freq){
  plus1<-plus+as.numeric(i*0.5)
  plus<-plus+i
  plus1.vector<-c(plus1.vector,plus1)
  plus.vector<-c(plus.vector,plus)
}

pdf("GDSC/51_GDSC_bar_per1.pdf",width=3.8,height=9)
ggplot(GDSC_yy2,aes(x=1,y=rev(Freq),fill=factor(Var1)))+
  geom_bar(stat="identity",position="stack",width=0.001)+
  scale_fill_manual(limits=seq(1,19),values =rev(col_select),guide=F)+
  scale_y_continuous(breaks=plus1.vector,labels=GDSC_yy2$TARGET.PATHWAY_per)+
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(colour=NA),
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=14,colour = rev(col_select)),
        axis.ticks=element_blank(),
        legend.title=element_blank(),legend.text=element_text(size=16))
dev.off()


GDSC_screen_cpd.m <- GDSC_screen_cpd.m[1:30,]
pdf("GDSC/top30_circadianGenes_correlared_to_GDSC_Drug_pathway.pdf",width=10,height=14)
par(mai=c(1,2.5, 0.2, 0.3))
b<-barplot(rev(GDSC_screen_cpd.m$cpd.f),col=adjustcolor("blue", alpha.f = 0.4),main="",horiz=T,space=0.3,names.arg=rev(GDSC_screen_cpd.m$DRUG.NAME),xlim=c(0,14),las=1,axes=FALSE,cex.names =1.3)
axis(2,at=b,labels=rev(GDSC_screen_cpd.m$TARGET.PATHWAY),hadj=0,las=2,line=-1.3,tick=F,cex.axis=1.6)
axis(1,labels=seq(0,14,length.out = 8),at=seq(0,14,length.out = 8),cex.axis=2,lwd=3,mgp=c(3,1.5,0))
#axis(2,labels=c(rep("",dim(b)[1])),at=b[,1],lwd.ticks = 0)
abline(v=0,lwd=5)
par(xpd=T)
mtext("The number of circadian genes",,side=1,line=3.5,cex=2)
dev.off()


CTRP.cpd.fm <- data.frame(CTRP.cpd.f)
CTRP.cpd.fm$DRUG.NAME <- rownames(CTRP.cpd.fm)
CTRP_screen_cpd <- read.csv("CTRP/CTRP_Drug.target.genes.pathway.csv",header=T)
CTRP_screen_cpd.m <- merge(CTRP.cpd.fm,CTRP_screen_cpd,by.x="DRUG.NAME",by.y="cpd_name")
CTRP_screen_cpd.m <- CTRP_screen_cpd.m[rev(order(CTRP_screen_cpd.m$CTRP.cpd.f)),]
CTRP_screen_cpd.m <- unique(CTRP_screen_cpd.m[which(CTRP_screen_cpd.m$target_or_activity_of_compound != "other"),c("DRUG.NAME","target_or_activity_of_compound","CTRP.cpd.f")])

CTRP_DRUG_circadian_pathway_Count <- table(CTRP_screen_cpd.m$target_or_activity_of_compound)[table(CTRP_screen_cpd.m$target_or_activity_of_compound) >=1]
CTRP_DRUG_circadian_pathway_Count <- data.frame(CTRP_DRUG_circadian_pathway_Count)
write.csv(CTRP_DRUG_circadian_pathway_Count,file="CTRP/CTRP_DRUG_circadian_pathway_Count.csv",quote = F)

CTRP_screen_cpd.m <- CTRP_screen_cpd.m[1:30,]
pdf("CTRP/top30_circadianGenes_correlared_to_CTRP_Drug_pathway.pdf",width=16,height=14)
par(mai=c(1,2.8, 0.2, 0.3))
b<-barplot(rev(CTRP_screen_cpd.m$CTRP.cpd.f),col=adjustcolor("blue", alpha.f = 0.4),main="",horiz=T,space=0.3,names.arg=rev(CTRP_screen_cpd.m$DRUG.NAME),xlim=c(0,14),las=1,axes=FALSE,cex.names =1.3)
axis(2,at=b,labels=rev(CTRP_screen_cpd.m$target_or_activity_of_compound),hadj=0,las=2,line=-1.3,tick=F,cex.axis=1.6)
axis(1,labels=seq(0,14,length.out = 8),at=seq(0,14,length.out = 8),cex.axis=2,lwd=3,mgp=c(3,1.5,0))
#axis(2,labels=c(rep("",dim(b)[1])),at=b[,1],lwd.ticks = 0)
abline(v=0,lwd=5)
par(xpd=T)
mtext("The number of circadian genes",,side=1,line=3.5,cex=2)
dev.off()

###The drugs of CTRP and GDSC database correlated to Actional genes
actionable.genes <-read.delim("~/Circadian/expression/correlation/pearson/CRY1_CRY2_cor0.3/actionable.genes.txt",header=T)
CTRP_cpd.target.actionable <- CTRP_cpd.target[which(CTRP_cpd.target$target.genes %in% actionable.genes$Gene),]
CTRP_cpd.target.actionable <- CTRP_cpd.target.actionable[,c(2,1)]
colnames(CTRP_cpd.target.actionable) <- c("DRUG.NAME","TARGET")
GDSC_cpd.target.actionable <- GDSC_cpd.target[which(GDSC_cpd.target$TARGET %in% actionable.genes$Gene),]
cpd.target.actionable <- rbind(CTRP_cpd.target.actionable,GDSC_cpd.target.actionable[,1:2])
cpd.target.actionable <- data.frame(table(cpd.target.actionable$TARGET)[table(cpd.target.actionable$TARGET)>0])
colnames(cpd.target.actionable) <- "DrugNum"
cpd.target.actionable$ActionableGene <- rownames(cpd.target.actionable)
cpd.target.actionable <- merge(cpd.target.actionable,actionable.genes,by.x="ActionableGene",by.y="Gene",all.y=T)
cpd.target.actionable[is.na(cpd.target.actionable)==T] <- 0
write.table(cpd.target.actionable,file="~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/the.number.of.cpd.target.actionable.txt",quote = F,row.names = F,sep='\t')
##indirectly target actional genes

indirect.CTRP_cpd.target.actionable <- Sign.CTRP_Drug_Exp.cor[which(Sign.CTRP_Drug_Exp.cor$Gene_Symbol %in% actionable.genes$Gene),]
indirect.CTRP_cpd.target.actionable <- indirect.CTRP_cpd.target.actionable[,c(3,1)]
colnames(indirect.CTRP_cpd.target.actionable) <- c("DRUG.NAME","TARGET")
indirect.cpd.target.actionable <- rbind(indirect.CTRP_cpd.target.actionable,indirect.GDSC_cpd.target.actionable)
indirect.cpd.target.actionable <- data.frame(table(indirect.cpd.target.actionable$TARGET)[table(indirect.cpd.target.actionable$TARGET)>0])
colnames(indirect.cpd.target.actionable) <- "DrugNum"
indirect.cpd.target.actionable$ActionableGene <- rownames(indirect.cpd.target.actionable)
indirect.cpd.target.actionable <- merge(indirect.cpd.target.actionable,actionable.genes,by.x="ActionableGene",by.y="Gene",all.y=T)
indirect.cpd.target.actionable[is.na(indirect.cpd.target.actionable)] <- 0
write.table(indirect.cpd.target.actionable,file="~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/51_the.number.of.indirect.cpd.target.actionable.txt",quote = F,row.names = F,sep='\t')


####compound directly targeted genes
corfiles <- list.files(path="~/Circadian/expression/correlation/pearson/",pattern = "*_gene_cor.pearson.txt")
pair_tumor <- read.delim("/home/yye1/Circadian/expression/tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)
corfiles <- corfiles[gsub("_gene_cor.pearson.txt","",corfiles) %in% pair_tumor$Type]
corfiles_tumor <- gsub("_gene_cor.pearson.txt","",corfiles)


target_sig <- unique(GDSC_cpd.target.cor[which(GDSC_cpd.target.cor$Compound_name %in% intersect(GDSC_cpd.target.cor$Compound_name,Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound_name)),1:2])
target_sig <- rbind(target_sig,data.frame(Gene_Symbol="NAMPT",Compound_name="FK866"))
Circadian_compound_sign <- Circadian_Sign.GDSC_Drug_Exp.cor.mmm[which(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound_name %in% target_sig$Compound_name),]
for(j in 1:length(corfiles)){
  cor_data <- read.delim(paste("~/Circadian/expression/correlation/pearson/",corfiles[j],sep = ""),header=T)
  for( i in unique(Circadian_compound_sign$Compound_name)){
    cir_select <- as.character(unique(Circadian_compound_sign[which(Circadian_compound_sign$Compound_name == i),"Gene_Symbol"]))
    target_select <- target_sig[which(target_sig$Compound_name==i),"Gene_Symbol"]
    cor_select <- cor_data[which(cor_data$gene %in% target_select),c("gene",cir_select)]
    cor_select.m <- melt(cor_select,id.vars="gene", measure.vars=cir_select)
    colnames(cor_select.m) <- c("Compound_TargetGene","Cir_Gene","Cor")
    cor_select.m$Compound <- rep(i,times=nrow(cor_select.m))
    if (i == "Afatinib"){
      cor_select.mm <- cor_select.m
    }else{
      cor_select.mm <- rbind(cor_select.mm,cor_select.m)
    }
  }
  cor_select.mm$tumor <-  rep(corfiles_tumor[j],times=nrow(cor_select.mm))
  if(j == 1){
    cor_selectAll <- cor_select.mm
  }else{
    cor_selectAll <- rbind(cor_selectAll,cor_select.mm)
  }
}
Kidey_select <- cor_selectAll#[which(cor_selectAll$tumor %in% c("KICH","KIRC","KIRP")),]
Kidey_select$class <- paste(Kidey_select$Compound_TargetGene,Kidey_select$Cir_Gene,Kidey_select$Compound,sep=",")
class <- sapply(split(Kidey_select[,"Cor"], factor(Kidey_select$class)),function(x){ifelse(abs(max(x))>abs(min(x)),max(x),min(x))})
class <- data.frame(cor=class,class=names(class))
class$Compound_TargetGene <- data.frame(do.call(rbind, strsplit(as.character(class$class),',')))$X1
class$Cir_Gene <- data.frame(do.call(rbind, strsplit(as.character(class$class),',')))$X2
class$Compound <- data.frame(do.call(rbind, strsplit(as.character(class$class),',')))$X3
class.M <- unique(merge(class,Circadian_compound_sign[,c(1,9)],by.x="Compound",by.y="Compound_name"))
class.M <- merge(class.M,GDSC_yy2[,1:2],by="TARGET.PATHWAY")
class_y <- unique(class.M[order(class.M$Var1),c(2,5,7)])
EGFR <- data.frame(Compound="Afatinib/Cetuximab/Gefitinib",Compound_TargetGene="EGFR",Var1=7)
class_y <- rbind(class_y[!(class_y$Compound %in% c("Afatinib","Cetuximab","Gefitinib") & class_y$Compound_TargetGene == "EGFR"),],EGFR)
class_y <- class_y[order(class_y$Var1),]
class_y_label <- paste(class_y$Compound,"-->",class_y$Compound_TargetGene)
class_y_col <- rep(brewer.pal(9, "Set1")[c(1:5,7:9)],times=3)[as.numeric(unfactor(class_y$Var1))]
class_y2 <- unique(class.M[order(class.M$Var1),c(1,5,7)])
class_y2 <- merge(class_y,class_y2[,1:2],by="Compound_TargetGene")
class_y2 <- class_y2[order(class_y2$Var1),]
class_y_label2 <- paste(class_y2$Compound,"-->",class_y2$TARGET.PATHWAY,": ",class_y2$Compound_TargetGene,sep="")

FK866  <- cor_selectAll[which(cor_selectAll$Compound== "FK866" & cor_selectAll$tumor %in% c("KICH","KIRC","KIRP")),]
pdf("GDSC/51Circadian_Sign.GDSC_Drug_correlation with compound target genes2.pdf",width=11,height = 3)#,width=4000,height=1500,res=400)
ggplot(class.M[order(class.M$TARGET.PATHWAY),],aes(x=Cir_Gene,y=Compound_TargetGene))+
  geom_tile(aes(fill=cor))+
  scale_fill_gradient2(low="darkgreen",mid="white", high = "magenta",midpoint=0,name="Cor Coeff.")+
  scale_x_discrete(limits=GDSC_x,expand=c(0.002,0.002))+
  scale_y_discrete(limits=class_y$Compound_TargetGene,labels=class_y_label2,expand=c(0.005,0.005))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor = element_blank(),#line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,color=class_y_col),
        axis.text.x=element_text(size=10,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.key.width = unit(0.2,"cm"),
        legend.key.heigh = unit(0.8,"cm"),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()

###scatter 
KIRC <- read.delim("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/KIRC_mRNA_each_exp_20160513",header=T)
KIRC$gene <- data.frame(do.call(rbind, strsplit(as.character(KIRC$gene),'\\|')))$X1
for( i in c("SIRT1","PER3","GNB2L1","PPP5C")){
  KIRC_sub <- KIRC[which(KIRC$gene %in% c("ATM",i)),]
  filter <- colnames(KIRC_sub)[as.numeric(substr(colnames(KIRC_sub),14,15)) ==1]
  filter <- filter[!is.na(filter)]
  KIRC_sub <- KIRC_sub[,filter]
  KIRC_sub.m <- data.frame(t(KIRC_sub))
  colnames(KIRC_sub.m) <- c("ATM",i)
  KIRC_sub.m <- KIRC_sub.m[which(KIRC_sub.m$ATM !=41.6069),]
  PT_cor <- signif(cor(KIRC_sub.m$ATM,KIRC_sub.m[i]),digits=2)
  #P_val <-  signif(cor.test(KIRC_sub.m$ATM,KIRC_sub.m[,i])$p.value,digits = 4)
  pdf(paste("GDSC/KIRC_ATM_",i,".expression.pdf",sep=""),width=5,height=5)
 p <- ggplot(KIRC_sub.m,aes(x=log2(as.numeric(as.vector(KIRC_sub.m[,i]))),y=log2(as.numeric(as.vector(ATM)))))+
    geom_point(color="black")+
    stat_smooth(method = "lm",color="red",size=1)+
    labs(x=paste("ATM (R = ",PT_cor,")",sep=""),
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



cor_selectAll[which(cor_selectAll$Compound== "CP466722" & cor_selectAll$Cir_Gene == "PPP5C"),]


###Top 10 gene correlation "darkgreen","white","magenta"
rgb.palette_inver <- colorRampPalette(c("green2","white","magenta"), space="rgb")
###GDSC
GDSC.Cir <- GDSC.exp[which(GDSC.exp$GENE_SYMBOL %in% circadian.genes),]
rownames(GDSC.Cir) <- GDSC.Cir$GENE_SYMBOL
GDSC.Cir.cor <- cor(t(GDSC.Cir[3:ncol(GDSC.Cir)]))
inverse_genes <-data.frame(GeneSymbol = GDSC_x[c(1:10,(length(GDSC_x)-9):length(GDSC_x))],group=rep(c("green2","magenta"),times=c(10,10)))
GDSC.Cir.cor.top10 <- GDSC.Cir.cor[which(rownames(GDSC.Cir.cor) %in% inverse_genes$GeneSymbol),which(colnames(GDSC.Cir.cor) %in% inverse_genes$GeneSymbol)]
GDSC.Cir.cor.top10[GDSC.Cir.cor.top10==1] <- 0
inverse_genes <- merge(data.frame(GeneSymbol=rownames(GDSC.Cir.cor.top10),order=1:nrow(GDSC.Cir.cor.top10)),inverse_genes,by="GeneSymbol")
inverse_genes <- inverse_genes[order(inverse_genes$order),]
pdf(paste("~/Circadian/Drug/GDSC/inverse_genes_cor/GDSC_CCLE_top10gene correlation.pdf",sep=""),width = 6,height = 6)
heatmap.2(GDSC.Cir.cor.top10,col = rgb.palette_inver(50),zlim=c(-max(abs(GDSC.Cir.cor.top10),max(GDSC.Cir.cor.top10))),dendrogram="both",ColSideColors =as.character(inverse_genes$group),RowSideColors =as.character(inverse_genes$group),margin=c(7,7),
          trace="none",key=T,key.title = "none",key.ylab = "",key.xlab ="Cor Coeff R",density.info="none", keysize=1.15,cexCol = 1.5,cexRow = 1.5,  key.xtickfun=function() {
            breaks <- parent.frame()$breaks
            return(list(
              at=parent.frame()$scale01(c(signif(breaks[1],digits = 1),0,
                                          signif(breaks[length(breaks)],digits = 1))),
              labels=c(as.character(signif(breaks[1],digits = 1)),0,
                       as.character(signif(breaks[length(breaks)],digits = 1)))
            ))
          }) #ColSideColors = genes.list$col,RowSideColors = RowSideColors,breaks=seq(-0.57,0.57,length= 101),zlim=c(-0.57,0.57),
dev.off()
GDSC.Cir.cor[GDSC.Cir.cor==1] <- 0
pdf(paste("~/Circadian/Drug/GDSC/inverse_genes_cor/GDSC_CCLE_top10gene correlationAll.pdf",sep=""),width = 12,height = 12)
heatmap.2(GDSC.Cir.cor,col = rgb.palette_inver(50),zlim=c(-max(abs(GDSC.Cir.cor),max(GDSC.Cir.cor))),dendrogram="both",margin=c(7,7),
          trace="none",key=T,key.title = "none",key.ylab = "",key.xlab ="Cor Coeff R",density.info="none", keysize=0.5,cexCol = 1.5,cexRow = 1.5)
dev.off()

#Compound target correlation with circadian genes

target_sig <- unique(GDSC_cpd.target.cor[which(GDSC_cpd.target.cor$Compound_name %in% intersect(GDSC_cpd.target.cor$Compound_name,Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound_name)),1:2])
#target_sig <- rbind(target_sig,data.frame(Gene_Symbol="NAMPT",Compound_name="FK866"))

Circadian_compound_sign <- Circadian_Sign.GDSC_Drug_Exp.cor.mmm[which(Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound_name %in% target_sig$Compound_name),]

GDSC.Cir.target <- GDSC.exp[which(GDSC.exp$GENE_SYMBOL %in% c(as.character(circadian.genes),as.character(unique(target_sig$Gene_Symbol)))),]
rownames(GDSC.Cir.target) <- GDSC.Cir.target$GENE_SYMBOLS
GDSC.Cir.target.cor <- cor(t(GDSC.Cir.target[3:ncol(GDSC.Cir.target)]))
GDSC.Cir.target.cor <- data.frame(GDSC.Cir.target.cor)
GDSC.Cir.target.cor$gene <- rownames(GDSC.Cir.target.cor)
for( i in unique(Circadian_compound_sign$Compound_name)){
    cir_select <- as.character(unique(Circadian_compound_sign[which(Circadian_compound_sign$Compound_name == i),"Gene_Symbol"]))
    target_select <- target_sig[which(target_sig$Compound_name==i),"Gene_Symbol"]
    
    cor_select <- GDSC.Cir.target.cor[which(rownames(GDSC.Cir.target.cor) %in% target_select),c("gene",cir_select)]
    cor_select.m <- melt(cor_select,id.vars="gene", measure.vars=cir_select)
    colnames(cor_select.m) <- c("Compound_TargetGene","Cir_Gene","Cor")
    cor_select.m$Compound <- rep(i,times=nrow(cor_select.m))
    if (i == "Afatinib"){
      cor_select.mm <- cor_select.m
    }else{
      cor_select.mm <- rbind(cor_select.mm,cor_select.m)
    }
}

cor_select.mm.M <- unique(merge(cor_select.mm,Circadian_compound_sign[,c(1,9)],by.x="Compound",by.y="Compound_name"))
cor_select.mm.M <- merge(cor_select.mm.M,GDSC_yy2[,1:2],by="TARGET.PATHWAY")
cor_select.mm.M <- merge(cor_select.mm.M,GDSC_cpd.target[,c(1,3)],by.x="Compound",by.y="Compound_name")
cor_select.mm.M <- unique(cor_select.mm.M)
cor_select.mm.M <- cor_select.mm.M[!(cor_select.mm.M$Compound_TargetGene %in% c("NAMPT","AKT1","CDK7","PLK3","CAMK1","DDR1","MET","EPHB4")),]
cor_select.mm_y <- unique(cor_select.mm.M[order(cor_select.mm.M$Var1),c(1,2,3,6,7)])
a <- unique(cor_select.mm.M[,c(1:3,7)])
EGFR <- data.frame(Compound="Afatinib/EKB-569/Gefitinib",TARGET.PATHWAY="EGFR signaling",Compound_TargetGene="EGFR",Var1=9,DRUG.NAME="Afatinib/EKB-569/Gefitinib")
TOP1 <- data.frame(Compound="Camptothecin/SN-38",TARGET.PATHWAY="DNA replication",Compound_TargetGene="TOP1",Var1=7,DRUG.NAME="Camptothecin/SN-38")
CDK9 <- data.frame(Compound="AT-7519/KIN001-270/THZ-2-49",TARGET.PATHWAY="cell cycle",Compound_TargetGene="CDK9",Var1=5,DRUG.NAME="AT-7519/KIN001-270/THZ-2-49")
HDAC6 <- data.frame(Compound="CAY10603/Tubastatin A",TARGET.PATHWAY="chromain  histone acetylation",Compound_TargetGene="HDAC6",Var1=2,DRUG.NAME="CAY10603/Tubastatin A")
FLT3 <- data.frame(Compound="CEP-701/WZ3105",TARGET.PATHWAY="RTK signaling",Compound_TargetGene="FLT3",Var1=2,DRUG.NAME="CEP.701/WZ3105")
MTOR <- data.frame(Compound="AZD8055/OSI-027/QL-X-138",TARGET.PATHWAY="RTK signaling",Compound_TargetGene="MTOR",Var1=2,DRUG.NAME="AZD8055/OSI.027/QL.X.138")
JAK2 <- data.frame(Compound="CEP-701/TG101348",TARGET.PATHWAY="RTK signaling",Compound_TargetGene="JAK2",Var1=2,DRUG.NAME="CEP.701/TG101348")
KIT <- data.frame(Compound="Masitinib/OSI-930",TARGET.PATHWAY="RTK signaling",Compound_TargetGene="KIT",Var1=2,DRUG.NAME="Masitinib/OSI.930")
PDK1 <- data.frame(Compound="BX-912/KIN001-244",TARGET.PATHWAY="PI3K signaling",Compound_TargetGene="PDK1",Var1=2,DRUG.NAME="BX.912/KIN001.244")
BCL2 <- data.frame(Compound="Navitoclax/TW 37",TARGET.PATHWAY="Apoptosis",Compound_TargetGene="PDK1",Var1=2,DRUG.NAME="BNavitoclax/TW.37")


cor_select.mm_y <- rbind(cor_select.mm_y[!(cor_select.mm_y$Compound_TargetGene %in% c("EGFR","TOP1","CDK9","HDAC6","BCL2","MTOR","JAK2","KIT","PDK1","FLT3")),],EGFR,TOP1,CDK9,HDAC6)
cor_select.mm_y <- cor_select.mm_y[order(cor_select.mm_y$Var1),]
cor_select.mm_y_label <- paste(cor_select.mm_y$Compound,"-->",cor_select.mm_y$Compound_TargetGene)
cor_select.mm_y_label2 <- paste(cor_select.mm_y$DRUG.NAME,"-->",cor_select.mm_y$TARGET.PATHWAY,": ",cor_select.mm_y$Compound_TargetGene,sep="")
cor_select.mm.M["Cor"][cor_select.mm.M["Cor"] > 0.4] <- 0.4
cor_select.mm.M["Cor"][cor_select.mm.M["Cor"] < -0.4] <- -0.4
write.table(cor_select.mm.M,file="~/Circadian/Drug/GDSC/cor_select.mm.M.txt",quote = F,sep="\t",row.names = F)
pdf("~/Circadian/Drug/GDSC/51_circadian target genes correlation1.pdf",width=10,height = 5)#,width=4000,height=1500,res=400)
ggplot(cor_select.mm.M,aes(x=Cir_Gene,y=Compound_TargetGene))+
  geom_tile(aes(fill=Cor))+
  scale_fill_gradientn(colours=colorRampPalette(c("green2","white","magenta"),space="rgb")(100),name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
 # scale_fill_gradient2(limits=c(-0.4,0.4),low="green2",mid="white", high = "magenta",midpoint=0,name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
  scale_x_discrete(limits=GDSC_x,expand=c(0.002,0.002))+
  scale_y_discrete(limits=cor_select.mm_y$Compound_TargetGene,labels=cor_select.mm_y_label2,expand=c(0.005,0.005))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor = element_blank(),#line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,color=col_select[as.numeric(unfactor(cor_select.mm_y$Var1))]),
        axis.text.x=element_text(size=10,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.width = unit(1,"cm"),
        legend.key.heigh = unit(0.3,"cm"),
        legend.key = element_rect(fill="white",colour = "black"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()

pdf("~/Circadian/Drug/GDSC/51_circadian target genes correlation_ytargetgenes.pdf",width=7,height =5)#,width=4000,height=1500,res=400)
ggplot(cor_select.mm.M,aes(x=Cir_Gene,y=Compound_TargetGene))+
  geom_tile(aes(fill=Cor))+
  scale_fill_gradientn(colours=colorRampPalette(c("green2","white","magenta"),space="rgb")(100),name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
  
 # scale_fill_gradient2(limits=c(-0.4,0.4),low="green2",mid="white", high = "magenta",midpoint=0,name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
  scale_x_discrete(limits=GDSC_x,expand=c(0.002,0.002))+
  scale_y_discrete(limits=cor_select.mm_y$Compound_TargetGene,labels=cor_select.mm_y$Compound_TargetGene,expand=c(0.005,0.005))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor = element_blank(),#line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,color=col_select[as.numeric(unfactor(cor_select.mm_y$Var1))]),
        axis.text.x=element_text(size=10,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.width = unit(1,"cm"),
        legend.key.heigh = unit(0.3,"cm"),
        legend.key = element_rect(fill="white",colour = "black"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()
####Link Circadian genes and drug target genes
cor_select.mm.MM <- cor_select.mm.M[which(abs(cor_select.mm.M$Cor)>=0.25),]
yy_order <- data.frame(Compound_TargetGene=cor_select.mm_y$Compound_TargetGene,locate.targets = seq(0,100,length.out = 28),color=col_select[as.numeric(unfactor(cor_select.mm_y$Var1))])
cir <- sapply(split(cor_select.mm.MM[,"Cor"],cor_select.mm.MM$Cir_Gene),sum)
cir <- cir[cir != 0]
cir_order <- data.frame(Cir_Gene=names(cir[order(cir)]),locate.genes = seq(2,96,length.out = 25))
cir_acOverlap <- merge(merge(cor_select.mm.MM,yy_order,by="Compound_TargetGene"),cir_order,by="Cir_Gene")
cir_acOverlap$colorline <- rep("magenta",nrow(cir_acOverlap))
cir_acOverlap["colorline"][cir_acOverlap["Cor"] < 0] <- "green2" 



pdf("~/Circadian/Drug/GDSC/target genes circadian genes lines.pdf",width=7,height =5)#,width=4000,height=1500,res=400)
ggplot(yy_order)+
  geom_text(aes(x=3,y=locate.targets,label=Compound_TargetGene),color=yy_order$color,hjust = 1,size=4)+
  geom_text(data=cir_order,aes(x=4.5,y=locate.genes,label=Cir_Gene),hjust = 0,size=4)+
  scale_x_continuous(limit=c(1,6))+
  geom_segment(aes(x = 3.05, y = locate.targets, xend = 4.45, yend = locate.genes),color = cir_acOverlap$colorline, data = cir_acOverlap)+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
dev.off()

pdf("~/Circadian/Drug/GDSC/51_circadian target genes correlation_yunderline.pdf",width=7,height = 5)#,width=4000,height=1500,res=400)
ggplot(cor_select.mm.M,aes(x=Cir_Gene,y=Compound_TargetGene))+
  geom_tile(aes(fill=Cor))+
  scale_fill_gradientn(colours=colorRampPalette(c("green2","white","magenta"),space="rgb")(100),name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
  #scale_fill_gradient2(limits=c(-0.4,0.4),low="green2",mid="white", high = "magenta",midpoint=0,name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
  scale_x_discrete(limits=GDSC_x,expand=c(0.002,0.002))+
  scale_y_discrete(limits=cor_select.mm_y$Compound_TargetGene,labels=rep("--------|",times=length(cor_select.mm_y$Compound_TargetGene)),expand=c(0.005,0.005))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor = element_blank(),#line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,color=col_select[as.numeric(unfactor(cor_select.mm_y$Var1))]),
        axis.text.x=element_text(size=10,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.width = unit(1,"cm"),
        legend.key.heigh = unit(0.3,"cm"),
        legend.key = element_rect(fill="white",colour = "black"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()




pdf("~/Circadian/Drug/GDSC/51_circadian target genes correlation_ycompound.pdf",width=7,height = 5)#,width=4000,height=1500,res=400)
ggplot(cor_select.mm.M,aes(x=Cir_Gene,y=Compound_TargetGene))+
  geom_tile(aes(fill=Cor))+
  scale_fill_gradient2(limits=c(-0.4,0.4),low="green2",mid="white", high = "magenta",midpoint=0,name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
  scale_x_discrete(limits=GDSC_x,expand=c(0.002,0.002))+
  scale_y_discrete(limits=cor_select.mm_y$Compound_TargetGene,labels=cor_select.mm_y$DRUG.NAME,expand=c(0.005,0.005))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor = element_blank(),#line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,color=col_select[as.numeric(unfactor(cor_select.mm_y$Var1))]),
        axis.text.x=element_text(size=10,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.width = unit(1,"cm"),
        legend.key.heigh = unit(0.3,"cm"),
        legend.key = element_rect(fill="white",colour = "black"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()

pdf("~/Circadian/Drug/GDSC/51_circadian target genes correlation_ysignaling.pdf",width=7,height = 5)#,width=4000,height=1500,res=400)
ggplot(cor_select.mm.M,aes(x=Cir_Gene,y=Compound_TargetGene))+
  geom_tile(aes(fill=Cor))+
  scale_fill_gradient2(limits=c(-0.4,0.4),low="green2",mid="white", high = "magenta",midpoint=0,name="Cor Coeff. R",breaks=c(-0.25,0.25),labels=c("Neg","Pos"))+
  scale_x_discrete(limits=GDSC_x,expand=c(0.002,0.002))+
  scale_y_discrete(limits=cor_select.mm_y$Compound_TargetGene,labels=cor_select.mm_y$TARGET.PATHWAY,expand=c(0.005,0.005))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        panel.grid.minor = element_blank(),#line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,color=col_select[as.numeric(unfactor(cor_select.mm_y$Var1))]),
        axis.text.x=element_text(size=10,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.width = unit(1,"cm"),
        legend.key.heigh = unit(0.3,"cm"),
        legend.key = element_rect(fill="white",colour = "black"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()

###
target_cir <- unique(cor_select.mm[,1:2])
target_cir$P_val <- rep(NA,times=nrow(target_cir))
for(i in 1:nrow(target_cir)){
  target_cir_each <- target_cir[1,]
  p_val <- cor.test(as.numeric(GDSC.exp[which(as.character(GDSC.exp$GENE_SYMBOLS) == target_cir_each[,1]),2:ncol(GDSC.exp)]),as.numeric(GDSC.exp[which(as.character(GDSC.exp$GENE_SYMBOLS) == target_cir_each[,2]),2:ncol(GDSC.exp)]))$p.value
  target_cir[i,"P_val"] <- p_val
}










for( Compound in unique(GDSC.Drug$DRUG.NAME)){
  EachDrug <-GDSC.Drug[which(GDSC.Drug$DRUG.NAME == Compound & (GDSC.Drug$COSMIC_ID %in% colnames(GDSC.exp))),]
  if(length(EachDrug$COSMIC_ID) > length(unique(EachDrug$COSMIC_ID))){
    EachDrug<- t(sapply(split(EachDrug[,c("LN_IC50","AUC")],EachDrug$COSMIC_ID),colMeans))
    EachDrug<-data.frame(EachDrug)
    EachDrug$COSMIC_ID <- rownames(EachDrug)
  }
  gene.cor <- apply( as.matrix(GDSC.exp[, EachDrug$COSMIC_ID]) , 1 , cor , y = EachDrug$AUC) 
  zscore <- log((1+gene.cor)/(1-gene.cor))/2*sqrt(nrow(EachDrug)-3)
  gene.cor <- data.frame(gene.cor)
  colnames(gene.cor) <- Compound
  zscore <- data.frame(zscore)
  colnames(zscore) <- paste(Compound,"_zscore",sep="")
  p_value <- apply(as.matrix(GDSC.exp[, EachDrug$COSMIC_ID]),1,function(x){cor.test(x, EachDrug$AUC)$p.value})
  #FDR_value <- p.adjust(p_value,method="fdr")
  p_value <- data.frame(p_value)
  colnames(p_value) <- paste(Compound,"_p",sep = "")
  GDSC_Drug_Exp.cor <- cbind(GDSC_Drug_Exp.cor,gene.cor,zscore,p_value)
}
Circadian_Sign.GDSC_Drug_Exp.cor.mmm["Compound_name"][Circadian_Sign.GDSC_Drug_Exp.cor.mmm["Compound_name"] == "X5.Fluorouracil"] <- "5.Fluorouracil"
Circadian_Sign.GDSC_Drug_Exp.cor.mmm["Compound_name"][Circadian_Sign.GDSC_Drug_Exp.cor.mmm["Compound_name"] == "X.5Z..7.Oxozeaenol"] <- "5Z..7.Oxozeaenol"
Circadian_Sign.GDSC_Drug_Exp.cor.mmm["Compound_name"][Circadian_Sign.GDSC_Drug_Exp.cor.mmm["Compound_name"] == "X17.AAG"] <- "17.AAG"


### JQ1
EachDrug <-GDSC.Drug[which(GDSC.Drug$DRUG.NAME == "JQ1" & (GDSC.Drug$COSMIC_ID %in% colnames(GDSC.exp))),]
if(length(EachDrug$COSMIC_ID) > length(unique(EachDrug$COSMIC_ID))){
  EachDrug<- t(sapply(split(EachDrug[,c("LN_IC50","AUC")],EachDrug$COSMIC_ID),colMeans))
  EachDrug<-data.frame(EachDrug)
  EachDrug$COSMIC_ID <- rownames(EachDrug)
}
JQ1_targetGeneExp <- GDSC.exp[which(GDSC.exp$GENE_SYMBOLS %in% c("CRY1","NR1D2","PPARG","PRKAA2")),c("GENE_SYMBOLS", EachDrug$COSMIC_ID)]
JQ1_targetGeneExp.m <- melt(JQ1_targetGeneExp,id.vars ="GENE_SYMBOLS",measure.vars = EachDrug$COSMIC_ID)
JQ1_targetGeneExp.m$AUC <- rep(EachDrug$AUC,each=4)
gene.cor <- apply( as.matrix(JQ1_targetGeneExp[, EachDrug$COSMIC_ID]) , 1 , cor , y = EachDrug$AUC) 
JQ1_cor <- data.frame(gene.cor=signif(gene.cor,digits = 2),zscore=Circadian_Sign.GDSC_Drug_Exp.cor[which(Circadian_Sign.GDSC_Drug_Exp.cor$Compound_name == "JQ1"),"zscore"],GENE_SYMBOLS=unique(JQ1_targetGeneExp.m$GENE_SYMBOLS))
JQ1_targetGeneExp.m$GENE_SYMBOLS <- factor(JQ1_targetGeneExp.m$GENE_SYMBOLS,levels = unique(JQ1_targetGeneExp.m$GENE_SYMBOLS))
levels(JQ1_targetGeneExp.m$GENE_SYMBOLS) <- paste(JQ1_cor$GENE_SYMBOLS, " ( z=",signif(as.numeric(JQ1_cor$zscore),digits=2),")",sep="")

pdf("~/Circadian/Drug/GDSC/the correlation of JQ1 AUC with circadian genes in GDSC.pdf",width = 7,height = 7)
ggplot(data = JQ1_targetGeneExp.m, aes(x = AUC, 
                                     y = value))+
  geom_point(size=0.9)+
  stat_smooth(method = "lm",color="red",size=1)+
  facet_wrap(~GENE_SYMBOLS,ncol=2,scales = "free_y")+
  #geom_text(data=JQ1_cor,aes(x=0.4,y=10,label=paste("z = ",signif(zscore,digits = 3),sep="")),size=5)+
  ylab('Gene Expression')+xlab("JQ1 AUC")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=14,color="black"),axis.title=element_text(size=16,color="black"),
        axis.text.x=element_text(size=14,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=14),
        strip.background = element_rect(fill=NA))
dev.off()

EachDrug <- CTRP.Drug[which(CTRP.Drug$cpd_name == "JQ-1"),]
EachDrug <- EachDrug[which(EachDrug$master_ccl_id %in% CTRP_CCLE_CellLine$master_ccl_id),]
filter <- as.character(EachDrug$ccl_name)
JQ1_CTRP_TargetExp <- CCLE_exp[which(CCLE_exp$Gene_Symbol %in% as.character(Circadian_Sign.CTRP_Drug_Exp.cor[which(tolower(Circadian_Sign.CTRP_Drug_Exp.cor$Compound_name) == "jq.1"),"Gene_Symbol"])),c("Gene_Symbol",sort(filter))]
JQ1_CTRP_TargetExp.m <- melt(JQ1_CTRP_TargetExp,id.vars = "Gene_Symbol",measure.vars = sort(filter))
JQ1_CTRP_TargetExp.m$AUC <- rep( EachDrug[order(EachDrug$ccl_name),"area_under_curve"],each=nrow(JQ1_CTRP_TargetExp))

gene.cor <- apply( as.matrix(2^CCLE_exp[which(CCLE_exp$Gene_Symbol %in% as.character(Circadian_Sign.CTRP_Drug_Exp.cor[which(tolower(Circadian_Sign.CTRP_Drug_Exp.cor$Compound_name) == "jq.1"),"Gene_Symbol"])),sort(filter)]) , 1 , cor , y = EachDrug[order(EachDrug$ccl_name),"area_under_curve"])
JQ1_CTRP_cor <- data.frame(gene.cor=signif(gene.cor,digits = 2),GENE_SYMBOLS=CCLE_exp[which(CCLE_exp$Gene_Symbol %in% as.character(Circadian_Sign.CTRP_Drug_Exp.cor[which(tolower(Circadian_Sign.CTRP_Drug_Exp.cor$Compound_name) == "jq.1"),"Gene_Symbol"])),2])
zscore <- Circadian_Sign.CTRP_Drug_Exp.cor[which(Circadian_Sign.CTRP_Drug_Exp.cor$Compound_name == "JQ.1"),c("Gene_Symbol","zscore")]
JQ1_CTRP_cor <- merge(zscore,JQ1_CTRP_cor,by.x="Gene_Symbol",by.y="GENE_SYMBOLS")
JQ1_CTRP_TargetExp.m$Gene_Symbol <- factor(JQ1_CTRP_TargetExp.m$Gene_Symbol,levels = unique(JQ1_CTRP_TargetExp.m$Gene_Symbol))
levels(JQ1_CTRP_TargetExp.m$Gene_Symbol) <- paste(JQ1_CTRP_cor$Gene_Symbol, " ( z=",signif(as.numeric(JQ1_CTRP_cor$zscore),digits=2),")",sep="")

pdf("~/Circadian/Drug/CTRP/the correlation of JQ1 AUC with circadian genes in CTRP.pdf",width = 12,height = 12)
ggplot(data = JQ1_CTRP_TargetExp.m, aes(x = AUC, 
                                        y = value))+
  geom_point(size=0.9)+
  stat_smooth(method = "lm",color="red",size=1)+
  facet_wrap(~Gene_Symbol,ncol=5,scale="free_y")+
  #geom_text(data=JQ1_CTRP_cor,aes(x=0.4,y=14,label=paste("z = ",signif(zscore,digits = 3),sep="")),size=5)+
  ylab('Gene Expression')+xlab("JQ1 AUC")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=12,color="black"),axis.title=element_text(size=16,color="black"),
        axis.text.x=element_text(size=12,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=14),
        strip.background = element_rect(fill=NA))
dev.off() 


###Correlation larger than 0.2
GDSC_Drug_Exp.cor <- read.delim("/home/yye1/Circadian/Drug/GDSC/GDSC_Drug_Exp_correlation.txt",header=T)
GDSC_Drug_Exp.cor.only <- melt(GDSC_Drug_Exp.cor,id.vars="Gene_Symbol",measure.vars=gsub("_zscore","",grep("_zscore",colnames(GDSC_Drug_Exp.cor)[2:ncol((GDSC_Drug_Exp.cor))],value=T)))
colnames(GDSC_Drug_Exp.cor.only) <- c("Gene_Symbol","Compound_name","cor")
GDSC_Drug_Exp.cor.onlySign <- GDSC_Drug_Exp.cor.only[which(GDSC_Drug_Exp.cor.only$Compound %in% GDSC_yyyy$Compound_name),]
GDSC_Drug_Exp.cor.onlySign <- GDSC_Drug_Exp.cor.onlySign[which(GDSC_Drug_Exp.cor.onlySign$Gene_Symbol %in% circadian.genes),]
GDSC_Drug_Exp.cor.onlySign <- GDSC_Drug_Exp.cor.onlySign[which(abs(GDSC_Drug_Exp.cor.onlySign$cor) >= 0.2),]
GDSC_Drug_Exp.cor.onlySign <-  merge(GDSC_Drug_Exp.cor.onlySign, GDSC_Drug_Exp.cor.p,by=c("Gene_Symbol","Compound_name"))
onlySignPvalue <- -log10(GDSC_Drug_Exp.cor.onlySign$p)
onlySignPvalue[is.infinite(onlySignPvalue)] <- 30
onlySignPvalue[onlySignPvalue>=30] <- 30
GDSC_Drug_Exp.cor.onlySign$pvalue <- onlySignPvalue

pdf("~/Circadian/Drug/GDSC/51Circadian_Sign.GDSC_Drug_Exp_correlation morethan 0.2.pdf",width=9,height = 12)#,width=4000,height=1500,res=400)
#tiff("GDSC/51Circadian_Sign.GDSC_Drug_Exp.cor.morethan_10cpds_pathway_ycolors2.tiff",width=3000,height = 4300)#,width=4000,height=1500,res=400)
ggplot(GDSC_Drug_Exp.cor.onlySign,aes(y=Compound_name,x=Gene_Symbol))+
  geom_point(aes(size=pvalue,col=cor))+
  scale_color_gradient2(low = "green2",mid="white", high = "magenta",midpoint = 0,na.value="white",name="cor")+ #breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),
  #scale_color_gradient2(&quot;red&quot;, mid = &quot;gray&quot;, high=&quot;blue&quot,midpoint = 0,na.value="white",name="z-score")+ #breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),
  scale_size_continuous(limit=c(-log10(0.05),30),range = c(1, 4.5),breaks=c(-log10(0.05),10,20,30),labels=c("0.05","1e-10","1e-20","<1e-30"),name="p-value")+
  scale_x_discrete(limit=GDSC_x,expand = c(0.02,0.02))+
  scale_y_discrete(limit=GDSC_yyyy$Compound_name,labels=GDSC_yyyy$DRUG.NAME,expand = c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = GDSC_ycol),
        axis.text.x=element_text(size=12,colour = GDSC_xcolor,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"),
        #legend.key.width = unit(1,"cm"),
        #  legend.key.heigh = unit(0.2,"cm"),
        legend.position = "bottom",legend.direction = "horizontal")
dev.off()

###CTRP and GDSC overlap
Circadian_Sign.CTRP_Drug_Exp.cor$Compound_withoutann <- tolower(gsub(" |-|,|\\/|\\.","",Circadian_Sign.CTRP_Drug_Exp.cor$Compound))
Circadian_Sign.GDSC_Drug_Exp.cor$Compound_withoutann <- tolower(gsub(" |-|,|\\/|\\.","",Circadian_Sign.GDSC_Drug_Exp.cor$Compound))
Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_withoutann <- tolower(gsub(" |-|,|\\/|\\.","",Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound))
Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound_withoutann <- tolower(gsub(" |-|,|\\/|\\.","",Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound))

CTRP_GDSC_sign_overlap <- intersect(Circadian_Sign.CTRP_Drug_Exp.cor$Compound_withoutann,Circadian_Sign.GDSC_Drug_Exp.cor$Compound_withoutann)
CTRP_GDSC_sign10_overlap <- intersect(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_withoutann,Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound_withoutann)

drug_ct<- read.delim("Anticancer drugs with  circadian timing.txt",header=F)
drug_ct$Compound_withoutann <-  tolower(gsub(" |-|,|\\/|\\.","",drug_ct$V1))
GDSC_knownDrug <- intersect(drug_ct$Compound_withoutann,Circadian_Sign.GDSC_Drug_Exp.cor$Compound_withoutann)
GDSC_knownDrug_cor <- Circadian_Sign.GDSC_Drug_Exp.cor[which(Circadian_Sign.GDSC_Drug_Exp.cor$Compound_withoutann %in% GDSC_knownDrug),]
for( i in unique(GDSC_knownDrug_cor$Compound_withoutann)){
  a <- GDSC_knownDrug_cor[which(GDSC_knownDrug_cor$Compound_withoutann == i),]
  print(c("pos",i,a[(a$zscore>0),"Gene_Symbol"]))
  print(c("neg",i,a[(a$zscore<0),"Gene_Symbol"]))
}

CTRP_knownDrug <- intersect(drug_ct$Compound_withoutann,Circadian_Sign.CTRP_Drug_Exp.cor$Compound_withoutann)
CTRP_knownDrug_cor <- Circadian_Sign.CTRP_Drug_Exp.cor[which(Circadian_Sign.CTRP_Drug_Exp.cor$Compound_withoutann %in% CTRP_knownDrug),]
for( i in unique(CTRP_knownDrug_cor$Compound_withoutann)){
  a <- CTRP_knownDrug_cor[which(CTRP_knownDrug_cor$Compound_withoutann == i),]
  print(c("pos",i,as.character(a[(a$zscore>0),"Gene_Symbol"])))
  print(c("neg",i,as.character(a[(a$zscore<0),"Gene_Symbol"])))
}
(intersect(GDSC_cpd.target$DRUG.NAME,CTRP.Drug.ann$cpd_name))
GDSC_CTRP_Drug.ann <- read.csv("GDSC_CCLE_CTRP_Drug_conversion.csv")[,c("GDSC.name","CTRP.name","CTRP.master.cpd.id")]
GDSC_CTRP_Drug.ann <- GDSC_CTRP_Drug.ann[which(GDSC_CTRP_Drug.ann$CTRP.name != "N/A"),]
GDSC_CTRP_Drug.ann$CTRP.name <- sub(" |-|,|\\/",".",GDSC_CTRP_Drug.ann$CTRP.name)
intersect(GDSC_CTRP_Drug.ann$CTRP.name,Circadian_Sign.CTRP_Drug_Exp.cor$Compound_name)
setdiff(tolower(GDSC_CTRP_Drug.ann$GDSC.name),union(intersect(tolower(GDSC_CTRP_Drug.ann$GDSC.name),tolower(Circadian_Sign.GDSC_Drug_Exp.cor$Compound_name)),intersect(tolower(GDSC_CTRP_Drug.ann$CTRP.name),tolower(Circadian_Sign.GDSC_Drug_Exp.cor$Compound_name))))


a<- union(intersect(tolower(GDSC_CTRP_Drug.ann$GDSC.name),tolower(Circadian_Sign.GDSC_Drug_Exp.cor$Compound_name)),intersect(tolower(GDSC_CTRP_Drug.ann$CTRP.name),tolower(Circadian_Sign.GDSC_Drug_Exp.cor$Compound_name)))
b <- union(intersect(tolower(GDSC_CTRP_Drug.ann$GDSC.name),tolower(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_name)),intersect(tolower(GDSC_CTRP_Drug.ann$CTRP.name),tolower(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_name)))

c <- Circadian_Sign.CTRP_Drug_Exp.cor.mmm[which(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_name %in% intersect(a,b)),]
b <- union(intersect(tolower(GDSC_CTRP_Drug.ann$GDSC.name),tolower(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_name)),intersect(tolower(GDSC_CTRP_Drug.ann$CTRP.name),tolower(Circadian_Sign.CTRP_Drug_Exp.cor.mmm$Compound_name)))
Circadian_Sign.CTRP_Drug_Exp.cor[which(tolower(Circadian_Sign.CTRP_Drug_Exp.cor$Compound_name) == "jq.1"),]

overlapdrug_than10_CTRP <- CTRP_screen_cpd.m[which(CTRP_screen_cpd.m$master_cpd_id %in% GDSC_CTRP_Drug.ann$CTRP.master.cpd.id),c(1,3)]
intersect(overlapdrug_than10_CTRP$DRUG.NAME,Circadian_Sign.GDSC_Drug_Exp.cor.mmm$Compound_name)



######

GDSC.exp <- read.delim("/home/yye1/Circadian/Drug/GDSC/Cell_line_RMA_proc_basalExp.txt",header=T)
corCircadian <- c("PER1","PER2","PER3","CRY1","CRY2","ARNTL","ARNTL2","RORA","RORB","RORC","NPAS2","CLOCK","NR1D1","NR1D2")
GDSC.core <- GDSC.exp[which(GDSC.exp$GENE_SYMBOLS %in% corCircadian),]

GDSC.core.MedianZ <- t(apply(GDSC.core[,intersect(unique(GDSC.Drug$COSMIC_ID),colnames(GDSC.core.MedianZ))],1,function(x){(x-median(x))/(sd(x))}))
GDSC.core.MedianZ <- data.frame(GeneSymbol=GDSC.core$GENE_SYMBOLS,GDSC.core.MedianZ)

matrix_Median <- as.matrix(GDSC.core.MedianZ[,2:ncol(GDSC.core.MedianZ)])
row.names(matrix_Median) <- GDSC.core$GENE_SYMBOLS
matrix_Median[matrix_Median > 3] <- 3
matrix_Median[matrix_Median < -3] <- -3
heatmap(matrix_Median[,order(SumGDSC.core.MedianZ)], Colv = NA,col=greenred(100),labCol = NA,cexRow = 2)

SumGDSC.core.MedianZ <- apply(GDSC.core.MedianZ[,2:ncol(GDSC.core.MedianZ)],2,sum)

GDSC.Drug <- read.csv("/home/yye1/Circadian/Drug/GDSC/v17_fitted_dose_response.csv",header=T)
Drug.cellIC50 <- data.frame(CellLine = names(SumGDSC.core.MedianZ))
for( Compound in unique(GDSC.Drug$DRUG.NAME)){
  EachDrug <-GDSC.Drug[which(GDSC.Drug$DRUG.NAME == Compound & (GDSC.Drug$COSMIC_ID %in% colnames(GDSC.core.MedianZ))),]
  if(length(EachDrug$COSMIC_ID) > length(unique(EachDrug$COSMIC_ID))){
    EachDrug<- t(sapply(split(EachDrug[,c("LN_IC50","AUC")],EachDrug$COSMIC_ID),colMeans))
    EachDrug<-data.frame(EachDrug)
    EachDrug$COSMIC_ID <- rownames(EachDrug)
  }
  Drugline <- data.frame(CellLine = EachDrug$COSMIC_ID,Compound = EachDrug$LN_IC50)
  Drug.cellIC50 <- merge(Drug.cellIC50,Drugline,by="CellLine",all.x=T)
  gene.cor <- cor.test(EachDrug$LN_IC50,SumGDSC.core.MedianZ[match(EachDrug$COSMIC_ID,names(SumGDSC.core.MedianZ))],method = "spearman")[c("estimate","p.value")]
  gene.corA <- data.frame(Compound,Cor=gene.cor$estimate,Pval = gene.cor$p.value)
  if(Compound == "Erlotinib"){
    gene.corAll <- gene.corA
  }else{
    gene.corAll <- rbind(gene.corAll,gene.corA)
  }
} 

colnames(Drug.cellIC50) <- c("CellLine",as.character(unique(GDSC.Drug$DRUG.NAME)))
gene.corAll.sign <- gene.corAll[which(gene.corAll$Pval < 0.05),]
Drug.cellIC50.m <- as.matrix(Drug.cellIC50[,gene.corAll.sign$Compound])
rownames(Drug.cellIC50.m) <- Drug.cellIC50$CellLine
Drug.cellIC50.m <- t(Drug.cellIC50.m[match(rownames(Drug.cellIC50.m),names(sort(SumGDSC.core.MedianZ))),])
Drug.cellIC50.mm <- t(apply(Drug.cellIC50.m,1,scale))
Drug.cellIC50.mm[Drug.cellIC50.mm>3] <- 3
Drug.cellIC50.mm[Drug.cellIC50.mm < -3] <- -3
heatmap(Drug.cellIC50.mm, na.color = "black",Colv = NA,scale=c("column"),col=greenred(100),labCol = NA,cexRow = 2,na.rm=T)





#### find circadian genes CNV
setwd("/home/yye1/Circadian/CNV")
library(ggplot2)
library(reshape2)
files.names <- list.files(path="/extraspace/TCGA/TCGA_CNV/Gene/",pattern=".txt")
pair_tumor <- read.delim("/home/yye1/Circadian/expression/tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)
files.names <- files.names[gsub(".txt|cnv_gene_","",files.names) %in% pair_tumor$Type]
files.Abs <- gsub(".txt","",files.names)
files.Abs <- gsub("cnv_gene_","",files.Abs)
circadian.genes <-  read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
#circadian.genes <- data.frame(do.call(rbind, strsplit(as.character(circadian.genes$V1),'\\|')))$X1
circadian.genes.CNV <- data.frame(circadian.genes$V1)
colnames(circadian.genes.CNV) <- "gene"
sample.num <- c("Tumor.type","Number")
CNV_count <- function(x){
  gain <- signif(length(x[x > log2(4/2)])/length(x) * 100,digits = 3)
  loss <- signif(length(x[x < log2(0.5)])/length(x) * 100,digits = 3)
  CNV <- c(gain,loss)
}
for(m in 1:length(files.names)){
  BLCA <- read.delim(paste("/extraspace/TCGA/TCGA_CNV/Gene/",files.names[m],sep=""),header=T)
  BLCA <- BLCA[which(BLCA$Gene.Symbol %in% circadian.genes$V1),] ##get circadian genes
  BLCA.names <- colnames(BLCA)[2:ncol(BLCA)]
  count <- c(files.Abs[m],length(BLCA.names))
  sample.num <- rbind(sample.num,count)
  CNV <- apply(BLCA[2:ncol(BLCA)],1,CNV_count)
  CNV <- data.frame(t(CNV))
  colnames(CNV) <- c(paste(files.Abs[m],".gain",sep=""),paste(files.Abs[m],".loss",sep=""))
  CNV$gene <- BLCA$Gene.Symbol
  circadian.genes.CNV <- merge(circadian.genes.CNV,CNV,by="gene")
}
circadian.genes.CNV[is.na(circadian.genes.CNV)] <- 0
write.table(circadian.genes.CNV,file="circadian.genes.CNV_percentage.txt",quote = F,row.names = F,sep="\t")
write.csv(sample.num,file="tumor.number.use.to.CNV_percentage.csv",quote = F)
circadian.genes.CNV <- read.delim("circadian.genes.CNV_percentage.txt",header=T)
cir.CNV.gain <- melt(circadian.genes.CNV,id.vars="gene",measure.vars=colnames(circadian.genes.CNV)[grep(".gain",colnames(circadian.genes.CNV))])
cir.CNV.lose <- melt(circadian.genes.CNV,id.vars="gene",measure.vars=colnames(circadian.genes.CNV)[grep(".loss",colnames(circadian.genes.CNV))])
cir.CNV.gain["value"][cir.CNV.gain["value"] < 5] <- 0
cir.CNV.lose["value"][cir.CNV.lose["value"] < 5 ] <- 0
cir.CNV <- rbind(cir.CNV.gain,cir.CNV.lose)
cir.CNV$class <- factor(c(rep("Gain",times=nrow(cir.CNV.gain)),rep("Loss",times=nrow(cir.CNV.lose))))
colnames(cir.CNV) <- c("geneSymbol","Type","Percentage","Class")
cir.CNV$Type <- gsub(".loss","",cir.CNV$Type)
cir.CNV$Type <- gsub(".gain","",cir.CNV$Type)

#gain_index <- grep(".gain",colnames(circadian.genes.CNV))
#lose_index <- grep(".loss",colnames(circadian.genes.CNV))
#cir.order <- circadian.genes.CNV
#cir.order[gain_index][(cir.order[gain_index ] < 5)] <- 0
#cir.order[gain_index][(cir.order[gain_index ] >= 5)] <- 1
#cir.order[lose_index][(cir.order[lose_index ] < 5)] <- 0
#cir.order[lose_index][(cir.order[lose_index ] >= 5)] <- -1
x_label_color <- read.table( "~/Circadian/expression/x_label_color.txt",sep="\t",header=T)
y_label <- read.table("~/Circadian/expression/expression_y_label.txt",header=T)
cir.CNV$geneSymbol <- factor(cir.CNV$geneSymbol,levels = x_label_color$x_label)
cir.CNV$Type <- factor(cir.CNV$Type,levels=y_label$x)
rgb.palette <- colorRampPalette(c("blue","white","red"), space="rgb")
cir.CNV.s <- cir.CNV
pdf("circadian.genes.CNV_percentage1.pdf",width=13,height = 5)#,width=4000,height=1500,res=400)
par(mai=c(1,1,0.5,1.5))
ggplot(cir.CNV,aes(x=geneSymbol,y=Type))+
  geom_point(aes(color=Class,size=Percentage))+
  scale_color_manual(limits = c("Gain","Loss"),values=c("Gain"=rgb.palette(100)[90],"Loss"=rgb.palette(100)[10]),name="CNV % ")+
  scale_y_discrete(limit=y_label.m$x,expand = c(0.05,0.05))+
  scale_x_discrete(limit=x_label_color$x_label,expand = c(0.02,0.5))+
  scale_size_continuous(limit=c(1,max(cir.CNV$Percentage)),range = c(1, 8),breaks=c(5,10,20,30),labels=c(5,10,20,">=30"),name = "Pvalue")+
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


###
cir.CNV.gain <- melt(circadian.genes.CNV,id.vars="gene",measure.vars=colnames(circadian.genes.CNV)[grep(".gain",colnames(circadian.genes.CNV))])
cir.CNV.lose <- melt(circadian.genes.CNV,id.vars="gene",measure.vars=colnames(circadian.genes.CNV)[grep(".loss",colnames(circadian.genes.CNV))])
cir.CNV.lose$value <- -cir.CNV.lose$value
cir.CNV.s <- rbind(cir.CNV.gain,cir.CNV.lose)
colnames(cir.CNV.s) <- c("geneSymbol","Type","Percentage")
cir.CNV.s$Type <- gsub(".loss","",cir.CNV.s$Type)
cir.CNV.s$Type <- gsub(".gain","",cir.CNV.s$Type)
cir.CNV.s["Percentage"][cir.CNV.s["Percentage"] >= 40] <- 40
cir.CNV.s["Percentage"][cir.CNV.s["Percentage"] < -40] <- -40
pdf("circadian.genes.CNV.FC1.5.percentage_tile.pdf",width=12,height = 5)#,width=4000,height=1500,res=400)
ggplot(cir.CNV.s,aes(x=geneSymbol,y=Type,fill=Percentage))+
  geom_tile(color="black")+
  # scale_fill_gradient2(low = "steelblue3",mid="white", high = "firebrick3",midpoint = 0,na.value="white",breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),name="log2 FC")+
  scale_fill_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,breaks=seq(-30,30,length.out = 5),labels=c(" -60","-30","0","30","60"),name="Percentage")+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit=y_label$x,expand=c(0,0))+
  scale_x_discrete(limit=x_label_color$x_label,expand = c(0,0))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = x_label_color$x_color,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()

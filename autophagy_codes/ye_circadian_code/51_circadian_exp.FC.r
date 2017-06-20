setwd("/home/yye1/Circadian/expression")
###mRNA expression data from TCGA mRNA expression
files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
pair_tumor <- read.delim("tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)
files.names <- files.names[which(gsub("_mRNA_each_exp_20160513","",files.names) %in% pair_tumor$Type)]
files.Abs <- gsub("_mRNA_each_exp_20160513","",files.names)
circadian.genes <- read.table("../circadian.genes.txt",header=F) 
circadian.genes.exp <- circadian.genes
colnames(circadian.genes.exp) <- "gene"
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
filter.mm <- c("Type","Tumor","Normal")
for(m in 1:length(files.names)){
  BLCA <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",files.names[m],sep=""),header=T)
  BLCA$gene <- data.frame(do.call(rbind,strsplit(as.character(BLCA$gene),"\\|")))$X1
  BLCA <- BLCA[which(BLCA$gene %in% circadian.genes$V1),] ##get circadian genes
  ###consider tumor type contain Tumor and normal genes
  BLCA.names <- colnames(BLCA)[2:ncol(BLCA)] 
  BLCA.names <- data.frame(BLCA.names)
  BLCA.names$BLCA.names <- as.character(BLCA.names$BLCA.names)
  BLCA.names <- data.frame(do.call(rbind, strsplit(BLCA.names$BLCA.names,'\\.')))
  BLCA.names <- BLCA.names[which(BLCA.names$X1=="TCGA"),]
  filter <- c()
  #BLCA.pair.FC <- data.frame(gene)
  for(i in unique(BLCA.names$X3)){
    participant <- BLCA.names[which(BLCA.names$X3 == i),]
    sample <- as.numeric(substr(BLCA.names[which(BLCA.names$X3 == i),"X4"],1,2))
    if ((nrow(participant) >=2)){
      sub <- apply( BLCA.names[grep(i,BLCA.names$X3),] , 1 , paste , collapse = "." )
      filter <- c(filter,sub)
    }
  }
  
  if(length(filter) >= 10){
    BLCA.f <- BLCA[,filter]
    t <- colnames(BLCA.f)[as.numeric(substr(colnames(BLCA.f),14,15)) == 1]
    n <- colnames(BLCA.f)[as.numeric(substr(colnames(BLCA.f),14,15)) ==11]
    if((length(t)>=10) & (length(n)>=10)){
      BLCA.genes.exp <- data.frame(BLCA$gene)
      colnames(BLCA.genes.exp) <- "gene"
      p_value <- apply(as.matrix(BLCA.f),1,function(x){t.test(x[t],x[n])$p.value})
      p_value <- data.frame(p_value)
      colnames(p_value) <- paste(files.Abs[m],"_pval",sep = "")
      FC <- log2((apply(BLCA.f[,t],1,mean)+0.1)/(apply(BLCA.f[,n],1,mean)+0.1))
      FC <- data.frame(FC)
      colnames(FC) <- files.Abs[m]
      #FC.p <- rbind(FC,p_value)
      BLCA.genes.exp <- cbind(BLCA.genes.exp,FC,p_value)
      circadian.genes.exp <- merge(circadian.genes.exp,BLCA.genes.exp,by="gene")
      filter.m <- c(files.Abs[m],length(t),length(n))
      filter.mm <- rbind(filter.mm,filter.m)
    }
  }
}
###
#write.table(filter.mm,file="tumor.normal.sample.calculate_morethan_10pairs.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(circadian.genes.exp,file="51_circadian.genes.exp.FC_morethan_10pairs.txt",quote=F,row.names=F,sep="\t")
filter.mm <- read.delim("tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)
library(ggplot2)
library(reshape2)
circadian.genes.exp <- read.delim("~/Circadian/expression/51_circadian.genes.exp.FC_morethan_10pairs.txt",header=T)
circadian.genes.exp[30:43] <- apply(circadian.genes.exp[,grep("_pval",colnames(circadian.genes.exp),value=T)],2,function(x){p.adjust(x,method="fdr")})
circadian.genes.exp[2:ncol(circadian.genes.exp)] <- signif(circadian.genes.exp[2:ncol(circadian.genes.exp)],digits = 2)

colnames(circadian.genes.exp) <- c(colnames(circadian.genes.exp)[1:29],gsub("_pval","_FDR",grep("_pval",colnames(circadian.genes.exp),value=T)))
pval_names <- grep("_pval",colnames(circadian.genes.exp),value=T)
FC_names <- colnames(circadian.genes.exp)[!(colnames(circadian.genes.exp) %in% c("gene","geneSymbol",pval_names,FDR_names))]
FDR_names <- gsub("_pval","_FDR",grep("_pval",colnames(circadian.genes.exp),value=T))
#circadian.genes.exp$geneSymbol <- data.frame(do.call(rbind, strsplit(as.character(circadian.genes.exp$gene),'\\|')))$X1
cir.exp.m <- melt(circadian.genes.exp,id.vars="gene",measure.vars=FC_names)
cir.pval.m <- melt(circadian.genes.exp,id.vars="gene",measure.vars=pval_names)
cir.fdr.m <- melt(circadian.genes.exp,id.vars = "gene",measure.vars = FDR_names)
cir.exp.pval <- cbind(cir.exp.m,cir.pval.m,cir.fdr.m)[c(1:3,6,9)]
colnames(cir.exp.pval) <- c("geneSymbol","Type","Fold","Pvalue","FDR")
cir.exp.pval$Fold_label <- cir.exp.pval$Fold

cir.exp.pval[,"Fold_label"][as.numeric(cir.exp.pval[,"Fold_label"]) >= log2(1.5)] <- 2^(cir.exp.pval[,"Fold_label"][cir.exp.pval[,"Fold_label"] >=log2(1.5)])
cir.exp.pval["Fold_label"][cir.exp.pval["Fold_label"] <= log2(2/3)] <- -1/2^(cir.exp.pval["Fold_label"][cir.exp.pval["Fold_label"] <= log2(2/3)])
cir.exp.pval <- cir.exp.pval[which(abs(cir.exp.pval$Fold) >= log2(1.5) & (cir.exp.pval$FDR <= 0.05)),]
#write.table(cir.exp.pval,file="~/Circadian/expression/cir.exp.pval_FC1.5.txt",sep="\t",quote = F,row.names = F)
#cir.exp.pval$Fold_label <- as.character(signif(cir.exp.pval$Fold_label,digits = 2))
cir.exp.pval[3][(cir.exp.pval[3]< log2(1.5))&(cir.exp.pval[3])> -log2(1.5) ] <- 0
cir.exp.pval[3][cir.exp.pval[3] > 2 ] <- 2
cir.exp.pval[3][cir.exp.pval[3] < -2 ] <- -2

cir.exp.pval["Fold_label"][abs(cir.exp.pval["Fold_label"]) < 1.5 ] <- 0
cir.exp.pval["Fold_label"][cir.exp.pval["Fold_label"] > 4 ] <- 4
cir.exp.pval["Fold_label"][cir.exp.pval["Fold_label"] < -4 ] <- -4


cir.exp.pval$Pvalue <- -log10(cir.exp.pval$Pvalue)
cir.exp.pval[4][cir.exp.pval[4] > 15] <- 15
#cir.exp.pval$class <- c(rep("exp",times=nrow(cir.exp.m)),rep("pval",times=nrow(cir.pval.m)))
cir.order <- circadian.genes.exp[,c("gene",FC_names)]
cir.order[2:ncol(cir.order)][(cir.order[2:ncol(cir.order)]<= -log2(1.5))&(circadian.genes.exp[pval_names] <= 0.05)] <- -1
cir.order[2:ncol(cir.order)][(cir.order[2:ncol(cir.order)]>= log2(1.5))&(circadian.genes.exp[pval_names] <= 0.05)] <- 1
cir.order[2:ncol(cir.order)][(((cir.order[2:ncol(cir.order)]< log2(1.5))&(cir.order[2:ncol(cir.order)]> -log2(1.5)))|(circadian.genes.exp[pval_names] > 0.05))] <- 0
# ###
colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
y_label <- names(table(cir.exp.pval$Type)[order(table(cir.exp.pval$Type))])
write.table(y_label,file="~/Circadian/expression/expression_y_label.txt",quote = F,row.names = F,sep="\t")
x_label <- as.character(cir.order[order(apply(cir.order[,FC_names], 1,sum)),"gene"])
x_color <- rep("black",times=length(x_label))
x_color[x_label %in% core.circadian] <- "red"
x_label_color <- data.frame(x_label,x_color)
write.table(x_label_color,file = "~/Circadian/expression/x_label_color.txt",sep="\t",quote = F,row.names = F)
x_label_color <- read.delim("~/Circadian/expression/x_label_color.txt")

pdf("51_circadian.genes.tumor.normal.morethan_10pairs_final_F1.5.pdf",width=12,height = 5)#,width=4000,height=1500,res=400)
ggplot(cir.exp.pval,aes(x=geneSymbol,y=Type,fill=Fold))+
  geom_tile(color="black")+
  # scale_fill_gradient2(low = "steelblue3",mid="white", high = "firebrick3",midpoint = 0,na.value="white",breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),name="log2 FC")+
  scale_fill_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,na.value="white",breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),name="log2 FC")+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit= y_label)+
  scale_x_discrete(limit= x_label,expand=c(0,0))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = x_color,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()

###point
pdf("~/Circadian/expression/51_circadian.genes.tumor.normal.morethan_10pairs_final_F1.5_point_foldlabel.pdf",width=11,height = 4)#,width=4000,height=1500,res=400)
ggplot(cir.exp.pval,aes(x=geneSymbol,y=Type))+
  geom_point(aes(size=Pvalue,col=Fold_label))+
  scale_color_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,na.value="white",breaks=seq(-4,4,length.out = 5),labels=c("<= -4","-1","0","1",">= 4"),name="Fold change")+
  scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 6),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit= y_label,expand = c(0.05,0.05))+
  scale_x_discrete(limit= x_label,expand = c(0.02,0.5))+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = x_color,angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()

cir.exp.pval["alteration"] <- rep("Upregulated",times=nrow(cir.exp.pval))
cir.exp.pval["alteration"][cir.exp.pval["Fold"] < 0] <- "Downregulated"
pdf("51_circadian.genes.tumor.normal.morethan_10pairs_count_FC1.5.pdf",width=10,height = 2)#,width=4000,height=1500,res=400)
ggplot(cir.exp.pval,aes(x=geneSymbol,fill=factor(alteration)))+
  geom_bar(color=NA,width = 0.5)+
  scale_fill_manual(limit=c("Downregulated","Upregulated"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(limit= c(-0.1,12.5),expand = c(0,0),breaks = seq(0,12,length.out = 5))+
  scale_x_discrete(limit= x_label,expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_line(linetype="dashed",color="lightgray"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.x = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()
#######




files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
pair_tumor <- read.delim("tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)
files.names <- files.names[which(gsub("_mRNA_each_exp_20160513","",files.names) %in% pair_tumor$Type)]
files.Abs <- gsub("_mRNA_each_exp_20160513","",files.names)
circadian.genes <- read.table("../circadian.genes.txt",header=F) 

for(m in 1:length(files.names)){
  BLCA <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",files.names[m],sep=""),header=T)
  
  BLCA.names <- colnames(BLCA)[2:ncol(BLCA)] 
  BLCA.names <- data.frame(BLCA.names)
  BLCA.names$BLCA.names <- as.character(BLCA.names$BLCA.names)
  BLCA.names <- data.frame(do.call(rbind, strsplit(BLCA.names$BLCA.names,'\\.')))
  BLCA.names <- BLCA.names[which(BLCA.names$X1=="TCGA"),]
  filter <- c()
  #BLCA.pair.FC <- data.frame(gene)
  for(i in unique(BLCA.names$X3)){
    participant <- BLCA.names[which(BLCA.names$X3 == i),]
    sample <- as.numeric(substr(BLCA.names[which(BLCA.names$X3 == i),"X4"],1,2))
    if ((nrow(participant) >=2)){
      sub <- apply( BLCA.names[grep(i,BLCA.names$X3),] , 1 , paste , collapse = "." )
      filter <- c(filter,sub)
    }
  }
  
  t <-filter[as.numeric(substr(filter,14,15)) == 1]
  n <- filter[as.numeric(substr(filter,14,15)) == 11]
  pair_sample <- data.frame(barcode=c(t,n),class=rep(c("tumor","normal"),times=c(length(t),length(n))),type=rep(files.Abs[m],times=length(c(t,n))))
  if(m==1){
    pair_sampleAll <- pair_sample
  }else{
    pair_sampleAll <- rbind(pair_sampleAll,pair_sample)
  }
}
pair_sampleAll$barcode <- gsub("\\.","\\-",pair_sampleAll$barcode)
write.table(pair_sampleAll,file="pair_sampleAll.txt",quote = F,row.names = F,sep="\t")
paste0(data.frame(do.call(rbind,strsplit(as.character(cir.exp.pval$Pvalue),"\\|")))$X1)
cir.exp.pval$labels_Rsq <- paste(gsub("e","ghf10^",cir.exp.pval$Pvalue))
  #paste(data.frame(do.call(rbind,strsplit(as.character(cir.exp.pval$Pvalue),"e")))$X1,"*","10^", format(data.frame(do.call(rbind,strsplit(as.character(cir.exp.pval$Pvalue),"e")))$X2, digits=2),sep="")

ggplot(cir.exp.pval,aes(x=geneSymbol,y=Type,fill=Fold))+
  geom_tile(color="black")+
  # scale_fill_gradient2(low = "steelblue3",mid="white", high = "firebrick3",midpoint = 0,na.value="white",breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),name="log2 FC")+
  scale_fill_gradient2(low = "blue",mid="white", high = "red",midpoint = 0,na.value="white",breaks=seq(-3,3,length.out = 5),labels=c("<= -3","-1.5","0","1.5",">= 3"),name="log2 FC")+
  geom_text(aes(label=as.character(labels_Rsq)), parse = TRUE, size=4)



Core.cir.count <- t(sapply(split(cir.exp.pval[,"Fold"],list(cir.exp.pval$geneSymbol)),function(x){c(length(x[x < 0]),length(x[x > 0]),(14-length(x[x < 0])-length(x[x > 0])))}))
Core.cir.count <- data.frame(Core.cir.count)
colnames(Core.cir.count) <- c("Down","Up","None")
Core.cir.count$Gene <- rownames(Core.cir.count)
#Core.cir.count <-Core.cir.count[which(Core.cir.count$Gene %in% core.circadian),]

Core.cir.count.m <- melt(Core.cir.count,id.vars = c("Gene"),measure.vars =c("Down","Up","None"))
Core.cir.count.m$per <- Core.cir.count.m$value/14

strip.x.color <- rep("black",times=length(unique(Core.cir.count.m$Gene)))
strip.x.color[sort(unique(Core.cir.count.m$Gene)) %in% core.circadian] <- "red"
pdf("/home/yye1/Circadian/expression/Percentage of cancer types in circadian genes expression FC1.pdf",width = 9,height=9)
ggplot(data = Core.cir.count.m, aes(x = factor(1), y = per, fill = factor(variable))) +
  geom_bar(stat = "identity", position = "stack", color = NA,alpha = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  coord_polar("y") +
  facet_grid(~Gene) +
  # geom_text(aes(x = factor(1), y= .5, label = val_mod, vjust = 4.5)) +
  theme(axis.text=element_blank(),axis.title = element_blank(), 
        panel.background = element_blank(),
        legend.title=element_blank(), axis.ticks = element_blank(),
        #strip.text.y = element_text(angle =0,hjust=0,color="black",size=11),strip.background = element_blank(),
        strip.text.x = element_text(color=strip.x.color,size=11,angle = 90,vjust = 0),legend.text = element_text(size=14),
        legend.position = "bottom",panel.spacing  = unit(0.02, "lines")) +
  scale_fill_manual(limits=c("Down","Up","None"),values=c("blue","red","lightgray"))
dev.off()

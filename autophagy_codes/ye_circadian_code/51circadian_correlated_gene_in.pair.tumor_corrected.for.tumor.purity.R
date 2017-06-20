#### find circadian genes correlated genes R > 0.5
setwd("/extraspace/yye1/analysis/Circadian/expression/correlation")
###mRNA expression data from TCGA mRNA expression
files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
pair_tumor <- read.delim("/home/yye1/Circadian/expression/tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)

circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
spearman_cor <- function(x,y){cor(x, y, use = "everything",method = "spearman")}
circadian.genes.pos.t <- circadian.genes
colnames(circadian.genes.pos.t) <- "gene"
circadian.genes.neg.t <- circadian.genes
colnames(circadian.genes.neg.t) <- "gene"
circadian.genes.pos.n <- circadian.genes.pos.t
circadian.genes.neg.n <- circadian.genes.neg.t
sample.num <- c("Tumor.type","Number")
purity.files <-list.files(path="/extraspace/yye1/analysis/Immunotherapy/bundance/tumor_purity_summary/",pattern=".txt")
purity.files <- purity.files[(gsub("AGP-|.txt","",purity.files) %in% tolower(pair_tumor$Type))]
files.names <- files.names[gsub("_mRNA_each_exp_20160513","",files.names) %in% toupper(gsub("AGP-|.txt","",purity.files))]
files.Abs <- gsub("_mRNA_each_exp_20160513","",files.names)

for(m in 1:length(files.names)){
  BLCA <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",files.names[m],sep=""),header=T)
  # BLCA <- BLCA[which(BLCA$gene %in% circadian.genes$V1),] ##get circadian genes
  ###consider tumor type contain Tumor and normal genes
  BLCA$gene <- data.frame(do.call(rbind,strsplit(as.character(BLCA$gene),"\\|")))$X1
  BLCA.names <- colnames(BLCA)[2:ncol(BLCA)] 
  BLCA.names <- data.frame(BLCA.names)
  BLCA.names$BLCA.names <- as.character(BLCA.names$BLCA.names)
  out <- strsplit(BLCA.names$BLCA.names,'\\.') 
  BLCA.names <- do.call(rbind, out)
  BLCA.names <- data.frame(BLCA.names)
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
  purity.data <- read.delim(paste("/extraspace/yye1/analysis/Immunotherapy/bundance/tumor_purity_summary/AGP-",tolower(files.Abs[m]),".txt",sep=""),header=T)
  a <- intersect(substr(purity.data$sampleid,1,12),substr(filter,1,12))
  filter <- filter[substr(filter,1,12) %in% a]
  if(length(filter) >= 10){
    BLCA.f <- BLCA[,filter]
    t <- colnames(BLCA.f)[as.numeric(substr(colnames(BLCA.f),14,15)) == 1]
    t <- t[!duplicated(substr(t,1,12))]
    t <- t[order(substr(t,1,12))]
    n <- colnames(BLCA.f)[as.numeric(substr(colnames(BLCA.f),14,15)) == 11]
    if((length(t)>=10) & (length(n)>=10)){
      BLCA.filter.t <- BLCA[,c("gene",t)]
      #  BLCA.filter.t <- BLCA.filter.t[(apply(BLCA.filter.t[,2:ncol(BLCA.filter.t)],1,sum) > 1),]
      BLCA.filter.n <- BLCA[,c("gene",n)]
      #  BLCA.filter.n <- BLCA.filter.n[(apply(BLCA.filter.n[,2:ncol(BLCA.filter.n)],1,sum) > 1),]
      BLCA.cor.t <- data.frame(BLCA.filter.t$gene)
      BLCA.cor.n <- data.frame(BLCA.filter.n$gene)
      colnames(BLCA.cor.t) <- "gene"
      colnames(BLCA.cor.n) <- "gene"
      pos_cor.t.all <- c()
      neg_cor.t.all <- c()
      pos_cor.n.all <- c()
      neg_cor.n.all <- c()
      for(j in circadian.genes$V1){
        cir.exp.t <- log2(as.matrix(BLCA.filter.t[which(BLCA.filter.t$gene==as.character(j)),2:ncol(BLCA.filter.t)])+1)
        cir.exp.n <- log2(as.matrix(BLCA.filter.n[which(BLCA.filter.n$gene==as.character(j)),2:ncol(BLCA.filter.n)])+1)
        purity.per <- purity.data[which(purity.data$sampleid %in% substr(t,1,16)),]
        purity.per <- purity.per[order(purity.per$sampleid),"purity"]
        p_value.t <- apply(log2(as.matrix(BLCA.filter.t[,2:ncol(BLCA.filter.t)]+1)),1,function(x){pcor.test(x,as.vector(as.numeric(cir.exp.t)),purity.per,method="pearson")$p.value})
        FDR_value.t <- p.adjust(p_value.t,method="fdr")
        FDR_value.t <- data.frame(FDR_value.t)
        colnames(FDR_value.t) <- paste(files.Abs[m],"_FDR",sep = "")
        p_value.n <- apply(log2(as.matrix(BLCA.filter.n[,2:ncol(BLCA.filter.n)])+1),1,function(x){cor.test(x,as.vector(as.numeric(cir.exp.n)),method="pearson")$p.value})
        FDR_value.n <- p.adjust(p_value.n,method="fdr")
        FDR_value.n <- data.frame(FDR_value.n)
        colnames(FDR_value.n) <- paste(files.Abs[m],"_FDR",sep = "")
        gene.cor.t <- apply(log2(as.matrix(BLCA.filter.t[,2:ncol(BLCA.filter.t)])+1) , 1 ,function(x){pcor.test(x,as.vector(as.numeric(cir.exp.t)),purity.per,method="pearson")$estimate})
        gene.cor.t <- data.frame(gene.cor.t)
        colnames(gene.cor.t) <- j
        BLCA.cor.t <- cbind(BLCA.cor.t,gene.cor.t,FDR_value.t)
        gene.cor.n <- apply(log2(as.matrix(BLCA.filter.n[,2:ncol(BLCA.filter.n)])+1), 1, function(x){cor.test(x,as.vector(as.numeric(cir.exp.n)),method="pearson")$estimate})
        gene.cor.n <- data.frame(gene.cor.n)
        colnames(gene.cor.n) <- j
        BLCA.cor.n <- cbind(BLCA.cor.n,gene.cor.n,FDR_value.n)
        # calculate correlated genes number (R > 0.3 or < -0.3,FDR < 0.05)
      }
      BLCA.cor.t[is.na(BLCA.cor.t)] <- 0
      BLCA.cor.n[is.na(BLCA.cor.n)] <- 0
      write.table(BLCA.cor.t,file=paste("/extraspace/yye1/analysis/Circadian/expression/correlation/",files.Abs[m],"_gene_cor.t_pearson.txt",sep = ""),quote = F,row.names = F,sep="\t")
      write.table(BLCA.cor.n,file=paste("/extraspace/yye1/analysis/Circadian/expression/correlation/",files.Abs[m],"_gene_cor.n_pearson.txt",sep = ""),quote = F,row.names = F,sep="\t")
  }
  filter.m <- c(files.Abs[m],length(t),length(n))
  print(filter.m)
}}


###Correlation larger than 0.5
setwd("/extraspace/yye1/analysis/Circadian/expression/correlation")
##circadian.genes.tumor
corfiles.t <- list.files(pattern = "*t_pearson.txt")
corfiles_tumor <- gsub("_gene_cor.t_pearson.txt","",corfiles.t)
cir_names <- gsub("\\|",".",circadian.genes$V1)
neg_cor.all.t <- data.frame(cir_names)
pos_cor.all.t <- data.frame(cir_names)
for(i in 1:length(corfiles.t)){
  cor_data <- read.delim(corfiles.t[i],sep = "",header=T)
  pos <- c()
  neg <- c()
  cir_names.FDR <-grep("_FDR",colnames(cor_data),value=T)
  for(j in 1:length(cir_names)){
    pos.num <- nrow(cor_data[which(cor_data[cir_names[j]] > 0.5 & cor_data[cir_names.FDR[j]] < 0.05),])
    neg.num <- nrow(cor_data[which(cor_data[cir_names[j]] < -0.5 & cor_data[cir_names.FDR[j]] < 0.05),])
    pos <- c(pos,pos.num)
    neg <- c(neg,neg.num)
  }
  neg_cor.all.t <- cbind(neg_cor.all.t,neg)
  pos_cor.all.t <- cbind(pos_cor.all.t,pos)
}
colnames(pos_cor.all.t) <- c("Gene",corfiles_tumor)
colnames(neg_cor.all.t) <- c("Gene",corfiles_tumor)

##circadian.genes.normal
corfiles.n <- list.files(pattern = "*n_pearson.txt")
corfiles_normal <- gsub("_gene_cor.n_pearson.txt","",corfiles.n)
cir_names <- gsub("\\|",".",circadian.genes$V1)
neg_cor.all.n <- data.frame(cir_names)
pos_cor.all.n <- data.frame(cir_names)
for(i in 1:length(corfiles.n)){
  cor_data.n <- read.delim(corfiles.n[i],sep = "",header=T)
  pos <- c()
  neg <- c()
  cir_names.FDR <-grep("_FDR",colnames(cor_data),value=T)
  for(j in 1:length(cir_names)){
    pos.num <- nrow(cor_data[which(cor_data[cir_names[j]] > 0.5 & cor_data[cir_names.FDR[j]] < 0.05),])
    neg.num <- nrow(cor_data[which(cor_data[cir_names[j]] < -0.5 & cor_data[cir_names.FDR[j]] < 0.05),])
    pos <- c(pos,pos.num)
    neg <- c(neg,neg.num)
  }
  neg_cor.all.n <- cbind(neg_cor.all.n,neg)
  pos_cor.all.n <- cbind(pos_cor.all.n,pos)
}
colnames(pos_cor.all.n) <- c("Gene",corfiles_normal)
colnames(neg_cor.all.n) <- c("Gene",corfiles_normal)

#####

pos_cor.all.t.m <- melt(pos_cor.all.t,id.vars="Gene",measure.vars=colnames(pos_cor.all.t)[2:(ncol(pos_cor.all.t)-1)])
colnames(pos_cor.all.t.m) <- c("GeneSymbol","Type","Num")
pos.t.index <- 2:(ncol(pos_cor.all.t)-1)
pos_cor.all.n.m <- melt(pos_cor.all.n,id.vars="Gene",measure.vars=colnames(pos_cor.all.n)[2:(ncol(pos_cor.all.n)-1)])
colnames(pos_cor.all.n.m) <- c("GeneSymbol","Type","Num")
pos.n.index <- 2:(ncol(pos_cor.all.n)-1)

pos_Yorder1 <- pos_cor.all.t[order(apply(pos_cor.all.t[,pos.t.index], 1,sum)),"Gene"]
pos_Xorder1 <- names(sort(apply(pos_cor.all.t[,pos.t.index], 2,sum)))
pos_fill_limit1 <- c(0,max(c(pos_cor.all.n.m$Num,pos_cor.all.t.m$Num)))
pos_ycolor <- rep("black",times=length(pos_Yorder1))
pos_ycolor[pos_Yorder1 %in% core.circadian] <- "red"

neg_cor.all.t.m <- melt(neg_cor.all.t,id.vars="Gene",measure.vars=colnames(neg_cor.all.t)[2:(ncol(neg_cor.all.t)-1)])
colnames(neg_cor.all.t.m) <- c("GeneSymbol","Type","Num")
neg.t.index <- 2:(ncol(neg_cor.all.t)-1)
neg_cor.all.n.m <- melt(neg_cor.all.n,id.vars="Gene",measure.vars=colnames(neg_cor.all.n)[2:(ncol(neg_cor.all.n)-1)])
colnames(neg_cor.all.n.m) <- c("GeneSymbol","Type","Num")
neg.n.index <- 2:(ncol(neg_cor.all.n)-1)
neg_Yorder1 <- neg_cor.all.t[order(apply(neg_cor.all.t[,neg.t.index], 1,sum)),"Gene"]
neg_Xorder1 <- names(sort(apply(neg_cor.all.t[,neg.t.index], 2,sum)))
neg_fill_limit1 <- c(0,max(c(neg_cor.all.n.m$Num,neg_cor.all.t.m$Num)))
neg_ycolor <- rep("black",times=length(neg_Yorder1))
neg_ycolor[neg_Yorder1 %in% core.circadian] <- "red"

pdf("pos_cor.all.t.pearson.R0.5.pdf",width=10,height = 12)#,width=4000,height=1500,res=400)
ggplot(pos_cor.all.t.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=as.numeric(Num)),col="white")+
  scale_fill_gradient(low = "white",high = "red",na.value="white",limit=pos_fill_limit1,name="Number")+
  scale_y_discrete(limit=pos_Yorder1)+
  scale_x_discrete(limit=pos_Xorder1)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = pos_ycolor),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5),
    axis.ticks=element_line("black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    #legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank(),
    legend.key = element_rect(colour = "black"))+
  geom_text(aes(label=Num),size=6)
dev.off()
pdf("pos_cor.all.n.pearson.R0.5.pdf",width=10,height = 12)#,width=4000,height=1500,res=400)
ggplot(pos_cor.all.n.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=as.numeric(Num)),col="white")+
  scale_fill_gradient(low = "white",high = "red",na.value="white",limit=pos_fill_limit1,name="Number")+
  scale_y_discrete(limit=pos_Yorder1)+
  scale_x_discrete(limit=pos_Xorder1)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = pos_ycolor),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5),
    axis.ticks=element_line("black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    #legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank(),
    legend.key = element_rect(colour = "black"))+
  geom_text(aes(label=Num),size=6)
dev.off()

pdf("neg_cor.all.t.pearson.R0.5.pdf",width=10,height = 12)#,width=4000,height=1500,res=400)
ggplot(neg_cor.all.t.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=as.numeric(Num)),col="white")+
  scale_fill_gradient(low = "white",high = "blue",na.value="white",limit=neg_fill_limit1,name="Number")+
  scale_y_discrete(limit=neg_Yorder1)+
  scale_x_discrete(limit=neg_Xorder1)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = neg_ycolor),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5),
    axis.ticks=element_line("black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    #legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank(),
    legend.key = element_rect(colour = "black"))+
  geom_text(aes(label=Num),size=6)
dev.off()


pdf("neg_cor.all.t.pearson.R0.5.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(neg_cor.all.t.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient(low = "white",high = "blue",na.value="white",limits=neg_fill_limit1,name="Number")+
  scale_y_discrete(limit=neg_Yorder1)+
  scale_x_discrete(limit=neg_Xorder1)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_rect(colour="black",fill=NA,size=2),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = "black",face = "bold"),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5,face = "bold"),
    axis.ticks=element_line(color="black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank())+
  geom_text(aes(label=Num),size=6)
dev.off()


pdf("neg_cor.all.n.pearson.R0.5.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(neg_cor.all.n.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient(low = "white",high = "blue",na.value="white",limits=neg_fill_limit1,name="Number")+
  scale_y_discrete(limit=neg_Yorder1)+
  scale_x_discrete(limit=neg_Xorder1)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_rect(colour="black",fill=NA,size=2),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = "black",face = "bold"),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5,face = "bold"),
    axis.ticks=element_line(color="black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank())+
  geom_text(aes(label=Num),size=6)
dev.off()




###minus
Neg.minus <- neg_cor.all.t[2:14] -neg_cor.all.n[2:14]
Neg.minus$GeneSymbol <- neg_cor.all.t$Gene
Neg.minus <- Neg.minus[!(Neg.minus$GeneSymbol %in% c("FOXO1","FOXM1","FOXO3","FOXA2")),]
Neg.minus$GeneSymbol <- factor(Neg.minus$GeneSymbol,levels = Neg.minus$GeneSymbol)
Neg.minus.m <- melt(Neg.minus,id.vars = "GeneSymbol",measure.vars = colnames(Neg.minus)[1:13])
colnames(Neg.minus.m) <- c("GeneSymbol","Type","Num")
neg.t.index <- 1:13
Neg.minus_Yorder <- Neg.minus[order(apply(abs(Neg.minus[,neg.t.index]), 1,sum)),"GeneSymbol"]
Neg.minus_Xorder <- names(sort(apply(abs(Neg.minus[,neg.t.index]), 2,sum)))
Neg.minus_fill_limit <- c(min(Neg.minus.m$Num),max(Neg.minus.m$Num))
limits <- max(abs(Neg.minus_fill_limit))
negminus_ycolor <- rep("black",times=length(Neg.minus_Yorder))
negminus_ycolor[Neg.minus_Yorder %in% core.circadian] <- "red"
pdf("neg_cor.all.tumor.minus.normal.pearson.R0.5_bluered.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(Neg.minus.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient2(low = "blue", mid="white",high = "red",na.value="white",limits = c(-limits,limits),name="Number")+
  scale_y_discrete(limit=Neg.minus_Yorder)+
  scale_x_discrete(limit=Neg.minus_Xorder)+
  theme(panel.background =element_blank(),#rect(colour="lightgray",fill = "white",size=2), 
        panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
        #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = negminus_ycolor),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        axis.line=element_blank())+
  geom_text(aes(label=Num),size=4)
dev.off()


pos.minus <- pos_cor.all.t[2:14] -pos_cor.all.n[2:14]
pos.minus$GeneSymbol <- pos_cor.all.t$Gene
pos.minus$GeneSymbol <- factor(pos.minus$GeneSymbol,levels = pos.minus$GeneSymbol)
pos.minus.m <- melt(pos.minus,id.vars = "GeneSymbol",measure.vars = colnames(pos.minus)[1:13])
colnames(pos.minus.m) <- c("GeneSymbol","Type","Num")
pos.t.index <- 1:13
pos.minus_Yorder <- pos.minus[order(apply(abs(pos.minus[,pos.t.index]), 1,sum)),"GeneSymbol"]
pos.minus_Xorder <- names(sort(apply(abs(pos.minus[,pos.t.index]), 2,sum)))
pos.minus_fill_limit <- c(min(pos.minus.m$Num),max(pos.minus.m$Num))
limits <- max(abs(pos.minus_fill_limit))
posminus_ycolor <- rep("black",times=length(pos.minus_Yorder))
posminus_ycolor[pos.minus_Yorder %in% core.circadian] <- "red"

pdf("pos_cor.all.tumor.minus.normal.pearson.R0.5_bluered.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(pos.minus.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient2(low = "blue", mid="white",high = "red",na.value="white",limits = c(-limits,limits),name="Number")+
  scale_y_discrete(limit=pos.minus_Yorder)+
  scale_x_discrete(limit=pos.minus_Xorder)+
  theme(panel.background =element_blank(),#rect(colour="lightgray",fill = "white",size=2), 
        panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
        #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = posminus_ycolor),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        axis.line=element_blank())+
  geom_text(aes(label=Num),size=4)
dev.off()
###only consider core circadian genes
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
Neg.minus <- neg_cor.all.t[2:14] -neg_cor.all.n[2:14]
Neg.minus$GeneSymbol <- neg_cor.all.t$Gene
Neg.minus <- Neg.minus[(Neg.minus$GeneSymbol %in% core.circadian),]
Neg.minus$GeneSymbol <- factor(Neg.minus$GeneSymbol,levels = Neg.minus$GeneSymbol)
Neg.minus.m <- melt(Neg.minus,id.vars = "GeneSymbol",measure.vars = colnames(Neg.minus)[1:13])
colnames(Neg.minus.m) <- c("GeneSymbol","Type","Num")
neg.t.index <- 1:13
Neg.minus_Yorder <- Neg.minus[order(apply(abs(Neg.minus[,neg.t.index]), 1,sum)),"GeneSymbol"]
Neg.minus_Xorder <- names(sort(apply(abs(Neg.minus[,neg.t.index]), 2,sum)))
Neg.minus_fill_limit <- c(min(Neg.minus.m$Num),max(Neg.minus.m$Num))
limits <- max(abs(Neg.minus_fill_limit))
pdf("core.neg_cor.all.tumor.minus.normal.pearson.R0.5_bluered.pdf",width=9,height = 5.5)#,width=4000,height=1500,res=400)
ggplot(Neg.minus.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient2(low = "blue", mid="white",high = "red",na.value="white",limits = c(-limits,limits),name="Number")+
  scale_y_discrete(limit=Neg.minus_Yorder)+
  scale_x_discrete(limit=Neg.minus_Xorder)+
  theme(panel.background =element_blank(),#rect(colour="lightgray",fill = "white",size=2), 
        panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
        #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        axis.line=element_blank())+
  geom_text(aes(label=Num),size=4)
dev.off()


pos.minus <- pos_cor.all.t[2:14] -pos_cor.all.n[2:14]
pos.minus$GeneSymbol <- pos_cor.all.t$Gene
pos.minus <- pos.minus[(pos.minus$GeneSymbol %in% core.circadian),]
pos.minus$GeneSymbol <- factor(pos.minus$GeneSymbol,levels = pos.minus$GeneSymbol)
pos.minus.m <- melt(pos.minus,id.vars = "GeneSymbol",measure.vars = colnames(pos.minus)[1:13])
colnames(pos.minus.m) <- c("GeneSymbol","Type","Num")
pos.t.index <- 1:13
pos.minus_Yorder <- pos.minus[order(apply(abs(pos.minus[,pos.t.index]), 1,sum)),"GeneSymbol"]
pos.minus_Xorder <- names(sort(apply(abs(pos.minus[,pos.t.index]), 2,sum)))
pos.minus_fill_limit <- c(min(pos.minus.m$Num),max(pos.minus.m$Num))
limits <- max(abs(pos.minus_fill_limit))
pdf("core.pos_cor.all.tumor.minus.normal.pearson.R0.5_bluered.pdf",width=9,height = 5.5)
ggplot(pos.minus.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient2(low = "blue", mid="white",high = "red",na.value="white",limits = c(-limits,limits),name="Number")+
  scale_y_discrete(limit=pos.minus_Yorder)+
  scale_x_discrete(limit=pos.minus_Xorder)+
  theme(panel.background =element_blank(),#rect(colour="lightgray",fill = "white",size=2), 
        panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
        #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        axis.line=element_blank())+
  geom_text(aes(label=Num),size=4)
dev.off()


###percentage
pos_max_cor <- pos_cor.all.t[1:22,2:14]
pos_max_corn <- pos_cor.all.n[1:22,2:14]
pos.minus.per <- pos.minus
for(i in 1:nrow(pos_max_cor)){
  for(j in 1:ncol(pos_max_cor)){
    pos_max_cor[i,j] <- ifelse(pos_max_cor[i,j] > pos_max_corn[i,j],pos_max_cor[i,j],pos_max_corn[i,j])
    pos.minus.per[i,j] <- signif(pos.minus[i,j]/pos_max_cor[i,j]*100,digits = 2)
  }
}
pos.minus.per.m <- melt(pos.minus.per,id.vars = "GeneSymbol",measure.vars = colnames(pos.minus.per)[1:13])
colnames(pos.minus.per.m) <- c("GeneSymbol","Type","Num")
pos.t.index <- 1:13
pos.minus.per_Yorder <- pos.minus.per[order(apply(abs(pos.minus.per[,pos.t.index]), 1,sum)),"GeneSymbol"]
pos.minus.per_Xorder <- names(sort(apply(abs(pos.minus.per[,pos.t.index]), 2,sum)))
pos.minus.per_fill_limit <- c(min(pos.minus.per.m$Num),max(pos.minus.per.m$Num))
limits <- max(abs(pos.minus.per_fill_limit))
pdf("core.pos_cor.all.tumor.minus.normal.pearson.R0.5_percentage_greenred.pdf",width=10,height = 8)#,width=4000,height=1500,res=400)
ggplot(pos.minus.per.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient2(low = "green", mid="white",high = "red",na.value="white",limits = c(-limits,limits),name="Percentage")+
  scale_y_discrete(limit=pos.minus.per_Yorder)+
  scale_x_discrete(limit=pos.minus.per_Xorder)+
  theme(panel.background =element_blank(),#rect(colour="lightgray",fill = "white",size=2), 
        panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
        #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        axis.line=element_blank())+
  geom_text(aes(label=Num),size=4)
dev.off()

Neg_max_cor <- neg_cor.all.t[1:22,2:14]
Neg_max_corn <- neg_cor.all.n[1:22,2:14]
Neg.minus.per <- Neg.minus
for(i in 1:nrow(Neg_max_cor)){
  for(j in 1:ncol(Neg_max_cor)){
    Neg_max_cor[i,j] <- ifelse(Neg_max_cor[i,j] > Neg_max_corn[i,j],Neg_max_cor[i,j],Neg_max_corn[i,j])
    Neg.minus.per[i,j] <- signif(Neg.minus[i,j]/Neg_max_cor[i,j]*100,digits = 2)
  }
}
Neg.minus.per.m <- melt(Neg.minus.per,id.vars = "GeneSymbol",measure.vars = colnames(Neg.minus.per)[1:13])
colnames(Neg.minus.per.m) <- c("GeneSymbol","Type","Num")
Neg.t.index <- 1:13
Neg.minus.per_Yorder <- Neg.minus.per[order(apply(abs(Neg.minus.per[,Neg.t.index]), 1,sum)),"GeneSymbol"]
Neg.minus.per_Xorder <- names(sort(apply(abs(Neg.minus.per[,Neg.t.index]), 2,sum)))
Neg.minus.per_fill_limit <- c(min(Neg.minus.per.m$Num),max(Neg.minus.per.m$Num))
limits <- max(abs(Neg.minus.per_fill_limit))
pdf("core.Neg_cor.all.tumor.minus.normal.pearson.R0.5_percentage_greenred.pdf",width=10,height = 8)#,width=4000,height=1500,res=400)
ggplot(Neg.minus.per.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient2(low = "green", mid="white",high = "red",na.value="white",limits = c(-limits,limits),name="Percentage")+
  scale_y_discrete(limit=Neg.minus.per_Yorder)+
  scale_x_discrete(limit=Neg.minus.per_Xorder)+
  theme(panel.background =element_blank(),#rect(colour="lightgray",fill = "white",size=2), 
        panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
        #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(2,"cm"),
        legend.key.width = unit(0.3,"cm"),
        axis.line=element_blank())+
  geom_text(aes(label=Num),size=4)
dev.off()

pdf("pos_cor.all.t.pearson.R0.5_as_percentage_order.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(pos_cor.all.t.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=as.numeric(Num)),col="white")+
  scale_fill_gradient(low = "white",high = "red",na.value="white",limit=pos_fill_limit1,name="Number")+
  scale_y_discrete(limit=pos.minus.per_Yorder)+
  scale_x_discrete(limit=pos.minus.per_Xorder)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_rect(colour="black",fill=NA,size=2),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = "black",face = "bold"),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5,face = "bold"),
    axis.ticks=element_line("black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank(),
    legend.key = element_rect(colour = "black"))+
  geom_text(aes(label=Num),size=6)
dev.off()

pdf("neg_cor.all.t.pearson.R0.5_as_percentage_order.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(neg_cor.all.t.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient(low = "white",high = "blue",na.value="white",limits=neg_fill_limit1,name="Number")+
  scale_y_discrete(limit=Neg.minus.per_Yorder)+
  scale_x_discrete(limit=Neg.minus.per_Xorder)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_rect(colour="black",fill=NA,size=2),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = "black",face = "bold"),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5,face = "bold"),
    axis.ticks=element_line(color="black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank())+
  geom_text(aes(label=Num),size=6)
dev.off()

pdf("pos_cor.all.n.pearson.R0.5_as_percentage_order.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(pos_cor.all.n.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient(low = "white",high = "red",na.value="white",limits=pos_fill_limit1,name="Number")+
  scale_y_discrete(limit=pos.minus.per_Yorder)+
  scale_x_discrete(limit=pos.minus.per_Xorder)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_rect(colour="black",fill=NA,size=2),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = "black",face = "bold"),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5,face = "bold"),
    axis.ticks=element_line(color="black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank())+
  geom_text(aes(label=Num),size=6)
dev.off()
pdf("neg_cor.all.n.pearson.R0.5_as_percentage_order.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(neg_cor.all.n.m,aes(x=Type,y=GeneSymbol))+
  geom_tile(aes(fill=Num),col="white")+
  scale_fill_gradient(low = "white",high = "blue",na.value="white",limits=neg_fill_limit1,name="Number")+
  scale_y_discrete(limit=Neg.minus.per_Yorder)+
  scale_x_discrete(limit=Neg.minus.per_Xorder)+
  theme(#panel.grid=element_blank(),#line(colour="grey",linetype="dashed"),
    #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
    panel.background=element_rect(colour="black",fill=NA,size=2),
    axis.title=element_blank(),
    axis.text.y=element_text(size=16,colour = "black",face = "bold"),
    axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5,face = "bold"),
    axis.ticks=element_line(color="black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.key = element_rect(colour = "black"),
    legend.key.height = unit(2,"cm"),
    legend.key.width = unit(0.3,"cm"),
    axis.line=element_blank())+
  geom_text(aes(label=Num),size=6)
dev.off()

core.n <- cor_data.n[which(cor_data.n$gene =="ARNTL"),c("CRY2",colnames(cor_data.t)[match("CRY2", colnames(cor_data.n))+1])]
core.t <- cor_data.t[which(cor_data.t$gene %in% core.circadian),c(core.circadian,colnames(cor_data.t)[match(core.circadian, colnames(cor_data.t))+1])]
######core circadian genes correlation in pair tumor
for(i in 1:length(corfiles.t)){
  cor_data.n <- read.delim(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/",corfiles.n[i],sep = ""),header=T)
  cor_data.t <- read.delim(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/",corfiles.t[i],sep = ""),header=T)
  core.n <- cor_data.n[which(cor_data.n$gene %in% core.circadian),core.circadian]
  core.t <- cor_data.t[which(cor_data.t$gene %in% core.circadian),core.circadian]
  rownames(core.t) <- cor_data.t[which(cor_data.t$gene %in% core.circadian),"gene"]
  rownames(core.n) <- cor_data.n[which(cor_data.n$gene %in% core.circadian),"gene"]
  
  upper.core.t <- core.t
  upper.core.t <- upper.core.t[order(rownames(upper.core.t)),order(colnames(upper.core.t))]
  upper.core.t[lower.tri(upper.core.t,diag = T)] <- NA
  upper.core.t$gene <- cor_data.t[which(cor_data.t$gene %in% core.circadian),"gene"]
  upper.core.t$gene <- upper.core.t$gene[order(upper.core.t$gene)]
  upper.core.tM <- melt(upper.core.t,id.vars = "gene",measure.vars = core.circadian,na.rm=T)
  
  lower.core.n <- core.n
  lower.core.n <- lower.core.n[order(rownames(lower.core.n)),order(colnames(lower.core.n))]
  
  lower.core.n[upper.tri(lower.core.n,diag = T)] <- NA
  lower.core.n$gene <- cor_data.n[which(cor_data.n$gene %in% core.circadian),"gene"]
  lower.core.n$gene <- lower.core.n$gene[order(lower.core.n$gene)]
  
  lower.core.nM <- melt(lower.core.n,id.vars = "gene",measure.vars = core.circadian,na.rm=T)
  core <- rbind(upper.core.tM,lower.core.nM)
  
  minus <- core.t-core.n
  minus$gene <- rownames(minus)
  minus.M <- melt(minus,id.vars="gene",measure.vars=core.circadian,na.rm=T)
  colnames(minus.M) <- c("gene","variable","minus")
  upper.core.tM <- merge(upper.core.tM,minus.M,by=c("gene","variable"))
  core.minus <- merge(core,minus.M,by=c("gene","variable"))
  
  core.minus <- merge(upper.core.tM,lower.core.nM,by.x=c("gene","variable"),by.y=c("variable","gene"))
 # core.minus <- core.minus[which(abs(core.minus$minus) >=0.2 & (abs(core.minus$value.x)  >= 0.3 | abs(core.minus$value.y) >= 0.3)),]
  core.minus.m <- melt(core.minus,id.vars = c("gene","variable"),measure.vars = c("value.x","value.y"))
  core.minus.t <- core.minus[,c(1:4)]
  core.minus.n <- core.minus[,c(2,1,5,4)]
  colnames(core.minus.n) <- c("gene","variable","value","minus")
  colnames(core.minus.t) <- c("gene","variable","value","minus")
  core.minus.mm <- rbind(core.minus.t,core.minus.n)
  #core.f <- core.minus.mm #core.minus.mm[which(abs(core.minus.mm$minus) >= 0.2),]
  #core.minus.mm["value"][abs(core.minus.mm["minus"]) < 0.2] <- 0
  mis <- data.frame(gene=rep(core.circadian,times=14),variable=rep(core.circadian,each=14))
  core.f <- merge(mis,core.minus.mm,by=c("gene","variable"),all.x=T)
  core.f["value"][core.f["value"] > 0.8] <- 0.8
  core.f["value"][core.f["value"] < -0.8] <- -0.8
  core.f$positive <- rep(0,times=nrow(core.f))
  core.f$negative <- rep(0,times=nrow(core.f))
  core.f["positive"][core.f["value"] > 0.3] <- 1
  core.f["negative"][core.f["value"] < -0.3] <- -1
  if(i==1){
    All.coreCor <- core.f[,c("gene","variable","positive","negative")]
  }else{
    All.coreCor$positive <- All.coreCor$positive + core.f$positive
    All.coreCor$negative <- All.coreCor$negative + core.f$negative
  }
  
  #upper.core.tM$cor <- rep(0,times=nrow(upper.core.tM))
  #upper.core.tM["cor"][upper.core.tM["value"] > 0] <- 1
  pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/core.circadian.cor.change/All/",substr(corfiles.t[i],1,4),"_core_circadian.genes.correlation.pdf",sep=""),width=8,height = 8)#,width=4000,height=1500,res=400)
  p<-ggplot(core.f,aes(x=gene,y=variable))+
    geom_tile(aes(fill=value),col="black")+
    scale_fill_gradient2(limit=c(-0.8,0.8),low = "blue",mid="white",high = "red",midpoint = 0,na.value="white",name="Cor Coeff.")+
    scale_y_discrete(limit=upper.core.t$gene,expand = c(0.002,0.002))+
    scale_x_discrete(limit=colnames(upper.core.t)[1:14],expand = c(0.002,0.002))+coord_fixed(ratio = 1)+
  #  geom_abline(slope=-1,intercept=1)+
    theme(panel.background=element_rect(colour="black",fill="white"),
          panel.grid=element_line(colour=NA),
          panel.grid.major=element_line(colour=NA),
          axis.title=element_blank(),
          axis.text.y=element_text(size=16,colour = "black"),
          axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5),
          axis.ticks=element_line(color="black"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14),
          legend.key = element_rect(colour = "black"),
          legend.key.height = unit(0.8,"cm"),
          legend.key.width = unit(2,"cm"),
          legend.position="bottom",
          legend.direction = "horizontal")
  print(p)
  dev.off()
  
}
###
All.coreCor.dup <- All.coreCor[which(All.coreCor$positive != 0 & All.coreCor$negative != 0),]
All.coreCor.uniq <- All.coreCor[!(rownames(All.coreCor) %in% rownames(All.coreCor.dup)),]
All.coreCor$merge <- rep(0,times=nrow(All.coreCor))
All.coreCor[rownames(All.coreCor.dup),"merge"] <- All.coreCor.dup$positive
All.coreCor[rownames(All.coreCor.uniq),"merge"] <- All.coreCor.uniq$positive + All.coreCor.uniq$negative
All.coreCor.uniq$negative <- abs(All.coreCor.uniq$negative)
All.coreCor.uniq[3:4][All.coreCor.uniq[3:4]==0] <- ""
#####expression correlation 
sss_x_label <- colnames(upper.core.t)[1:14]
sss_y_label <- upper.core.t$gene
All.coreCor.dup$x = sapply(All.coreCor.dup$gene,function(x){match(x,sss_x_label)}) ## add the geom_polygon ID, convert the factor to numeric 
All.coreCor.dup$y = sapply(All.coreCor.dup$variable,function(x){match(x,sss_y_label)})#as.numeric(dup.data$geneSymbol)
All.coreCor.dup$negative <- abs(All.coreCor.dup$negative)


poly_x = c()
poly_y = c()
for ( i in 1: nrow(All.coreCor.dup)){
  poly_x = c(poly_x,All.coreCor.dup[i,"x"]-0.48,All.coreCor.dup[i,'x'],All.coreCor.dup[i,"x"],All.coreCor.dup[i,"x"]-0.48 )
  poly_y = c(poly_y,All.coreCor.dup[i,"y"]-0.48,All.coreCor.dup[i,"y"]-0.48,All.coreCor.dup[i,"y"]+0.48,All.coreCor.dup[i,"y"]+0.48)
}
polygonID_negative <- data.frame( group = rep(seq(1:nrow(All.coreCor.dup)), each = 4),poly_x, poly_y )
polygonID_negative$negative <- rep(All.coreCor.dup$negative,each=4)

pdf("/extraspace/yye1/analysis/Circadian/expression/correlation/core.circadian.cor.change/All/core circadian genes correlation.pdf",width=8,height = 8)#,width=4000,height=1500,res=400)
ggplot(All.coreCor,aes(x=gene,y=variable))+
  geom_tile(aes(fill=merge),col="black")+
  scale_fill_gradient2(limit=c(-13,13),low = "blue",mid="white",high = "red",midpoint = 0,na.value="white",name="Cor Coeff.")+
  scale_y_discrete(limit=upper.core.t$gene,expand = c(0.002,0.002))+
  scale_x_discrete(limit=colnames(upper.core.t)[1:14],expand = c(0.002,0.002))+coord_fixed(ratio = 1)+
  geom_polygon(data=polygonID_negative,aes(x=poly_x,y=poly_y,group=group,fill=negative))+
  #  geom_abline(slope=-1,intercept=1)+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid=element_line(colour=NA),
        panel.grid.major=element_line(colour=NA),
        axis.title=element_blank(),
        axis.text.y=element_text(size=16,colour = "black"),
        axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(0.8,"cm"),
        legend.key.width = unit(2,"cm"),
        legend.position="bottom",
        legend.direction = "horizontal")+
  geom_text(data = All.coreCor.uniq,aes(x=gene,y=variable,label=positive))+
  geom_text(data = All.coreCor.uniq,aes(x=gene,y=variable,label=abs(as.numeric(negative))))+
  geom_text(data = All.coreCor.dup,aes(x=x+0.25,y=y,label=positive))+
  geom_text(data = All.coreCor.dup,aes(x=x-0.25,y=y,label=abs(as.numeric(negative))))
dev.off()



ggplot(All.coreCor,aes(x=tumor,y=geneSymbol))+
  geom_tile(aes(fill=factor(value)),col="darkgray")+
  scale_fill_manual(limits=c("High","Low","Stage","Subtype"),values = c("red","blue","green","gold"),na.value="white",labels=c("High_level_Worse","Low_level_Worse","Stage","Subtype"),name="")+
  scale_y_discrete(limit=sss_y_label)+
  scale_x_discrete(limit=sss_x_label)+
  geom_polygon(data=polygonID_stage,aes(x=poly_x,y=poly_y,group=group),fill="green")+
  geom_polygon(data=polygonID_subtype,aes(x=poly_x,y=poly_y,group=group),fill="gold")+
  geom_polygon(data=polygonID_threeOverlap_subtype,aes(x=poly_x,y=poly_y,group=group),fill="gold")+
  geom_polygon(data=polygonID_threeOverlap_stage,aes(x=poly_x,y=poly_y,group=group),fill="green")+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=16,colour = sss_y_color),
        axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=14),legend.position="bottom",legend.direction="horizontal")



pair_sampleAll <- read.delim("/home/yye1/Circadian/expression/pair_sampleAll.txt",header=T)
for(m in unique(pair_sampleAll$type)){
  BLCA <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",m,"_mRNA_each_exp_20160513",sep=""),header=T)
  # BLCA <- BLCA[which(BLCA$gene %in% circadian.genes$V1),] ##get circadian genes
  ###consider tumor type contain Tumor and normal genes
  BLCA$gene <- data.frame(do.call(rbind,strsplit(as.character(BLCA$gene),"\\|")))$X1
  
  colnames(BLCA) <- gsub("\\.","\\-",colnames(BLCA))
  BLCA <- BLCA[which(BLCA$gene %in% core.circadian),c("gene",colnames(BLCA)[colnames(BLCA) %in% pair_sampleAll[which(pair_sampleAll$type==m),"barcode"]])]
  BLCA.m <- melt(BLCA,id.vars = "gene",measure.vars = colnames(BLCA)[2:ncol(BLCA)])
  if(m=="BLCA"){
    pair_sample.exp <- BLCA.m 
  }else{
    pair_sample.exp <- rbind(pair_sample.exp,BLCA.m)
  }
}
pair_sample.expAll <- merge(pair_sampleAll,pair_sample.exp,by.x=c("barcode"),by.y="variable")
pair_sample.expAll$barcode <- substr(pair_sample.expAll$barcode,1,15)
pair_sample.expAll <- pair_sample.expAll[!duplicated(pair_sample.expAll[,c("barcode","gene")]),]
pairs <- substr(pair_sample.expAll$barcode,1,12)[duplicated(substr(pair_sample.expAll$barcode,1,12))]
pair_sample.expAll <- pair_sample.expAll[which(substr(pair_sample.expAll$barcode,1,12) %in% pairs),]

lm.out1 <- lm(log2(pair_sample.expAll[which(pair_sample.expAll$class=="tumor" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="PER3"),"value"])~log2(pair_sample.expAll[which(pair_sample.expAll$class=="tumor" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="ARNTL"),"value"])) 
abline(lm.out1, col="red")

plot(log2(pair_sample.expAll[which(pair_sample.expAll$class=="tumor" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="PER2"),"value"]),log2(pair_sample.expAll[which(pair_sample.expAll$class=="tumor" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="CRY2"),"value"]),xlim=c(5,13),ylim=c(5,13),xlab="PER3",ylab="CRY1",col="red")
par(new=T)
plot(log2(pair_sample.expAll[which(pair_sample.expAll$class=="normal" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="PER2"),"value"]),log2(pair_sample.expAll[which(pair_sample.expAll$class=="normal" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="CRY2"),"value"]),col="blue",xlim=c(5,13),ylim=c(5,13),xlab=NA,ylab=NA)
lm.out2 <- lm(log2(pair_sample.expAll[which(pair_sample.expAll$class=="normal" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="CRY1"),"value"])~log2(pair_sample.expAll[which(pair_sample.expAll$class=="normal" & pair_sample.expAll$type=="BRCA" & pair_sample.expAll$gene=="PER3"),"value"])) 
abline(lm.out2, col="blue")
dev.off()

####PER3 with "ARNTL","NPAS2","NR1D1","NR1D2
LUSC_pair_sample.exp <- pair_sample.expAll[which(pair_sample.expAll$type=="LUSC"),]
LUSC_pair_sample.exp_other <- Reduce(function(x,y){merge(x,y,by="barcode")},list(LUSC_pair_sample.exp[which(LUSC_pair_sample.exp$gene=="ARNTL"),],
                                                                                 LUSC_pair_sample.exp[which(LUSC_pair_sample.exp$gene=="NPAS2"),c(1,5)],
                                                                                 LUSC_pair_sample.exp[which(LUSC_pair_sample.exp$gene=="NR1D1"),c(1,5)],
                                                                                 LUSC_pair_sample.exp[which(LUSC_pair_sample.exp$gene=="NR1D2"),c(1,5)]))

colnames(LUSC_pair_sample.exp_other) <- c("barcode","class","type","gene","ARNTL","NPAS2","NR1D1","NR1D2")
LUSC_pair_sample.exp_other.m <- melt(LUSC_pair_sample.exp_other,id.vars=colnames(LUSC_pair_sample.exp_other)[1:3],measure.vars=colnames(LUSC_pair_sample.exp_other)[5:8])
LUSC_pair_sample.expM <- merge(LUSC_pair_sample.exp[which(LUSC_pair_sample.exp$gene=="PER3"),],LUSC_pair_sample.exp_other.m,by=c("barcode","class","type"))

colnames(LUSC_pair_sample.expM) <- c("barcode","class","type","PER3","PER3_RPKM","Other","Other_RPKM")
core.n_cor <- signif(core.n[which(rownames(core.n) %in% c("ARNTL","NPAS2","NR1D1","NR1D2")),"PER3"],digits = 2)
core.t_cor <- signif(core.t[which(rownames(core.t) %in% c("ARNTL","NPAS2","NR1D1","NR1D2")),"PER3"],digits = 2)
core.n_cor <- data.frame(core.n_cor)
core.t_cor <- data.frame(core.t_cor)
core.n_cor$Other <- c("ARNTL","NPAS2","NR1D1","NR1D2")
core.t_cor$Other <- c("ARNTL","NPAS2","NR1D1","NR1D2")


pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/core.circadian.cor.change/PER3_other.expressioncorrelation.pdf",sep=""),width=10,height=3)
ggplot(LUSC_pair_sample.expM,aes(x=log2(Other_RPKM),y=log2(PER3_RPKM)))+
  geom_point(aes(color=factor(class)))+
  geom_text(data = core.t_cor,aes(y = 5.7,x=12.2, label = paste("R = ",core.t_cor,sep="")),color="red",size=4,hjust=1)+
  geom_text(data = core.n_cor,aes(y = 5.1,x=12.2, label = paste("R = ",core.n_cor,sep="")),color="blue",size=4,hjust=1)+
  facet_wrap(~Other,ncol=4)+
  scale_color_manual(limits=c("normal","tumor"),values=c("blue","red"),labels=c("Normal","Tumor"),name="")+
 # stat_smooth(method = "lm",color="red",size=1)+
  labs(x="Gene Expression",
       y="PER3 Gene Expression")+
  #annotate("text",x=-2,y=12,label=paste("n=",NB_num,sep=""),size=8,color="darkorange2")+
  # annotate("text",x=10.1,y=-2,hjust=0,label=paste("n=",n_num,sep=""),size=8,color="steelblue2")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=20),
        strip.background = element_rect(fill=NA),legend.background = element_rect(fill=NA),legend.text=element_text(size=14))
dev.off()

###KIRP CRY2 with "ARNTL2","ARNTL","NPAS2","RORB","RORC","PER1","PER2","PER3"
KIRP_pair_sample.exp <- pair_sample.expAll[which(pair_sample.expAll$type=="KIRP"),]
KIRP_pair_sample.exp_other <- Reduce(function(x,y){merge(x,y,by="barcode")},list(KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="ARNTL2"),],
                                                                                 KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="ARNTL"),c(1,5)],
                                                                                 KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="NPAS2"),c(1,5)],
                                                                                 KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="RORB"),c(1,5)],
                                                                                 KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="RORC"),c(1,5)],
                                                                                 KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="PER1"),c(1,5)],
                                                                                 KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="PER2"),c(1,5)],
                                                                                 KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="PER3"),c(1,5)]))

colnames(KIRP_pair_sample.exp_other) <- c("barcode","class","type","gene","ARNTL2","ARNTL","NPAS2","RORB","RORC","PER1","PER2","PER3")
KIRP_pair_sample.exp_other.m <- melt(KIRP_pair_sample.exp_other,id.vars=colnames(KIRP_pair_sample.exp_other)[1:3],measure.vars=colnames(KIRP_pair_sample.exp_other)[5:12])
KIRP_pair_sample.expM <- merge(KIRP_pair_sample.exp[which(KIRP_pair_sample.exp$gene=="CRY2"),],KIRP_pair_sample.exp_other.m,by=c("barcode","class","type"))

colnames(KIRP_pair_sample.expM) <- c("barcode","class","type","CRY2","CRY2_RPKM","Other","Other_RPKM")
core.n$Other <- rownames(core.n)
core.t$Other <- rownames(core.t)

core.n_cor <- core.n[which(rownames(core.n) %in% c("ARNTL2","ARNTL","NPAS2","RORB","RORC","PER1","PER2","PER3")),c("CRY2","Other")]
core.t_cor <- core.t[which(rownames(core.t) %in% c("ARNTL2","ARNTL","NPAS2","RORB","RORC","PER1","PER2","PER3")),c("CRY2","Other")]
colnames(core.t_cor) <- c("core.t_cor","Other")
colnames(core.n_cor) <- c("core.n_cor","Other")
core.t_cor$core.t_cor <- signif(core.t_cor$core.t_cor, digits = 2)
core.n_cor$core.n_cor <- signif(core.n_cor$core.n_cor, digits = 2)
KIRP_pair_sample.expM$Other_RPKM <- log2(as.numeric(KIRP_pair_sample.expM$Other_RPKM)+1)
KIRP_pair_sample.expM$CRY2_RPKM <- log2(as.numeric(KIRP_pair_sample.expM$CRY2_RPKM)+1)

pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/core.circadian.cor.change/KIRP_CRY2_other.expressioncorrelation.pdf",sep=""),width=9,height=6)
ggplot(KIRP_pair_sample.expM,aes(x=Other_RPKM,y=CRY2_RPKM))+
  geom_point(aes(color=factor(class)))+
  geom_text(data = core.t_cor,aes(y = 12.5,x=0.2, label = paste("R = ",core.t_cor,sep="")),color="red",size=4,hjust=0)+
  geom_text(data = core.n_cor,aes(y = 12.1,x=0.2, label = paste("R = ",core.n_cor,sep="")),color="blue",size=4,hjust=0)+
 # scale_x_continuous(breaks=seq(0,12,by=4))+
  facet_wrap(~Other,ncol=4)+
  scale_color_manual(limits=c("normal","tumor"),values=c("blue","red"),labels=c("Normal","Tumor"),name="")+
  # stat_smooth(method = "lm",color="red",size=1)+
  labs(x="Gene Expression",
       y="CRY2 Gene Expression")+
  #annotate("text",x=-2,y=12,label=paste("n=",NB_num,sep=""),size=8,color="darkorange2")+
  # annotate("text",x=10.1,y=-2,hjust=0,label=paste("n=",n_num,sep=""),size=8,color="steelblue2")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=20),
        strip.background = element_rect(fill=NA),legend.background = element_rect(fill=NA),legend.text=element_text(size=14),
        legend.position = "bottom",legend.direction="horizontal",legend.key=element_blank())
dev.off()



#########
setwd("/extraspace/yye1/analysis/Circadian/expression/correlation")
##circadian.genes.tumor
corfiles.t <- list.files(pattern = "*t_pearson.txt")
corfiles.n <- list.files(pattern = "*n_pearson.txt")
corfiles_tumor <- gsub("_gene_cor.t_pearson.txt","",corfiles.t)
for(i in 1:length(corfiles.t)){
  cor_data.t <- read.delim(corfiles.t[i],sep = "",header=T)
  cor_data.n <- read.delim(corfiles.n[i],sep = "",header=T)
  cor_data.t <- cor_data.t[which(cor_data.t$gene != "?"),]
  cor_data.n <- cor_data.n[which(cor_data.n$gene != "?"),]
  cor.t <- cor_data.t[,grep("FDR",colnames(cor_data.t))-1]
  cor.t[cor.t> 0.5] <- 1
  cor.t[cor.t< -0.5] <- -1
  corFDR.t <- cor_data.t[,grep("FDR",colnames(cor_data.t))]
  corFDR.t[corFDR.t < 0.05] <- 0
  corsign.t <- cor.t+corFDR.t
  corsign.t[abs(corsign.t) != 1] <- 0
  cor.n <- cor_data.n[,grep("FDR",colnames(cor_data.n))-1]
  cor.n[cor.n> 0.5] <- 3
  cor.n[cor.n< -0.5] <- -3
  corFDR.n <- cor_data.n[,grep("FDR",colnames(cor_data.n))]
  corFDR.n[corFDR.n < 0.05] <- 0
  corsign.n <- cor.n+corFDR.n
  corsign.n[abs(corsign.n) != 3] <- 0
  corSum <- corsign.t + corsign.n
  corSum <- data.frame(gene = cor_data.t$gene,corSum)
  if(i==1){
    AllcorSum <- corSum
  }else{
    AllcorSum <-  cbind(AllcorSum,corSum[,2:ncol(corSum)])
  }
}
write.csv(AllcorSum,file="/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/AllcorSum.csv",quote = F,row.names = F)

setwd("/extraspace/yye1/analysis/Circadian/expression/correlation/")
##circadian.genes.tumor
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
library(gplots)
AllcorSum <- read.csv("/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/AllcorSum.csv",header=T)

library(gdata)
library(doParallel)
library(doMC)
registerDoMC()
library(foreach)
library(methods)

foreach(i = core.circadian) %dopar% {
  subCor <- AllcorSum[,c(1,grep(i,colnames(AllcorSum)))]
  colnames(subCor) <- c("gene",corfiles_tumor)
  negtspCor <- subCor
  negtspCor[2:ncol(negtspCor)][negtspCor[2:ncol(negtspCor)] != -1] <- 0
  negtspCor <- negtspCor[(apply(negtspCor[2:ncol(negtspCor)],1,sum)< 0),]
  negtspCorNum <- apply(negtspCor[,2:ncol(negtspCor)],2,sum)
  negtspCorNum <- sort(negtspCorNum)
  negtspCor.m <- negtspCor[do.call(order, as.data.frame(negtspCor[match(names(negtspCorNum),colnames(negtspCor))])),names(negtspCorNum)]
  rownames(negtspCor.m) <- negtspCor$gene
  pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/",i,".tumorspecific.gain.negative.correlation.pdf",sep=""),width = 6,height = 8)
  heatmap.2(t(as.matrix(negtspCor.m)),Colv = NA,Rowv = NA,col=c("blue","gray"),breaks=c(-0.9,-0.1,0),trace="none",labCol = NA,scale = "none",key=F)
  dev.off()
  
  postspCor <- subCor
  postspCor[2:ncol(postspCor)][postspCor[2:ncol(postspCor)] != 1] <- 0
  postspCor <- postspCor[(apply(postspCor[2:ncol(postspCor)],1,sum) > 0),]
  postspCorNum <- apply(postspCor[,2:ncol(postspCor)],2,sum)
  postspCorNum <- rev(sort(postspCorNum))
  postspCor.m <- postspCor[do.call(order, as.data.frame(postspCor[match(names(postspCorNum),colnames(postspCor))])),names(postspCorNum)]
  rownames(postspCor.m) <- postspCor$gene
  postspCor.m <- postspCor.m[nrow(postspCor.m):1,]
  
  pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/",i,".tumorspecific.gain.postive.correlation.pdf",sep=""),width = 6,height = 8)
  heatmap.2(t(as.matrix(postspCor.m)),Colv = NA,Rowv = NA,col=c("gray","red"),breaks=c(0,0.9,1),trace="none",labCol = NA,scale = "none",key=F)
  dev.off()
  
  negnspCor <- subCor
  negnspCor[2:ncol(negnspCor)][negnspCor[2:ncol(negnspCor)] != -3] <- 0
  negnspCor <- negnspCor[(apply(negnspCor[2:ncol(negnspCor)],1,sum)< 0),]
  negnspCorNum <- apply(negnspCor[,2:ncol(negnspCor)],2,sum)
  negnspCorNum <- sort(negnspCorNum)
  negnspCor.m <- negnspCor[do.call(order, as.data.frame(negnspCor[match(names(negnspCorNum),colnames(negnspCor))])),names(negnspCorNum)]
  rownames(negnspCor.m) <- negnspCor$gene
  pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/",i,".tumorspecific.loss.negative.correlation.pdf",sep=""),width = 6,height = 8)
  heatmap.2(t(as.matrix(negnspCor.m)),Colv = NA,Rowv = NA,col=c("blue","gray"),breaks=c(-0.29,-0.1,0),trace="none",labCol = NA,scale = "none",key=F)
  dev.off()
  
  posnspCor <- subCor
  posnspCor[2:ncol(posnspCor)][posnspCor[2:ncol(posnspCor)] != 3] <- 0
  posnspCor <- posnspCor[(apply(posnspCor[2:ncol(posnspCor)],1,sum) > 0),]
  posnspCorNum <- apply(posnspCor[,2:ncol(posnspCor)],2,sum)
  posnspCorNum <- rev(sort(posnspCorNum))
  posnspCor.m <- posnspCor[do.call(order, as.data.frame(posnspCor[match(names(posnspCorNum),colnames(posnspCor))])),names(posnspCorNum)]
  rownames(posnspCor.m) <- posnspCor$gene
  posnspCor.m <- posnspCor.m[nrow(posnspCor.m):1,]
  
  pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/",i,".tumorspecific.loss.postive.correlation.pdf",sep=""),width = 6,height = 8)
  heatmap.2(t(as.matrix(posnspCor.m)),Colv = NA,Rowv = NA,col=c("gray","red"),breaks=c(0,0.29,3),trace="none",labCol = NA,scale = "none",key=F)
  dev.off()
  
  negcomCor <- subCor
  negcomCor[2:ncol(negcomCor)][negcomCor[2:ncol(negcomCor)] != -4] <- 0
  negcomCor <- negcomCor[(apply(negcomCor[2:ncol(negcomCor)],1,sum)< 0),]
  negcomCorNum <- apply(negcomCor[,2:ncol(negcomCor)],2,sum)
  negcomCorNum <- sort(negcomCorNum)
  negcomCor.m <- negcomCor[do.call(order, as.data.frame(negcomCor[match(names(negcomCorNum),colnames(negcomCor))])),names(negcomCorNum)]
  rownames(negcomCor.m) <- negcomCor$gene
  pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/",i,".tumors_normal.negative.correlation.pdf",sep=""),width = 6,height = 8)
  heatmap.2(t(as.matrix(negcomCor.m)),Colv = NA,Rowv = NA,col=c("blue","gray"),breaks=c(-0.39,-0.3,0),trace="none",labCol = NA,scale = "none",key=F)
  dev.off()
  
  poscomCor <- subCor
  poscomCor[2:ncol(poscomCor)][poscomCor[2:ncol(poscomCor)] != 4] <- 0
  poscomCor <- poscomCor[(apply(poscomCor[2:ncol(poscomCor)],1,sum) > 0),]
  poscomCorNum <- apply(poscomCor[,2:ncol(poscomCor)],2,sum)
  poscomCorNum <- rev(sort(poscomCorNum))
  poscomCor.m <- poscomCor[do.call(order, as.data.frame(poscomCor[match(names(poscomCorNum),colnames(poscomCor))])),names(poscomCorNum)]
  rownames(poscomCor.m) <- poscomCor$gene
  poscomCor.m <- poscomCor.m[nrow(poscomCor.m):1,]
  pdf(paste("/extraspace/yye1/analysis/Circadian/expression/correlation/correlationGenes/",i,".tumor_normal.postive.correlation.pdf",sep=""),width = 6,height = 8)
  heatmap.2(t(as.matrix(poscomCor.m)),Colv = NA,Rowv = NA,col=c("gray","red"),breaks=c(0,0.39,4),trace="none",labCol = NA,scale = "none",key=F)
  dev.off()
  
}


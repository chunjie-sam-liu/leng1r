setwd("~/Circadian/expression/correlation/")
library(ggplot2)
library(reshape2)
files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
files.Abs <- gsub("_mRNA_each_exp_20160513","",files.names)
circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
spearman_cor <- function(x,y){cor(x, y, use = "everything",method = "spearman")}
sample.num <- c("Tumor.type","Number")
for(m in 1:length(files.names)){
  BLCA <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",files.Abs[m],"_mRNA_each_exp_20160513",sep=""),header=T)
  BLCA$gene <- data.frame(do.call(rbind,strsplit(as.character(BLCA$gene),"\\|")))$X1
  ###consider tumor type contain Tumor and normal genes
  BLCA.names <- colnames(BLCA)[2:ncol(BLCA)] 
  BLCA.names <- BLCA.names[grep("TCGA",BLCA.names)]
  filter.names <- c("gene",BLCA.names[(as.numeric(substr(BLCA.names,14,15)) ==1)])
  filter.names <- filter.names[!is.na(filter.names)]
  if(length(filter.names)>10){
    count <- c(files.Abs[m],length(filter.names))
    sample.num <- rbind(sample.num,count)
    BLCA.filter <- BLCA[,filter.names]
    BLCA.cor <- data.frame(BLCA.filter$gene)
    colnames(BLCA.cor) <- "gene"
    for(i in circadian.genes$V1){
      cir.exp <- as.matrix(BLCA.filter[which(BLCA.filter$gene==i),2:ncol(BLCA.filter)])
      gene.cor <- apply( as.matrix(BLCA.filter[,2:ncol(BLCA.filter)]) , 1 ,cor, y = as.vector(as.numeric(cir.exp)))
      gene.cor <- data.frame(gene.cor)
      colnames(gene.cor) <- i
      gene.cor <- data.frame(gene.cor)
      colnames(gene.cor) <- i
      p_value <- apply(as.matrix(BLCA.filter[,2:ncol(BLCA.filter)]),1,function(x){cor.test(x,as.vector(as.numeric(cir.exp)))$p.value})
      FDR_value <- p.adjust(p_value,method="fdr")
      FDR_value <- data.frame(FDR_value)
      colnames(FDR_value) <- paste(i,"_FDR",sep = "")
      BLCA.cor <- cbind(BLCA.cor,gene.cor,FDR_value)
    }
    BLCA.cor[is.na(BLCA.cor)] <- 0
    write.table(BLCA.cor,file=paste("pearson/",files.Abs[m],"_gene_cor.pearson.txt",sep = ""),quote = F,row.names = F,sep="\t")
  }
}
setwd("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes")
actional.genes <- read.delim("~/Circadian/expression/correlation/pearson/CRY1_CRY2_cor0.3/actionable.genes.txt",header=T)
actional.genes$Gene <- as.character(actional.genes$Gene)
circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
corfiles <- list.files(path="../",pattern = "*_gene_cor.pearson.txt")
#pair_tumor <- read.delim("/home/yye1/Circadian/expression/tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)
#corfiles <- corfiles[gsub("_gene_cor.pearson.txt","",corfiles) %in% pair_tumor$Type]
corfiles_tumor <- gsub("_gene_cor.pearson.txt","",corfiles)
cir_names <- gsub("\\|",".",circadian.genes$V1)
neg_action_All <- actional.genes
for( i in cir_names){
  negGeneAll <- c()
  for(j in 1:length(corfiles)){
    cor_data <- read.delim(paste("../",corfiles[j],sep = ""),header=T)
    cor_data <- cor_data[colnames(cor_data)[c(1,grep(i,colnames(cor_data)))]]
    negGene <- as.character(cor_data[which(cor_data[,2] <= -0.3 & cor_data[,3] < 0.05),"gene"])
    negGene <- intersect(negGene,actional.genes$Gene)
    if(length(negGeneAll)==0){
      negGeneAll <- negGene
    }else{
      negGeneAll <- c(negGeneAll,negGene)
    }
    
  }
  neg_action.num <- table(as.factor(negGeneAll))
  if(length(neg_action.num) >=1){
    neg_action.num <- data.frame(neg_action.num)
    colnames(neg_action.num) <- c("Gene",i)
  }else{
    neg_action.num <- actional.genes
    neg_action.num[i] <- rep(0,times=nrow(actional.genes))
  }
  #write.table(neg_action,file=paste("all_cancer_type/",i,"overlap.with.actional.genes.txt",sep=""),quote = F,row.names = F,sep="\t")
  # neg_action_All[i] <- neg_action_All$Gene %in% neg_action
  neg_action_All <- merge(neg_action_All,neg_action.num,by="Gene",all.x=T)
}

########positive correlated to actional genes
pos_action_All <- actional.genes
for( i in cir_names){
  posGeneAll <- c()
  for(j in 1:length(corfiles)){
    cor_data <- read.delim(paste("../",corfiles[j],sep = ""),header=T)
    cor_data <- cor_data[colnames(cor_data)[c(1,grep(i,colnames(cor_data)))]]
    posGene <- as.character(cor_data[which(cor_data[,2] >= 0.3 & cor_data[,3] < 0.05),"gene"])
    posGene <- intersect(posGene,actional.genes$Gene)
    if(length(posGeneAll)==0){
      posGeneAll <- posGene
    }else{
      posGeneAll <- c(posGeneAll,posGene)
    }
    
  }
  pos_action.num <- table(as.factor(posGeneAll))
  if(length(pos_action.num) >=1){
    pos_action.num <- data.frame(pos_action.num)
    colnames(pos_action.num) <- c("Gene",i)
  }else{
    pos_action.num <- actional.genes
    pos_action.num[i] <- rep(0,times=nrow(actional.genes))
  }
  #write.table(pos_action,file=paste("all_cancer_type/",i,"overlap.with.actional.genes.txt",sep=""),quote = F,row.names = F,sep="\t")
  # pos_action_All[i] <- pos_action_All$Gene %in% pos_action
  pos_action_All <- merge(pos_action_All,pos_action.num,by="Gene",all.x=T)
  print(i)
}


action_All <- actional.genes
for( i in cir_names){
  for(j in 1:length(corfiles)){
    cor_data <- read.delim(paste("/home/yye1/Circadian/expression/correlation/pearson/",corfiles[j],sep = ""),header=T)
    cor_data <- cor_data[colnames(cor_data)[c(1,grep(i,colnames(cor_data)))]]
    posGene <- as.character(cor_data[which(cor_data[,2] >= 0.3 & cor_data[,3] < 0.05),"gene"])
    posGene <- intersect(posGene,actional.genes$Gene)
    negGene <- as.character(cor_data[which(cor_data[,2] <= -0.3 & cor_data[,3] < 0.05),"gene"])
    negGene <- intersect(negGene,actional.genes$Gene)
    pos_Ac_cir <- data.frame(Actionable = posGene,Type = rep(j,length(posGene)),Cir = rep(i,length(posGene)))
    neg_Ac_cir <- data.frame(Actionable = negGene,Type = rep(j,length(negGene)),Cir = rep(i,length(negGene)))
    if( j ==1){
      pos_Ac_cirA <- pos_Ac_cir
      neg_Ac_cirA <- neg_Ac_cir
    }else{
      pos_Ac_cirA <- rbind(pos_Ac_cirA,pos_Ac_cir)
      neg_Ac_cirA <- rbind(neg_Ac_cirA,neg_Ac_cir)
    }
  }
  if(i == cir_names[1]){
    pos_Ac_cirAll <- pos_Ac_cirA
    neg_Ac_cirAll <- neg_Ac_cirA
  }else{
    pos_Ac_cirAll <- rbind(pos_Ac_cirAll,pos_Ac_cirA)
    neg_Ac_cirAll <- rbind(neg_Ac_cirAll,neg_Ac_cirA)
  }
}





neg_action_All[is.na(neg_action_All)] <- 0
#colnames(neg_action_All) <-  data.frame(do.call(rbind, strsplit(as.character(colnames(neg_action_All)),'\\.')))$X1
pos_action_All[is.na(pos_action_All)] <- 0
#colnames(pos_action_All) <-  data.frame(do.call(rbind, strsplit(as.character(colnames(pos_action_All)),'\\.')))$X1


neg_action_All <- neg_action_All[which(apply(neg_action_All[2:ncol(neg_action_All)],1,sum)>=1),]
write.table(neg_action_All,file="all_cancer_type/circadian.correlated.genes_overlap.with.actional.genes.txt",quote = F,row.names = F,sep="\t")
neg_action_All <- read.delim("all_cancer_type/circadian.correlated.genes_overlap.with.actional.genes.txt",header=T)
#set in more than 5 cancer types
neg_action_All[2:52][neg_action_All[2:52] <= 4] <- 0
neg_action_All <- neg_action_All[which(apply(neg_action_All[2:52],1,sum)>=5),]
neg_action_All$Gene <- as.character(neg_action_All$Gene)
###
neg_action_All.m <- melt(neg_action_All,id.vars = "Gene",measure.vars = colnames(neg_action_All)[2:ncol(neg_action_All)])

neg_action_All.m["Overlap"] <- rep(0,times=nrow(neg_action_All.m))
neg_action_All.m["Overlap"][neg_action_All.m["value"] >=1 ] <- -1
colnames(neg_action_All.m) <- c("Actional","Circadian","Number","Overlap")
neg.ann <- neg_action_All.m[which(neg_action_All.m$Number >=1),]
neg_xlimit <- neg_action_All[order(apply(neg_action_All[,2:ncol(neg_action_All)], 1,sum)),"Gene"]
ac_neg_y <- names(rev(sort(apply(neg_action_All[,2:ncol(neg_action_All)], 2,sum))))
ac_neg_ycol <- rep("black",times=length(ac_neg_y))
ac_neg_ycol[ac_neg_y %in% core.circadian] <- "red"
pdf("all_cancer_type/circadian.negatively.correlated.genes_overlap.with.actional.genes1_5cancerTypes.pdf",width=18,height = 10)#,width=4000,height=1500,res=400)
ggplot(neg_action_All.m,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=factor(Overlap)),col="white")+
  scale_fill_manual(limit=c(0,-1),values = c("gray","lightblue"),labels=c("Non-Cor","Negatively-Cor"),name="",guide=F)+
  scale_x_discrete(limit=neg_action_All[order(apply(neg_action_All[,2:ncol(neg_action_All)], 1,sum)),"Gene"])+
  scale_y_discrete(limit=ac_neg_y)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = ac_neg_ycol),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(1,"cm"),
        legend.key.width = unit(1,"cm"))+
  annotate("text",y=neg.ann$Circadian,x=neg.ann$Actional,label=neg.ann$Number,size=4)
dev.off()

pos_action_All <- pos_action_All[which(apply(pos_action_All[2:ncol(pos_action_All)],1,sum)>=1),]
write.table(pos_action_All,file="all_cancer_type/circadian.positively.correlated.genes_overlap.with.actional.genes.txt",quote = F,row.names = F,sep="\t")
pos_action_All <- read.delim("all_cancer_type/circadian.positively.correlated.genes_overlap.with.actional.genes.txt",header = T)
#set in more than 5 cancer types
pos_action_All[2:52][pos_action_All[2:52] <= 4] <- 0

pos_action_All <- pos_action_All[which(apply(pos_action_All[2:52],1,sum)>=5),]
pos_action_All$Gene <- as.character(pos_action_All$Gene)
###
pos_action_All.m <- melt(pos_action_All,id.vars = "Gene",measure.vars = colnames(pos_action_All)[2:ncol(pos_action_All)])
pos_action_All.m["Overlap"] <- rep(0,times=nrow(pos_action_All.m))
pos_action_All.m["Overlap"][pos_action_All.m["value"] >=1 ] <- 1
colnames(pos_action_All.m) <- c("Actional","Circadian","Number","Overlap")
pos.ann <- pos_action_All.m[which(pos_action_All.m$Number >=1),]
pos_xlimit <- pos_action_All[order(apply(pos_action_All[,2:ncol(pos_action_All)], 1,sum)),"Gene"]
ac_pos_y <- names(rev(sort(apply(pos_action_All[,2:ncol(pos_action_All)], 2,sum))))
ac_pos_ycol <- rep("black",times=length(ac_pos_y))
ac_pos_ycol[ac_pos_y %in% core.circadian] <- "red"
pdf("all_cancer_type/circadian.positively.correlated.genes_overlap.with.actional.genes1_5cancerTypes.pdf",width=24,height = 10)#,width=4000,height=1500,res=400)
ggplot(pos_action_All.m,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=factor(Overlap)),col="white")+
  scale_fill_manual(limit=c(0,1),values = c("gray","lightpink"),labels=c("Non-Cor","Positively-Cor"),name="",guide=F)+
  scale_x_discrete(limit=pos_action_All[order(apply(pos_action_All[,2:ncol(pos_action_All)], 1,sum)),"Gene"])+
  scale_y_discrete(limit=ac_pos_y)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = ac_pos_ycol),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(1,"cm"),
        legend.key.width = unit(1,"cm"))+
  annotate("text",y=pos.ann$Circadian,x=pos.ann$Actional,label=pos.ann$Number,size=4)
dev.off()

#intersect(neg_action_All[which(neg_action_All$PER1 > 0),c("Gene","PER1")]$Gene,pos_action_All[which(pos_action_All$PER1 > 0),c("Gene","PER1")]$Gene)

neg_action_All.f <- neg_action_All.m[which(neg_action_All.m$Actional != "JAK2" & neg_action_All.m$Circadian != "FOXO1"),]
pos_action_All.f <- pos_action_All.m[which(pos_action_All.m$Actional != "JAK2" & pos_action_All.m$Circadian != "FOXO1"),]
non_cor <- merge(pos_action_All.m[which(pos_action_All.m$Number ==0),],neg_action_All.m[which(neg_action_All.m$Number == 0),],by=c("Actional","Circadian"),all=T)[1:4]
colnames(non_cor) <- c("Actional","Circadian","Number","Overlap")
non_cor[is.na(non_cor)] <- 0
action_All <- rbind(pos_action_All.m[which(pos_action_All.m$Number >=1),],neg_action_All.m[which(neg_action_All.m$Number >=1),])
o <- merge(non_cor,action_All,by=c("Actional","Circadian"))[1:4]
fi <- paste(o$Actional,o$Circadian,sep="")
non_cor$fi <- paste(non_cor$Actional,non_cor$Circadian,sep="")
non_cor.f <- non_cor[!(non_cor$fi %in% fi),][1:4]
action_All.f <- rbind(pos_action_All.m[which(pos_action_All.m$Number >=1),],neg_action_All.m[which(neg_action_All.m$Number >=1),],non_cor.f)
All_xlimit <- sapply(split(action_All.f[,"Overlap"],action_All.f$Actional),sum)
#All_xlimit <- All_xlimit[All_xlimit !=0]
All_xlimit_label <- names(All_xlimit)[order(All_xlimit)]
All_ylimit <- sapply(split(action_All.f[,"Number"],action_All.f$Circadian),sum)
#All_ylimit <- All_ylimit[All_ylimit !=0]
All_ylimit_label <- names(All_ylimit)[order(All_ylimit)]
All_ann <- action_All.f[which(action_All.f$Number >1),]
All_ylimit_col <- rep("black",times=length(All_ylimit_label))
All_ylimit_col[All_ylimit_label %in% core.circadian] <- "red"
pdf("all_cancer_type/51circadian.correlated.genes_overlap.with.actional.genes1_5cancerTypes.pdf",width=28,height = 12)#,width=4000,height=1500,res=400)
ggplot(action_All.f,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=factor(Overlap)),col="white")+
  scale_fill_manual(limit=c(-1,0,1),values = c("lightblue","gray","lightpink"),na.value="gray",labels=c("Negatively-Cor","Non-Cor","Positively-Cor"),name="")+
  scale_x_discrete(limit=All_xlimit_label)+
  scale_y_discrete(limit=All_ylimit_label)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = All_ylimit_col ),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=14),
        # legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(1,"cm"),
        legend.position = "bottom",
        legend.direction ="horizontal",
        legend.key.width = unit(1,"cm"))+
  annotate("text",y=All_ann$Circadian,x=All_ann$Actional,label=All_ann$Number,size=4)#+
# geom_polygon(data= data.frame(group=rep(1,times=4),poly_x=c()),aes(x=poly_x,y=poly_y,fill="lightblue", group=group)))
dev.off()
write.table(All_xlimit_label,file = "all_cancer_type/Actional.genes.correlated.circadian.txt",quote = F,sep="\t",row.names = F)

action_All.ff <- action_All.f
action_All.ff['Number'][action_All.ff["Overlap"]==-1] <- -action_All.ff['Number'][action_All.ff["Overlap"]==-1]
write.table(action_All.ff,file="/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/action_All.ff.txt",sep="\t",row.names = F,quote = F)
pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/51circadian.correlated.genes_overlap.with.actional.genes1_5cancerTypes_colorful_withoutann_sigPPI.pdf",width=28,height = 14)#,width=4000,height=1500,res=400)
ggplot(action_All.ff,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=Number),col="lightgray")+
  #scale_fill_gradientn(colours=colorRampPalette(c("blue","white","red"),space="rgb")(100))+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",midpoint = 0,name="Cancer types")+
  # scale_fill_manual(limit=c(-1,0,1),values = c("lightblue","gray","lightpink"),na.value="gray",labels=c("Negatively-Cor","Non-Cor","Positively-Cor"),name="")+
  scale_x_discrete(limit=All_xlimit_label)+
  scale_y_discrete(limit=All_ylimit_label)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=18,colour = All_ylimit_col ),
        axis.text.x=element_text(size=18,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(0.8,"cm"),
        legend.position = "bottom",
        legend.direction ="horizontal",
        legend.key.width = unit(2,"cm"))+
  geom_tile(data = Ac_cirPPI_uniq_heatmap,aes(x=Actional,y=Circadian),alpha = 0,color="black",size=.4)
  
  #geom_tile(data = Ac_cirPPI_uniq_heatmap[which(Ac_cirPPI_uniq_heatmap$class=="Dir"),],aes(x=Actional,y=Circadian),alpha = 0,color="red",size=.8)+
 # geom_tile(data = Ac_cirPPI_uniq_heatmap[which(Ac_cirPPI_uniq_heatmap$class=="NonDir"),],aes(x=Actional,y=Circadian),alpha = 0,color="black",size=.8)

# annotate("text",y=All_ann$Circadian,x=All_ann$Actional,label=All_ann$Number,size=4)#+
# geom_polygon(data= data.frame(group=rep(1,times=4),poly_x=c()),aes(x=poly_x,y=poly_y,fill="lightblue", group=group)))
dev.off()
###Only show core circadian genes
colnames(Ac_cirPPI_uniq) <- c("Actional","Circadian","class")
Ac_cirPPI_uniq_heatmap <- merge(Ac_cirPPI_uniq,action_All.ff[which(action_All.ff$Overlap !=0),],by=c("Actional","Circadian"))

####
merge(action_All.ff[which(action_All.ff$Overlap %in% c(-1,1)),],actionable_ChIP.m[which(actionable_ChIP.m$value == "target"),],by.x=c("Actional","Circadian"),by.y=c("actionable","variable"))

####Only show 22 Circadian genes' pattern
FOXO <- c("FOXO1","FOXO3","FOXA2","FOXM1")
pos_action_All.m <- pos_action_All.m[which(!(pos_action_All.m$Circadian %in% FOXO)),]
pos_action_All.m$Circadian <- factor(pos_action_All.m$Circadian,levels = names(table(pos_action_All.m$Circadian))[1:22])
pos_ac <- sapply(split(pos_action_All.m[,3],pos_action_All.m$Actional),sum)
pos_action_All.m <- pos_action_All.m[which(pos_action_All.m$Actional %in% names(pos_ac)[pos_ac > 0]),]
neg_action_All.m <- neg_action_All.m[which(!(neg_action_All.m$Circadian %in% FOXO)),]
neg_action_All.m$Circadian <- factor(neg_action_All.m$Circadian,levels = names(table(neg_action_All.m$Circadian))[1:22])
neg_ac <- sapply(split(neg_action_All.m[,3],neg_action_All.m$Actional),sum)
neg_action_All.m <- neg_action_All.m[which(neg_action_All.m$Actional %in% names(neg_ac)[neg_ac > 0]),]

non_cor <- merge(pos_action_All.m[which(pos_action_All.m$Number ==0),],neg_action_All.m[which(neg_action_All.m$Number == 0),],by=c("Actional","Circadian"),all=T)[1:4]
colnames(non_cor) <- c("Actional","Circadian","Number","Overlap")
non_cor[is.na(non_cor)] <- 0
action_All <- rbind(pos_action_All.m[which(pos_action_All.m$Number >=1),],neg_action_All.m[which(neg_action_All.m$Number >=1),])
o <- merge(non_cor,action_All,by=c("Actional","Circadian"))[1:4]
fi <- paste(o$Actional,o$Circadian,sep="")
non_cor$fi <- paste(non_cor$Actional,non_cor$Circadian,sep="")
non_cor.f <- non_cor[!(non_cor$fi %in% fi),][1:4]
action_All.f <- rbind(pos_action_All.m[which(pos_action_All.m$Number >=1),],neg_action_All.m[which(neg_action_All.m$Number >=1),],non_cor.f)
All_xlimit <- sapply(split(action_All.f[,"Overlap"],action_All.f$Actional),sum)
#All_xlimit <- All_xlimit[All_xlimit !=0]
All_xlimit_label <- names(All_xlimit)[order(All_xlimit)]
write.table(All_xlimit_label,file = "all_cancer_type/Actional.genes.correlated.circadian.txt",quote = F,sep="\t",row.names = F)
All_ylimit <- sapply(split(action_All.f[,"Number"],action_All.f$Circadian),sum)
#All_ylimit <- All_ylimit[All_ylimit !=0]
All_ylimit_label <- names(All_ylimit)[order(All_ylimit)]
All_ann <- action_All.f[which(action_All.f$Number > 1),]

pdf("all_cancer_type/Only.circadian.correlated.genes_overlap.with.actional.genes1_5cancerTypes.pdf",width=18,height = 8)#,width=4000,height=1500,res=400)
ggplot(action_All.f,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=factor(Overlap)),col="white")+
  scale_fill_manual(limit=c(-1,0,1),values = c("lightblue","lightgray","lightpink"),na.value="gray",labels=c("Negatively-Cor","Non-Cor","Positively-Cor"),name="")+
  scale_x_discrete(limit=All_xlimit_label)+
  scale_y_discrete(limit=All_ylimit_label)+
  theme(panel.background=element_rect(colour="white",fill="white",size=2),
        panel.grid=element_blank(),
        #  panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=12,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        # legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(0.6,"cm"),
        legend.position = "bottom",
        legend.direction ="horizontal",
        legend.key.width = unit(0.6,"cm"))+
  annotate("text",y=All_ann$Circadian,x=All_ann$Actional,label=All_ann$Number,size=3)#+
# geom_polygon(data= data.frame(group=rep(1,times=4),poly_x=c()),aes(x=poly_x,y=poly_y,fill="lightblue", group=group)))
dev.off()

##
cpd.target.actionable <- read.delim("all_cancer_type/the.number.of.cpd.target.actionable.txt",header=T)
indirect.cpd.target.actionable<- read.delim("all_cancer_type/the.number.of.indirect.cpd.target.actionable.txt",header=T)
All_together <- rbind(cpd.target.actionable,indirect.cpd.target.actionable)
All_together$Class <- rep(c("Target drug","Correlated drug"),times=c(nrow(cpd.target.actionable),nrow(indirect.cpd.target.actionable)))
pdf("all_cancer_type/the.number.of.target.and.correlated.cpd.directly.target.actionable.pdf",width=18,height=4)
ggplot(All_together,aes(x=ActionableGene,y=DrugNum,fill=factor(Class)))+
  geom_bar(color="black",stat = "identity")+
  scale_y_continuous(expand=c(0,0),limit=c(-1,(max(sapply(split(indirect.cpd.target.actionable[,"DrugNum"],indirect.cpd.target.actionable$ActionableGene),sum))+30)))+
  scale_x_discrete(limits = All_xlimit_label)+
  scale_fill_manual(limit=c("Target drug","Correlated drug"),values = c("red","blue"),labels=c("Target drug","Correlated drug"),name="")+
  ylab("Compounds Count")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        panel.grid.minor=element_line(colour=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=12,color="black",angle=90,vjust=0.5,hjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank(),legend.text=element_text(size=10),
        legend.position=c(0.05,0.9),legend.background=element_blank())#,legend.direction="horizontal")#+
dev.off()
All_together$Class <- factor(All_together$Class,levels=c("Target drug","Correlated drug"))
pdf("all_cancer_type/the.number.of.target.and.correlated.cpd.directly.target.actionable_separate.pdf",width=18,height=5)
ggplot(All_together,aes(x=ActionableGene,y=DrugNum,fill=factor(Class)))+
  geom_bar(color="black",stat = "identity")+
  scale_y_continuous(expand=c(0,0.2))+#,limit=c(-1,(max(sapply(split(indirect.cpd.target.actionable[,"DrugNum"],indirect.cpd.target.actionable$ActionableGene),sum))+30)))+
  scale_x_discrete(limits = All_xlimit_label)+
  scale_fill_manual(limit=c("Target drug","Correlated drug"),values = c("red","blue"),labels=c("Target drug","Correlated drug"),guide=F)+
  facet_wrap(~Class,nrow=2,scale="free_y")+
  ylab("Compounds Count")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        panel.grid.minor=element_line(colour=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.y=element_text(size=12,color="black",angle=90,vjust=0.5,hjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),strip.text=element_text(size=16),
        legend.title=element_blank(),legend.text=element_text(size=10),
        legend.position=c(0.05,0.9),legend.background=element_blank())#,legend.direction="horizontal")#+
dev.off()



neg_action_All <- neg_action_All[,c("Gene",core.circadian)]
neg_action_All[2:15][neg_action_All[2:15] < 5] <- 0
neg_action_All <- neg_action_All[which(apply(neg_action_All[2:15],1,sum)>=5),]
neg_action_All$Gene <- as.character(neg_action_All$Gene)
###
neg_action_All.m <- melt(neg_action_All,id.vars = "Gene",measure.vars = colnames(neg_action_All)[2:ncol(neg_action_All)])

neg_action_All.m["Overlap"] <- rep(0,times=nrow(neg_action_All.m))
neg_action_All.m["Overlap"][neg_action_All.m["value"] >=1 ] <- -1
colnames(neg_action_All.m) <- c("Actional","Circadian","Number","Overlap")
neg.ann <- neg_action_All.m[which(neg_action_All.m$Number >=1),]
neg_xlimit <- neg_action_All[order(apply(neg_action_All[,2:ncol(neg_action_All)], 1,sum)),"Gene"]
ac_neg_y <- names(rev(sort(apply(neg_action_All[,2:ncol(neg_action_All)], 2,sum))))
ac_neg_ycol <- rep("black",times=length(ac_neg_y))
ac_neg_ycol[ac_neg_y %in% core.circadian] <- "red"
pdf("all_cancer_type/core_circadian.negatively.correlated.genes_overlap.with.actional.genes1_5cancerTypes.pdf",width=10,height = 6)#,width=4000,height=1500,res=400)
ggplot(neg_action_All.m,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=factor(Overlap)),col="white")+
  scale_fill_manual(limit=c(0,-1),values = c("gray","lightblue"),labels=c("Non-Cor","Negatively-Cor"),name="",guide=F)+
  scale_x_discrete(limit=neg_action_All[order(apply(neg_action_All[,2:ncol(neg_action_All)], 1,sum)),"Gene"])+
  scale_y_discrete(limit=ac_neg_y)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(1,"cm"),
        legend.key.width = unit(1,"cm"))+
  annotate("text",y=neg.ann$Circadian,x=neg.ann$Actional,label=neg.ann$Number,size=4)
dev.off()

pos_action_All <- read.delim("all_cancer_type/circadian.positively.correlated.genes_overlap.with.actional.genes.txt",header = T)
#set in more than 5 cancer types
pos_action_All <- pos_action_All[,c("Gene",core.circadian)]
pos_action_All[2:15][pos_action_All[2:15] < 5] <- 0

pos_action_All <- pos_action_All[which(apply(pos_action_All[2:15],1,sum)>=5),]
pos_action_All$Gene <- as.character(pos_action_All$Gene)
###
pos_action_All.m <- melt(pos_action_All,id.vars = "Gene",measure.vars = colnames(pos_action_All)[2:ncol(pos_action_All)])
pos_action_All.m["Overlap"] <- rep(0,times=nrow(pos_action_All.m))
pos_action_All.m["Overlap"][pos_action_All.m["value"] >=1 ] <- 1
colnames(pos_action_All.m) <- c("Actional","Circadian","Number","Overlap")
pos.ann <- pos_action_All.m[which(pos_action_All.m$Number >=1),]
pos_xlimit <- pos_action_All[order(apply(pos_action_All[,2:ncol(pos_action_All)], 1,sum)),"Gene"]
ac_pos_y <- names(rev(sort(apply(pos_action_All[,2:ncol(pos_action_All)], 2,sum))))
ac_pos_ycol <- rep("black",times=length(ac_pos_y))
ac_pos_ycol[ac_pos_y %in% core.circadian] <- "red"
pdf("all_cancer_type/core.circadian.positively.correlated.genes_overlap.with.actional.genes1_5cancerTypes.pdf",width=18,height = 6)#,width=4000,height=1500,res=400)
ggplot(pos_action_All.m,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=factor(Overlap)),col="white")+
  scale_fill_manual(limit=c(0,1),values = c("gray","lightpink"),labels=c("Non-Cor","Positively-Cor"),name="",guide=F)+
  scale_x_discrete(limit=pos_action_All[order(apply(pos_action_All[,2:ncol(pos_action_All)], 1,sum)),"Gene"])+
  scale_y_discrete(limit=ac_pos_y)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour ="black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(1,"cm"),
        legend.key.width = unit(1,"cm"))+
  annotate("text",y=pos.ann$Circadian,x=pos.ann$Actional,label=pos.ann$Number,size=4)
dev.off()
non_cor <- merge(pos_action_All.m[which(pos_action_All.m$Number ==0),],neg_action_All.m[which(neg_action_All.m$Number == 0),],by=c("Actional","Circadian"),all=T)[1:4]
colnames(non_cor) <- c("Actional","Circadian","Number","Overlap")
non_cor[is.na(non_cor)] <- 0
action_All <- rbind(pos_action_All.m[which(pos_action_All.m$Number >=1),],neg_action_All.m[which(neg_action_All.m$Number >=1),])
o <- merge(non_cor,action_All,by=c("Actional","Circadian"))[1:4]
fi <- paste(o$Actional,o$Circadian,sep="")
non_cor$fi <- paste(non_cor$Actional,non_cor$Circadian,sep="")
non_cor.f <- non_cor[!(non_cor$fi %in% fi),][1:4]
action_All.f <- rbind(pos_action_All.m[which(pos_action_All.m$Number >=1),],neg_action_All.m[which(neg_action_All.m$Number >=1),],non_cor.f)
All_xlimit <- sapply(split(action_All.f[,"Overlap"],action_All.f$Actional),sum)
#All_xlimit <- All_xlimit[All_xlimit !=0]
All_xlimit_label <- names(All_xlimit)[order(All_xlimit)]
All_ylimit <- sapply(split(action_All.f[,"Number"],action_All.f$Circadian),sum)
#All_ylimit <- All_ylimit[All_ylimit !=0]
All_ylimit_label <- names(All_ylimit)[order(All_ylimit)]
All_ann <- action_All.f[which(action_All.f$Number >1),]
All_ylimit_col <- rep("black",times=length(All_ylimit_label))
All_ylimit_col[All_ylimit_label %in% core.circadian] <- "red"
pdf("all_cancer_type/core_circadian.correlated.genes_overlap.with.actional.genes1_5cancerTypes.pdf",width=24,height = 6)#,width=4000,height=1500,res=400)
ggplot(action_All.f,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=factor(Overlap)),col="white")+
  scale_fill_manual(limit=c(-1,0,1),values = c("lightblue","gray","lightpink"),na.value="gray",labels=c("Negatively-Cor","Non-Cor","Positively-Cor"),name="")+
  scale_x_discrete(limit=All_xlimit_label)+
  scale_y_discrete(limit=All_ylimit_label)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_line(colour="grey",linetype="dashed"),
        panel.grid.major=element_line(colour="grey",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black" ),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=14),
        # legend.title=element_text(size=14),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(1,"cm"),
        legend.position = "bottom",
        legend.direction ="horizontal",
        legend.key.width = unit(1,"cm"))+
  annotate("text",y=All_ann$Circadian,x=All_ann$Actional,label=All_ann$Number,size=4)#+
# geom_polygon(data= data.frame(group=rep(1,times=4),poly_x=c()),aes(x=poly_x,y=poly_y,fill="lightblue", group=group)))
dev.off()
write.table(All_xlimit_label,file = "all_cancer_type/Actional.genes.correlated.core_circadian.txt",quote = F,sep="\t",row.names = F)
##directly target actionable genes
cpd.target.actionable <- read.delim("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/the.number.of.cpd.target.actionable.txt",header = T)
indirect.cpd.target.actionable <- read.delim("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/51_the.number.of.indirect.cpd.target.actionable.txt",header=T)

pdf("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/51_the.number.of.cpd.directly.target.actionable.pdf",width=24,height=2.5)
ggplot(cpd.target.actionable,aes(x=ActionableGene,y=DrugNum))+
  geom_bar(color=NA,fill="blue",stat = "identity")+
  scale_y_continuous(expand=c(0,0),limit=c(-0.05,(max(cpd.target.actionable$DrugNum)+2)))+
  scale_x_discrete(limits = All_xlimit_label)+
  ylab("Count")+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.1),
        panel.grid.minor=element_line(colour=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18,color="black",angle=90,vjust=0.5,hjust=0.5),
        legend.title=element_blank(),legend.text=element_text(size=14),
        #axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(color="black"),
       # axis.ticks=element_line(color="black"),
        legend.position="bottom",legend.direction="horizontal")#+
dev.off()

##indirectly target actional genes
pdf("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/indirect.cpd.target.actionable51_the.number.of.indirect.cpd.directly.target.actionable1.pdf",width=24,height=2.5)
ggplot(indirect.cpd.target.actionable,aes(x=ActionableGene,y=DrugNum))+
  geom_bar(color=NA,fill="red",stat = "identity")+
  scale_y_continuous(expand=c(0,0),limit=c(-0.5,(max(indirect.cpd.target.actionable$DrugNum)+20)))+
  scale_x_discrete(limits = All_xlimit_label)+
  ylab("Count")+
  theme(panel.background=element_rect(colour="black",fill="white"),
        panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.1),
        panel.grid.minor=element_line(colour=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=18,color="black",angle=90,vjust=0.5,hjust=0.5),
        legend.title=element_blank(),legend.text=element_text(size=14),
        #axis.text.x=element_text(size=16,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(color="black"),
        # axis.ticks=element_line(color="black"),
        legend.position="bottom",legend.direction="horizontal")#+
dev.off()


for( i in cir_names){
  for(j in 1:length(corfiles)){
    cor_data <- read.delim(paste("/home/yye1/Circadian/expression/correlation/pearson/",corfiles[j],sep = ""),header=T)
    cor_data <- cor_data[colnames(cor_data)[c(1,grep(i,colnames(cor_data)))]]
    posGene <- as.character(cor_data[which(cor_data[,2] >= 0.3 & cor_data[,3] < 0.05),"gene"])
    posGene <- intersect(posGene,actional.genes$Gene)
    negGene <- as.character(cor_data[which(cor_data[,2] <= -0.3 & cor_data[,3] < 0.05),"gene"])
    negGene <- intersect(negGene,actional.genes$Gene)
    pos_Ac_cir <- data.frame(Actionable = posGene,Type = rep(j,length(posGene)),Cir = rep(i,length(posGene)))
    neg_Ac_cir <- data.frame(Actionable = negGene,Type = rep(j,length(negGene)),Cir = rep(i,length(negGene)))
    if( j ==1){
      pos_Ac_cirA <- pos_Ac_cir
      neg_Ac_cirA <- neg_Ac_cir
    }else{
      pos_Ac_cirA <- rbind(pos_Ac_cirA,pos_Ac_cir)
      neg_Ac_cirA <- rbind(neg_Ac_cirA,neg_Ac_cir)
    }
  }
  if(i == cir_names[1]){
    pos_Ac_cirAll <- pos_Ac_cirA
    neg_Ac_cirAll <- neg_Ac_cirA
  }else{
    pos_Ac_cirAll <- rbind(pos_Ac_cirAll,pos_Ac_cirA)
    neg_Ac_cirAll <- rbind(neg_Ac_cirAll,neg_Ac_cirA)
  }
}

corfilesAnn <- data.frame(CancerType=gsub("_gene_cor.pearson.txt","",corfiles),Type=seq(1,32))
Ac_cirAll <- rbind(pos_Ac_cirAll,neg_Ac_cirAll)
Ac_cirAll$Class <- rep(c("pos","neg"),times=c(nrow(pos_Ac_cirAll),nrow(neg_Ac_cirAll)))
Ac_cirAll <- merge(Ac_cirAll,corfilesAnn,by="Type")

All_AcCir_pair <- data.frame(sapply(split(Ac_cirAll[,"Class"],Ac_cirAll$CancerType),table))
All_AcCir_pair$Class <- c("neg","pos")
All_AcCir_pair.m <- melt(All_AcCir_pair,id.var="Class",measure.vars=colnames(All_AcCir_pair)[1:(ncol(All_AcCir_pair)-1)])
pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/cir_actionable_correlation pair in each cancer type.pdf",width=8,height = 2.8)#,width=4000,height=1500,res=400)
ggplot(Ac_cirAll,aes(x=CancerType,fill=factor(Class)))+
  geom_bar(color=NA,width = 0.6)+
  scale_fill_manual(limit=c("neg","pos"),values = c("blue","red"),labels=c("Negative","Positive"),name="")+
  ylab("Correlation pairs")+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(expand = c(0.02,0),breaks = seq(0,2500,length.out = 6))+
  scale_x_discrete(limit= names(table(Ac_cirAll$CancerType))[order(table(Ac_cirAll$CancerType))],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_blank(),#line(linetype="dashed",color="lightgray"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=14,color="black"),
        axis.text.y=element_text(size=12,colour = "black"),
        axis.text.x=element_text(size=12,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
       legend.position = c(0.15,0.86),legend.direction ="horizontal",
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_hline(yintercept = mean(All_actionabel_pairs),linetype="dashed")
dev.off()
sapply(split(Ac_cirAll[,"Actionable"],list(Ac_cirAll$CancerType,Ac_cirAll$Class)),function(x){length(unique(x))})

All_actionabel_pairs <- sapply(split(All_AcCir_pair.m[,"value"],list(All_AcCir_pair.m$variable)),sum)



All_actionabel_core_pairs <- sapply(split(Ac_corecirAll[,"value"],Ac_corecirAll$variable),sum)

Ac_corecirAll <- Ac_cirAll[which(Ac_cirAll$Cir %in% core.circadian),]
pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/core cir_actionable_correlation pair in each cancer type.pdf",width=8,height = 2.8)#,width=4000,height=1500,res=400)
ggplot(Ac_corecirAll,aes(x=CancerType,fill=factor(Class)))+
  geom_bar(color=NA,width = 0.5)+
  scale_fill_manual(limit=c("neg","pos"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(expand = c(0.02,0),breaks = seq(0,600,length.out = 4))+
  scale_x_discrete(limit= names(table(Ac_corecirAll$CancerType))[order(table(Ac_corecirAll$CancerType))],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_blank(),#line(linetype="dashed",color="lightgray"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))
dev.off()
actionable_pos_neg_num_eachCancer <- sapply(split(Ac_cirAll[,"Actionable"],list(Ac_cirAll$CancerType,Ac_cirAll$Class)),function(x){length(unique(x))})
actionable_pos_neg_num_eachCancer <- data.frame(actionable_pos_neg_num_eachCancer)
colnames(actionable_pos_neg_num_eachCancer) <- "value"
actionable_pos_neg_num_eachCancer$Type <-  data.frame(do.call(rbind,strsplit(as.character(rownames(actionable_pos_neg_num_eachCancer)),"\\.")))$X1
actionable_pos_neg_num_eachCancer$Class <-  data.frame(do.call(rbind,strsplit(as.character(rownames(actionable_pos_neg_num_eachCancer)),"\\.")))$X2
function(x){
  a <- as.numeric(data.frame(do.call(rbind,strsplit(as.character(x),"e")))$X1)
  b <- as.numeric(data.frame(do.call(rbind,strsplit(as.character(x),"e")))$X2)
  text(4,4,substitute(x*10^y, list(x=a,y=b)))
}


pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/actionable_pos_neg in each cancer type.pdf",width=8,height = 2.8)#,width=4000,height=1500,res=400)
ggplot(actionable_pos_neg_num_eachCancer,aes(x=Type,y=value,fill=factor(Class)))+
  geom_bar(color=NA,stat = "identity",width = 0.5)+
  scale_fill_manual(limit=c("neg","pos"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(expand = c(0.02,0),breaks = seq(0,200,length.out = 5))+
  ylab("Correlation pairs")+
  scale_x_discrete(limit=names(sapply(split(actionable_pos_neg_num_eachCancer[,"value"],actionable_pos_neg_num_eachCancer$Type),sum))[order(sapply(split(actionable_pos_neg_num_eachCancer[,"value"],actionable_pos_neg_num_eachCancer$Type),sum))],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_blank(),#line(linetype="dashed",color="lightgray"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=15,color="black"),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_hline(yintercept =mean(sapply(split(actionable_pos_neg_num_eachCancer[,"value"],actionable_pos_neg_num_eachCancer$Type),sum)),linetype="dashed")
dev.off()

actionable_num_eachCancer <- sapply(split(Ac_cirAll[,"Actionable"],list(Ac_cirAll$CancerType)),function(x){length(unique(x))})


core.actionable_num_eachCancer <- sapply(split(Ac_corecirAll[,"Actionable"],list(Ac_corecirAll$CancerType)),function(x){length(unique(x))})
core.actionable_pos_neg_num_eachCancer <- sapply(split(Ac_corecirAll[,"Actionable"],list(Ac_corecirAll$CancerType,Ac_corecirAll$Class)),function(x){length(unique(x))})
core.actionable_pos_neg_num_eachCancer <- data.frame(core.actionable_pos_neg_num_eachCancer)
colnames(core.actionable_pos_neg_num_eachCancer) <- "value"
core.actionable_pos_neg_num_eachCancer$Type <-  data.frame(do.call(rbind,strsplit(as.character(rownames(core.actionable_pos_neg_num_eachCancer)),"\\.")))$X1
core.actionable_pos_neg_num_eachCancer$Class <-  data.frame(do.call(rbind,strsplit(as.character(rownames(core.actionable_pos_neg_num_eachCancer)),"\\.")))$X2

pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/core circadian actionable_pos_neg in each cancer type.pdf",width=8,height = 2.8)#,width=4000,height=1500,res=400)
ggplot(core.actionable_pos_neg_num_eachCancer,aes(x=Type,y=value,fill=factor(Class)))+
  geom_bar(color=NA,stat = "identity",width = 0.5)+
  scale_fill_manual(limit=c("neg","pos"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(expand = c(0.02,0),breaks = seq(0,200,length.out = 5))+
  scale_x_discrete(limit=names(sapply(split(core.actionable_pos_neg_num_eachCancer[,"value"],core.actionable_pos_neg_num_eachCancer$Type),sum))[order(sapply(split(core.actionable_pos_neg_num_eachCancer[,"value"],core.actionable_pos_neg_num_eachCancer$Type),sum))],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_blank(),#line(linetype="dashed",color="lightgray"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=14,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_hline(yintercept =mean(sapply(split(core.actionable_pos_neg_num_eachCancer[,"value"],core.actionable_pos_neg_num_eachCancer$Type),sum)),linetype="dashed")
dev.off()
core.actionable_num_eachCancer <- data.frame(core.actionable_num_eachCancer)
colnames(core.actionable_num_eachCancer) <- "value"
core.actionable_num_eachCancer$Type <- rownames(core.actionable_num_eachCancer)
pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/core circadian actionable_unique in each cancer type.pdf",width=8,height = 2.8)#,width=4000,height=1500,res=400)
ggplot(core.actionable_num_eachCancer,aes(x=Type,y=value))+
  geom_bar(color=NA,stat = "identity",width = 0.5,fill="darkgreen")+
  #scale_fill_manual(limit=c("neg","pos"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_continuous(expand = c(0.02,0),breaks = seq(0,200,length.out = 5))+
  ylab("Number of Actionable genes")+
  scale_x_discrete(limit=core.actionable_num_eachCancer[order(core.actionable_num_eachCancer$value),"Type"],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_blank(),#line(linetype="dashed",color="lightgray"),
        axis.title.y=element_text(size=14,colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size=14,colour = "black"),
        axis.text.x=element_text(size=12,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_hline(yintercept =mean(core.actionable_num_eachCancer$value),linetype="dashed")
dev.off()

actionable_num_eachCancer <- data.frame(actionable_num_eachCancer)
colnames(actionable_num_eachCancer) <- "value"
actionable_num_eachCancer$Type <- rownames(actionable_num_eachCancer)
pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/circadian actionable_unique in each cancer type.pdf",width=5,height = 2.8)#,width=4000,height=1500,res=400)
ggplot(actionable_num_eachCancer,aes(x=Type,y=value))+
  geom_bar(color=NA,stat = "identity",width = 0.5,fill="darkgreen")+
  #scale_fill_manual(limit=c("neg","pos"),values = c("blue","red"),guide=FALSE)+
  # scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  ylab("Actionable genes")+
  scale_y_continuous(expand = c(0.02,0),breaks = seq(0,200,length.out = 5))+
  scale_x_discrete(limit=actionable_num_eachCancer[order(actionable_num_eachCancer$value),"Type"],expand=c(0.01,0.01))+
  theme(panel.background=element_rect(colour="black",fill="white",size=1),
        panel.grid.major=element_blank(),#line(linetype="dashed",color="lightgray"),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size = 12,color = "black"),
        axis.text.y=element_text(size=10,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.key = element_rect(fill="white",colour = "black"))+
  geom_hline(yintercept =mean(actionable_num_eachCancer$value),linetype="dashed")
dev.off()

###
PPI <- read.table("/extraspace/liucj/reference/PPI/all_PPI_data.tsv",sep="\t",header=T)
Ac_to_cir <- PPI[which((PPI$P1_symbol %in% actionable.genes$x) & (PPI$P2_symbol %in% circadian.genes)),]
cir_to_Ac <- PPI[which((PPI$P2_symbol %in% actionable.genes$x) & (PPI$P1_symbol %in% circadian.genes)),]

Ac_cirPPIAll <- unique(data.frame(Actionable = c(Ac_to_cir$P1_symbol,cir_to_Ac$P2_symbol),Circadian = c(Ac_to_cir$P2_symbol,cir_to_Ac$P1_symbol),evidence=c(Ac_to_cir$evidence,cir_to_Ac$evidence),pubmed=c(Ac_to_cir$pubmed,cir_to_Ac$pubmed),source=c(Ac_to_cir$source,cir_to_Ac$source)))
Ac_cirPPI <- unique(Ac_cirPPIAll[,1:3])
write.csv(Ac_cirPPI,file="/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/Ac_cirPPI.csv",quote=F)
Ac_cirPPI$class <- rep("Dir",nrow(Ac_cirPPI))
Ac_cirPPI[,"class"][Ac_cirPPI[,"evidence"] %in% c("Dosage Lethality","Synthetic Lethality","Synthetic Growth Defect","in vitro","in vitro;in vivo","in vivo")] <- "NonDir"
Ac_cirPPI <- unique(Ac_cirPPI[,c(1,2,4)])
Ac_cirPPI_uniq <- data.frame()
for( i in 1:nrow(unique(Ac_cirPPI[,1:2]))){
  x <- unique(Ac_cirPPI[,1:2])[i,]
  sub <- Ac_cirPPI[which(Ac_cirPPI$Actionable == x[,1] & Ac_cirPPI$Circadian == x[,2]),]
  if(is.na(match("Dir",sub$class)) == FALSE){
    sub$class <- rep("Dir",nrow(sub))
  }
  if(nrow(Ac_cirPPI_uniq)==0){
    Ac_cirPPI_uniq <- sub
  }else(
    Ac_cirPPI_uniq <- rbind(Ac_cirPPI_uniq,sub)
  )
}
Ac_cirPPI_uniq <- unique(Ac_cirPPI_uniq)

Ac_corecirPPI <- Ac_cirPPI_uniq[which(Ac_cirPPI_uniq$Circadian %in% core.circadian),]

###
df <- graph.data.frame(Ac_corecirPPI,directed=T)
df1 <- graph.data.frame(data.frame(First = unique(Ac_corecirPPI[,1]), second = unique(Ac_corecirPPI[,1])))
tt1 <- layout.circle(df1)
df2 <- graph.data.frame(data.frame(First = unique(Ac_corecirPPI[,2]), second = unique(Ac_corecirPPI[,2])))
tt2 <- layout.circle(df2)
tt <- layout.circle(df)
#########  rescale two layout
n1=length(unique(Ac_corecirPPI[,1]))
n2=length(unique(Ac_corecirPPI[,2]))
tt[1:n1,1] <- tt1[,1]*1.7
tt[1:n1,2] <- tt1[,2]*1.7
tt[(n1+1):(n1+n2),1] <- tt2[,1]
tt[(n1+1):(n1+n2),2] <- tt2[,2]
#sig.cir <- as.vector(unique(BRCA.nodes$Gene)[unique(BRCA.nodes$Gene) %in% core.circadian])
ttname <- V(df)$name
V(df)$color <- "tomato"
V(df)[1:n1]$color <- "yellowgreen"
#V(df)[ttname %in% sig.cir]$color <- "tomato"
E(df)$arrow.mode = 0
pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/Interaction of actionable and core circadian genes.pdf",width = 6,height = 6)
plot(df,
     layout=tt,vertex.label.color="black",vertex.size=15,vertex.label.cex=0.8)
dev.off()

Ac_cirPPI_uniq <- Ac_cirPPI_uniq[which(Ac_cirPPI_uniq$Circadian != "CREBBP"),]
df <- graph.data.frame(Ac_cirPPI_uniq,directed=T)
df1 <- graph.data.frame(data.frame(First = unique(Ac_cirPPI_uniq[,1]), second = unique(Ac_cirPPI_uniq[,1])))
tt1 <- layout.circle(df1)
df2 <- graph.data.frame(data.frame(First = unique(Ac_cirPPI_uniq[,2]), second = unique(Ac_cirPPI_uniq[,2])))
tt2 <- layout.circle(df2)
tt <- layout.circle(df)
#########  rescale two layout
n1=length(unique(Ac_cirPPI_uniq[,1]))
n2=length(unique(Ac_cirPPI_uniq[,2]))
tt[1:n1,1] <- tt1[,1]*1.7
tt[1:n1,2] <- tt1[,2]*1.7
tt[(n1+1):(n1+n2),1] <- tt2[,1]
tt[(n1+1):(n1+n2),2] <- tt2[,2]
#sig.cir <- as.vector(unique(BRCA.nodes$Gene)[unique(BRCA.nodes$Gene) %in% core.circadian])
ttname <- V(df)$name
V(df)$color <- "magenta"
V(df)[1:n1]$color <- "green2"
#V(df)[ttname %in% sig.cir]$color <- "magenta2"
E(df)$arrow.mode = 0
pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/Interaction of actionable and circadian genes.pdf",width = 10,height = 10)
plot(df,#edge.color=c("gold","black")[(E(df)$class=="NonDir")+1],
     layout=tt,vertex.label.color="black",vertex.size=9,vertex.label.cex=0.8)
dev.off()



CRY2_sp <- setdiff(action_All.ff[which(action_All.ff$Circadian == "CRY2" & abs(action_All.ff$Number) > 0),"Actional"],action_All.ff[which(action_All.ff$Circadian == "CRY1" & abs(action_All.ff$Number) > 0),"Actional"])
CRY1_sp <- setdiff(action_All.ff[which(action_All.ff$Circadian == "CRY1" & abs(action_All.ff$Number) > 0),"Actional"],action_All.ff[which(action_All.ff$Circadian == "CRY2" & abs(action_All.ff$Number) > 0),"Actional"])
CRY1_2_overlap <- 
write.table(CRY2_sp,file = "/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/CRY2_sp_interact_actionable.genes.txt",quote = F,sep="\t",row.names = F,col.names = F)
write.table(CRY1_sp,file = "/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/CRY1_sp_interact_actionable.genes.txt",quote = F,sep="\t",row.names = F,col.names = F)

write.table(action_All[which(action_All$Circadian == "CRY2" & action_All$Overlap==1),"Actional"],file ="/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/CRY2_positive_actionableGenes.txt",quote = F,row.names = F,col.names = F )
write.table(action_All[which(action_All$Circadian == "CRY2" & action_All$Overlap==-1),"Actional"],file ="/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/CRY2_negative_actionableGenes.txt",quote = F,row.names = F,col.names = F )



#######circadian genes ChIP-seq regulated actionable genes
setwd("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/circadian.genes.ChIPseq")
actional.genes <- read.delim("~/Circadian/expression/correlation/pearson/CRY1_CRY2_cor0.3/actionable.genes.txt",header=T)
GSE70686 <- read.delim("GSE70686_ClockBmal.txt",header=F)
library(VennDiagram)
group<-list(GSE70686$V1,actional.genes$Gene)
 names(group)<-c("CLOCK/ARNTL","Actionable genes")
# pdf("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/circadian.genes.ChIPseq/GSE70686 and actionable genes overlap.pdf",width=5,height=5)
par(mar=c(2,4,1,10),xpd=T)
ven.plot<-venn.diagram(group, fill=c("red3","blue3","green3"),alpha=c(0.8,0.8),cat.cex=2.2, cex=2,filename=NULL)
grid.draw(ven.plot)
dev.off()

homology <- read.csv("mouse.human.homology.genes.csv",header=T)
actionable_ChIP <- data.frame(actionable = actional.genes$Gene)
Sheetnames <- c("ARNTL","CLOCK","NPAS2","PER1","PER2","CRY1","CRY2")
for(i in 1:7){
  sub <- readxl::read_xlsx("Koike_TableS2_MasterPeakList_072612.xlsx",sheet= i)
  homoloysub <- unique(homology[which(homology$mgi_symbol %in% unique(sub$Symbol)),"hgnc_symbol"])
  homoloysub <- data.frame(homoloysub)
  homoloysub$ChIP <- rep(1,times=nrow(homoloysub))
  colnames(homoloysub) <- c("actionable",Sheetnames[i])
  actionable_ChIP <- merge(actionable_ChIP,homoloysub,by="actionable",all.x=T)}
actionable_ChIP[is.na(actionable_ChIP)] <- 0
actionable_ChIP <- actionable_ChIP[which(apply(actionable_ChIP[,2:8],1,sum)>0),]
actionable_ChIP.m <- melt(actionable_ChIP,id.vars = "actionable",measure.vars = Sheetnames)
cir_order <- sapply(split(actionable_ChIP.m[,"value"],actionable_ChIP.m$variable),sum)
actionable_order <- sapply(split(actionable_ChIP.m[,"value"],actionable_ChIP.m$actionable),sum)
actionable_ChIP.m["value"][actionable_ChIP.m["value"]==1] <- "target"
actionable_ChIP.m["value"][actionable_ChIP.m["value"]==0] <- "none"

pdf("~/Circadian/expression/correlation/pearson/Circadian_Actional_genes/circadian.genes.ChIPseq/circadian targeted actionable genes ChIP.pdf",width=10,height = 4)#,width=4000,height=1500,res=400)
ggplot(actionable_ChIP.m,aes(x=actionable,y=variable,fill=factor(value)))+
  geom_tile(color="white",show.legend = F)+
   scale_fill_manual(limit=c("target","none"),values=c("red","gray"))+
# scale_size_continuous(limit=c(-log10(0.05),15),range = c(1, 10),breaks=c(-log10(0.05),5,10,15),labels=c("0.05","1e-5","1e-10","<1e-15"))+
  scale_y_discrete(limit= names(cir_order[order(cir_order )]))+
  scale_x_discrete(limit= names(actionable_order[order(actionable_order)]))+
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size=10,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"))
dev.off()

NG <- read.delim("NG_ADARB1_list.txt",header = T)
NGhomo <- homology[which(homology$mgi_symbol %in% NG$gene),]

Science <- actionable_ChIP.m[which(actionable_ChIP.m$value == "target"),1:2]

BMAL1_GSE69100 <- read.delim("GSE69100_BMAL1_U2OS.cell",header=F)
BMAL1_GSE69100 <- BMAL1_GSE69100[which(BMAL1_GSE69100$V1 %in% actional.genes$Gene),]
BMAL1_GSE69100 <- data.frame(actionable=BMAL1_GSE69100,variable=rep("ARNTL",length(BMAL1_GSE69100)))
GSE34019_NR1D1 <- read.delim("GSM840528_NR1D1.uniq.Gene",header=F)
GSE34019_NR1D1 <- unique(homology[which(homology$mouse_ensembl_id %in% GSE34019_NR1D1$V1),"hgnc_symbol"])
GSE34019_NR1D1 <- GSE34019_NR1D1[GSE34019_NR1D1 %in% actional.genes$Gene]
GSE34019_NR1D1 <- data.frame(actionable=GSE34019_NR1D1,variable=rep("NR1D1",length(GSE34019_NR1D1)))
GSE34019_NR1D2 <- read.delim("GSM840529_NR1D2.uniq.Gene",header=F)
GSE34019_NR1D2 <- unique(homology[which(homology$mouse_ensembl_id %in% GSE34019_NR1D2$V1),"hgnc_symbol"])
GSE34019_NR1D2 <- GSE34019_NR1D2[GSE34019_NR1D2 %in% actional.genes$Gene]
GSE34019_NR1D2 <- data.frame(actionable=GSE34019_NR1D2,variable=rep("NR1D2",length(GSE34019_NR1D2)))
GSE67962_RORA <- read.delim("GSE67962_RORA.tsv",header=T)
GSE67962_RORA <- unique(homology[which(homology$mgi_symbol %in% GSE67962_RORA$Target_genes),"hgnc_symbol"])
GSE67962_RORA <- GSE67962_RORA[GSE67962_RORA %in% actional.genes$Gene]
GSE67962_RORA <- data.frame(actionable=GSE67962_RORA,variable=rep("NR1D2",length(GSE67962_RORA)))
GSE67962_RORC <- read.delim("GSE67962_RORC.tsv",header=T)
GSE67962_RORC <- unique(homology[which(homology$mgi_symbol %in% GSE67962_RORC$Target_genes),"hgnc_symbol"])
GSE67962_RORC <- GSE67962_RORC[GSE67962_RORC %in% actional.genes$Gene]
GSE67962_RORC <- data.frame(actionable=GSE67962_RORC,variable=rep("NR1D2",length(GSE67962_RORC)))


CoreCircadianChIP <- rbind(Science,BMAL1_GSE69100,GSE34019_NR1D1,GSE34019_NR1D2,GSE67962_RORA,GSE67962_RORC)
CoreCircadianChIP <- unique(CoreCircadianChIP )
CoreCircadianChIPCor <- merge(action_All.ff[which(action_All.ff$Overlap %in% c(-1,1)),],CoreCircadianChIP,by.x=c("Actional","Circadian"),by.y=c("actionable","variable"))

pdf("/home/yye1/Circadian/expression/correlation/pearson/Circadian_Actional_genes/all_cancer_type/51circadian.correlated.genes_overlap.with.actional.genes1_5cancerTypes_colorful_withoutann_sigPPI_ChIP.pdf",width=28,height = 14)#,width=4000,height=1500,res=400)
ggplot(action_All.ff,aes(y=Circadian,x=Actional))+
  geom_tile(aes(fill=Number),col="lightgray")+
  #scale_fill_gradientn(colours=colorRampPalette(c("blue","white","red"),space="rgb")(100))+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",midpoint = 0,name="Cancer types")+
  # scale_fill_manual(limit=c(-1,0,1),values = c("lightblue","gray","lightpink"),na.value="gray",labels=c("Negatively-Cor","Non-Cor","Positively-Cor"),name="")+
  scale_x_discrete(limit=All_xlimit_label)+
  scale_y_discrete(limit=All_ylimit_label)+
  theme(panel.background=element_rect(colour="black",fill="white",size=2),
        panel.grid=element_blank(),
        panel.grid.major=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size=18,colour = All_ylimit_col ),
        axis.text.x=element_text(size=18,colour = "black",angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.key = element_rect(colour = "black"),
        legend.key.height = unit(0.8,"cm"),
        legend.position = "bottom",
        legend.direction ="horizontal",
        legend.key.width = unit(2,"cm"))+
  geom_tile(data = Ac_cirPPI_uniq_heatmap,aes(x=Actional,y=Circadian),alpha = 0,color="black",size=.4)+
  geom_point(data = CoreCircadianChIPCor,aes(x=Actional,y=Circadian,size=1.5),shape=4,color="black")
dev.off()
###Only show core circadian genes
colnames(Ac_cirPPI_uniq) <- c("Actional","Circadian","class")
Ac_cirPPI_uniq_heatmap <- merge(Ac_cirPPI_uniq,action_All.ff[which(action_All.ff$Overlap !=0),],by=c("Actional","Circadian"))

####
merge(action_All.ff[which(action_All.ff$Overlap %in% c(-1,1)),],actionable_ChIP.m[which(actionable_ChIP.m$value == "target"),],by.x=c("Actional","Circadian"),by.y=c("actionable","variable"))


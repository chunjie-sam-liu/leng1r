###Merge surival, subtype, stage data 
setwd("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype")
stage <- read.csv("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.stage/significant.DE.clock.genes.in.different.stage.csv",header=T)
subtype <- read.csv("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/significant.DE.clock.genes.in.different.subtype.csv",header = T)
survival <- read.delim("/home/yye1/Circadian/survival/circadian.genes.survival.count.txt",header=T)
survival <- survival[which(!is.na(survival$survival)==T),]
stage$stage <- rep("Stage",times=nrow(stage))
stage <- stage[which(abs(stage$FC_log2) >= log2(1.5)),]
subtype$subtype <- rep("Subtype",times=nrow(subtype))
subtype <- subtype[which(abs(subtype$FC_log2) >= log2(1.5)),]
union_cancer <- union(union(stage$tumor,subtype$tumor),survival$tumor)
#union_cancer <- union_cancer[!(union_cancer %in% c("DLBC","READ","CESC"))]
union_data <- data.frame(geneSymbol=rep(circadian.genes,each=length(union_cancer)),tumor=rep(union_cancer,times=51))
union_data <- merge(merge(merge(union_data,survival,by=c("tumor","geneSymbol"),all.x=T),stage[,c("tumor","geneSymbol","stage")],by=c("tumor","geneSymbol"),all.x=T),subtype[,c("tumor","geneSymbol","subtype")],by=c("tumor","geneSymbol"),all.x=T)
union_data.s <- union_data
three_overlap <- union_data.s[which(!is.na(union_data.s$survival)==T & !is.na(union_data.s$stage)==T & !is.na(union_data.s$subtype)==T),]

survival_stage <- union_data.s[which(!is.na(union_data.s$survival)==T & !is.na(union_data.s$stage)==T & !is.na(union_data.s$subtype)==F),]
survival_subtype <- union_data.s[which(!is.na(union_data.s$survival)==T & !is.na(union_data.s$stage)==F & !is.na(union_data.s$subtype)==T),]
stage_subtype <- union_data.s[which(!is.na(union_data.s$survival)==F & !is.na(union_data.s$stage)==T & !is.na(union_data.s$subtype)==T),]



union_data.m <- melt(union_data,id.vars=c("tumor","geneSymbol"),measure.vars=c("survival","stage","subtype"))
union_data.m <- union_data.m[which(!is.na(union_data.m$value==T)),]
unique_data.m <- union_data.m[!duplicated(union_data.m[,1:2]),]
unique_data.m$geneSymbol <- factor(unique_data.m$geneSymbol,levels=unique(unique_data.m$geneSymbol))
dup <- union_data.m[duplicated(union_data.m[,1:2]),]
unique_dup <- dup[!duplicated(dup[,1:2]),]
unique_dup_1 <- unique_dup[!(unique_dup$tumor == "KIRC" & unique_dup$geneSymbol %in% c("HLF","PRKAA2","TEF","PER3")),]
dup_dup_stage <- unique_dup[which(unique_dup$tumor == "KIRC" & unique_dup$geneSymbol %in% c("HLF","PRKAA2","TEF","PER3")),]
dup_dup_subtype <- dup[duplicated(dup[,1:2]),]

sss_x_label <- names(table(union_data.m$tumor))#[order(table(union_data.m$tumor))] #union_cancer[order(union_cancer)] #
sss_y_label <- table(union_data.m$geneSymbol)[table(union_data.m$geneSymbol)>0]
sss_y_label <- names(sss_y_label)[order(sss_y_label)]
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
sss_y_color <- rep("black",times=length(sss_y_label))
sss_y_color[sss_y_label %in% core.circadian] <- "red"


unique_data.mm <- merge(data.frame(geneSymbol=rep(sss_y_label,each=length(union_cancer)),tumor=rep(union_cancer,times=length(sss_y_label))),unique_data.m,by=c("tumor","geneSymbol"),all.x=T)


unique_dup$x = sapply(unique_dup$tumor,function(x){grep(x,sss_x_label)}) ## add the geom_polygon ID, convert the factor to numeric 
unique_dup$y = sapply(unique_dup$geneSymbol,function(x){grep(x,sss_y_label)})#as.numeric(dup.data$geneSymbol)

###stage duplicate
unique_dup_stage <- unique_dup[which(unique_dup$variable=="stage"),]
poly_x = c()
poly_y = c()
for ( i in 1: nrow(unique_dup_stage)){
  poly_x = c(poly_x, unique_dup_stage[i,"x"]-0.498, unique_dup_stage[i,'x'],unique_dup_stage[i,"x"],unique_dup_stage[i,"x"]-0.498 )
  poly_y = c(poly_y, unique_dup_stage[i,"y"]-0.498,unique_dup_stage[i,"y"]-0.498,unique_dup_stage[i,"y"]+0.498,unique_dup_stage[i,"y"]+0.498)
}
polygonID_stage <- data.frame( group = rep(seq(1:nrow(unique_dup_stage)), each = 4),poly_x, poly_y )

###subtype duplicate
unique_dup_subtype <- unique_dup[which(unique_dup$variable=="subtype"),]

poly_x = c()
poly_y = c()
for ( i in 1: nrow(unique_dup_subtype)){
  poly_x = c(poly_x, unique_dup_subtype[i,"x"]-0.498, unique_dup_subtype[i,'x'],unique_dup_subtype[i,"x"],unique_dup_subtype[i,"x"]-0.498 )
  poly_y = c(poly_y, unique_dup_subtype[i,"y"]-0.49,unique_dup_subtype[i,"y"]-0.49,unique_dup_subtype[i,"y"]+0.49,unique_dup_subtype[i,"y"]+0.49)
}
polygonID_subtype <- data.frame( group = rep(seq(1:nrow(unique_dup_subtype)), each = 4),poly_x, poly_y )

###three ovelap ---stage
dup_dup_stage$x <- match(dup_dup_stage$tumor,sss_x_label)
dup_dup_stage$y <- match(dup_dup_stage$geneSymbol,sss_y_label)

poly_x = c(dup_dup_stage[1,"x"]-0.498+1/3, dup_dup_stage[1,'x']-0.498+2/3,dup_dup_stage[1,"x"]-0.498+2/3,dup_dup_stage[1,"x"]-0.498+1/3 )
poly_y = c(dup_dup_stage[1,"y"]-0.498,dup_dup_stage[1,"y"]-0.498,dup_dup_stage[1,"y"]+0.498,dup_dup_stage[1,"y"]+0.498)
polygonID_threeOverlap_stage <- data.frame(group=rep(1,times=4),poly_x,poly_y)

###three ovelap ---subtype
dup_dup_subtype$x <- match(dup_dup_subtype$tumor,sss_x_label)
dup_dup_subtype$y <- match(dup_dup_subtype$geneSymbol,sss_y_label)

poly_x = c(dup_dup_subtype[1,"x"]-0.5, dup_dup_subtype[1,'x']-0.498+1/3,dup_dup_subtype[1,"x"]-0.498+1/3,dup_dup_subtype[1,"x"]-0.5 )
poly_y = c(dup_dup_subtype[1,"y"]-0.498,dup_dup_subtype[1,"y"]-0.498,dup_dup_subtype[1,"y"]+0.498,dup_dup_subtype[1,"y"]+0.498)
polygonID_threeOverlap_subtype <- data.frame(group=rep(1,times=4),poly_x,poly_y)



pdf("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype/survival_stage_subtype_merge_atleast_1gene.pdf",width=10,height = 11)#,width=4000,height=1500,res=400)
ggplot(unique_data.mm,aes(x=tumor,y=geneSymbol))+
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
dev.off()

###Add age data
age <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/clock.sig.txt",header=T)


###Example PER3
clock.exp.stage <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.stage/clock.exp.stage.txt",header=T)
PER3_stage <- clock.exp.stage[which(clock.exp.stage$GeneSymbol=="PER3" & clock.exp.stage$type=="KIRC"),]
library(RColorBrewer)
fill_colors <- brewer.pal(4,"Set1")
pdf("PER3_KIRC_ExpForStage_withTitle.pdf",width = 4.5,height = 5)
ggplot(PER3_stage, aes( x = pathologic_stage, y = log2(value), fill = pathologic_stage))+ stat_boxplot(geom = "errorbar", width = 0.15 )+
     geom_boxplot(width = 0.4,outlier.shape = NA)+
     scale_fill_manual( limits= c("Stage I","Stage II","Stage III","Stage IV"),values = rainbow(4),labels=c("I","II","III","IV"),guide=F)+
     scale_x_discrete(limits= c("Stage I","Stage II","Stage III","Stage IV"),labels=c("I","II","III","IV"))+
     scale_y_continuous(limits=c(6.8,11.8),breaks=seq(7,11))+
     labs(list(title =paste("PER3",' (FDR = ',stage[which(stage$tumor=="KIRC" & stage$geneSymbol=="PER3"),"FDR"],")",sep="") , x = "Stages", y = 'Gene Expression'))+
     #annotate("text",x=1.3,y=7,label=paste('FDR = ',stage[which(stage$tumor=="KIRC" & stage$geneSymbol=="PER3"),"FDR"],sep=""),size=6)+
      theme(panel.background = element_rect(fill=NA,color="black",size=1),
            panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
         axis.text.x = element_text(size = 16,color = "black"),
         axis.text.y = element_text(size = 16,color = "black"),
         axis.title = element_text(size = 18),
         axis.ticks = element_line(color = "black"),
         axis.ticks.length = unit(0.25, "cm"),
        plot.title=element_text(size = 18,hjust=0.5))
dev.off()


clock.exp.subtype <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/clock.exp.in.different.subtype.txt")
PER3_subtype <- clock.exp.subtype[which(clock.exp.subtype$GeneSymbol=="PER3" & clock.exp.subtype$type=="KIRC"),]
subtype_sign <- read.csv("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/significant.DE.clock.genes.in.different.subtype.csv",header=T)
library(RColorBrewer)
#fill_colors <- brewer.pal(4,"Set1")
pdf("PER3_KIRC_ExpForsubtype_withTitle.pdf",width = 4.5,height = 5)
ggplot(PER3_subtype, aes( x = subtype, y = log2(Exp_val), fill = subtype))+ stat_boxplot(geom = "errorbar", width = 0.15 )+
  geom_boxplot(width = 0.4,outlier.shape = NA)+
  scale_fill_manual( limits= c("Normal",1,2,3,4),values = rainbow(5),labels= c("Normal",1,2,3,4),guide=F)+
  scale_x_discrete(limits=  c("Normal",1,2,3,4))+
  scale_y_continuous(limits=c(6.8,11.8),breaks=seq(7,11))+
  labs(list(title =paste("PER3",' (FDR = ',subtype_sign[which(subtype_sign$tumor=="KIRC" &subtype_sign$geneSymbol=="PER3"),"FDR"],")",sep="") , x = "Subtypes", y = 'Gene Expression'))+
  #annotate("text",x=1.3,y=7,label=paste('FDR = ',subtype[which(subtype$tumor=="KIRC" & subtype$geneSymbol=="PER3"),"FDR"],sep=""),size=6)+
  theme(panel.background = element_rect(fill=NA,color="black",size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 16,color = "black"),
        axis.text.y = element_text(size = 16,color = "black"),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        plot.title=element_text(size = 18,hjust=0.5))
dev.off()
###gene expression of CRY2 in different stages for KIRP
CRY2_stage <- clock.exp.stage[which(clock.exp.stage$GeneSymbol=="CRY2" & clock.exp.stage$type %in% c("KIRP")),]
CRY2_FDR <- data.frame(type=c("KIRP","LIHC"),FDR=stage[which(stage$tumor %in% c("KIRP") & stage$geneSymbol=="CRY2"),"FDR"])
library(RColorBrewer)
fill_colors <- brewer.pal(4,"Set1")
pdf("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype/CRY2_stage_KIRP.pdf",width = 4,height = 3.5)
ggplot(CRY2_stage, aes( x = pathologic_stage, y = log2(value), fill = pathologic_stage))+ stat_boxplot(geom = "errorbar", width = 0.15 )+
  geom_boxplot(width = 0.6,outlier.shape = NA)+
  scale_fill_manual( limits= c("Stage I","Stage II","Stage III","Stage IV"),values = brewer.pal(4,"Set1"),labels=c("I","II","III","IV"),guide=F)+
  scale_x_discrete(limits= c("Stage I","Stage II","Stage III","Stage IV"),labels=c("I","II","III","IV"))+
 # scale_y_continuous(limits=c(6.8,11.8),breaks=seq(7,11))+
  labs(list(title =paste("KIRP-CRY2",' (FDR = ',stage[which(stage$tumor=="KIRP" & stage$geneSymbol=="CRY2"),"FDR"],")",sep="") , x = "Stages", y = 'Gene Expression'))+
  # facet_wrap(~type,ncol=2)+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=16),
        plot.title=element_text(size=18,hjust=0.5),
        strip.background = element_rect(fill=NA))
dev.off()
pdf("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype/CRY2_subtype_KIRC.pdf",width = 4,height = 3.5)
ggplot(CRY2_subtype, aes( x = subtype, y = log2(Exp_val), fill = subtype))+ stat_boxplot(geom = "errorbar", width = 0.15 )+
  geom_boxplot(width = 0.4,outlier.shape = NA)+
  scale_fill_manual( limits= c(1,2,3,4),values = brewer.pal(5,"Set1"),labels= c(1,2,3,4),guide=F)+
  scale_x_discrete(limits=  c(1,2,3,4))+
  #scale_y_continuous(limits=c(6.8,11.8),breaks=seq(7,11))+
  labs(list(title =paste("KRIC-CRY2",' (FDR = ',subtype_sign[which(subtype_sign$tumor=="KIRC" &subtype_sign$geneSymbol=="CRY2"),"FDR"],")",sep="") , x = "Subtypes", y = 'Gene Expression'))+
  #annotate("text",x=1.3,y=7,label=paste('FDR = ',subtype[which(subtype$tumor=="KIRC" & subtype$geneSymbol=="CRY2"),"FDR"],sep=""),size=6)+
  # facet_wrap(~type,ncol=2)+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.text.x=element_text(size=16,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=16),
        plot.title=element_text(size=18,hjust=0.5),
        strip.background = element_rect(fill=NA))
dev.off()
##
for ( i in unique(clock.exp.stage$type)){
  sub_stage <- clock.exp.stage[which(clock.exp.stage$type == i),]
  for(j in unique(sub_stage$GeneSymbol)){
    sub_stageGene <-  sub_stage[which( sub_stage$GeneSymbol == j),]
    GeneExp <- sapply(split(sub_stageGene[,"value"], sub_stageGene$pathologic_stage),mean)
    GeneExp <- data.frame(t(GeneExp))
    GeneExp$GeneSymbol <- j
    if(GeneExp[,1] >  1.5 * GeneExp[,4]  ){
      GeneExp$Trend <- "down"
    }else if( GeneExp[,1] *1.5 < GeneExp[,4]){
      GeneExp$Trend <- "up"
    }else{
      GeneExp$Trend <- "none"
    }
    if(j ==  unique(sub_stage$GeneSymbol)){
      GeneExpA <- GeneExp
    }else{
      GeneExpA <- rbind(GeneExpA,GeneExp)
    }
  }
  GeneExpA$type <- rep(i,times=nrow(GeneExpA))
  if(i == unique(clock.exp.stage$type)[1]){
    GeneExpAll <-  GeneExpA
  }else{
    GeneExpAll <- rbind(GeneExpAll,GeneExpA)
  }
}

down <- GeneExpAll[which(GeneExpAll$Trend=="down"),]
up <- GeneExpAll[which(GeneExpAll$Trend=="up"),]
down <- merge(down,stage,by.x=c("GeneSymbol","type"),by.y=c("geneSymbol","tumor"))
up <- merge(up,stage,by.x=c("GeneSymbol","type"),by.y=c("geneSymbol","tumor"))



CRY2_subtype <- clock.exp.subtype[which(clock.exp.subtype$GeneSymbol=="CRY2" & clock.exp.subtype$type=="KIRC"),]
#subtype_sign <- read.csv("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.subtype/significant.DE.clock.genes.in.different.subtype.csv",header=T)
library(RColorBrewer)
#fill_colors <- brewer.pal(4,"Set1")


###Subtype for BRCA
BRCA_subtype <- subtype_sign[which(subtype_sign$tumor=="BRCA"),]
for(i in as.character(BRCA_subtype$geneSymbol)){
  gene_subtype <- clock.exp.subtype[which(clock.exp.subtype$GeneSymbol==i & clock.exp.subtype$type=="BRCA"),]
  pdf(paste("subtype_BRCA/",i,"_BRCA_ExpForsubtype_withTitle.pdf",sep=""),width = 4.5,height = 5)
  P <-  ggplot(gene_subtype, aes( x = subtype, y = log2(Exp_val), fill = subtype))+ stat_boxplot(geom = "errorbar", width = 0.15 )+
    geom_boxplot(width = 0.4,outlier.shape = NA)+
    scale_fill_manual( limits= as.character(unique(gene_subtype$subtype)),values = rainbow(length(unique(gene_subtype$subtype))),labels= c("Normal",1,2,3,4),guide=F)+
    scale_x_discrete(limits=  unique(gene_subtype$subtype))+
    #scale_y_continuous(limits=c(6.8,11.8),breaks=seq(7,11))+
    labs(list(title =paste(i,' (FDR = ',subtype_sign[which(subtype_sign$tumor=="BRCA" &subtype_sign$geneSymbol==i),"FDR"],")",sep="") , x = "Subtypes", y = 'Gene Expression'))+
    #annotate("text",x=1.3,y=7,label=paste('FDR = ',subtype[which(subtype$tumor=="KIRC" & subtype$geneSymbol=="PER3"),"FDR"],sep=""),size=6)+
    theme(panel.background = element_rect(fill=NA,color="black",size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 16,color = "black",angle=90,vjust=0.5,hjust=1),
          axis.text.y = element_text(size = 16,color = "black"),
          axis.title = element_text(size = 18),
          axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(0.25, "cm"),
          plot.title=element_text(size = 18,hjust=0.5))
  print(P)
 dev.off()
}
BRCA_subtype.exp <- clock.exp.subtype[which(clock.exp.subtype$type=="BRCA" & clock.exp.subtype$subtype %in% c("Basal","Her2","LumA","LumB")),]
BRCA_subtype.exp.Basal <- BRCA_subtype.exp[which(BRCA_subtype.exp$subtype == "Basal"),]
BRCA_subtype.exp.Basal.mean <- sapply(split(BRCA_subtype.exp.Basal[,"Exp_val"],BRCA_subtype.exp.Basal$GeneSymbol),mean)
BRCA_subtype.exp.Other <- BRCA_subtype.exp[which(BRCA_subtype.exp$subtype != "Basal"),]
BRCA_subtype.exp.Other.mean <- sapply(split(BRCA_subtype.exp.Other[,"Exp_val"],BRCA_subtype.exp.Other$GeneSymbol),mean)
Basal.sp <-  names(BRCA_subtype.exp.Other.mean)[abs(log2((BRCA_subtype.exp.Other.mean+1)/(BRCA_subtype.exp.Basal.mean+1))) >= log2(1.5)]
Basal.sp.down <-  names(BRCA_subtype.exp.Other.mean)[log2((BRCA_subtype.exp.Other.mean)/(BRCA_subtype.exp.Basal.mean)) >= log2(1.5)]
Basal.sp.up <-  names(BRCA_subtype.exp.Other.mean)[log2((BRCA_subtype.exp.Other.mean)/(BRCA_subtype.exp.Basal.mean)) <= -log2(1.5)]
Basal.sp.up.exp <- BRCA_subtype.exp[which(BRCA_subtype.exp$GeneSymbol %in% Basal.sp.up),]
Basal.sp.up.exp <- Basal.sp.up.exp [which(Basal.sp.up.exp$GeneSymbol != "FBXL21"),]
pdf("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype/subtype_BRCA/Basal.sp.up.expression.pdf",width = 4,height = 5)
ggplot(Basal.sp.up.exp, aes( x = subtype, y = log2(Exp_val+1), fill = subtype))+ stat_boxplot(geom = "errorbar", width = 0.2 )+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  scale_fill_manual( limits= c("Basal","LumA","LumB","Her2"),values = c("red","darkgray","gray","lightgray"),guide=F)+
  scale_x_discrete(limits= c("Basal","LumA","LumB","Her2"),labels= c("TNBC","LumA","LumB","Her2"))+
 # geom_text(data = LGG_sig_Pval,aes(y = 11.5,x=1, label = paste("p = ",t.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
  scale_y_continuous(breaks=seq(0,12,by=4))+
  facet_wrap(~GeneSymbol,ncol=3,scale="free_y")+
  ylab('Gene Expression')+xlab("")+
  theme(panel.background=element_rect(colour="black",fill="white",size=0.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=10,color="black"),axis.title=element_text(size=14,color="black"),
        axis.text.x=element_text(size=12,color="black",angle=45,vjust=1,hjust=1),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=12),
        strip.background = element_rect(fill=NA))
dev.off()


Basal.sp.down.exp <- BRCA_subtype.exp[which(BRCA_subtype.exp$GeneSymbol %in% Basal.sp.down),]
#Basal.sp.down.exp <- Basal.sp.down.exp [which(Basal.sp.down.exp$GeneSymbol != "FBXL21"),]
pdf("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype/subtype_BRCA/Basal.sp.down.expression.pdf",width = 4,height = 5)
ggplot(Basal.sp.down.exp, aes( x = subtype, y = log2(Exp_val+1), fill = subtype))+ stat_boxplot(geom = "errorbar", width = 0.2 )+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  scale_fill_manual( limits= c("Basal","LumA","LumB","Her2"),values = c("blue","darkgray","gray","lightgray"),guide=F)+
  scale_x_discrete(limits= c("Basal","LumA","LumB","Her2"),labels= c("TNBC","LumA","LumB","Her2"))+
  # geom_text(data = LGG_sig_Pval,aes(y = 11.5,x=1, label = paste("p = ",t.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
  scale_y_continuous(breaks=seq(0,12,by=4))+
  facet_wrap(~GeneSymbol,ncol=3,scale="free_y")+
  ylab('Gene Expression')+xlab("")+
  theme(panel.background=element_rect(colour="black",fill="white",size=0.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=10,color="black"),axis.title=element_text(size=14,color="black"),
        axis.text.x=element_text(size=12,color="black",angle=45,vjust=1,hjust=1),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=12),
        strip.background = element_rect(fill=NA))
dev.off()


sub.exp <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/","BRCA_mRNA_each_exp_20160513",sep=""),header=T)
sub.exp$gene <- data.frame(do.call(rbind, strsplit(as.character(sub.exp$gene),'\\|')))$X1
sub.exp <- sub.exp[which(sub.exp$gene %in% Basal.sp),] ##get target genes as.character(BRCA_subtype$geneSymbol)
a <- c("PER2","CRY2","BHLHE40","RORB","BTRC","RCRC","CSNK2A2","NR1I3","NFIL3","ARNTL2","CSNK2B","TNF","NPAS2")
sub.tumor.names <- colnames(sub.exp)[as.numeric(substr(colnames(sub.exp),14,15)) ==1]
sub.tumor.names <- sub.tumor.names[!is.na(sub.tumor.names)]

BRCA_barcode <- unique(clock.exp.subtype[which(clock.exp.subtype$type=="BRCA" & clock.exp.subtype$subtype %in% c("Basal","Her2","LumA","LumB")),])
sub.tumor.names <- sub.tumor.names[gsub("\\.","\\-",substr(sub.tumor.names,1,12)) %in% BRCA_barcode$barcode]
sub.exp.m <- sub.exp[,sub.tumor.names]
colnames(sub.exp.m) <- gsub("\\.","\\-",substr(colnames(sub.exp.m),1,12)) 
rownames(sub.exp.m) <- sub.exp$gene
sub_matrix <- t(apply(log2(sub.exp.m+1),1,scale))
colnames(sub_matrix) <- gsub("\\.","\\-",substr(colnames(sub.exp.m),1,12)) 

sub_col <- unique(BRCA_barcode[,c("barcode","subtype")])
sub_col <- sub_col[order(sub_col$subtype),]
sub_col$subtype <- factor(sub_col$subtype,levels=unique(sub_col$subtype))
sub_matrix <- sub_matrix[,as.character(sub_col$barcode)]
sub_matrix[sub_matrix > 5] <- 5
sub_matrix[sub_matrix < -5] <- -5

sub_col$colors <- rep(brewer.pal(4,"Set1"),times=table(sub_col$subtype))


heatmap.2(sub_matrix[,order(sub_col$subtype)],col = rgb.palette(100),breaks=c(seq(min(sub_matrix),1,length.out = 50),seq(0.01,max(sub_matrix),length.out = 51)),cexRow=1.1,margin=c(0.55,6),distfun = distfunc,
          hclustfun = hclustfunc, Colv = F,lwid=c(0.5,2), lhei=c(0.3,1.2), 
          ColSideColors=sub_col$colors,labCol = FALSE,trace = "none",key=F)
legend(0,1.1,legend=names(table(sub_col$subtype)),
       fill=brewer.pal(6,"Set1"), border=FALSE, bty="n", y.intersp = 1.2, cex=0.9,xpd=T)





core_BRCA <- clock.exp.subtype[which(clock.exp.subtype$GeneSymbol %in% core.circadian & clock.exp.subtype$type=="BRCA"),]
sapply(split(core_BRCA[,"Exp_val"],core_BRCA$subtype),function(x){cor(matrix(unlist(x)),method = "spearman")})
for( i in unique(core_BRCA$subtype)){
  sub <- core_BRCA[which(core_BRCA$subtype==i),]
  sub <- sub[order(sub$GeneSymbol),]
  sub_exp <- matrix(sub$Exp_val,ncol = length(unique(sub$GeneSymbol)))
  colnames(sub_exp) <- unique(sub$GeneSymbol)
  sub_cor<- cor(sub_exp)
  sub_cor <- data.frame(sub_cor)
  sub_cor$subtype <- rep(i,times=nrow(sub_cor))  
  sub_cor$GeneSymbol <- rownames(sub_cor)
  if(i == "Normal"){
    BRCA_cor <- sub_cor
  }else{
    BRCA_cor <- rbind(BRCA_cor,sub_cor)
  }
}

a <- BRCA_cor[which(BRCA_cor$GeneSymbol %in% c("CLOCK","NPAS2","ARNTL","ARNTL2")),c("PER1","PER2","PER3","CRY2","RORA","RORB","RORC","NR1D1","NR1D2","GeneSymbol","subtype")]
image(t(as.matrix(a[order(a$GeneSymbol,a$subtype),1:9])),col=rgb.palette(100),zlim=c(-0.6,0.6))






###Add age data
age <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/clock.sig.txt",header=T)
age <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/clock.sig.txt",header=T)

union_cancer <- Reduce(function(x,y)union(x,y),list(stage$tumor,subtype$tumor,survival$tumor,age$type))
#union_cancer <- union_cancer[!(union_cancer %in% c("DLBC","READ","CESC"))]
union_data <- data.frame(geneSymbol=rep(circadian.genes$V1,each=length(union_cancer)),tumor=rep(union_cancer,times=51))
union_data <- merge(merge(merge(merge(union_data,survival,by=c("tumor","geneSymbol"),all.x=T),
                                stage[,c("tumor","geneSymbol","stage")],by=c("tumor","geneSymbol"),all.x=T),
                          subtype[,c("tumor","geneSymbol","subtype")],by=c("tumor","geneSymbol"),all.x=T),age[,c("tumor","geneSymbol","age")],by=c("tumor","geneSymbol"),all.x=T)
union_data.s <- union_data

union_data.m <- melt(union_data,id.vars=c("tumor","geneSymbol"),measure.vars=c("survival","stage","subtype","age"))
union_data.m <- union_data.m[which(!is.na(union_data.m$value==T)),]
unique_data.m <- union_data.m[!duplicated(union_data.m[,1:2]),]
unique_data.m$geneSymbol <- factor(unique_data.m$geneSymbol,levels=unique(unique_data.m$geneSymbol))
dup <- union_data.m[duplicated(union_data.m[,1:2]),]
unique_dup <- dup[!duplicated(dup[,1:2]),]

unique_dup_1 <- unique_dup[!(unique_dup$tumor == "KIRC" & unique_dup$geneSymbol %in% c("HLF","PRKAA2","TEF","PER3")),]
dup_dup_stage <- unique_dup[which(unique_dup$tumor == "KIRC" & unique_dup$geneSymbol %in% c("HLF","PRKAA2","TEF","PER3")),]
dup_dup_subtype <- dup[duplicated(dup[,1:2]),]

sss_x_label <- union_cancer[order(union_cancer)] #names(table(union_data.m$tumor))[order(table(union_data.m$tumor))] #
sss_y_label <- table(union_data.m$geneSymbol)[table(union_data.m$geneSymbol)>0]
sss_y_label <- names(sss_y_label)[order(sss_y_label)]
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
sss_y_color <- rep("black",times=length(sss_y_label))
sss_y_color[sss_y_label %in% core.circadian] <- "red"


unique_data.mm <- merge(data.frame(geneSymbol=rep(sss_y_label,each=length(union_cancer)),tumor=rep(union_cancer,times=length(sss_y_label))),unique_data.m,by=c("tumor","geneSymbol"),all.x=T)


unique_dup$x = sapply(unique_dup$tumor,function(x){grep(x,sss_x_label)}) ## add the geom_polygon ID, convert the factor to numeric 
unique_dup$y = sapply(unique_dup$geneSymbol,function(x){grep(x,sss_y_label)})#as.numeric(dup.data$geneSymbol)

###stage duplicate
unique_dup_stage <- unique_dup[which(unique_dup$variable=="stage"),]
poly_x = c()
poly_y = c()
for ( i in 1: nrow(unique_dup_stage)){
  poly_x = c(poly_x, unique_dup_stage[i,"x"]-0.498, unique_dup_stage[i,'x'],unique_dup_stage[i,"x"],unique_dup_stage[i,"x"]-0.498 )
  poly_y = c(poly_y, unique_dup_stage[i,"y"]-0.498,unique_dup_stage[i,"y"]-0.498,unique_dup_stage[i,"y"]+0.498,unique_dup_stage[i,"y"]+0.498)
}
polygonID_stage <- data.frame( group = rep(seq(1:nrow(unique_dup_stage)), each = 4),poly_x, poly_y )

###subtype duplicate
unique_dup_subtype <- unique_dup[which(unique_dup$variable=="subtype"),]

poly_x = c()
poly_y = c()
for ( i in 1: nrow(unique_dup_subtype)){
  poly_x = c(poly_x, unique_dup_subtype[i,"x"]-0.498, unique_dup_subtype[i,'x'],unique_dup_subtype[i,"x"],unique_dup_subtype[i,"x"]-0.498 )
  poly_y = c(poly_y, unique_dup_subtype[i,"y"]-0.49,unique_dup_subtype[i,"y"]-0.49,unique_dup_subtype[i,"y"]+0.49,unique_dup_subtype[i,"y"]+0.49)
}
polygonID_subtype <- data.frame( group = rep(seq(1:nrow(unique_dup_subtype)), each = 4),poly_x, poly_y )
###age duplicate
unique_dup_age <- unique_dup[which(unique_dup$variable=="age"),]

poly_x = c()
poly_y = c()
for ( i in 1: nrow(unique_dup_age)){
  poly_x = c(poly_x, unique_dup_age[i,"x"]-0.498, unique_dup_age[i,'x'],unique_dup_age[i,"x"],unique_dup_age[i,"x"]-0.498 )
  poly_y = c(poly_y, unique_dup_age[i,"y"]-0.49,unique_dup_age[i,"y"]-0.49,unique_dup_age[i,"y"]+0.49,unique_dup_age[i,"y"]+0.49)
}
polygonID_age <- data.frame( group = rep(seq(1:nrow(unique_dup_age)), each = 4),poly_x, poly_y )


###three ovelap ---stage
dup_dup_stage$x <- grep(dup_dup_stage$tumor,sss_x_label)
dup_dup_stage$y <- grep(dup_dup_stage$geneSymbol,sss_y_label)

poly_x = c(dup_dup_stage[1,"x"]-0.498+1/3, dup_dup_stage[1,'x']-0.498+2/3,dup_dup_stage[1,"x"]-0.498+2/3,dup_dup_stage[1,"x"]-0.498+1/3 )
poly_y = c(dup_dup_stage[1,"y"]-0.498,dup_dup_stage[1,"y"]-0.498,dup_dup_stage[1,"y"]+0.498,dup_dup_stage[1,"y"]+0.498)
polygonID_threeOverlap_stage <- data.frame(group=rep(1,times=4),poly_x,poly_y)

###three ovelap ---subtype
dup_dup_subtype$x <- grep(dup_dup_subtype$tumor,sss_x_label)
dup_dup_subtype$y <- grep(dup_dup_subtype$geneSymbol,sss_y_label)

poly_x = c(dup_dup_subtype[1,"x"]-0.5, dup_dup_subtype[1,'x']-0.498+1/3,dup_dup_subtype[1,"x"]-0.498+1/3,dup_dup_subtype[1,"x"]-0.5 )
poly_y = c(dup_dup_subtype[1,"y"]-0.498,dup_dup_subtype[1,"y"]-0.498,dup_dup_subtype[1,"y"]+0.498,dup_dup_subtype[1,"y"]+0.498)
polygonID_threeOverlap_subtype <- data.frame(group=rep(1,times=4),poly_x,poly_y)



pdf("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype/survival_stage_subtype_age_merge_atleast_1gene.pdf",width=10,height = 11)#,width=4000,height=1500,res=400)
ggplot(unique_data.mm,aes(x=tumor,y=geneSymbol))+
  geom_tile(aes(fill=factor(value)),col="darkgray")+
  scale_fill_manual(limits=c("High","Low","Stage","Subtype","Young","Old"),values = c("red","blue","green","gold","cyan","magenta"),na.value="white",labels=c("High_level_Worse","Low_level_Worse","Stage","Subtype","Young","Old"),name="")+
  scale_y_discrete(limit=sss_y_label)+
  scale_x_discrete(limit=sss_x_label)+
  geom_polygon(data=polygonID_stage,aes(x=poly_x,y=poly_y,group=group),fill="green")+
  geom_polygon(data=polygonID_subtype,aes(x=poly_x,y=poly_y,group=group),fill="gold")+
  geom_polygon(data=polygonID_age,aes(x=poly_x,y=poly_y,group=group),fill="magenta")+
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
dev.off()


#######Breast cancer pathway score distribution
###
BRCA.pathwayScore <- read.csv("/extraspace/yye1/analysis/Circadian/RPPA/RPPARelated/BRCA-Pathway-Scores.csv")
BRCA.clinical <- clock.exp.subtype[which(clock.exp.subtype$type == "BRCA"),]
BRCA.pathwayScore_subtype <- merge(unique(BRCA.clinical[,c(1,4)]),BRCA.pathwayScore,by.x="barcode",by.y="Sample")
BRCA.pathwayScore_subtype <- BRCA.pathwayScore_subtype[which(BRCA.pathwayScore_subtype$subtype %in%  c("Basal","LumA","LumB","Her2")),]
BRCA.pathwayScore_subtype.m <- melt(BRCA.pathwayScore_subtype,id.vars = "subtype",measure.vars = colnames(BRCA.pathwayScore_subtype)[4:13])
BRCA.pathwayScore_subtype.m$variable <- gsub("Score","",BRCA.pathwayScore_subtype.m$variable)
pdf("/extraspace/yye1/analysis/Circadian/expression/sur_stage_subtype/subtype_BRCA/BRCA_subtype_pathway_scores.pdf",width = 9,height = 4.8)
ggplot(BRCA.pathwayScore_subtype.m, aes( x = subtype, y =  value, color = subtype))+ stat_boxplot(geom = "errorbar", width = 0.2 )+
  geom_boxplot(width = 0.6,fill=NA,outlier.shape = NA)+#)+
  scale_color_manual( limits= c("Basal","LumA","LumB","Her2"),values = rainbow(4),guide=F)+
  scale_x_discrete(limits= c("Basal","LumA","LumB","Her2"),labels= c("Basal","LumA","LumB","Her2"))+
  # geom_text(data = LGG_sig_Pval,aes(y = 11.5,x=1, label = paste("p = ",t.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
  # scale_y_continuous(breaks=seq(0,12,by=4))+
  facet_wrap(~variable,ncol=5)+
  ylab('Gene Expression')+xlab("")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text.y=element_text(size=14,color="black"),axis.title=element_text(size=18,color="black"),
        axis.text.x=element_text(size=14,color="black",angle=30,vjust=1,hjust=1),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=14),
        strip.background = element_rect(fill=NA))
dev.off()

####TNBC VS other samples for correlation with other genes
BRCA <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/","BRCA_mRNA_each_exp_20160513",sep=""),header=T)
BRCA$gene <- data.frame(do.call(rbind,strsplit(as.character(BRCA$gene),"\\|")))$X1
TumorSampleNames <- colnames(BRCA)[as.numeric(substr(colnames(BRCA),14,15))==1]
TumorSampleNames <- TumorSampleNames[!is.na(TumorSampleNames)]
BRCA <- BRCA[,c("gene",TumorSampleNames )]
colnames(BRCA) <- gsub("\\.","\\-",substr(colnames(BRCA),1,12))
TNBC_sampleNames <- as.character(unique(BRCA_subtype.exp[which(BRCA_subtype.exp$subtype=="Basal"),"barcode"]))
OtherType_sampleNames <- as.character(unique(BRCA_subtype.exp[which(BRCA_subtype.exp$subtype %in% c("LumA","LumB","Her2")),"barcode"]))
OtherType_sampleNames_samples <- sample(OtherType_sampleNames,length(TNBC_sampleNames))
TNBC.corAll <- data.frame(gene=BRCA$gene)
OtherType.corAll <- data.frame(gene=BRCA$gene)
for(i in circadian.genes){
  TNBC.exp <- BRCA[which(BRCA$gene == i),TNBC_sampleNames]
  OtherType.exp <- BRCA[which(BRCA$gene == i),OtherType_sampleNames_samples  ]
  TNBC.cor <- apply(log2(as.matrix(BRCA[,TNBC_sampleNames])+1) , 1 ,function(x){cor.test(x,as.vector(as.numeric(TNBC.exp)),method="pearson")$estimate})
  OtherType.cor <- apply(log2(as.matrix(BRCA[,OtherType_sampleNames_samples ])+1) , 1 ,function(x){cor.test(x,as.vector(as.numeric(OtherType.exp)),method="pearson")$estimate})
  TNBC.pvalue <- apply(log2(as.matrix(BRCA[,TNBC_sampleNames])+1) , 1 ,function(x){cor.test(x,as.vector(as.numeric(TNBC.exp)),method="pearson")$p.value})
  OtherType.pvalue <- apply(log2(as.matrix(BRCA[,OtherType_sampleNames_samples ])+1) , 1 ,function(x){cor.test(x,as.vector(as.numeric(OtherType.exp)),method="pearson")$p.value})
  TNBC.FDR <- p.adjust(TNBC.pvalue,method="fdr")
  OtherType.FDR <- p.adjust(OtherType.pvalue,method="fdr")
  TNBC.cor_sub <- data.frame(TNBC.cor,TNBC.pvalue,TNBC.FDR)
  colnames(TNBC.cor_sub) <- paste(i,c("","_pvalue","_FDR"),sep="")
  OtherType.cor_sub <- data.frame(OtherType.cor,OtherType.pvalue,OtherType.FDR)
  colnames(OtherType.cor_sub) <- paste(i,c("","_pvalue","_FDR"),sep="")
  TNBC.corAll <- cbind(TNBC.corAll,TNBC.cor_sub)
  OtherType.corAll <- cbind(OtherType.corAll,OtherType.cor_sub)
}


for(j in circadian.genes){
  OtherType.corNum <- nrow(OtherType.corAll[which(OtherType.corAll[j] > 0.3 & OtherType.corAll[paste(j,"_FDR",sep="")] < 0.05),])
  TNBC.corNum <- nrow(TNBC.corAll[which(TNBC.corAll[j] > 0.3 & TNBC.corAll[paste(j,"_FDR",sep="")] < 0.05),])
  corNum <- data.frame(gene=j,TBNC_CorNum=TNBC.corNum, OtherType_CorNum=OtherType.corNum)
   if(j ==circadian.genes[1]){
    corNumAll <- corNum
   
  }else{
    corNumAll <- rbind(corNumAll,corNum)
  }
}
for( m in 1:100){
  set.seed(m)
  OtherType_sampleNames_samples <- sample(OtherType_sampleNames,length(TNBC_sampleNames))
  OtherType.corAll <- data.frame(gene=BRCA$gene)
  for(i in circadian.genes){
    OtherType.exp <- BRCA[which(BRCA$gene == i),OtherType_sampleNames_samples  ]
     OtherType.cor <- apply(log2(as.matrix(BRCA[,OtherType_sampleNames_samples ])+1) , 1 ,function(x){cor.test(x,as.vector(as.numeric(OtherType.exp)),method="pearson")$estimate})
    OtherType.pvalue <- apply(log2(as.matrix(BRCA[,OtherType_sampleNames_samples ])+1) , 1 ,function(x){cor.test(x,as.vector(as.numeric(OtherType.exp)),method="pearson")$p.value})
    OtherType.FDR <- p.adjust(OtherType.pvalue,method="fdr")
     OtherType.cor_sub <- data.frame(OtherType.cor,OtherType.pvalue,OtherType.FDR)
    colnames(OtherType.cor_sub) <- paste(i,c("","_pvalue","_FDR"),sep="")
    OtherType.corAll <- cbind(OtherType.corAll,OtherType.cor_sub)
  }
  for(j in circadian.genes){
    OtherType.corNum <- nrow(OtherType.corAll[which(OtherType.corAll[j] > 0.3 & OtherType.corAll[paste(j,"_FDR",sep="")] < 0.05),])
    corNum <- data.frame(gene=j, OtherType_CorNum=OtherType.corNum)
    if(j ==circadian.genes[1]){
      OthercorNumAll <- corNum
      
    }else{
      OthercorNumAll <- rbind(OthercorNumAll,corNum)
    }
  }
  if(m==1){
    OthercorNumAllSum <- OthercorNumAll
  }else{
    OthercorNumAllSum <- merge(OthercorNumAllSum,OthercorNumAll,by="gene")
  }
}

#### circadian genes expression across different tumor stage
setwd("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.stage")
###mRNA expression data from TCGA mRNA expression
exp.files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
clinic.files.names <- list.files(path="/extraspace/TCGA/TCGA_clinical/",pattern = "_clinical_clean.txt")
commonsamples <- intersect(gsub("_clinical_clean.txt","",clinic.files.names),gsub("_mRNA_each_exp_20160513","",exp.files.names))
library(RColorBrewer)
fill_colors <- brewer.pal(4,"Set1")
i=0
clock.exp.stage <- data.frame()
aov.PvalAll <- data.frame()
for(cancer in commonsamples){
  clinic <- read.delim(paste("/extraspace/TCGA/TCGA_clinical/",cancer, "_clinical_clean.txt",sep=""),header=T)
  if(!is.na(match("pathologic_stage",colnames(clinic)))){
    clinic <- clinic[which(clinic$pathologic_stage %in% c("Stage I","Stage II","Stage III","Stage IV")),]
    clinic$pathologic_stage <- factor(clinic$pathologic_stage,levels=c("Stage I","Stage II","Stage III","Stage IV"))
    if(sum(table(clinic$pathologic_stage)) >=40){
      exp <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",cancer,"_mRNA_each_exp_20160513",sep=""),header=T)
      exp <- exp[which(data.frame(do.call(rbind, strsplit(as.character(exp$gene),'\\|')))$X1 %in% circadian.genes$V1),]
      commonbarcode <- intersect(clinic$barcode,gsub("\\.","\\-",substr(colnames(exp),1,12)))
      commonexpnames <- colnames(exp)[gsub("\\.","\\-",substr(colnames(exp),1,12)) %in% commonbarcode]
      commonexpnames <- commonexpnames[!duplicated(gsub("\\.","\\-",substr(commonexpnames,1,12)))]
      exp <- exp[,c("gene",commonexpnames)]
      colnames(exp) <- c("gene",gsub("\\.","\\-",substr(commonexpnames,1,12)))
      exp$GeneSymbol <- data.frame(do.call(rbind,strsplit(as.character(exp$gene),'\\|')))$X1
      exp.m <- melt(exp,id.vars = "GeneSymbol",measure.vars =gsub("\\.","\\-",substr(commonexpnames,1,12)))
      exp.mm <- merge(exp.m,clinic[,c("barcode","pathologic_stage")],by.x="variable",by.y="barcode")
      aov.Pval <- signif(sapply(split(exp.mm[,3:4],exp.mm$GeneSymbol),function(x){oneway.test(log2(value+1)~pathologic_stage,data = x)$p.value}),digits = 2)
      aov.Pval <- data.frame(aov.Pval)
      aov.Pval$GeneSymbol <- rownames(aov.Pval)
      
      #pdf(paste(cancer,".clock.genes.expression.across.different.stage.pdf",sep=""),width = 16,height = 12)
      #P1 <- ggplot(exp.mm, aes( x = pathologic_stage, y = log2(value+1), fill = pathologic_stage))+ stat_boxplot(geom = "errorbar", width = 0.2 )+
     #   geom_boxplot(width = 0.6,outlier.shape = NA)+
     #   scale_fill_manual( limits= c("Stage I","Stage II","Stage III","Stage IV"),values = fill_colors,labels=c("I","II","III","IV"),guide=F)+
     #   scale_x_discrete(limits= c("Stage I","Stage II","Stage III","Stage IV"),labels=c("I","II","III","IV"))+
     #   geom_text(data = aov.Pval,aes(y = 14,x=1.2, label = paste("p=",aov.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
     #   facet_wrap(~GeneSymbol,nrow=4)+
    #    #  geom_line( aes(group = cancer), linetype = 'longdash',position=position_jitter(w=0, h=0), size = 0.3)+
        #  geom_jitter(position = position_jitter(width = 0.05, height = 0), size = 0.6)+
     #   ylab('Gene Expression')+xlab("")+
        # ggtitle(paste0('Wilcox Pvalue: ',pvalue ))+
    #    theme(panel.background = element_rect(fill=NA,color="black",size=2),
    #          panel.grid.major = element_blank(),
    #          panel.grid.minor = element_blank(),
     #         axis.text.x = element_text(size = 16),
     #         axis.text.y = element_text(size = 16),
     #         axis.title.y = element_text(size = 18),
    #          axis.text.x = element_blank(),
    #          # plot.title=element_text(size = 12),
     #         strip.text=element_text(size=16))
     # print(P1)
    #  dev.off()
      exp.mm$type <- rep(cancer,times=nrow(exp.mm))
      aov.Pval$type <- rep(cancer,times=nrow(aov.Pval))
      if(nrow(clock.exp.stage)==0){
        clock.exp.stage <- exp.mm
        
      }else{
        clock.exp.stage <- rbind(clock.exp.stage,exp.mm)
      }
      
      if(nrow(aov.PvalAll)==0){
        aov.PvalAll <- aov.Pval
      }else{aov.PvalAll <- rbind(aov.PvalAll,aov.Pval)}
    }
    
  }
}

write.table(clock.exp.stage,file="clock.exp.stage.txt",sep="\t",quote = F,row.names = F) 
write.table(aov.PvalAll,file="clock.exp.stage.aov.Pvalue.txt",sep="\t",quote = F,row.names = F)


sur_all.p <- read.delim("~/Circadian/survival/circadian.genes.survival.count.txt",header=T)
sur_all.p <-  sur_all.p[!( sur_all.p$geneSymbol %in% c("FOXO1","FOXM1","FOXA2","FOXO3")),]

sur_all.p <- sur_all.p[!is.na(sur_all.p$survival),]
sur_all.p$geneSymbol <- factor(sur_all.p$geneSymbol,levels = unique(sur_all.p$geneSymbol))
sur_all.p$tumor <- factor(sur_all.p$tumor,levels=unique(sur_all.p$tumor))
x_label <- names(table(sur_all.p$tumor))[order(table(sur_all.p$tumor))]
y_label <- names(table(sur_all.p$geneSymbol))[order(table(sur_all.p$geneSymbol))]
miss <- data.frame(geneSymbol=rep(unique(sur_all.p$geneSymbol),times=length(unique(sur_all.p$tumor))),tumor=rep(unique(sur_all.p$tumor),each=length(unique(sur_all.p$geneSymbol))),survival=rep(0,times=length(unique(sur_all.p$geneSymbol))*length(unique(sur_all.p$tumor))))
miss <- merge(miss,sur_all.p,by=c("geneSymbol","tumor"),all.x=T) 
miss <- miss[(is.na(miss$survival.y)==T),1:3]
colnames(miss) <- colnames(sur_all.p)
sur_all.pp <- rbind(sur_all.p,miss)
pdf("~/Circadian/survival/circadian.genes.survival.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(sur_all.pp,aes(x=tumor,y=geneSymbol))+
  geom_tile(aes(fill=factor(survival)),col="white")+
  scale_fill_manual(values = c("High"="red","Low" = "blue"),labels=c("H_Poor_Sur","L_Poor_Sur"),na.value="lightgray",name="Survival")+
  scale_y_discrete(limit=y_label)+
  scale_x_discrete(limit=x_label)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=16,colour = "black"),
        axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))

dev.off()
aov.PvalAll <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.stage/clock.exp.stage.aov.Pvalue.txt",header=T)
aov.PvalAll <- aov.PvalAll[,c(2,3,1)]
colnames(aov.PvalAll) <- c("geneSymbol","tumor","Pval")
aov.PvalAll$geneSymbol <- factor(aov.PvalAll$geneSymbol,levels = unique(aov.PvalAll$geneSymbol))
aov.PvalAll$tumor <- factor(aov.PvalAll$tumor,levels=unique(aov.PvalAll$tumor))
#pmiss <- data.frame(geneSymbol=rep(unique(aov.PvalAll$geneSymbol),times=length(unique(aov.PvalAll$tumor))),tumor=rep(unique(aov.PvalAll$tumor),each=length(unique(aov.PvalAll$geneSymbol))),Pval=rep(0,times=length(unique(aov.PvalAll$geneSymbol))*length(unique(aov.PvalAll$tumor))))
#pmiss <- merge(pmiss,aov.PvalAll,by=c("geneSymbol","tumor"),all.x=T) 
#pmiss <- pmiss[(is.na(pmiss$Pval.y)==T),1:3]
#colnames(pmiss) <- colnames(aov.PvalAll)
#aov.PvalAllp <- rbind(aov.PvalAll,pmiss)
#aov.PvalAllp["Pval"][aov.PvalAllp["Pval"] > 0] <- 1
aov.PvalAll["Pval"][aov.PvalAll["Pval"] <= 0.05] <- 0
aov.PvalAll["Pval"][aov.PvalAll["Pval"] > 0.05] <- 1
aov.PvalAll.sig <- aov.PvalAll[which(aov.PvalAll$Pval==0),]
px_label <- names(table(aov.PvalAll.sig$tumor))[order(table(aov.PvalAll.sig$tumor))]
py_label <- names(table(aov.PvalAll.sig$geneSymbol))[order(table(aov.PvalAll.sig$geneSymbol))]
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
y_color <- rep("black",times=length(py_label))
y_color[py_label %in% core.circadian] <- "red"
pdf("stage.significant.Pval.pdf",width=9,height = 10)#,width=4000,height=1500,res=400)
ggplot(aov.PvalAll,aes(x=tumor,y=geneSymbol))+
  geom_tile(aes(fill=factor(Pval)),col="lightgray")+
  scale_fill_manual(limits=c(1,0),values = c("white","red"),guide=F)+
  scale_y_discrete(limit=py_label)+
  scale_x_discrete(limit=px_label)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=16,colour = y_color),
        axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))
dev.off()
intersect(names(table(aov.PvalAll$type)),names(table(sur_all.p$tumor)))
union(names(table(aov.PvalAll$type)),names(table(sur_all.p$tumor)))

##
clock.exp.stage <- read.delim("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.stage/clock.exp.stage.txt",header=T)


head(clock.exp.stage)
trend <- c()
for( i in unique(aov.PvalAll$tumor)){
  sub <- aov.PvalAll[which(aov.PvalAll$tumor == i),]
  for(j in unique(sub$geneSymbol)){
    sub2 <- clock.exp.stage[which(clock.exp.stage$GeneSymbol==as.character(j) & clock.exp.stage$type==as.character(i)),]
    mean_val <- sapply(split(sub2[,"value"],sub2$pathologic_stage),mean)
    FC <- log2((max(mean_val)+1)/(min(mean_val)+1))
  # print(c(j,FC,mean_val))}
    pval <- t.test(sub2[which(sub2$pathologic_stage == "Stage I"),"value"],sub2[which(sub2$pathologic_stage == "Stage IV"),"value"])$p.value
    pval[is.nan(pval)==T] <- 1
    aovPval <- oneway.test(value~pathologic_stage,data = sub2)$p.value
    aovPval[is.nan(aovPval) == T] <- 1
   # aovPval[is.nan(aovPval)==T] <- 1
    
    if(aovPval <= 0.05 & abs(FC) >= 0.5){
      if(log2((mean_val["Stage I"]+1)/(mean_val["Stage IV"]+1)) >= 0.5 & pval <= 0.05){
        st <- c(i,j,"Down")
      }
      else if(log2((mean_val["Stage I"]+1)/(mean_val["Stage IV"]+1)) <= -0.5 & pval <= 0.05){
        st <- c(i,j,"Up")
      }else{
        st <- c(i,j,"Mix")
      }
      
      trend <- rbind(trend,st)
    }
    }
}



library(varhandle)
trend <- data.frame( trend)
colnames(trend) <- c("tumor","geneSymbol","change")
aov.PvalAllpp <- merge(aov.PvalAll,trend,by=c("geneSymbol","tumor"),all.x=T)
write.table(aov.PvalAllpp,file="Exp.dif.stage.log2F_0.5_pval0.5.txt",quote = F,sep="\t",row.names = F)
#aov.PvalAllpp$change <- unfactor(aov.PvalAllpp$change)
#aov.PvalAllpp["change"][is.na(aov.PvalAllpp["change"])==T] <- "None"
pdf("stage.significant.Pval_with.expression.trend.pdf",width=9,height = 10)#,width=4000,height=1500,res=400)
ggplot(aov.PvalAllpp,aes(x=tumor,y=geneSymbol))+
  geom_tile(aes(fill=factor(change)),col="lightgray")+
  scale_fill_manual(limits=c("Up","Down","Mix"),values = c("red","blue","gold"),labels=c("Late","Early","Differential expression"),na.value="white",name="")+
  scale_y_discrete(limit=py_label)+
  scale_x_discrete(limit=px_label)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=16,colour = y_color),
        axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),legend.position="bottom",legend.direction="horizontal")

dev.off()



####
trend_stage <- c()
for( i in as.character(unique(clock.exp.stage$type))){
  if(nrow(clock.exp.stage[which( clock.exp.stage$type==i),]) >= 51*20){
    for(j in unique(clock.exp.stage$GeneSymbol)){
      sub <- clock.exp.stage[which(clock.exp.stage$GeneSymbol==as.character(j) & clock.exp.stage$type==i),]
      sub$pathologic_stage <- factor(sub$pathologic_stage,levels = unique(sub$pathologic_stage))
      mean_val <- sapply(split(sub[,"value"],sub$pathologic_stage),mean)
      FC <- log2((max(mean_val)+1)/(min(mean_val)+1))
      aovPval <- oneway.test(value~pathologic_stage,data = sub)$p.value
      aovPval[is.nan(aovPval) == T] <- 1
      #  aovPval <- signif(aovPval,digits = 2)
      trend_stage <- rbind(trend_stage,c(i,j,FC,aovPval))
      
    }
  }
}
library(varhandle)
trend_stage <- data.frame( trend_stage)
colnames(trend_stage) <- c("tumor","geneSymbol","FC_log2","p_val")
trend_stage$FDR <- signif(p.adjust(as.numeric(unfactor(trend_stage$p_val)),method = "fdr"),digits = 2)
trend_stage_sign <- trend_stage[which(abs(as.numeric(unfactor(trend_stage$FC_log2)))>=0.5 & trend_stage$FDR <= 0.05),]
write.csv(trend_stage_sign,file="significant.DE.clock.genes.in.different.stage.csv",quote = F,row.names = F)

tumor.stage <- unique(clock.exp.stage[,c("variable","pathologic_stage","type")])
tumor.stage_sum <- sapply(split(tumor.stage[,"pathologic_stage"],tumor.stage$type),table)
tumor.stage_sum <- data.frame(tumor.stage_sum)
tumor.stage_sum$stage <- rownames(tumor.stage_sum)
tumor.stage_sum.m <- melt(tumor.stage_sum,id.vars = "stage",measure.vars=as.character(unique(trend_stage$tumor)))
sapply(split(tumor.stage_sum.m[,"value"],tumor.stage_sum.m$variable),sum)
tumor.stage_sum.m <- tumor.stage_sum.m[which(tumor.stage_sum.m$value !=0),]
colnames(tumor.stage_sum.m) <- c("Stage","Cancer","Samples Number")
tumor.stage_sum.m <- tumor.stage_sum.m[!(tumor.stage_sum.m$Cancer %in% c("SKCM","READ")),c(2,1,3)]
write.csv(tumor.stage_sum.m,file="tumor stage statistics.csv",quote = F,row.names = F)

stage.sig <- stage[which(stage$tumor != "READ"),]
stage.sig <- merge(clock.exp.stage,stage.sig,by.y=c("tumor","geneSymbol"),by.x=c("type","GeneSymbol"))
stage.sig$type_symbol <- paste(stage.sig$type,stage.sig$GeneSymbol,sep=" - ")
stage_sig_Pval <- unique(stage.sig[,c("type_symbol","p_val","FDR")])

pdf("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.stage/expression of Significant DE circadian genes.pdf",width = 16,height = 14)
ggplot(stage.sig, aes( x = pathologic_stage, y = log2(value+1), fill = pathologic_stage))+ stat_boxplot(geom = "errorbar", width = 0.2 )+
  geom_boxplot(width = 0.6,outlier.shape = NA)+
  scale_fill_manual( limits= c("Stage I","Stage II","Stage III","Stage IV"),values = brewer.pal(4,"Set1"),guide=F)+
  scale_x_discrete(limits= c("Stage I","Stage II","Stage III","Stage IV"),labels= c("I","II","III","IV"))+
  geom_text(data = stage_sig_Pval,aes(y = 13.5,x=2, label = paste("p = ",signif(p_val,digits = 2),sep=""),fill=NA),size=4,hjust=0.5)+
  scale_y_continuous(breaks=seq(0,12,by=4))+
  facet_wrap(~type_symbol,ncol=8)+
  ylab('Gene Expression')+xlab("")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=16),
        strip.background = element_rect(fill=NA))
dev.off()


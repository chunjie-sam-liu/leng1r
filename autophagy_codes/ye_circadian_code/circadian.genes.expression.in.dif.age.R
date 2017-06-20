#### circadian genes expression across different tumor age
setwd("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age")
###mRNA expression data from TCGA mRNA expression
exp.files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
clinic.files.names <- list.files(path="/extraspace/TCGA/TCGA_clinical/",pattern = "_clinical_clean.txt")
commonsamples <- intersect(gsub("_clinical_clean.txt","",clinic.files.names),gsub("_mRNA_each_exp_20160513","",exp.files.names))
library(RColorBrewer)
fill_colors <- brewer.pal(4,"Set1")
i=0
clock.exp.age <- data.frame()
t.PvalAll <- data.frame()
for(cancer in commonsamples){
  clinic <- read.delim(paste("/extraspace/TCGA/TCGA_clinical/",cancer, "_clinical_clean.txt",sep=""),header=T)
  if(!is.na(match("age_at_initial_pathologic_diagnosis",colnames(clinic)))){
    clinic <- clinic[which(is.na(clinic$age_at_initial_pathologic_diagnosis)==F),]
#    clinic$pathologic_stage <- factor(clinic$pathologic_stage,levels=c("Stage I","Stage II","Stage III","Stage IV"))
    if(nrow(clinic[which(clinic$age_at_initial_pathologic_diagnosis >= 60),]) >= 10 & nrow(clinic[which(clinic$age_at_initial_pathologic_diagnosis <= 50),]) >= 10 ){
      exp <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",cancer,"_mRNA_each_exp_20160513",sep=""),header=T)
      exp$GeneSymbol <- data.frame(do.call(rbind,strsplit(as.character(exp$gene),'\\|')))$X1
      exp <- exp[which(exp$GeneSymbol %in% circadian.genes$V1),]
      commonbarcode <- intersect(clinic$barcode,gsub("\\.","\\-",substr(colnames(exp),1,12)))
      commonexpnames <- colnames(exp)[gsub("\\.","\\-",substr(colnames(exp),1,12)) %in% commonbarcode]
      commonexpnames <- commonexpnames[!duplicated(gsub("\\.","\\-",substr(commonexpnames,1,12)))]
      exp <- exp[,c("GeneSymbol",commonexpnames)]
      colnames(exp) <- c("GeneSymbol",gsub("\\.","\\-",substr(commonexpnames,1,12)))
     
      exp.m <- melt(exp,id.vars = "GeneSymbol",measure.vars =gsub("\\.","\\-",substr(commonexpnames,1,12)))
      exp.mm <- merge(exp.m,clinic[,c("barcode","age_at_initial_pathologic_diagnosis")],by.x="variable",by.y="barcode")
      colnames(exp.mm) <- c("barcode","GeneSymbol","value","age")
      exp.mm$GeneSymbol <- factor(exp.mm$GeneSymbol,levels=unique(exp.mm$GeneSymbol))
      #sapply(split(exp.mm[,c("value","age")],exp.mm$GeneSymbol),function(x){cor(log2(x[,1]+1),x[2],method="spearman")})
      exp.mm["age"][exp.mm["age"] >= 60] <- "Old"
      exp.mm["age"][exp.mm["age"] <= 50] <- "Young"
      
      exp.mm <- exp.mm[which(exp.mm$age %in% c("Young","Old")),]
     # size <- nrow(unique(exp.mm[which(exp.mm$age=="Young"),c("barcode","age")]))
     # random_old <- sample(unique(exp.mm[which(exp.mm$age=="Old"),c("barcode")]),size)
     # exp.mm <- exp.mm[which(exp.mm$barcode %in% c(as.character(unique(exp.mm[which(exp.mm$age=="Young"),"barcode"])),as.character(random_old))),]
      exp.mm$age <- factor(exp.mm$age,levels=unique(exp.mm$age))
      t.Pval <- sapply(split(exp.mm[,3:4],exp.mm$GeneSymbol),function(x){t.test(log2(value+1)~age,data = x)$p.value})
      t.Pval <- signif(p.adjust(t.Pval,method = "fdr"),digits = 2)
      t.Pval <- data.frame(t.Pval)
      t.Pval$GeneSymbol <- rownames(t.Pval)
      
      exp.mm$type <- rep(cancer,times=nrow(exp.mm))
      t.Pval$type <- rep(cancer,times=nrow(t.Pval))
      if(nrow(clock.exp.stage)==0){
        clock.exp.age <- exp.mm
        
      }else{
        clock.exp.age <- rbind(clock.exp.age,exp.mm)
      }
      
      if(nrow(t.PvalAll)==0){
        t.PvalAll <- t.Pval
      }else{t.PvalAll <- rbind(t.PvalAll,t.Pval)}
    }
    
  }
}
clock.exp.age_50 <- clock.exp.age
age_patient <- unique(clock.exp.age[,c("barcode","age","type")])
age_patient_count <- data.frame(t(sapply(split(age_patient[,"age"],age_patient$type),table)))
write.csv(age_patient_count,file="/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/age_patient_summary.csv",quote = F,row.names = T)
write.table(clock.exp.age,file="/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/clock.exp.age.txt",quote = F,row.names = F,sep="\t")

clock.exp.age.mean <- c()
for(i in unique(clock.exp.age$type)){
  for( j in unique(clock.exp.age$GeneSymbol)){
    sub <- clock.exp.age[which(clock.exp.age$type == i & clock.exp.age$GeneSymbol==j),]
    sub.exp <- sapply(split(sub[,"value"],sub$age),mean)
    sub.exp <- c(as.numeric(sub.exp),j,i)
    clock.exp.age.mean <- rbind(clock.exp.age.mean,sub.exp)
  }
}
clock.exp.age.mean <- data.frame(clock.exp.age.mean)
colnames(clock.exp.age.mean) <- c("young","old","GeneSymbol","type")
clock.exp.age.mean <- merge(clock.exp.age.mean,t.PvalAll,by=c("GeneSymbol","type"))

library(varhandle)
clock.exp.age.mean$young <- unfactor(clock.exp.age.mean$young)
clock.exp.age.mean$old <- unfactor(clock.exp.age.mean$old)
clock.exp.age.mean$class <- rep("None",times=nrow(clock.exp.age.mean))
clock.exp.age.mean["class"][log2(clock.exp.age.mean["young"]/ clock.exp.age.mean["old"]) >= log2(1.45) & clock.exp.age.mean["t.Pval"] <= 0.05] <- "Young"
clock.exp.age.mean["class"][log2(clock.exp.age.mean["young"]/ clock.exp.age.mean["old"]) <= -log2(1.45) & clock.exp.age.mean["t.Pval"] <= 0.05] <- "Old"
types.sig <- names(table(clock.exp.age.mean[which(clock.exp.age.mean$class != "None"),"type"])[table(clock.exp.age.mean[which(clock.exp.age.mean$class != "None"),"type"]) > 0])
clock.exp.age.mean <- clock.exp.age.mean[which(clock.exp.age.mean$type %in% types.sig),]
clock.exp.age.mean$type <- factor(clock.exp.age.mean$type,levels = types.sig)
clock.sig <- clock.exp.age.mean[which(clock.exp.age.mean$class != "None"),]
write.table(clock.sig,file="/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/clock.sig.txt",quote = F,row.names = F,sep="\t")
px_label <- names(table(clock.sig$type))[order(table(clock.sig$type))]
py_label <- names(table(clock.sig$GeneSymbol))[order(table(clock.sig$GeneSymbol))]
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
y_color <- rep("black",times=length(py_label))
y_color[py_label %in% core.circadian] <- "red"
write.table(clock.exp.age.mean,file="/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/clock.exp.age.mean.txt",quote = F,row.names = F,col.names = T,sep="\t")
pdf("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/age.significant.Pval_with.expression.trend_50_60_FC1.3.pdf",width=6,height = 10)#,width=4000,height=1500,res=400)
ggplot(clock.exp.age.mean,aes(x=type,y=GeneSymbol))+
  geom_tile(aes(fill=factor(class)),col="gray")+
  scale_fill_manual(limits=c("Young","Old","None"),values = c("blue","red","white"),name="")+
  scale_y_discrete(limit=py_label)+
  scale_x_discrete(limit=px_label)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=14,colour = y_color),
        axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),legend.position="bottom",legend.direction="horizontal")

dev.off()


LGG_sig <- clock.exp.age[which(clock.exp.age$type == "LGG" & clock.exp.age$GeneSymbol %in% c("RORC","NPAS2","ARNTL2")),]
LGG_sig_Pval <- clock.exp.age.mean[which(clock.exp.age.mean$type == "LGG" & clock.exp.age.mean$GeneSymbol %in% c("RORC","NPAS2","ARNTL2")),c(1,5)]
pdf("/extraspace/yye1/analysis/Circadian/expression/circadian.in.dif.age/LGG_RORC_NPAS2_ARNTL2_expression.pdf",width = 6,height = 2.5)
ggplot(LGG_sig, aes( x = age, y = log2(value+1), fill = age))+ stat_boxplot(geom = "errorbar", width = 0.2 )+
    geom_boxplot(width = 0.6,outlier.shape = NA)+
    scale_fill_manual( limits= c("Young","Old"),values = c("cyan","magenta"),guide=F)+
    scale_x_discrete(limits= c("Young","Old"),labels= c("Young","Old"))+
    geom_text(data = LGG_sig_Pval,aes(y = 11.5,x=1, label = paste("p = ",t.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
    scale_y_continuous(breaks=seq(0,12,by=4))+
    facet_wrap(~GeneSymbol,ncol=3)+
    ylab('Gene Expression')+xlab("")+
  theme(panel.background=element_rect(colour="black",fill="white",size=1.5),
        panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA),
        axis.text=element_text(size=16,color="black"),axis.title=element_text(size=18,color="black"),
        axis.ticks.length = unit(.25, "cm"),strip.text=element_text(size=16),
        strip.background = element_rect(fill=NA))
 dev.off()


#####
pdf(paste(cancer,".clock.genes.expression.across.different.age.pdf",sep=""),width = 16,height = 12)
P1 <- ggplot(exp.mm, aes( x = factor(age), y = log2(value+1), fill =  factor(age)))+ stat_boxplot(geom = "errorbar", width = 0.2 )+
  geom_boxplot(width = 0.6,outlier.shape = NA)+
  scale_fill_manual( limits= c("< 60",">= 60"),values = fill_colors,labels=c("<60",">=60"),guide=F)+
  scale_x_discrete(limits= c("< 60",">= 60"),labels=c("<60",">=60"))+
  geom_text(data = t.Pval,aes(y = 14,x=1.2, label = paste("p=",t.Pval,sep=""),fill=NA),size=4,hjust=0.5)+
  facet_wrap(~GeneSymbol,nrow=4)+
  #  geom_line( aes(group = cancer), linetype = 'longdash',position=position_jitter(w=0, h=0), size = 0.3)+
  #  geom_jitter(position = position_jitter(width = 0.05, height = 0), size = 0.6)+
  ylab('Gene Expression')+xlab("")+
  # ggtitle(paste0('Wilcox Pvalue: ',pvalue ))+
  theme(panel.background = element_rect(fill=NA,color="black",size=2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        # plot.title=element_text(size = 12),
        strip.text=element_text(size=16))
print(P1)
dev.off()


####correlation

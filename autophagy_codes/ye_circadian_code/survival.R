#!usr/bin/Rscript
setwd("~/Circadian/survival")
library(gdata)
library(doParallel)
library(doMC)
registerDoMC()
TCGAPath = '/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/'
filenameall = list.files(path =  TCGAPath, pattern = '_mRNA_each_exp_20160513')
circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
library(foreach)
library(methods)
library(survival)

folder <- "Survival_new"
if (!file.exists(folder)) { dir.create(folder) }
nn1 <- 40
foreach (kkk1 = c(1:length(filenameall))) %dopar% {
  filename <- filenameall[kkk1]  ## read apa result file one by one for each cancer type
  Cond <- gsub("_mRNA_each_exp_20160513", "", filename) ## get the cancer name
  exp_data <- read.table(file.path(TCGAPath, filename), header = T, row.names = 1)  ## using the first column as the row names
 
  exp_data <- exp_data[which(data.frame(do.call(rbind,strsplit(as.character(rownames(exp_data)),"\\|")))$X1 %in% circadian.genes$V1),]
  filter.names <- colnames(exp_data)[substr(colnames(exp_data),14,15) == "01"] ## get tumor samples 01
  exp_data <- exp_data[,filter.names]
  colnames(exp_data) <- substr(colnames(exp_data),1,12)
  colnames(exp_data) <- gsub("\\.", "-", colnames(exp_data))
  
  #########  proceed if tumor have at least nn1 samples exist
  ## tumor column selection
  if (ncol(exp_data) >= 40) {# at lest 20 samples required
    nameOfClinicFile <- paste( "/extraspace/TCGA/TCGA_clinical/",Cond, "_clinical_clean.txt", sep = "")
    if (file.exists(nameOfClinicFile)) {
      clinic_data <- read.table(nameOfClinicFile, sep = "\t", header = TRUE, comment.char = "", quote = "", fill = TRUE)
      clinic_ID <- as.vector(clinic_data[,1])
      commonsample <- intersect(colnames(exp_data), clinic_ID) ## get the sample with clinc info
      exp_tumor_W_clinicID <- exp_data[,match(commonsample, colnames(exp_data))]
      clinic_data <- clinic_data[match(commonsample, clinic_ID),]
      
      write.table(clinic_data, file = file.path(folder, paste(Cond, "-Match-Subtype-Info.txt", sep = "")), sep = "\t", quote=FALSE, row.names=FALSE)
      ###########  survival output
      Result <- matrix(NA, nrow(exp_tumor_W_clinicID), 8, byrow=TRUE)
      rownames(Result) <- rownames(exp_tumor_W_clinicID)
      colnames(Result) <- c("coef", "Exp(coef)", "Coxp", "KMp", "N", "low", "high", "FC")
      for (i1 in 1:nrow(Result)){
        
        mRNA_exp <- as.vector(t(exp_tumor_W_clinicID[i1,]))
        keeplink <- which(!is.na(clinic_data[,"os_status"]) & !is.na(clinic_data[,"os_days"]) & (clinic_data[,"os_days"]) >= 0 & !is.na(mRNA_exp))
        Time.dfs <- as.numeric(as.vector(clinic_data[keeplink,"os_days"]))
        cen.status <- ifelse(as.vector(clinic_data[keeplink,"os_status"]) == "Dead", 1,0)
        mRNA_exp_refined <-  mRNA_exp[keeplink]	
      #  keeplink_quantile <- c(1:length(mRNA_exp_refined))[mRNA_exp_refined >= quantile(mRNA_exp_refined,seq(0,1,0.2))[4] | mRNA_exp_refined <= quantile(mRNA_exp_refined,seq(0,1,0.2))[3]]
       # Time.dfs <- Time.dfs[keeplink_quantile]
       # cen.status <- cen.status[keeplink_quantile]
       # mRNA_exp_refined <- mRNA_exp_refined[keeplink_quantile]
        if(length(mRNA_exp_refined) >= nn1){
          Result[i1, "N"] <- length(mRNA_exp_refined)
          test.data1 <- list(time     = Time.dfs,
                             status   = cen.status,
                             group    = mRNA_exp_refined)
          tempresult<-try(model1 <- coxph(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude),silent=TRUE)
          if(!is(tempresult, "try-error")){
            Result[i1, c("coef", "Exp(coef)", "Coxp")] <- summary(model1)$coefficients[1,c("coef", "exp(coef)", "Pr(>|z|)" )]
            cutgroup <- ifelse(as.vector(mRNA_exp_refined) <= median(as.vector(mRNA_exp_refined)), "Low","High")
            if(length(unique(cutgroup)) >1){
              test.data1 <- list(time     = Time.dfs,
                                 status   = cen.status,
                                 group    = as.factor(cutgroup))
              model1 <- survdiff(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
              Result[i1, c("KMp")] <- 1-pchisq(model1$chisq, df=length(levels(factor(cutgroup)))-1)
              Result[i1, "low"] <- mean(mRNA_exp_refined[cutgroup == "Low"])
              Result[i1, "high"] <- mean(mRNA_exp_refined[cutgroup == "High"])
              Result[i1, "FC"] <- log2((mean(mRNA_exp_refined[cutgroup == "High"])+1)/(mean(mRNA_exp_refined[cutgroup == "Low"])+1))
              fit <- survfit(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
              if(Result[i1, c("KMp")] < 0.05){
                pdf(paste(folder,"/",Cond,"_",row.names(Result)[i1],".pdf",sep=""),width=6,height = 6)
                #pdf("22.pdf",width=6,height = 6)
                plot(fit,col=c("red","blue"),lty=1,lwd=2,mark.time=TRUE,main=paste("Kaplan-Meier Curves ",rownames(Result)[i1],sep=""),xlab = "Survival in days",cex.lab=1.5,cex.axis=1.2)
                legend("topright", attributes(as.factor(test.data1$group))$levels, col=c("red","blue"),lty=1,box.lwd=2,lwd=2)
                text(max(Time.dfs)/6,0.1,paste("p = ",signif(as.numeric(Result[i1, c("KMp")]),digits = 2),sep=""),cex=1.2)
                dev.off()
              
              }
            }
          }
        }######### end of if each marker has at least 20 no NA samples
      }######### end of each marker
      
      tmp <- as.numeric(as.vector(Result[, "Coxp"]))
      CoxFDRAdjustPvalue <- p.adjust(tmp, method="fdr")
      CoxOtherAdjustPvalue <- p.adjust(tmp, method="bonferroni")
      tmp <- as.numeric(as.vector(Result[, "KMp"]))
      KMpFDRAdjustPvalue <- p.adjust(tmp, method="fdr")
      KMpOtherAdjustPvalue <- p.adjust(tmp, method="bonferroni")
      Result <- cbind(Marker = rownames(Result), Result,  FDR = CoxFDRAdjustPvalue, Bonferroni = CoxOtherAdjustPvalue,
                      KMFDR = KMpFDRAdjustPvalue, KMBonferroni = KMpOtherAdjustPvalue)
      nameused <- paste(Cond, "-Survival.txt", sep = "")
      write.table(Result, file =file.path(folder,nameused), sep="\t", quote=FALSE, row.names=FALSE)
      keeplink <- which((as.numeric(Result[,"KMFDR"]) < 0.054 & !is.na(KMpFDRAdjustPvalue)))
      if(length(keeplink) > 0){
        nameused <- paste(Cond, "-Survival-Select.txt", sep = "")
        if(length(keeplink) == 1){
          
          write.table(data.frame(t(as.matrix(Result[keeplink,]))), file =file.path(folder,nameused), sep="\t", quote=FALSE, row.names=FALSE)
        }else{
          write.table(Result[keeplink,], file =file.path(folder,nameused), sep="\t", quote=FALSE, row.names=FALSE)
        }

        
              } ###########  select survival with FDR < 0.05 markers
      
    } ##########  end of survival
    
  }#########  end of proceed if tumor have at least nn1 samples exist
  
#  write.table("Done", file =file.path("Survival", paste("Cond", Cond, ".txt", sep = "")), row.names=F, col.names=F, quote=F)
  
}

folder <- "~/Circadian/survival/Survival_new/"
list.files(path =  folder, pattern = '-Subtype-Info.txt')
outfileall <- list.files(path =  folder, pattern = '-Survival-Select.txt')
outfileallnames <- gsub("-Survival-Select.txt","",outfileall)
for(i in 1:length(outfileall)){
  sur_select <- read.delim(paste(folder,"/",outfileall[i],sep=""), header=T)
  sur_select$survival <- ifelse(sur_select$Exp.coef. > 1,"High","Low")
  sur_select$tumor <- rep(outfileallnames[i],times=nrow(sur_select))
  sur_all.m <- sur_select[which(sur_select$KMp < 0.05),c("Marker","survival")]
  colnames(sur_all.m) <-c("Marker",outfileallnames[i])
  if(i==1){
    sur_all <- sur_select
    sur_all.mm  <- sur_all.m
  }else{
    sur_all <- rbind(sur_all,sur_select)
    sur_all.mm <- merge(sur_all.mm,sur_all.m,by="Marker",all=T)
  }
}
sur_all$geneSymbol <- data.frame(do.call(rbind, strsplit(as.character(sur_all$Marker),'\\|')))$X1
sur_all.mm$geneSymbol <- data.frame(do.call(rbind, strsplit(as.character(sur_all.mm$Marker),'\\|')))$X1
sur_all.p <- melt(sur_all.mm,id.vars="geneSymbol",measure.vars=colnames(sur_all.mm)[2:7])
colnames(sur_all.p) <- c("geneSymbol","tumor","survival")
write.table(sur_all.p,file="~/Circadian/survival/circadian.genes.survival.count.txt",quote = ,row.names = F,sep="\t")
#sur_all.p <- read.delim("~/Circadian/survival/circadian.genes.survival.count.txt",header=T)
sur_all.p <-  sur_all.p[!( sur_all.p$geneSymbol %in% c("FOXO1","FOXM1","FOXA2","FOXO3")),]

sur_all.p <- sur_all.p[!is.na(sur_all.p$survival),]
sur_all.p$geneSymbol <- factor(sur_all.p$geneSymbol,levels = unique(sur_all.p$geneSymbol))
sur_all.p$tumor <- factor(sur_all.p$tumor,levels=unique(sur_all.p$tumor))
x_label <- names(table(sur_all.p$tumor))[order(table(sur_all.p$tumor))]
y_label <- names(table(sur_all.p$geneSymbol))[order(table(sur_all.p$geneSymbol))]
core.circadian <- c("ARNTL2","ARNTL","CLOCK","NPAS2","NR1D1","NR1D2","RORA","RORB","RORC","PER1","PER2","PER3","CRY1","CRY2")
y_color <- rep("black",times=length(y_label))
y_color[y_label %in% core.circadian] <- "red"
miss <- data.frame(geneSymbol=rep(unique(sur_all.p$geneSymbol),times=length(unique(sur_all.p$tumor))),tumor=rep(unique(sur_all.p$tumor),each=length(unique(sur_all.p$geneSymbol))),survival=rep(0,times=length(unique(sur_all.p$geneSymbol))*length(unique(sur_all.p$tumor))))
miss <- merge(miss,sur_all.p,by=c("geneSymbol","tumor"),all.x=T) 
miss <- miss[(is.na(miss$survival.y)==T),1:3]
colnames(miss) <- colnames(sur_all.p)
sur_all.pp <- rbind(sur_all.p,miss)
pdf("~/Circadian/survival/circadian.genes.survival1.pdf",width=9,height = 10)#,width=4000,height=1500,res=400)
ggplot(sur_all.pp,aes(x=tumor,y=geneSymbol))+
  geom_tile(aes(fill=factor(survival)),col="lightgray")+
  scale_fill_manual(values = c("High"="red","Low" = "blue"),labels=c("H_Poor_Sur","L_Poor_Sur","None"),na.value="white",name="Survival")+
  scale_y_discrete(limit=y_label)+
  scale_x_discrete(limit=x_label)+
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


cir.exp.pval <- read.delim("~/Circadian/expression/cir.exp.pval.txt",header=T)
cir.exp.pval <- cir.exp.pval[!(cir.exp.pval$geneSymbol %in% c("FOXO1","FOXO3","FOXM1","FOXA2")),]
cir.exp.pval$geneSymbol <- factor(cir.exp.pval$geneSymbol,levels=unique(cir.exp.pval$geneSymbol))
cir.exp.pval$class <- rep(0,times=nrow(cir.exp.pval))
cir.exp.pval["class"][cir.exp.pval["Fold"] > 0] <- "Tumor-Up"
cir.exp.pval["class"][cir.exp.pval["Fold"] < 0] <- "Tumor-Down"
cir.exp.pval.m <- cir.exp.pval[,c(1,2,6)]
colnames(cir.exp.pval.m)<- c("geneSymbol","Type","Class")
sur_all.p.m <- sur_all.p
colnames(sur_all.p.m) <- c("geneSymbol","Type","Class")
sur_exp <- rbind(cir.exp.pval.m,sur_all.p.m)
dup <- apply(sur_exp[duplicated(sur_exp[,c("geneSymbol","Type")]),c("geneSymbol","Type")],1,function(x){paste(x[1],x[2],sep="\t")})
sur_exp.sp <- sur_exp[!(duplicated(sur_exp[,c("geneSymbol","Type")])),]#sur_exp[!( apply(sur_exp[,c("geneSymbol","Type")],1,function(x){paste(x[1],x[2],sep="\t")}) %in% dup),]
#dup.data <- sur_exp[( apply(sur_exp[,c("geneSymbol","Type")],1,function(x){paste(x[1],x[2],sep="\t")}) %in% dup),] 
dup.data <- dup.data[which(dup.data$Class %in% c("High","Low")),]
miss <- data.frame(geneSymbol=rep(unique(sur_exp$geneSymbol),times=length(unique(sur_exp$Type))),Type=rep(unique(sur_exp$Type),each=length(unique(sur_exp$geneSymbol))),Class=rep(0,times=length(unique(sur_exp$geneSymbol))*length(unique(sur_exp$Type))))
miss <- merge(miss,sur_exp.sp,by=c("geneSymbol","Type"),all.x=T) 
miss <- miss[(is.na(miss$Class.y)==T),1:3]
colnames(miss) <- colnames(sur_exp.sp)
sur_exp.spp <- rbind(sur_exp.sp,miss)
sur_exp.spp["Class"][sur_exp.spp["Class"] == 0] <- "None"

library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
x_label <- names(table(sur_exp.sp$Type))[order(table(sur_exp.sp$Type))]
y_label <- names(table(sur_exp.sp$geneSymbol))[order(table(sur_exp.sp$geneSymbol))]
dup.data$x = sapply(unfactor(dup.data$Type),function(x){grep(x,x_label)}) ## add the geom_polygon ID, convert the factor to numeric 
dup.data$y = sapply(unfactor(dup.data$geneSymbol),function(x){grep(x,y_label)})#as.numeric(dup.data$geneSymbol)
poly_x = c()
poly_y = c()
for ( i in 1: nrow(dup.data)){
  poly_x = c(poly_x, dup.data[i,"x"]-0.48, dup.data[i,'x'],dup.data[i,"x"],dup.data[i,"x"]-0.48 )
  poly_y = c(poly_y, dup.data[i,"y"]-0.48,dup.data[i,"y"]-0.48,dup.data[i,"y"]+0.48,dup.data[i,"y"]+0.48)
}
polygonID_2 <- data.frame( group = rep(seq(1:nrow(dup.data)), each = 4),poly_x, poly_y )


pdf("~/Circadian/survival/circadian.genes.survival_exp1.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
ggplot(sur_exp.spp,aes(x=Type,y=geneSymbol))+
  geom_tile(aes(fill=factor(Class)),col="white")+
  scale_fill_manual(limits=c("Tumor-Up","Tumor-Down","High","Low","None"),values = c("red","blue","magenta","cyan","lightgray"),labels=c("Tumor-Up","Tumor-Down","H_Poor_Sur","L_Poor_Sur","None"),na.value="lightgray",name="")+
  scale_y_discrete(limit=y_label)+
  scale_x_discrete(limit=x_label)+
  geom_polygon(data=polygonID_2,aes(x=poly_x,y=poly_y,group=group),fill="magenta")+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=16,colour = "black"),
        axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.position="bottom",legend.direction="horizontal")
 
dev.off()

as.numeric(sur_exp.sp$geneSymbol)













outfileall <- list.files(path =  "Survival_new/", pattern = '-Survival.txt')
outfileallnames <- gsub("-Survival.txt","",outfileall)
for(i in 1:length(outfileall)){
  sur_select <- read.delim(paste("Survival_new/",outfileall[i],sep=""), header=T)
  sur_select <- sur_select[which(sur_select$FDR   < 0.05),]
  sur_select$survival <- ifelse(sur_select$Exp.coef. > 1,"High","Low")
  sur_select$tumor <- rep(outfileallnames[i],times=nrow(sur_select))
  sur_all.m <- sur_select[,c("Marker","survival")]
  colnames(sur_all.m) <-c("Marker",outfileallnames[i])
  if(i==1){
    sur_all <- sur_select
    sur_all.mm  <- sur_all.m
  }else{
    sur_all <- rbind(sur_all,sur_select)
    sur_all.mm <- merge(sur_all.mm,sur_all.m,by="Marker",all=T)
  }
}
sur_all$geneSymbol <- data.frame(do.call(rbind, strsplit(as.character(sur_all$Marker),'\\|')))$X1
sur_all.mm$geneSymbol <- data.frame(do.call(rbind, strsplit(as.character(sur_all.mm$Marker),'\\|')))$X1
sur_all.p <- melt(sur_all.mm,id.vars="geneSymbol",measure.vars=colnames(sur_all.mm)[2:29])
colnames(sur_all.p) <- c("geneSymbol","tumor","survival")
#write.table(sur_all.p,file="~/Circadian/survival/circadian.genes_pvalue.survival.count.txt",quote = ,row.names = F,sep="\t")
#sur_all.p <- read.delim("~/Circadian/survival/circadian.genes_pvalue.survival.count.txt",header=T)
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
pdf("~/Circadian/survival/circadian.genes_pvalue.survival.pdf",width=10,height = 10)#,width=4000,height=1500,res=400)
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

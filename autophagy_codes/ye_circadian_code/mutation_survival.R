#### find circadian genes mutation
setwd("/home/yye1/Circadian/mutation")
library(ggplot2)
library(reshape2)
files.names <- list.files(path="/extraspace/TCGA/Mutation/",pattern=".txt")
#pair_tumor <- read.delim("/home/yye1/Circadian/expression/tumor.normal.sample.calculate_morethan_10pairs.txt",header=T)
#files.names <- files.names[gsub(".txt|mutation_","",files.names) %in% pair_tumor$Type]
files.Abs <- gsub(".txt","",files.names)
files.Abs <- gsub("mutation_","",files.Abs)
circadian.genes <- read.table("/home/yye1/Circadian/circadian.genes.txt",header=F)
circadian.genes <- circadian.genes$V1#data.frame(do.call(rbind, strsplit(as.character(circadian.genes$V1),'\\|')))$X1
sample.num <- c("Tumor.type","Number")
for(m in 1:length(files.names)){
  BLCA <- read.delim(paste("/extraspace/TCGA/Mutation/",files.names[m],sep=""),header=T)
  AllMutCount <- apply(BLCA[,2:ncol(BLCA)],2,sum)
  AllMutCountNames <- names(AllMutCount)[AllMutCount < 1000]
  AllMutCountNames <- AllMutCountNames[!is.na(AllMutCountNames)]
  BLCA <- BLCA[,c(colnames(BLCA)[1], AllMutCountNames)]
  BLCA <- BLCA[which(BLCA$X %in% core.circadian),] ##get circadian genes
  mutSample <- apply(BLCA[,2:ncol(BLCA)],2,sum)
  mutSample <- data.frame(mutSample)
  mutSample$barcode <- gsub("\\.","\\-",rownames(mutSample))
  if(nrow(mutSample[which(mutSample$mutSample >= 1),]) >= 10){
    nameOfClinicFile <- paste( "/extraspace/TCGA/TCGA_clinical/",files.Abs[m], "_clinical_clean.txt", sep = "")
    clinic_data <- read.table(nameOfClinicFile, sep = "\t", header = TRUE, comment.char = "", quote = "", fill = TRUE)
    clinic_ID <- as.vector(clinic_data[,1])
    clinic_data <- merge(mutSample,clinic_data,by="barcode")
    keeplink <- which(!is.na(clinic_data[,"os_status"]) & !is.na(clinic_data[,"os_days"]) & (clinic_data[,"os_days"]) >= 0 )
    Time.dfs <- as.numeric(as.vector(clinic_data[keeplink,"os_days"]))
    cen.status <- ifelse(as.vector(clinic_data[keeplink,"os_status"]) == "Dead", 1,0)
    mutexp <- clinic_data$mutSample[keeplink]
    test.data1 <- list(time     = Time.dfs,
                       status   = cen.status,
                       group    = mutexp)
    tempresult<-try(model1 <- coxph(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude),silent=TRUE)
    cutgroup <- ifelse(as.vector(mutexp) <= 0, "Non","Mutation")
     test.data1 <- list(time     = Time.dfs,
                       status   = cen.status,
                       group    = as.factor(cutgroup))
    model1 <- survdiff(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
    fit <- survfit(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(cutgroup)))-1)
    if(KMP <= 0.06){
      pdf(paste("~/Circadian/mutation/core circadian genes mutation_survival in",files.Abs[m],"_filter_ultramutation.pdf",sep=""),width=5,height = 5)
      plot(fit,col=c("red","blue"),lty=1,lwd=2,mark.time=TRUE,main=paste("KM Curves ","core circadian genes",sep=""),xlab = "Survival in days",cex.lab=1.5,cex.axis=1.2)
      legend("topright", attributes(as.factor(test.data1$group))$levels, col=c("red","blue"),lty=1,box.lwd=2,lwd=2)
      text(max(Time.dfs)/6,0.1,paste("p = ",signif(as.numeric(KMP),digits = 2),sep=""),cex=1.2)
      dev.off()
    }
    print(c(files.Abs[m],KMP))
  }
  
}
clinic_files <- list.files("/extraspace/TCGA/TCGA_clinical/",pattern= "_clinical_clean.txt")
commonfiles <- intersect(gsub("mutation_|\\.txt","",files.names),gsub("_clinical_clean.txt","",clinic_files))
clinic_dataAll <- data.frame()
for(m in commonfiles){
  BLCA <- read.delim(paste("/extraspace/TCGA/Mutation/",files.names[grep(m,files.names)],sep=""),header=T)
  AllMutCount <- apply(BLCA[,2:ncol(BLCA)],2,sum)
  AllMutCountNames <- names(AllMutCount)[AllMutCount < 1000]
  AllMutCountNames <- AllMutCountNames[!is.na(AllMutCountNames)]
  BLCA <- BLCA[which(BLCA$X %in% c("PER1","PER2","PER3")),] ##get circadian genes
  mutSample <- apply(BLCA[,2:ncol(BLCA)],2,sum)
  mutSample <- data.frame(mutSample)
  mutSample$barcode <- gsub("\\.","\\-",rownames(mutSample))
  if(nrow(mutSample[which(mutSample$mutSample >= 1),]) >=5){
    nameOfClinicFile <- paste( "/extraspace/TCGA/TCGA_clinical/",m, "_clinical_clean.txt", sep = "")
    clinic_data <- read.table(nameOfClinicFile, sep = "\t", header = TRUE, comment.char = "", quote = "", fill = TRUE)
    clinic_ID <- as.vector(clinic_data[,1])
    clinic_data <- merge(mutSample,clinic_data,by="barcode")
    clinic_data <- clinic_data[,c("barcode","mutSample","os_days","os_status")]
    if(nrow(clinic_dataAll)==0){
      clinic_dataAll <- clinic_data
    }else{
      clinic_dataAll <- rbind(clinic_dataAll,clinic_data)
    }

    keeplink <- which(!is.na(clinic_dataAll[,"os_status"]) & !is.na(clinic_dataAll[,"os_days"]) & (clinic_dataAll[,"os_days"]) >= 0 )
    Time.dfs <- as.numeric(as.vector(clinic_dataAll[keeplink,"os_days"]))
    cen.status <- ifelse(as.vector(clinic_dataAll[keeplink,"os_status"]) == "Dead", 1,0)
    mutexp <- clinic_dataAll$mutSample[keeplink]
    test.data1 <- list(time     = Time.dfs,
                       status   = cen.status,
                       group    = mutexp)
    tempresult<-try(model1 <- coxph(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude),silent=TRUE)
    cutgroup <- ifelse(as.vector(mutexp) <= 0, "Non","Mutation")
    test.data1 <- list(time     = Time.dfs,
                       status   = cen.status,
                       group    = as.factor(cutgroup))
    model1 <- survdiff(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
    fit <- survfit(Surv(time, status) ~ group, data=test.data1, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(cutgroup)))-1)
    
    if(KMP <= 0.05){
      pdf("~/Circadian/mutation/PER mutation_survival across all cancer samples_filter_ultramutation.pdf",width=5,height = 5)
      plot(fit,col=c("red","blue"),lty=1,lwd=2,mark.time=TRUE,main=paste("Kaplan-Meier Curves ","PER mutation",sep=""),xlab = "Survival in days",cex.lab=1.5,cex.axis=1.2)
      legend("topright", attributes(as.factor(test.data1$group))$levels, col=c("red","blue"),lty=1,box.lwd=2,lwd=2)
      text(max(Time.dfs)/6,0.1,paste("p = ",signif(as.numeric(KMP),digits = 2),sep=""),cex=1.2)
      dev.off()
    }
    print(c(files.Abs[m],KMP))
  }
}
###

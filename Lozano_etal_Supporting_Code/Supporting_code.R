#This code is a supporting resource for Lozano, Chaudhuri, Nene et al., Nature Medicine 2021
#Scripts for training and validation of composite models, evaluating freedom from severe toxicity, and generating related figures are included.
#Contact: Aaron M. Newman (amnewman@stanford.edu)
#Last modified: 11-02-2021
#Version 1.0

#   Requirements
#1. Ensure the following R packages are installed: pROC, survminer, ggplot2
#2. Place 'model_input.txt' and this R script in the same directory, or specify the full path to  'model_input.txt' in the INPUT variable below
#3. From the command line, run: "Rscript Supporting_code.R"

suppressMessages(library(survival))
suppressMessages(library(pROC))
suppressMessages(library(survminer))
suppressMessages(library(ggplot2))

options(scipen = 999)

INPUT <- "model_input.txt" #Path to model_input.txt

actCD4mem_ceiling <- 0.1 #All activated CD4 T cell values identified as outliers by the ROUT test in Prism with 10% FDR have frequency >0.1, regardless of training subset in the analyses below (Methods)

x <- read.table(INPUT,sep="\t", header=T, row.names=1)

#================Functions===================

#Build and run composite model on training and validation sets
makeCompositeModel <- function(input, index, actCD4mem_ceiling, LOOCV){
    
    if(!LOOCV){
        ##### replace outliers with max value of non-outliers
        AMCD4 <- input$Activated_CD4_memory_T_CSx
        mx <- max((AMCD4[index])[which(AMCD4[index]<actCD4mem_ceiling)])
        AMCD4[which(AMCD4 > mx)] <- mx
        input$AMCD4 <- AMCD4
        #####

        #Train logistic regression model ('composite model') on activated CD4 memory T cell abundance and TCR Shannon entropy using patient subset defined by 'index'
        rf <- glm(formula=as.factor(Severe_irAE)~AMCD4 + TCR_clonotype_diversity_ShannonEntropy,data=input[index,],family="binomial")
        #Run composite model on all pts
        out <- predict(rf, input,type=c("response"))
        modelPred <- out
        input$CompositeModel <- modelPred #Continuous model results
        
    }else{
        
        #Leave-one-out cross-validatoin (LOOCV)
        input2 <- input
        
        for(i in 1:nrow(input)){
            
            ##### replace outliers with max value of non-outliers
            AMCD4 <- input$Activated_CD4_memory_T_CSx
            mx <- max((AMCD4[-i])[which(AMCD4[-i]<actCD4mem_ceiling)])
            AMCD4[which(AMCD4 > mx)] <- mx
            input2$AMCD4 <- AMCD4
            #####

            #Train logistic regression model ('composite model') on activated CD4 memory T cell abundance and TCR Shannon entropy using patient subset using LOOCV
            rf <- glm(formula=as.factor(Severe_irAE)~AMCD4 + TCR_clonotype_diversity_ShannonEntropy,data=input2[-i,],family="binomial")
            #Run composite model on all pts
            out <- predict(rf, input2,type=c("response"))
            modelPred <- out
            input$CompositeModelLOOCV[i] <- modelPred[i] #Continuous model results
            
            #Optimize split using Youden's J statistic
            suppressMessages(g <- pROC::roc(as.factor(input2$Severe_irAE[-i]) ~ modelPred[-i], data=input2))
            cutp <- g$thresholds[which.max(g$sensitivities+g$specificities)]
            
            if(modelPred[i] > cutp)  input$BinaryModelLOOCV[i] <- 1
            else input$BinaryModelLOOCV[i] <- 0
        }
    }
    
    input #return input
}

#Calculate AUC of composite model in patient subset 'index'
calculateAUC <- function(input, index, LOOCV){
    
    if(!LOOCV){
        suppressMessages(g <- pROC::roc(as.factor(input$Severe_irAE[index]) ~ CompositeModel[index], data=input))
    }else{
        suppressMessages(g <- pROC::roc(as.factor(input$Severe_irAE[index]) ~ CompositeModelLOOCV[index], data=input))
    }
    g #return AUC model

}

#Optimize split using Youden's J statistic and generate binary composite model
binarizeModel <- function(input, g){

    cutp <- g$thresholds[which.max(g$sensitivities+g$specificities)]
    q <- rep(0,nrow(input))
    q[which(input$CompositeModel>=cutp)] <- 1
    input$BinaryModel <- q
    input #return input
}

#Get hazard ratio of binary composite model
getCox <- function(input, index, LOOCV){
    
    if(!LOOCV){
        mod <- suppressWarnings(summary(coxph(Surv(Time_to_severe_irAE_months, Severe_irAE) ~ BinaryModel, data = input[index,])))
    }else{
        mod <- suppressWarnings(summary(coxph(Surv(Time_to_severe_irAE_months, Severe_irAE) ~ BinaryModelLOOCV, data = input[index,])))
    }
    c(mod$coefficients[1,2],mod$sctest[3])
}

makeKMplot <- function(input, index, label, LOOCV){
    
    pdf(file=paste(label,".pdf", sep=""), height=5, width = 4.5, onefile=FALSE)
    if(!LOOCV){
        fit <- survfit(Surv(Time_to_severe_irAE_months, Severe_irAE) ~ BinaryModel, data = input[index,])
    }else{
        fit <- survfit(Surv(Time_to_severe_irAE_months, Severe_irAE) ~ BinaryModelLOOCV, data = input[index,])
    }
    print(ggsurvplot(fit, data = input[index,], risk.table = TRUE, pval = TRUE, xlim=c(0,4),break.time.by = 1, risk.table.y.text = FALSE, palette=c("#1B76FF","#BD111D"), xlab="Months after treatment initiation", ylab="Fraction free from severe irAE", censor=F))
    q <- dev.off()
    
}

############################################################################

#Figure 4d (left)

#Train on Bulk cohort 1 (all patients)
train_lab <- "Bulk cohort 1"

train <- which(x$Cohort==train_lab)
x <- makeCompositeModel(x, train, actCD4mem_ceiling, FALSE)
AUCtrain <- calculateAUC(x, train, FALSE)
x <- binarizeModel(x, AUCtrain)
#Write composite model to disk
write.table(data.frame(PID=rownames(x),CompositeModel=x$CompositeModel,BinaryModel=x$BinaryModel),paste("CompositeModel_Train_",gsub(" ","",train_lab),".txt", sep=""),sep="\t",quote=F)

#Test on Bulk cohort 2 (all patients)
test_lab <- "Bulk cohort 2"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab), sep="")
test <- which(x$Cohort==test_lab)

#Calculate AUC on test cohort (all patients)
AUCtest <- calculateAUC(x, test, FALSE)
cat(paste("Figure 4d (left) | Train:",train_lab,", Test: ", test_lab,", AUC: ",round(AUCtest$auc,3), "\n", sep=""))

#Test on Bulk cohort 2 (combination therapy patients)
test_lab <- "Bulk cohort 2"
test_lab2 <- "Combination"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab), gsub(" ","",test_lab2), sep="")
test <- which(x$Cohort==test_lab & x$Therapy==test_lab2)

#Calculate AUC on test cohort (combination therapy patients)
AUCtest2 <- calculateAUC(x, test, FALSE)
cat(paste("Figure 4d (left) | Train:",train_lab,", Test: ", test_lab," ", test_lab2, ", AUC: ",round(AUCtest2$auc,3),"\n", sep=""))

#Test on Bulk cohort 2 (PD1 therapy patients)
test_lab <- "Bulk cohort 2"
test_lab2 <- "Anti-PD1"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab), gsub(" ","",test_lab2), sep="")
test <- which(x$Cohort==test_lab & x$Therapy==test_lab2)

#Calculate AUC on test cohort (PD1 therapy patients)
AUCtest3 <- calculateAUC(x, test, FALSE)
cat(paste("Figure 4d (left) | Train:",train_lab,", Test: ", test_lab," ", test_lab2, ", AUC: ",round(AUCtest3$auc,3),"\n", sep=""))

#Combined AUC plot
pdf("Figure4d_left.pdf",height=4,width=4)
plot(AUCtest, las=1, col="#FFBF19",lwd=2.5, cex.lab=1.2, cex.axis=1.2)
plot(AUCtest3, col="#BDBDBD",lwd=2.5, add=T)
plot(AUCtest2, col="#1F4DE6",lwd=2.5, add=T)
legend("bottomright", legend=c(paste("Both therapies\nAUC = ",round(AUCtest$auc,2),"\n",sep=""), paste("Combination\nAUC = ",round(AUCtest2$auc,2),"\n",sep=""),paste("PD1\nAUC = ",round(AUCtest3$auc,2),"\n",sep="")),
col=c("#FFBF19", "#1F4DE6", "#BDBDBD"), lty=1, cex=0.8, bty = "n")
q <- dev.off()

######################################

#Figure 4d (right)

#Train on Combination therapy patients
train_lab <- "Combination"

train <- which(x$Therapy==train_lab)
x <- makeCompositeModel(x, train, actCD4mem_ceiling, FALSE)
AUCtrain <- calculateAUC(x, train, FALSE)
x <- binarizeModel(x, AUCtrain)
#Write composite model to disk
write.table(data.frame(PID=rownames(x),CompositeModel=x$CompositeModel,BinaryModel=x$BinaryModel),paste("CompositeModel_Train_",gsub(" ","",train_lab),".txt", sep=""),sep="\t",quote=F)

#Test on PD1 monotherapy patients
test_lab <- "Anti-PD1"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab), sep="")
test <- which(x$Therapy==test_lab)

#Calculate AUC on test cohort (all patients)
AUCtest <- calculateAUC(x, test, FALSE)
cat(paste("Figure 4d (right) | Train:",train_lab,", Test: ", test_lab,", AUC: ",round(AUCtest$auc,3), "\n", sep=""))

#Train on PD1 therapy patients
train_lab <- "Anti-PD1"

train <- which(x$Therapy==train_lab)
x <- makeCompositeModel(x, train, actCD4mem_ceiling, FALSE)
AUCtrain <- calculateAUC(x, train, FALSE)
x <- binarizeModel(x, AUCtrain)
#Write composite model to disk
write.table(data.frame(PID=rownames(x),CompositeModel=x$CompositeModel,BinaryModel=x$BinaryModel),paste("CompositeModel_Train_",gsub(" ","",train_lab),".txt", sep=""),sep="\t",quote=F)

#Test on Combination therapy patients
test_lab <- "Combination"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab), sep="")
test <- which(x$Therapy==test_lab)

#Calculate AUC on test cohort (all patients)
AUCtest2 <- calculateAUC(x, test, FALSE)
cat(paste("Figure 4d (right) | Train:",train_lab,", Test: ", test_lab,", AUC: ",round(AUCtest2$auc,3), "\n", sep=""))

#Combined AUC plot
pdf("Figure4d_right.pdf",height=4,width=4)
plot(AUCtest, las=1, col="#BDBDBD",lwd=2.5, cex.lab=1.2, cex.axis=1.2)
plot(AUCtest2, col="#1F4DE6",lwd=2.5, add=T)
legend("bottomright", legend=c(paste("Train on PD1,\ntest on combination\nCombination AUC = ",round(AUCtest2$auc,2),"\n",sep=""), paste("Train on combination,\ntest on PD1\nPD1 AUC = ",round(AUCtest$auc,3),"\n",sep="")),
col=c("#1F4DE6","#BDBDBD"), lty=1, cex=0.8, bty = "n")
q <- dev.off()

############################################################################

#Extended Data Figure 8

#Freedom from severe irAE analysis (remove two patients [YUMEDIC, YUDIME] that changed therapy due to progression within 3 months of ICI treatment initiation)

#Train on Bulk cohort 1
train_lab <- "Bulk cohort 1"

train <- which(x$Cohort==train_lab & !(rownames(x) %in% c("YUMEDIC", "YUDIME")))
x <- makeCompositeModel(x, train, actCD4mem_ceiling, FALSE)
AUCtrain <- calculateAUC(x, train, FALSE)
x <- binarizeModel(x, AUCtrain)
#Write composite model to disk
write.table(data.frame(PID=rownames(x),CompositeModel=x$CompositeModel,BinaryModel=x$BinaryModel),paste("CompositeModel_RmEarlyProgressors_Train_",gsub(" ","",train_lab),".txt", sep=""),sep="\t",quote=F)

#Test on Bulk cohort 2 (all patients)
test_lab <- "Bulk cohort 2"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab),"_RmEarlyProgressors", sep="")
test <- which(x$Cohort==test_lab)

#Calculate AUC and Hazard Ratio on test cohort
AUCtest <- calculateAUC(x, test, FALSE)
HRtest <- getCox(x, test, FALSE)
cat(paste("Extended Data Figure 8a | Train:",train_lab,", Test: ", test_lab,", HR: ",round(HRtest[1],3), ", Logrank p-value: ",round(HRtest[2],4), "\n", sep=""))
makeKMplot(x, test, "Extended_Data_Figure8a", FALSE)

#Test on Bulk cohort 2 (combination therapy patients)
test_lab <- "Bulk cohort 2"
test_lab2 <- "Combination"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab),test_lab2,"_RmEarlyProgressors", sep="")
test <- which(x$Cohort==test_lab & x$Therapy==test_lab2)

#Calculate AUC on test cohort (combination therapy patients)
AUCtest <- calculateAUC(x, test, FALSE)
HRtest <- getCox(x, test, FALSE)
cat(paste("Extended Data Figure 8b | Train:",train_lab,", Test: ", test_lab, " ",test_lab2,", HR: ",round(HRtest[1],3), ", Logrank p-value: ",round(HRtest[2],4), "\n", sep=""))
makeKMplot(x, test, "Extended_Data_Figure8b", FALSE)

#Test on Bulk cohort 2 (PD1 therapy patients)
test_lab <- "Bulk cohort 2"
test_lab2 <- "Anti-PD1"
label2 <- paste("Train_",gsub(" ","",train_lab),"_Test_", gsub(" ","",test_lab),test_lab2,"_RmEarlyProgressors", sep="")
test <- which(x$Cohort==test_lab & x$Therapy==test_lab2)

#Calculate AUC on test cohort (PD1 therapy patients)
AUCtest <- calculateAUC(x, test, FALSE)
HRtest <- getCox(x, test, FALSE)
cat(paste("Extended Data Figure 8c | Train:",train_lab,", Test: ", test_lab, " ", test_lab2,", HR: ",round(HRtest[1],3), ", Logrank p-value: ",round(HRtest[2],4), "\n", sep=""))
makeKMplot(x, test, "Extended_Data_Figure8c", FALSE)

############################################################################

#Extended Data Figure 7a

#Train on all patients using LOOCV
x <- makeCompositeModel(x, "", actCD4mem_ceiling, TRUE)

#Write LOOCV composite model to disk
write.table(data.frame(PID=rownames(x),CompositeModelLOOCV=x$CompositeModelLOOCV,BinaryModelLOOCV=x$BinaryModelLOOCV),paste("CompositeModelLOOCV_Train_",gsub(" ","",train_lab),".txt", sep=""),sep="\t",quote=F)

#Calculate AUC on test cohort (all patients)
AUCtest <- calculateAUC(x, 1:nrow(x), TRUE)
cat(paste("Extended Data Figure 7a | Train:LOOCV, Test: All, AUC: ",round(AUCtest$auc,3), "\n", sep=""))

#Test on Combination therapy patients
test_lab <- "Combination"

label2 <- paste("Train_LOOCV_Test_", gsub(" ","",test_lab), sep="")
test <- which(x$Therapy==test_lab)

#Calculate AUC on test cohort (combination therapy patients)
AUCtest2 <- calculateAUC(x, test, TRUE)
cat(paste("Extended Data Figure 7a | Train:LOOCV, Test: ", test_lab, ", AUC: ",round(AUCtest2$auc,3),"\n", sep=""))

#Test on PD1 monotherapy patients
test_lab <- "Anti-PD1"

label2 <- paste("Train_LOOCV_Test_", gsub(" ","",test_lab), sep="")
test <- which(x$Therapy==test_lab)

#Calculate AUC on test cohort (combination therapy patients)
AUCtest3 <- calculateAUC(x, test, TRUE)
cat(paste("Extended Data Figure 7a | Train:LOOCV, Test: ", test_lab, ", AUC: ",round(AUCtest3$auc,3),"\n", sep=""))

#Combined AUC plot
pdf("Extended_Data_Figure7a.pdf",height=4,width=4)
plot(AUCtest, las=1, col="#FFBF19",lwd=2.5, cex.lab=1.2, cex.axis=1.2)
plot(AUCtest3, col="#BDBDBD",lwd=2.5, add=T)
plot(AUCtest2, col="#1F4DE6",lwd=2.5, add=T)
legend("bottomright", legend=c(paste("Both therapies\nAUC = ",round(AUCtest$auc,2),"\n",sep=""), paste("Combination\nAUC = ",round(AUCtest2$auc,2),"\n",sep=""),paste("PD1\nAUC = ",round(AUCtest3$auc,2),"\n",sep="")),
col=c("#FFBF19", "#1F4DE6", "#BDBDBD"), lty=1, cex=0.8, bty = "n")
q <- dev.off()

############################################################################

#Figure 5a

#Freedom from severe irAE analysis in combination therapy patients via LOOCV (remove two patients [YUMEDIC, YUDIME] that changed therapy due to progression within 3 months of ICI treatment initiation)

#Train on all patients using LOOCV
x2 <- makeCompositeModel(x[which(!(rownames(x) %in% c("YUMEDIC", "YUDIME"))),], "", actCD4mem_ceiling, TRUE)

#Write LOOCV composite model to disk
write.table(data.frame(PID=rownames(x2),CompositeModelLOOCV=x2$CompositeModelLOOCV,BinaryModelLOOCV=x2$BinaryModelLOOCV),paste("CompositeModelLOOCV_RmEarlyProgressors_Train_",gsub(" ","",train_lab),".txt", sep=""),sep="\t",quote=F)

#Combination therapy patients
test_lab <- "Combination"
label2 <- paste("Train_LOOCV_Test_", gsub(" ","",test_lab),"_RmEarlyProgressors", sep="")
test <- which(x2$Therapy==test_lab)

AUCtest <- calculateAUC(x2, test, TRUE)
HRtest <- getCox(x2, test, TRUE)
cat(paste("Figure 5a | Train:LOOCV, Test: ", test_lab,", HR: ",round(HRtest[1],3), ", Logrank p-value: ",round(HRtest[2],4), "\n", sep=""))
makeKMplot(x2, test, "Figure5a", TRUE)

#PD1 therapy patients
test_lab <- "Anti-PD1"
label2 <- paste("Train_LOOCV_Test_", gsub(" ","",test_lab),"_RmEarlyProgressors", sep="")
test <- which(x2$Therapy==test_lab)

AUCtest <- calculateAUC(x2, test, TRUE)
HRtest <- getCox(x2, test, TRUE)
cat(paste("Main Text | Train:LOOCV, Test: ", test_lab,", HR: ",round(HRtest[1],3), ", Logrank p-value: ",round(HRtest[2],4), "\n", sep=""))

############################################################################

#Figure 4e (composite model LOOCV by grade)

x <- makeCompositeModel(x, "", actCD4mem_ceiling, TRUE)
d <- data.frame(model=x$CompositeModelLOOCV, grade=as.factor(x$irAE_grade), therapy=x$Therapy)
d$grade[which(d$grade==0)] <- "1"
d$grade[which(d$grade==1)] <- "1"
p <- ggplot(d, aes(x=grade, y=model, fill=grade)) +
  geom_boxplot(outlier.shape=NA, lwd=.75)
p <- p + scale_fill_manual(values=c("#E4EEFF",  "#ACCFFF","#DDA5A4", "#FF9381"), guide="none") + theme_classic()
p <- p + geom_point(pch=16, size=4,position=position_jitterdodge(0.5),  aes(color=grade))+ scale_color_manual(values=c("#B6C4CF",  "#68ABFF","#C56769", "#FC5C4A")) + theme(text = element_text(size=20), axis.text.x=element_text(colour="black", size=20),axis.text.y=element_text(colour="black", size=20),legend.position="bottom", axis.ticks.length = unit(.25, "cm"))+ scale_x_discrete(name="irAE grade") +
  scale_y_continuous(name="Composite model score", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1))
pdf(file="Figure4e.pdf", height=5.5,width=5)
p
q <- dev.off()
cat(paste("Figure 4e | Kruskal-Wallis p-value: ",round(kruskal.test(d$model~d$grade,data=d)$p.value,6),"\n",sep=""))

############################################################################

#Extended Data Figure 7e (composite model LOOCV by grade for combination therapy only)

x <- makeCompositeModel(x, "", actCD4mem_ceiling, TRUE)
x2 <- x[which(x$Therapy=="Combination"),]
d <- data.frame(model=x2$CompositeModelLOOCV, grade=as.factor(x2$irAE_grade), therapy=x2$Therapy)
d$grade[which(d$grade==0)] <- "1"
d$grade[which(d$grade==1)] <- "1"
p <- ggplot(d, aes(x=grade, y=model, fill=grade)) +
  geom_boxplot(outlier.shape=NA, lwd=.75)
p <- p + scale_fill_manual(values=c("#E4EEFF",  "#ACCFFF","#DDA5A4", "#FF9381"), guide="none") + theme_classic()
p <- p + geom_point(pch=16, size=4,position=position_jitterdodge(0.5),  aes(color=grade))+ scale_color_manual(values=c("#B6C4CF",  "#68ABFF","#C56769", "#FC5C4A")) + theme(text = element_text(size=20), axis.text.x=element_text(colour="black", size=20),axis.text.y=element_text(colour="black", size=20),legend.position="bottom", axis.ticks.length = unit(.25, "cm"))+ scale_x_discrete(name="irAE grade") +
  scale_y_continuous(name="Composite model score", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1))
pdf(file="Extended_Data_Figure7e.pdf", height=5.5,width=5)
p
q <- dev.off()
cat(paste("Extended Data Figure 7e | Kruskal-Wallis p-value: ",round(kruskal.test(d$model~d$grade,data=d)$p.value,5),"\n",sep=""))

############################################################################

#Extended Data Figure 7g (composite model LOOCV vs. No. symptomatic irAEs)

x <- makeCompositeModel(x, "", actCD4mem_ceiling, TRUE)

d <- data.frame(model=x$CompositeModelLOOCV, gr2plus=as.factor(x$No_grade_2plus_irAE), therapy=x$Therapy)
p <- ggplot(d, aes(x=gr2plus, y=model, fill=gr2plus)) +
  geom_boxplot(outlier.shape=NA, lwd=.75)
p <- p + scale_fill_manual(values=c("#E4EEFF",  "#9BACFA","#DDA5A4", "#FF9381"), guide="none") + theme_classic()
p <- p + geom_point(pch=16, size=4,position=position_jitterdodge(0.5),  aes(color=gr2plus))+ scale_color_manual(values=c("#B6C4CF",  "#3557F5","#C56769", "#FC5C4A")) + theme(text = element_text(size=20), axis.text.x=element_text(colour="black", size=20),axis.text.y=element_text(colour="black", size=20),legend.position="bottom", axis.ticks.length = unit(.25, "cm"))+ scale_x_discrete(name="No. of symptomatic irAEs") +
  scale_y_continuous(name="Composite model score", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1))
pdf(file="Extended_Data_Figure7g.pdf", height=5.5,width=5)
p
q <- dev.off()
cat(paste("Extended Data Figure 7g | Kruskal-Wallis p-value: ",round(kruskal.test(d$model~d$gr2plus,data=d)$p.value,5),"\n",sep=""))

############################################################################

#Extended Data Figure 7h (composite model LOOCV vs. No. organ system toxicities)

x <- makeCompositeModel(x,"", actCD4mem_ceiling, TRUE)

d <- data.frame(model=x$CompositeModelLOOCV, No_Organs=as.factor(x$No_organ_systems_irAE), therapy=x$Therapy)
d$No_Organs[which(d$No_Organs==4)] <- "3"
d$No_Organs[which(d$No_Organs==3)] <- "3"
p <- ggplot(d, aes(x=No_Organs, y=model, fill=No_Organs)) +
  geom_boxplot(outlier.shape=NA, lwd=.75)
p <- p + scale_fill_manual(values=c("#E4EEFF",  "#9BACFA","#DDA5A4", "#FF9381"), guide="none") + theme_classic()
p <- p + geom_point(pch=16, size=4,position=position_jitterdodge(0.5),  aes(color=No_Organs))+ scale_color_manual(values=c("#B6C4CF",  "#3557F5","#C56769", "#FC5C4A")) + theme(text = element_text(size=20), axis.text.x=element_text(colour="black", size=20),axis.text.y=element_text(colour="black", size=20),legend.position="bottom", axis.ticks.length = unit(.25, "cm"))+ scale_x_discrete(name="No. of organ system toxicities") +
  scale_y_continuous(name="Composite model score", limits=c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1))
pdf(file="Extended_Data_Figure7h.pdf", height=5.5,width=5)
p
q <- dev.off()
cat(paste("Extended Data Figure 7h | Kruskal-Wallis p-value: ",round(kruskal.test(d$model~d$No_Organs,data=d)$p.value,5),"\n",sep=""))

############################################################################

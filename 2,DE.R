library(ggplot2)
library(gridExtra)
library(ggpubr)

setwd("/home/wang/Documents/2024/2024_10/HE4_COPD/Data/")

symbole_ENTREZ <- read.table("hg19_symbol_entrez.txt", header = F, sep = "\t", stringsAsFactors = F)

need_symbole_ENTREZ <- symbole_ENTREZ[grep("WFDC",symbole_ENTREZ[,1]),][-1,]


GSE239897_exp <- read.table("GSE239897_TissueCoreProject_fpkm.txt", header = T, sep = ",", stringsAsFactors = F)
need_index <-  grep("WFDC",GSE239897_exp[,1])[-8]



need_GSE239897_exp <- GSE239897_exp[need_index, -1]
rownames(need_GSE239897_exp) <- GSE239897_exp[need_index,1]

GSE239897_info <- rep("COPD", 82)
GSE239897_info[grep("CTRL", colnames(need_GSE239897_exp))] <- "Control"


count <<- 1
DE_stat <- t(apply(
  as.matrix(need_GSE239897_exp),
  1,
  function(x){
    cat(count, "\n")
    count <<- count +1
    COPD_exp = x[which(GSE239897_info == "COPD")]
    nc_exp = x[which(GSE239897_info == "Control")]
    
    if(mean(COPD_exp) == 0 | mean(nc_exp) == 0) {
      c(mean(nc_exp), mean(COPD_exp), 1, 1)
    } else{
      FC = mean(COPD_exp) /  mean(nc_exp)
      
      p = wilcox.test(COPD_exp, nc_exp)$p.value
      
      return(c(mean(nc_exp), mean(COPD_exp), FC, p))
    }
    
  }
))

write.table(DE_stat, "WFDC_DE.txt", row.names = T, col.names = F, sep = "\t", quote = F)

###
library(ggpubr)
plot_data <- data.frame(
  
  WFDC2 = unlist(need_GSE239897_exp["WFDC2", ]),
  WFDC3 = unlist(need_GSE239897_exp["WFDC3", ]), 
  WFDC5 = unlist(need_GSE239897_exp["WFDC5", ]),
  WFDC10A = unlist(need_GSE239897_exp["WFDC10A", ]),
  WFDC12 = unlist(need_GSE239897_exp["WFDC12", ]),
  WFDC13 = unlist(need_GSE239897_exp["WFDC13", ]),
  type = GSE239897_info
)

mycompare_list <- list(c("COPD", "Control"))
p6 <- ggplot(data = plot_data, aes(type, WFDC13)) +
  geom_violin(aes(fill = type), colour = "white",show.legend = F) +
  scale_fill_manual(values = c("#cabbe9", "#e0c45c")) +
  geom_boxplot(colour = "gray", width = 0.1) +
  ylab("Expression of WFDC13") +
  xlab("") +
  #scale_x_discrete(labels =  c("Metastasis tumor", "Primary tumor")) +
  stat_compare_means(comparisons = mycompare_list, method = "wilcox.test") +
  theme_classic()






library(ComplexHeatmap)
library(circlize)


mat <- as.matrix(need_GSE239897_exp[, c(which(GSE239897_info == "Control"), which(GSE239897_info == "COPD"))])
mat <- t(scale(t(mat)))
ha_column = HeatmapAnnotation(final_data = data.frame(type = rep(c("Control", "COPD"), c(43, 39))), 
                              col = list(type = c("COPD" = "#FF6666", "Control" = "#99CC66")))
Heatmap(mat,top_annotation = ha_column,name="legend", cluster_columns = F, 
        show_column_names = F, 
        show_row_names = T,
        row_names_gp = gpar(fontsize = 9),
        col = colorRamp2(c(-4, 0, 4), c("#003399", "#CCCCCC", "#CC0033")))


plot_data <- data.frame(
  GSVA = unlist(need_GSE239897_exp["WFDC2", -c(59:82)]),
  #GSVA = unlist(GSE162955_exp["TUBGCP3",]),
  type = GSE239897_info[-c(59:82)]
)

mycompare_list <- list(c("Control", "COPD"))
p1 <- ggplot(data = plot_data, aes(type, GSVA)) +
  geom_violin(aes(fill = type), colour = "white",show.legend = F) +
  scale_fill_manual(values = c("#e3c4a8", "#4592af")) +
  geom_boxplot(colour = "gray", width = 0.1) +
  ylab("Expression of WFDC2") +
  xlab("GSE239897") +
  #scale_x_discrete(labels =  c("Metastasis tumor", "Primary tumor")) +
  stat_compare_means(comparisons = mycompare_list, method = "wilcox.test") +
  theme_classic()

###
set.seed(123)
library(fbroc)

map_table <- c("COPD" = 0, "Control" = 1)

type <- map_table[GSE239897_info]

final_data <- data.frame("WFDC12" = unlist(need_GSE239897_exp["WFDC12",]) ,
                    "WFDC3" = unlist(need_GSE239897_exp["WFDC3",]),
                    "WFDC2"= unlist(need_GSE239897_exp["WFDC2",]),
                    "type" = type)
library(pROC)
library(caret)

cutoff <- createDataPartition(final_data$type, p=0.50, list=FALSE)
# select 15% of the data for validation
testfinal_data <- final_data[-cutoff,]
# use the remaining 85% of data to training and testing the models
trainfinal_data <- final_data[cutoff,]


control <- trainControl(method="cv",number=10, classProbs = TRUE)
metric <- "Accuracy"
model <- train(type ~ WFDC2, data = trainfinal_data, method = "svmRadial",
               tuneLength = 8,preProc = c("center","scale"), 
               trControl=control)

predict <- predict(model, newdata = testfinal_data)

roc_result <- boot.roc(predict, as.logical(final_data[-cutoff,2]), n.boot = 100)
plot(roc_result)

roc1 <- roc(as.logical(final_data[-cutoff,2]),predict, levels = c(TRUE, FALSE))
a <- ci(roc1)
plot(ci.thresholds(roc1))

reportROC::reportROC(as.logical(final_data[-cutoff,2]), predict)


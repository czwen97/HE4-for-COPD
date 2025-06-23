setwd("D://work//2024//2024_08//维生素D与老年糖尿病患者心理状态的相关性研究//Data")

library(forestplot)
library(dplyr)
library(randomForest)
general <- read.csv("genal.csv", stringsAsFactors = F)


n_anx_index <- which(general$SAS水平 < 50)
y_anx_index <- which(general$SAS水平 >= 50)


n_dep_index <- which(general$SDS水平 < 53)
y_dep_index <- which(general$SDS水平 >= 53)

anx_type <- rep(0, 1139)
anx_type[y_anx_index] <- 1


dep_type <- rep(0, 1139)
dep_type[y_dep_index] <- 1


set.seed(123)
model_data = cbind(general[,c(3,6,10, 12,16,17,20)], dep_type)

otu_train.forest <- randomForest(dep_type ~ ., 
                                 data = model_data, importance = TRUE)
otu_train.forest


importance_otu <- otu_train.forest$importance
head(importance_otu)

#或者使用函数 importance()
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)

write.table(importance_otu, "dep_importance_otu.txt", col.names = F, row.names = T, sep = "\t", quote = F) ###输出的结果

#作图展示 top30 重要的 OTUs
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'variable importance')
###

model_data = cbind(general[,c(3)], anx_type)

otu_train.forest <- randomForest(anx_type ~ ., 
                                 data = model_data, importance = TRUE)


roc_result <- boot.roc(predict(otu_train.forest), as.logical(anx_type), n.boot = 500)
plot(roc_result)

roc1 <- roc(as.logical(otu_train.forest,anx_type), levels = c(TRUE, FALSE))
a <- ci(roc1)
plot(ci.thresholds(roc1))

reportROC::reportROC(as.logical(anx_type), predict(otu_train.forest))






###
上面为随机森林
下面为罗杰斯特回归模型
###


model_data =  genal[,c(5:7,14:16,19,23:26, 1)] ###选择要纳入模型的数据列


model <- glm( group ~ `The.number.of.lobes.involved`,
              data = model_data, family = binomial(link = "logit")) ###单因素逐个分析

summary(model)
OR = round(exp(coef(model))[2], digits = 3)
CI = round(exp (confint(model, level = 0.95))[2,], digits = 3)

str =  paste0(OR, "[", CI[1], "-", CI[2],"]")
str
###

model_data = cbind(general[,c(3,6,10, 12,16,17,20)], anx_type)

model <- glm( anx_type ~., 
              data = model_data, family = binomial(link = "logit"))

summary(model)
OR = round(exp(coef(model))[-1], digits = 3)
CI = round(exp (confint(model, level = 0.95))[-1,], digits = 3)

str = apply(
  cbind(OR, CI),
  1,
  function(x){
    return(paste0(x[1], "[", x[2], "-", x[3],"]"))
    
  }
)
str
info <- cbind(str, summary(model)$coefficients[-1,c(1,4)])
write.table(info, "anx_m_model_info.txt", row.names = F, col.names = F, sep = "\t", quote = F) 

cochrane_from_rmeta <- structure(list(mean  = c(NA, NA, OR), 
                                      lower = c(NA, NA, CI[,1]),
                                      upper = c(NA, NA, CI[,2])),
                                 .Names = c("mean", "lower", "upper"), 
                                 row.names = c(NA, -4L), 
                                 class = "data.frame")

tabletext <- cbind(c("Multivariate Logistic regression model of Anxiety", "Index", names(OR)),
                   c("", "OR[95%CI]", str),
                   c("", "p-value", round(as.numeric(info[,3]), digits = 4)))

p1 <- cochrane_from_rmeta %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(rep(TRUE, 2), rep(FALSE, 8), F),
             clip = c(0.1, 20), 
             xlog = TRUE, 
             col = fpColors(box = "#F27457",
                            line = "#03A696",
                            summary = "royalblue"))




###


model_data = cbind(general[,c(3,6,10, 12,16,17,20)], dep_type)

model <- glm( dep_type ~., 
              data = model_data, family = binomial(link = "logit"))

summary(model)
OR = round(exp(coef(model))[-1], digits = 3)
CI = round(exp (confint(model, level = 0.95))[-1,], digits = 3)

str = apply(
  cbind(OR, CI),
  1,
  function(x){
    return(paste0(x[1], "[", x[2], "-", x[3],"]"))
    
  }
)
str
info <- cbind(str, summary(model)$coefficients[-1,c(1,4)])
write.table(info, "dep_m_model_info.txt", row.names = F, col.names = F, sep = "\t", quote = F) 
# forest

cochrane_from_rmeta <- structure(list(mean  = c(NA, NA, OR), 
                                      lower = c(NA, NA, CI[,1]),
                                      upper = c(NA, NA, CI[,2])),
                                 .Names = c("mean", "lower", "upper"), 
                                 row.names = c(NA, -4L), 
                                 class = "data.frame")

tabletext <- cbind(c("Multivariate Logistic regression model of Depression", "Index", names(OR)),
                   c("", "OR[95%CI]", str),
                   c("", "p-value", round(as.numeric(info[,3]), digits = 4)))

p2 <- cochrane_from_rmeta %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(rep(TRUE, 2), rep(FALSE, 8), F),
             clip = c(0.1, 20), 
             xlog = TRUE, 
             col = fpColors(box = "#F27457",
                            line = "#03A696",
                            summary = "royalblue"))


###
library(fbroc)

final_data <-data.frame(general[,3], anx_type)
model <- glm( anx_type ~ . , 
              data = final_data, family = binomial(link = "logit"))



summary(model)
exp(coef(model))
exp (confint(model, level = 0.95))

need_coef = coef(model)

ResourceSelection::hoslem.test(model$y, fitted(model))

###

glm_function <- function(data) {
  data = as.numeric(data)
  pred = need_coef[1] +  data[1] * need_coef[2]
  return(c(pred, as.numeric(data[2])))
}

pred_result <-  t(apply(
  matrix(sample(1:1139, 1000), ncol = 1),
  1,
  function(index){
    result = glm_function(final_data[index,])
    result[which(result == 1)] <- TRUE
    result[which(result == 0)] <- FALSE
    return(result)
  }
))
library(pROC)


roc_result <- boot.roc(pred_result[,1], as.logical(pred_result[,2]), n.boot = 500)
plot(roc_result)

roc1 <- roc(as.logical(pred_result[,2]), pred_result[,1], levels = c(TRUE, FALSE))
a <- ci(roc1)
plot(ci.thresholds(roc1))

reportROC::reportROC(as.logical(pred_result[,2]), pred_result[,1])

### SVM
library(caret)

df <- data.frame(general[,3], dep_type)
cutoff <- createDataPartition(df$dep_type, p=0.50, list=FALSE)
# select 15% of the data for validation
testdf <- df[-cutoff,]
# use the remaining 85% of data to training and testing the models
traindf <- df[cutoff,]


control <- trainControl(method="cv",number=10, classProbs = TRUE)
metric <- "Accuracy"
model <- train(dep_type ~., data = traindf, method = "svmRadial",
               tuneLength = 8,preProc = c("center","scale"), 
              trControl=control)

predict <- predict(model, newdata = testdf)

roc_result <- boot.roc(predict, as.logical(df[-cutoff,2]), n.boot = 500)
plot(roc_result)

roc1 <- roc(as.logical(df[-cutoff,2]),predict, levels = c(TRUE, FALSE))
a <- ci(roc1)
plot(ci.thresholds(roc1))

reportROC::reportROC(as.logical(df[-cutoff,2]), predict)


rm(list = ls())
options(stringsAsFactors = F)
load("20240701GSE221921.Rdata")
rm(list = setdiff(ls(), c("filter_count","group")))
var=c("ARHGEF37","DYRK3","RGS17")
exp1=filter_count[rownames(filter_count) %in% var,]

##xgboost
library(shapviz)
library(xgboost)
library(ROCit)
library(tibble)
library(caret)
library(pROC)
xgboost_data=exp1
xgboost_data=as.data.frame(t(xgboost_data))
xgboost_data=xgboost_data[match(group$Sample,rownames(xgboost_data)),]
xgboost_data$group=group$Etiology
xgboost_data$group=as.factor(xgboost_data$group)
xgboost_data$group=as.numeric(xgboost_data$group)-1


set.seed(1234)
param <- list(objective = "binary:logistic",
              booster = "gbtree",
              eval_metric = "error",
              eta = 0.3,
              max_depth = 3,
              subsample = 1,
              colsample_bytree = 1,
              gamma = 0.5)

colnames(xgboost_data)
x <- as.matrix(xgboost_data[, c(1:3)])
y <- xgboost_data$group

mat <- xgb.DMatrix(data = x, label = y)
mat

library(caret)
grid <- expand.grid(nrounds = c(75, 100),
                    colsample_bytree = 1,
                    min_child_weight = 1,
                    eta = c(0.01, 0.1),
                    gamma = c(0.5, 0.25),
                    subsample = 0.5,
                    max_depth = c(2, 3))
cntrl <- trainControl(method = "cv",
                      number = 5,
                      verboseIter = F,
                      returnData = F,
                      returnResamp = "final")

set.seed(1)
train.xgb <- train(x = x,
                   y = y,
                   trControl = cntrl,
                   tuneGrid = grid,
                   method = "xgbTree")
train.xgb

param <- list(objective = "binary:logistic",
              booster = "gbtree",
              eval_metric = "error",
              eta = 0.01,
              max_depth = 2,
              subsample = 0.5,
              colsample_bytree = 1,
              gamma = 0.25,
              min_child_weight = 1,
              subsample = 0.5)


set.seed(1)
xgb.fit <- xgb.train(params = param, 
                     data = mat, 
                     nrounds = 100)
pred <- predict(xgb.fit, newdata = mat)
head(pred)
library(ROCR)
pred <- prediction(pred, xgboost_data$group)
perf <- performance(pred, "tpr", "fpr")
auc <- round(performance(pred, "auc")@y.values[[1]],digits = 4)

plot(perf, 
     main = paste("ROC curve (", "AUC = ",auc,")"), 
     col = 2, 
     lwd = 2)
abline(0,1, lty = 2, lwd = 2)

save.image(file="train.Rdata")

#Test
rm(list = ls())
options(stringsAsFactors = F)
load(file="step1output.Rdata")
test=exp
library(stringr)
Group=ifelse(str_detect(pd$`source_name_ch1`,"Control"),
             "control",
             "FM")
Group = factor(Group,levels = c("control","FM"))
Group
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL11532')
##8109161=ARHGEF37
#7909225=DYRK3
#8130394=RGS17
test=test[rownames(test) %in% c(8109161,7909225,8130394),]
rownames(test)=c("ARHGEF37","DYRK3","RGS17")
test=as.data.frame(t(test))
test$group=Group
test$group=as.numeric(test$group)-1

x <- as.matrix(test[, c(1:3)])
y <- test$group

test.mat <- xgb.DMatrix(data = x, label = y)


library(caret)
param <- list(objective = "binary:logistic",
              booster = "gbtree",
              eval_metric = "error",
              eta = 0.01,
              max_depth = 2,
              subsample = 0.5,
              colsample_bytree = 1,
              gamma = 0.25,
              min_child_weight = 1,
              subsample = 0.5)



set.seed(1)
xgb.fit <- xgb.train(params = param, 
                     data = test.mat, 
                     nrounds = 100)

pred <- predict(xgb.fit, newdata = test.mat)
head(pred)
library(ROCR)
pred <- prediction(pred, test$group)
perf <- performance(pred, "tpr", "fpr")
auc <- round(performance(pred, "auc")@y.values[[1]],digits = 4)

plot(perf, 
     main = paste("ROC curve (", "AUC = ",auc,")"), 
     col = 2, 
     lwd = 2)
abline(0,1, lty = 2, lwd = 2)



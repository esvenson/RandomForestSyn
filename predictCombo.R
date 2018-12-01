library(randomForest)
library(rsample)
library(pROC)

combination.tags = read.table("Gayvertetal2017PLOSCB_BRAF_synergy_labels.txt", header = FALSE, col.names = c('Combination', 'Type'))
single.drug.screen = read.table("Gayvertetal2017PLOSCB_singleagent_GI50.txt",
                                quote = "", 
                                sep = "\t",
                                header = TRUE)

## Prepare Training Set Matrix ##
cell.lines <- colnames(single.drug.screen)
drugs <- rownames(single.drug.screen)
train.data <- sapply(combination.tags$Combination, function(x){
  local.drugs <- unlist(strsplit(as.character(x), split = '_'))
  local.drugs.idx <- match(local.drugs, table = drugs)
  if(length(local.drugs.idx) == 2){
    local.dat <- sapply(1:length(cell.lines), simplify = FALSE, function(line.idx){
      return(setNames(c(mean(single.drug.screen[local.drugs.idx, line.idx]),
                        single.drug.screen[local.drugs.idx[1], line.idx] - single.drug.screen[local.drugs.idx[2], line.idx]), 
                      nm = c(paste(cell.lines[line.idx], '_mean', sep = ''),
                             paste(cell.lines[line.idx], '_diff', sep = ''))))
    })
    local.dat <- do.call(c, local.dat)
    return(local.dat)
  }else{
    return(rep(NA, length(cell.lines)))
  }
})
colnames(train.data) <- combination.tags$Combination
train.data <- na.roughfix(t(train.data))

## Create K-Fold ##
nRepeat <- 100
split.data <- vfold_cv(train.data, 
                       V = 10, 
                       repeats = nRepeat)

## Run Test ##
for(i in 1:nRepeat){
  local.train <- training(split.data$splits[[i]])
  local.train.labels <- combination.tags$Type[match(rownames(local.train),table = combination.tags$Combination)]
  local.test <- testing(split.data$splits[[i]])
  local.test.labels <- combination.tags$Type[match(rownames(local.test), table = combination.tags$Combination)]
  
  gmodel = randomForest(na.roughfix(local.train),
                        local.train.labels,
                        ntree = 1000,
                        mtry = sqrt(ncol(local.train)),
                        importance = TRUE)
  res <- predict(gmodel, newdata = na.roughfix(local.test))
  auc.res <- auc(roc(response = as.factor(as.numeric(local.test.labels == 'Synergistic')), predictor = as.numeric(res == 'Synergistic')))
  cat(sprintf('AUC: %.2f\n', auc.res[1]))
}

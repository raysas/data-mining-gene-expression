## ----installing, echo=FALSE, message=FALSE------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install("edgeR")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

if (!requireNamespace("caret", quietly = TRUE)) {
  install.packages("caret")
}

if (!requireNamespace("leaps", quietly = TRUE)) {
  install.packages("leaps")
}


## ----loading, message=FALSE---------------------------------------------------
library(tidyverse) #data anlaysis
library(GEOquery) #connecion to GEO data
library(biomaRt) #connection to ensembl
library(edgeR) #normalize seq counts
library(pROC) #plot the roc curve
library(caret) #for cv
library(leaps) #for subset selection
library(MASS)


#loading data accessed from GEO
gene_exp <- read.csv("../data/GSE231409_BRAVE_RNASeq_counts.csv")
rownames(gene_exp) <- gene_exp$target_id
gene_exp <- gene_exp %>%
  dplyr::select(-target_id)

#loading the metadata from GEO directly
gse <- getGEO("GSE231409", GSEMatrix =TRUE, getGPL=FALSE)
metadata <- pData(gse[[1]])

get_stat <- function(df) {
  cat("1st 5 columns:", head(colnames(df), n = 5), "\n")
  cat("1st 5 rows:", head(rownames(df), n = 5), "\n")
}
get_sample <- function(df){
  submatrix <- df[1:3, 1:3]
  print(submatrix)
}

##get_stat(gene_exp)
get_sample(gene_exp)


## ----normalization------------------------------------------------------------

dge <- DGEList(counts = gene_exp)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
gene_exp <- as.data.frame(normalized_counts)

#get_stat(gene_exp)
get_sample(gene_exp)



## ----var_genes----------------------------------------------------------------

sel_indices <- order(apply(gene_exp, 1, var), decreasing = TRUE)[1:100]
gene_exp <- gene_exp[sel_indices, , drop = FALSE]


## ----accessing_gene_names_ensembl---------------------------------------------
#Fixing row names: changing the ensemble accession IDs into meaningful gene names
rownames(gene_exp) <- sub("\\..*", "", rownames(gene_exp)) #removing the .* from the ensemble ids since that dows not affect the gene names
#dont worry about duplicates they will be handled ;)

ensembl_ids <- rownames(gene_exp)
gene_exp <- gene_exp %>% 
  rownames_to_column(var = "ens_accession") %>%
  rename(ensembl_gene_id = ens_accession) #we're in the process of changing the rownames

ensembl <- useEnsembl(biomart = "genes") #choosing the database genes 
connect <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #building a connection with ensembl to the database for human genes


#att <- listAttributes(connect) #listing the attributes to see which map for ids and names, just exploring around :p
#filters <- listFilters(connect) #listing the filters to see which map for ids and names
#retrieves with filter for the ensemble ids and returns the gene names only for the values that we provide - in this case the ensembl_ids

gene_names <- getBM(attributes = c("ensembl_gene_id","external_gene_name"), 
                    filters ="ensembl_gene_id", 
                    values = ensembl_ids, 
                    mart = connect)
gene_names$ensembl_gene_id <- as.character(gene_names$ensembl_gene_id)
gene_names$external_gene_name <- as.character(gene_names$external_gene_name)


#wanna get the list of ensembl ids that are not in gene_names
ensembl_ids_not_in_gene_names <- setdiff(ensembl_ids, gene_names$ensembl_gene_id)


#join tables and remove the ids that do not map to gene names
gene_exp <- merge(gene_exp, gene_names, by = "ensembl_gene_id", all.x = TRUE)
gene_exp <- gene_exp %>% 
  filter(!is.na(gene_exp$external_gene_name))


#adding a prefix to duplicate external_gene_name
gene_exp <- gene_exp %>%
  mutate(external_gene_name = ifelse(duplicated(external_gene_name), 
                                     paste(external_gene_name, row_number(), sep = "_"), 
                                     external_gene_name))


rownames(gene_exp) <- gene_exp$external_gene_name
gene_exp <- gene_exp[,-1] %>%
  dplyr::select(-external_gene_name)

#get_stat(gene_exp)
get_sample(gene_exp)



## ----transpose----------------------------------------------------------------
gene_exp <- t(gene_exp) %>%
  as.data.frame()

#get_stat(gene_exp)
get_sample(gene_exp)


## ----add_labels---------------------------------------------------------------

#get the labels of the metadata in a separate column: data frame contains title, accession number and covid result
response <- metadata %>%
  dplyr::select(1,16) %>%
  rename(covid = characteristics_ch1.6) %>%
  mutate(covid = gsub("covid: ", "" ,covid))

gene_exp <- gene_exp %>%
  rownames_to_column(var = "title") %>% #convert the row names into a column
  left_join(., response, by = "title")  
rownames(gene_exp) <- gene_exp$title
gene_exp <- gene_exp[,-1]

gene_exp <- gene_exp[!is.na(gene_exp$covid), ] #to remove unlabeled samples, NA as response



## ----factor_response----------------------------------------------------------
gene_exp$covid <- factor(gene_exp$covid, labels=c("No","Yes"))
#colnames(gene_exp)[18383]


## ----useful_functions---------------------------------------------------------

measures = function(m){
  spe <- (m[2,2] / (m[1,2]+m[2,2]))
  sen <- (m[1,1] / (m[1,1]+m[2,1]))
  a <- (sum(diag(m)) / sum(m))
  
  cat("Specificity: ", spe, "\n")
  cat("Sensitivity: ", sen, "\n")
  cat("Accuracy: ", a, "\n")
  
  return (c(se=sen, sp=spe, acc=a))
}#helper method to use in the next ne

print_confusion_and_measures <- function(model, test_data) {
  #only prints does not return !!
  pred <- predict(model, newdata = test_data)
  confusion <- table(pred$class, test_data$Outcome)
  print(confusion)
  measures(confusion)
}

plot_roc_curves <- function(model1, model2, test_data) {
  #this function is intended to compare different models (qda, lda) models curve on the same plot: 
  #run it like this plot_roc_curves(model1, model2, model3, testing)
  
  model1_name <- deparse(substitute(model1))
  model2_name <- deparse(substitute(model2))
  
  pred1 <- predict(model1, newdata = test_data)
  pred2 <- predict(model2, newdata = test_data)
  
  roc1 <- roc(response = test_data$covid, predictor = pred1$posterior[,2])
  roc2 <- roc(response = test_data$covid, predictor = pred2$posterior[,2])
  
  auc1 <- auc(roc1)
  auc2 <- auc(roc2)
  
  print(paste(model1_name, "AUC:", auc1))
  print(paste(model2_name, "AUC:", auc2))
  
  ggroc(list(Model1 = roc1, Model2 = roc2))
}



split_data <- function(df, p=0.7) {
  #this function takes a data frame and proportion of split and splits it to training and testing - if proportion not provided ir'll assume it is 0.7
  df <- df %>%
    mutate(id=row_number())

  tr <- df %>%
    slice_sample(prop=p)
  te <- anti_join(df, tr, by='id')

  tr <- tr %>% dplyr::select(-id)
  te <- te %>% dplyr::select(-id)
  
  return (list(train = tr, test = te))
}



## ----forward_selection--------------------------------------------------------


forward_selection <- regsubsets(covid ~ ., data= gene_exp, nvmax = ncol(gene_exp) - 1, method = "forward")

cp_summary <- summary(forward_selection)$cp
min_cp_index <- which.min(cp_summary)
min_cp_value <- cp_summary[min_cp_index]
plot(cp_summary, type='l')
points(min_cp_index, min_cp_value, col='red', pch=19)
cat("Number of features at the lowest Cp:", min_cp_index - 1, "\n")

adjr2_summary <- summary(forward_selection)$adjr2
max_adjr2_index <- which.max(adjr2_summary)
max_adjr2_value <- adjr2_summary[max_adjr2_index]
plot(adjr2_summary, type='l', main='Adjusted R-squared')
points(max_adjr2_index, max_adjr2_value, col='red', pch=19)
cat("Number of features at the highest Adjusted R-squared:", max_adjr2_index - 1, "\n")


bic_summary <- summary(forward_selection)$bic
min_bic_index <- which.min(bic_summary)
min_bic_value <- bic_summary[min_bic_index]
plot(bic_summary, type='l', main='BIC')
points(min_bic_index, min_bic_value, col='blue', pch=19)
cat("Number of features at the lowest BIC:", min_bic_index, "\n")



## ----seven_genes--------------------------------------------------------------
coef(forward_selection, 6)



## ----val_set_approach---------------------------------------------------------

set.seed(42)
my_split <- split_data(gene_exp, 0.8)
train_set <- my_split$train
test_set <- my_split$test
#View(train_set)
#View(test_set)


## ----one_gene-----------------------------------------------------------------

coef(forward_selection, 1)

logistic_model <- glm(covid ~ IFI27, data = train_set, family = "binomial")
#Perform the logistic model on the 6 features that we chose using forward selection
summary(logistic_model)

pred <- predict(logistic_model, newdata = test_set)
threshold <- 0.5  #probability
predicted_labels <- ifelse(pred > threshold, 1, 0)
confusion_matrix <- table(predicted_labels, test_set$covid)
print(confusion_matrix)


measures(confusion_matrix)



## ----logistic_regression_forward_val------------------------------------------

logistic_model <- glm(covid ~ RPS2+IFI27+EEF2+`MT-CO2` + `MT-CO3` + IGHM, data = train_set, family = "binomial")
#Perform the logistic model on the 6 features that we chose using forward selection
summary(logistic_model)

pred <- predict(logistic_model, newdata = test_set)
threshold <- 0.5  #probability
predicted_labels <- ifelse(pred > threshold, 1, 0)
confusion_matrix <- table(predicted_labels, test_set$covid)
print(confusion_matrix)


measures(confusion_matrix)



## ----lda_forward_val----------------------------------------------------------

lda_model <- lda(covid ~ RPS2+IFI27+EEF2+`MT-CO2` + `MT-CO3` + IGHM, data = train_set)
#summary(lda_model)

pred <- predict(lda_model, newdata = test_set)
confusion <- table(pred$class, test_set$covid)
confusion
measures(confusion)



## ----qda_model_forward_val----------------------------------------------------

qda_model <- qda(covid ~ RPS2+IFI27+EEF2+`MT-CO2` + `MT-CO3` + IGHM, data = train_set)
#summary(lda_model)

pred <- predict(qda_model, newdata = test_set)
confusion <- table(pred$class, test_set$covid)
confusion
measures(confusion)



## ----lda_qda_roc_-------------------------------------------------------------
plot_roc_curves(lda_model, qda_model, test_set)


## -----------------------------------------------------------------------------
train_control <- trainControl(method = "cv", number = 5)
#View(train_control)
logistic_cv <- train(covid ~RPS2+IFI27+EEF2+`MT-CO2`+`MT-CO3`+IGHG1, data = gene_exp, method = "glm", family=binomial, trControl = train_control, metric="Accuracy")
logistic_cv


## -----------------------------------------------------------------------------

lda_cv <- train(covid ~RPS2+IFI27+EEF2+`MT-CO2`+`MT-CO3`+IGHG1, data = gene_exp, method = "lda", trControl = train_control, metric="Accuracy")
lda_cv



## ----qda_cv-------------------------------------------------------------------

qda_cv <- train(covid ~RPS2+IFI27+EEF2+`MT-CO2`+`MT-CO3`+IGHG1, data = gene_exp, method = "lda", trControl = train_control, metric="Accuracy")
qda_cv



## ----loocv--------------------------------------------------------------------

loocv <- trainControl(method = "LOOCV")



## ----log_reg_loocv------------------------------------------------------------

logistic_model <- train( covid ~ RPS2+IFI27+EEF2+`MT-CO2` + `MT-CO3` + IGHG1, data=gene_exp, method="glm", family=binomial(link=logit), metric="Accuracy" , trControl=loocv)

logistic_model



## ----lda_loocv----------------------------------------------------------------

lda_model  <- train( covid ~ RPS2+IFI27+EEF2+`MT-CO2` + `MT-CO3` + IGHG1, data=gene_exp, method="lda", metric="Accuracy" , trControl=loocv)

lda_model



## ----qda_loocv----------------------------------------------------------------

qda_model  <- train( covid ~ RPS2+IFI27+EEF2+`MT-CO2` + `MT-CO3` + IGHG1, data=gene_exp, method="qda", metric="Accuracy" , trControl=loocv)

qda_model
#plot_roc_curves(lda_model, qda_model, test_set)



## ----back_sub-----------------------------------------------------------------

backward_sel <- regsubsets(covid ~ ., data= gene_exp, nvmax = ncol(gene_exp) - 1, method="backward")

bic_summary <- summary(backward_sel)$bic
min_bic_index <- which.min(bic_summary)
min_bic_value <- bic_summary[min_bic_index]
plot(bic_summary, type='l', main='BIC')
points(min_bic_index, min_bic_value, col='red', pch=19)



## ----imp_genes_back-----------------------------------------------------------

p <- min_bic_index
coef(backward_sel, p)



## ----backward_resampling_val--------------------------------------------------

#the data split is still available as train_set and test_set
logistic_model <- glm(covid ~ CD74+MYH9+COTL1+FKBP8+ OAS3 +  LCP1 +TLN1 + RPS11 + RPL13A + EEF1A1 + EEF2 + IFIT1 + `MT-ND4` + `MT-CO3` + UBA52, data = train_set, family = "binomial")
summary(logistic_model)

pred <- predict(logistic_model, newdata = test_set)
threshold <- 0.5  #probability
predicted_labels <- ifelse(pred > threshold, 1, 0)
confusion_matrix <- table(predicted_labels, test_set$covid)
measures(confusion_matrix)




## ----log_reg_cv_back----------------------------------------------------------

cv <- trainControl(method = "cv", n=5)

logistic_model <- train( covid ~ CD74 +MYH9+COTL1+FKBP8+ OAS3 +  LCP1 +TLN1 + RPS11 + RPL13A + EEF1A1 + EEF2 + IFIT1 + `MT-ND4` + `MT-CO3` + UBA52, data=gene_exp, method="glm", family=binomial(link=logit), metric="Accuracy" , trControl=cv)

logistic_model



## ----log_reg_loocv_back-------------------------------------------------------

logistic_model <- train( covid ~ CD74+MYH9+COTL1+FKBP8+ OAS3 +  LCP1 +TLN1 + RPS11 + RPL13A + EEF1A1 + EEF2 + IFIT1 + `MT-ND4` + `MT-CO3` + UBA52, data=gene_exp, method="glm", family=binomial(link=logit), metric="Accuracy" , trControl=loocv)

logistic_model



## -----------------------------------------------------------------------------
#best_selection <- regsubsets(covid ~ ., data= gene_exp, nvmax = ncol(gene_exp) - 1)


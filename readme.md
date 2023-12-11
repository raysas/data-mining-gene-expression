# Data Mining Project

This is a mutliphase projet for the Data Mining course at LAU in Fall 2023. The project is divided into 3 phases:

- [Phase 1](https://github.com/raysas/CO2-emissions-prediction): Regression
- [Phase 2](https://github.com/raysas/data-mining-gene-expression/phase-2): Classification
- [Phase 3](https://github.com/raysas/data-mining-gene-expression/phase-3): Clustering

Phase 1 is in a seperate [repo](https://github.com/raysas/CO2-emissions-prediction). Phase 2 and 3 are done using the same dataset, withdrwan from GEO. The dataset is a collection of _RNA-seq raw counts_ for 455 samples, involving 60705 genes, Out of the 455, 417 have  clear covid status (positive/negative). The dataset is available in the data folder. 

## Phase 2: Classification

In this phase, we used preprocessing steps to clean the data, and then used 3 different classification algorithms to classify the samples into positive and negative covid status. The algorithms used are:
- logistic regression
- Linear Discriminant Analysis
- Quadratic Discriminant Analysis

We used different resampling techniques to evaluate the performance of the models and compare test errors. The resampling techniques used are:
- Validation set approach 80/20
- 5-fold cross validation (empirical k=5)
- Leave one out cross validation (LOOCV)

Feature selection was applied as well to reduce the number of features and improve the performance of the models. The feature selection techniques used are:
- Forward selection
- Backward selection

Note that out of the 60K genes, ~18K was left after normalization and preprocessing. Then we applied the highest variance technique - i.e. got the top 100 genes that showed between samples variance. Then we applied the feature selection techniques on the 100 genes. In an ideal situation we would've applied dimensionality reduction (PCA) which bases on the same idea. However, this will be applied in phase 3.

Details on this phase [here](https://github.com/raysas/data-mining-gene-expression/phase-2).

## Phase 3: Classification and Clustering

In this phase we extended classification to include Decision Trees (DT) and moved to unsupervided learning approaches including Dimensionality Reduction (DR) and Clustering.
% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/BRLasso.R
\name{BRLasso}
\alias{BRLasso}
\title{Bootstrap ranking LASSO method.}
\usage{
BRLasso(x, y, B = 5, Boots = 100, kfold = 10)
}
\arguments{
\item{x}{the predictor matrix}

\item{y}{the response variable, a factor object with values of 0 and 1}

\item{B}{the external loop for intersection operation, with the default value 5}

\item{Boots}{the internal loop for bootstrap sampling, with the default value 100}

\item{kfold}{the K-fold cross validation, with the default value 10}
}
\description{
This function performs a bootstrap ranking LASSO Logistic regression model for variable selection.
}
\examples{
library(datasets)
head(iris)
X <- as.matrix(subset(iris,iris$Species!="setosa")[,-5])
Y <- as.factor(ifelse(subset(iris,iris$Species!="setosa")[,5]=='versicolor',0,1))
# Fitting a BRLasso Logistic regression model
BRLasso.fit <- BRLasso(x=X, y=Y, B=2, Boots=10, kfold=10)
# Variables selected
BRLasso.fit$Var_selected
# Coefficients of the selected variables
BRLasso.fit$Asso
}

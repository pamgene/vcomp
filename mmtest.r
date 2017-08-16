library(bnshiny)
library(bnutil)
library(lme4)
library(dplyr)
library(ggplot2)
library(data.table)
library(vcomp2)
library(pgMulticore)


bndata = getData()
mdl = formula("value ~ 1+FFactor1 +(FFactor1|SFactor6) + (1|CubeID)")
result = vcomp2::modelOperator(bndata$data, model = mdl, nomfac = c("SFactor6", "CubeID") )



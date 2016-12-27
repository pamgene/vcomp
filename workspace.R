library(bnutil)
library(vcomp)
library(shiny)
library(lme4)
library(dplyr)
library(data.table)
library(ggplot2)

getData = function() AnnotatedData$new(data=vcdata, metadata=vcmeta)

getFolder = function() file.path(getwd(), 'run')

setResult = function(annotatedResult){

}

bnMessageHandler = bnshiny::BNMessageHandler$new()
bnMessageHandler$getFolderHandler = getFolder
bnMessageHandler$getRunFolderHandler = getFolder
bnMessageHandler$getDataHandler = getData
bnMessageHandler$setResultHandler = setResult

bnshiny::startBNTestShiny('vcomp', sessionType='run', bnMessageHandler=bnMessageHandler)

#' @import plyr
#' @import dplyr
#' @import bnutil
#' @import ggplot2
#' @import lme4
#' @import data.table
#' @export
shinyServerRun = function(input, output, session, context) {

  getFolderReactive = context$getRunFolder()
  getDataReactive = context$getData()
  getSettingsReactive = context$getFolder()

  output$body = renderUI({
    sidebarLayout(
      sidebarPanel(
        tags$div(HTML("<strong><font color = blue>Variance Component Analysis</font></strong>")),
        tags$hr(),
        actionButton("start", "Run model"),
        tags$hr(),
        checkboxInput("reml", "REML", value = TRUE),
        textInput("model", "value ~", value = "1"),
        actionButton("add", "Add term for:"),
        selectInput("terms", "", choices = list(), multiple = TRUE),
        tags$hr(),
        selectInput("nominal", "Categorical factors:", choices = list(), multiple = TRUE),
        tags$hr(),
        actionButton("done", "Done")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Model Output",
                   selectInput("pid", "Select peptide", choices = list()),
                   verbatimTextOutput("modelOutput"),
                   tableOutput("ic"),
                   tableOutput("dataInput")
          ),
          tabPanel("Random Factors",
                   selectInput("graphtype", "Graph Type", choices = c("CV", "relative", "SD")),
                   checkboxInput("logx" , "Logarithmic x-axis", value = FALSE),
                   splitLayout(
                     textInput("xmin",    label = "x-axis lower limit", value = "0"),
                     textInput("xmax",    label=  "x-axis upper limit", value = "auto"),
                     textInput('ymin',    label = "y-axis lower limit", value = "0"),
                     textInput('ymax',    label = "y-axis upper limit", value = "1")
                   ),
                   plotOutput("graphout", height = 700),
                   actionLink("png", "Save graph"),
                   verbatimTextOutput("status")
          ),

          tabPanel("Fixed Factors",
                   plotOutput("fxdout", height = 700),
                   actionLink("fxdpng", "Save graph"),
                   verbatimTextOutput("fxdstatus")

          )
        )
      )
    )
  })

  observe({
    getFolder = getFolderReactive$value
    if (is.null(getFolder)) return()
    folder = getFolder()

    getData=getDataReactive$value
    if (is.null(getData)) return()
    bndata = getData()

    getSettings = getSettingsReactive$value
    if (is.null(getSettings)) return()
    settingsFile = file.path(getSettings(), "settings.RData")

    if(file.exists(settingsFile)){
      settings = load(settingsFile)
      if(exists("aModelSpec")){
        updateTextInput(session, "model", value =aModelSpec)
      }
      if(exists("aNomFac")){
        updateSelectInput(session, "nominal", selected = aNomFac)
      }
      if(exists("aReml")){
        updateCheckboxInput(session, "reml", value = aReml)
      }
    }

    factorColumnLabels = bndata$arrayColumnNames

    if (bndata$hasXAxis){
      factorColumnLabels = c(factorColumnLabels, bndata$xAxisColumnName)
    }

    updateSelectInput(session, "terms", choices = factorColumnLabels, selected = factorColumnLabels[1])

    bNum = vector(length = length(factorColumnLabels))
    for(i in 1:length(bNum)){
      bNum[i] = class(bndata$data[[ factorColumnLabels[i] ]]) == "numeric"
    }
    updateSelectInput(session, "nominal", choices = factorColumnLabels, selected = factorColumnLabels[!bNum])

    observeEvent(input$add, {
      newterms = paste( "(1|", input$terms, ")", sep = "", collapse = "+")
      newmodel = paste(input$model, newterms, sep = "+")
      updateTextInput(session, "model", value = newmodel)
    })

    observeEvent(input$done, {
      res = modelReactive()
      cRes = getdfout(res)
      cRes = cRes[c("rowSeq", "colSeq", "cValue", "residuals")]
      meta = data.frame(labelDescription = c("rowSeq","colSeq", "cValue", "residuals"),
                        groupingType = c("rowSeq", "colSeq", "QuantitationType", "QuantitationType"))

      result = AnnotatedData$new(data = cRes, metadata = meta)
      context$setResult(result)
    }, ignoreInit = TRUE)


    modelReactive = reactive({
      input$start
      isolate({
         aModelSpec = input$model
         aNomFac    = input$nominal
         aReml      = input$reml
      })
      save(file = settingsFile, aModelSpec, aNomFac, aReml)
      models = modelOperator(bndata$data, model = formula(paste("value ~",aModelSpec) ), nomfac = aNomFac, reml =  aReml)
      pidList = models$rowSeq
      names(pidList) = models$ID
      updateSelectInput(session, "pid", choices = pidList, selected = pidList[1])
      return(models)
    })

    vcReactive = reactive({
      vcdf = getVarComp(modelReactive())
    })

    fxdReactive = reactive({
      fxddf = getFxdComp(modelReactive())
    })

    output$modelOutput = renderPrint({
      if (input$start == 0) return(".")
      res = modelReactive()
      if(!is.null(res$rowSeq)){
        this = res %>% filter(rowSeq == input$pid)
        return(this$aLme)
      } else {
        return(".")
      }

    })

     output$ic = renderTable({
      if (input$start == 0) return(NULL)
      res = modelReactive()
      this = res %>% filter(rowSeq == input$pid)
      if (length(this$aLme)>0){
        df = data.frame("Parameter" = c("AIC", "BIC"), value = c(AIC(this$aLme[[1]]), BIC(this$aLme[[1]])))
        return(df)
      } else {
        return(NULL)
      }
    })

    output$dataInput = renderTable({
      if (input$start == 0) return(NULL)
      res = modelReactive()
      this = res %>% filter(rowSeq == input$pid)
      if (length(this$aLme)>0){
        idf = this$df[[1]]
        isolate({
          bShowCol = apply(as.matrix(colnames(idf)), MARGIN = 1, function(x)grepl(x, input$model))
          bShowCol = bShowCol | colnames(idf) %in% c("value","rowSeq", "colSeq")
        })
        return(idf[,bShowCol])
      } else {
        return(NULL)
      }
    })

    plotReactive = reactive({
      vcdf = vcReactive()
      if(input$graphtype == "SD"){
        aPlot = ggplot(vcdf, aes(x = y0, y = s, colour = comp)) + ylab("SD")
      } else if (input$graphtype == "relative"){
        aPlot = ggplot(vcdf, aes(x = y0, y = r, colour = comp)) + ylab("Relative Var")
      } else{
        aPlot = ggplot(vcdf, aes(x = y0, y = cv, colour = comp)) + ylab("CV")
      }

      xLim = as.numeric(c(input$xmin, input$xmax))
      yLim = as.numeric(c(input$ymin, input$ymax))

      aPlot = aPlot + ylim(yLim)
      aPlot = aPlot + xlim(xLim)

      if(input$logx){
        if (!any(is.na(xLim)) & all(xLim > 0)){
          aPlot = aPlot + scale_x_log10(limits = xLim)
        } else {
          aPlot = aPlot + scale_x_log10()
        }
      }
      aPlot = aPlot + geom_point() + facet_wrap(~comp)

    })

    fxdPlotReactive = reactive({
      df = fxdReactive()
      aPlot = ggplot(df, aes(x = y0, y = fxd, colour = comp)) + geom_point()
      aPlot = aPlot + facet_wrap(~comp, scales = "free_y")
    })

    save.png = reactive({
      filename = file.path(getFolder(), paste("VC", format(Sys.time(), "%Y%m%d-%H%M.png")) )
      aPlot = plotReactive()
      ggsave(aPlot, file= filename)
    })

    save.fxdpng = reactive({
      filename = file.path(getFolder(), paste("FXD", format(Sys.time(), "%Y%m%d-%H%M.png")) )
      aPlot = fxdPlotReactive()
      ggsave(aPlot, file= filename)
    })

    output$graphout = renderPlot({
      if (input$start == 0) return()
      aPlot = plotReactive()
      print(aPlot)
    })

    output$fxdout = renderPlot({
      if (input$start == 0) return()
      aPlot = fxdPlotReactive()
      print(aPlot)
    })

    output$status = renderText({
      if (input$png > 0){
        save.png()
        return("Plot saved as image in results folder")
      }
      return(".")
    })

    output$fxdstatus = renderText({
      if (input$fxdpng > 0){
        save.fxdpng()
        return("Plot saved as image in results folder")
      }
      return(".")
    })

  })
}

#' @export
shinyServerShowResults = function(input, output, session, context){
  getFolderReactive = context$getRunFolder()

  output$body = renderUI({
    mainPanel(
      tags$hr(),
      tags$div(HTML("<strong><font color = blue>Find your saved plots in the results folder</font></strong>")),
      tags$hr(),
      actionLink("resfolder", "Open results folder")
    )
  })

  observe({
    getFolder = getFolderReactive$value
    if (is.null(getFolder)) return()
    observeEvent(input$resfolder, {
      shell.exec(getFolder())
    })
  })
}

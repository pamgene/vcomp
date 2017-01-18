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
        textInput("model", "value ~", value = "1"),
        actionButton("add", "Add term for:"),
        selectInput("terms", "", choices = list(), multiple = TRUE),
        tags$hr(),
        actionButton("done", "Done")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Model Output",
                   selectInput("pid", "Select peptide", choices = list()),
                   verbatimTextOutput("modelOutput")
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
      updateTextInput(session, "model", value = modelSpec)
    }

    arrayColumnLabels = bndata$arrayColumnNames
    updateSelectInput(session, "terms", choices = arrayColumnLabels, selected = arrayColumnLabels[1])


    observeEvent(input$add, {
      newterms = paste( "(1|", input$terms, ")", sep = "", collapse = "+")
      newmodel = paste(input$model, newterms, sep = "+")
      updateTextInput(session, "model", value = newmodel)
    })

    observeEvent(input$done, {
      # dummy return
      meta = data.frame(labelDescription = c("rowSeq","colSeq", "dummy"),
                        groupingType = c("rowSeq", "colSeq", "QuantitationType"))

      result = AnnotatedData$new(data = data.frame(rowSeq = 1, colSeq = 1, dummy = NaN), metadata = meta)
      context$setResult(result)
    })
    modelReactive = reactive({
      input$start
      isolate({
        modelSpec = input$model
        save(file = settingsFile, modelSpec)
        models = modelOperator(bndata$data, model = formula(paste("value ~",input$model) ))
        pidList = paste(models$ID, " (",models$rowSeq,")", sep = "")
        updateSelectInput(session, "pid", choices = pidList, selected = pidList[1])
        result = list(models, pidList)
      })
      return(result)
    })

    vcReactive = reactive({
      input$start
      res = modelReactive()
      vcdf = getVarComp(res[[1]])
    })

    fxdReactive = reactive({
      input$start
      res = modelReactive()
      fxddf = getFxdComp(res[[1]])
    })

    output$modelOutput = renderPrint({
      if (input$start == 0) return(".")
      res = modelReactive()
      bIdx = input$pid == res[[2]]
      models = res[[1]]
      return(models[bIdx,]$V1[[1]])
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

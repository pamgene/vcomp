#' @import dplyr
#' @import bnutil
#' @import ggplot2
#' @import lme4
#' @import data.table
#' @export
shinyServerRun = function(input, output, session, context) {

  getFolderReactive = context$getFolder()
  getDataReactive = context$getData()

  output$body = renderUI({
    sidebarLayout(
      sidebarPanel(
        tags$div(HTML("<strong><font color = blue>Variance Component Analysis</font></strong>")),
        tags$hr(),
        textInput("model", "value ~", value = "1"),
        actionButton("add", "Add term for:"),
        selectInput("terms", "", choices = list(), multiple = TRUE),
        actionButton("start", "Start")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Model Output",
                   selectInput("pid", "Select peptide", choices = list()),
                   verbatimTextOutput("modelOutput")
          ),
          tabPanel("Graph",
                   selectInput("graphtype", "Graph Type", choices = c("CV", "relative", "SD")),
                   checkboxInput("logx" , "Logarithmic x-axis", value = FALSE),
                   splitLayout(
                     textInput("xmin",    label = "x-axis lower limit", value = "0"),
                     textInput("xmax",    label=  "x-axis upper limit", value = "auto"),
                     textInput('ymin',    label = "y-axis lower limit", value = "0"),
                     textInput('ymax',    label = "y-axis upper limit", value = "1")
                   ),
                   plotOutput("graphout", height = 700)
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

    arrayColumnLabels = bndata$arrayColumnNames
    updateSelectInput(session, "terms", choices = arrayColumnLabels, selected = arrayColumnLabels[1])

    observeEvent(input$add, {
      newterms = paste( "(1|", input$terms, ")", sep = "", collapse = "+")
      newmodel = paste(input$model, newterms, sep = "+")
      updateTextInput(session, "model", value = newmodel)
    })

    modelReactive = reactive({
      input$start
      isolate({
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

    output$modelOutput = renderPrint({
      if (input$start == 0) return(".")
      res = modelReactive()
      bIdx = input$pid == res[[2]]
      models = res[[1]]
      return(models[bIdx,]$V1[[1]])
    })

    output$graphout = renderPlot({
      if (input$start == 0) return()
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
      print(aPlot)
    })

  })
}


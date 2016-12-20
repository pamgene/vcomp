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
        textInput("model", "value ~", value = "1 + (1|SFactor6)"),
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
              selectInput("graphtype", "Graph Type", choices = c("SD", "relative", "CV")),
              plotOutput("graphout")
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
      vcdf = vcReactive()
      if(input$graphtype == "SD"){
        aPlot = ggplot(vcdf, aes(x = y0, y = s, colour = comp)) + ylab("SD")
      } else if (input$graphtype == "relative"){
        aPlot = ggplot(vcdf, aes(x = y0, y = r, colour = comp)) + ylab("Relative Var")
      } else{
        aPlot = ggplot(vcdf, aes(x = y0, y = r, colour = comp)) + ylab("CV")
      }
      aPlot = aPlot + geom_point() + ylim(c(0, 0.5)) + facet_wrap(~comp)
      print(aPlot)
    })

  })
}


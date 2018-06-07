library(shiny)
library(markdown)
library(dplyr)
library(ggplot2)
library(mzID)
library(shinyjs)
source("helper.R")

shinyUI(fluidPage(
  useShinyjs(),
  tags$style(appCSS),
           sidebarLayout(
             sidebarPanel(width = 4,
                          HTML(markdownToHTML(text =
'
Import mzid file and assess decoy quality.
'
                          )), hr(),
                          fileInput('file1','Choose mzid File',accept = c('.mzid')),
                          withBusyIndicatorUI(actionButton("processFileBtn","Process mzid file",class = "btn-primary")),
                          checkboxInput("log", '-log10 transform score?', TRUE),
 			  numericInput("nBreaks", "# histogram bins:", 50, min = 1, max = 200),
        numericInput("alpha", "FDR significance level", 0.01, min = 0, max = 1,step=.01),
                          HTML(markdownToHTML(text =
'')),
        uiOutput("ColumnSelector")
),
mainPanel(width = 8,  plotOutput('histPlot'),plotOutput('ppPlot'))
           ))
)

library(shiny)
library(markdown)
library(dplyr)
library(ggplot2)
library(mzID)
library(shinyjs)
source("helper.R") # Load all the code needed to show feedback on a button click

shinyServer(function(input, output, session) {
#######################
  ## read data from given path
  dataFlat <- eventReactive(input$processFileBtn,{
  withBusyIndicatorServer("processFileBtn", {
  flatten(mzID(req(input$file1$datapath)))
  })
  })
  output$ColumnSelector <- renderUI({
	  tagList(
    selectInput("decoys","select the decoy column", choices = colnames(dataFlat())),

    selectInput("score","select the score column", choices = colnames(dataFlat())) )
  }
  )
  data <- reactive({
      data <- dataFlat()[,c(input$decoys,input$score)]
      data <- na.exclude(data)
      names(data)<-c("decoy","score")
      data$score<-as.double(data$score)
      if (input$log) data$score<--log10(data$score+1e-100)
      data <- data %>% arrange(desc(score)) %>% mutate(fdr=cumsum(decoy)/cumsum(!decoy))
      data$fdr[data$decoy] <- NA
      data
  })
####################
  ##plots
#######################
  output$histPlot <- renderPlot({
	binwidth <- diff(range(data()$score,na.rm=TRUE))/input$nBreaks
  nPSMsig <- sum(data()$fdr<=input$alpha,na.rm=TRUE)
  scoreThresh <- min(data()$score[data()$fdr<=input$alpha],na.rm=TRUE)
    if (nPSMsig>0){
    return(ggplot(data(), aes(score, fill = decoy, col=I("black")))+ geom_histogram(alpha = 0.5, binwidth=binwidth, position = 'identity') +  labs(x = 'Score', y = 'Counts' ,title = paste0('Histogram of targets and decoys\n',nPSMsig,' PSMs significant at ',input$alpha*100,'% FDR level')) +
    geom_vline(xintercept=scoreThresh,colour="red", size=1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.2)),
      axis.text = element_text(size = rel(1.2)),
      axis.title.y = element_text(angle = 0)))
      } else return(ggplot(data(), aes(score, fill = decoy, col=I("black")))+ geom_histogram(alpha = 0.5, binwidth=binwidth, position = 'identity') +  labs(x = 'Score', y = 'Counts' ,title = paste0('Histogram of targets and decoys\n',nPSMsig,' PSMs significant at ',input$alpha*100,'% FDR level')) +
      theme_bw() +
      theme(
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        axis.title.y = element_text(angle = 0)))
  })
  output$ppPlot <- renderPlot({
	if (class(data()$decoy)=="logical"&class(data()$score)=="numeric"){
    pi0 <- sum(data()$decoy)/sum(!data()$decoy)
    ppPlot <- ggplot()  +
    geom_abline(slope = pi0,color = 'black') +
    labs(x = 'Decoy Percentile', y = 'Target\nPercentile' ,title = 'PP plot of target PSMs') +
    xlim(0,1) + ylim(0,1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.2)),
      axis.text = element_text(size = rel(1.2)),
      axis.title.y = element_text(angle = 0))

  x <- data()$score[!data()$decoy]
  Ft <- ecdf(x)
  Fd <- ecdf(data()$score[data()$decoy])
  df <- data_frame(Fdp = Fd(x), Ftp = Ft(x))

  ppPlot + geom_point(data = df,aes(Fdp,Ftp),color = 'dark grey')
   }
  })

}
)

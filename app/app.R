library(shiny)
library(shinydashboard)
library(httr)
library(globals)
library(future)

source("analysis.R")

size = 450

ui<-dashboardPage(
  dashboardHeader(title = "CAD Macrophage ATAC-seq"),
  dashboardSidebar(textInput(inputId = "Gene",
                             label = "Gene Symbol",
                             value = "STAT1"),
                   sliderInput(inputId = "Flank",
                               label = "Zoom",min = 0,max = 100000,value = 10000)
                   # ,actionButton(inputId = "neg10K", label = "< 10k")
                   # ,actionButton(inputId = "pos10K", label = "> 10k")
  ),
  dashboardBody(
    fluidRow(
    box(title = "Healthy Vs CAD",
    # the MA-plot
    plotOutput("plotma", click="plotma_click"),# width=size, height=size),
    collapsible = TRUE),
    box(title = "Normalized Counts",
    # # the counts plot for a single gene
    plotOutput("plotcounts"),# width=size, height=size),
    #
      collapsible = TRUE,
    )
  ),
  fluidRow(box(width = 12,
    tableOutput(outputId = 'clickedPoints')
  )
  ),
  fluidRow(box(width = 12, #height = 900,
    plotOutput(outputId = "Track")
  )
  )
  )
)

server<-function(input,output,session){

  observeEvent(input$neg10K,{
    print(as.numeric(input$neg10K))
  })
  
  observeEvent(input$pos10K,{
    print(as.numeric(input$pos10K))
  })
  
  # MA-plot of all genes
  output$plotma     <-  renderPlot({
    # par(mar=c(5,5,3,2),cex.lab=2,cex.main=2,cex.axis=1.5)
    
    ggplot(tops,aes(x=AveExpr, y=logFC, color=Significance,ymin = -ymax, ymax = ymax)) +
      geom_point() +
      coord_cartesian() +
      ylab("log2 FC") +
      xlab("log2 Normalized Reads")
    
    # message("Up=",nrow(tops_sig_up))
    # message("Down=",nrow(tops_sig_down))
    # text(xloc,yloc,paste0("Up=",nrow(tops_sig_up)))
    # text(xloc,-yloc,paste0("Down=",nrow(tops_sig_down)))
    # 
    # idx = current$idx
    # 
    # # add circle for the selected point
    # if (!is.null(idx)) points( data2[idx,1], data2[idx,2],
    #                            col="dodgerblue", cex=3, lwd=3 )
  })
  
  #get the clicked points!
  clicked <- reactive({
    # We need to tell it what the x and y variables are:
    nearPoints(tops, input$plotma_click, xvar = "AveExpr", yvar = "logFC")
    # updateTextInput(session, "Gene", value = paste(rownames(clicked()[1,])))
  })
  
  #output those points into a table
  output$clickedPoints <- renderTable({
    if ((nrow(clicked())>0)) clicked()[1,]
  }, rownames = T)
  
  # This will change the value of input$Gene, based on x
  observe({
    if(nrow(clicked())>0){
        updateTextInput(session, "Gene", value = rownames(clicked()[1,]))
    }
  })
  
  # counts plot
  output$plotcounts <-  renderPlot({
    # par(mar=c(5,5,3,2),cex.lab=2,cex.main=2,cex.axis=1.5)
    # plot the counts for the selected gene
    # idx = current$idx
    # if (!is.null(idx)) plotCounts( dds, idx ,"Group")
    tryCatch(
      if ((nrow(clicked())>0)){
        # plotCounts( dds, input$Gene ,"Group")
        ggplot() +
          geom_point(data = vstNormalizedCounts_Macrophage[vstNormalizedCounts_Macrophage$Peaks %in% input$Gene,],
                     aes(x=Group, y=vstNormalizedCounts, color = Group )) +
          scale_color_manual(values= c("blue","red"))
      }
      else{
        textOutput("Please select a region in MA plot")
      },
      error=function(e) e
    )
  })
  
  output$hover_info <- renderPrint({
    if(!is.null(input$plot_hover)){
      hover=input$plot_hover
      dist=sqrt((hover$x-vstNormalizedCounts_Macrophage$Group)^2+(hover$y-vstNormalizedCounts_Macrophage$vstNormalizedCounts)^2)
      cat("SampleID: \n")
      if(min(dist) < 3)
        vstNormalizedCounts_Macrophage$SampleID[which.min(dist)]
    }
  })
  
  output$Track      <-  renderPlot({
    Track_list<-list(
      HC_1 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_1.bw",
      HC_2 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_2.bw",
      HC_3 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_3.bw",
      HC_4 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_4.bw",
      CAD_1 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_1.bw",
      CAD_2 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_2.bw",
      CAD_3 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_3.bw",
      CAD_4 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_4.bw"
    )
    
    Track_cols<-c("blue","blue","blue","blue","red","red","red","red")
    
    url <- "https://atac-tracks-api-e2gbey6dba-uw.a.run.app"
    
    body <- list(
      distflank = input$Flank,
      Genes = input$Gene,
      Track_list = Track_list,
      Track_cols = Track_cols
    )
    
    path <- 'atacTracks'
    
    raw.result %<-% POST(url = url, path = path, body = body, encode = 'json')
    
    apiIn %<-% base::unserialize(httr::content(raw.result))
    
    vp %<-% viewTracks(apiIn$tracks, gr=apiIn$geneRegion, viewerStyle=apiIn$view)
    
  })
}

shinyApp(ui = ui,server = server)
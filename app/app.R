library(shiny)
library(shinydashboard)
library(httr)
library(globals)
library(future)

source("analysis.R")
source("plot_exception.R")

size = 450

ui<-dashboardPage(
  dashboardHeader(title = "CAD Macrophage ATAC-seq"),
  dashboardSidebar(selectInput(inputId = "Comparison",
                               label = "Differential Comparisons",
                               choices = c("HC vs CAD" = "HCvsCAD")),
                   textInput(inputId = "Gene",
                             label = "Gene Symbol",
                             value = "STAT1"),
                   sliderInput(inputId = "Flank",
                               label = "Zoom",
                               min = 0,max = 50000,
                               value = 10000,step = 10000,
                               pre = "x", sep = ",",
                               animate = F)
                   # ,actionButton(inputId = "neg10K", label = "< 10k")
                   # ,actionButton(inputId = "pos10K", label = "> 10k")
  ),
  dashboardBody(
    fluidRow(
    box(title = "MA Plot",
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

  # observeEvent(input$neg10K,{
  #   print(as.numeric(input$neg10K))
  # })
  # 
  # observeEvent(input$pos10K,{
  #   print(as.numeric(input$pos10K))
  # })
  # 
  # MA-plot of all genes
  output$plotma     <-  renderPlot({
    comparison <- input$Comparison
    ymax <- 8
    xloc <- 6 # px
    yloc <- 6 # py
    # tops<-topTable(fitqwr2.cqn, coef = comparison, number = Inf,sort.by = "none")
    # write.table(x = tops,file = "./HCvsCAD.txt",row.names = T,col.names = T,quote = F,sep = '\t')
    tops <- read.table(paste0(File_location,"/",comparison,".txt"),header = T)
    tops$Significance<-"NA"
    tops[tops$adj.P.Val <= 0.05,"Significance"]<-"< 0.05"
    tops[tops$adj.P.Val > 0.05,"Significance"]<-"> 0.05"
    tops_sig<-subset(tops,adj.P.Val<0.05)
    tops_sig_up<-subset(tops_sig,logFC>0)
    tops_sig_down<-subset(tops_sig,logFC < 0)
    tops_not_sig<-subset(tops,adj.P.Val>=0.05)
    PvsU_up<-tops_sig_up
    PvsU_down<-tops_sig_down
    PvsU_not_sig<-tops_not_sig
    tops_sig_mod<-tops_sig_up
    tops_sig_mod<-rbind(tops_sig_mod,tops_sig_down)
    
    ggplot(tops,aes(x=AveExpr, y=logFC, color=Significance,ymin = -ymax, ymax = ymax)) +
      geom_point() +
      coord_cartesian() +
      ylab("log2 FC") +
      xlab("log2 Normalized Reads") +
      ggtitle(gsub("vs"," vs ",input$Comparison))
    
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
    comparison <- input$Comparison
    ymax <- 8
    xloc <- 6 # px
    yloc <- 6 # py
    tops <- read.table(paste0(File_location,"/",comparison,".txt"),header = T)
    tops$Significance<-"NA"
    tops[tops$adj.P.Val <= 0.05,"Significance"]<-"< 0.05"
    tops[tops$adj.P.Val > 0.05,"Significance"]<-"> 0.05"
    
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
          scale_color_manual(values= c("red","blue"))
      }
      else{
        textOutput("Please select a region in MA plot")
      },
      error=function(e) e
    )
  })
  
  output$hover_info <- renderPrint({
    if(!is.null(input$plotcounts)){
      hover=input$plotcounts
      dist=sqrt((hover$x-vstNormalizedCounts_Macrophage$Group)^2+(hover$y-vstNormalizedCounts_Macrophage$vstNormalizedCounts)^2)
      cat("SampleID: \n")
      if(min(dist) < 3)
        vstNormalizedCounts_Macrophage$SampleID[which.min(dist)]
    }
  })
  
  output$Track      <-  renderPlot({
    if((grepl(":", input$Gene)==F) && (input$Gene %notin% UCSC.hg19.genes$V1)){
      return(plot_exception("Please enter a valid Gene Symbol or region."))
    } else {
      # showModal(modalDialog("Loading Tracks", footer=NULL))
      tryCatch({
        Track_list<-list(
          CAD_4 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_4.bw",
          CAD_3 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_3.bw",
          CAD_2 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_2.bw",
          CAD_1 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/CAD_1.bw",
          HC_4 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_4.bw",
          HC_3 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_3.bw",
          HC_2 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_2.bw",
          HC_1 = "http://storage.googleapis.com/gbsc-gcp-lab-jgoronzy_group/Rohit/Tracks/macrophage/HC_1.bw"
        )
        
        Track_cols<-c("red","red","red","red","blue","blue","blue","blue")
        
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
        
        numTracks<-length(Track_list)
        numAll<-length(names(apiIn$tracks))
        j=1
        for(i in ((numAll-numTracks)+1):numAll){
          names(apiIn$tracks)[i]<-gsub("_","-",names(apiIn$tracks)[i])
          names(apiIn$tracks)[i]<-paste0("\t",names(apiIn$tracks)[i])
          setTrackStyleParam(apiIn$tracks[[i]], "ylabpos", "bottomleft")
          setTrackStyleParam(apiIn$tracks[[i]], "ylabgp", list(cex=1, col=Track_cols[j]))
          j=j+1
        }
        
        
        vp %<-% viewTracks(apiIn$tracks, gr=apiIn$geneRegion, viewerStyle=apiIn$view)
      }, warning = function(war) {
        
        # warning handler picks up where error was generated
        message("Tracks could not be loaded!")
        
      }, error = function(err) {
        message("Tracks could not be loaded!")
        
      }
      )
      # removeModal()
    }
  })
}
shinyApp(ui = ui,server = server)
library(shiny)
library(shinydashboard)

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
    distflank<-input$Flank
    if(grepl(":", input$Gene)){
      message("This is region")
      clicked<-paste0(input$Gene,":-")
      message(clicked)
      Genes_toplot_gr<-grFromLocationString(clicked)
    } 
    else {
      message("This is gene")
      Genes<-input$Gene
      Genes_toplot <- Genes
      message(Genes_toplot)
      entrezIDforGenes_toplot <- get(Genes_toplot, org.Hs.egSYMBOL2EG)
      # Genes_toplot <- geneTrack(entrezIDforGenes_toplot,TxDb.Hsapiens.UCSC.hg19.knownGene,type = "gene")[[1]]
      hg19_genes<-genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
      Genes_toplot_gr<-hg19_genes[hg19_genes$gene_id %in% entrezIDforGenes_toplot]
      seqlevels(Genes_toplot_gr) <- seqlevelsInUse(Genes_toplot_gr)
    }

    Genes_toplot_gr<-Genes_toplot_gr + distflank
    
    seqlevelsStyle(Genes_toplot_gr) <- "UCSC"
    
    HC_1 <- importScore(file = HC_1, format="BigWig",ranges = Genes_toplot_gr)
    HC_2 <- importScore(file = HC_2, format="BigWig",ranges = Genes_toplot_gr)
    HC_3 <- importScore(file = HC_3, format="BigWig",ranges = Genes_toplot_gr)
    HC_4 <- importScore(file = HC_4, format="BigWig",ranges = Genes_toplot_gr)
    
    CAD_1 <- importScore(file = CAD_1, format="BigWig",ranges = Genes_toplot_gr)
    CAD_2 <- importScore(file = CAD_2, format="BigWig",ranges = Genes_toplot_gr)
    CAD_3 <- importScore(file = CAD_3, format="BigWig",ranges = Genes_toplot_gr)
    CAD_4 <- importScore(file = CAD_4, format="BigWig",ranges = Genes_toplot_gr)
    
    setTrackStyleParam(HC_1, "color", c("blue","blue"))
    setTrackStyleParam(HC_2, "color", c("blue","blue"))
    setTrackStyleParam(HC_3, "color", c("blue","blue"))
    setTrackStyleParam(HC_4, "color", c("blue","blue"))
    
    setTrackStyleParam(CAD_1, "color", c("red","red"))
    setTrackStyleParam(CAD_2, "color", c("red","red"))
    setTrackStyleParam(CAD_3, "color", c("red","red"))
    setTrackStyleParam(CAD_4, "color", c("red","red"))
    
    y_max<-ceiling(max(c(HC_1$dat$score,HC_2$dat$score,HC_3$dat$score,HC_4$dat$score,
                         CAD_1$dat$score,CAD_2$dat$score,CAD_3$dat$score,CAD_4$dat$score)))
    
    setTrackStyleParam(HC_1, "ylim", c(0,y_max))
    setTrackStyleParam(HC_2, "ylim", c(0,y_max))
    setTrackStyleParam(HC_3, "ylim", c(0,y_max))
    setTrackStyleParam(HC_4, "ylim", c(0,y_max))
    
    setTrackStyleParam(CAD_1, "ylim", c(0,y_max))
    setTrackStyleParam(CAD_2, "ylim", c(0,y_max))
    setTrackStyleParam(CAD_3, "ylim", c(0,y_max))
    setTrackStyleParam(CAD_4, "ylim", c(0,y_max))
    
    # setTrackStyleParam(Genes_toplot, "ylabpos", 'upstream')
    # setTrackStyleParam(Genes_toplot, "color", 'black')
    # setTrackYaxisParam(Genes_toplot, "gp", list(col = "black", lty = "solid", lwd = 3, fontsize = 32))
    # 
    t <- try(Refseq_Genes <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      org.Hs.eg.db,
                                      gr=Genes_toplot_gr))
    if ("try-error" %in% class(t)){
      optSty <- optimizeStyle(trackList(HC_1, HC_2, HC_3, HC_4, CAD_1, CAD_2, CAD_3, CAD_4), theme=NULL)
      trackList <- optSty$tracks
      viewerStyle <- optSty$style
      vp <- viewTracks(trackList, gr=Genes_toplot_gr, viewerStyle=viewerStyle)
    }
    else {
      Refseq_Genes_names<-c()
    for (i in 1:length(Refseq_Genes)){
      setTrackStyleParam(Refseq_Genes[[i]], "ylabpos", "upstream")
      setTrackStyleParam(Refseq_Genes[[i]], "ylabgp", list(cex=.6))
      setTrackStyleParam(Refseq_Genes[[i]], "color", 'black')
      setTrackYaxisParam(Refseq_Genes[[i]], "gp", list(col = "black", lty = "solid", lwd = 3, fontsize = 32))
      Refseq_Genes_names<-c(Refseq_Genes_names,paste0(Refseq_Genes[[i]]$dat$symbol,"::",Refseq_Genes[[i]]$dat$transcript)[1])
    }
    names(Refseq_Genes)<-Refseq_Genes_names

    # eval(parse(text=(paste("obj<-list(",Genes,"=Genes_toplot,Transcripts=Refseq_Genes)",sep=""))))
    
    eval(parse(text=(paste("obj<-list(","Transcripts=Refseq_Genes)",sep=""))))
    
    # eval(parse(text=(paste("obj<-list(",Genes,"=Genes_toplot)",sep=""))))
    
    optSty <- optimizeStyle(trackList(obj[1], HC_1, HC_2, HC_3, HC_4, CAD_1, CAD_2, CAD_3, CAD_4), theme=NULL)
    trackList <- optSty$tracks
    viewerStyle <- optSty$style
    vp <- viewTracks(trackList, gr=Genes_toplot_gr, viewerStyle=viewerStyle)
    }
  })
}

shinyApp(ui = ui,server = server)
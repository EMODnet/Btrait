
# ----------------------------------------------
# ----------------------------------------------
# SHINY server for the macroalgal model
# ----------------------------------------------
# ----------------------------------------------

# for debugging only
INPUT       <<- NULL  
DD <<- NULL
#  INPUT <<-  reactiveValuesToList(input) # TO SEE VALUES OF input, as a list

shinyServer(function(input, output, session) {

#===========================================
# Create lists on the fly
#===========================================
  
#---------------------------  
# taxa 
#---------------------------  
  output$SettaxonList <- renderUI({
    DAT      <- getDataset()
    taxLevel <- getTaxLevel()
    label    <- paste("Select", taxLevel)
    
    if (taxLevel == "taxa"){
      if (input$all)
        LIST <- sort(unique(DAT$density$taxon))
      else
        LIST <- sort(names(sort(table(DAT$density$taxon), 
                     decreasing=TRUE))[1:30])
    } else if (taxLevel == "abiotics"){ 
       LIST <- colnames(DAT$abiotics)[-1] 
       LIST <- LIST[!LIST=="EUNIScode"]
    } else if (taxLevel == "traits"){
       meta <- metadata(Traits_nioz)
       LIST <- unique(meta$trait)
    } else if (taxLevel == "groups"){
       meta <- metadata(Groups)
       LIST <- unique(meta$description)
    } else if (taxLevel == "summary")
       LIST <- c("density", "taxa", "bioturbation", "bioirrigation")
    
    isolate({
      selectInput(inputId="taxon", label=label, choices=LIST, selectize=FALSE)
    })   
  })

#---------------------------  
# abiotic
#---------------------------  
  output$SetabioticsList <- renderUI({
    DAT      <- getDataset()
    label    <- "Select abiotic x-variable"
    
    LIST <- colnames(MWTL$abiotics)[-(1:2)]
    LIST <- LIST[!LIST=="EUNIScode"]

    isolate({
      selectInput(inputId="abioticVar", label=label, choices=LIST, selectize=FALSE, 
        selected = "mud")
    })   
  })

#===========================================
# Get data
#===========================================

  getDataset <- reactive({
    DAT      <- MWTL
    DAT$main <- "MWTL"
    DAT
  })

  getAbiotics <- reactive({
    input$abioticVar
  })

#---------------------------  
  getTraitData <- reactive({
#    tr <- input$settrait
#    if (tr =="Nioz")               Dat <- Traits_nioz
#    else if (tr == "Cefas")        Dat <- Traits_cefas 
#    else if (tr == "Bioturbation") Dat <- Traits_Db
#    else if (tr == "Irrigation")   Dat <- Traits_irr
    Dat <- Traits_nioz
    Dat
  })

#---------------------------  
  getWhat <- reactive({
    what     <- input$what
    what
  })
  
#---------------------------  
  getTaxLevel <- reactive({
    Var     <- input$settaxon
#    Var <- "species"
    Var
  })
  
#---------------------------  
  getSpData <- reactive({  # get the taxon data that are selected
    INPUT    <<-  reactiveValuesToList(input)
    DAT      <- getDataset()
    what     <- getWhat()
    taxLevel <- getTaxLevel()

    getSpeciesData(DATA    = DAT, 
                   taxon   = input$taxon, 
                   settaxon= taxLevel, 
                   what    = what)
  })

#---------------------------  
  getSpYear <- reactive({  # keep year
    INPUT    <<-  reactiveValuesToList(input)
    DAT      <- getDataset()
    what     <- getWhat()

    getSpeciesYear(DATA=DAT, taxon=input$taxon, 
                   settaxon=getTaxLevel(), what=what)  
  })

#---------------------------  
  getUnits <- reactive({  # get the label of the y-value 
    DAT <- getDataset()
    taxLevel <- getTaxLevel()

    if (taxLevel == "traits"){
        units  <- subset(attributes(Traits_nioz)$description, 
                         trait==input$taxon) [1,"units"]
    } else if (taxLevel == "abiotics"){
        units  <- subset(attributes(DAT$abiotics)$description, 
                         name==input$taxon) [1,"units"]
    } else  
        units <- subset(attributes(DAT$density)$description, 
                        name==getWhat() )$units
  })
  
#---------------------------  
  getlog <- reactive({
      LA <- input$logaxes
      
      Log <- ""
      if      (LA == "x-axis") Log <- "x"
      else if (LA == "y-axis") Log <- "y"
      else if (LA == "both")   Log <- "xy" 
      Log
    })

#---------------------------  
  getylog <- reactive({
      LA <- input$logaxes
      
      Log <- ""
      if  (LA %in% c("y-axis", "both"))  Log <- "y" 
      Log
    })
#---------------------------  
  getclog <- reactive({
    
      Log <- ""
      if  (input$clog) Log <- "c"
      Log
    })

#===========================================
# plotting and printing
#===========================================

  output$printTraits <- renderPrint({
    taxLevel <- getTaxLevel()

    if (taxLevel == "taxa")
      getTraitVals(Trait=getTraitData(), taxon=input$taxon)
    else if (taxLevel == "traits")
      getTraitSpecs(tr=input$taxon)
    else if (taxLevel == "groups")
      getGroupSpecs(tr=input$taxon)
    else if (taxLevel == "summary"){
     try(Data  <- getSpYear(), silent=TRUE)
      DD <- Data
        if (! is.null(Data))  summary(Data[,5])
    }
    else if (taxLevel == "abiotics")
      subset(metadata(getDataset()$abiotics), name==input$taxon)
  })
  
  output$printTaxonomy <- renderPrint({
    if  (input$settaxon == "taxa")
      getTaxon(taxon=input$taxon)
    else{
     cat(" ")
    }
  })
  
  getColorScheme <- reactive({
    if (input$color == "green/blue") 
        Col <- ramp.col(c("lightgreen","darkblue"),100) 
    else
        Col <- jet2.col(100)
    Col
  })

#--------------------------
# species/trait maps
#--------------------------

  plotmap <- reactive({
        
     DAT  <- getDataset()
     Log  <- getclog()
     Data <- NULL   # to prevent an error
     Col  <- getColorScheme()
     try(Data  <- getSpData(), silent=TRUE)

     #     Data <- na.omit(Data)
     if (! is.null(Data)){
       Main  <- input$taxon
       taxLevel <- getTaxLevel()

       if (taxLevel == "traits"){
         Clab <- c("",getUnits())
         Main <- paste(Main, " (mean based on ", getWhat(),")", sep="")
       
       } else if (taxLevel == "abiotics"){
         Clab <- c("", getUnits())
       
       } else  
         Clab <- c(getWhat(), getUnits())
       
       par(las=1, mar=c(3,3,3,0))
       nc <- ncol(Data)
       cex <- 2
       
       mapBtrait(contours=DAT$contour, x=Data$x, y=Data$y, 
                 colvar=Data[,nc], log=Log, pch=15, cex=cex, 
                 main=Main, clab=Clab, draw.levels=TRUE, col=Col,
                 colkey=list(length=0.2, width=0.5, 
                   cex.axis=0.8, shift=-0.25, dist=-0.25,         
                   cex.clab=par("cex.lab"))) 
       
        NSnot <- DAT$stations[!DAT$stations$station %in%Data$station,]
        
        points2D(NSnot$x, NSnot$y, colvar=NULL, col="darkgrey", 
                 pch="+", add=TRUE, cex=0.75)
     } 

  })  
  
#--------------------------

  output$spMap <- renderPlot({
    plotmap()
  })

#--------------------------
# give mean value (not used)
#--------------------------
  output$meanValue <- renderText({
    Data  <- getSpData()
    nc    <- ncol(Data)
    Mean  <- mean(Data[,nc])
    prettyNum(Mean)
  })

#--------------------------
# species/trait vs abiotic
#--------------------------
  output$xyPlot <- renderPlot({
     Log   <- getlog()
     Data  <- NULL   # to prevent an error
     DAT   <- getDataset()
     try(Data  <- getSpData(), silent=TRUE)
     
     if (! is.null(Data)){
       nc <- ncol(Data)
       colnames(Data)[nc] <- "value"
           taxLevel <- getTaxLevel()
      Main  <- paste(input$taxon)
      if (taxLevel == "taxa")
        Main  <- paste(Main, getWhat(), sep=", ") 
       taxLevel <- getTaxLevel()

       sabi  <- getAbiotics()
       Datx  <- merge(DAT$station,DAT$abiotics[c("station", "depth", sabi)] )    
       data  <- merge(Datx, Data)
       
       if (taxLevel != "abiotics")
         ylab <- paste(" [",getUnits(),"]", sep="")
       else
         ylab <- paste("[", getUnits(),"]", sep="")
       
       xlab2  <- paste(sabi, ",", subset(attributes(DAT$abiotics)$description, 
                         name==sabi) [1,"units"])

       par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(4,4,2,2))      
       plot(data$depth, data$value, main="", xlab="depth, m", ylab=ylab, 
          log=Log, pch=16, col="lightblue", cex=2, las=1)
       plot(data[,sabi], data$value, main="", xlab=xlab2, ylab=ylab, 
          log=Log, pch=16, col="lightblue", cex=2, las=1)
       mtext(side=3, outer=TRUE, line=-1.5, Main, cex=1.5)
     }
  }) 

#--------------------------
# time series
#--------------------------
  output$tsPlot <- renderPlot({
     Log   <- getylog()
     Data  <- NULL   # to prevent an error
     try(Data  <- getSpYear(), silent=TRUE)
     DD <<- Data
     
     if (! is.null(Data)){
       taxLevel <- getTaxLevel()
       
       nc <- ncol(Data)-1
       Main  <- paste(input$taxon, ",", getWhat() )
       if (taxLevel != "abiotics")
         ylab <- paste(" [", getUnits(),"]", sep="")
       else
         ylab <- paste("[", getUnits(),"]", sep="")

       if (nrow(Data) == 0) 
         plotNone(main=Main, xlab="year", 
            ylab=ylab)
       else
         boxplot(Data[,nc]~ Data$year, main=Main, xlab="year", 
           ylab=ylab, log=Log, col="lightblue", las=1, 
           at=sort(unique(Data$year)))
     }
  }) 

})


# Load functions
source('www/xGlobalFunctions.R')[[1]]
source('www/uploadFunctions.R')[[1]]
source('www/processFunctions.R')[[1]]
source('www/fitFunctions.R')[[1]]
source('www/logFunctions.R')[[1]]

# Define model types and default variables
.GlobalEnv$modelTypes <- c("MW / Single ideal species",
                           "Kd / Monomer-Nmer")
.GlobalEnv$fitPlotWidth <- 290
.GlobalEnv$fitPlotHeight <- 300
defineVars()

# Define server logic required to draw a histogram
function(input, output, session) {

  # ------------ STARTUP

  
  # Load upload tab
  output$upload <- renderUI({
    source('www/upload.R')[[1]]
  })  
  
  
  # Load fit tab
  output$fit <- renderUI({
    source('www/fit.R')[[1]]
  })  

  
  # Render process tab
  output$process <- renderUI({
    source('www/process.R')[[1]]
  })
  
  
  # Reset variables on refresh
  observe({
    # Load files
    defineVars()
    resetUI(output)
    updateLog(output,'UltraAnalysis loaded')
  })
  
  
  # ------------ ON TAB SWITCH
  
  
  # Observe a tab switch
  observeEvent(input$tabSwitch,{
    
    if(input$tabSwitch=='upload'){
      renderPsvInput(output)
      renderSdInput(output)
    }
    
    if(input$tabSwitch=='process'){
      renderScanSelection(output)
      renderScanPlot(input,output)
      renderSavedSelections(output)
    }
    
    if(input$tabSwitch=='fit'){
      modelTypeToShow <- ifelse(length(dataList$selectedModel)==0,modelTypes[1],dataList$selectedModel)
      renderModelTypeSelection(output,modelTypeToShow)
    }
    
  })
  
  
  # ------------ UPLOAD TAB EVENTS

  
  # When user uploads data
  observeEvent(input$uploadInput,{
    
    stop <- ''
    
    # Reset vars
    defineVars()
    updateLog(output,'UltraAnalysis loaded')
    
    # Reset ui
    resetUI(output)
    
    # Get uploaded file
    uploadInput <- input$uploadInput
    file <- uploadInput$datapath
    filenames <- uploadInput$name
    
    # Update log
    updateLog(output,paste("Uploaded ",length(file)," files",sep=""))
    
    # Create progress bar
    progress <- shiny::Progress$new()
    progress$set(message='Loading uploaded file...',value=0.1)
    
    if(str_detect(file[1],'.Rdata')){
      
      if(length(file)>1){
        
        stop <- "You cannot upload multiple .Rdata files at once"
        
      } else{
        
        # Load to global environment
        load(file,envir=.GlobalEnv)
        
        # Generate hooks for saved sectors
        savedSectors <- dataList$savedSectors
        lapply(savedSectors,function(x){
          # Create hook for editing 
          generateSectorHooks(input,output,session,x)
        })
        
        # Render UIs
        defaultScansToAnalyze <- dataList$scansToAnalyze
        
      }
      
    } else{
      # Read each scan
      progress$set(message='Reading scan data...',value=0.5)
      scanData <- readScans(file,filenames)
      
      # Update the datalist with scan data
      updateDataList('scanData',scanData)
      defaultScansToAnalyze <- NULL
    }
    
    if(stop==''){
      # Format scan plot output
      renderScanPlots(input,output,defaultScansToAnalyze)
      
      # Render global inputs
      renderPsvInput(output)
      renderSdInput(output)
      
      # Close progress bar
      progress$close()
      
      #defaultScansToAnalyze <- NULL
    } else{
      updateLog(output,stop)
    }
  })
  

  # ------------ PROCESS TAB EVENTS
  

  # Observe scaling y axis
  observeEvent(input$scaleSelectionPlot,{
    renderScanPlot(input,output)
  })
  
  
  # When user selects a scan to see up close
  observeEvent(input$selectScan,{
    
    selectedScan <- as.numeric(input$selectScan)
    renderPlot <- ifelse(length(dataList$skipPlotRerender)>0,dataList$skipPlotRerender,FALSE)
    
    if(!is.na(selectedScan)&&!renderPlot){
      # Update data list
      updateDataList('selectedScan',selectedScan)
      # Plot
      resetBrush(output,session)
      renderScanPlot(input,output)
    } else{
      updateDataList('skipPlotRerender',FALSE)
    }
    
  })
  
  
  # When user brushes an area on the plot (selects a sector)
  observeEvent(input$selectingSector,{
    
    currentScanData <- retrieveScanData()
    selectedSector <- input$selectingSector
    
    # If selectedSector is alreayd equal to the selected sector, then return instead of running 2x
    if(length(dataList$selectedSector)>0){
      if(selectedSector$xmax==dataList$selectedSector$xmax&&selectedSector$xmin==dataList$selectedSector$xmin){
        return()
      }
    }
    
    # Update datalist
    updateDataList('selectedSector',selectedSector)
    
    # Update plots
    renderSelectedSectorPlots(input,output,currentScanData,selectedSector)
    
    # Stop skipping rerender
    updateDataList('skipPlotRerender',FALSE)
    
  })
  
  
  # When a user clicks the plot, reset to normal view
  observeEvent(input$resetPlotView,{
    resetBrush(output,session)
    renderScanPlot(input,output)
  })
  
  
  # When user saves a sector
  observeEvent(input$saveSector,{
    
    # Get receptor concentration
    receptorConcentration <- tryCatch(ifelse(length(input$receptorInput)==0,0,
                                             ifelse(input$receptorInput=='',0,input$receptorInput)
                                             ),error=function(e)return(0))
    
    # Get scan data
    currentScanData <- retrieveScanData()
    
    # Get sector data
    selectedSector <- dataList$selectedSector
    
    # Get current sector df
    currentSectorData <- subset(currentScanData,r>selectedSector$xmin&r<selectedSector$xmax)
    
    # Get current sector df
    currentSectorList <- createSector(dataList$selectedScan,currentSectorData,receptor=receptorConcentration)
    currentSectorDf <- currentSectorList$sector
    sectorID <- currentSectorList$id
    
    # Get saved sectors
    savedSectors <- dataList$savedSectors
    
    # Check if sector ID already exists
    otherIDs <- unlist(lapply(savedSectors,'[[',6))
    
    if(!(sectorID%in%otherIDs)){
      currentSectorDf$ID <- sectorID
      
      # Define sector color
      currentSectorDf$Color <- getRandomColor()
      
      newSectorLog <- paste("Added sector: Scan ",currentSectorDf$Scan," Min ",currentSectorDf$Min,' Max ',currentSectorDf$Max, " n = ",currentSectorDf$n,sep="")
      tryCatch(updateLog(output,newSectorLog),error=function(e)return())
      
      # Save as list
      savedSectors[[length(savedSectors)+1]] <- currentSectorDf
      
      # Update dataList
      updateDataList('savedSectors',savedSectors)
      
      # Create hook for editing 
      generateSectorHooks(input,output,session,currentSectorDf)
      
      # Update UI
      renderSavedSelections(output)
      
      # Update plot
      resetBrush(output,session)
      renderScanPlot(input,output)
    }
    
    
  })

  
  # For auto sector finding
  observeEvent(input$autoSector,{
    
    if(length(dataList$scansToAnalyze)==0){
      return()
    }
    
    if(length(dataList$scanData)==0){
      return()
    }
    
    # Create progress bar
    .GlobalEnv$progress <- shiny::Progress$new()
    progress$set(message='Loading data...',value=0)
    
    # For each scan, run find miniscus
    allScanData <- dataList$scanData[dataList$scansToAnalyze]
   
    # Find menisci
    identifiedSectors <- lapply(1:length(allScanData),findMenisci,input,output,session,allScanData) %>% bind_rows()
    
    # Covnert each row to a list
    identifiedSectorsList <- lapply(1:nrow(identifiedSectors),function(x){return(identifiedSectors[x,])})
    
    # Update saved ranges
    # Save as list
    savedSectors <- dataList$savedSectors
    newSavedSectors <- c(savedSectors,identifiedSectorsList)
    
    # Update dataList
    updateDataList('savedSectors',newSavedSectors)
    
    # Update UI
    renderSavedSelections(output)
    
    # Update plot
    resetBrush(output,session)
    renderScanPlot(input,output)
    
    progress$close()
    
  })
  
  
  # ------------ FIT TAB EVENTS
  
  
  # When user selects a model type
  observeEvent(input$modelType,{
  
    # Render UI for model parms
    selectedModel <- input$modelType
    updateDataList('selectedModel',selectedModel)
    
    output$modelParms <- renderUI({
      
      if(selectedModel=='MW / Single ideal species'){
        div(class='column',
          actionButton('fitData','🏃 Run fit',style='align-self:center')
        )
      } else if(selectedModel=='Kd / Monomer-Nmer'){
        div(class='column',
            textInput('mwInput',NULL,dataList$mwValue,300,'MW (Da)'),
            textInput('eCoefInput',NULL,dataList$eCoefValue,300,'Extinction coef. (M^-1cm^-1)'),
            textInput('nInput',NULL,dataList$nValue,300,'N (mers)'),
            actionButton('fitData','🏃 Run fit',style='align-self:center')
          )
      }
      
    })
    
    output$fitDescription <- renderText({
      paste("Global fit for: ",dataList$selectedModel,sep="")
    })
    
  })
  
  
  # When user tries to fit the data
  observeEvent(input$fitData,{
    
    checkInputs <- tryCatch(checkForInputs(input,output),error=function(e)return(FALSE))
  
    if(checkInputs==FALSE){
      return()
    }
    
    .GlobalEnv$progress <- shiny::Progress$new()
    progress$set(message='Fitting data...',value=0.5)
    
    
    # Get saved sectors
    savedSectors <- dataList$savedSectors 
    
    # Get all absolute scans
    allScans <- unique(lapply(dataList$scanData,'[[','AbsoluteScan') %>% unlist())
    
    if(length(dataList$selectedModel)==0){
      updateDataList('selectedModel','MW / Single ideal species')
    }
    
    # Define parameters of fit
      # Returns fit function among other things depending on selected model to fit
    .GlobalEnv$fitParms <- defineFitParms(input,output)
    
    # Processes sector data per fit function, adding parms accordingly
    processedSectors <- lapply(1:length(savedSectors),
                               calculateReferenceRadius,
                               savedSectors=savedSectors,
                               allScans=allScans,
                               fitParms=fitParms)
    
    if(length(which(is.na(processedSectors)))>0){
      print("Error in fitting baseline")
      tryCatch(updateLog(output,"Error fitting baselines"),error=function(e)return())
      return()
    }
    
    # Combine into 1 global dataset
    globalData <- processedSectors %>% bind_rows()

    # Dynamically render fit function from parms
    processedFitFunction <- processFitFunction(fitParms,globalData)
   
    # Do fit
    startTime <- Sys.time()
    globalFit <- eval(parse(text=processedFitFunction$fitString))
    endTime <- Sys.time()
    
    # Report processing time
    timeDiff <- as.numeric(difftime(endTime,startTime,units='secs'))
    
    timeString <- paste("Processed ",globalFit$convInfo$finIter," iterations fitting ",length(globalFit$m$getPars()),' parameters in ',round(timeDiff,1),' seconds with ',globalFit$convInfo$trsName,'. Stop message: ',globalFit$convInfo$stopMessage,sep="")
    tryCatch(updateLog(output,timeString),error=function(e)return())
    
    # Process fit results for plotting
    fit <- processFitResults(globalFit,globalData,processedFitFunction,fitParms)
    dev.off()
    
    # Render plots as outputs
    output$nonlinearPlot <- renderPlot(fit$nonlinear,width=fitPlotWidth,height=fitPlotHeight)
    output$residualPlot <- renderPlot(fit$residuals,width=fitPlotWidth,height=fitPlotHeight)
    #output$linearPlot <- renderPlot(fit$linear,width=320,height=300)
    
    # Save fit to global env
    .GlobalEnv$fit <- fit
    
    # Render ui
    renderFitOutput(output)
    progress$close()
  })
  
  
  # When user clicks residual plot
  observeEvent(input$clickFitResultNonlinearResiduals,{
    
    .GlobalEnv$userClick <- input$clickFitResultNonlinearResiduals
    
    # Get residual plot data
    residualPlot <- fit$residuals
    residualPlotData <- residualPlot$data
    
    # Get x and y
    x <- residualPlotData$r-residualPlotData$r0
    y <- residualPlotData$Ar-residualPlotData$pAr
    
    # Get errs
    errs <- sqrt((userClick$x-x)^2+(userClick$y-y)^2)
    nearestPoint <- residualPlotData[which(errs==min(errs,na.rm=TRUE)),]
    
    # Get sector data
    sectorData <- subset(residualPlotData,ID==nearestPoint$ID)
    
    # Rerender plot with nearest point
    newPlot <- residualPlot+
      geom_point(data=nearestPoint,aes(x=r-r0,
                                       y=Ar-pAr),
                 col='blue',size=4)+
      geom_label(data=nearestPoint,aes(x=-Inf,
                                       y=Inf),
                 label=paste("Scan ",
                             nearestPoint$AbsoluteScan,
                             " Sector ",paste(round(min(sectorData$r),3),' - ',round(max(sectorData$r),3),sep=""),
                             "\n(",round(nearestPoint$r,3),
                             ", ",round(nearestPoint$Ar-nearestPoint$offset,3),")",
                             sep=""),col='black',hjust=0,vjust=1,fontface='bold')
    
    output$residualPlot <- renderPlot(newPlot,width=fitPlotWidth,height=fitPlotHeight)
    
  })
  
  
  # When user clicks nonlinear plot
  observeEvent(input$clickFitResultNonlinear,{
    
    .GlobalEnv$userClick <- input$clickFitResultNonlinear
    
    # Get nonlinear plot data
    nonlinearPlot <- fit$nonlinear
    nonlinearPlotData <- nonlinearPlot$data
    
    # Get x and y
    x <- nonlinearPlotData$r
    y <- nonlinearPlotData$Ar-nonlinearPlotData$offset
    
    # Get errs
    errs <- sqrt((userClick$x-x)^2+(userClick$y-y)^2)
    nearestPoint <- nonlinearPlotData[which(errs==min(errs,na.rm=TRUE)),]
    
    # Get sector data
    sectorData <- subset(nonlinearPlotData,ID==nearestPoint$ID)
    
    # Rerender plot with nearest point
    newPlot <- nonlinearPlot+
      geom_point(data=nearestPoint,aes(x=r,
                                       y=Ar-offset),
                 col='blue',size=4)+
      geom_label(data=nearestPoint,aes(x=-Inf,
                                      y=Inf),
                label=paste("Scan ",
                            nearestPoint$AbsoluteScan,
                            " Sector ",paste(round(min(sectorData$r),3),' - ',round(max(sectorData$r),3),sep=""),
                            "\n(",round(nearestPoint$r,3),
                            ", ",round(nearestPoint$Ar-nearestPoint$offset,3),")",
                            sep=""),col='black',hjust=0,vjust=1,fontface='bold')
    
    output$nonlinearPlot <- renderPlot(newPlot,width=fitPlotWidth,height=fitPlotHeight)
    
  })
  
  
  # Reset plots
  observeEvent(input$resetFitPlots,{
    output$nonlinearPlot <- renderPlot(fit$nonlinear,width=fitPlotWidth,height=fitPlotHeight)
    output$residualPlot <- renderPlot(fit$residuals,width=fitPlotWidth,height=fitPlotHeight)
    #output$linearPlot <- renderPlot(fit$linear,width=325,height=300)
  })
  
  
}

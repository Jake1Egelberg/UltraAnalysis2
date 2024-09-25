
# Define fitting function
.GlobalEnv$singleIdealSpecies <- function(r, w, R, temp, r0, A0, Mb, offset){
  Ar <- offset + A0 * exp(1)^( ( (Mb*w^2)/(2*R*temp) )*(r^2 - r0^2) )
  return(Ar)
}

# file <- c(
#   "C:/1_Documents/Ferguson Lab/UltraAnalysis2/Data/00027.RA4",
#   "C:/1_Documents/Ferguson Lab/UltraAnalysis2/Data/00026.RA4",
#   "C:/1_Documents/Ferguson Lab/UltraAnalysis2/Data/00025.RA4",
#   "C:/1_Documents/Ferguson Lab/UltraAnalysis2/Data/00024.RA4"
# )
# filenames <-  c(
#   "00027.RA4",
#   "00026.RA4",
#   "00025.RA4",
#   "00024.RA4"
# )

# Read each scan
.GlobalEnv$readScans <- function(file,filenames){
  
  # Read each scan
  readAScan <- function(x,file,filenames){
    
    row <- file[x]
    filename <- filenames[x]
    
    # Read data
    scanLines <- read.table(row,header=TRUE,sep="\t")
    
    # Get scan name
    scanName <- names(scanLines)
    
    # Extract metadata
    rawMeta <- scanLines[1,]
    metaVec <- unlist(strsplit(rawMeta," "))
    
    cellNumber <- as.numeric(metaVec[2])
    temp <- as.numeric(metaVec[3])+273.15
    rpm <- as.numeric(metaVec[4])
    w <- (rpm/60)*2*pi
    
    # Get only values
    valueLines <- scanLines[-1,]
    
    # Split values
    splitValues <- strsplit(valueLines," ")
    
    # Remove blanks from split values
    # IMPORTANT, negative values and positive values have differnt indices before blanks removed
    noBlanks <- lapply(splitValues,function(x){
      blankInds <- which(x=="")
      if(length(blankInds)>0){
        return(x[-blankInds])
      } else{
        return(x)
      }
    })
    
    # Extract values
    radialPosition <- unlist(lapply(noBlanks,'[[',1))
    absorbance <- unlist(lapply(noBlanks,'[[',2))
    noise <- unlist(lapply(noBlanks,'[[',3))
    
    R <- 8.3144 * 1000 * 100 * 100 # (g * cm^2) / s^2 * mol * K
    
    # Get scan number
    scanNum <- as.numeric(str_remove_all(str_sub(filename,end=str_locate(filename,'\\.')[[1]]-1),'[:alpha:]|[:punct:]'))
    
    # Create df
    dfValues <- data.frame(
      Name = NA,
      Cell = cellNumber,
      Speed = rpm,
      Scan = x,
      AbsoluteScan = scanNum,
      CellScanSpeed = NA,
      r = radialPosition,
      Ar = absorbance,
      r0=NA,
      w = w,
      temp = temp,
      R = R,
      theta = w^2 / (2*R*temp)
    ) %>% mutate_all(as.numeric)
    dfValues$Name <- scanName
    dfValues$File <- file[x]
    dfValues$CellScanSpeed <- paste("Cell ",cellNumber," Scan ",scanNum," ",rpm," rpm",sep="")
    
    naInds <- which(is.na(dfValues$Ar))
    if(length(naInds)>0){
      curatedValues <- dfValues[-naInds,]
    } else{
      curatedValues <- dfValues
    }
    
    return(curatedValues)
    
  }
  out <- lapply(1:length(file),readAScan,file=file,filenames=filenames)
  
  # Get scan names
  scanNames <- lapply(out,'[[','AbsoluteScan')%>%unlist()%>%unique()
  names(out) <- scanNames

  return(out)
  
}

# function to get random character
.GlobalEnv$getRandomCharacter <- function(){
  id <- paste('r',sample(1:10000)[1],sep="")
  return(id)
}

# Function to get random number
.GlobalEnv$getRandomNumber <- function(){
  id <- sample(1:100000)[1]
  return(id)
}

# Function to get random color
.GlobalEnv$getRandomColor <- function(){
  col <- sample(colorVec,1)
  return(col)
}

# Summarize a nls as a string
.GlobalEnv$summarizeNLS <- function(fit,type){
  
  # Get fit summary
  fitSummary <- summary(fit)
  
  # Get coef
  fitCoefs <- as.data.frame(coef(fitSummary))
  
  # Format summary string
  parmsFit <- fitSummary$df[1]
  coefInfoSpec <- 2
  coefInfo <- paste(paste(rownames(fitCoefs),": ",round(fitCoefs$Estimate,coefInfoSpec)," +/- ",round(fitCoefs$`Std. Error`,coefInfoSpec),sep=""),collapse=", ")
  itToConverge <- fitSummary$convInfo$finIter
  residualSE <- round(fitSummary$sigma,5)
  
  # Format summary string
  summaryString <- paste(
    "Fitting for: ",type," | ",itToConverge," iterations to convergence | Residual SE: ",residualSE,
    " | ",parmsFit," parameters fit: ",coefInfo,sep=""
    
  )
  
  # Save coefficients and summary string
  out <- list()
  out$summary <- summaryString
  out$coefs <- fitCoefs
  return(out)
  
}

# Define colros to use
.GlobalEnv$colorVec <- rainbow(10)

# Define model types
.GlobalEnv$modelTypes <- c("MW / Single ideal species",
                           "Kd / Monomer-Nmer")
psv <- 0.71
sd <- 1.003

# Define starting vars
defineVars <- function(){
  #.GlobalEnv$logText <- NULL
  #.GlobalEnv$scanData <- NULL #testData[[1]]
  #.GlobalEnv$savedSectors <- NULL #list(data.frame(Scan=1,Min=5.887,Max=6.046,n=68,Receptor=0,ID='r1726670506',Color='#00FFFF'))
  #.GlobalEnv$selectedModel <- modelTypes[1]
  
  # Define data list
    # selectedModel - Type of model to fit
    # logText - text of log
    # scanData - all loaded scans
    # scansToAnalyze - scans selected for inclusion in analysis
    # selectedScan - scan selected for preview
    # selectedSector - currently brushed sector on the plot
    # savedSectors - sectors saved by the user for fitting
    # psvValue - partial specific volume
    # sdValue - solvent density
  dataList <- list()
  # dataList$selectedModel <- modelTypes[1]
  # dataList$selectedScan <- 1
  # dataList$savedSectors <- list(data.frame(Scan=1,Min=5.887,Max=6.046,n=68,Receptor=0,ID='r1726670506',Color='#00FFFF'))
  # dataList$scanData <- testData[[1]]
  # dataList$scansToAnalyze <- 1
  .GlobalEnv$dataList <- dataList
  
}
defineVars()

testplot <- ggplot(data.frame(x=1:10,y=1:10),aes(x=x,y=y))+
  geom_point()

# function to update data list
.GlobalEnv$updateDataList <- function(dataType,newValue){
  
  dataList[[dataType]] <- newValue
  .GlobalEnv$dataList <- dataList
  
  
}


# Define server logic required to draw a histogram
function(input, output, session) {

  # ------------ STARTUP

  # Update the log text
  updateLog <- function(update){
    
    logText <- c(update,dataList$logText)
    updateDataList('logText',logText)
    
    output$log <- renderUI({
      source('www/log.R')[[1]]
    })
    
  }
  
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
  
  # Reset the ui
  resetUI <- function(){
    
    # Process
    output$savedSectorUI <- NULL
    output$sectorPlotFrame <- NULL
    output$currentScanPlotOutput <- NULL
    
    # Fit
    output$fitResult <- NULL
    
    
  }
  
  # Reset variables on refresh
  observe({
    defineVars()
    resetUI()
    updateLog('UltraAnalysis loaded')
  })
  
  # ------------ ON TAB SWITCH
  
  
  # Observe a tab switch
  observeEvent(input$tabSwitch,{
    
    if(input$tabSwitch=='upload'){
      renderPsvInput()
      renderSdInput()
    }
    
    if(input$tabSwitch=='process'){
      renderScanSelection()
      renderScanPlot()
      renderSavedSelections()
    }
    
    if(input$tabSwitch=='fit'){
      modelTypeToShow <- ifelse(length(dataList$selectedModel)==0,modelTypes[1],dataList$selectedModel)
      renderModelTypeSelection(modelTypeToShow)
    }
    
  })
  
  
  # ------------ UPLOAD
  
  # Function when someone checks a scan to include
  includeScan <- function(scanDataInd,scanNumber){
    
    observeEvent(eval(parse(text=paste('input$scan',scanDataInd,sep=""))),{
      
      # Retrieve value and scan
      value <- eval(parse(text=paste('input$scan',scanDataInd,sep="")))
      
      # Update scans to analyze
      if(value){
        
        if(!(scanNumber%in%dataList$scansToAnalyze)){
          scansToAnalyze <- c(scanNumber,dataList$scansToAnalyze)
          scansToAnalyze <- scansToAnalyze[order(scansToAnalyze)]
          
        } else{
          scansToAnalyze <- dataList$scansToAnalyze
        }
        
      } else{
        if(scanNumber%in%dataList$scansToAnalyze){
          scansToAnalyze <- dataList$scansToAnalyze[-which(dataList$scansToAnalyze==scanNumber)] 
        } else{
          scansToAnalyze <- dataList$scansToAnalyze
        }
      }
      
      updateDataList('scansToAnalyze',scansToAnalyze)
      
      # Get absolute scan numbers for scans to analyze
      absoluteScanNumbers <- lapply(dataList$scanData[dataList$scansToAnalyze],'[[','AbsoluteScan')%>%unlist()%>%unique()
      updateDataList('absoluteScanNumbers',absoluteScanNumbers)
      
    })
    
  }
  
  # Render preview of scan data within scanPlotPreviews
  generateScanPlotPreview <- function(currentScan,scanDataInd,defaultScansToAnalyze){
    
    scanPlot <- ggplot(currentScan,aes(x=r,y=Ar))+
      geom_line()+
      xlab("r (cm)")+
      geom_hline(yintercept=1,col='gray50',lty='dashed')+
      geom_hline(yintercept=0,col='gray50',lty='dashed')+
      theme_prism()
    
    absoluteScanNumber <- unique(currentScan$AbsoluteScan)
    scanNumber <- unique(currentScan$Scan)
    scanSpeed <- unique(currentScan$Speed)
    scanInput1 <- paste("scan",scanDataInd,sep="")
    
    # Get value
    if(length(defaultScansToAnalyze)==0){
      checkValue <- 0
    } else{
      if(scanNumber%in%defaultScansToAnalyze){
        checkValue <- 1
      } else{
        checkValue <- 0
      } 
    }
    
    # Create observe event when someone selects a scan to include
    includeScan(scanDataInd,scanNumber)
    
    firstPlot <- div(
      class='column',style='width:fit-content;',
      div(class='row',style='margin:10px;',
          div(
            class='column',style='width:fit-content;',
            renderPlot(scanPlot,width=300,height=200)
          ),
          div(
            class='column',style='width:100px;align-items:flex-start;margin-left:10px;',
            p(HTML(paste("<b>Scan ",absoluteScanNumber,"</b>",sep="")),class='body'),
            p(HTML(paste(scanSpeed," rpm",sep="")),class='body'),
            checkboxInput(scanInput1,"Include",checkValue,width=300)
          )
      )
    )
    return(firstPlot)
    
  }
  
  # Render scan plot previews
  renderScanPlots <- function(defaultScansToAnalyze=NULL){
    
    output$dataPreview <- renderUI({
      
      dataLength <- length(dataList$scanData)
      lengthVec <- 1:dataLength
      
      # Get row indices
      rowInds <- lengthVec%%2
      
      # Calculate row as sum of previous rows
      rowVec <- cumsum(rowInds)
      
      # Get unique rows
      uniqueRows <- unique(rowVec)
      
      div(
        lapply(1:length(uniqueRows), function(i) { #generate row for each plot
          
          # Get current row plotting to
          currentRow <- uniqueRows[i]
          
          # Get scanData indices
          scanDataInds <- which(rowVec==currentRow)
          
          # Generate first scan plot
          currentScan <- dataList$scanData[[scanDataInds[1]]]
          firstPlot <- generateScanPlotPreview(currentScan,scanDataInd=scanDataInds[1],defaultScansToAnalyze)
          
          if(length(scanDataInds)>1){
            
            # Generate second scan plot
            nextScan <- dataList$scanData[[scanDataInds[2]]]
            secondPlot <- generateScanPlotPreview(nextScan,scanDataInd=scanDataInds[2],defaultScansToAnalyze)
            
            # Return ui
            div(style='display:flex;flex-direction:column',
                div(class='row',style='margin:10px;',
                    firstPlot,
                    secondPlot
                )
            ) # parent div 
            
          } else{
            
            # Return ui
            div(style='display:flex;flex-direction:column',
                div(class='row',style='margin:10px;',
                    firstPlot
                )
            ) # parent div
            
          }
          
        })
      )
    })
    
  }
  
  # When user uploads data
  observeEvent(input$uploadInput,{
    
    stop <- ''
    
    # Reset vars
    defineVars()
    updateLog('UltraAnalysis loaded')
    
    # Reset ui
    resetUI()
    
    # Get uploaded file
    uploadInput <- input$uploadInput
    file <- uploadInput$datapath
    filenames <- uploadInput$name
    
    # Update log
    updateLog(paste("Uploaded ",length(file)," files",sep=""))
    
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
          generateSectorHooks(x)
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
      renderScanPlots(defaultScansToAnalyze)
      
      # Render global inputs
      renderPsvInput()
      renderSdInput()
      
      # Close progress bar
      progress$close()
      
      #defaultScansToAnalyze <- NULL
    } else{
      updateLog(stop)
    }
    
  })
  
  # Render psv input
  renderPsvInput <- function(){
    output$psvInputArea <- renderUI({
      textInput('psvInput',NULL,dataList$psvValue,300,'Partial specific volume (mL/g)')
    })
  }
  renderSdInput <- function(){
    output$sdInputArea <- renderUI({
      textInput('sdInput',NULL,dataList$sdValue,300,'Solvent density (g/mL)')
    })
  }
  
  
  # ------------ PROCESS
  
  # Function to define a sector
  createSector <- function(scanNumber,data,receptor=0){
    
    randNum <- getRandomNumber()
    sectorID <- paste('r',randNum,sep="")
    rangeID <- paste('a',randNum,sep="")
    
    sector <- data.frame(
      Scan=scanNumber,
      Min=min(data$r),
      Max=max(data$r),
      n=nrow(data),
      Receptor=receptor,
      ID=sectorID,
      RangeID = rangeID,
      Color=getRandomColor()
    )
    
    out <- list()
    out$sector <- sector
    out$id <- sectorID
    return(out)
    
  }
  
  # Render hooks for saved sector
  generateSectorHooks <- function(sector){
    
    observeEvent(eval(parse(text=paste('input$',sector$ID,sep=""))),{
      editSelectedSector(sector$ID)
    })
    
    observeEvent(eval(parse(text=paste('input$',sector$RangeID,sep=""))),{
      zoomToSector(sector$ID)
    })
    
    
  }
  
  # Automatically find the meniscus
  findMenisci <- function(x,allScanData){
    
    tmp <- allScanData[[x]]
    
    sectorNum <- 3
    
    # Get scan number
    scanNumber <- unique(tmp$Scan)
    absoluteScanNumber <- unique(tmp$AbsoluteScan)
    
    updateMessage <- paste("Finding sectors for Scan ",absoluteScanNumber," (",x," of ",length(allScanData),")",sep="")
    tryCatch(updateLog(updateMessage),error=function(e)return())
    progress$set(message=updateMessage,
                 value=(x-1)/length(allScanData))
    
    # Get derivative for clustering
    tmp$Dr <- c(0,diff(tmp$Ar))
    
    # Get mean derivative
    meanDr <- mean(tmp$Dr)
    tmp$meanDr <- meanDr
    
    # Cluster into 3 groups by x,y position
    clustVec <- sqrt(tmp$Dr^2+tmp$r^2)
    kclust <- kmeans(clustVec,sectorNum)
    clusterVector <- kclust$cluster
    
    # ORder cluster vector from small to large to keep conssitent
    clusterVectorOrdered <- clusterVector[order(clusterVector)]
    
    tmp$Sector <- clusterVectorOrdered
    
    # Plot clusters
    clusterPlot <- ggplot(tmp,aes(x=r,y=Ar))+
      #geom_line()+
      geom_point(data=tmp,aes(y=Dr,col=as.character(Sector)))+
      geom_hline(yintercept=meanDr,col='red')+
      theme_prism()
    
    # Get sectors
    autoRanges <- lapply(1:sectorNum,function(sector){
      
      rawSectorData <- subset(tmp,Sector==sector)
      
      # Transform sector data to be linear
      sectorData <- data.frame(
        x = (rawSectorData$r^2),
        y = suppressWarnings(log(rawSectorData$Ar))
      )

      # Split data into groups of 10 datapoints
      n <- ceiling(nrow(sectorData)*0.25)
      errList <- lapply(1:(nrow(sectorData)-n),function(y){
     
        rawSub <- rawSectorData[y:(y+n),]
        sub <- sectorData[y:(y+n),]
        
        # Skip fits for data with NA or -Inf
        badInds <- which(is.na(sub$y)|sub$y==-Inf)
        if(length(badInds)>0){
          out <- list()
          out$metrics <- NA
          return(out)
        }
        
        # Format for lm.fit
        X <- cbind(1, as.matrix(sub$x))

        # Get linear fit for group 3
        subFit <- .lm.fit(X,sub$y)

        # # Get predictions
        # preds <- predict(tempFit,data.frame(r=sub$r,w=sub$w,R=sub$R,temp=sub$temp,r0=sub$r0))
        # 
        # Get proxy for normality of residuals as mean over median
        fitResiduals <- (subFit$residuals)
        
        # Get if normally distributed
        normalTest <- shapiro.test(fitResiduals)
        
        # Calcualte rsqr
        meanErr <- sum((sub$y-mean(sub$y,na.rm=TRUE))^2)
        pErr <- sum((fitResiduals^2))
        rsqr <- 1-(pErr/meanErr)
        
        # Calculate metric weighting distance to 1 by GoF
        outMetrics <-  (1-rsqr) * (1-normalTest$statistic[[1]])
        
        # Save output in list
        outlist <- list()
        outlist$data <- rawSub
        outlist$metrics <- outMetrics
        outlist$rsqr <- rsqr
        outlist$W <- normalTest$statistic[[1]]
        outlist$coef <- subFit$coefficients
        
        return(outlist)
        
      }) 
      errs <- unlist(lapply(errList,'[[','metrics'))
      
      # Get which minimizes err
      minInd <- which(errs==min(errs,na.rm=TRUE))[1]
      
      # Get best set of data
      bestData <- errList[[minInd]]
      
      # Plot
      bestFit <- function(x,slope,intercept){
        y <- x*slope+intercept
        return(r)
      }
      meniscusData <- bestData$data
      plotForBugfixing <- ggplot(meniscusData,
             aes(x=r^2,y=log(Ar)))+
        geom_point()+
        geom_abline(slope=bestData$coef[[2]],intercept=bestData$coef[[1]],col='red')+
        theme_prism()
      
      # Format as range
      rangeFormatList <- createSector(absoluteScanNumber,meniscusData)
      rangeFormat <- rangeFormatList$sector
      sectorID <- rangeFormatList$id
      
      # Create hook for editing 
      generateSectorHooks(rangeFormat)
      
      return(rangeFormat)
      
      
    }) %>% bind_rows()
    
    tryCatch(updateLog(paste("Sector finding complete for scan ",scanNumber,sep="")),error=function(e)return())
    
    # Return ranges
    return(autoRanges)
    
  }
  
  # Render selection of scan
  renderScanSelection <- function(selectedScan=NULL){
    
    if(length(selectedScan)==0){
      defaultSelection <- dataList$selectedScan
    } else{
      defaultSelection <- selectedScan
    }
    
    # Update saved sectors when changing which scans are included
    if(length(dataList$savedSectors)>0){
      sectorScans <- lapply(dataList$savedSectors,'[[',1)%>%unlist()
      savedSectorsToInclude <- which(sectorScans%in%dataList$absoluteScanNumbers)
      newSectors <- dataList$savedSectors[savedSectorsToInclude]
      updateDataList('savedSectors',newSectors)
    }
    
    output$scanSelectionInput <- renderUI({
      selectInput('selectScan',NULL,choices=dataList$absoluteScanNumbers,selected=defaultSelection,multiple=FALSE,width=300)
    })
    
  }
  
  # Function to render sector preview
  renderSectorPreview <- function(sectorLabel,bothPlots){

    # Render output
    output$sectorPlotFrame <- renderUI(
      
      div(
        class='column',
        p(sectorLabel),
        renderPlot(
          bothPlots,
          width=350,
          height=230
        ),
        div(style='margin:10px auto;',
            #textInput('receptorInput',NULL,NULL,width=300,'Receptor concentration (uM)')
        ),
        div(style='margin:-10px auto;',
            actionButton('saveSector',label="💾 Save sector")
        )
      )
      
    )
  }
  
  # Function to retrieve scan data given data list and a selected scan
  retrieveScanData <- function(){
    
    scanNames <- names(dataList$scanData)
    selectedScanInd <- which(dataList$selectedScan==scanNames)
    selectedScanData <- dataList$scanData[[selectedScanInd]]
    
    return(selectedScanData)
    
  }
  
  # Function to render scan plot
  renderScanPlot <- function(selectedSector=NULL,sectorLabel=NULL,bothPlots=NULL){

    if(length(dataList$selectedScan)==0){
      return()
    }
    
    if(length(dataList$scanData)==0){
      return()
    }
    
    if(length(dataList$scansToAnalyze)==0){
      return()
    }
    
    currentScanData <- retrieveScanData()
    
    if(length(currentScanData)==0){
      return()
    }
    
    if(!is.null(selectedSector)){
      
      scanDataToPlot <- subset(currentScanData,r>(as.numeric(selectedSector$xmin)*0.99)&r<(as.numeric(selectedSector$xmax)*1.01))
    
      renderSectorPreview(sectorLabel,bothPlots)
      
    } else{
      scanDataToPlot <- currentScanData
    }
 
    currentScanPlot <- ggplot(scanDataToPlot,aes(x=r,y=Ar))+
      geom_point()+
      geom_hline(yintercept=0,col='gray50',lty='dashed')+
      geom_hline(yintercept=1,col='gray50',lty='dashed')+
      theme_prism()
    
    # Check which selections are in the current scan
    savedSectorsDf <- bind_rows(dataList$savedSectors)
    if(nrow(savedSectorsDf)>0){
      
      inScan <- subset(savedSectorsDf,Scan==dataList$selectedScan)
      
      # Get data within ranges
      dataInRange <- lapply(1:nrow(inScan),function(x){
        
        row <- inScan[x,]
   
        sub <- subset(scanDataToPlot,r>=row$Min&r<=row$Max)
        
        if(nrow(sub)==0){
          return(NULL)
        }
        
        sub$Color <- row$Color
        return(sub)
      
      }) %>% bind_rows()
  
      if(nrow(dataInRange)>0){
        currentScanPlot <- currentScanPlot + 
          geom_point(data=dataInRange,aes(x=r,y=Ar,col=Color))+
          scale_color_identity()
      }
      
    }
    
    scaleValue <- ifelse(length(input$scaleSelectionPlot)>0,input$scaleSelectionPlot,FALSE)
    if(scaleValue){
      currentScanPlot <- currentScanPlot +
        scale_y_continuous(limits=c(0,1))
    }
    
    #Render output
    output$currentScanPlotOutput <- renderPlot(
      currentScanPlot,width=600,height=300
    )
    
  }
  
  # Observe scaling y axis
  observeEvent(input$scaleSelectionPlot,{
    renderScanPlot()
  })
  
  # When user selects a scan to see up close
  observeEvent(input$selectScan,{
    
    selectedScan <- as.numeric(input$selectScan)
    renderPlot <- ifelse(length(dataList$skipPlotRerender)>0,dataList$skipPlotRerender,FALSE)
    
    if(!is.na(selectedScan)&&!renderPlot){
      # Update data list
      updateDataList('selectedScan',selectedScan)
      # Plot
      resetBrush()
      renderScanPlot()
    } else{
      updateDataList('skipPlotRerender',FALSE)
    }
    
  })
  
  # Render plots from selected sector
  renderSelectedSectorPlots <- function(currentScanData,selectedSector){
    
    # Get sector data from scan data
    sectorData <- subset(currentScanData,r>selectedSector$xmin&r<selectedSector$xmax)
    
    if(nrow(sectorData)==0){
      return()
    }
    
    # Transform
    transy <- tryCatch(
      log(sectorData$Ar),
      warning=function(w){
        warningString <- paste("<b>Warning</b> in ln(Ar) of selected sector: ",w[[1]],sep="")
        updateLog(warningString)
        return(log(sectorData$Ar))
      }
    )
    transformed <- data.frame(
      x = sectorData$r^2,
      y = transy
    )
    
    # Get linear fit
    linearFit <- tryCatch(lm(y~x,data=transformed),error=function(e)return(NA))
    linearFitSummary <- tryCatch(summary(linearFit),error=function(e)return(NA))
    rsqr <- tryCatch(linearFitSummary$r.squared,error=function(e)return(NA))
    coefs <- tryCatch(coef(linearFitSummary),error=function(e)return(c(0,0)))
    
    # Predict from fit
    preds <- tryCatch(predict(linearFit,transformed),error=function(e)return(NA))
    transformed$pAr <- preds
    
    # Generate df to match with currently selected area
    tmpData <- data.frame(
      Scan = dataList$selectedScan,
      Min = min(sectorData$r),
      Max = max(sectorData$r),
      n = nrow(sectorData)
    )
    if(nrow(tmpData)==0){
      return()
    }
    
    # Prevents rendering of plots and ui if selected range is the same as what is already displayed
    checkForSector <- tryCatch(print(currentSectorDf),error=function(e)return(NA))
    if(length(checkForSector)>1){
      if(identical(tmpData,currentSectorDf)){
        return()
      }
    }
    
    # Format current sector dataframe
    currentSectorDf <- data.frame(
      Scan = dataList$selectedScan,
      Min = min(sectorData$r),
      Max = max(sectorData$r),
      n = nrow(sectorData)
    )
    if(nrow(currentSectorDf)==0){
      return()
    }
    
    # Format sector label
    sectorLabel <- paste(currentSectorDf$Min," - ",currentSectorDf$Max,
                         " | n = ",currentSectorDf$n,
                         " | R^2 = ",round(rsqr,3),sep="")
    
    
    # Plot
    currentSectorPlot <- ggplot(transformed,aes(x=x,y=y))+
      geom_point(size=1)+
      geom_abline(slope=coefs[2],intercept=coefs[1],col='red')+
      scale_y_continuous(n.breaks=3)+
      ylab("ln(Ar)")+
      xlab("r^2")+
      theme_prism()+
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=11)
      )
    
    # Plot residuals
    residualPlot <- ggplot(transformed,aes(x=x,y=y-pAr))+
      geom_point(size=1)+
      geom_hline(yintercept=0,col='gray20')+
      scale_y_continuous(n.breaks=3)+
      ylab("y - ŷ")+
      xlab('r^2')+
      theme_prism()+
      theme(axis.text = element_text(size=10),
            axis.title = element_text(size=11)
      )
    
    # combine plots
    bothPlots <- plot_grid(currentSectorPlot,residualPlot,
                           rel_heights=c(1,1),ncol=1,align='v')
    
    # Update scan plot to zoom
    renderScanPlot(selectedSector,sectorLabel,bothPlots)
    
  }
  
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
    renderSelectedSectorPlots(currentScanData,selectedSector)
    
    # Stop skipping rerender
    updateDataList('skipPlotRerender',FALSE)
    
  })
  
  # Function to reset brush
  resetBrush <- function(){
    session$resetBrush('selectingSector')
    output$sectorPlotFrame <- NULL
  }
  
  # When a user clicks the plot, reset to normal view
  observeEvent(input$resetPlotView,{
    resetBrush()
    renderScanPlot()
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
      tryCatch(updateLog(newSectorLog),error=function(e)return())
      
      # Save as list
      savedSectors[[length(savedSectors)+1]] <- currentSectorDf
      
      # Update dataList
      updateDataList('savedSectors',savedSectors)
      
      # Create hook for editing 
      generateSectorHooks(currentSectorDf)
      
      # Update UI
      renderSavedSelections()
      
      # Update plot
      resetBrush()
      renderScanPlot()
    }
    
    
  })

  # Function to render saved selections
  renderSavedSelections <- function(){

    if(length(dataList$savedSectors)==0){
      output$savedSectorUI <- NULL
      return()
    }
    
    # Update UI
    output$savedSectorUI <- renderUI({
      
      div(class='column',style='align-items:center;width;100%;',
        lapply(1:length(dataList$savedSectors),function(x){
          
          tmpSector <- dataList$savedSectors[[x]]
          sectorID <- tmpSector$ID
          rangeID <- tmpSector$RangeID
          colorString <- paste(as.numeric(col2rgb(tmpSector$Color)[,1]),collapse=',')
          
          div(
            class='row',style=paste('margin:5px;width:100%;border-radius:10px;background-color:rgba(',colorString,',0.2)',sep=""),
            div(
              class='column',style='align-items:center;width:10%;margin-left:5%',
              p(HTML("<b>Scan</b>"),style='font-size:12px;'),
              p(tmpSector$Scan,style='font-size:10px;')
            ),
            div(
              class='column',style='align-items:center;width:45%',
              p(HTML("<b>Range</b>"),style='font-size:12px;'),
              actionButton(rangeID,
                           paste(tmpSector$Min,' - ',tmpSector$Max,sep=""),
                           style='font-size:10px;padding:0px;background-color:transparent;border:0px;color:black;font-weight:normal;text-decoration:underline')
            ),
            div(
              class='column',style='align-items:center;width:15%',
              p(HTML("<b>uM</b>"),style='font-size:12px;'),
              p(tmpSector$Receptor,style='font-size:10px;')
            ),
            div(class='column',style='align-items:center;width:15%;justify-content:center;margin-left:10px',
              actionButton(sectorID,label="🗑️",style='font-size:12px;background-color:transparent;border: 0px solid black')
            )
          ) # End row
        })
      )
      
    })
  }
  
  # Define function to remove a sector
  editSelectedSector <- function(sectorID){

    savedSectors <- dataList$savedSectors
    
    ids <- lapply(savedSectors,'[[',6) %>% unlist()
    idToRemove <- which(ids==sectorID)[1]
    
    sectorToRemove <- savedSectors[[idToRemove]]
    # Update log
    removeSectorLog <- paste("Removed sector: Scan ",sectorToRemove$Scan," Min ",sectorToRemove$Min,' Max ',sectorToRemove$Max, " n = ",sectorToRemove$n,sep="")
    updateLog(removeSectorLog)
    
    savedSectors <- savedSectors[-idToRemove]
    updateDataList('savedSectors',savedSectors)
    
    #Update UI
    renderSavedSelections()
    renderScanPlot()
  }
  
  # Define function to zoom to a sector
  zoomToSector <- function(sectorID){
    
    savedSectors <- dataList$savedSectors
    
    ids <- lapply(savedSectors,'[[',6) %>% unlist()
    idToZoom <- which(ids==sectorID)[1]
    
    sectorToZoom <- savedSectors[[idToZoom]]
    
    sectorList <-list()
    sectorList$xmin <- sectorToZoom$Min
    sectorList$xmax <- sectorToZoom$Max
    
    # Switch scan data
    if(sectorToZoom$Scan!=dataList$selectedScan){
      updateDataList('skipPlotRerender',TRUE)
      renderScanSelection(sectorToZoom$Scan)
    }
    
    updateDataList('selectedSector',sectorList)
    updateDataList('selectedScan',sectorToZoom$Scan)
    
    # Render plot
    resetBrush()

    # Get current scan data
    currentScanData <- retrieveScanData()
    
    # Render sector plots
    renderSelectedSectorPlots(currentScanData,sectorList)
    
  }
  
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
    identifiedSectors <- lapply(1:length(allScanData),findMenisci,allScanData) %>% bind_rows()
    
    # Covnert each row to a list
    identifiedSectorsList <- lapply(1:nrow(identifiedSectors),function(x){return(identifiedSectors[x,])})
    
    # Update saved ranges
    # Save as list
    savedSectors <- dataList$savedSectors
    newSavedSectors <- c(savedSectors,identifiedSectorsList)
    
    # Update dataList
    updateDataList('savedSectors',newSavedSectors)
    
    # Update UI
    renderSavedSelections()
    
    # Update plot
    resetBrush()
    renderScanPlot()
    
    progress$close()
    
  })
  
  # ------------ FIT
  
  # Render the model type to display
  renderModelTypeSelection<-function(model){
    output$modelTypeSelection <- renderUI({
      selectInput('modelType',NULL,choices=modelTypes,selected=model,multiple=FALSE,width=300)
    })
  }
  
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
  
  # Create function to check for proper inputs before fitting
  checkForInputs <- function(){
    
    if(length(dataList$selectedModel)==0){
      selectedModel <- input$modelType
      updateDataList('selectedModel',selectedModel)
    }
    
    # Get psv and sd inputs
    
    # Format check vector
    checkVec <- c(input$psvInput,input$sdInput)
    if(length(checkVec)==0){
      updateLog("Invalid partial specific volume (PSV) input. Defaulting to 0.71 mL/g")
      psv <- 0.71
      updateLog("Invalid solvent density (SD). Defaulting to 1.003 g/mL")
      sd <- 1.003
    } else{
      
      if(input$psvInput==""){
        updateLog("Invalid partial specific volume (PSV). Defaulting to 0.71 mL/g")
        psv <- 0.71
      } else{
        psv <- as.numeric(input$psvInput)
      }
      
      if(input$sdInput==""){
        updateLog("Invalid solvent density (SD) input. Defaulting to 1.003 g/mL")
        sd <- 1.003
      } else{
        sd <- as.numeric(input$sdInput)
      }
    }
    
    # Update datalist
    updateDataList('psvValue',psv)
    updateDataList('sdValue',sd)
    
    # Get saved sectors
    savedSectors <- dataList$savedSectors 
    
    # Loop through saved sectors
    if(length(savedSectors)==0){
      updateLog("No sectors to fit. Save sectors to fit a model.")
      return(FALSE)
    }
    
    return(TRUE)
    
  }
  
    # Define function parameters
    .GlobalEnv$defineFitParms <- function(){
      
      if(dataList$selectedModel=='MW / Single ideal species'){
        
        dataCols <- c("r","r0","w","temp","Ar","R","ID")
        localParms <- c("A0","offset")
        globalParms <- c("Mb")
        addedParms <- NULL
        
        fitFunction <- "
          function(DATA_COLUMNS,GLOBAL_PARAMETERS,LOCAL_PARAMETERS){
            LOCAL_FIT_BOILERPLATE
            
            Ar <- (offset) + ( (A0) * exp(1)^( ( (Mb*w^2)/(R*temp*2) )*(r^2 - r0^2) ) )
            
            return(Ar)
          }
        "
        
      } else if(dataList$selectedModel=='Kd / Monomer-Nmer'){
        
        # Get fit-specific inputs
        mwValue <- tryCatch(as.numeric(input$mwInput),error=function(e)return(80000))
        mwValue <- ifelse(is.na(mwValue),80000,mwValue)
        eCoefValue <- tryCatch(as.numeric(input$eCoefInput),error=function(e)return(76000))
        eCoefValue <- ifelse(is.na(eCoefValue),76000,eCoefValue)
        nValue <- tryCatch(as.numeric(input$nInput),error=function(e)return(2))
        nValue <- ifelse(is.na(nValue),2,nValue)

        # Update datalist
        updateDataList('mwValue',mwValue)
        updateDataList('eCoefValue',eCoefValue)
        updateDataList('nValue',nValue)
        
        dataCols <- c("r","r0","w","temp","Ar","R","Mb","N","ID")
        localParms <- c("A0","offset")
        globalParms <- c("K")
        addedParms <- c("Mb" = mwValue*(1-(dataList$psvValue*dataList$sdValue)), # covnert MW to Mb
                        "N" = nValue)
        
        fitFunction <- "
          function(DATA_COLUMNS,GLOBAL_PARAMETERS,LOCAL_PARAMETERS){
            LOCAL_FIT_BOILERPLATE
          
            theta <- w^2 / (2*R*temp)
            monomerTerm <- A0 * exp(1)^( Mb * theta * (r^2 - r0^2) )
            dimerTerm <-  N * K * A0^N * exp(1)^( N * Mb * theta * (r^2 - r0^2) )
            Ar <- offset + monomerTerm + dimerTerm
            
            return(Ar)
          }
        "
        
      }
      
      
      out <- list()
      out$dataCols <- dataCols
      out$localParms <- localParms
      out$globalParms <- globalParms
      out$addedParms <- addedParms
      out$fitFunction <- fitFunction
      
      return(out)
    }
    
    # Calculate reference radius for each sector
    .GlobalEnv$calculateReferenceRadius <- function(x,savedSectors,allScans,fitParms){
      
      dataCols <- fitParms$dataCols
      addedParms <- fitParms$addedParms
      
      # Get sector
      sector <- savedSectors[[x]]
      
      # Get data from scan
      scanInd <- which(allScans==as.numeric(sector$Scan))
      sectorScan <- dataList$scanData[[scanInd]]
      
      # Subset out sector
      data <- subset(sectorScan,r>=sector$Min&r<=sector$Max)
      
      # Check that scan and sector match
      checkScan <- data$AbsoluteScan[1]==sector$Scan
      checkMinr <- min(data$r)==sector$Min
      checkMaxr <- max(data$r)==sector$Max
      
      if(!checkScan||!checkMinr||!checkMaxr){
        print(paste("In fit for baseline: data not subset for x = ",x,sep=""))
        tryCatch(updateLog("<b>Fatal error:</b> data mismatch. Sector not properly subset."),error=function(e)return())
        return(NA)
      }
      
      # Add parms in added parms
      if(length(addedParms)>0){
        for(i in 1:length(addedParms)){
          parm <- addedParms[i]
          data[,names(parm)] <- as.numeric(parm)
        } 
      }
      
      # Define r0 reference radius as 2/3 down sector per hetero documentation
      r0 <- (diff(range(data$r))*(2/3))+min(data$r)
      data$r0 <- r0
      
      # Assign ID
      data$ID <- x
      
      # Select only data that will be used in fits
      fitData <- data %>% select(all_of(dataCols))
      
      # add absolute scan column
      fitData$AbsoluteScan <- data$AbsoluteScan
      
      # Return data
      return(fitData)
      
  }
    
    # Function to dynamically generate function string from parms
    .GlobalEnv$processFitFunction <- function(fitParms,globalData){
      
      dataCols <- fitParms$dataCols
      globalParms <- fitParms$globalParms
      localParms <- fitParms$localParms
      
      # Define local fit boilerplate
      localFitBoilerplate <- "
        A0 <- lapply(paste('A0_',ID,sep=''),function(x){get(x)}) %>% unlist()
        offset <- lapply(paste('offset_',ID,sep=''),function(x){get(x)}) %>% unlist()
      "
      
      # Process inputs
      dataColumnList <- paste(dataCols,collapse=',')
      
      globalParmListStart <- paste(globalParms,'=1',sep='',collapse=',')
      
      parmArrayLength <- length(dataList$savedSectors)
      localParmList <- lapply(localParms,function(x){
        vec<-rep(x,parmArrayLength)
        new <- paste(vec,'_',1:parmArrayLength,sep="")
        return(new)
      }) %>% unlist()
      localParmListCollapsed <- paste(localParmList,collapse=",")
      localParmListStart <- paste(localParmList,'=1',sep='',collapse=',')
      
      startString <- paste('c(',globalParmListStart,',',localParmListStart,")",sep='')
      
      # Substitute parameters into fit function
      dynamicFitFunction <- gsub("LOCAL_FIT_BOILERPLATE",localFitBoilerplate,fitParms$fitFunction)
      dynamicFitFunction <- gsub('DATA_COLUMNS',dataColumnList,dynamicFitFunction)
      dynamicFitFunction <- gsub('LOCAL_PARAMETERS',localParmListCollapsed,dynamicFitFunction)
      dynamicFitFunction <- gsub('GLOBAL_PARAMETERS',globalParms,dynamicFitFunction)
      
      # Load into real function
      loadedFunction <- eval(parse(text=dynamicFitFunction))
      .GlobalEnv$loadedFunction <- loadedFunction
      
      # Format the fit string
      fitString <- "
        gsl_nls(
          Ar~loadedFunction(DATA_COLUMNS,GLOBAL_PARAMETERS,LOCAL_PARAMETERS),
          data=globalData,
          algorithm='lm',
          start=START_STRING
        )
      "
      
      # Configure parameters
      dynamicFitString <- gsub('DATA_COLUMNS',dataColumnList,fitString)
      dynamicFitString <- gsub('GLOBAL_PARAMETERS',globalParms,dynamicFitString)
      dynamicFitString <- gsub('LOCAL_PARAMETERS',localParmListCollapsed,dynamicFitString)
      dynamicFitString <- gsub('START_STRING',startString,dynamicFitString)
      
      # Format output list
      out <- list()
      out$fitFunction <- loadedFunction
      out$fitString <- dynamicFitString
      out$localParmList <- localParmList
      return(out)
    }
    
    # Create function to plot nonlinear fit
    .GlobalEnv$plotNonlinearFit <- function(globalFitDup){
      
      nonlinear <- ggplot(globalFitDup,aes(x=r))+
        geom_point(aes(y=Ar-offset))+
        geom_line(aes(y=pAr-offset,group=ID),col='red',lwd=1)+
        theme_prism()
      
      tryCatch(print(nonlinear),
               warning=function(w){
                 warningString <- paste("<b>Warning</b> in nonlinear plotting: ",w[[1]],sep="")
                 updateLog(warningString)
                 return()
               })
      
      return(nonlinear)
      
    }
    
    # Create function to plot residuals
    .GlobalEnv$plotNonlinearResiduals <- function(globalFitDup){
      
      
      residualPlot <- ggplot(globalFitDup,aes(x=r-r0,y=Ar-pAr,col=as.character(ID)))+
        geom_point()+
        geom_hline(yintercept=0,col='red')+
        theme_prism()+
        theme(legend.position = 'none')
      
      return(residualPlot)
      
    }
    
    # Create function to plot linear fit
    .GlobalEnv$plotLinearFit <- function(globalFitDup,showSimulation=FALSE){
      
      linear <- ggplot(globalFitDup,aes(x=r2_minus_r0_over_2))+
        geom_point(aes(y=ln_A_minus_offset_over_A0))+
        geom_line(aes(x=r2_minus_r0_over_2,y=linearpAr),col='red',lwd=1)+
        xlab("(r^2 - r0^2) / 2")+
        ylab("ln( (Ar - offset) / A0 )")+
        theme_prism()
      
      # Predict if monomer
      if(showSimulation==TRUE){
        
        monomerData <- singleIdealSpecies(r=globalFitDup$r,
                                          w=globalFitDup$w,
                                          R=globalFitDup$R,
                                          temp=globalFitDup$temp,
                                          r0=globalFitDup$r0,
                                          A0=globalFitDup$A0,
                                          Mb=globalFitDup$Mb,
                                          offset=globalFitDup$offset)
        dimerData <- singleIdealSpecies(r=globalFitDup$r,
                                        w=globalFitDup$w,
                                        R=globalFitDup$R,
                                        temp=globalFitDup$temp,
                                        r0=globalFitDup$r0,
                                        A0=globalFitDup$A0,
                                        Mb=(globalFitDup$Mb*2),
                                        offset=globalFitDup$offset)
        
        globalFitDup$monomerSimulation <- monomerData
        globalFitDup$dimerSimulation <- dimerData
        
        linear <- linear +
          geom_line(data=globalFitDup,aes(x=r2_minus_r0_over_2,y=log((monomerSimulation-offset)/A0),color='Simulated Monomer',group=1),lwd=1)+
          geom_line(data=globalFitDup,aes(x=r2_minus_r0_over_2,y=log((dimerSimulation-offset)/A0),color='Simulated Dimer',group=1),lwd=1)+
          scale_color_manual(
            values=c(
              "Simulated Dimer"='goldenrod',
              "Simulated Monomer"='firebrick'
            )
          )+
          theme(legend.position='top')
        
      }
      
      tryCatch(print(linear),
               warning=function(w){
                 warningString <- paste("<b>Warning</b> in linear plotting: ",w[[1]],sep="")
                 updateLog(warningString)
                 return()
               })
      
      return(linear)
      
    }
    
    # Process fit results
    .GlobalEnv$processFitResults <- function(globalFit,globalData,processedFitFunction,fitParms){
      
      localParms <- fitParms$localParms
      globalParms <- fitParms$globalParms
      
      # Extract fitted Mb
      fitSummary <- summarizeNLS(globalFit,dataList$selectedModel)
      parmRownames <- rownames(fitSummary$coefs)
      
      # Create copy of df to avoid modification
      globalFitDup <- globalData
      
      # Update log
      tryCatch(updateLog(fitSummary$summary),error=function(e){return()})
      
      # Add local and global parm columns to data
      parmVector <- c(localParms,globalParms)
      for(i in seq_along(parmVector)){
        
        parm <- parmVector[i]
        globalFitDup[,parm] <- NA
      }
      
      # Fill columns with data
      allParms <- c(processedFitFunction$localParmList,globalParms)
      for(i in seq_along(allParms)){
        
        # Get parm
        parm <- allParms[i]
        
        # Extract coefficient
        estimateInd <- which(parmRownames==parm)
        
        # Get relevant subset of data
        parmSplit <- strsplit(parm,'_')[[1]]
        
        if(length(parmSplit)>1){
          dataInd <- which(globalFitDup$ID==as.numeric(parmSplit[2])) 
        } else{
          dataInd <- 1:nrow(globalFitDup)
        }
        parmType <- parmSplit[1]
        
        # Add parm to df
        globalFitDup[dataInd,parmType] <- fitSummary$coef$Estimate[estimateInd]
        
      }
      
      # predict
      preds <- predict(globalFit,globalFitDup)
      globalFitDup$pAr <- preds
      
      # Calculate rsqr
      meanErr <- sum((globalFitDup$Ar-mean(globalFitDup$Ar))^2)
      pErr <- sum((globalFitDup$Ar-globalFitDup$pAr)^2)
      rsqr <- 1-(pErr/meanErr)
      
      # Plot non linear fit
      nonlinear <- plotNonlinearFit(globalFitDup)
      residuals <- plotNonlinearResiduals(globalFitDup)
      
      # Format output based on model
      if(dataList$selectedModel=='MW / Single ideal species'){
        coefRow <- fitSummary$coefs[which(rownames(fitSummary$coefs)=='Mb'),]
        displayedMW <- (coefRow$Estimate / (1-(dataList$psvValue*dataList$sdValue)))/1000
        displayedMWErr <- (coefRow$`Std. Error` / (1-(dataList$psvValue*dataList$sdValue)))/1000
        displayedText <- paste(
          "MW = ",round(displayedMW,1)," +/- ",round(displayedMWErr,1)," kDa",sep=""
        )
      } else if(dataList$selectedModel=='Kd / Monomer-Nmer'){
        coefRow <- fitSummary$coefs[which(rownames(fitSummary$coefs)=='K'),]
        displayedKd <- ((1/coefRow$Estimate)*(dataList$nValue/dataList$eCoefValue))*1000000
        displayedKdErr <- (coefRow$`Std. Error`/coefRow$Estimate)*displayedKd
        displayedText <- paste(
          "Kd = ",round(displayedKd,2)," +/- ",round(displayedKdErr,2)," uM",sep=""
        )
      }
      
      parmRow <- div(
        p(HTML(displayedText),style='padding-left:10px;'),
        renderTable(fitSummary$coefs,rownames = TRUE)
      )
      
      # Return out list
      outlist <- list()
      outlist$nonlinear <- nonlinear
      outlist$data <- globalFitDup
      outlist$description <- parmRow
      outlist$residuals <- residuals
      
      return(outlist)
      
    }

    .GlobalEnv$fitPlotWidth <- 290
    .GlobalEnv$fitPlotHeight <- 300
  
  # When user tries to fit the data
  observeEvent(input$fitData,{
    
    checkInputs <- tryCatch(checkForInputs(),error=function(e)return(FALSE))
    
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
    .GlobalEnv$fitParms <- defineFitParms()
    
    # Processes sector data per fit function, adding parms accordingly
    processedSectors <- lapply(1:length(savedSectors),
                               calculateReferenceRadius,
                               savedSectors=savedSectors,
                               allScans=allScans,
                               fitParms=fitParms)
    
    if(length(which(is.na(processedSectors)))>0){
      print("Error in fitting baseline")
      tryCatch(updateLog("Error fitting baselines"),error=function(e)return())
      return()
    }
    
    # Combine into 1 global dataset
    globalData <- processedSectors %>% bind_rows()

    # Dynamically render fit function from parms
    processedFitFunction <- processFitFunction(fitParms,globalData)
   
    # Do fit
    globalFit <- eval(parse(text=processedFitFunction$fitString))
    
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
    renderFitOutput()
    progress$close()
  })
  
  # Function to render the fit output
  renderFitOutput <- function(){
    output$fitResult <- renderUI({
      
      #Define download handler
      output$downloadFit <- downloadHandler(
        filename=function(){
          paste("UA_FIT_",paste(unlist(str_extract_all(Sys.time(),"[:digit:]")),collapse=""),".csv",sep="")
        },
        content=function(file){
          write.csv(fit$data,file=file,row.names=FALSE)
        }
      )#download handler
      
      output$downloadExperiment <- downloadHandler(
        filename=function(){
          paste("UA_EXP_",paste(unlist(str_extract_all(Sys.time(),"[:digit:]")),collapse=""),".Rdata",sep="")
        },
        content=function(file){
          save(dataList,file=file)
        }
      )
      
      # UI
      div(
        class='row',
        div(
          class='column',style='width:fit-contents;margin-left:15px;',
          plotOutput('nonlinearPlot',click='clickFitResultNonlinear',dblclick='resetFitPlots',height = paste(fitPlotHeight,'px',sep="")),
        ),
        div(
          class='column',style='width:fit-contents;margin-left:5px;',
          plotOutput('residualPlot',click='clickFitResultNonlinearResiduals',dblclick='resetFitPlots',height = paste(fitPlotHeight,'px',sep=""))
        ),
        div(
          class='column',style='margin-left:15px',
          div(class='row',style=paste('overflow-y: auto;height:',fitPlotHeight/2,'px;',sep=""),
              fit$description
          ),
          div(
            class='row',style='margin-top:15px;align-self:flex-end',
            p('Download',class='header')
          ),
          div(
            class='row',style='align-self:flex-end;',
            downloadButton("downloadFit","Fit data",icon=icon('database'))
          ),
          div(
            class='row',style='align-self:flex-end;margin-top:10px',
            downloadButton("downloadExperiment","Experiment file",icon=icon('flask'))
          ),
        )
      )
    })
  }
  
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

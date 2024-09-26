
# Function to retrieve scan data given data list and a selected scan
retrieveScanData <- function(){
  
  scanNames <- names(dataList$scanData)
  selectedScanInd <- which(dataList$selectedScan==scanNames)
  selectedScanData <- dataList$scanData[[selectedScanInd]]
  
  return(selectedScanData)
  
}

# Automatically find the meniscus
findMenisci <- function(x,input,output,session,allScanData){
  
  tmp <- allScanData[[x]]
  
  sectorNum <- 3
  
  # Get scan number
  scanNumber <- unique(tmp$Scan)
  absoluteScanNumber <- unique(tmp$AbsoluteScan)
  
  updateMessage <- paste("Finding sectors for Scan ",absoluteScanNumber," (",x," of ",length(allScanData),")",sep="")
  tryCatch(updateLog(output,updateMessage),error=function(e)return())
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
    generateSectorHooks(input,output,session,rangeFormat)
    
    return(rangeFormat)
    
    
  }) %>% bind_rows()
  
  tryCatch(updateLog(output,paste("Sector finding complete for scan ",scanNumber,sep="")),error=function(e)return())
  
  # Return ranges
  return(autoRanges)
  
}

# Render selection of scan
renderScanSelection <- function(output,selectedScan=NULL){
  
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
renderSectorPreview <- function(output,sectorLabel,bothPlots){
  
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

# Function to render scan plot
renderScanPlot <- function(input,output,selectedSector=NULL,sectorLabel=NULL,bothPlots=NULL){
  
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
    
    renderSectorPreview(output,sectorLabel,bothPlots)
    
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
    
    if(nrow(inScan)>0){
      
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

# Render plots from selected sector
renderSelectedSectorPlots <- function(input,output,currentScanData,selectedSector){
  
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
      updateLog(output,warningString)
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
  renderScanPlot(input,output,selectedSector,sectorLabel,bothPlots)
  
}

# Function to define a sector
createSector <- function(scanNumber,data,receptor=0){
  
  randNum <- getRandomNumber()
  sectorID <- paste('r',randNum,sep="")
  rangeID <- paste('a',randNum,sep="")
  
  sector <- tibble(
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

# Function to render saved selections
renderSavedSelections <- function(output){
  
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
editSelectedSector <- function(input,output,sectorID){
  
  savedSectors <- dataList$savedSectors
  
  ids <- lapply(savedSectors,'[[',6) %>% unlist()
  idToRemove <- which(ids==sectorID)[1]
  
  sectorToRemove <- savedSectors[[idToRemove]]
  # Update log
  removeSectorLog <- paste("Removed sector: Scan ",sectorToRemove$Scan," Min ",sectorToRemove$Min,' Max ',sectorToRemove$Max, " n = ",sectorToRemove$n,sep="")
  updateLog(output,removeSectorLog)
  
  savedSectors <- savedSectors[-idToRemove]
  updateDataList('savedSectors',savedSectors)
 
  #Update UI
  renderSavedSelections(output)
  renderScanPlot(input,output)
}

# Define function to zoom to a sector
zoomToSector <- function(input,output,session,sectorID){
  
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
  resetBrush(output,session)
  
  # Get current scan data
  currentScanData <- retrieveScanData()
  
  # Render sector plots
  renderSelectedSectorPlots(input,output,currentScanData,sectorList)
  
}

# Render hooks for saved sector
generateSectorHooks <- function(input,output,session,sector){
  
  observeEvent(eval(parse(text=paste('input$',sector$ID,sep=""))),{
    editSelectedSector(input,output,sector$ID)
  })
  
  observeEvent(eval(parse(text=paste('input$',sector$RangeID,sep=""))),{
    zoomToSector(input,output,session,sector$ID)
  })
  
  
}

# Function to reset brush
resetBrush <- function(output,session){
  session$resetBrush('selectingSector')
  output$sectorPlotFrame <- NULL
}

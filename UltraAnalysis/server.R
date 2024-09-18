
# Define fitting function
.GlobalEnv$singleIdealSpecies <- function(r, w, R, temp, r0, A0, Mb, offset){
  Ar <- offset + A0 * exp(1)^( ( (Mb*w^2)/(R*temp*2) )*(r^2 - r0^2) )
  return(Ar)
}

# Read each scan
.GlobalEnv$readScans <- function(file){
  
  # Read each scan
  readAScan <- function(x,file){
    
    row <- file[x]
    
    # Read data
    .GlobalEnv$scanLines <- read.table(row,header=TRUE,sep="\t")
    
    # Get scan name
    scanName <- names(scanLines)
    
    # Extract metadata
    rawMeta <- scanLines[1,]
    metaVec <- unlist(strsplit(rawMeta," "))
    
    cellNumber <- as.numeric(metaVec[2])
    temp <- as.numeric(metaVec[3])+273.15
    rpm <- as.numeric(metaVec[4])
    w <- rpm * 0.1047 
    
    # Get only values
    valueLines <- scanLines[-c(1:2),]
    
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
    
    R <- 8.314 * 1000 * 100 * 100 # (g * cm^2) / s^2 * mol * K
    
    # Create df
    dfValues <- data.frame(
      Name = NA,
      Cell = cellNumber,
      Speed = rpm,
      Scan = x,
      r = radialPosition,
      r0=NA,
      w = w,
      R = R,
      temp = temp,
      A0 = NA,
      offset = NA,
      Ar = absorbance,
      A0Err = NA,
      OffsetErr = NA
    ) %>% mutate_all(as.numeric)
    dfValues$Name <- scanName
    dfValues$File <- file[x]
    
    
    naInds <- which(is.na(dfValues$Ar))
    if(length(naInds)>0){
      curatedValues <- dfValues[-naInds,]
    } else{
      curatedValues <- dfValues
    }
    
    return(curatedValues)
    
  }
  out <- lapply(1:length(file),readAScan,file=file)
  
  # Name scan list with condition
  cell <- unique(out[[1]]$Cell)
  scan <- unique(out[[1]]$Scan)
  speed <- unique(out[[1]]$Speed)
  cellScanSpeed <- paste("Cell ",cell," Scan ",scan," ",speed," rpm",sep="")
  names(out) <- cellScanSpeed

  return(out)
  
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
  convegeTolerance <- fitSummary$convInfo$finTol
  
  # Format summary string
  summaryString <- paste(
    "Fitting for: ",type," | ",itToConverge," iterations to convergence | Convergence tolerance: ",convegeTolerance,
    " | ",parmsFit," parameters fit: ",coefInfo,sep=""
    
  )
  
  # Save coefficients and summary string
  out <- list()
  out$summary <- summaryString
  out$coefs <- fitCoefs
  return(out)
  
}

# Read test data 
testFiles <- list.files('www/testData/',full.names=TRUE)
.GlobalEnv$testData <- lapply(testFiles,readScans)

# Plot
bothPlots <- ggplot(testData[[1]][[1]],aes(x=r,y=Ar))+
  geom_point()+
  theme_prism()+
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )
sectorLabel <- "testing"

# Define colros to use
colorVec <- rainbow(10)

# Define model types
.GlobalEnv$modelTypes <- c("MW / Single ideal species")
psv <- 0.71
sd <- 1.003

# Define starting vars
defineVars <- function(){
  .GlobalEnv$logText <- NULL
  .GlobalEnv$psvValue <- NULL
  .GlobalEnv$sdValue <- NULL
  .GlobalEnv$scansToAnalyze <- 1
  .GlobalEnv$scanData <- NULL #testData[[1]]
  .GlobalEnv$savedSectors <- NULL #list(data.frame(Scan=1,Min=5.887,Max=6.046,n=68,Receptor=0,ID='r1726670506',Color='#00FFFF'))
  .GlobalEnv$selectedModel <- modelTypes[1]
}
defineVars()

testplot <- ggplot(data.frame(x=1:10,y=1:10),aes(x=x,y=y))+
  geom_point()


# Define server logic required to draw a histogram
function(input, output, session) {

  # Update the log text
  updateLog <- function(update){
    
    logText <- c(update,logText)
    .GlobalEnv$logText <- logText
    
    output$log <- renderUI({
      source('www/log.R')[[1]]
    })
    
  }
  
  # Reset variables on refresh
  observe({
    defineVars()
    updateLog('UltraAnalysis loaded')
  })
  
  output$upload <- renderUI({
    source('www/upload.R')[[1]]
  })  
  
  output$fit <- renderUI({
    source('www/fit.R')[[1]]
  })  
  
  #Render output
  output$currentScanPlotOutput <- NULL
  
  renderProcess <- function(){
    output$process <- renderUI({
      source('www/process.R')[[1]]
    })
  }
  renderProcess()

  
  # ------------ UPLOAD
  
  
  # When user uploads data
  observeEvent(input$uploadInput,{
    
    # Reset vars
    defineVars()
    updateLog('UltraAnalysis loaded')
    
    # Get uploaded file
    .GlobalEnv$uploadInput <- input$uploadInput
    
    file <- uploadInput$datapath
    
    # Update log
    updateLog(paste("Uploaded ",length(file)," files",sep=""))
    
    # Create progress bar
    progress <- shiny::Progress$new()
    progress$set(message='Loading uploaded file...',value=0.1)

    # Read each scan
    progress$set(message='Reading scan data...',value=0.5)
    .GlobalEnv$scanData <- readScans(file)
    
    # Format scan plot output
    output$dataPreview <- renderUI({
      
      dataLength <- length(scanData)
      lengthVec <- 1:dataLength
      
      # Get row indices
      rowInds <- lengthVec%%2
      
      # Calculate row as sum of previous rows
      rowVec <- cumsum(rowInds)
      
      # Get unique rows
      uniqueRows <- unique(rowVec)
      
      div(
        lapply(1:length(uniqueRows), function(i) { #generate row for each plot
          
          currentRow <- uniqueRows[i]
          
          # Get scanData indices
          scanDataInds <- which(rowVec==currentRow)
          
          currentScan <- scanData[[scanDataInds[1]]]
          scanPlot <- ggplot(currentScan,aes(x=r,y=Ar))+
            geom_line()+
            xlab("r (cm)")+
            theme_prism()
          
          scanNumber <- unique(currentScan$Scan)
          scanSpeed <- unique(currentScan$Speed)
          scanInput1 <- paste("scan",scanDataInds[1],sep="")
          
          # Create observe event when someone selects a scan to include
          observeEvent(eval(parse(text=paste('input$scan',scanDataInds[1],sep=""))),{
            # Retrieve value and scan
            value <- eval(parse(text=paste('input$scan',scanDataInds[1],sep="")))
            
            if(value){
              scansToAnalyze <- c(scanNumber,scansToAnalyze)
              scansToAnalyze <- scansToAnalyze[order(scansToAnalyze)]
            } else{
              scansToAnalyze <- scansToAnalyze[-which(scansToAnalyze==scanNumber)]
            }

            .GlobalEnv$scansToAnalyze <- scansToAnalyze
            
            if(length(scansToAnalyze)>0){
              updateLog(paste("Scans included: ",paste(scansToAnalyze,collapse=", "),sep=""))
            } else{
              if(!is.null(scansToAnalyze)){
                updateLog("No scans included")
              }
            }
            
            # Update process tab
            renderProcess()
            
          })
          
          firstPlot <- div(
            class='column',style='width:fit-content;',
            div(class='row',style='margin:10px;',
                div(
                  class='column',style='width:fit-content;',
                  renderPlot(scanPlot,width=300,height=200)
                ),
                div(
                  class='column',style='width:100px;align-items:flex-start;margin-left:10px;',
                  p(HTML(paste("<b>Scan ",scanNumber,"</b>",sep="")),class='body'),
                  p(HTML(paste(scanSpeed," rpm",sep="")),class='body'),
                  checkboxInput(scanInput1,"Include",0,width=300)
                )
            )
          )
          
          
          
          if(length(scanDataInds)>1){
            
            nextScan <- scanData[[scanDataInds[2]]]
            scanPlot2 <- ggplot(nextScan,aes(x=r,y=Ar))+
              geom_line()+
              xlab("r (cm)")+
              theme_prism()
            
            nextScanNumber <- unique(nextScan$Scan)
            nextScanSpeed <- unique(nextScan$Speed)
            scanInput2 <- paste("scan",scanDataInds[2],sep="")
            
            # Create second observe event
            observeEvent(eval(parse(text=paste('input$scan',scanDataInds[2],sep=""))),{
              # Retrieve value and scan
              value <- eval(parse(text=paste('input$scan',scanDataInds[2],sep="")))
        
              if(value){
                scansToAnalyze <- c(nextScanNumber,scansToAnalyze)
                scansToAnalyze <- scansToAnalyze[order(scansToAnalyze)]
              } else{
                scansToAnalyze <- scansToAnalyze[-which(scansToAnalyze==nextScanNumber)]
              }
              .GlobalEnv$scansToAnalyze <- scansToAnalyze
              
              if(length(scansToAnalyze)>0){
                updateLog(paste("Scans included: ",paste(scansToAnalyze,collapse=", "),sep="")) 
              } else{
                if(!is.null(scansToAnalyze)){
                  updateLog("No scans included")
                }
              }
              
              # Update process tab
              renderProcess()
              
            })
            
            secondPlot <- div(
              class='column',style='width:fit-content;',
              div(class='row',style='margin:10px;',
                  div(
                    class='column',style='width:fit-content;',
                    renderPlot(scanPlot2,300,height=200)
                  ),
                  div(
                    class='column',style='width:100px;align-items:flex-start;margin-left:10px;',
                    p(HTML(paste("<b>Scan ",nextScanNumber,"</b>",sep="")),class='body'),
                    p(HTML(paste(nextScanSpeed," rpm",sep="")),class='body'),
                    checkboxInput(scanInput2,"Include",0,width=300)
                  )
              )
            )
            
            div(style='display:flex;flex-direction:column',
                div(class='row',style='margin:10px;',
                    firstPlot,
                    secondPlot
                )
            ) # parent div 
          } else{
            div(style='display:flex;flex-direction:column',
                div(class='row',style='margin:10px;',
                    firstPlot
                )
            ) # parent div
          }
          
        })
      )
    })

    progress$close()
    
  })
  
  # Function to render scan plot
  output$sectorPlotFrame <- NULL
  renderScanPlot <- function(selectedSector=NULL,sectorLabel=NULL,bothPlots=NULL){
    
    if(length(currentScanData)==0){
      return()
    }
    
    if(!is.null(selectedSector)){
      
      scanDataToPlot <- subset(currentScanData,r>(as.numeric(selectedSector$xmin)*0.99)&r<(as.numeric(selectedSector$xmax)*1.01))
    
      # Render output
      output$sectorPlotFrame <- renderUI(
        
        div(
          class='column',
          p(sectorLabel),
          renderPlot(
            bothPlots,
            width=350,
            height=180
          ),
          div(style='margin:10px auto;',
              textInput('receptorInput',NULL,NULL,width=300,'Receptor concentration (uM)')
          ),
          div(style='margin:-10px auto;',
              actionButton('saveSector',label="💾 Save sector")
          )
        )
        
      )
      
    } else{
      scanDataToPlot <- currentScanData
    }
 
    currentScanPlot <- ggplot(scanDataToPlot,aes(x=r,y=Ar))+
      geom_point()+
      theme_prism()
    
    # Check which selections are in the current scan
    savedSectors <- bind_rows(savedSectors)
    if(nrow(savedSectors)>0){
      
      inScan <- subset(savedSectors,Scan==selectedScan)
      
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
    
    #Render output
    output$currentScanPlotOutput <- renderPlot(
      currentScanPlot,width=600,height=300
    )
    
  }
  
  
  # ------------ PROCESS
  
  
  # When user selects a scan to see up close
  observeEvent(input$selectScan,{
    .GlobalEnv$selectedScan <- input$selectScan
    
    if(selectedScan!=""){
      # Get scan data
      .GlobalEnv$currentScanData <- scanData[[as.numeric(selectedScan)]]
      
      # Plot
      renderScanPlot()

    }
    
  })
  
  # When user selects a brush
  observeEvent(input$selectingSector,{
    
    selectedSector <- input$selectingSector

    # Get sector data from scan data
    sectorData <- subset(currentScanData,r>selectedSector$xmin&r<selectedSector$xmax)
    
    # Get linear fit
    linearFit <- tryCatch(lm(Ar~r,data=sectorData),error=function(e)return(NA))
    linearFitSummary <- tryCatch(summary(linearFit),error=function(e)return(NA))
    rsqr <- tryCatch(linearFitSummary$r.squared,error=function(e)return(NA))
    coefs <- tryCatch(coef(linearFitSummary),error=function(e)return(c(0,0)))
    
    # Predict from fit
    preds <- tryCatch(predict(linearFit,sectorData),error=function(e)return(NA))
    sectorData$pAr <- preds
    
    tmpData <- data.frame(
      Scan = selectedScan,
      Min = min(sectorData$r),
      Max = max(sectorData$r),
      n = nrow(sectorData)
    )
    
    # Prevents rendering of plots and ui if selected range is the same as what is already displayed
    checkForSector <- tryCatch(print(currentSectorDf),error=function(e)return(NA))
    if(length(checkForSector)>1){
      if(identical(tmpData,currentSectorDf)){
        return()
      }
    }
    
    # Format current sector dataframe
    .GlobalEnv$currentSectorDf <- data.frame(
      Scan = selectedScan,
      Min = min(sectorData$r),
      Max = max(sectorData$r),
      n = nrow(sectorData)
    )
    
    # Format sector label
    sectorLabel <- paste(currentSectorDf$Min," - ",currentSectorDf$Max,
                         " | n = ",currentSectorDf$n,
                         " | R^2 = ",round(rsqr,3),sep="")
    
    
    # Plot
    currentSectorPlot <- ggplot(sectorData,aes(x=r,y=Ar))+
      geom_point(size=1)+
      geom_abline(slope=coefs[2],intercept=coefs[1],col='red')+
      scale_y_continuous(n.breaks=3)+
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
    residualPlot <- ggplot(sectorData,aes(x=r,y=Ar-pAr))+
      geom_point(size=1)+
      geom_hline(yintercept=0,col='gray20')+
      scale_y_continuous(n.breaks=3)+
      ylab("Ar - pAr")+
      theme_prism()+
      theme(axis.text = element_text(size=10),
            axis.title = element_text(size=11)
            )
    
    # combine plots
    bothPlots <- plot_grid(currentSectorPlot,residualPlot,
                           rel_heights=c(1,1),ncol=1,align='v')
    
    # Update scan plot to zoom
    renderScanPlot(selectedSector,sectorLabel,bothPlots)
    
  })
  
  # When a user clicks the plot, reset to normal view
  observeEvent(input$resetPlotView,{
    renderScanPlot()
  })

  # Function to render saved selections
  renderSavedSelections <- function(){

    if(length(savedSectors)==0){
      output$savedSectorUI <- NULL
      return()
    }
    
    # Update UI
    output$savedSectorUI <- renderUI({
      
      div(class='column',style='align-items:center;width;100%;',
        lapply(1:length(savedSectors),function(x){
          
          tmpSector <- savedSectors[[x]]
          sectorID <- tmpSector$ID
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
              p(paste(tmpSector$Min,' - ',tmpSector$Max,sep=""),style='font-size:10px;')
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
  renderSavedSelections()
  
  # Define function to edit a secotr
  editSelectedSector <- function(receptorID){
    ids <- lapply(savedSectors,'[[',6) %>% unlist()
    idToRemove <- which(ids==receptorID)
    
    sectorToRemove <- savedSectors[[idToRemove]]
    # Update log
    removeSectorLog <- paste("Removed sector: Scan ",sectorToRemove$Scan," Min ",sectorToRemove$Min,' Max ',sectorToRemove$Max, " n = ",sectorToRemove$n,sep="")
    updateLog(removeSectorLog)
    
    savedSectors <- savedSectors[-idToRemove]
    .GlobalEnv$savedSectors <- savedSectors
    
    #Update UI
    renderSavedSelections()
    renderScanPlot()
  }
  
  # When user saves a selection
  observeEvent(input$saveSector,{
    
    
    # Get receptor concentration
    receptorConcentration <- ifelse(input$receptorInput=="",0,input$receptorInput)
    currentSectorDf$Receptor <- receptorConcentration
    
    # Define sector ID
    sectorID <- paste('r',as.character(round(as.numeric(Sys.time()))),sep="")
    
    # Check if sector ID already exists
    otherIDs <- unlist(lapply(savedSectors,'[[',6))
    
    if(!(sectorID%in%otherIDs)){
      currentSectorDf$ID <- sectorID
      
      # Define sector color
      currentSectorDf$Color <- sample(colorVec,1)
      
      newSectorLog <- paste("Added sector: Scan ",currentSectorDf$Scan," Min ",currentSectorDf$Min,' Max ',currentSectorDf$Max, " n = ",currentSectorDf$n,sep="")
      updateLog(newSectorLog)
      
      # Save as list
      savedSectors[[length(savedSectors)+1]] <- currentSectorDf
      .GlobalEnv$savedSectors <- savedSectors
    
      # Create hook for editing 
      observeEvent(eval(parse(text=paste('input$',sectorID,sep=""))),{
        editSelectedSector(sectorID)
      })
  
      # Update UI
      renderSavedSelections()
      
      # Update plot
      session$resetBrush('selectingSector')
      renderScanPlot()
    }
    
    
  })
  
  
  # ------------ FIT
  
  # When user selects a model type
  observeEvent(input$modelType,{
  
    # Render UI for model parms
    .GlobalEnv$selectedModel <- input$modelType
    
  
    output$modelParms <- renderUI({
      
      if(selectedModel=='MW / Single ideal species'){
        
        div(class='column',
          actionButton('fitData','🏃 Run fit',style='align-self:center')
        )
      }

    })
    
  })
  
  
  # When user tries to fit the data
  observeEvent(input$fitData,{
    
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
        psv <- input$psvInput
      }
      
      if(input$sdInput==""){
        updateLog("Invalid solvent density (SD) input. Defaulting to 1.003 g/mL")
        sd <- 1.003
      } else{
        sd <- input$sdInput
      }
    }
    
   
    
    # Loop through saved sectors
    if(length(savedSectors)==0){
      updateLog("No sectors to fit. Save sectors to fit a model.")
      return()
    }
    
    # Fit baselines
    baselineFits <- lapply(1:length(savedSectors),function(x){
      
      # Get sector
      sector <- savedSectors[[x]]
      
      # Get data from scan
      sectorScan <- scanData[[as.numeric(sector$Scan)]]
      
      # Subset out sector
      data <- subset(sectorScan,r>=sector$Min&r<=sector$Max)
      
      # Define r0 reference radius
      r0 <- median(data$r)
      data$r0 <- r0
      
      # Fit baseline and reference absorbance
      baselineFit <- gsl_nls(
        Ar~singleIdealSpecies(r,w,R,temp,r0,A0, Mb, offset),
        data=data,
        algorithm='lm',
        start=c(A0=1,offset=1,Mb=1),
        lower=c(A0=0,offset=NULL,Mb=0)
      )
      
      # Get fit info
      fitSummary <- summarizeNLS(baselineFit,"baseline")
      
      # Update log
      tryCatch(updateLog(fitSummary$summary),error=function(e)return())
      
      # Get coefficients
      coefs <- fitSummary$coefs
      
      # Extract A0 and offset
      A0 <- coefs$Estimate[which(rownames(coefs)=="A0")]
      A0Err <- coefs$`Std. Error`[which(rownames(coefs)=="A0")]
      offset <- coefs$Estimate[which(rownames(coefs)=="offset")]
      offsetErr <- coefs$`Std. Error`[which(rownames(coefs)=="offset")]
      
      # Add to data
      data$A0 <- A0
      data$A0Err <- A0Err
      data$offset <- offset
      data$OffsetErr <- offsetErr
      data$ID <- paste("GROUP",x,sep="")
      
      # Return data
      return(data)
      
      
    }) %>% bind_rows()
    
    if(selectedModel=='MW / Single ideal species'){
      
      # Fit for MW
      mwFit <-  gsl_nls(
        Ar~singleIdealSpecies(r,w,R,temp,r0,A0,Mb,offset),
        data=baselineFits,
        algorithm='lm',
        start=c(Mb=1)
      )
      
      # Extract fitted Mb
      mwFitSummary <- summarizeNLS(mwFit,selectedModel)
      
      # Update log
      tryCatch(updateLog(mwFitSummary$summary),error=function(e){return()})
      
      # Get coef
      Mb <- mwFitSummary$coefs[which(rownames(mwFitSummary$coefs)=="Mb"),]
      
      # Convert Mb to MW
      mw <- (Mb$Estimate/(1-(psv*sd)))/1000
      mwErr <- (Mb$`Std. Error`/(1-(psv*sd)))/1000
      
      # predict
      preds <- predict(mwFit,baselineFits)
      baselineFits$pAr <- preds
      
      # Calculate rsqr
      meanErr <- sum((baselineFits$Ar-mean(baselineFits$Ar))^2)
      pErr <- sum((baselineFits$Ar-baselineFits$pAr)^2)
      rsqr <- 1-(pErr/meanErr)
      
      # Plot non linear fit
      nonlinear <- ggplot(baselineFits,aes(x=r))+
        geom_point(aes(y=Ar-offset))+
        geom_line(aes(y=pAr-offset,group=ID),col='red')+
        theme_prism()
      
      # Plot linear
      linear <- ggplot(baselineFits,aes(x=((r^2-r0^2)/2)))+
        geom_point(aes(y=log((Ar-offset)/A0)))+
        geom_line(aes(y=log((pAr-offset)/A0)),col='red')+
        xlab("(r^2 - r0^2) / 2")+
        ylab("ln( (Ar - offset) / A0 )")+
        theme_prism()
      
      # Format plot title
      plotTitle <- paste("MW = ",round(mw,1)," +/- ",round(mwErr,1)," kDa",
                         "\n Buoyant MW = ",round(Mb$Estimate/1000,1)," +/- ",round(Mb$`Std. Error`/1000,1)," kDa, R^2 = ",round(rsqr,3),sep="")
      
      # Format grid
      plotgrid <- plot_grid(
        nonlinear,
        linear,
        ncol=2,
        nrow=1,
        align='h'
      )
   
      
    }
    
    # Render ui
    output$fitResult <- renderUI({
      
      #Define download handler
      output$downloadFit<-downloadHandler(
        
        filename=function(){
          paste("ULTRAANALYSIS_",paste(unlist(str_extract_all(Sys.time(),"[:digit:]")),collapse=""),".csv",sep="")
        },
        content=function(file){
          write.csv(baselineFits,file=file,row.names=FALSE)
        }
        
      )#download handler
      
    
      div(
        class='row',
        div(
          class='column',style='width:fit-contents;margin-left:15px;',
          renderPlot(plotgrid,width=750,height=300)
        ),
        div(
          class='column',style='width:fit-contents;margin-left:30px',
          div(class='row',
            div(
              class='column',
              p(HTML(paste("<b>MW</b> = ",round(mw,1)," +/- ",round(mwErr,1)," kDa",sep="")),class='body'),
              p(HTML(paste("<b>Buoyant MW</b> = ",round(Mb$Estimate/1000,1)," +/- ",round(Mb$`Std. Error`/1000,1)," kDa",sep="")),class='body'),
              p(HTML(paste("<b>R^2</b> = ",round(rsqr,3),sep="")),class='body') 
            )
          ),
          div(
            class='row',style='align-self:center;margin-top;10px',
            downloadButton("downloadFit","Download")
          )
        )
      )
      
      
    })
    
    
    
    
  })
  
  

}

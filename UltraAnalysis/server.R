

# Read each scan
readScans <- function(file){
  
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
    
    # Create df
    dfValues <- data.frame(
      Name = NA,
      Cell = cellNumber,
      Speed = rpm,
      Scan = x,
      r = radialPosition,
      r0=NA,
      w = w,
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


# Define starting vars
defineVars <- function(){
  .GlobalEnv$logText <- NULL
  .GlobalEnv$psvValue <- NULL
  .GlobalEnv$sdValue <- NULL
  .GlobalEnv$scansToAnalyze <- NULL
  .GlobalEnv$scanData <- NULL
  .GlobalEnv$savedSectors <- list()
}
defineVars()

testplot <- ggplot(data.frame(x=1:10,y=1:10),aes(x=x,y=y))+
  geom_point()


# Define server logic required to draw a histogram
function(input, output, session) {

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
  
  # Update the log text
  updateLog <- function(update){
    
    logText <- c(update,logText)
    .GlobalEnv$logText <- logText
    
    output$log <- renderUI({
      source('www/log.R')[[1]]
    })
    
  }
  updateLog('UltraAnalysis loaded')

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
  
  # When user selects a scan to see up close
  observeEvent(input$selectScan,{
    .GlobalEnv$selectedScan <- input$selectScan
    
    if(selectedScan!=""){
      # Get scan data
      .GlobalEnv$currentScanData <- scanData[[as.numeric(selectedScan)]]
      
      # Plot
      currentScanPlot <- ggplot(currentScanData,aes(x=r,y=Ar))+
        geom_point()+
        theme_prism()
      
      #Render output
      output$currentScanPlotOutput <- renderPlot(
        currentScanPlot,width=600,height=300
      )
    }
    
  })
  
  # When user selects a brush
  observeEvent(input$selectingSector,{
    
    .GlobalEnv$selectedSector <- input$selectingSector
    
    # Get sector data from scan data
    sectorData <- subset(currentScanData,r>selectedSector$xmin&r<selectedSector$xmax)
    
    # Get linear fit
    linearFit <- lm(Ar~r,data=sectorData)
    linearFitSummary <- summary(linearFit)
    rsqr <- linearFitSummary$r.squared
    coefs <- coef(linearFitSummary)
    
    # Predict from fit
    preds <- predict(linearFit,sectorData)
    sectorData$pAr <- preds
    
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
      geom_point()+
      geom_abline(slope=coefs[2],intercept=coefs[1],col='red')+
      theme_prism()+
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
    # Plot residuals
    residualPlot <- ggplot(sectorData,aes(x=r,y=Ar-pAr))+
      geom_point()+
      geom_hline(yintercept=0,col='gray20')+
      scale_y_continuous(n.breaks=3)+
      ylab("Ar - pAr")+
      theme_prism()
    
    # combine plots
    bothPlots <- plot_grid(currentSectorPlot,residualPlot,
                           rel_heights=c(1,1),ncol=1,align='v')
    
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
          removeSectorID <- paste('removeSector',x,sep="")
          
          div(
            class='row',style='margin:5px;width:100%;',
            div(
              class='column',style='align-items:center;width:25%',
              p(HTML("<b>Scan</b>"),style='font-size:12px;'),
              p(tmpSector$Scan,style='font-size:10px;')
            ),
            div(
              class='column',style='align-items:center;width:50%',
              p(HTML("<b>Range</b>"),style='font-size:12px;'),
              p(paste(tmpSector$Min,' - ',tmpSector$Max,sep=""),style='font-size:10px;')
            ),
            div(
              class='column',style='align-items:center;width:25%',
              p(HTML("<b>uM</b>"),style='font-size:12px;'),
              p(ifelse(tmpSector$Receptor=="",0,tmpSector$Receptor),style='font-size:10px;')
            )
          ) # End row
        })
      )
      
    })
  }
  
  # When a user undos a selection
  observeEvent(input$undoSector,{
    
    
    removeSectorLog <- paste("Removed sector: Scan ",savedSectors[[length(savedSectors)]]$Scan," Min ",savedSectors[[length(savedSectors)]]$Min,' Max ',savedSectors[[length(savedSectors)]]$Max, " n = ",savedSectors[[length(savedSectors)]]$n,sep="")
    updateLog(removeSectorLog)
    
    
    savedSectors <- savedSectors[-length(savedSectors)]
    .GlobalEnv$savedSectors <- savedSectors
    renderSavedSelections()
    
  })
  
  # When user saves a selection
  observeEvent(input$saveSector,{
    
    
    # Get receptor concentration
    receptorConcentration <- input$receptorInput
    currentSectorDf$Receptor <- receptorConcentration
    
    newSectorLog <- paste("Added sector: Scan ",currentSectorDf$Scan," Min ",currentSectorDf$Min,' Max ',currentSectorDf$Max, " n = ",currentSectorDf$n,sep="")
    updateLog(newSectorLog)
    
    # Save as list
    savedSectors[[length(savedSectors)+1]] <- currentSectorDf
    .GlobalEnv$savedSectors <- savedSectors
    
    print(savedSectors)
    
    # Update UI
    renderSavedSelections()
    
    
  })
  


}


# Read uploaded scans
readScans <- function(file,filenames){
  
  # Read each scan
  readAScan <- function(x,file,filenames){
    
    tmpFile <- file[x]
    filename <- filenames[x]
    
    # Extract metadata
    metaData <- read_lines(tmpFile,n_max=2)
    scanName <- metaData[1]
    metaDf <- metaData[2] %>%
      strsplit(' ') %>%
      unlist() %>%
      as.matrix() %>% 
      t() %>%
      `colnames<-`(c("uk1","Cell","Celcius","Rpm","Seconds","w_squared_times_t","uk2",'uk3')) %>%
      as.data.frame() %>%
      mutate(across(-1,as.numeric)) %>%
      as_tibble()
    cellNumber <- metaDf$Cell
    temp <- metaDf$Celcius+273.15
    rpm <- metaDf$Rpm
    w <- (rpm/60)*2*pi
    
    # Read data as tibble
    scanData <- read_lines(tmpFile,skip=2)
    scanDf <- scanData %>%
      str_squish() %>%
      read_table(col_names=c('r',"Ar",'noise'))
    
    # Define gas constant with proper units
    R <- 8.3144 * 1000 * 100 * 100 # (g * cm^2) / s^2 * mol * K
    
    # Get scan number
    scanNum <- filename %>%
      str_sub(end=str_locate(filename,'\\.')[[1]]-1) %>%
      str_remove_all('[:alpha:]|[:punct:]') %>%
      as.numeric()
    
    # Create df
    dfValues <- tibble(
      Name = scanName,
      Cell = cellNumber,
      Speed = rpm,
      Scan = x,
      AbsoluteScan = scanNum,
      CellScanSpeed = paste("Cell ",cellNumber," Scan ",scanNum," ",rpm," rpm",sep=""),
      r = scanDf$r,
      Ar = scanDf$Ar,
      r0=NA,
      w = w,
      temp = temp,
      R = R,
      theta = w^2 / (2*R*temp),
      Noise=scanDf$noise,
      File = file[x]
    ) %>% mutate(across(-c("Name","CellScanSpeed","File"),as.numeric))
    
    return(dfValues)
  }
  out <- lapply(1:length(file),readAScan,file=file,filenames=filenames)
  
  # Get scan names
  scanNames <- lapply(out,'[[','AbsoluteScan')%>%
    unlist()%>%
    unique()
  names(out) <- scanNames
  
  return(out)
  
}

# Render psv input
renderSdInput <- function(output){
  output$sdInputArea <- renderUI({
    textInput('sdInput',NULL,dataList$sdValue,300,'Solvent density (g/mL)')
  })
}
renderPsvInput <- function(output){
  output$psvInputArea <- renderUI({
    textInput('psvInput',NULL,dataList$psvValue,300,'Partial specific volume (mL/g)')
  })
}

# Function when someone checks a scan to include
includeScan <- function(input,scanDataInd,scanNumber){
  
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
generateScanPlotPreview <- function(input,currentScan,scanDataInd,defaultScansToAnalyze){
  
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
  includeScan(input,scanDataInd,scanNumber)
  
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
renderScanPlots <- function(input,output,defaultScansToAnalyze=NULL){
  
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
        firstPlot <- generateScanPlotPreview(input,currentScan,scanDataInd=scanDataInds[1],defaultScansToAnalyze)
        
        if(length(scanDataInds)>1){
          
          # Generate second scan plot
          nextScan <- dataList$scanData[[scanDataInds[2]]]
          secondPlot <- generateScanPlotPreview(input,nextScan,scanDataInd=scanDataInds[2],defaultScansToAnalyze)
          
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




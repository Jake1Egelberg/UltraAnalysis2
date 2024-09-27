# Read default fit inputs
.GlobalEnv$userInputs <- read_csv('www/xInputs.csv',show_col_types =FALSE)

# Define fit function boilerplate
fitFunction <- "
  function(DATA_COLUMNS,GLOBAL_PARAMETERS,LOCAL_PARAMETERS){
    LOCAL_FIT_BOILERPLATE
    
    MODEL
    
    return(Ar)
  }
"

# Define local fit boilerplate
localFitBoilerplate <- "
  A0 <- lapply(paste('A0_',ID,sep=''),function(x){get(x)}) %>% unlist()
  offset <- lapply(paste('offset_',ID,sep=''),function(x){get(x)}) %>% unlist()
"

# Render the model type to display
renderModelTypeSelection<-function(output,model){
  output$modelTypeSelection <- renderUI({
    selectInput('modelType',NULL,choices=modelTypes,selected=model,multiple=FALSE,width=300)
  })
}

# Define function to get model file from selected model type
getModelInfo <- function(modelFile=NULL){
  
  if(length(modelFile)==0){
    modelFile <- modelFiles[str_detect(modelFiles,paste(dataList$selectedModel,'.txt',sep=""))]
  } 
  
  loadedModel <- read_lines(modelFile)
  modelTibble <- data.frame(loadedModel) %>% 
    t() %>%
    `colnames<-`(c("Model","Local","Global","Inputs","Mutations","Displayed","DisplayedUnits")) %>%
    as_tibble()
  return(modelTibble)
}

# Render model parms
renderModelParms <- function(input,output){
  
  # Render UI for model parms
  selectedModel <- input$modelType
  updateDataList('selectedModel',selectedModel)
  
  # Get model tibble 
  modelTibble <- getModelInfo()
  # Get inputs
  inputs <- strsplit(modelTibble$Inputs,',') %>% unlist()
  
  output$modelParms <- renderUI({
    

    div(class='column',style='height:250px;overflow-y:auto;',
        actionButton('fitData','🏃 Run fit',style='align-self:center;margin-bottom:10px;margin-top:10px;'),
        div(
          lapply(seq_along(inputs),function(x){
            tmpInput <- inputs[x]
            tmpCorrespondingRow <- subset(userInputs,Key==tmpInput)
            tmpCorrespondingValue <- as.numeric(dataList[str_detect(tmpInput,names(dataList))])
            
            textInput(tmpCorrespondingRow$Input,NULL,tmpCorrespondingValue,300,tmpCorrespondingRow$Placeholder)
          })
        )
    )
    
  })
  
  output$fitDescription <- renderText({
    paste("Global fit for: ",dataList$selectedModel,sep="")
  })
  
}

# Create generic check for input function
  # If input not formatted properly, then updates data list with default and returns FALSE
  # If input formatted properly, updates data list with input and returns true
checkNumericInput <- function(key,input,output){
 
  # Update console
  updateString <- paste("Validating input for: ",key,sep="")
  #updateLog(output,updateString)
  
  # Get input for key
  correspondingInput <- subset(userInputs,Key==key)
  
  if(nrow(correspondingInput)==0){
    errorString <- paste("Error checking numeric input: key/input mismatch for: ",key,sep="")
    tryCatch(updateLog(output,errorString),error=function(e)return())
    return()
  }
  
  # Get user input
  newInput <- input[[correspondingInput$Input]]
  
  # Check for null or empty string
  failedCheck <- length(newInput)==0 || newInput=="" || is.na(as.numeric(newInput))
  
  
  if(failedCheck){
    warningString <- paste("<b>Warning</b>: Invalid ",key," input. Defaulting to ",correspondingInput$Default,sep="")
    updateLog(output,warningString)
    updateDataList(key,correspondingInput$Default)
    return(FALSE)
  } else{
    validatedString <- paste(key," validated. Applying value of ",as.numeric(newInput),sep="")
    updateLog(output,validatedString)
    updateDataList(key,as.numeric(newInput))
    return(TRUE)
  }
  
}

# Create function to check for proper inputs before fitting
checkForInputs <- function(input,output){
  
  if(length(dataList$selectedModel)==0){
    selectedModel <- input$modelType
    updateDataList('selectedModel',selectedModel)
  }
  
  # Get saved sectors
  savedSectors <- dataList$savedSectors 
  
  # Loop through saved sectors
  if(length(savedSectors)==0){
    updateLog(output,"No sectors to fit. Save sectors to fit a model.")
    return(FALSE)
  }
  
  return(TRUE)
  
}

# Check fit-specific inputs
checkMonomerNmer <- function(input,output){
  
  # Check fit-specific inputs
  checkNumericInput(key='mwValue',input$mwInput,output,'mwValue',80000)
  checkNumericInput(key='nValue',input$nInput,output,'nValue',2)
  checkNumericInput(key='eCoefValue',input$eCoefInput,output,'eCoefValue',76000)
  
}

# Define function to calculate columns based on user inputs as additional fit parms
calculateAdditionColumns <- function(x,coefs=NULL){
  
  # Get index of equals
  eqtStartInd <- str_locate(x,'=')[[1]]
  
  # Get equation after  equals
  equation <- str_sub(x,start=eqtStartInd+1)
  
  # Get dataList values
  dataListValues <- dataList[str_detect(x,names(dataList))]
  
  # Get all identified variables
  identifiedVariables <- equation %>%
    str_extract_all("\\b[a-zA-Z]+") %>% unlist()
  
  # Check if all identified variables are accounted for, if they are not, then pull from coefs
  matchInds <- which(identifiedVariables%in%names(dataListValues))
  lengthOfMatch <- length(matchInds)
  lengthOfAllVars <- length(identifiedVariables)
  if(lengthOfMatch!=lengthOfAllVars){
    
    # Check that coefs is not null
    if(length(coefs)==0){
      tryCatch(updateLog(output,'<b>Fatal error</b>: attempted to pull coef from empty fit'),error=function(e)return())
      return()
    }
    
    # Get unidentified variables
    unidentifiedVars <- identifiedVariables[-matchInds]
    for(i in 1:length(unidentifiedVars)){
      
      # Get match from coefs
      coefMatch <- coefs[which(rownames(coefs)==unidentifiedVars[i]),]
      
      if(nrow(coefMatch)==0){
        tryCatch(updateLog(output,'<b>Fatal error</b>: attempted to pull coef that doesnt exist in fit'),error=function(e)return())
        return()
      }
      
      # Add to data list values
      dataListValues[unidentifiedVars[i]] <- coefMatch$Estimate
    }
    
  }
  
  # Substitute into equation
  for(y in 1:length(dataListValues)){equation <- gsub(names(dataListValues[y]),as.numeric(dataListValues[y]),equation)}
  
  # calculate value
  calculatedValue <- equation %>%
    parse(text=.) %>%
    eval()
  
  # Get value name
  valueName <- str_sub(x,end=eqtStartInd-1)
  
  # Format vector output
  vectorOutput <- c(calculatedValue)
  names(vectorOutput) <- valueName
  
  return(vectorOutput)
}

# Define function parameters
defineFitParms <- function(input,output){
  
  # Identify model file
  modelFile <- modelFiles[str_detect(modelFiles,paste(dataList$selectedModel,'.txt',sep=""))]

  # Load model
  #modelFile <- "C:/Users/Jake/Documents/Code/UltraAnalysis2/UltraAnalysis/www/Models/Single ideal species.txt"
  #modelFile <- "C:/Users/Jake/Documents/Code/UltraAnalysis2/UltraAnalysis/www/Models/Monomer-Nmer.txt"
  #modelFile <- "C:/1_Documents/Ferguson Lab/UltraAnalysis2/UltraAnalysis/www/Models/Monomer-Nmer.txt"
  #modelFile <- "C:/1_Documents/Ferguson Lab/UltraAnalysis2/UltraAnalysis/www/Models/Single ideal species.txt"
  #modelFile <- "C:/1_Documents/Ferguson Lab/UltraAnalysis2/UltraAnalysis/www/Models/Monomer-Nmer-Mmer.txt"
  modelTibble <- getModelInfo(modelFile)
  .GlobalEnv$modelTibble <- modelTibble
  

  # Extract model variables
  modelVariables <- str_replace_all(modelTibble$Model,'[^[:alpha:]|0]|exp'," ") %>%
    strsplit(" ") %>%
    unlist() %>%
    unique() %>%
    .[.!=""]

  # Define as data, local, or global per file
  localParms <- strsplit(modelTibble$Local,',') %>% unlist()
  globalParms <- strsplit(modelTibble$Global,',') %>% unlist()
  dataCols <- c(modelVariables[!modelVariables%in%c(localParms,globalParms)],"ID")

  # Check necessary inputs
  inputstoCheck <- strsplit(modelTibble$Inputs,",") %>% 
    unlist() %>%
    lapply(checkNumericInput,input=input,output=output)
  
  # Re-render input text boxes with new values
  renderModelParms(input,output)
  
  # Calculate any necessary additional columns
  if(modelTibble$Mutations!='NA'){
    # Format mutations to perform
    addedParms <- strsplit(modelTibble$Mutations,',')[[1]] %>%
      lapply(calculateAdditionColumns) %>%
      unlist()
  } else{
    addedParms <- NULL
  }
  
  # Substitute model into fit function
  fitFunctionWithModel <- gsub("MODEL",modelTibble$Model,fitFunction)
  
  # Return
  out <- list()
  out$dataCols <- dataCols
  out$localParms <- localParms
  out$globalParms <- globalParms
  out$addedParms <- addedParms
  out$fitFunction <- fitFunctionWithModel
  out$displayedValue <- modelTibble$Displayed
  out$displayedUnits <- modelTibble$DisplayedUnits
  
  return(out)
}

# Calculate reference radius for each sector
calculateReferenceRadius <- function(x,savedSectors,allScans,fitParms){
  
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
    tryCatch(updateLog(output,"<b>Fatal error:</b> data mismatch. Sector not properly subset."),error=function(e)return())
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
processFitFunction <- function(fitParms,globalData){
  
  # Get parameters from fit parms
  dataCols <- fitParms$dataCols
  globalParms <- fitParms$globalParms
  localParms <- fitParms$localParms
  
  # Collapse global parms
  globalParmsList <- paste(globalParms,collapse=',')
  
  # Collapse data columns into a vector for insertion into dynamic fit string
  dataColumnList <- paste(dataCols,collapse=',')
  
  # Collapse global parms into vector for start 
  globalParmListStart <- paste(globalParms,'=1',sep='',collapse=',')
  
  # Get length of files to determine length of local variable vectors
  parmArrayLength <- length(dataList$savedSectors)
  
  # Create local variables and collapse into vectors for data and start
  localParmList <- lapply(localParms,function(x){
    vec<-rep(x,parmArrayLength)
    new <- paste(vec,'_',1:parmArrayLength,sep="")
    return(new)
  }) %>% unlist()
  localParmListCollapsed <- paste(localParmList,collapse=",")
  localParmListStart <- paste(localParmList,'=1',sep='',collapse=',')
  
  # Format start string
  startString <- paste('c(',globalParmListStart,',',localParmListStart,")",sep='')
  
  # Substitute parameters into fit function
  dynamicFitFunction <- gsub("LOCAL_FIT_BOILERPLATE",localFitBoilerplate,fitParms$fitFunction)
  dynamicFitFunction <- gsub('DATA_COLUMNS',dataColumnList,dynamicFitFunction)
  dynamicFitFunction <- gsub('LOCAL_PARAMETERS',localParmListCollapsed,dynamicFitFunction)
  dynamicFitFunction <- gsub('GLOBAL_PARAMETERS',globalParmsList,dynamicFitFunction)
  
  # Load into real function
  loadedFunction <- eval(parse(text=dynamicFitFunction))
  .GlobalEnv$loadedFunction <- loadedFunction
  
  # Format command to run the fit
  fitString <- "
        gsl_nls(
          Ar~loadedFunction(DATA_COLUMNS,GLOBAL_PARAMETERS,LOCAL_PARAMETERS),
          data=globalData,
          algorithm='lm',
          start=START_STRING,
          control=gsl_nls_control(maxiter=1000)
        )
      "
  
  # Insert parms into command
  dynamicFitString <- gsub('DATA_COLUMNS',dataColumnList,fitString)
  dynamicFitString <- gsub('GLOBAL_PARAMETERS',globalParmsList,dynamicFitString)
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
plotNonlinearFit <- function(output,globalFitDup){
  
  nonlinear <- ggplot(globalFitDup,aes(x=r))+
    geom_point(aes(y=Ar-offset))+
    geom_line(aes(y=pAr-offset,group=ID),col='red',lwd=1)+
    theme_prism()
  
  tryCatch(print(nonlinear),
           warning=function(w){
             warningString <- paste("<b>Warning</b> in nonlinear plotting: ",w[[1]],sep="")
             updateLog(output,warningString)
             return()
           })
  
  return(nonlinear)
  
}

# Create function to plot residuals
plotNonlinearResiduals <- function(output,globalFitDup){
  
  # Calculate residuals
  globalFitDup$residuals <- globalFitDup$Ar - globalFitDup$pAr
  
  # Get residuals vector ordered by r-r0 to mirror plot
  resVec <- globalFitDup$residuals
  
  # Check residuals for outliars
  resH <- IQR(resVec)*1.5
  resMin <- quantile(resVec,0.25)[[1]] - resH
  resMax <- quantile(resVec,0.75)[[1]] + resH
  
  # Label outliars
  outliarInds <- which(resVec>resMax|resVec<resMin)
  globalFitDup$Outliar <- "F"
  if(length(outliarInds)>0){
    globalFitDup[outliarInds,]$Outliar <-'T' 
    
    # Warn user
    outliarString <- paste("<b>Warning</b>: ",length(outliarInds)," outliar residuals detected. Consider removing outliars to improve fitting.",sep="")
    tryCatch(updateLog(output,outliarString),error=function(e)return())
  }
  
  # Order resvec for distribution test
  resVecOrdered <- resVec[order((globalFitDup$r-globalFitDup$r0),decreasing=FALSE)]
  
  # Split residuals into groups of 10
  chunkLength <- 10
  groups <- split(resVecOrdered,ceiling(seq_along(resVecOrdered) / chunkLength) )
  
  # Remove groups with too few points
  lengths <- lapply(groups,length) %>% unlist()
  groupsCur <- groups[!lengths<chunkLength]
  
  # For each group, get number of points above 0 (should be 0.5)
  pointsAboveZero <- lapply(groupsCur,function(x){length(which(x>0))}) %>% unlist()
  fractionAboveZero <- pointsAboveZero/chunkLength
  
  # Simulate what value should be with uniform distirbution
  # uniformRes <- runif(10000,-1,1)
  # uniformResGroups <- split(uniformRes,ceiling(seq_along(uniformRes)/chunkLength))
  # unifAboveZero <- lapply(uniformResGroups,function(x){length(which(x>0))}) %>% unlist()
  # unifFractionAboveZero <- unifAboveZero/chunkLength
  # # Calculate threshholds from simulation
  # zScore <- 1.96
  # maxThresh <- mean(unifFractionAboveZero) + zScore*sd(unifFractionAboveZero) # about 0.8
  # minThresh <- mean(unifFractionAboveZero) - zScore*sd(unifFractionAboveZero) # about 0.2
  
  # Calculate how many points are significantly greater or less than 0.5
  fractionHigh <- length(which(fractionAboveZero>0.8))/length(fractionAboveZero)
  fractionLow <- length(which(fractionAboveZero<0.2))/length(fractionAboveZero)
  
  # Sum should be less than 0.05, doing 0.1 for less stringent thresh
  sumOfFractions <- fractionLow + fractionHigh
  
  if(sumOfFractions>0.1){
    improperResidualDist<- "<b>Warning</b>: Residuals do not appear randomly distributed around the fit. This suggests an invalid fit."
    tryCatch(updateLog(output,improperResidualDist),error=function(e)return())
  }
  
  residualPlot <- ggplot(globalFitDup,aes(x=r-r0,y=residuals,col=as.character(ID)))+
    geom_point(shape=ifelse(globalFitDup$Outliar=="T",1,16))+
    geom_hline(yintercept=0,col='red')+
    ylab("Ar - pAr")+
    theme_prism()+
    theme(legend.position = 'none')
  
  return(residualPlot)
  
}

# Process fit results
processFitResults <- function(globalFit,globalData,processedFitFunction,fitParms){
  
  localParms <- fitParms$localParms
  globalParms <- fitParms$globalParms
  
  # Extract fit summary
  fitSummary <- summarizeNLS(globalFit,dataList$selectedModel)
  parmRownames <- rownames(fitSummary$coefs)
  
  # Create copy of df to avoid modification
  globalFitDup <- globalData
  
  # Update log
  tryCatch(updateLog(output,fitSummary$summary),error=function(e){return()})
  
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
  
  # Save to global env for debugging
  .GlobalEnv$globalFitDup <- globalFitDup
  
  # Plot fits
  nonlinear <- plotNonlinearFit(output,globalFitDup)
  residuals <- plotNonlinearResiduals(output,globalFitDup)
  
  # Get global parm fit coefs
  globalCoef <- fitSummary$coefs[which(rownames(fitSummary$coefs)==fitParms$globalParms),]
  
  # Calculate displayed value
  displayedValue <- fitParms$displayedValue %>%
    strsplit(',') %>%
    unlist() %>%
    lapply(calculateAdditionColumns,coefs=fitSummary$coefs) %>%
    unlist()
  
  # Propagate error
  displayedValueError <- (globalCoef$`Std. Error`/globalCoef$Estimate)*as.numeric(displayedValue)
    
  # Format displayed text
  displayedText <- paste(names(displayedValue),' = ',round(as.numeric(displayedValue),1)," +/- ",round(displayedValueError,1),sep="")

  # Generate row to display fit coefficients
  parmRow <- div(
    lapply(
      displayedText,
      function(x){
        p(HTML(x))
      }
    ),
    renderUI({
      div(
        style='background-color:white',
        HTML(tab_model(globalFit,p.style='stars',title=dataList$selectedModel,dv.labels = '',string.pred="Parameters",show.ci=FALSE,show.se=TRUE)$knitr)
      )
      })
    #renderTable(fitSummary$coefs,rownames = TRUE)
  )
  
  
  
  
  
  # Return out list
  outlist <- list()
  outlist$nonlinear <- nonlinear
  outlist$data <- globalFitDup
  outlist$description <- parmRow
  outlist$residuals <- residuals
  
  return(outlist)
  
}

# Function to render the fit output
renderFitOutput <- function(output){
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
        div(class='row',style=paste('overflow-y: auto;height:',fitPlotHeight*0.7,'px;margin-left:5px',sep=""),
            fit$description
        ),
        div(
          class='row',style='margin-top:15px;align-self:flex-end',
          p('Download',class='header')
        ),
        div(
          class='row',style='align-self:flex-end;',
          div(style='margin-right:10px;',
            downloadButton("downloadFit","Fit data",icon=icon('database'))
          ),
          downloadButton("downloadExperiment","Experiment file",icon=icon('flask'))
        ),
      )
    )
  })
}


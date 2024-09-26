
# Render the model type to display
renderModelTypeSelection<-function(output,model){
  output$modelTypeSelection <- renderUI({
    selectInput('modelType',NULL,choices=modelTypes,selected=model,multiple=FALSE,width=300)
  })
}

# Create function to check for proper inputs before fitting
checkForInputs <- function(input,output){
  
  if(length(dataList$selectedModel)==0){
    selectedModel <- input$modelType
    updateDataList('selectedModel',selectedModel)
  }
  
  # Get psv and sd inputs
  
  # Format check vector
  checkVec <- c(input$psvInput,input$sdInput)
  if(length(checkVec)==0){
    updateLog(output,"Invalid partial specific volume (PSV) input. Defaulting to 0.71 mL/g")
    psv <- 0.71
    updateLog(output,"Invalid solvent density (SD). Defaulting to 1.003 g/mL")
    sd <- 1.003
  } else{
    
    if(input$psvInput==""){
      updateLog(output,"Invalid partial specific volume (PSV). Defaulting to 0.71 mL/g")
      psv <- 0.71
    } else{
      psv <- as.numeric(input$psvInput)
    }
    
    if(input$sdInput==""){
      updateLog(output,"Invalid solvent density (SD) input. Defaulting to 1.003 g/mL")
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
    updateLog(output,"No sectors to fit. Save sectors to fit a model.")
    return(FALSE)
  }
  
  return(TRUE)
  
}

# Define function parameters
defineFitParms <- function(input){
  
  if(dataList$selectedModel=='MW / Single ideal species'){
    
    dataCols <- c("r","r0","w","temp","Ar","R","ID")
    localParms <- c("A0","offset")
    globalParms <- c("Mb")
    lowerParms <- c('Mb'=0)
    upperParms <- c("Mb"=10000)
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
    mwValue <- tryCatch(as.numeric(input$mwInput),error=function(e)return(as.numeric(dataList$mwValue)))
    mwValue <- ifelse(is.na(mwValue),80000,mwValue)
    eCoefValue <- tryCatch(as.numeric(input$eCoefInput),error=function(e)return(as.numeric(dataList$eCoefValue)))
    eCoefValue <- ifelse(is.na(eCoefValue),76000,eCoefValue)
    nValue <- tryCatch(as.numeric(input$nInput),error=function(e)return(as.numeric(dataList$nValue)))
    nValue <- ifelse(is.na(nValue),2,nValue)
    
    # Update datalist
    updateDataList('mwValue',mwValue)
    updateDataList('eCoefValue',eCoefValue)
    updateDataList('nValue',nValue)
    
    dataCols <- c("r","r0","w","temp","Ar","R","Mb","N","ID")
    localParms <- c("A0","offset")
    globalParms <- c("K")
    lowerParms <- c("K"=0)
    upperParms <- c("K"=100)
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
  out$lowerParms <- lowerParms
  out$upperParms <- upperParms
  
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
  
  lowerString <- paste('c("',paste(names(fitParms$lowerParms),'"=',as.numeric(fitParms$lowerParms),sep=""),")",sep="")
  
  upperString <- paste('c("',paste(names(fitParms$upperParms),'"=',as.numeric(fitParms$upperParms),sep=""),")",sep="")
  
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
          algorithm='lmaccel',
          start=START_STRING,
          control=gsl_nls_control(maxiter=1000)
        )
      "
  
  # Configure parameters
  dynamicFitString <- gsub('DATA_COLUMNS',dataColumnList,fitString)
  dynamicFitString <- gsub('GLOBAL_PARAMETERS',globalParms,dynamicFitString)
  dynamicFitString <- gsub('LOCAL_PARAMETERS',localParmListCollapsed,dynamicFitString)
  dynamicFitString <- gsub('START_STRING',startString,dynamicFitString)
  dynamicFitString <- gsub('LOWER_STRING',lowerString,dynamicFitString)
  dynamicFitString <- gsub('UPPER_STRING',upperString,dynamicFitString)
  
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

# Create function to plot linear fit
plotLinearFit <- function(output,globalFitDup,showSimulation=FALSE){
  
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
             updateLog(output,warningString)
             return()
           })
  
  return(linear)
  
}

# Process fit results
processFitResults <- function(globalFit,globalData,processedFitFunction,fitParms){
  
  localParms <- fitParms$localParms
  globalParms <- fitParms$globalParms
  
  # Extract fitted Mb
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
  
  .GlobalEnv$globalFitDup <- globalFitDup
  
  # Plot non linear fit
  nonlinear <- plotNonlinearFit(output,globalFitDup)
  residuals <- plotNonlinearResiduals(output,globalFitDup)
  
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


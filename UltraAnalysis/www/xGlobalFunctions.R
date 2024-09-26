
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

# Function to get random character
getRandomCharacter <- function(){
  id <- paste('r',sample(1:10000)[1],sep="")
  return(id)
}

# Function to get random number
getRandomNumber <- function(){
  id <- sample(1:100000)[1]
  return(id)
}

# Function to get random color
getRandomColor <- function(){
  col <- sample(colorVec,1)
  return(col)
}

# Summarize a nls as a string
summarizeNLS <- function(fit,type){
  
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

# Define colors to use
colorVec <- rainbow(10)

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

# Define fitting function
singleIdealSpecies <- function(r, w, R, temp, r0, A0, Mb, offset){
  Ar <- offset + A0 * exp(1)^( ( (Mb*w^2)/(2*R*temp) )*(r^2 - r0^2) )
  return(Ar)
}

# function to update data list
updateDataList <- function(dataType,newValue){
  
  dataList[[dataType]] <- newValue
  .GlobalEnv$dataList <- dataList
  
  
}

# Reset the ui
resetUI <- function(output){
  
  # Process
  output$savedSectorUI <- NULL
  output$sectorPlotFrame <- NULL
  output$currentScanPlotOutput <- NULL
  
  # Fit
  output$fitResult <- NULL
  
  
}


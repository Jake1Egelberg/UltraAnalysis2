library(zoo)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggprism)
library(gslnls)

# Define MW function
mwFunction <- function(r, w, temp, r0, R, A0, Mb, offset){
  Ar <- offset + A0 * exp(1)^( ( (Mb*w^2)/(R*temp*2) )*(r^2 - r0^2) )
  return(Ar)
}

# Define universal gas constant in proper units
# Correct gas constant to units of (g * cm^2) / s^2 * mol * K
R <- 8.3144 * 1000 * 100 * 100 # (g * cm^2) / s^2 * mol * K

selectedCells <- 'C:/1_Documents/Ferguson Lab/UltraAnalysis/Data/618 WT EGFR/cell4abs.log'
selectedCells <- 'C:/Users/Jake/Documents/Code/UltraAnalysis2/Data/cell4abs.log'

# Read selected log file
readLogFile <- function(selectedCells){
  
  currentLog <- selectedCells
  
  # Get log directory
  logDirectory <- str_sub(currentLog,end=max(str_locate_all(currentLog,'\\/')[[1]]))
  
  # Read log file
  logLines <- readLines(currentLog,warn=FALSE)
  
  # Remove blanks
  blankInds <- which(logLines=="")
  if(length(blankInds)>0){
    logLinesNoBlanks <- logLines[-blankInds]
  } else{
    logLinesNoBlanks <- logLines
  }
  
  out <- list()
  out$Log <- logLinesNoBlanks
  out$Directory <- logDirectory
  
  return(out)
  
}
logFile <- readLogFile(selectedCells)

# Extract information for each scan
extractScanInformation <- function(logFile){
  
  log <- logFile$Log
  dir <- logFile$Directory
  
  # Function to get log information
  parseScan <- function(x,dir){
    
    # Spit data into components
    splitData <- trimws(strsplit(log[x],",")[[1]])
    
    if(splitData %>% length() < 4 ){
      return(NULL)
    }
    
    # Get file name from prev line
    prevLine <- log[x-1]
    fileName <- trimws(str_sub(prevLine,start=str_locate(prevLine,'file name:')[[2]]+1))
    
    if(is.na(fileName)==TRUE){
      return(NULL)
    }
    
    # Curate to full file name assuming in same directory as log
    fullFileName <- paste(dir,fileName,sep="")
    
    # Extract indices of colon
    colonInds <- str_locate(splitData,":")[,1]
    
    # Extract column names and data values
    dataNames <- str_sub(splitData,end=colonInds-1)
    dataValues <- as.numeric(str_sub(splitData,start=colonInds+1))
    
    # Format as a matrix
    mat <- matrix(ncol=length(dataNames),nrow=1)
    colnames(mat) <- dataNames
    mat[1,] <- dataValues
    
    # Covnert matrix to df
    df <- as.data.frame(mat)
    df$file <- fullFileName
    return(df)
  }
  # Get log info
  out <- lapply(1:length(log),parseScan,dir=dir)
  
  # Combine to df
  logDf <- bind_rows(out)
  
  # Add cell number, scan number, calculate w from w2t
  cellNumber <- as.numeric(paste(str_extract_all(str_extract(log[1],'Cell [:digit:]'),'[:digit:]')[[1]],collapse=""))
  logDf$cell <- cellNumber
  logDf$scan <- 1:nrow(logDf)
  
  # Calculate omega from omega^2 * t
  logDf$w <- sqrt((logDf$w2t/logDf$time))
  
  # Convert temperature from C to K
  logDf$temp <- as.numeric(logDf$temp)+273.15
  
  # Reformat output
  newFormat <- data.frame(
    Cell = logDf$cell,
    Scan = logDf$scan,
    Speed_rpm = logDf$speed,
    W_radPers = logDf$w,
    Temp_K = logDf$temp,
    File = logDf$file
  )
  
  # Return
  return(newFormat)
  
}
extractedScanInformation <- extractScanInformation(logFile)

# Read each scan
readScans <- function(extractedScanInformation){
  
  # Read each scan
  readAScan <- function(x,extractedScanInformation){
    
    row <- extractedScanInformation[x,]
    
    # Read data
    scanLines <- read.table(row$File,header=TRUE,sep="\t")
    
    # Get scan name
    scanName <- names(scanLines)
    
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
    radialPosition <- as.numeric(unlist(lapply(noBlanks,'[[',1)))
    absorbance <- unlist(lapply(noBlanks,'[[',2))
    noise <- unlist(lapply(noBlanks,'[[',3))
    
    # Create df
    dfValues <- data.frame(
      Name = NA,
      Cell = row$Cell,
      Scan = row$Scan,
      r = radialPosition,
      r0 = NA,
      w = (row$Speed_rpm/60)*2*pi,
      temp = row$Temp_K,
      Ar = absorbance,
      R=R
    ) %>% mutate_all(as.numeric)
    dfValues$Name <- scanName
    dfValues$File <- row$File
    
    
    naInds <- which(is.na(dfValues$Ar))
    if(length(naInds)>0){
      curatedValues <- dfValues[-naInds,]
    } else{
      curatedValues <- dfValues
    }
    
    return(dfValues)
    
  }
  out <- lapply(1:nrow(extractedScanInformation),readAScan,extractedScanInformation=extractedScanInformation)
  
  # Name scan list with condition
  cellScanSpeed <- paste("Cell ",extractedScanInformation$Cell," Scan ",extractedScanInformation$Scan," ",extractedScanInformation$Speed_rpm," rpm",sep="")
  names(out) <- cellScanSpeed
  
  return(out)
  
}
scanData <- readScans(extractedScanInformation)

# Create a list of scans to analyze
s1 <- scanData$`Cell 4 Scan 1 6000 rpm`
s1sub <- subset(s1,r>5.9&r<6.1)
s2 <- scanData$`Cell 4 Scan 2 6000 rpm`
s2sub <- subset(s2,r>5.9&r<6.1)


# Format list
rawscanList <- list(s1sub,
                    s2sub
                    )

# calculate reference radii
scanList <- lapply(rawscanList,function(x){
  tmp <- x
  tmp$r0 <- (diff(range(tmp$r))*(2/3))+min(tmp$r)
  return(tmp)
})

# Using list length, create parameter array
parmArrayLength <- length(scanList)

# Number each analysis and collapse list to df
globalData <- lapply(1:parmArrayLength,function(x){
  tmp <- scanList[[x]]
  tmp$ID <- x
  return(tmp)
}) %>% bind_rows()
plot(globalData$r,globalData$Ar)

# Select fit data
fitData <- globalData %>% select(
  c(
    "r",
    "r0",
    "w",
    "temp",
    "Ar",
    "R",
    "ID"
  )
)

modelType <- c("MW / Single ideal species",
               "Kd / A + A <-> AA")[2]

# Define local fit boilerplate
localFitBoilerplate <- "
  A0 <- lapply(paste('A0',ID,sep=''),function(x){get(x)}) %>% unlist()
  offset <- lapply(paste('offset',ID,sep=''),function(x){get(x)}) %>% unlist()
"


if(modelType=='MW / Single ideal species'){
  
  dataCols <- c("r","r0","w","temp","Ar","R","ID")
  localParms <- c("A0","offset")
  globalParms <- c("Mb")
  
  fitFunction <- "
    function(DATA_COLUMNS,GLOBAL_PARAMETERS,LOCAL_PARAMETERS){
      LOCAL_FIT_BOILERPLATE
      
      Ar <- (offset) + ( (A0) * exp(1)^( ( (Mb*w^2)/(R*temp*2) )*(r^2 - r0^2) ) )
      
      return(Ar)
    }
  "
  
  fitData$Mb <- NA
  
} else if(modelType=='Kd / A + A <-> AA'){
  
  dataCols <- c("r","r0","w","temp","Ar","R","Mb","N","ID")
  localParms <- c("A0","offset")
  globalParms <- c("K")
  
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
  
  fitData$Mb <- 24469
  fitData$N <- 2
  
}

# Process parms

dataColumnList <- paste(dataCols,collapse=',')

globalParmList <- globalParms
globalParmListStart <- paste(globalParmList,'=1',sep='',collapse=',')

localParmList <- lapply(localParms,function(x){
  vec<-rep(x,parmArrayLength)
  new <- paste(vec,1:length(vec),sep="")
  return(new)
}) %>% unlist()
localParmListCollapsed <- paste(localParmList,collapse=",")
localParmListStart <- paste(localParmList,'=1',sep='',collapse=',')

startString <- paste('c(',globalParmListStart,',',localParmListStart,")",sep='')

# Substitute parameters into fit function
dynamicFitFunction <- gsub("LOCAL_FIT_BOILERPLATE",localFitBoilerplate,fitFunction)
dynamicFitFunction <- gsub('DATA_COLUMNS',dataColumnList,dynamicFitFunction)
dynamicFitFunction <- gsub('LOCAL_PARAMETERS',localParmListCollapsed,dynamicFitFunction)
dynamicFitFunction <- gsub('GLOBAL_PARAMETERS',globalParmList,dynamicFitFunction)

# Load into real function

loadedFunction <- eval(parse(text=dynamicFitFunction))

# Format the fit string
fitString <- "
  gsl_nls(
    Ar~loadedFunction(DATA_COLUMNS,GLOBAL_PARAMETERS,LOCAL_PARAMETERS),
    data=fitData,
    algorithm='lm',
    start=START_STRING,
    control=gsl_nls_control(maxiter=1000)
  )
"

# Configure parameters
dynamicFitString <- gsub('DATA_COLUMNS',dataColumnList,fitString)
dynamicFitString <- gsub('GLOBAL_PARAMETERS',globalParmList,dynamicFitString)
dynamicFitString <- gsub('LOCAL_PARAMETERS',localParmListCollapsed,dynamicFitString)
dynamicFitString <- gsub('START_STRING',startString,dynamicFitString)

# Do fit
globalFit <- eval(parse(text=dynamicFitString))

# Get summary
sum <- summary(globalFit)

# Plot predictions
fitData$pAr <- predict(globalFit,fitData)
ggplot(fitData,aes(x=r))+
  geom_point(aes(y=Ar))+
  geom_line(aes(y=pAr,group=ID),col='red',lwd=1)+
  theme_prism()

sum




readScans <- function(file,filenames){
  
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
      theta = w^2 / (2*R*temp),
      Noise=noise
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
# Process kyles data
paths <- c(
  'C:/1_Documents/Ferguson Lab/AUC Data/TGFa - Copy/00005.RA2',
  'C:/1_Documents/Ferguson Lab/AUC Data/TGFa - Copy/00011.RA2',
  'C:/1_Documents/Ferguson Lab/AUC Data/TGFa - Copy/00015.RA2'
)

scan5 <- readScans(paths[1],'00005.RA2')[[1]]
scan11 <- readScans(paths[2],'00011.RA2')[[1]]
scan15 <- readScans(paths[3],'00015.RA2')[[1]]

# Subset out ranges
scan5sector1 <- subset(scan5,r>=5.924&r<=6.057)
scan5sector2 <- subset(scan5,r>=6.404&r<=6.513)
scan5sector3 <- subset(scan5,r>=6.94&r<=7.068)

scan11sector1 <- subset(scan11,r>=5.967&r<=6.094)
scan11sector2 <- subset(scan11,r>=6.409&r<=6.526)
scan11sector3 <- subset(scan11,r>=6.918&r<=7.059)

scan15sector1 <- subset(scan15,r>=6.367&r<=6.499)
scan15sector2 <- subset(scan15,r>=6.985&r<=7.115)

# Create secotr list
rawscanList <- list(
  scan5sector1,
  scan5sector2,
  scan5sector3,
  scan11sector1,
  scan11sector2,
  scan11sector3,
  scan15sector1,
  scan15sector2
)

# Save each file as a .RA2 file for loading in hetero to check
lapply(rawscanList,function(x){
  
  name <- "TGFa"
  
  if(unique(x$AbsoluteScan)==5){
    headerString <- 'R 2 20.0 04000 0067906 1.1970E10 496 5' 
  } else if(unique(x$AbsoluteScan)==11){
    headerString <- 'R 2 20.0 06000 0131548 3.5150E10 496 5'
  } else if(unique(x$AbsoluteScan)==15){
    headerString <- 'R 2 20.0 09000 0186174 8.0320E10 480 5'
  }
  
  # Get data
  dataToSave <- select(x,c('r','Ar','Noise'))
  dataToSave$Blank <- ''
  dataToSave <- dataToSave[,c(4,1,2,3)]
  
  # Save
  newFile <- gsub('.RA2','_sub.RA2',unique(x$File))
  
  write.table(
    name,
    newFile,
    sep='\t',
    quote=FALSE,
    row.names = FALSE,
    col.names=FALSE
  )
  
  write.table(
    headerString,
    newFile,
    sep='\t',
    quote=FALSE,
    row.names = FALSE,
    col.names=FALSE,
    append=TRUE
  )
  
  write.table(
    dataToSave,
    newFile,
    sep='\t',
    quote=FALSE,
    row.names = FALSE,
    col.names=FALSE,
    append=TRUE
  )
  
})






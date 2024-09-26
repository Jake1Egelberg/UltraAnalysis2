
# Function to update the log
updateLog <- function(output,update){
  
  logText <- c(update,dataList$logText)
  updateDataList('logText',logText)
  
  output$log <- renderUI({
    source('www/log.R')[[1]]
  })
  
}
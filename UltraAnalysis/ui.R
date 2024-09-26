library(shiny)

# Data handling
library(readr)
library(stringr)
library(dplyr)

# Plotting
library(ggplot2)
library(ggprism)
library(cowplot)

# Nonlinear fits
library(gslnls)

# Define UI for application that draws a histogram
fluidPage(

  tags$body(
    style='background-color:#DCEAF7;overflow-y:hidden'
  ),
  tags$link(rel = "stylesheet", type = "text/css", href = "xStyles.css"), # link external stylesheet 
  p('UltraAnalysis V1.0.0',style='position:absolute;right:2vw;top:3vh;z-index:100;font-size:20px;',class='header'),
  tabsetPanel(
    id='tabSwitch',
    tabPanel(
      value='upload',
      title='📂 Upload',
      div(
        style='margin-top:12vh;',
        uiOutput('upload')
      )
    ),
    tabPanel(
      value='process',
      title='📝 Process',
      div(
        style='margin-top:12vh',
        uiOutput('process')
      )
    ),
    tabPanel(
      value='fit',
      title="📈 Fit",
      div(
        style='margin-top:12vh',
        uiOutput('fit')
      )
    ),
    selected='upload'
  ),
  uiOutput('log')
)


div(
  class='row',style='height:65vh;',
  div(
    class='column',style='width:19%;height:100%;padding:10px;',
    p('Included scans',class='header'),
    selectInput('selectScan',NULL,choices=scansToAnalyze,multiple=FALSE,width=300),
    
    div(class='row',style='margin:2px',
      div(
        class='column',
        p('Saved sectors',class='header')
      ),
      div(
        class='column',style='margin:0 0 0 auto',
        actionButton('undoSector',"⬅️ Undo",style='font-size:10px')
      )
    ),
    
    div(
      style='margin:5px;height:200px;overflow-y:scroll',
      uiOutput('savedSectorUI')
    )

  ), # Column 1
  div(
    class='column',style='width:1%;height:90%;border-left: 1.5px solid #002060;margin-left:1%;align-self:center',
  ), # divider column
  div(
    class='column',style='width:fit-content;height:100%;padding:10px;',
    p('Select a sector',class='header'),
    plotOutput('currentScanPlotOutput',brush = 'selectingSector')
  ), #column 2
  div(
    class='column',style='padding:10px;',
    p('Sector preview',class='header'),
    uiOutput('sectorPlotFrame')
  ) #column 3
  
  
) # Parent div

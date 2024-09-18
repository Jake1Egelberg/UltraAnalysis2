div(
  class='row',style='height:65vh;',
  div(
    class='column',style='width:19%;height:100%;padding:10px;',
    p('Included scans',class='header'),
    uiOutput('scanSelectionInput'),
    
    p('Saved sectors',class='header'),
    
    div(
      style='margin:5px;height:200px;overflow-y:scroll',
      uiOutput('savedSectorUI')
    )

  ), # Column 1
  div(
    class='column',style='width:1%;height:90%;border-left: 1.5px solid #002060;margin-left:1%;align-self:center',
  ), # divider column
  div(
    class='column',style='width:fit-content;height:100%;padding:10px',
    div(
      class='row',style='padding-right:15px;padding-left:15px',
      p('Select a sector',class='header'),
      actionButton('autoSector',label='🔍 Auto',style='font-size:10px;margin-right:10px;margin-left:10px;padding:5px;height:fit-content'),
      p('Drag over an area to select it, double click to reset the zoom',class='subbody',style='align-self:flex-end; margin:0 0 10px auto'),
    ),
    plotOutput('currentScanPlotOutput',brush = 'selectingSector', dblclick='resetPlotView')
  ), #column 2
  div(
    class='column',style='padding:10px;',
    p('Sector preview',class='header'),
    uiOutput('sectorPlotFrame')
  ) #column 3
  
  
) # Parent div

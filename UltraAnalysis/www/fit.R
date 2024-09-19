
div(
  class='row',style='height:65vh;',
  div(
    class='column',style='width:19%;height:100%;padding:10px;',
    p('Model to fit',class='header'),
    uiOutput('modelTypeSelection'),
    uiOutput('modelParms')
  ), # Column 1
  div(
    class='column',style='width:1%;height:90%;border-left: 1.5px solid #002060;margin-left:1%;align-self:center',
  ), # divider column
  div(
    class='column',style='padding:10px;',
    textOutput('fitDescription')%>%tagAppendAttributes(class = 'header'),
    uiOutput('fitResult')
  ), #column 2
  
  
) # Parent div

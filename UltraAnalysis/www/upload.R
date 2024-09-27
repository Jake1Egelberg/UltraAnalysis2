
div(
  class='row',style='height:65vh;',
  div(
    class='column',style='width:19%;height:100%;padding:10px;',
    p('Upload scans (.RA) or a UA experiment file (.Rdata)',class='header'),
    fileInput('uploadInput',NULL,TRUE,accept=c(paste(".RA",1:10,sep=""),'.Rdata'),width=300,buttonLabel = '🗉 Select file(s)',placeholder = '0 file(s)'),
    uiOutput('uploadInformation')
  ), # Column 1
  div(
    class='column',style='width:1%;height:90%;border-left: 1.5px solid #002060;margin-left:1%;align-self:center',
  ), # divider column
  div(
    class='column',style='width:80%;height:100%;padding:10px;',
    p('Scan previews',class='header'),
    div(style='height:100%;overflow-y:scroll;margin-right:5px;',
      uiOutput('dataPreview')
    ) # inner box
  ) #column 2
  
  
) # Parent div

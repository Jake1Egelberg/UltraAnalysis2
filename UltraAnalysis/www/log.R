
div(
  class='column',style='
  background-color:#d6d6d6;
  width:100%;
  padding:10px;
  border-top: 1.5px solid black;
  overflow-y:scroll',id='logDiv',
  lapply(dataList$logText,function(text){
    p(HTML(text), class='log')
  })
)
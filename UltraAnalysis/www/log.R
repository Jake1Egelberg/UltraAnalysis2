
div(
  class='column',style='background-color:#d6d6d6;width:100%;height:20vh;padding:10px;border-top: 1.5px solid black;overflow-y:scroll',
  lapply(logText,function(text){
    p(text, class='log')
  })
)
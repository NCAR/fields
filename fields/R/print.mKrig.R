print.mKrig<- function (x, digits = 4,...){
  object<- summary.mKrig( x,...)
  print.mKrigSummary(object, digits = digits)
}

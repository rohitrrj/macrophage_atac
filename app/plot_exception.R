#' plot message for exception
#' 
#' Useful to display messages in \code{shiny} reports
#'
#' Typically call \code{return(plot_exception(...))} where you would have called \code{stop(...)}
#' @param ... text to display, concatenated with sep
#' @param sep separator used for concatenation
#' @param type function to use to print in console
#' @param color text color, by default red for message and warning else black
#' @param console if TRUE print in console, if FALSE just plot
#' @param size text size
#' @examples
#' plot_exception("no data for current filter selection")
#' plot_exception("NO","WAY","!!!",color="blue",size=12,console=FALSE)
#' @export
plot_exception <-function(
  ...,
  sep=" ",
  type=c("message","warning","cat","print"),
  color="auto",
  console=TRUE,
  size = 6){      
  type=match.arg(type)
  txt = paste(...,collapse=sep)
  if(console){
    if(type == "message") message(txt)
    if(type == "warning") warning(txt)
    if(type == "cat") cat(txt)
    if(type == "print") print(txt)
  }
  if(color =="auto") color <- if(type == "cat") "black" else "red"
  if(txt == "warning") txt <- paste("warning:",txt)
  print(ggplot2::ggplot() +
          ggplot2::geom_text(ggplot2::aes(x=0,y=0,label=txt),color=color,size=size) + 
          ggplot2::theme_void())
  invisible(NULL)
}
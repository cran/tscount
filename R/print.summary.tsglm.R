print.summary.tsglm <- function(x, ...){ 
  if(length(coef(x)) > 0){
    cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    cat("Coefficients:\n")
    print(format.data.frame(as.data.frame(coef(x)), digits=3), print.gap=2, quote=FALSE, na.print="")
    if(!is.null(coef(x)$"Std.Error")){
      if(x$se.type == "normapprox") cat("Standard errors and confidence intervals (level = ", x$level*100, "%) obtained\nby normal approximation.\n")
      if(x$se.type == "bootstrap") cat("Standard errors and confidence intervals (level = ", x$level*100, "%) obtained\nby parametric bootstrap with", x$se.bootstrapsamples, "replications.\n")
    }
    cat(
      "\nLink function:", x$link,
      "\nDistribution family:", x$distr, if(x$distr=="nbinom"){"(with overdispersion coefficient 'sigmasq')"}else{NULL},
      "\nNumber of coefficients:", x$number.coef,
      "\nLog-likelihood:", x$logLik,
      "\nAIC:", x$AIC,
      "\nBIC:", x$BIC,
      "\nQIC:", x$QIC,
    "\n\n")
  }else{ 
    if(length(x$init)>0){
      print(x, ...)
    }else{
      cat("No coefficients\n")
    }
  }
  invisible(x)  
}

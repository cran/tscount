#Link function:
g <- function(x, link=c("identity", "log")){
  link <- match.arg(link)
  result <- if(is.null(x)) NULL else switch(link,
    "identity" = x,
    "log" = log(x)
  )
  return(result)
}

#Transformation function:
trafo <- function(x, link=c("identity", "log")){
  link <- match.arg(link)
  result <- if(is.null(x)) NULL else switch(link,
    "identity" = x,
    "log" = log(x+1)
  )
  return(result)
}

#Inverse of transformation function:
trafo_inv <- function(x, link=c("identity", "log")){
  link <- match.arg(link)
  result <- if(is.null(x)) NULL else switch(link,
    "identity" = x,
    "log" = exp(x)-1
  )
  return(result)
}

#Inverse of link function:
g_inv <- function(x, link=c("identity", "log")){
  link <- match.arg(link)
  result <- if(is.null(x)) NULL else switch(link,
    "identity" = x,
    "log" = exp(x)
  )
  return(result)
}

#1st derivative of inverse link function:
g_inv_1st <- function(x, link=c("identity", "log")){
  link <- match.arg(link)
  result <- if(is.null(x)) NULL else switch(link,
    "identity" = 1,
    "log" = exp(x)
  )
  return(result)
}

#2nd derivative of inverse link function:
g_inv_2nd <- function(x, link=c("identity", "log")){
  link <- match.arg(link)
  result <- if(is.null(x)) NULL else switch(link,
    "identity" = 0,
    "log" = exp(x)
  )
  return(result)
}

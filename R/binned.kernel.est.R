#' @title Kernel Density Estimation
#' @description This function implements an interface to the kernel density estimation functions provided by the \pkg{KernSmooth} package.
#' @param data a numeric vector containing the sample on which the kernel density estimate is to be constructed.
#' @param kernel 
#' @param bandwidth
#' @param canonical
#' @param scalest
#' @param level
#' @param gridsize
#' @param range.data
#' @param truncate
#' @author Hajk-Georg Drost
#' 
#' @references
#' 
#' Matt Wand (2015). KernSmooth: Functions for Kernel Smoothing Supporting Wand & Jones (1995). R package version 2.23-14. \url{http://CRAN.R-project.org/package=KernSmooth}
#'
#' Henry Deng and Hadley Wickham (2011). Density estimation in R. \url{http://vita.had.co.nz/papers/density-estimation.pdf}.
#' @export

binned.kernel.est <- function(data, 
                              kernel     = "normal", 
                              bandwidth  = NULL, 
                              canonical  = FALSE,
                              scalest    = "minim",
                              level      = 2L,
                              gridsize   = 401L, 
                              range.data = range(data), 
                              truncate   = TRUE) {
        
        if(is.null(bandwidth)){
                bandwidth.estimate <- KernSmooth::dpik( x        = data,
                                                        scalest  = scalest,
                                                        level    = level,
                                                        kernel   = kernel,
                                                        gridsize = gridsize,
                                                        range.x  = range.data,
                                                        truncate = truncate )
        } else {
                
                bandwidth.estimate <- bandwidth
        }
        
        est <- KernSmooth::bkde(x         = data,
                                kernel    = kernel,
                                canonical = canonical,
                                bandwidth = bandwidth.estimate,
                                gridsize  = gridsize,
                                range.x   = range.data,
                                truncate  = truncate
                                )
        
        
        return(est)
}

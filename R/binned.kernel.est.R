#' @title Kernel Density Estimation
#' @description This function implements an interface to the kernel density estimation functions provided by the \pkg{KernSmooth} package.
#' @param data a numeric vector containing the sample on which the kernel density estimate is to be constructed.
#' @param kernel character string specifying the smoothing kernel
#' @param bandwidth the kernel bandwidth smoothing parameter.
#' @param canonical a logical value indicating whether canonically scaled kernels should be used
#' @param scalest estimate of scale. 
#' \itemize{
#'  \item \code{"stdev"} - standard deviation is used.
#'  \item \code{"iqr"} - inter-quartile range divided by 1.349 is used.
#'  \item \code{"minim"} - minimum of \code{"stdev"} and \code{"iqr"} is used.
#' }
#' @param level number of levels of functional estimation used in the plug-in rule.
#' @param gridsize the number of equally-spaced points over which binning is performed to obtain kernel functional approximation.
#' @param range.data vector containing the minimum and maximum values of \code{data} at which to compute the estimate. The default is the minimum and maximum data values.
#' @param truncate logical value indicating whether data with x values outside the range specified by \code{range.data} should be ignored.
#' @author Hajk-Georg Drost
#' 
#' @references
#' 
#' Matt Wand (2015). KernSmooth: Functions for Kernel Smoothing Supporting Wand & Jones (1995). R package version 2.23-14.
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

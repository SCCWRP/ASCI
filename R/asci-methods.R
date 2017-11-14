#' @include class-asci.R
NULL

#' @rdname asci-class
setMethod('show', 
          signature = 'asci', 
          definition = function(object){
            
            # defined methods from below
            meths <- .S4methods(class = 'asci')
            meths <-  attr(meths, 'info')$generic
            meths <- meths[!meths %in% 'show']

            # on screen
            cat('An object of class', class(object), '\n')
            cat('Scores calculated for', paste(object@taxa, collapse = ', ') , '\n')
            cat('Use these functions for access:', paste(meths, collapse = ', '), '\n')
            invisible(NULL)
            
          })

#' @param object \code{asci} object created with \code{\link{ASCI}}
#' @rdname asci-class
#' 
#' @export
setGeneric('scores', function(object) standardGeneric('scores'))

#' @rdname asci-class
setMethod('scores', 'asci', function(object) object@scores)
          
#' @rdname asci-class
#' @export
setGeneric('Supp1_mmi', function(object) standardGeneric('Supp1_mmi'))

#' @rdname asci-class
setMethod('Supp1_mmi', 'asci', function(object) object@Supp1_mmi)

#' @rdname asci-class
#' @export
setGeneric('Supp1_OE', function(object) standardGeneric('Supp1_OE'))

#' @rdname asci-class
setMethod('Supp1_OE', 'asci', function(object) object@Supp1_OE)

#' @rdname asci-class
#' @export
setGeneric('Supp2_OE', function(object) standardGeneric('Supp2_OE'))

#' @rdname asci-class
setMethod('Supp2_OE', 'asci', function(object) object@Supp2_OE)


setMethod('show', 
          signature = 'asci', 
          definition = function(object){
            
            # defined methods from below
            meths <- .S4methods(class = 'asci')
            meths <-  attr(meths, 'info')$generic
            meths <- meths[!meths %in% 'show']
            
            # unique samples 
            n <- scores(object)
            n <- unique(n$SampleID)
            n <- length(n)
            
            # on screen
            cat('An object of class', class(object), '\n')
            cat('Scores calculated for', paste(object@taxa, collapse = ', ') , 'indices for', n, 'unique samples\n')
            cat('Use these functions for access:', paste(meths, collapse = ', '), '\n')
            invisible(NULL)
            
          })


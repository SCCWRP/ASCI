# methods for asci class

# show method
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

##
# scores

# generic
setGeneric('scores', function(object) standardGeneric('scores'))

# method
setMethod('scores', 'asci', function(object) object@scores)

##
# Supp1_mmi

# generic
setGeneric('Supp1_mmi', function(object) standardGeneric('Supp1_mmi'))

# method
setMethod('Supp1_mmi', 'asci', function(object) object@Supp1_mmi)

##
# Supp1_OE

# generic
setGeneric('Supp1_OE', function(object) standardGeneric('Supp1_OE'))

# method
setMethod('Supp1_OE', 'asci', function(object) object@Supp1_OE)

##
# Supp2_OE

# generic
setGeneric('Supp2_OE', function(object) standardGeneric('Supp2_OE'))

# method
setMethod('Supp2_OE', 'asci', function(object) object@Supp2_OE)

##
# Supp3_OE

# generic
setGeneric('Supp3_OE', function(object) standardGeneric('Supp3_OE'))

# method
setMethod('Supp3_OE', 'asci', function(object) object@Supp3_OE)
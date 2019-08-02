##- interactionMatrix
setGeneric(name = "interactionMatrix",
            def = function(object) {
                standardGeneric("interactionMatrix")
            }
)

##- DIR
setGeneric(name = "DIR",
            def = function(object, pvalue) {
                standardGeneric("DIR")
            }
)

##- compartments
setGeneric(name = "compartments",
            def = function(object) {
                standardGeneric("compartments")
            }
)

##- concordances
setGeneric(name = "concordances",
            def = function(object) {
                standardGeneric("concordances")
            }
)

##- parameters
setGeneric(name = "parameters",
            def = function(object) {
                standardGeneric("parameters")
            }
)

setGeneric(name = "parameters<-",
            def = function(object, value) {
                standardGeneric("parameters<-")
            }
)

##- parameters
setGeneric(name = "printHiCDOCParameters",
            def = function(object) {
                standardGeneric("printHiCDOCParameters")
            }
)

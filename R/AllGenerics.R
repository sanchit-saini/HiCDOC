## - chromosomes
setGeneric(
    name = "chromosomes",
    def = function(object) {
        standardGeneric("chromosomes")
    }
)

## - positions
setGeneric(
    name = "positions",
    def = function(object) {
        standardGeneric("positions")
    }
)

## - conditions
setGeneric(
    name = "conditions",
    def = function(object) {
        standardGeneric("conditions")
    }
)

## - replicates
setGeneric(
    name = "replicates",
    def = function(object) {
        standardGeneric("replicates")
    }
)

## - interactions
setGeneric(
    name = "interactions",
    def = function(object) {
        standardGeneric("interactions")
    }
)

## - resolution
setGeneric(
    name = "resolution",
    def = function(object) {
        standardGeneric("resolution")
    }
)

## - compartments
setGeneric(
    name = "compartments",
    def = function(object) {
        standardGeneric("compartments")
    }
)

## - differences
setGeneric(
    name = "differences",
    def = function(object, threshold = 0.05) {
        standardGeneric("differences")
    }
)

## - concordances
setGeneric(
    name = "concordances",
    def = function(object) {
        standardGeneric("concordances")
    }
)

## - parameters
setGeneric(
    name = "parameters",
    def = function(object) {
        standardGeneric("parameters")
    }
)

setGeneric(
    name = "parameters<-",
    def = function(object, value) {
        standardGeneric("parameters<-")
    }
)

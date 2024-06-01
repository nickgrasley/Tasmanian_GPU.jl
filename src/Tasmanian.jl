module Tasmanian

import Base: show, read, run, write
using Plots
using Random
using SparseArrays
using Tasmanian_jll

const global TASlib = libtasmaniansparsegrid

# includes
include("libTasmanian.jl")
include("TSG.jl")
include("../examples/examples.jl")

export TasmanianSG, LocalRules, copyGrid, evaluateBatch, evaluateBatch!, getDims,
getLoadedPoints, getNeededPoints, getNout, getNumDimensions, getNumLoaded,
getNumNeeded, getNumOutputs, getNumPoints, getOrder, getPoints, isFourier,
isGlobal, isLocalPolynomial, isSequence, isWavelet, loadNeededPoints!,
makeLocalPolynomialGrid!, setDomainTransform!, setSurplusRefinement!, differentiate!

end # module

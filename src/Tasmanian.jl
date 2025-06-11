module Tasmanian

import Base: show, read, run, write
using LinearAlgebra
using Plots
using Random
using SparseArrays
using Tasmanian_jll

# Check for custom library path, fallback to jll
# If want to use custom GPU installation of Tasmanian, link the path via ENV's TASMANIAN_LIBRARY_PATH
const global TASlib = get(ENV, "TASMANIAN_LIBRARY_PATH", "") != "" ? 
    ENV["TASMANIAN_LIBRARY_PATH"] : libtasmaniansparsegrid

# includes
include("libTasmanian.jl")
include("TSG.jl")
include("../examples/examples.jl")

export TasmanianSG, LocalRules, compareGrids, copyGrid, differentiate!, enableAcceleration,
    estimateAnisotropicCoefficients, evaluate, evaluateBatch,
    evaluateBatch!, evaluateHierarchicalFunctions, evaluateSparseHierarchicalFunctions,
    evaluateThreadSafe, getAlpha, getBeta, getConformalTransformASIN, getDims,
    getDomainTransform, getHierarchicalCoefficients,
    getInterpolationWeights, getInterpolationWeightsBatch, getGPUMemory, getLevelLimits, getLoadedPoints,
    getNeededPoints, getNout, getNumDimensions, getGPUName, getNumGPUs, getNumLoaded,
    getNumNeeded, getNumOutputs, getNumPoints, getOrder, getPoints,
    getQuadratureWeights, integrate, isAccelerationAvailable, isFourier, isGlobal,
    isLocalPolynomial, isSequence, isSetDomainTransform, isSetConformalTransformASIN,
    isWavelet,
    loadNeededPoints!, loadNeededValues!, makeFourierGrid!, makeGlobalGrid!,
    makeLocalPolynomialGrid!, makeSequenceGrid!, removePointsByHierarchicalCoefficient!, setAnisotropicRefinement!,
    setConformalTransformASIN!, makeWaveletGrid!,
    setDomainTransform!, setGPUID!, setHierarchicalCoefficients!, setSurplusRefinement!, updateFourierGrid!, updateGlobalGrid!,
    updateLocalPolynomialGrid!, updateSequenceGrid!, updateWaveletGrid!


end # module

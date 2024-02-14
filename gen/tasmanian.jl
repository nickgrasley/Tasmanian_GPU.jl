module Tasmanian

using Tasmanian_jll
export Tasmanian_jll

using CEnum

# no prototype is found for this function at TasmanianSparseGrid.h:38:7, please use with caution
function tsgConstructTasmanianSparseGrid()
    ccall((:tsgConstructTasmanianSparseGrid, tasmanian), Ptr{Cvoid}, ())
end

function tsgDestructTasmanianSparseGrid(grid)
    ccall((:tsgDestructTasmanianSparseGrid, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgCopyGrid(destination, source)
    ccall((:tsgCopyGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), destination, source)
end

function tsgCopySubGrid(destination, source, outputs_begin, outputs_end)
    ccall((:tsgCopySubGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cint), destination, source, outputs_begin, outputs_end)
end

# no prototype is found for this function at TasmanianSparseGrid.h:42:13, please use with caution
function tsgGetVersion()
    ccall((:tsgGetVersion, tasmanian), Ptr{Cchar}, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:43:13, please use with caution
function tsgGetLicense()
    ccall((:tsgGetLicense, tasmanian), Ptr{Cchar}, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:44:5, please use with caution
function tsgGetVersionMajor()
    ccall((:tsgGetVersionMajor, tasmanian), Cint, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:45:5, please use with caution
function tsgGetVersionMinor()
    ccall((:tsgGetVersionMinor, tasmanian), Cint, ())
end

# no prototype is found for this function at TasmanianSparseGrid.h:46:5, please use with caution
function tsgIsOpenMPEnabled()
    ccall((:tsgIsOpenMPEnabled, tasmanian), Cint, ())
end

function tsgWrite(grid, filename)
    ccall((:tsgWrite, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}), grid, filename)
end

function tsgWriteBinary(grid, filename)
    ccall((:tsgWriteBinary, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}), grid, filename)
end

function tsgRead(grid, filename)
    ccall((:tsgRead, tasmanian), Cint, (Ptr{Cvoid}, Ptr{Cchar}), grid, filename)
end

function tsgMakeGlobalGrid(grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, alpha, beta, custom_filename, limit_levels)
    ccall((:tsgMakeGlobalGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Cdouble, Cdouble, Ptr{Cchar}, Ptr{Cint}), grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, alpha, beta, custom_filename, limit_levels)
end

function tsgMakeSequenceGrid(grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, limit_levels)
    ccall((:tsgMakeSequenceGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, dimensions, outputs, depth, sType, sRule, anisotropic_weights, limit_levels)
end

function tsgMakeLocalPolynomialGrid(grid, dimensions, outputs, depth, order, sRule, limit_levels)
    ccall((:tsgMakeLocalPolynomialGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cint}), grid, dimensions, outputs, depth, order, sRule, limit_levels)
end

function tsgMakeWaveletGrid(grid, dimensions, outputs, depth, order, limit_levels)
    ccall((:tsgMakeWaveletGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Cint, Ptr{Cint}), grid, dimensions, outputs, depth, order, limit_levels)
end

function tsgMakeFourierGrid(grid, dimensions, outputs, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgMakeFourierGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, dimensions, outputs, depth, sType, anisotropic_weights, limit_levels)
end

function tsgUpdateGlobalGrid(grid, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgUpdateGlobalGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, depth, sType, anisotropic_weights, limit_levels)
end

function tsgUpdateSequenceGrid(grid, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgUpdateSequenceGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, depth, sType, anisotropic_weights, limit_levels)
end

function tsgUpdateFourierGrid(grid, depth, sType, anisotropic_weights, limit_levels)
    ccall((:tsgUpdateFourierGrid, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}), grid, depth, sType, anisotropic_weights, limit_levels)
end

function tsgGetAlpha(grid)
    ccall((:tsgGetAlpha, tasmanian), Cdouble, (Ptr{Cvoid},), grid)
end

function tsgGetBeta(grid)
    ccall((:tsgGetBeta, tasmanian), Cdouble, (Ptr{Cvoid},), grid)
end

function tsgGetOrder(grid)
    ccall((:tsgGetOrder, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumDimensions(grid)
    ccall((:tsgGetNumDimensions, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumOutputs(grid)
    ccall((:tsgGetNumOutputs, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetRule(grid)
    ccall((:tsgGetRule, tasmanian), Ptr{Cchar}, (Ptr{Cvoid},), grid)
end

function tsgCopyRuleChars(grid, buffer_size, name, num_actual)
    ccall((:tsgCopyRuleChars, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cchar}, Ptr{Cint}), grid, buffer_size, name, num_actual)
end

function tsgGetCustomRuleDescription(grid)
    ccall((:tsgGetCustomRuleDescription, tasmanian), Ptr{Cchar}, (Ptr{Cvoid},), grid)
end

function tsgGetNumLoaded(grid)
    ccall((:tsgGetNumLoaded, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumNeeded(grid)
    ccall((:tsgGetNumNeeded, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetNumPoints(grid)
    ccall((:tsgGetNumPoints, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetLoadedPointsStatic(grid, x)
    ccall((:tsgGetLoadedPointsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgGetLoadedPoints(grid)
    ccall((:tsgGetLoadedPoints, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetNeededPointsStatic(grid, x)
    ccall((:tsgGetNeededPointsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgGetNeededPoints(grid)
    ccall((:tsgGetNeededPoints, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetPointsStatic(grid, x)
    ccall((:tsgGetPointsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgGetPoints(grid)
    ccall((:tsgGetPoints, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetQuadratureWeightsStatic(grid, weights)
    ccall((:tsgGetQuadratureWeightsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, weights)
end

function tsgGetQuadratureWeights(grid)
    ccall((:tsgGetQuadratureWeights, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetInterpolationWeightsStatic(grid, x, weights)
    ccall((:tsgGetInterpolationWeightsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, x, weights)
end

function tsgGetInterpolationWeights(grid, x)
    ccall((:tsgGetInterpolationWeights, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid}, Ptr{Cdouble}), grid, x)
end

function tsgLoadNeededPoints(grid, vals)
    ccall((:tsgLoadNeededPoints, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, vals)
end

function tsgGetLoadedValues(grid)
    ccall((:tsgGetLoadedValues, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetLoadedValuesStatic(grid, values)
    ccall((:tsgGetLoadedValuesStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, values)
end

function tsgEvaluate(grid, x, y)
    ccall((:tsgEvaluate, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, x, y)
end

function tsgEvaluateFast(grid, x, y)
    ccall((:tsgEvaluateFast, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, x, y)
end

function tsgIntegrate(grid, q)
    ccall((:tsgIntegrate, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, q)
end

function tsgEvaluateBatch(grid, x, num_x, y)
    ccall((:tsgEvaluateBatch, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, num_x, y)
end

function tsgBatchGetInterpolationWeightsStatic(grid, x, num_x, weights)
    ccall((:tsgBatchGetInterpolationWeightsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, num_x, weights)
end

function tsgBatchGetInterpolationWeights(grid, x, num_x)
    ccall((:tsgBatchGetInterpolationWeights, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid}, Ptr{Cdouble}, Cint), grid, x, num_x)
end

function tsgIsGlobal(grid)
    ccall((:tsgIsGlobal, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsSequence(grid)
    ccall((:tsgIsSequence, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsLocalPolynomial(grid)
    ccall((:tsgIsLocalPolynomial, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsWavelet(grid)
    ccall((:tsgIsWavelet, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgIsFourier(grid)
    ccall((:tsgIsFourier, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgSetDomainTransform(grid, a, b)
    ccall((:tsgSetDomainTransform, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, a, b)
end

function tsgIsSetDomainTransfrom(grid)
    ccall((:tsgIsSetDomainTransfrom, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgClearDomainTransform(grid)
    ccall((:tsgClearDomainTransform, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgGetDomainTransform(grid, a, b)
    ccall((:tsgGetDomainTransform, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), grid, a, b)
end

function tsgSetConformalTransformASIN(grid, truncation)
    ccall((:tsgSetConformalTransformASIN, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), grid, truncation)
end

function tsgIsSetConformalTransformASIN(grid)
    ccall((:tsgIsSetConformalTransformASIN, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgClearConformalTransform(grid)
    ccall((:tsgClearConformalTransform, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgGetConformalTransformASIN(grid, truncation)
    ccall((:tsgGetConformalTransformASIN, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), grid, truncation)
end

function tsgClearLevelLimits(grid)
    ccall((:tsgClearLevelLimits, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgGetLevelLimits(grid, limits)
    ccall((:tsgGetLevelLimits, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), grid, limits)
end

function tsgSetAnisotropicRefinement(grid, sType, min_growth, output, level_limits)
    ccall((:tsgSetAnisotropicRefinement, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}), grid, sType, min_growth, output, level_limits)
end

function tsgEstimateAnisotropicCoefficients(grid, sType, output, num_coefficients)
    ccall((:tsgEstimateAnisotropicCoefficients, tasmanian), Ptr{Cint}, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}), grid, sType, output, num_coefficients)
end

function tsgEstimateAnisotropicCoefficientsStatic(grid, sType, output, coefficients)
    ccall((:tsgEstimateAnisotropicCoefficientsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}), grid, sType, output, coefficients)
end

function tsgSetGlobalSurplusRefinement(grid, tolerance, output, level_limits)
    ccall((:tsgSetGlobalSurplusRefinement, tasmanian), Cvoid, (Ptr{Cvoid}, Cdouble, Cint, Ptr{Cint}), grid, tolerance, output, level_limits)
end

function tsgSetLocalSurplusRefinement(grid, tolerance, sRefinementType, output, level_limits, scale_correction)
    ccall((:tsgSetLocalSurplusRefinement, tasmanian), Cvoid, (Ptr{Cvoid}, Cdouble, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}), grid, tolerance, sRefinementType, output, level_limits, scale_correction)
end

function tsgClearRefinement(grid)
    ccall((:tsgClearRefinement, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgMergeRefinement(grid)
    ccall((:tsgMergeRefinement, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgBeginConstruction(grid)
    ccall((:tsgBeginConstruction, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgIsUsingConstruction(grid)
    ccall((:tsgIsUsingConstruction, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

function tsgGetCandidateConstructionPointsVoidPntr(grid, sType, output, anisotropic_weights, limit_levels)
    ccall((:tsgGetCandidateConstructionPointsVoidPntr, tasmanian), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}), grid, sType, output, anisotropic_weights, limit_levels)
end

function tsgGetCandidateConstructionPointsSurplusVoidPntr(grid, tolerance, sRefType, output, limit_levels, scale_correction)
    ccall((:tsgGetCandidateConstructionPointsSurplusVoidPntr, tasmanian), Ptr{Cvoid}, (Ptr{Cvoid}, Cdouble, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}), grid, tolerance, sRefType, output, limit_levels, scale_correction)
end

function tsgGetCandidateConstructionPoints(grid, sType, output, anisotropic_weights, limit_levels, num_points, x)
    ccall((:tsgGetCandidateConstructionPoints, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cdouble}}), grid, sType, output, anisotropic_weights, limit_levels, num_points, x)
end

function tsgGetCandidateConstructionSurplusPoints(grid, tolerance, sRefType, output, limit_levels, scale_correction, num_points, x)
    ccall((:tsgGetCandidateConstructionSurplusPoints, tasmanian), Cvoid, (Ptr{Cvoid}, Cdouble, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Ptr{Cdouble}}), grid, tolerance, sRefType, output, limit_levels, scale_correction, num_points, x)
end

function tsgGetCandidateConstructionPointsPythonGetNP(grid, vecx)
    ccall((:tsgGetCandidateConstructionPointsPythonGetNP, tasmanian), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), grid, vecx)
end

function tsgGetCandidateConstructionPointsPythonStatic(vecx, x)
    ccall((:tsgGetCandidateConstructionPointsPythonStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), vecx, x)
end

function tsgGetCandidateConstructionPointsPythonDeleteVect(vecx)
    ccall((:tsgGetCandidateConstructionPointsPythonDeleteVect, tasmanian), Cvoid, (Ptr{Cvoid},), vecx)
end

function tsgLoadConstructedPoint(grid, x, numx, y)
    ccall((:tsgLoadConstructedPoint, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, numx, y)
end

function tsgFinishConstruction(grid)
    ccall((:tsgFinishConstruction, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgRemovePointsByHierarchicalCoefficient(grid, tolerance, output, scale_correction)
    ccall((:tsgRemovePointsByHierarchicalCoefficient, tasmanian), Cvoid, (Ptr{Cvoid}, Cdouble, Cint, Ptr{Cdouble}), grid, tolerance, output, scale_correction)
end

function tsgEvaluateHierarchicalFunctions(grid, x, num_x, y)
    ccall((:tsgEvaluateHierarchicalFunctions, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), grid, x, num_x, y)
end

function tsgEvaluateSparseHierarchicalFunctions(grid, x, num_x, pntr, indx, vals)
    ccall((:tsgEvaluateSparseHierarchicalFunctions, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), grid, x, num_x, pntr, indx, vals)
end

function tsgEvaluateSparseHierarchicalFunctionsGetNZ(grid, x, num_x)
    ccall((:tsgEvaluateSparseHierarchicalFunctionsGetNZ, tasmanian), Cint, (Ptr{Cvoid}, Ptr{Cdouble}, Cint), grid, x, num_x)
end

function tsgEvaluateSparseHierarchicalFunctionsStatic(grid, x, num_x, pntr, indx, vals)
    ccall((:tsgEvaluateSparseHierarchicalFunctionsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), grid, x, num_x, pntr, indx, vals)
end

function tsgGetHierarchicalSupportStatic(grid, support)
    ccall((:tsgGetHierarchicalSupportStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, support)
end

function tsgGetHierarchicalCoefficients(grid)
    ccall((:tsgGetHierarchicalCoefficients, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgGetHierarchicalCoefficientsStatic(grid, coeff)
    ccall((:tsgGetHierarchicalCoefficientsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, coeff)
end

function tsgSetHierarchicalCoefficients(grid, c)
    ccall((:tsgSetHierarchicalCoefficients, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, c)
end

function tsgIntegrateHierarchicalFunctions(grid)
    ccall((:tsgIntegrateHierarchicalFunctions, tasmanian), Ptr{Cdouble}, (Ptr{Cvoid},), grid)
end

function tsgIntegrateHierarchicalFunctionsStatic(grid, integrals)
    ccall((:tsgIntegrateHierarchicalFunctionsStatic, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}), grid, integrals)
end

function tsgPythonGetGlobalPolynomialSpace(grid, interpolation, num_indexes)
    ccall((:tsgPythonGetGlobalPolynomialSpace, tasmanian), Ptr{Cint}, (Ptr{Cvoid}, Cint, Ptr{Cint}), grid, interpolation, num_indexes)
end

function tsgGetGlobalPolynomialSpace(grid, interpolation, num_indexes, indexes)
    ccall((:tsgGetGlobalPolynomialSpace, tasmanian), Cvoid, (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}), grid, interpolation, num_indexes, indexes)
end

function tsgPrintStats(grid)
    ccall((:tsgPrintStats, tasmanian), Cvoid, (Ptr{Cvoid},), grid)
end

function tsgEnableAcceleration(grid, accel)
    ccall((:tsgEnableAcceleration, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}), grid, accel)
end

function tsgEnableAccelerationGPU(grid, accel, gpu)
    ccall((:tsgEnableAccelerationGPU, tasmanian), Cvoid, (Ptr{Cvoid}, Ptr{Cchar}, Cint), grid, accel, gpu)
end

function tsgGetAccelerationType(grid)
    ccall((:tsgGetAccelerationType, tasmanian), Ptr{Cchar}, (Ptr{Cvoid},), grid)
end

function tsgSetGPUID(grid, gpuID)
    ccall((:tsgSetGPUID, tasmanian), Cvoid, (Ptr{Cvoid}, Cint), grid, gpuID)
end

function tsgGetGPUID(grid)
    ccall((:tsgGetGPUID, tasmanian), Cint, (Ptr{Cvoid},), grid)
end

# no prototype is found for this function at TasmanianSparseGrid.h:143:5, please use with caution
function tsgGetNumGPUs()
    ccall((:tsgGetNumGPUs, tasmanian), Cint, ())
end

function tsgGetGPUMemory(gpu)
    ccall((:tsgGetGPUMemory, tasmanian), Cint, (Cint,), gpu)
end

function tsgIsAccelerationAvailable(accel)
    ccall((:tsgIsAccelerationAvailable, tasmanian), Cint, (Ptr{Cchar},), accel)
end

function tsgGetGPUName(gpu, num_buffer, buffer, num_actual)
    ccall((:tsgGetGPUName, tasmanian), Cvoid, (Cint, Cint, Ptr{Cchar}, Ptr{Cint}), gpu, num_buffer, buffer, num_actual)
end

function tsgDeleteInts(p)
    ccall((:tsgDeleteInts, tasmanian), Cvoid, (Ptr{Cint},), p)
end

end # module

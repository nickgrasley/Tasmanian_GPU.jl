using Tasmanian
using Test

function getSparseGridTests()
    # The list llTest contains list-of-lists (could make into tuples)
    # Each sub-list (tuple) has two entries, one is the command that should raise the exception,
    # the second is the "sError" variable in the exception class.
    # notError tests here are needed to ensure that the multi-statement commands fail for the correct function.
    return [["grid = makeGlobalGrid(-1, 1,  4, \"level\", \"clenshaw-curtis\")", "iDimension"],
            ["grid = makeGlobalGrid(2, -1,  4, \"level\", \"clenshaw-curtis\")", "iOutputs"],
            ["grid = makeGlobalGrid(2,  1, -4, \"level\", \"clenshaw-curtis\")", "iDepth"],
            ["grid = makeGlobalGrid(2,  1,  4, \"wrong\", \"clenshaw-curtis\")", "sType"],
            ["grid = makeGlobalGrid(2,  1,  4, \"level\", \"clenshaw-wrong\")", "sRule"],
            ["grid = makeGlobalGrid(2,  1,  4, \"level\", \"clenshaw-curtis\", anisotropic_weights=[1,2,3])", "liAnisotropicWeights"],
            ["grid = makeGlobalGrid(2,  1,  4, \"level\", \"clenshaw-curtis\", anisotropic_weights=[1,2], levelLimits = [1, 2, 3])", "liLevelLimits"],
            ["grid = makeGlobalGrid(2,  1,  4, \"level\", \"clenshaw-curtis\", anisotropic_weights=[1,2], levelLimits = [1, 2])", "notError"],
            ["grid = makeSequenceGrid(-1, 1,  4, \"level\", \"leja\")", "iDimension"],
            ["grid = makeSequenceGrid(2, -1,  4, \"level\", \"leja\")", "iOutputs"],
            ["grid = makeSequenceGrid(2,  1, -4, \"level\", \"leja\")", "iDepth"],
            ["grid = makeSequenceGrid(2,  1,  4, \"wrong\", \"leja\")", "sType"],
            ["grid = makeSequenceGrid(2,  1,  4, \"level\", \"weja\")", "sRule"],
            ["grid = makeSequenceGrid(2,  1,  4, \"level\", \"leja\", anisotropic_weights=[1, 2, 3])", "liAnisotropicWeights"],
            ["grid = makeSequenceGrid(2,  1,  4, \"level\", \"leja\", anisotropic_weights=[1, 2], levelLimits = [1, 2, 3])", "liLevelLimits"],
            ["grid = makeSequenceGrid(2,  1,  4, \"level\", \"leja\", levelLimits = [1, 2])", "notError"],
            ["grid = makeLocalPolynomialGrid(-1, 1,  4,  2, \"localp\")", "iDimension"],
            ["grid = makeLocalPolynomialGrid(2, -1,  4,  2, \"localp\")", "iOutputs"],
            ["grid = makeLocalPolynomialGrid(2,  1, -4,  2, \"localp\")", "iDepth"],
            ["grid = makeLocalPolynomialGrid(2,  1,  4, -2, \"localp\")", "iOrder"],
            ["grid = makeLocalPolynomialGrid(2,  1,  4,  2, \"lowrong\")", "sRule"],
            ["grid = makeLocalPolynomialGrid(2,  1,  4,  2, \"localp\", levelLimits=[1, 2, 3])", "liLevelLimits"],
            ["grid = makeLocalPolynomialGrid(2,  1,  4,  2, \"localp\", levelLimits=[1, 2])", "notError"],
            ["grid = makeWaveletGrid(-1, 1,  4,  1)", "iDimension"],
            ["grid = makeWaveletGrid(2, -1,  4,  1)", "iOutputs"],
            ["grid = makeWaveletGrid(2,  1, -4,  3)", "iDepth"],
            ["grid = makeWaveletGrid(2,  1,  4,  2)", "iOrder"],
            ["grid = makeWaveletGrid(2,  1,  4,  1, levelLimits=[1, 2, 3])", "liLevelLimits"],
            ["grid = makeWaveletGrid(2,  1,  4,  1, levelLimits=[2, 1])", "notError"],
            ["grid = makeFourierGrid(-1, 1,  4, \"level\")", "iDimension"],
            ["grid = makeFourierGrid(2, -1,  4, \"level\")", "iOutputs"],
            ["grid = makeFourierGrid(2,  1, -4, \"level\")", "iDepth"],
            ["grid = makeFourierGrid(2,  1,  4, \"wrong\")", "sType"],
            ["grid = makeFourierGrid(2,  1,  4, \"level\", anisotropic_weights=[1, 2, 3])", "liAnisotropicWeights"],
            ["grid = makeFourierGrid(2,  1,  4, \"level\", anisotropic_weights=[1, 2], levelLimits = [1, 2, 3])", "liLevelLimits"],
            ["grid = makeFourierGrid(2,  1,  4, \"level\", levelLimits = [1, 2])", "notError"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\")", "notError"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"chebyshev\")", "notError"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededValues!(grid, zeros(2, 6))", "notError"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); updateGlobalGrid!(grid, 1,\"iptotal\")", "updateGlobalGrid"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"chebyshev\"); updateGlobalGrid!(grid, -1,\"iptotal\")", "iDepth"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"chebyshev\"); updateGlobalGrid!(grid, 4,\"wrong\")", "sType"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"chebyshev\"); updateGlobalGrid!(grid, 4,\"iptotal\", anisotropic_weights=[1,2,3])", "liAnisotropicWeights"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"chebyshev\"); updateGlobalGrid!(grid, 4,\"iptotal\",levelLimits = [1,2,3])", "liLevelLimits"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"chebyshev\"); updateGlobalGrid!(grid, 4,\"iptotal\",levelLimits = [1,2])", "notError"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"rleja\"); updateSequenceGrid!(grid, 4,\"iptotal\")", "updateSequenceGrid"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); updateSequenceGrid!(grid, -1,\"iptotal\")", "iDepth"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); updateSequenceGrid!(grid, 4,\"wrong\")", "sType"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); updateSequenceGrid!(grid, 4,\"iptotal\",anisotropic_weights=[1,2,3])", "liAnisotropicWeights"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); updateSequenceGrid!(grid, 4,\"iptotal\",levelLimits = [1,2,3])", "liLevelLimits"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); updateSequenceGrid!(grid, 4,\"iptotal\",levelLimits = [2,3])", "notError"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); updateFourierGrid!(grid, 3,\"level\")", "updateFourierGrid"],
            ["grid = makeFourierGrid(2, 1, 2, \"level\"); updateFourierGrid!(grid, -3,\"level\")", "iDepth"],
            ["grid = makeFourierGrid(2, 1, 2, \"level\"); updateFourierGrid!(grid, 3,\"wrong\")", "sType"],
            ["grid = makeFourierGrid(2, 1, 2, \"level\"); updateFourierGrid!(grid, 3,\"iptotal\",anisotropic_weights=[1])", "liAnisotropicWeights"],
            ["grid = makeFourierGrid(2, 1, 2, \"level\"); updateFourierGrid!(grid, 3,\"iptotal\",levelLimits = [4])", "liLevelLimits"],
            ["grid = makeFourierGrid(2, 1, 2, \"level\"); updateFourierGrid!(grid, 3,\"iptotal\", anisotropic_weights=[4, 3], levelLimits=[4, 3])", "notError"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); aW = getInterpolationWeights(grid, [1,2,3])", "lfX"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); aW = getInterpolationWeightsBatch(grid, [1,2,3])", "llfX"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); aW = getInterpolationWeightsBatch(grid, [1 1; 2 2; 3 3])", "llfX"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(1))", "llfVals"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(3,6))", "llfVals"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,5))", "llfVals"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,6)); loadNeededPoints!(grid, ones(2,5))", "llfVals"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); evaluate(grid, zeros(2,1))", "evaluate"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); evaluateThreadSafe(grid, zeros(1,2))", "evaluateThreadSafe"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,6)); evaluate(grid, zeros(3,1))", "lfX"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,6)); evaluate(grid, zeros(3))", "lfX"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,6)); evaluateThreadSafe(grid, zeros(1,3))", "lfX"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,6)); evaluateThreadSafe(grid, zeros(3))", "lfX"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); evaluateBatch(grid, zeros(2,1))", "evaluateBatch"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,6)); evaluateBatch(grid, zeros(3,1))", "llfX"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); loadNeededPoints!(grid, zeros(2,6)); evaluateBatch(grid, zeros(1, 2))", "llfX"],
            ["grid = makeSequenceGrid(2, 2, 2, \"level\", \"rleja\"); integrate(grid)", "integrate"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"gauss-legendre\"); setDomainTransform!(grid, zeros(2))", "llfTransform"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"gauss-legendre\"); setDomainTransform!(grid, zeros(1,2))", "llfTransform"],
            ["grid = makeGlobalGrid(3, 1, 2, \"level\", \"gauss-legendre\"); setDomainTransform!(grid, zeros(2,2))", "llfTransform"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"gauss-legendre\")", "notError"],
            ["grid = makeGlobalGrid(3, 1, 2, \"level\", \"gauss-legendre\"); setDomainTransform!(grid, zeros(3,2))", "notError"],
            ["grid = makeGlobalGrid(3, 1, 2, \"level\", \"clenshaw-curtis\")", "notError"],
            ["grid = makeGlobalGrid(3, 1, 2, \"level\", \"clenshaw-curtis\"); setConformalTransformASIN!(grid, [0,2,4])", "notError"],
            ["grid = makeGlobalGrid(3, 1, 2, \"level\", \"clenshaw-curtis\"); setConformalTransformASIN!(grid, [0,2])", "liTruncation"],
            ["grid = makeGlobalGrid(3, 1, 2, \"level\", \"clenshaw-curtis\"); setConformalTransformASIN!(grid, [0 2 3; 1 2 3])", "liTruncation"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"rleja\"); setAnisotropicRefinement!(grid, \"iptotal\", 10, 1)", "setAnisotropicRefinement"],
            ["grid = makeGlobalGrid(2, 0, 2, \"level\", \"clenshaw-curtis\"); setAnisotropicRefinement!(grid, \"iptotal\", 10, 1)", "setAnisotropicRefinement"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"fejer2\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", -2, 1)", "iMinGrowth"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, -1)", "iOutput"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, 0)", "notError"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, 0, [2, 3])", "notError"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, 0, [2, 3, 3])", "liLevelLimits"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, -1)", "notError"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, -2)", "iOutput"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, 5)", "iOutput"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"wrong\", 10, 0)", "sType"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, -1, [3, 4])", "notError"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); setAnisotropicRefinement!(grid, \"iptotal\", 10, -1, [3, 4, 5])", "liLevelLimits"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"rleja\"); estimateAnisotropicCoefficients(grid, \"iptotal\", 1)", "estimateAnisotropicCoefficients"],
            ["grid = makeGlobalGrid(2, 0, 2, \"level\", \"clenshaw-curtis\"); estimateAnisotropicCoefficients(grid, \"iptotal\", 1)", "estimateAnisotropicCoefficients"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); loadExpN2(grid); estimateAnisotropicCoefficients(grid, \"iptotal\", -1)", "iOutput"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); estimateAnisotropicCoefficients(grid, \"iptotal\", -1)", "notError"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); estimateAnisotropicCoefficients(grid, \"ipcurved\", -1)", "notError"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); estimateAnisotropicCoefficients(grid, \"iptotal\", -2)", "iOutput"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); estimateAnisotropicCoefficients(grid, \"iptotal\", 5)", "iOutput"],
            ["grid = makeSequenceGrid(2, 1, 3, \"iptotal\", \"leja\"); loadExpN2(grid); estimateAnisotropicCoefficients(grid, \"wrong\", 0)", "sType"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\")", "setSurplusRefinement"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\")", "setSurplusRefinement"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); loadExpN2(grid); setSurplusRefinement!(grid, -1.E-4, output=0, refinement_type=\"classic\")", "fTolerance"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\")", "sCriteria"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0)", "notError"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"\", level_Limits=[2, 3, 4])", "liLevelLimits"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"\", level_Limits=[2, 3])", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0)", "sCriteria"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"class\")", "sCriteria"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\")", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\", level_Limits=[2, 3, 4])", "liLevelLimits"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\", level_Limits=[], scale_correction=ones(1, 3))", "llfScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\", level_Limits=[], scale_correction=ones(getNumPoints(grid) - 1, 1))", "llfScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\", level_Limits=[], scale_correction=ones(getNumPoints(grid), 2))", "llfScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 2, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=-1, refinement_type=\"classic\", level_Limits=[], scale_correction=ones(getNumPoints(grid), 2))", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 2, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=-1, refinement_type=\"classic\", level_Limits=[], scale_correction=ones(getNumPoints(grid), 3))", "llfScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); setSurplusRefinement!(grid, 1.E-4, output=0, refinement_type=\"classic\", level_Limits=[2, 3])", "notError"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); removePointsByHierarchicalCoefficient!(grid, 1.E-4, 0)", "removePointsByHierarchicalCoefficient"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); removePointsByHierarchicalCoefficient!(grid, -1.E-4, 0)", "fTolerance"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); removePointsByHierarchicalCoefficient!(grid, 1.E-4, -2)", "iOutput"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); removePointsByHierarchicalCoefficient!(grid, 1.E-4, 3)", "iOutput"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); removePointsByHierarchicalCoefficient!(grid, 1.E-4, 0)", "removePointsByHierarchicalCoefficient"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); loadExpN2(grid); removePointsByHierarchicalCoefficient!(grid, 1.E-4, 0)", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, -1, ones(3))", "aScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, -1, ones(10))", "aScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, -1, ones(2,11))", "aScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, -1, ones(2,13))", "aScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, -1, ones(3,13))", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, 0,  ones(3,13))", "aScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, 0,  ones(11))", "aScaleCorrection"],
            ["grid = makeLocalPolynomialGrid(2, 3, 2, 1, \"localp\"); loadNeededPoints!(grid, ones(getNumOutputs(grid), getNumNeeded(grid))); removePointsByHierarchicalCoefficient!(grid, 1.E-4, 0,  ones(1,13))", "notError"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); evaluateHierarchicalFunctions(grid, [1.0 1.0; 0.5 0.3])", "notError"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); evaluateHierarchicalFunctions(grid, [1.0; 0.5])", "llfX"],
            ["grid = makeGlobalGrid(2, 1, 2, \"level\", \"clenshaw-curtis\"); evaluateHierarchicalFunctions(grid, [1.0, 1.0])", "llfX"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); evaluateSparseHierarchicalFunctions(grid, [1.0 1.0; 0.5 0.3])", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); evaluateSparseHierarchicalFunctions(grid, [1.0 0.5])", "llfX"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); evaluateSparseHierarchicalFunctions(grid, [1.0, 1.0])", "llfX"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); setHierarchicalCoefficients!(grid, [1.0, 1.0])", "llfCoefficients"],
            ["grid = makeLocalPolynomialGrid(2, 1, 2, 1, \"localp\"); setHierarchicalCoefficients!(grid, [1.0 1.0])", "llfCoefficients"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); setHierarchicalCoefficients!(grid, [1.0 1.0; 1.0 1.0; 1.0 1.0; 1.0 1.0; 1.0 1.0])", "llfCoefficients"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); setHierarchicalCoefficients!(grid, ones(1, 5))", "notError"],
            ["grid = makeFourierGrid(2, 1, 1, \"level\"); setHierarchicalCoefficients!(grid, ones(1, 5))", "llfCoefficients"],
#=
            ["grid = makeGlobalGrid(2, 1, 3, \"level\", \"rleja\"); getCandidateConstructionPoints(grid, \"level\", -1)", "getCandidateConstructionPoints"],
            ["grid = makeGlobalGrid(2, 1, 3, \"level\", \"rleja\"); beginConstruction(grid); getCandidateConstructionPoints(grid, \"lev\", -1)", "sType"],
            ["grid = makeGlobalGrid(2, 1, 3, \"level\", \"rleja\"); beginConstruction(grid); getCandidateConstructionPoints(grid, \"level\", [0])", "liAnisotropicWeightsOrOutput"],
            ["grid = makeGlobalGrid(2, 1, 3, \"level\", \"rleja\"); beginConstruction(grid); getCandidateConstructionPoints(grid, \"level\", \"string\")", "liAnisotropicWeightsOrOutput"],
            ["grid = makeGlobalGrid(2, 1, 3, \"level\", \"rleja\"); beginConstruction(grid); getCandidateConstructionPoints(grid, \"level\", -1, [2])", "liLevelLimits"],
            ["grid = makeGlobalGrid(2, 1, 3, \"level\", \"rleja\"); beginConstruction(grid); getCandidateConstructionPoints(grid, \"level\", -1, [2, 1])", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); getCandidateConstructionPointsSurplus(grid, 1.E-5, \"classic\")", "getCandidateConstructionPointsSurplus"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); getCandidateConstructionPointsSurplus(grid, 1.E-5, \"classic\", 0, [2])", "liLevelLimits"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); getCandidateConstructionPointsSurplus(grid, 1.E-5, \"classic\", -1, [2, 3])", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [0.0, 0.0], [1.0])", "notError"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); loadConstructedPoint(grid, [0.0, 0.0], [1.0])", "loadConstructedPoint"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [0.0], [1.0])", "lfX"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [0.0, 0.0], [1.0, 2.0])", "lfY"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [[0.0, 0.0], [1.0, 0.0]], [1.0, 2.0])", "lfY"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [[0.0,], [1.0,]], [[1.0,], [2.0,]])", "lfX"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [[0.0, 0.0], [1.0, 0.0]], [[1.0, 2.0], [1.0, 2.0]])", "lfY"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [[0.0, 0.0], [1.0, 0.0]], [[1.0,], [2.0,], [3.0,]])", "lfY"],
            ["grid = makeLocalPolynomialGrid(2, 1, 1, 1, \"localp\"); beginConstruction(grid); loadConstructedPoint(grid, [[[0.0, 0.0], [1.0, 0.0]],], [[1.0,], [2.0,]])", "lfX"],
            =#
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); enableAcceleration(grid, \"gpu-wrong\")", "sAccelerationType"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); enableAcceleration(grid, \"gpu-default\")", "notError"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); enableAcceleration(grid, \"gpu-default\", GPUID=isAccelerationAvailable(grid, \"gpu-cuda\") ? 0 : nothing)", "notError"],
            ["grid = makeSequenceGrid(2, 1, 2, \"level\", \"leja\"); enableAcceleration(grid, \"gpu-default\", GPUID=-11)", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); isAccelerationAvailable(grid1, \"cpu-wrong\")", "sAccelerationType"],
            ["grid1 = TasmanianSG(0,0,0); isAccelerationAvailable(grid1, \"cpu-blas\")", "notError"],
            ["grid1 = TasmanianSG(0,0,0); getGPUMemory(grid1, -1)", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); getGPUMemory(grid1, 1000000)", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); getGPUMemory(grid1, getNumGPUs())", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); getGPUName(grid1, -1)", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); getGPUName(grid1, 1000000)", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); getGPUName(grid1, getNumGPUs())", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); setGPUID!(grid1, -1)", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); setGPUID!(grid1, 1000000)", "iGPUID"],
            ["grid1 = TasmanianSG(0,0,0); setGPUID!(grid1, getNumGPUs())", "iGPUID"],
            #=
            ["grid = makeLocalPolynomialGrid(1, 1, 1, 1, \"localp\"); Tasmanian.loadNeededPoints!(grid, lambda x, tid : x, grid, 1)", "notError"],
            ["grid = makeLocalPolynomialGrid(1, 1, 1, 1, \"localp\"); Tasmanian.loadNeededPoints!(grid, lambda x, tid : np.ones((2,)) * x, grid, 1)", "loadNeededValues"],
            =#
            ]
end

function makeGlobalGrid(dims, out, depth, sType, sRule; anisotropic_weights=Vector{Int32}(undef, 0), alpha=0.0, beta=0.0, custom_filename="", levelLimits=Vector{Int32}(undef, 0))
    tsg = TasmanianSG(dims, out, depth)
    makeGlobalGrid!(tsg, sType=sType, sRule=sRule, anisotropic_weights=anisotropic_weights, alpha=alpha, beta=beta, custom_filename=custom_filename, levelLimits=levelLimits)
    return tsg
end

function makeSequenceGrid(dims, out, depth, sType, sRule; anisotropic_weights=Vector{Int32}(undef, 0), levelLimits=Vector{Int32}(undef, 0))
    tsg = TasmanianSG(dims, out, depth)
    makeSequenceGrid!(tsg, sType=sType, sRule=sRule, anisotropic_weights=anisotropic_weights, levelLimits=levelLimits)
    return tsg
end

function makeLocalPolynomialGrid(dims, out, depth, order, sRule; levelLimits=Vector{Int32}(undef, 0))
    tsg = TasmanianSG(dims, out, depth)
    makeLocalPolynomialGrid!(tsg, order=order, sRule=sRule, levelLimits=levelLimits)
    return tsg
end

function makeWaveletGrid(dims, out, depth, order; levelLimits=Vector{Int32}(undef, 0))
    tsg = TasmanianSG(dims, out, depth)
    makeWaveletGrid!(tsg, order=order, levelLimits=levelLimits)
    return tsg
end

function makeFourierGrid(dims, out, depth, sType; anisotropic_weights=Vector{Int32}(undef, 0), levelLimits=Vector{Int32}(undef, 0))
    tsg = TasmanianSG(dims, out, depth)
    makeFourierGrid!(tsg; sType=sType, anisotropic_weights=anisotropic_weights, levelLimits=levelLimits)
    return tsg
end

function testListedExceptions(Tests)
    for test in Tests
        @show test
        if test[2] == "notError"
            eval(Meta.parse(test[1]))
        else
            @test_throws Tasmanian.TasmanianInputError eval(Meta.parse(test[1]))
        end
    end
end

testListedExceptions(getSparseGridTests())

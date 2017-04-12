__precompile__()
module JFVMM

import Base: +, -, *, /, .*, ./
export MeshStructure, BoundaryCondition, CellValue, FaceValue, CellVector,
       arithmeticMean, geometricMean, harmonicMean, upwindMean, linearMean,
       tvdMean, createBC, boundaryConditionTerm, cellBoundary, solvePDE,
       divergenceTerm, gradientTerm, convectionUpwindTerm, createCellVector,
       convectionTerm, convectionTvdTerm, diffusionTerm, createCellVariable,
       createFaceVariable, copyCell, fluxLimiter, createMesh1D,
       createMesh2D, createMesh3D, createMeshRadial2D, createMeshCylindrical2D,
       createMeshCylindrical3D, createMeshCylindrical1D, solveLinearPDE,
       linearSourceTerm, constantSourceTerm, transientTerm,
       solveMUMPSLinearPDE, faceEval, cellEval, permfieldlogrndg, permfieldlogrnde,
       # plot, imshow, xlabel, ylabel, figure, legend, pcolor, contour, colorbar,
       JFVM_test, solveExplicitPDE, reshapeCell,
       cellVolume, reshapeInternalCell, internalCells, domainInt, convectionTvdRHS

# Ignore visualization for now
# visualizeCellVectors
# visualizeCells

include("jfvmm_1types.jl")
include("jfvmm_2meshstructure.jl")
include("jfvmm_3boundarycondition.jl")
include("jfvmm_4domain_variables.jl")
include("jfvmm_5diffusion_terms.jl")
include("transientTerms.jl")
include("domainOperators.jl")
include("convectionTerms.jl")
include("averagingTerms.jl")
include("calculusTerms.jl")
include("sourceTerms.jl")
include("solveVisualizePDE.jl")
include("JFVMtools.jl")
include("jfvm_test.jl")

end # module

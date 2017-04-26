__precompile__()
module JFVMM

import Base: +, -, *, /
# export MeshStructure, BoundaryCondition, CellValue, FaceValue, CellVector,
#        arithmeticMean, geometricMean, harmonicMean, upwindMean, linearMean,
#        tvdMean, createBC, boundaryConditionTerm, cellBoundary, solvePDE,
#        divergenceTerm, gradientTerm, convectionUpwindTerm, createCellVector,
#        convectionTerm, convectionTvdTerm, diffusionTerm, createCellVariable,
#        createFaceVariable, copyCell, fluxLimiter, createMesh1D,
#        createMesh2D, createMesh3D, createMeshRadial2D, createMeshCylindrical2D,
#        createMeshCylindrical3D, createMeshCylindrical1D, solveLinearPDE,
#        linearSourceTerm, constantSourceTerm, transientTerm,
#        solveMUMPSLinearPDE, faceEval, cellEval, permfieldlogrndg, permfieldlogrnde,
#        # plot, imshow, xlabel, ylabel, figure, legend, pcolor, contour, colorbar,
#        JFVM_test, solveExplicitPDE, reshapeCell,
#        cellVolume, reshapeInternalCell, internalCells, domainInt, convectionTvdRHS

# Ignore visualization for now
# visualizeCellVectors
# visualizeCells

include("jfvmm_01_types.jl")
include("jfvmm_02_domain_variables.jl")
include("jfvmm_03_mesh_structure.jl")
include("jfvmm_04_boundary_condition.jl")
# include("jfvmm_05_diffusion_terms.jl")
# include("jfvmm_06_transient_terms.jl")
# include("jfvmm_07_domain_operators.jl")
# include("jfvmm_08_convection_terms.jl")
# include("jfvmm_09_averaging_terms.jl")
# include("jfvmm_10_calculus_terms.jl")
# include("jfvmm_11_source_terms.jl")
# include("jfvmm_12_solve_pde.jl")
# include("jfvmm_13_tools.jl")
# include("jfvm_test.jl")

end # module

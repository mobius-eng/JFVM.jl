using Base.Test

import JFVMM
# %% Import


# %% points array mesh
# m2d = FV.createMesh([FV.Point(i*0.1,j*0.5) for i in 0:4, j in 0:3])

# %% assign variables
FV = JFVMM
m2d = FV.createRectMesh(3.0, 2.0, 3, 4)
bc2d = FV.createBoundaryCondition(m2d)
φ = FV.createCellField(m2d, 1.0)
φf = FV.linearMean(φ)
dφ = FV.createCellField(m2d, FV.Point2D)
# %% Test simple meshing
# 2D
# test dimensions
@test size(m2d) == (3, 4)
# test (generated) vertices
@test size(FV.nodes(m2d)) == (4, 5)
@test FV.centroid(FV.nodes(m2d)[3,2]) ≈ FV.Point2D(2.0, 0.5)
@test FV.centroid(FV.nodes(m2d)[1,1]) ≈ FV.Point2D(0.0, 0.0)
# test cells
@test size(FV.cells(m2d)) == (3,4)
@test all(x -> x.volume ≈ 0.5, m2d |> FV.cells)
@test FV.centroid(FV.cells(m2d)[3,2]) ≈ FV.Point2D(2.5, 0.75)
#test faces
@test size(m2d |> FV.ifaces) == (4,4)
@test size(m2d |> FV.jfaces) == (3,5)
@test all(x -> x.area ≈ 0.5, m2d |> FV.ifaces)
@test all(x -> x.area ≈ 1, m2d |> FV.jfaces)
@test FV.centroid(FV.ifaces(m2d)[3,2]) ≈ FV.Point2D(2.0, 0.75)
@test FV.centroid(FV.jfaces(m2d)[3,2]) ≈ FV.Point2D(2.5, 0.5)
@test m2d.ifaces[3,2].Sf / m2d.ifaces[3,2].area ≈ FV.Point2D(-1.0, 0.0)
@test m2d.jfaces[3,2].Sf ≈ FV.Point2D(0.0, -1.0) * m2d.jfaces[3,2].area
# special attention to boundary faces on top-right side
@test FV.centroid(m2d.ifaces[4,2]) ≈ FV.Point2D(3.0, 0.75)
@test FV.centroid(m2d.jfaces[3, 5]) ≈ FV.Point2D(2.5, 2.0)
@test FV.facevector(FV.ifaces(m2d)[4,2]) / FV.facearea(FV.ifaces(m2d)[4,2]) ≈ FV.Point2D(1,0)
@test FV.facevector(FV.jfaces(m2d)[3,5]) / FV.facearea(FV.jfaces(m2d)[3,5]) ≈ FV.Point2D(0,1)
# %% Test boundary conditions
FV.cellGradient!(dφ, φf)


bc2d.left.value[:] = [1.0, 0.0, 0.0, 0.0]
bc2d.left.value
m2d.faces.ifaces[1, 4].area
FV.boundaryConditionTerm(bc2d, FV.FaceValue(m2d, [0.0], [0.0], [0.0]))
m2d.cells[2,4]
# FV.dimensions(m2d)

# %%
# map(x -> x.normal, m3d.faces.ifaces)
FV.createCellField(m2d, 3.0)

size(m2d)

FV.createCellField(m2d, [(i,j) for i in 1:3,
    j in 1:4])

FV.createCellField(m2d, p -> p)


FV.createCellValue(m3d, 3.0)
FV.dimensions(m2d)


FV.Fields.boundaryFaceValuesOfCellValue

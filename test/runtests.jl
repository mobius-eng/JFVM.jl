using Base.Test

# %% Import
import JFVMM
FV = JFVMM

# %% assign variables
m2d = FV.createMesh2D(3.0, 2.0, 3, 4)
bc2d = FV.createBC(m2d)
m3d = FV.createMesh3D(3.0, 2.0, 4.0, 3, 4, 2)
bc3d = FV.createBC(m3d)
# %% Test simple meshing
# 2D
# test dimensions
@test m2d.dims[1] == 3
@test m2d.dims[2] == 4
@test length(m2d.dims) == 2
# test (generated) vertices
@test size(m2d.vertices) == (4, 5)
@test m2d.vertices[1,1] ≈ FV.Point(0.0, 0.0)
@test m2d.vertices[3,2] ≈ FV.Point(2.0, 0.5)
# test cells
@test size(m2d.cells.volumes)[1:2] == (3,4)
@test size(m2d.cells.centroids)[1:2] == (3,4)
@test all(x -> x ≈ 0.5, m2d.cells.volumes)
@test m2d.cells.centroids[3,2] ≈ FV.Point(2.5, 0.75)
#test faces
@test size(m2d.faces.ifaces)[1:2] == (4,4)
@test size(m2d.faces.jfaces)[1:2] == (3,5)
@test all(x -> x.area ≈ 0.5, m2d.faces.ifaces)
@test all(x -> x.area ≈ 1, m2d.faces.jfaces)
@test m2d.faces.ifaces[3,2].loc ≈ FV.Point(2.0, 0.75)
@test m2d.faces.jfaces[3,2].loc ≈ FV.Point(2.5, 0.5)
@test m2d.faces.ifaces[3,2].normal ≈ FV.Point(1.0, 0.0)
@test m2d.faces.jfaces[3,2].normal ≈ FV.Point(0.0, 1.0)
# special attention to boundary faces on top-right side
@test m2d.faces.jfaces[3, 5].loc ≈ FV.Point(2.5, 2.0)
@test m2d.faces.ifaces[4,2].loc ≈ FV.Point(3.0, 0.75)
@test m2d.faces.jfaces[3, 5].normal ≈ FV.Point(0.0, 1.0)
@test m2d.faces.ifaces[4,2].normal ≈ FV.Point(1.0, 0.0)
# 3D
# test dimensions
@test m3d.dims[1] == 3
@test m3d.dims[2] == 4
@test m3d.dims[3] == 2
@test length(m3d.dims) == 3
# test (generated) vertices
@test size(m3d.vertices) == (4, 5, 3)
@test m3d.vertices[1,1,1] ≈ FV.Point(0.0, 0.0, 0.0)
@test m3d.vertices[3,2,3] ≈ FV.Point(2.0, 0.5, 4.0)
# test cells
@test size(m3d.cells.volumes) == (3,4,2)
@test size(m3d.cells.centroids) == (3,4,2)
@test all(x -> x ≈ 1.0, m3d.cells.volumes)
@test m3d.cells.centroids[3,2,2] ≈ FV.Point(2.5, 0.75, 3.0)
#test faces
@test size(m3d.faces.ifaces) == (4,4,2)
@test size(m3d.faces.jfaces) == (3,5,2)
@test size(m3d.faces.kfaces) == (3,4,3)
@test all(x -> x.area≈ 1.0, m3d.faces.ifaces)
@test all(x -> x.area ≈ 2.0, m3d.faces.jfaces)
@test all(x -> x.area ≈ 0.5, m3d.faces.kfaces)
@test m3d.faces.ifaces[3,2,2].loc ≈ FV.Point(2.0, 0.75, 3.0)
@test m3d.faces.jfaces[3,2,2].loc ≈ FV.Point(2.5, 0.5,3.0)
@test m3d.faces.kfaces[3,2,2].loc ≈ FV.Point(2.5, 0.75,2.0)
@test m3d.faces.ifaces[3,2,2].normal ≈ FV.Point(1.0, 0.0, 0.0)
@test m3d.faces.jfaces[3,2,2].normal ≈ FV.Point(0.0, 1.0, 0.0)
@test m3d.faces.kfaces[3,2,2].normal ≈ FV.Point(0.0, 0.0, 1.0)
# special attention to boundary faces on top-right side
@test m3d.faces.ifaces[4,2,2].loc ≈ FV.Point(3.0, 0.75,3.0)
@test m3d.faces.jfaces[3,5,2].loc ≈ FV.Point(2.5, 2.0,3.0)
@test m3d.faces.kfaces[3,2,3].loc ≈ FV.Point(2.5, 0.75,4.0)

@test m3d.faces.ifaces[4,2,2].normal ≈ FV.Point(1.0, 0.0, 0.0)
@test m3d.faces.jfaces[3,5,2].normal ≈ FV.Point(0.0, 1.0, 0.0)
@test m3d.faces.kfaces[3,2,3].normal ≈ FV.Point(0.0, 0.0, 1.0)
# %% Test boundary conditions
FV.FaceValue(m2d, [0.0], [0.0], [0.0])
bc2d.left.value[:] = [1.0, 0.0, 0.0, 0.0]
bc2d.left.value
m2d.faces.ifaces[1, 4].area
FV.boundaryConditionTerm(bc2d, FV.FaceValue(m2d, [0.0], [0.0], [0.0]))

# FV.dimensions(m2d)


# map(x -> x.normal, m3d.faces.ifaces)

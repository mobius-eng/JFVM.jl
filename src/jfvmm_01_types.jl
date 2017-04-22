"""
Stores cell (centres) locations
"""
immutable CellLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

CellLocation{T<:Real}(x :: Array{T,1}) = CellLocation(x, zeros(T,1), zeros(T,1))
CellLocation{T<:Real}(x:: Array{T,1}, y:: Array{T,1}) = CellLocation(x, y, zeros(T,1))

"""
Stores the size of each cell in each direction
"""
immutable CellSize{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

CellSize{T<:Real}(x :: Array{T,1}) = CellSize(x, zeros(T,1), zeros(T,1))
CellSize{T<:Real}(x :: Array{T,1}, y :: Array{T,1}) = CellSize(x, y, zeros(T,1))

# TODO: check what exactly is meant by the face location
"""
Stores face locations
"""
immutable FaceLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

FaceLocation{T<:Real}(x::Array{T,1}) = FaceLocation(x, zeros(T,1), zeros(T,1))
FaceLocation{T<:Real}(x::Array{T,1}, y::Array{T,1}) = FaceLocation(x, y, zeros(T,1))

"""
Abstract type shared among all meshes
"""
abstract MeshType

"""
Cartesian 1D mesh
"""
immutable Mesh1D <: MeshType end

const mesh1D = Mesh1D()

"""
Polar 1D mesh (symmetrical in angular direction)
"""
immutable Mesh1DPolar <: MeshType end

const mesh1DPolar = Mesh1DPolar()

"""
Cartesian 2D mesh
"""
immutable Mesh2D <: MeshType end

const mesh2D = Mesh2D()

"""
Cyllindrical 2D mesh (symmetrical in angular dimension)
"""
immutable Mesh2DCylindrical <: MeshType end

const mesh2DCylindrical = Mesh2DCylindrical()

"""
Polar 2D mesh
"""
immutable Mesh2DPolar <: MeshType end

const mesh2DPolar = Mesh2DPolar()

"""
Cartesian 3D mesh
"""
immutable Mesh3D <: MeshType end

const mesh3D = Mesh3D()

"""
Cyllindrical 3D mesh
"""
immutable Mesh3DCylindrical <: MeshType end

const mesh3DCylindrical = Mesh3DCylindrical()

"""
Mesh structure

Fields:

    + meshtype: determines the type of the mesh
    + dim: dimensions
    + cellsize, cellcenters, facecenters: mesh particularities
    + corner: array of indices indicating position of the corner
    + edge: array of indices for the edges

"""
immutable MeshStructure{T<:Real}
  meshtype :: MeshType
  dims :: Array{Int,1}
  cellsize :: CellSize{T}
  cellcenters :: CellLocation{T}
  facecenters :: FaceLocation{T}
  corner :: Array{Int,1}
  edge :: Array{Int,1}
end

"""
Representation of the variable defined on the cell
"""
immutable CellValue{T<:Real}
  domain::MeshStructure{T}
  value::Array{T}
end

"""
Representation of the vector variable defined on the cell
"""
immutable CellVector{T<:Real}
  domain::MeshStructure{T}
  xvalue::Array{T}
  yvalue::Array{T}
  zvalue::Array{T}
end

"""
Representation of the quantity on the face of the cell.

Fields:

- `domain` : the domain (mesh) on which the value is defined
- `xvalue` : the values on the faces with normal (±1,0,0)
- `yvalue` : the values on the faces with normal (0,±1,0)
- `zvalue` : the values on the faces with norma (0,0,±1)
"""
immutable FaceValue{T<:Real}
    domain::MeshStructure{T}
    xvalue::Array{T}
    yvalue::Array{T}
    zvalue::Array{T}
end

"""
Border value

    a ∇φ + b φ = c
"""
type BorderValue{T<:Real}
    a::Array{T}
    b::Array{T}
    c::Array{T}
    periodic::Bool
end

"""
Boundary conditions
"""
type BoundaryCondition{T<:Real}
    domain::MeshStructure{T}
    left::BorderValue{T}
    right::BorderValue{T}
    bottom::BorderValue{T}
    top::BorderValue{T}
    back::BorderValue{T}
    front::BorderValue{T}
end

"""
Stores cell (centres) locations
"""
immutable CellLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

"""
Stores the size of each cell in each direction
"""
immutable CellSize{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

# TODO: check what exactly is meant by the face location
"""
Stores face locations
"""
immutable FaceLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

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

const mesh2D = Mesh2D

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
type CellValue{T<:Real}
  domain::MeshStructure{T}
  value::Array{T}
end

"""
Representation of the vector variable defined on the cell
"""
type CellVector{T<:Real}
  domain::MeshStructure{T}
  xvalue::Array{T}
  yvalue::Array{T}
  zvalue::Array{T}
end

# TODO: check if can make it work for scalar more efficiently
"""
Representation of the quantity on the face of the cell.
"""
type FaceValue{T<:Real}
  domain::MeshStructure{T}
  xvalue::Array{T}
  yvalue::Array{T}
  zvalue::Array{T}
end

"""
"""
type BorderValue{T<:Real}
  a::Array{T}
  b::Array{T}
  c::Array{T}
  periodic::Bool
end

"""
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

# %% These definitions are just for readability in the future


type Cells{T<:Real}
    volumes :: Array{T, 3}
    centroids :: Array{Point{T},3}
end

function Cells{T<:Real}(::Type{T}, Nx :: Int64, Ny :: Int64, Nz :: Int64 = 1)
    volumes = zeros(T, Nx, Ny, Nz)
    c = [Point(T) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
    Cells(volumes, c)
end

type Face{T<:Real}
    area :: T
    normal :: Point{T}
    loc :: Point{T}
end

function createFace{T<:Real}(a :: Point{T}, b :: Point{T}, nd :: Point{T})
    loc = (a+b) ./ (2one(T))
    p = a - b
    area = norm(a - b)
    n = (Point(p.y, -p.x) ./ area)
    normal = n .* sign(n ⋅ nd)
    Face(area, normal, loc)
end

function createFace{T<:Real}(a :: Point{T}, b :: Point{T}, c :: Point{T}, d ::Point{T}, nd :: Point)
    c1 = centroid(d, a, b)
    c2 = centroid(b, c, d)
    a1 = parallelogramarea(b - a, d - a) / 2
    a2 = parallelogramarea(b - c, d - c) / 2
    area = a1 + a2
    loc = (a1/area) .* c1 + (a2/area) .* c2
    p = (b - a) × (d - a)
    n = p ./ norm(p)
    normal = n .* sign(n ⋅ nd)
    Face(area, normal, loc)
end

type Faces{T<:Real}
    ifaces :: Array{Face{T}}
    jfaces :: Array{Face{T}}
    kfaces :: Array{Face{T}}
end

function Faces{T<:Real}(::Type{T}, Ni, Nj, Nk=1)
    ifaces = [Face(zero(T), Point(T), Point(T)) for i in 1:(Ni+1), j in 1:Nj, k in 1:Nk]
    jfaces = [Face(zero(T), Point(T), Point(T)) for i in 1:Ni, j in 1:(Nj+1), k in 1:Nk]
    kfaces = [Face(zero(T), Point(T), Point(T)) for i in 1:Ni, j in 1:Nj, k in 1:(Nk+1)]
    Faces(ifaces, jfaces, kfaces)
end

# %% Mesh structure
"""
Mesh structure

Fields:
    + dim: vector dimensions for each index (length = N)
    + cellsize, cellcenters, facecenters: mesh particularities
"""
type MeshStructure{T<:Real}
    dims        :: Vector{Int64}
    vertices    :: Array{Point{T}}
    cells       :: Cells{T}
    faces       :: Faces{T}
end

dimensions(m :: MeshStructure) = m.dims

# %% Cell and Face values

"""
Representation of the variable defined on the cell
"""
type CellValue{T1<:Real, T2}
    domain :: MeshStructure{T1}
    value  :: Array{T2}
end

"""
Representation of the vector variable defined on the cell
"""
type CellVector{T<:Real}
    domain :: MeshStructure{T}
    vector :: Array{Point{T}}
end

"""
Representation of the quantity on the face of the cell.

Fields:

- `domain` : the domain (mesh) on which the value is defined
- `value` : the values on the faces with normal (±1,0,0)
"""
type FaceValue{T<:Real, T2}
    domain :: MeshStructure{T}
    ivalue :: Array{T2}
    jvalue :: Array{T2}
    kvalue :: Array{T2}
end

type FaceVector{T<:Real}
    domain :: MeshStructure{T}
    ivalue :: Array{Point{T}}
    jvalue :: Array{Point{T}}
    kvalue :: Array{Point{T}}
end

# %% Boundary condition representation
"""
Border value: either flux or value
- `isflux` contains the array of trues where flux is specified, otherwise -- falses
- `value` array either values or flux ⋅ (area norm)
"""
type BorderValue{T}
    isflux :: Array{Bool}
    value :: Array{T}
end

"""
Boundary conditions
"""
immutable BoundaryCondition{T<:Real}
    domain::MeshStructure{T}
    left :: BorderValue{T}
    right :: BorderValue{T}
    bottom :: BorderValue{T}
    top :: BorderValue{T}
    back :: BorderValue{T}
    front :: BorderValue{T}
end

function BoundaryCondition(;kwargs...)
    d = Dict(kwargs)
    BoundaryCondition(d[:domain], d[:left], d[:right], d[:bottom], d[:top], d[:back], d[:front])
end

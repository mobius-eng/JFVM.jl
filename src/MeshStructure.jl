__precompile__()
module MeshStructure

# %% Preambule
using JFVMM.Geometry

export Node, Face, Cell, Mesh
export meshindex, centroid, dimensions, createMesh, createRectMesh
export facevector, facearea, faceowneri, faceneighbouri, centroid_correction
export cellifacesi, celljfacesi
export cells, nodes, ifaces, jfaces, boundaryfaces, globalindex, localindex

# %% Node
"""
Single vertex of the mesh

Fields:
- `centroid` is the space position of the vertex
- `meshindex` is the position of the vertex within the mesh
"""
immutable Node
    meshindex :: CartesianIndex{2}
    centroid :: Point2D
end

"""
Returns the centroid of the mesh constituent (node, face, element)
"""
centroid(n :: Node) = n.centroid

"""
Return the index of the mesh constituent within the mesh
"""
meshindex(n :: Node) = n.meshindex

# %% Face
"""
Represents the face of the mesh

Fields:
- `meshindex` is the position of the face within the mesh
- `centroid` is the location of the centroid of the face
- `Sf` is face vector: Sf = Sx i + Sy j + Sz k, such as Sf || normal to face and |Sf| = area
- `area` is the area of the face
- `gf` is "geometric factor" of the face for linear interpolation
- `dCF` is the distance (vector) CF connecting centroids of owner cell (C) and the neighbour cell (F)
- `dCf` is the distance (vector) Cf connecting centroid of the owner cell (C) and the centroid of face
"""
immutable Face
    isboundary :: Bool
    meshindex :: CartesianIndex{2}
    direction :: CartesianIndex{2}
    Sf :: Point2D
    centroid :: Point2D
    falsecentroid :: Point2D
    ff :: Point2D
    CF :: Point2D
    Cf :: Point2D
    Ef :: Point2D
    gf :: Float64
    area :: Float64
    gDiff :: Float64
end

meshindex(f :: Face) = f.meshindex
centroid(f :: Face) = f.centroid

"""
Vector in the normal direction to the face such that its abs value = face area
"""
facevector(f :: Face) = f.Sf

facearea(f :: Face) = f.area

faceowneri(f :: Face) = meshindex(f)
faceneighbouri(f :: Face) = meshindex(f) - f.direction

centroid_correction(f :: Face) = f.ff

#%% Cell
"""
Cell (element) of the mesh

Fields:
- `meshindex` is the position of the cell in the mesh
- `neighbours` is the list of neighboring cells
- `faces` is the list of cell faces, ordered in the same way as `neighbours`, i.e. `faces[i]` is the face between this cell and `neighbours[i]` cell
- `nodes` is the list of nodes making up the cell.
- `volume` is the volume of the cell
- `faceSign` is the indicator if this cell is the owner of its faces: +1 if it is and -1 otherwise
"""
immutable Cell
    meshindex :: CartesianIndex{2}
    volume :: Float64
    centroid :: Point2D
    ifaceSign :: Dict{CartesianIndex{2},Float64}
    jfaceSign :: Dict{CartesianIndex{2},Float64}
end

meshindex(c :: Cell) = c.meshindex
centroid(c :: Cell) = c.centroid

function cellifacesi(c :: Cell)
    I = meshindex(c)
    Ii = CartesianIndex(1,0)
    (I, I+Ii)
end

function celljfacesi(c :: Cell)
    I = meshindex(c)
    Ij = CartesianIndex(0,1)
    (I, I+Ij)
end

# %% Mesh
type Mesh
    dims :: Tuple{Int64,Int64}
    cells :: Array{Cell,2}
    nodes :: Array{Node,2}
    ifaces :: Array{Face,2}
    jfaces :: Array{Face, 2}
end

dimensions(m :: Mesh) = m.dims
Base.size(m :: Mesh) = m.dims

function ifaces(m :: Mesh, internal = false)
    if internal
        @view m.ifaces[2:end-1,:]
    else
        @view m.ifaces[:,:]
    end
end

function jfaces(m :: Mesh, internal = false)
    if internal
        @view m.jfaces[:,2:end-1]
    else
        @view m.jfaces[:,:]
    end
end

function cells(m :: Mesh, internal = false)
    if internal
        @view m.cells[2:end-1,2:end-1]
    else
        @view m.cells[:,:]
    end
end

nodes(m :: Mesh) = m.nodes

"""
Returns global index (in the system matrix) for local indices (i,j)
or expressed as Cartesian index (`CartesianIndex(i,j)`).
"""
globalindex(m :: Mesh, i, j) = sub2ind(size(m), i, j)
globalindex(m :: Mesh, I :: CartesianIndex{2}) = sub2ind(size(m), I.I...)

"""
Returns local indices from global index `n` as `CartesianIndex(i,j)`
"""
localindex(m :: Mesh, n) = CartesianIndex(ind2sub(size(m), n)...)

function boundaryfaces(m :: Mesh)
    [(@view m.ifaces[1,:]), (@view m.jfaces[:,end]), (@view m.ifaces[end,:]), (@view m.jfaces[:,end])]
end

# %% Compute cell volume and centroid
# 2D case
function computeCellVolumeCentroid(a :: Point2D, b :: Point2D, c :: Point2D, d :: Point2D)
    a1 = triangleArea(a, b, c)
    a2 = triangleArea(a, c, d)
    c1 = triangleCentroid(a,b,c)
    c2 = triangleCentroid(a,c,d)
    a = a1+a2
    (a,(a1/a)*c1 + (a2/a)*c2)
end

# %% Compute face area
# 2D case
function computeFaceAreaCentroid(a :: Point2D, b :: Point2D)
    p = b-a
    area = norm(p)
    Sf = Point2D(p.y, -p.x)
    (area, Sf, (a+b)/2)
end

# %% Mesh from the array of points
#2D case
function createFace(points, nodes :: Array{Node,2}, cells :: Array{Cell, 2}, ind :: CartesianIndex{2}, facedir :: CartesianIndex{2})
    a, b = points
    Ni, Nj = size(cells)
    area, Sfu, f = computeFaceAreaCentroid(a, b)
    islowerboundary = false
    isupperboundary = false
    # determine if the face belong to boundary
    if any(x -> x == 0, (ind - facedir).I)
        islowerboundary = true
    end
    if any(x -> x > 0, (ind - CartesianIndex(Ni, Nj)).I)
        isupperboundary = true
    end
    # centroid of the owner
    if !isupperboundary
        # if not on upper boundary: face's owner has the same index
        C = centroid(cells[ind])
    else
        # otherwise, it is "previous" cell in the direction of face's normal
        C = centroid(cells[ind-facedir])
    end
    Cf = f - C
    # correct Sf sign
    sgn = sign(Cf ⋅ Sfu)
    Sf = sgn * Sfu
    # Neighbour information and face weigh factor
    if isupperboundary || islowerboundary
        # if on the boundary, define CF as Cf
        CF = Cf
        gf = 1.0
        falsecentroid = C
        ff = Point2D(0.0,0.0)
    else
        # neighbour for inner faces is "previous" cell
        indF = ind - facedir
        F = centroid(cells[indF])
        CF = F - C
        falsecentroid = C + ((Cf ⋅ CF) / (CF ⋅ CF)) * CF
        ff = f - falsecentroid
        # weighting factor for false centroid point f'
        gf = norm(falsecentroid - C) / norm(CF)
        # compure weighing factor (robust way)
        # ef = Sf / area
        # fF = F - f
        # Cf_ef = Cf ⋅ ef
        # gf = Cf_ef / (Cf_ef + (fF ⋅ ef))
    end
    # Compute Ef using over-relaxed method
    gDiff = area/(Sf ⋅ CF)
    Ef = gDiff * CF
    Face(islowerboundary || isupperboundary, ind, facedir, Sf, f, falsecentroid, ff, CF, Cf, Ef, gf, area, gDiff)
end

function cellifaceSign(I :: CartesianIndex{2}, cellmax)
    Ni, Nj = cellmax
    I2 = I + CartesianIndex(1,0)
    d = Dict(I => 1.0)
    if I.I[1] == Ni
        d[I2] = 1.0
    else
        d[I2] = -1.0
    end
    d
end

function celljfaceSign(I :: CartesianIndex{2}, cellmax)
    Ni, Nj = cellmax
    I2 = I + CartesianIndex(0,1)
    d = Dict(I => 1.0)
    if I.I[2] == Nj
        d[I2] = 1.0
    else
        d[I2] = -1.0
    end
    d
end


function getCellPointsIndices(I :: CartesianIndex{2})
    di = CartesianIndex(1,0)
    dj = CartesianIndex(0,1)
    [I, I+di, I+di+dj, I+dj]
end

function getFacePointsIndices(I :: CartesianIndex{2}, dir :: CartesianIndex{2})
    di = CartesianIndex(dir.I[2], dir.I[1])
    [I, I + di]
end

function cycletuple(t, n = 1)
    (t[end-n+1:end]..., t[1:end-n]...)
end
#
# function getFacePointsIndices(I :: CartesianIndex{3}, dir)
#     dir1 = CartesianIndex(cycletuple(dir.I,1)...)
#     dir2 = CartesianIndex(cycletuple(dir.I,2)...)
#     [I, I+dir1, I+dir1+dir2, I+dir2]
# end

function createMesh(points :: Array{Point2D,2})
    nmax = [size(points)...]
    cmax = nmax - 1
    dim = length(nmax)
    d = Tuple([zeros(Int64, dim-1); 1])
    dirs = [CartesianIndex(cycletuple(d, i)...) for i in 1:dim]
    cend = CartesianIndex(cmax...)
    cbeg = CartesianIndex(ones(Int64,dim)...)

    nodesRange = CartesianRange(size(points))
    cellsRange = CartesianRange(Tuple(cmax))
    nodes = [Node(I, points[I]) for I in nodesRange]
    cells = [Cell(I, (computeCellVolumeCentroid(points[getCellPointsIndices(I)]...)..., cellifaceSign(I, cmax), celljfaceSign(I,cmax))...) for I in cellsRange]

    faces = [
        [createFace(points[getFacePointsIndices(I, dirs[id])], nodes, cells, I, dirs[id]) for I in CartesianRange(cbeg, cend+dirs[id])]
        for id in 1:dim
        ]
    Mesh((cmax...), cells, nodes, faces...)
end

function createRectMesh(width :: Float64, height :: Float64, Nw, Nh)
    xs = linspace(0.0, width, Nw+1)
    ys = linspace(0.0, height, Nh+1)
    points = [Point2D(x,y) for x in xs, y in ys]
    createMesh(points)
end


end

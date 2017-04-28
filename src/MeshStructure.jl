__precompile__()
module MeshStructure

# %% Preambule
using JFVMM.Geometry

export Node, Face, Cell
export meshindex, centroid

# %% Node
"""
Single vertex of the mesh

Fields:
- `centroid` is the space position of the vertex
- `meshindex` is the position of the vertex within the mesh
"""
immutable Node{T<:Real}
    meshindex :: Tuple{Int64, Int64, Int64}
    centroid :: Point{T}
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
immutable Face{T<:Real}
    meshindex :: Tuple{Int64,Int64,Int64}
    direction :: Symbol
    area :: T
    Sf :: Point{T}
    centroid :: Point{T}
    gf :: T
    dCF :: Point{T}
    dCf :: Point{T}
end

meshindex(f :: Face) = f.meshindex
centroid(f :: Face) = f.centroid

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
immutable Cell{T<:Real}
    meshindex :: Tuple{Int64,Int64,Int64}
    volume :: T
    centroid :: Point{T}
end

meshindex(c :: Cell) = c.meshindex
centroid(c :: Cell) = c.centroid

# %% Compute cell volume and centroid
# 2D case
function computeCellVolumeCentroid{T<:Real}(a :: Point{T}, b, c, d)
    a1, Sf1 = triangleArea(a, b, c)
    a2, Sf2 = triangleArea(a, c, d)
    c1 = triangleCentroid(a,b,c)
    c2 = triangleCentroid(a,c,d)
    a = a1+s2
    (a,(a1/a)*c1 + (a2/a)*c2)
end

# 3D case
function computeCellVolumeCentroid{T<:Real}(a :: Point{T}, b, c, d, e, f, g, h)
    # geometric centre of the cell
    xg = (a+b+c+d+e+f+g+h) / 8
    # compute volume and the centroid of each pyramid
    # Pyramid abcd_xg
    v1, c1 = pyramidVolumeCentroid(a,b,c,d,xg)
    # Pyramid abfe_xg
    v2, c2 = pyramidVolumeCentroid(a,b,f,e,xg)
    # Pyramid efgh_xg
    v3, c3 = pyramidVolumeCentroid(e,f,g,h,xg)
    # Pyramid hgcd_xg
    v4, c4 = pyramidVolumeCentroid(h,g,c,d,xg)
    # Pyramid adhe_xg
    v5, c5 = pyramidVolumeCentroid(a,d,h,e,xg)
    # Pyramid cbfg_xg
    v6, c6 = pyramidVolumeCentroid(c,b,f,g,xg)
    # combine
    v = v1+v2+v3+v4+v5+v6
    (v, (v1/v)*c1 + (v2/v)*c2 + (v3/v)*c3 + (v4/v)*c4 + (v5/v)*c5 + (v6/v)*c6)
end

# %% Compute face area
# 2D case
function computeFaceAreaCentroid{T<:Real}(a :: Point{T}, b)
    p = b-a
    area = norm(p)
    Sf = Point(p.y, - p.x)
    (area, Sf, (a+b)/2)
end

# 3D case
function computeFaceAreaCentroid{T<:Real}(a :: Point{T}, b, c, d)
    # Triangle abc
    a1, Sf1 = triangleArea(a,b,c)
    c1 = triangleCentroid(a,b,c)
    # Triangle acd
    a2, Sf2 = triangleArea(a,c,d)
    c2 = triangleCentroid(a,c,d)
    # combine
    a = a1+a2
    e = Sf1 / a1
    (a, a*e, (a1/a)*c1 + (a2/a)*c2)
end

# %% Mesh from the array of points
#2D case
function createFace(a, b, nodes, cells, i, j, facedir)
    ind = (i,j,0)
    di, dj = facedir
    Npi, Npj = size(nodes)
    area, Sf, centroid = computeFaceAreaCentroid(a, b)
end

function createMesh{T<:Real}(points :: Array{Point{T},2})
    Npi, Npj = size(points)
    Ni = Npi-1
    Nj = Npj-1
    nodes = [Node((i,j,0), point[i,j]) for i ∈ 1:Npi, j ∈ 1:Npj]
    cells = [Cell((i,j,0), computeCellVolumeCentroid(points[i,j], points[i+1,j], points[i+1,j+1],points[i,j+1])...) for i ∈ 1:Ni, j ∈ 1:Nj]
    ifaces = [cretateFace(points[i,j], points[i,j+1], nodes, cells, i, j, (1,0)) for i ∈ 1:Npi, j ∈ 1:Nj]
end

end

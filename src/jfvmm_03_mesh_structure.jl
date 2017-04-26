# %% General 2D mesh
centroid{T<:Real}(a :: Point{T}, b, c) = (a .+ b .+ c) ./ (3one(T))
centroid{T<:Real}(a :: Point{T},b,c,d) = (a .+ b .+ c .+ d) ./ (4one(T))

vecvol(a, b, c) = abs(det([a.x a.y a.z; b.x b.y b.z; c.x c.y c.z]))
vecvol(a, b) = abs(a.x*b.y - a.y * b.x)


function parallelogramarea(a,b)
    norm(a × b)
end

function centroidvol(a, b, c, d)
    c1 = centroid(a,b,c)
    c2 = centroid(a,d,c)
    v1 = vecvol(b .- a, c .- a) / 2
    v2 = vecvol(d .- a,c .- a) / 2
    v = v1 + v2
    ((v1/v) .* c1 + (v2/v) .* c2, v)
end

function pyramidvol(a,b,c,d,p)
    v1 = vecvol(a - p, b - p, c - p) / 6
    v2 = vecvol(a -p, d - p, c - p) / 6
    (v1, v2)
end

function centroidvol(a, b, c, d, e, f, g, h)
    # cannot compute centroid and volume directly. Needs to split to tetrahedra
    # Split as follows: compute geometric mean and construct pyramids from al the faces and
    # the mean. Split the pyramid to 2 tetrahedra via base diagonal.
    # mean centre of cuboid
    p = (a + b + c + d + e + f + g + h) / 8
    # splitting cuboid to tetrahedra
    # back plane
    v_abcp, v_acdp = pyramidvol(a,b,c,d,p)
    c_abcp = centroid(a,b,c,p)
    c_acdp = centroid(a,c,d,p)
    # front plane
    v_efgp,  v_ehgp = pyramidvol(e, f, g, h, p)
    c_efgp = centroid(e,f,g,p)
    c_ehgp = centroid(e,h,g,p)
    # left plane
    v_adhp, v_ahep = pyramidvol(a,d,h,e,p)
    c_adhp = centroid(a,d,h,p)
    c_ahep = centroid(a,h,e,p)
    # right plane
    v_bcgp, v_bgfp = pyramidvol(b,c,g,f,p)
    c_bcgp = centroid(b,c,g,p)
    c_bgfp = centroid(b,g,f,p)
    # bottom plane
    v_abfp, v_afep = pyramidvol(a,b,f,e,p)
    c_abfp = centroid(a,b,f,p)
    c_afep = centroid(a,f,e,p)
    # top plane
    v_dcgp, v_dhgp = pyramidvol(d,c,g,h,p)
    c_dcgp = centroid(d,c,g,p)
    c_dhgp = centroid(d,h,g,p)
    v = v_abcp + v_acdp + v_efgp + v_ehgp + v_adhp + v_ahep + v_bcgp + v_bgfp + v_abfp + v_afep + v_dcgp + v_dhgp
    c = (v_abcp/v) .* c_abcp + (v_acdp/v) .* c_acdp + (v_efgp/v) .* c_efgp + (v_ehgp/v) .* c_ehgp + (v_adhp/v) .* c_adhp + (v_ahep/v) .* c_ahep + (v_bcgp/v) .* c_bcgp + (v_bgfp/v) .* c_bgfp + (v_abfp/v) .* c_abfp + (v_afep/v) .* c_afep + (v_dcgp/v) .* c_dcgp + (v_dhgp/v) .* c_dhgp
    (c,v)
end

function centroidarea(a, b, c, d)
    c1 = centroid(d, a, b)
    c2 = centroid(b, c, d)
    a1 = parallelogramarea(b .- a, d .- a) / 2
    a2 = parallelogramarea(b .- c, d .- c) / 2
    at = a1 + a2
    ((a1/at) .* c1 + (a2/at) .* c2, at)
end

function rectpoints(vertices, i, j, k=1)
    va = vertices[i,j,k]
    vb = vertices[i+1,j,k]
    vc = vertices[i+1,j+1,k]
    vd = vertices[i,j+1,k]
    (va, vb, vc, vd)
end

function createMesh2D{T<:Real}(vertices :: Array{Point{T}, 2})
    dims = [size(vertices)[1:2]...;] .- 1
    cells = Cells(T, dims...)
    faces = Faces(T, dims...)
    for j in 1:dims[2]
        for i in 1:dims[1]
            va, vb, vc, vd = rectpoints(vertices, i, j)
            c, v = centroidvol(va, vb, vc, vd)
            cells.centroids[i,j] = c
            cells.volumes[i,j] = v
            faces.ifaces[i,j] = createFace(va, vd, vc - vd)
            faces.jfaces[i,j] = createFace(va, vb, vd - va)
        end
    end
    # loop over right boundary for facelocs and faceareas
    for j in 1:dims[2]
        va = vertices[end,j]
        vd = vertices[end,j+1]
        faces.ifaces[end,j] = createFace(va, vd, faces.ifaces[end-1,j].normal)
    end
    #loop over top boundary
    for i in 1:dims[1]
        va = vertices[i,end]
        vb = vertices[i+1,end]
        faces.jfaces[i,end] = createFace(va, vb, faces.jfaces[i, end-1].normal)
    end
    MeshStructure(dims, vertices, cells, faces)
end

function createMesh3D{T<:Real}(vertices :: Array{Point{T}, 3})
    dims = [size(vertices)[1:3]...;] - 1
    cells = Cells(T, dims...)
    faces = Faces(T, dims...)
    for k in 1:dims[3]
        for j in 1:dims[2]
            for i in 1:dims[1]
                va, vb, vc, vd = rectpoints(vertices, i, j, k)
                ve, vf, vg, vh = rectpoints(vertices, i, j, k+1)
                c, v = centroidvol(va, vb, vc, vd, ve, vf, vg, vh)
                cells.centroids[i,j,k] = c
                cells.volumes[i,j,k] = v
                faces.ifaces[i,j,k] = createFace(va, ve, vh, vd, vg - vh)
                faces.jfaces[i,j,k] = createFace(va, vb, vf, ve, vd - va)
                faces.kfaces[i,j,k] = createFace(va, vb, vc, vd, ve - va)
            end
        end
    end
    # loop over front boundary k = Nz+1
    for j in 1:dims[2]
        for i in 1:dims[1]
            va, vb, vc, vd = rectpoints((@view vertices[:,:,end]), i, j)
            faces.kfaces[i,j,end] = createFace(va, vb, vc, vd, faces.kfaces[i,j,end-1].normal)
        end
    end
    # loop over top boundary j = Ny + 1
    for k in 1:dims[3]
        for i in 1:dims[1]
            va, vb, vc, vd = rectpoints((@view vertices[:, end, :]), i, k)
            faces.jfaces[i,end,k] = createFace(va, vb, vc, vd, faces.jfaces[i,end-1,k].normal)
        end
    end
    # loop over right boundary i = Nx + 1
    for k in 1:dims[3]
        for j in 1:dims[2]
            va, vb, vc, vd = rectpoints((@view vertices[end,:,:]), j, k)
            faces.ifaces[end,j,k] = createFace(va, vb, vc, vd, faces.ifaces[end-1,j,k].normal)
        end
    end
    MeshStructure(dims, vertices, cells, faces)
end

# %% ========================= 2D CARTESIAN MESH ==========================
"""
m = createMesh2D(Nx::Int, Ny::Int, Width::Real, Height::Real)
It creates a uniform mesh on a 2D Cartesian domain.

Inputs:

   + Nx: number of cells in the x-direction (integer)
   + Ny: number of cells in the y-direction (integer)
	 + Width: length of the domain in the x-direction (Real)
   + Height: length of the domain in the y-direction (Real)

Outputs:

   + m: a mesh structure

Usage:

Nx=10
Ny=15
Lx=1.0
Ly=2.5
m=createMesh2D(Nx, Ny, Lx, Ly)
"""
function createMesh2D{T<:Real}(width::T, height::T, Nw::Int, Nh::Int)
    verticesx = linspace(zero(T), width, Nw+1) |> collect
    verticesy = linspace(zero(T), height, Nh+1) |> collect
    vertices = [Point(x,y) for x in verticesx, y in verticesy]
    createMesh2D(vertices)
end

# %% ========================== 3D CARTESIAN MESH ============================
"""
Creates a uniform 3D mesh:

m = createMesh3D(Nx, Ny, Nz, Width, Height, Depth)

Inputs:

    + Nx (int): number of cells in horizontal direction
    + Ny (int): number of cells in vertical direction
    + Nz (int): number of cells in depth direction
    + Width, Height, Depth (reals): size of each dimension

Output:

    + m: a mesh structure

"""
function createMesh3D{T<:Real}(width::T, height::T, depth::T, Nw::Int64, Nh::Int64, Nd::Int64)
    verticesx = linspace(zero(T), width, Nw+1) |> collect
    verticesy = linspace(zero(T), height, Nh+1) |> collect
    verticesz = linspace(zero(T), depth, Nd+1) |> collect
    vertices = [Point(x,y,z) for x in verticesx, y in verticesy, z in verticesz]
    createMesh3D(vertices)
end

globalindex(mesh, i, j) = (i-1)*dimensions(mesh)[2] + j

function globalindex(mesh, i, j, k)
    Ni, Nj, Nk = dimensions(m)
    (i-1)*Nj*Nk + (j-1)*Nk + k
end

function globalneighbours(mesh, n, dirs = [(-1,1), (-1,1), (-1,1)])
    dims = dimensions(mesh)
    if length(dims) == 2
        # 2D case
        Nj = dims[2]
        [map(p -> n + p*Nj, dirs[1]), map(p -> n + p, dirs[2])]
    else
        # 3D
        Nj,Nk = dims[2:3]
        [map(p -> n + p*Nj*Nk, dirs[1]), map(p -> n + p*Nk, dirs[2]), map(p -> n + p, dirs[3])]
    end
end

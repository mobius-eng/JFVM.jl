# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# %% ============== Utility =========================================
"""
!!! INTERNAL !!!
Calculates cell sizes
"""
function cellsize{T<:Real}(faceloc :: Vector{T})
    [
        faceloc[2]-faceloc[1];
        faceloc[2:end]-faceloc[1:end-1];
        faceloc[end]-faceloc[end-1]
    ]
end

cellsize{T<:Real}(δ :: T, N :: Int64) = δ .* ones(T,N+2)

"""
!!! INTERNAL !!!
Calculates cell locations
"""
cellloc{T<:Real}(faceloc :: Vector{T}) = (faceloc[2:end] .+ faceloc[1:end-1]) ./ 2

cellloc{T<:Real}(δ :: T, N :: Int64) = [1:N;] .* δ .- δ ./ 2

faceloc{T<:Real}(δ :: T, N :: Int64) = [0:N] .* δ

# %% ========== Internal 1D mesh creation ===========================
"""
!!! INTERNAL !!!
Creates a unifrom 1D mesh
"""
function meshStructure1D{T<:Real}(meshtype :: MeshType, cellsizev :: T, cellnum :: Int64)
    MeshStructure(meshtype,
        [cellnum],
        CellSize(cellsize(cellsizev), cellnum),
        CellLocation(cellloc(cellsizev, cellnum)),
        FaceLocation(faceloc(cellsizev, cellnum)),
        [1], [1])
end

function meshStructure1D{T<:Real}(meshtype :: MeshType, faceloc :: Vector{T})
    MeshStructure(meshtype,
    		[length(faceloc)-1],
    		CellSize(cellsize(faceloc)),
    		CellLocation(cellloc(faceloc)),
    		FaceLocation(faceloc),
    		[1], [1])
end

# %% ====================== 1D CARTESIAN MESH =======================
"""
m = createMesh1D(Nx::Int, Width::Real)
It creates a uniform mesh on a 1D Cartesian domain.

Inputs:

   + Nx: number of cells in the domain (integer)
	 + Width: width or length of the domain (Real)

Outputs:

   + m: a mesh structure

Usage:

Nx=10
Lx=1.0
m=createMesh1D(Nx, Lx)
"""
createMesh1D{T<:Real}(Nx::Int, Width::T) = meshStructure1D(mesh1D, Width/Nx, Nx)

"""
m = createMesh1D(facelocationX::Array{Real,1})
It creates a non/uniform mesh on a 1D Cartesian domain.

Inputs:

   + facelocationX: location of the the cell boundaries on
	 the x-axis.

Outputs:

   + m: a mesh structure

Usage:

	x= [0.0, 1.0, 1.4, 2.5, 4.1, 6.0, 10.0]
	m=createMesh1D(x)
"""
createMesh1D{T<:Real}(facelocationX::Array{T,1}) = meshStructure1D(mesh1D, facelocationX)

# %% ================= 1D CYLINDRICAL MESH ==========================
"""
m = createMeshCylindrical1D(Nr::Int, Radius::Real)
It creates a uniform mesh on a 1D Radial domain.

Inputs:

   + Nr: number of cells in the domain (integer)
	 + Radius: Radius of the domain (Real)

Outputs:

   + m: a mesh structure

Usage:

Nr=10
R=1.0
m=createMesh1D(Nr, R)
"""
createMeshCylindrical1D{T<:Real}(Nr::Int, Radius::T) =
    meshStructure1D(mesh1DPolar, Radius/Nr, Nr)

"""
m = createMeshCylindrical1D(facelocationR::Array{Real,1})
It creates a non/uniform mesh on a 1D radial domain.

Inputs:

   + facelocationX: location of the the cell boundaries on
	 the r-axis.

Outputs:

   + m: a mesh structure

Usage:

	r= [0.1, 1.0, 1.4, 2.5, 4.1, 6.0, 10.0]
	m=createMeshCylindrical1D(r)
"""
createMeshCylindrical1D{T<:Real}(facelocationR::Array{T,1}) =
    meshStructure1D(mesh1DPolar, facelocationR)

# %% ================= Internal 2D mesh generation ========================
function meshStructure2D{T<:Real}(meshtype :: MeshType, dx :: T, dy :: T, Nx :: Int64, Ny :: Int64)
    G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
    MeshStructure(meshtype,
    	[Nx, Ny],
    	CellSize(cellsize(dx, Nx), cellsize(dy, Ny)),
    	CellLocation(cellloc(dx, Nx), cellloc(dy, Ny)),
    	FaceLocation(faceloc(dx, Nx),faceloc(dy, Ny)),
    	G[[1,end],[1,end]][:],
    	[1])
end

function meshStructure2D{T<:Real}(meshtype :: MeshType, faceloc1 :: Vector{T}, faceloc2 :: Vector{T})
    Nx = length(faceloc1)-1
    Ny = length(faceloc2)-1
    G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
    MeshStructure(meshtype,
    	[Nx, Ny],
    	CellSize(cellsize(faceloc1), cellsize(faceloc2)),
    	CellLocation(cellloc(faceloc1), cellloc(faceloc2)),
    	FaceLocation(faceloc1, faceloc2),
    	G[[1,end],[1,end]][:],
    	[1])
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
createMesh2D{T<:Real}(Nx::Int, Ny::Int, Width::T, Height::T) =
    meshStructure2D(mesh2D, Width / Nx, Height / Ny, Nx, Ny)

"""
Creates a mesh on a 2D Cartesian domain with specified face locations

m = createMesh2D(faceLocationX, faceLocationY)

Inputs:

   + faceLocationX (Vector of reals): face locations along horizontal axis
   + faceLocationY (Vector of reals): face locations along vertical axis

Outputs:

   + m: a mesh structure

Usage:

fX = [0.0; 0.05; 0.1:0.1:1.0;]
fY = [0.0; 0.05; 0.12; 0.2:0.2:3.0;]

m=createMesh2D(fx, fy)
"""
createMesh2D{T<:Real}(facelocationX::Array{T,1}, facelocationY::Array{T,1}) =
    meshStructure2D(mesh2D, facelocationX, facelocationY)

# %% ========================= 2D RADIAL MESH ==============================
"""
Builds a uniform radial (polar coordinates) 2D mesh:

m = createMeshRadial2D(Nr, Ntheta, Radius, Angle)

Inputs:

    + Nr (Integer): number of cells in r (radial) direction
    + Ntheta (Integer): number of cells in angular direction
    + Radius (Real): domain length in r (radial) direction
    + Angle (Real): angular domain size (in radians)

Outputs:

    + m: a mesh structure

"""
function createMeshRadial2D(Nr::Int, Ntheta::Int, Radius::Real, Angle::Real)
    if Angle>2*pi
    	Angle = 2*pi
    	warn("The domain size adjusted to match a maximum of 2π.")
    end
    dr = Radius/Nr
    dtheta = Angle/Ntheta
    meshStructure2D(mesh2DPolar, dr, dtheta, Nr, Ntheta)
end

"""
Creates a mesh on a 2D radial (polar coordinates) domain with specified face
locations

m = createMesh2D(faceLocationR, faceLocationTheta)

Inputs:

   + faceLocationR (Vector of reals): face locations along radial axis
   + faceLocationTheta (Vector of reals): face locations along angular axis

Outputs:

   + m: a mesh structure

Usage:

fR = [0.0; 0.05; 0.1:0.1:1.0;]
fΘ = [0.0; 0.05; 0.12; 0.2:0.2:π;]

m=createMesh2D(fR, fΘ)
"""
function createMeshRadial2D{T<:Real}(facelocationR::Array{T,1}, facelocationTheta::Array{T,1})
    if facelocationTheta[end]>2.0*pi
    	facelocationTheta = facelocationTheta/facelocationTheta[end]*(2.0*pi)
    	warn("The domain size adjusted to match a maximum of 2π.")
    end
    meshStructure2D(mesh2DPolar, facelocationR, facelocationTheta)
end

# %% ===================== 2D CYLINDRICAL MESH ==========================
"""
Builds a uniform cylindrical 2D mesh

m = createMeshCylindrical2D(Nr, Nz, Radius, Height)

Inputs:
    + Nr (Integer): number of cells in r (radial) direction
    + Nz (Integer): number of cells in z (axial, vertical) direction
    + Radius (Real): radius of the cylinder
    + Height (Real): height of the cylinder

Output:

    + m: a mesh structure
"""
createMeshCylindrical2D(Nr::Int, Nz::Int, Radius::Real, Height::Real) =
    meshStructure2D(mesh2DCylindrical, Radius / Nr, Height / Nz, Nr, Nz)

"""
Creates a mesh on a 2D cylindrical domain with specified face
locations

m = createMesh2D(faceLocationR, faceLocationZ)

Inputs:

   + faceLocationR (Vector of reals): face locations along radial axis
   + faceLocationZ (Vector of reals): face locations along vertical axis

Outputs:

   + m: a mesh structure

Usage:

fR = [0.0; 0.05; 0.1:0.1:1.0;]
fZ = [0.0; 0.05; 0.12; 0.2:0.2:3.0;]

m=createMesh2D(fR, fZ)
"""
createMeshCylindrical2D{T<:Real}(facelocationR::Array{T,1}, facelocationY::Array{T,1}) =
    meshStructure2D(mesh2DCylindrical, facelocationR, facelocationY)

# %% ======================= Internal 3D mesh generation ====================
function meshStructure3D{T<:Real}(meshtype :: MeshType, dx :: T, dy :: T, dz :: T, Nx, Ny, Nz)
    G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
    MeshStructure(meshtype,
    	[Nx, Ny, Nz],
    	CellSize(cellsize(dx, Nx), cellsize(dy, Ny), cellsize(dz, Nz)),
    	CellLocation(cellloc(dx, Nx), cellloc(dy, Ny), cellloc(dz, Nz)),
    	FaceLocation(faceloc(dx, Nx), faceloc(dy, Ny), faceloc(dz, Nz)),
    	G[[1,end],[1,end],[1,end]][:],                         # corners
    	[                                                      # edges
            G[[1, end], [1, end], 2:Nz+1][:];
        	G[[1, end], 2:Ny+1, [1, end]][:];
        	G[2:Nx+1, [1, end], [1, end]][:]
        ])
end

function meshStructure3D{T<:Real}(meshtype :: MeshType, facelocationX::Array{T,1}, facelocationY::Array{T,1}, facelocationZ::Array{T,1})
    Nx = length(facelocationX)-1
    Ny = length(facelocationY)-1
    Nz = length(facelocationZ)-1
    G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
    MeshStructure(meshtype,
    	[Nx, Ny, Nz],
    	CellSize(cellsize(facelocationX), cellsize(facelocationY), cellsize(facelocationZ)),
    	CellLocation(cellloc(facelocationX), cellloc(facelocationY), cellloc(facelocationZ)),
    	FaceLocation(facelocationX, facelocationY, facelocationZ),
    	G[[1,end],[1,end],[1,end]][:],
    	[
            G[[1, end], [1, end], 2:Nz+1][:];
    	    G[[1, end], 2:Ny+1, [1, end]][:];
        	G[2:Nx+1, [1, end], [1, end]][:]
        ])
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
function createMesh3D(Nx::Int, Ny::Int, Nz::Int, Width::Real, Height::Real, Depth::Real)
    dx = Width/Nx
    dy = Height/Ny
    dz = Depth/Nz
    meshStructure3D(mesh3D, dx, dy, dz, Nx, Ny, Nz)
end

"""
Creates a mesh on a 3D Cartesian domain with specified faces
locations

m = createMesh2D(faceLocationX, faceLocationY, faceLocationZ)

Inputs:

   + faceLocationX (Vector of reals): face locations along horizontal direction
   + faceLocationY (Vector of reals): face locations along vertical direction
   + faceLocationZ (Vector of reals): face locations along depth direction

Outputs:

   + m: a mesh structure

Usage:

```julia
fX = [0.0; 0.05; 0.1:0.1:1.0;]
fY = [0.0; 0.05; 0.12; 0.2:0.2:3.0;]
fz = [0.0:0.1:0.5;]

m=createMesh3D(fX, fY, fZ)
```
"""
function createMesh3D{T<:Real}(facelocationX::Array{T,1}, facelocationY::Array{T,1}, facelocationZ::Array{T,1})
    meshStructure3D(mesh3D, facelocationX, facelocationY, facelocationZ)
end

# %% ================== 3D CYLINDRICAL MESH ===================================
"""
Creates a uniform mesh on a 3D cylindrical domain

m = createMeshCylindrical3D(Nr, Ntheta, Nz, Radius, Angle, Height)

Inputs:

   + Nr (int): number of cells along radial direction
   + Ntheta (int): number of cells along angular direction
   + Nz (int): number of cells in vertical direction
   + Radius, Angle, Height: upper boundary of the domain, notice: Angle ≤ 2π

Outputs:

   + m: a mesh structure

Usage:

```julia
m=createMeshCylindrical3D(20, 30, 40, 1.0, π, 3.2)
```
"""
function createMeshCylindrical3D(Nr::Int, Ntheta::Int, Nz::Int, Radius::Real, Angle::Real, Height::Real)
    if Angle>2*pi
    	Angle = 2*pi
    	warn("The domain size adjusted to match a maximum of 2π")
    end
    dr = Radius/Nr
    dtheta = Angle/Ntheta
    dz = Height/Nz
    meshStructure3D(mesh3DCylindrical, dr, dtheta, dz, Nr, Ntheta, Nz)
end

"""
Creates a mesh on a 3D cylindrical domain with specified face locations

```
m = createMeshCylindrical3D(facelocationR, facelocationTheta, facelocationZ)
```

Inputs:

   + facelocation*: location of faces in each direction: radial, angular, vertical

Outputs:

   + m: a mesh structure

Usage:

```julia
fR = [0.0:0.02:0.08;0.1:0.1:1.0;]
fΘ = [0.0; 0.1; 0.2:0.2:1.0;]
fZ = [0.0:0.1:3.1;]
m=createMeshCylindrical3D(fR, fΘ, fZ)
```
"""
function createMeshCylindrical3D{T<:Real}(facelocationR::Array{T,1}, facelocationTheta::Array{T,1}, facelocationZ::Array{T,1})
    if facelocationTheta[end]>2*pi
    	facelocationTheta = facelocationTheta/facelocationTheta[end]*2.0*pi
    	warn("The domain size adjusted to match a maximum of 2π.")
    end
    meshStructure3D(mesh3DCylindrical, facelocationR, facelocationTheta, facelocationZ)
end

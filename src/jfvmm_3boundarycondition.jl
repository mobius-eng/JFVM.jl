# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 28 December, 2014
# ===============================

# ==========================================================
# Last changes:
#  * 2014-12-28: - support for cylindrical 3D and radial 2D
#  * 2014-12-30: - cell boundary functions added;
#                - Debugging and testing
# ==========================================================


# %% =================== CREATE BOUNDARY CONDITION =======================
"""
Creates boundary condition structure scaffolding for a mesh structure

```julia
bc = createBC(mesh)
```
"""
createBC(m::MeshStructure) = createBC(m.meshtype, m)

createBC(::Mesh1D, m) = create1DBC(m)
createBC(::Mesh1DPolar, m) = create1DBC(m)
createBC(::Mesh2D, m) = create2DBC(m)
createBC(::Mesh2DPolar, m) = create2DBC(m)
createBC(::Mesh2DCylindrical, m) = create2DBC(m)
createBC(::Mesh3D, m) = create3DBC(m)
createBC(::Mesh3DCylindrical, m) = create3DBC(m)

# All BC are non-periodic by default, hence `false`
function create1DBC{T<:Real}(m :: MeshStructure{T})
    onev = one(T)
    zerov = zero(T)
    BoundaryCondition(m,
        # Active: left and right
        BorderValue([onev], [zerov], [zerov], false),
        BorderValue([onev], [one], [one], false),
        # Nonactive
        BorderValue([zerov], [zerov], [zerov], false),
        BorderValue([zerov], [zerov], [zerov], false),
        BorderValue([zerov], [zerov], [zerov], false),
        BorderValue([zerov], [zerov], [zerov], false))
end

function create2DBC{T<:Real}(m :: MeshStructure{T})
    (Nx, Ny) = m.dims[1:2]
    BoundaryCondition(m,
        # Active: left, right, top, bottom
        BorderValue(ones(T,1,Ny), zeros(T,1,Ny), zeros(T,1,Ny), false),
        BorderValue(ones(T,1,Ny), zeros(T,1,Ny), zeros(T,1,Ny), false),
        BorderValue(ones(T,Nx,1), zeros(T,Nx,1), zeros(T,Nx,1), false),
        BorderValue(ones(T,Nx,1), zeros(T,Nx,1), zeros(T,Nx,1), false),
        # Inactive: front and back
        BorderValue([zero(T)], [zero(T)], [zero(T)], false),
        BorderValue([zero(T)], [zero(T)], [zero(T)], false))
end

function create3DBC{T<:Real}(m :: MeshStructure{T})
    (Nx, Ny, Nz) = m.dims
    BoundaryCondition(m,
        BorderValue(ones(T,1,Ny,Nz), zeros(T,1,Ny,Nz), zeros(T,1,Ny,Nz), false),
        BorderValue(ones(T,1,Ny,Nz), zeros(T,1,Ny,Nz), zeros(T,1,Ny,Nz), false),
        BorderValue(ones(T,Nx,1,Nz), zeros(T,Nx,1,Nz), zeros(T,Nx,1,Nz), false),
        BorderValue(ones(T,Nx,1,Nz), zeros(T,Nx,1,Nz), zeros(T,Nx,1,Nz), false),
        BorderValue(ones(T,Nx,Ny,1), zeros(T,Nx,Ny,1), zeros(T,Nx,Ny,1), false),
        BorderValue(ones(T,Nx,Ny,1), zeros(T,Nx,Ny,1), zeros(T,Nx,Ny,1), false))
end

"""
Creates the matrix of coefficients and RHS vector for a specified boundary
condition

```julia
(MBC, RHSBC) = boundaryConditionTerm(bc)
```
"""
boundaryConditionTerm(bc :: BoundaryCondition) =
    boundaryConditionTerm(bc.domain.meshType, bc)


boundaryConditionTerm(::Mesh1D, bc) = boundaryCondition1D(bc)
boundaryConditionTerm(::Mesh1DPolar, bc) = boundaryCondition1D(bc)
boundaryConditionTerm(::Mesh2D, bc) = boundaryCondition2D(bc)
boundaryConditionTerm(::Mesh2DCylindrical, bc) = boundaryCondition2D(bc)
boundaryConditionTerm(::Mesh2DPolar, bc) = boundaryConditionRadial2D(bc)
boundaryConditionTerm(::Mesh3D, bc) = boundaryCondition3D(bc)
boundaryConditionTerm(::Mesh3DCylindrical, bc) = boundaryConditionCylindrical3D(bc)

# %% ======================= BOUNDARY CARTESIAN 1D ===========================
"""
!!! INTERNAL !!!

Creates the matrix of coefficients and RHS vector for a specified boundary
condition in 1D domain (Cartesian and polar)
"""
function boundaryCondition1D{T<:Real}(BC::BoundaryCondition{T})
    Nx = BC.domain.dims[1]
    G = [1:Nx+2;]
    dx_1 = BC.domain.cellsize.x[1]
    dx_end = BC.domain.cellsize.x[end]
    # define the RHS column vector
    BCRHS = zeros(T,Nx+2)
    # initialize
    q = 0
    # Assign values to the boundary condition matrix and the RHS vector based
    # on the BC structure
    if !BC.right.periodic && !BC.left.periodic
        # non-periodic boundary condition
        # number of non-zero matrix entries, two on each side:
        nb = 4
        # define the vectors to be used for the creation of the sparse matrix
        # s - stores matrix elements in a vector
        # ii and jj provide the mapping from vector index q to matrix indices (i,j)
        # q -> (ii[q], jj[q])
        ii = zeros(Int64,nb)
        jj = zeros(Int64,nb)
        s = zeros(T, nb)
        # Right boundary
        i = Nx+2
        q=q+1
        ii[q] = G[i]
        jj[q] = G[i]
        # diagonal element
        s[q] = BC.right.b[1]/(2one(T)) + BC.right.a[1]/dx_end
        q=q+1
        ii[q] = G[i]
        jj[q] = G[i-1]
        # left-to diagonal
        s[q] = BC.right.b[1]/(2one(T)) - BC.right.a[1]/dx_end
        BCRHS[G[i]] = BC.right.c[1]
        # Left boundary
        i = 1
        q=q+1
        ii[q] = G[i]
        jj[q] = G[i+1]
        # diagonal
        s[q] = -(BC.left.b[1]/(2one(T)) + BC.left.a[1]/dx_1)
        q=q+1
        ii[q] = G[i]
        jj[q] = G[i]
        # right to diagonal
        s[q] = -(BC.left.b[1]/(2one(T)) - BC.left.a[1]/dx_1)
        BCRHS[G[i]] = -BC.left.c[1]
    elseif BC.right.periodic || BC.left.periodic
        # periodic boundary condition
        # why does it have 8 items?
        nb=8
        # define the vectors to be used for the creation of the sparse matrix
        # s - stores matrix elements in a vector
        # ii and jj provide the mapping from vector index q to matrix indices (i,j)
        # q -> (ii[q], jj[q])
        ii = zeros(Int64,nb)
        jj = zeros(Int64,nb)
        s = zeros(T, nb)
        # Right boundary
        i = Nx+2
        # diagonal
        q=q+1
        ii[q] = G[i]
        jj[q] = G[i]
        s[q] = one(T)
        # left to diagonal
        q=q+1
        ii[q] = G[i]
        jj[q] = G[i-1]
        s[q] = -one(T)
        # left-bottom corner of matrix
        q=q+1
        ii[q] = G[i]
        jj[q] = G[1]
        s[q] = dx_end/dx_1
        # left-bottom corner -> one right
        q=q+1
        ii[q] = G[i]
        jj[q] = G[2]
        s[q] = -dx_end/dx_1
        # RHS vector value
        BCRHS[G[i]] = zero(T)
        # Left boundary
        i = 1
        # top diagonal
        q=q+1
        ii[q] = G[1]
        jj[q] = G[1]
        s[q] = one(T)
        # top right to diagonal
        q=q+1
        ii[q] = G[1]
        jj[q] = G[2]
        s[q] = one(T)
        # top-right corner <- one left
        q=q+1
        ii[q] = G[1]
        jj[q] = G[Nx+1]
        s[q] = -one(T)
        # top-right corner
        q=q+1
        ii[q] = G[1]
        jj[q] = G[Nx+2]
        s[q] = -one(T);
        BCRHS[G[i]] = zero(T)
    end
    # Build the sparse matrix of the boundary conditions
    BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], Nx+2, Nx+2)
    (BCMatrix, BCRHS)
end

# %% ========================= BOUNDARY CARTESIAN 2D =========================
"""
!!! INTERNAL !!!

Creates the matrix of coefficients and RHS vector for a specified boundary
condition in 2D Cartesian and cylindrical domains
"""
function boundaryCondition2D{T<:Real}(BC::BoundaryCondition{T})
    Nx, Ny = BC.domain.dims[1:2]
    # correspondence of positions in domain and in system matrix
    G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
    dx_1 = BC.domain.cellsize.x[1]
    dx_end = BC.domain.cellsize.x[end]
    dy_1 = BC.domain.cellsize.y[1]
    dy_end = BC.domain.cellsize.y[end]
    # number of matrix entries:
    nb = 8*(Nx+Ny+2)
    # s stores the entries, ii and jj provide the mapping from vector entries
    # to matrix indices
    ii = zeros(Int64,nb)
    jj = zeros(Int64,nb)
    s = zeros(T, nb)
    # define the RHS column vector
    BCRHS = zeros(T,(Nx+2)*(Ny+2))
    # deal with phantom corner entries
    for q = 1:4
      ii[q] = BC.domain.corner[q]
      jj[q] = BC.domain.corner[q]
      s[q] = maximum(BC.top.b/2+BC.top.a/dy_end)
      BCRHS[BC.domain.corner[q]] = zero(T)
    end
    q = 4
    # top and bottom first
    if !BC.top.periodic && !BC.bottom.periodic
        # non-periodic boundary condition
        # top boundary
        j=Ny+2
        # iterate over all non-phantom indices
        for i=2:Nx+1
            q+=1
            # diagonal
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = BC.top.b[i-1]/2 + BC.top.a[i-1]/dy_end
            q+=1
            # next to diagonal, in 2D case it is not super clear where it will be
            ii[q] = G[i,j]
            jj[q] = G[i,j-1]
            s[q] = BC.top.b[i-1]/2 - BC.top.a[i-1]/dy_end
            # BC are only provided for real cells, hence need to -1 on c
            BCRHS[G[i,j]] = BC.top.c[i-1]
        end
        # bottom boundary
        j=1
        for i=2:Nx+1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j+1]
            s[q] = -(BC.bottom.b[i-1]/2 + BC.bottom.a[i-1]/dy_1)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = -(BC.bottom.b[i-1]/2 - BC.bottom.a[i-1]/dy_1)
            BCRHS[G[i,j]] = -(BC.bottom.c[i-1])
            end
    elseif BC.top.periodic || BC.bottom.periodic
        # periodic boundary condition
        # top boundary
        j=Ny+2
        for i=2:Nx+1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j-1]
            s[q] = -one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,1]
            s[q] = dy_end/dy_1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,2]
            s[q] = -dy_end/dy_1
            BCRHS[G[i,j]] = zero(T)
        end
        # bottom boundary
        j=1
            for i=2:Nx+1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j+1]
            s[q] = one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,Ny+1]
            s[q] = -one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,Ny+2]
            s[q] = -one(T)
            BCRHS[G[i,j]] = zero(T)
        end
    end
    # left and right
    if !BC.right.periodic && !BC.left.periodic
        # non-periodic boundary condition
        # Right boundary
        i=Nx+2
        for j=2:Ny+1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = BC.right.b[j-1]/2 + BC.right.a[j-1]/dx_end
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i-1,j]
            s[q] = BC.right.b[j-1]/2 - BC.right.a[j-1]/dx_end
            BCRHS[G[i,j]] = BC.right.c[j-1]
        end
        # Left boundary
        i = 1
        for j=2:Ny+1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i+1,j]
            s[q] = -(BC.left.b[j-1]/2 + BC.left.a[j-1]/dx_1)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = -(BC.left.b[j-1]/2 - BC.left.a[j-1]/dx_1)
            BCRHS[G[i,j]] = -(BC.left.c[j-1])
        end
    elseif BC.right.periodic || BC.left.periodic
        # periodic boundary condition
        # Right boundary
        i=Nx+2
        for j=2:Ny+1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i-1,j]
            s[q] = -one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[1,j]
            s[q] = dx_end/dx_1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[2,j]
            s[q] = -dx_end/dx_1
            BCRHS[G[i,j]] = zero(T)
        end
        # Left boundary
        i = 1;
        for j=2:Ny+1
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i,j]
            s[q] = one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[i+1,j]
            s[q] = one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[Nx+1,j]
            s[q] = -one(T)
            q+=1
            ii[q] = G[i,j]
            jj[q] = G[Nx+2,j]
            s[q] = -one(T)
            BCRHS[G[i,j]] = zero(T)
        end
    end
    # Build the sparse matrix of the boundary conditions
    BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
    (BCMatrix, BCRHS)
end

# TODO: continue changing types
# %% ========================= BOUNDARY Radial 2D ==========================
"""
!!! INTERNAL !!!

Creates the matrix of coefficients and RHS vector for a specified boundary
condition in 2D polar domain
"""
function boundaryConditionRadial2D(BC::BoundaryCondition)
    Nx, Ntheta = tuple(BC.domain.dims...)
    G=reshape([1:(Nx+2)*(Ntheta+2);], Nx+2, Ntheta+2)
    dx_1 = BC.domain.cellsize.x[1]
    dx_end = BC.domain.cellsize.x[end]
    dtheta_1 = BC.domain.cellsize.y[1]
    dtheta_end = BC.domain.cellsize.y[end]
    rp = BC.domain.cellcenters.x
    # number of boundary nodes:
    nb = 8*(Nx+Ntheta+2)

    # define the vectors to be used for the creation of the sparse matrix
    ii = zeros(Int64,nb)
    jj = zeros(Int64,nb)
    s = zeros(Float64, nb) # Float64 by default, but specify the type just in case

    # define the RHS column vector
    BCRHS = zeros((Nx+2)*(Ntheta+2))

    for q = 1:4
      ii[q] = BC.domain.corner[q]
      jj[q] = BC.domain.corner[q]
      s[q] = maximum(BC.top.b/2.0+BC.top.a./(dtheta_end*rp))
      BCRHS[BC.domain.corner[q]] = 0.0
    end
    q = 4
    # Assign values to the boundary condition matrix and the RHS vector based
    # on the BC structure
    if !BC.top.periodic && !BC.bottom.periodic # non-periodic boundary condition
      # top boundary
      j=Ntheta+2
      for i=2:Nx+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = BC.top.b[i-1]/2.0 + BC.top.a[i-1]/(dtheta_end*rp[i-1])
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j-1]
        s[q] = BC.top.b[i-1]/2.0 - BC.top.a[i-1]/(dtheta_end*rp[i-1])
        BCRHS[G[i,j]] = BC.top.c[i-1]
      end
      # bottom boundary
      j=1
      for i=2:Nx+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j+1]
        s[q] = -(BC.bottom.b[i-1]/2.0 + BC.bottom.a[i-1]/(dtheta_1*rp[i-1]))
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = -(BC.bottom.b[i-1]/2.0 - BC.bottom.a[i-1]/(dtheta_1*rp[i-1]))
        BCRHS[G[i,j]] = -(BC.bottom.c[i-1])
      end
    elseif BC.top.periodic || BC.bottom.periodic  # periodic boundary condition
      # top boundary
      j=Ntheta+2
      for i=2:Nx+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = 1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j-1]
        s[q] = -1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,1]
        s[q] = dtheta_end/dtheta_1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,2]
        s[q] = -dtheta_end/dtheta_1
        BCRHS[G[i,j]] = 0.0
      end
      # bottom boundary
      j=1
      for i=2:Nx+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = 1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j+1]
        s[q] = 1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,Ny+1]
        s[q] = -1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,Ny+2]
        s[q] = -1.0
        BCRHS[G[i,j]] = 0.0
      end
    end

    if !BC.right.periodic && !BC.left.periodic # non-periodic boundary condition
      # Right boundary
      i=Nx+2
      for j=2:Ntheta+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = BC.right.b[j-1]/2.0 + BC.right.a[j-1]/dx_end
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i-1,j]
        s[q] = BC.right.b[j-1]/2.0 - BC.right.a[j-1]/dx_end
        BCRHS[G[i,j]] = BC.right.c[j-1]
      end
      # Left boundary
      i = 1
      for j=2:Ntheta+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i+1,j]
        s[q] = -(BC.left.b[j-1]/2.0 + BC.left.a[j-1]/dx_1)
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = -(BC.left.b[j-1]/2.0 - BC.left.a[j-1]/dx_1)
        BCRHS[G[i,j]] = -(BC.left.c[j-1])
      end
    elseif BC.right.periodic || BC.left.periodic  # periodic boundary condition
      # Right boundary
      i=Nx+2
      for j=2:Ntheta+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = 1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i-1,j]
        s[q] = -1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[1,j]
        s[q] = dx_end/dx_1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[2,j]
        s[q] = -dx_end/dx_1
        BCRHS[G[i,j]] = 0.0
      end
      # Left boundary
      i = 1;
      for j=2:Ntheta+1
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i,j]
        s[q] = 1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[i+1,j]
        s[q] = 1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[Nx+1,j]
        s[q] = -1.0
        q+=1
        ii[q] = G[i,j]
        jj[q] = G[Nx+2,j]
        s[q] = -1.0
        BCRHS[G[i,j]] = 0.0
      end
    end

    # Build the sparse matrix of the boundary conditions
    BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], (Nx+2)*(Ntheta+2), (Nx+2)*(Ntheta+2))
    (BCMatrix, BCRHS)
end

# %% ================================ BOUNDARY 3D =============================
"""
!!! INTERNAL !!!

Creates the matrix of coefficients and RHS vector for a specified boundary
condition in 3D Cartesian domain
"""
function boundaryCondition3D(BC::BoundaryCondition)
    Nx, Ny, Nz = tuple(BC.domain.dims...)
    G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
    dx_1 = BC.domain.cellsize.x[1]
    dx_end = BC.domain.cellsize.x[end]
    dy_1 = BC.domain.cellsize.y[1]
    dy_end = BC.domain.cellsize.y[end]
    dz_1 = BC.domain.cellsize.z[1]
    dz_end = BC.domain.cellsize.z[end]

    # number of boundary nodes (exact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
    nb = 8*((Nx+1)*(Ny+1)+(Nx+1)*(Nz+1)+(Ny+1)*(Nz+1))

    # define the vectors to be used for the creation of the sparse matrix
    ii = zeros(Int64, nb)
    jj = zeros(Int64, nb)
    s = zeros(Float64, nb)

    # define the RHS column vector
    BCRHS = zeros(Float64, (Nx+2)*(Ny+2)*(Nz+2))

    # assign value to the corner nodes (useless cells)
    q = 1:8
    ii[q] = BC.domain.corner
    jj[q] = BC.domain.corner
    s[q] = 1.0
    BCRHS[BC.domain.corner] = 0.0

    # assign values to the edges (useless cells)
    q = q[end]+[1:length(BC.domain.edge);]
    ii[q] = BC.domain.edge
    jj[q] = BC.domain.edge
    s[q] = 1.0
    BCRHS[BC.domain.edge] = 0.0

    # Assign values to the boundary condition matrix and the RHS vector based
    # on the BC structure
    if !BC.top.periodic && !BC.bottom.periodic
        # top boundary
        j=Ny+2
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = BC.top.b/2.0 + BC.top.a/dy_end
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j-1,k]
        s[q] = BC.top.b/2.0 - BC.top.a/dy_end
        BCRHS[G[i,j,k]] = BC.top.c

        # Bottom boundary
        j=1
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j+1,k]
        s[q] = -(BC.bottom.b/2.0 + BC.bottom.a/dy_1) # consider the reverse direction of normal
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = -(BC.bottom.b/2.0 - BC.bottom.a/dy_1) # consider the reverse direction of normal
        BCRHS[G[i,j,k]] = -(BC.bottom.c)
    elseif BC.top.periodic || BC.bottom.periodic # periodic
        # top boundary
        j=Ny+2
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j-1,k]
        s[q] = -1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,1,k]
        s[q] = dy_end/dy_1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,2,k]
        s[q] = -dy_end/dy_1
        BCRHS[G[i,j,k]] = 0.0

        # Bottom boundary
        j=1
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j+1,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,Ny+1,k]
        s[q] = -1.0
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,Ny+2,k]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0
    end

    if !BC.right.periodic && !BC.left.periodic
        # Right boundary
        i=Nx+2
        j=2:Ny+1
        k=2:Nz+1
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = BC.right.b/2.0 + BC.right.a/dx_end
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i-1,j,k]
        s[q] = BC.right.b/2.0 - BC.right.a/dx_end
        BCRHS[G[i,j,k]] = BC.right.c

        # Left boundary
        i = 1
        j=2:Ny+1
        k=2:Nz+1
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i+1,j,k]
        s[q] = -(BC.left.b/2.0 + BC.left.a/dx_1)
        # consider the reverse direction of normal
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = -(BC.left.b/2.0 - BC.left.a/dx_1)  # consider the reverse direction of normal
        BCRHS[G[i,j,k]] = -(BC.left.c)
    elseif BC.right.periodic || BC.left.periodic # periodic
        # Right boundary
        i=Nx+2
        j=2:Ny+1
        k=2:Nz+1
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i-1,j,k]
        s[q] = -1.0
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[1,j,k]
        s[q] = dx_end/dx_1
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[2,j,k]
        s[q] = -dx_end/dx_1
        BCRHS[G[i,j,k]] = 0.0

        # Left boundary
        i = 1
        j=2:Ny+1
        k=2:Nz+1
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i+1,j,k]
        s[q] = 1.0
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[Nx+1,j,k]
        s[q] = -1.0
        q = q[end]+[1:Ny*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[Nx+2,j,k]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0
    end

    if !BC.front.periodic && !BC.back.periodic
        # Back boundary
        k=1
        i = 2:Nx+1
        j=2:Ny+1
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k+1]
        s[q] = -(BC.back.b/2.0 + BC.back.a/dz_1)  # consider the reverse direction of normal
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = -(BC.back.b/2.0 - BC.back.a/dz_1)  # consider the reverse direction of normal
        BCRHS[G[i,j,k]] = -(BC.back.c)

        # Front boundary
        k=Nz+2
        i = 2:Nx+1
        j=2:Ny+1
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = BC.front.b/2.0 + BC.front.a/dz_end
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k-1]
        s[q] = BC.front.b/2.0 - BC.front.a/dz_end
        BCRHS[G[i,j,k]] = BC.front.c
    elseif BC.front.periodic || BC.back.periodic  # periodic
        # Back boundary
        k=1
        i = 2:Nx+1
        j=2:Ny+1
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k+1]
        s[q] = 1.0
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,Nz+1]
        s[q] = -1.0
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,Nz+2]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0

        # Front boundary
        k=Nz+2
        i = 2:Nx+1
        j=2:Ny+1
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k-1]
        s[q] = -1.0
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,1]
        s[q] = dz_end/dz_1
        q = q[end]+[1:Nx*Ny;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,2]
        s[q] = -dz_end/dz_1
        BCRHS[G[i,j,k]] = 0.0
    end

    # Build the sparse matrix of the boundary conditions
    BCMatrix = sparse(ii[1:q[end]], jj[1:q[end]], s[1:q[end]],
        (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))

    (BCMatrix, BCRHS)
end

# %% ================== BOUNDARY Cylindrical 3D ==============================
"""
!!! INTERNAL !!!

Creates the matrix of coefficients and RHS vector for a specified boundary
condition in 2D cylindrical domain
"""
function boundaryConditionCylindrical3D(BC::BoundaryCondition)
    Nx, Ntheta, Nz = tuple(BC.domain.dims...)
    G=reshape([1:(Nx+2)*(Ntheta+2)*(Nz+2);], Nx+2, Ntheta+2, Nz+2)
    dx_1 = BC.domain.cellsize.x[1]
    dx_end = BC.domain.cellsize.x[end]
    dtheta_1 = BC.domain.cellsize.y[1]
    dtheta_end = BC.domain.cellsize.y[end]
    dz_1 = BC.domain.cellsize.z[1]
    dz_end = BC.domain.cellsize.z[end]

    rp = BC.domain.cellcenters.x
    #rp = repmat(BC.domain.cellcenters.x', 1, Nz)

    # number of boundary nodes (exact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
    nb = 8*((Nx+1)*(Ntheta+1)+(Nx+1)*(Nz+1)+(Ntheta+1)*(Nz+1))

    # define the vectors to be used for the creation of the sparse matrix
    ii = zeros(Int64, nb)
    jj = zeros(Int64, nb)
    s = zeros(Float64, nb)

    # define the RHS column vector
    BCRHS = zeros(Float64, (Nx+2)*(Ntheta+2)*(Nz+2))

    # assign value to the corner nodes (useless cells)
    q = 1:8
    ii[q] = BC.domain.corner
    jj[q] = BC.domain.corner
    s[q] = 1.0
    BCRHS[BC.domain.corner] = 0.0

    # assign values to the edges (useless cells)
    q = q[end]+[1:length(BC.domain.edge);]
    ii[q] = BC.domain.edge
    jj[q] = BC.domain.edge
    s[q] = 1.0
    BCRHS[BC.domain.edge] = 0.0

    # Assign values to the boundary condition matrix and the RHS vector based
    # on the BC structure
    if !BC.top.periodic && !BC.bottom.periodic
        # top boundary
        j=Ntheta+2
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = BC.top.b/2.0 + BC.top.a./(dtheta_end*rp)
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j-1,k]
        s[q] = BC.top.b/2.0 - BC.top.a./(dtheta_end*rp)
        BCRHS[G[i,j,k]] = BC.top.c

        # Bottom boundary
        j=1
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j+1,k]
        s[q] = -(BC.bottom.b/2.0 + BC.bottom.a./(dtheta_1*rp)) # consider the reverse direction of normal
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = -(BC.bottom.b/2.0 - BC.bottom.a./(dtheta_1*rp)) # consider the reverse direction of normal
        BCRHS[G[i,j,k]] = -(BC.bottom.c)
    elseif BC.top.periodic || BC.bottom.periodic # periodic
        # top boundary
        j=Ntheta+2
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,2,k]
        s[q] = -1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,2,k]
        s[q] = -1
        BCRHS[G[i,j,k]] = 0.0

        # Bottom boundary
        j=1
        i=2:Nx+1
        k=2:Nz+1
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,Ntheta+1,k]
        s[q] = -1.0
        q = q[end]+[1:Nx*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,Ntheta+1,k]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0
    end

    if !BC.right.periodic && !BC.left.periodic
        # Right boundary
        i=Nx+2
        j=2:Ntheta+1
        k=2:Nz+1
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = BC.right.b/2.0 + BC.right.a/dx_end
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i-1,j,k]
        s[q] = BC.right.b/2.0 - BC.right.a/dx_end
        BCRHS[G[i,j,k]] = BC.right.c

        # Left boundary
        i = 1
        j=2:Ntheta+1
        k=2:Nz+1
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i+1,j,k]
        s[q] = -(BC.left.b/2.0 + BC.left.a/dx_1)
        # consider the reverse direction of normal
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = -(BC.left.b/2.0 - BC.left.a/dx_1)  # consider the reverse direction of normal
        BCRHS[G[i,j,k]] = -(BC.left.c)
    elseif BC.right.periodic || BC.left.periodic # periodic
        # Right boundary
        i=Nx+2
        j=2:Ntheta+1
        k=2:Nz+1
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[2,j,k]
        s[q] = -1.0
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[2,j,k]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0

        # Left boundary
        i = 1
        j=2:Ntheta+1
        k=2:Nz+1
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[Nx+1,j,k]
        s[q] = -1.0
        q = q[end]+[1:Ntheta*Nz;]
        ii[q] = G[i,j,k]
        jj[q] = G[Nx+1,j,k]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0
    end

    if !BC.front.periodic && !BC.back.periodic
        # Back boundary
        k=1
        i = 2:Nx+1
        j=2:Ntheta+1
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k+1]
        s[q] = -(BC.back.b/2.0 + BC.back.a/dz_1)  # consider the reverse direction of normal
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = -(BC.back.b/2.0 - BC.back.a/dz_1)  # consider the reverse direction of normal
        BCRHS[G[i,j,k]] = -(BC.back.c)

        # Front boundary
        k=Nz+2
        i = 2:Nx+1
        j=2:Ntheta+1
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = BC.front.b/2.0 + BC.front.a/dz_end
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k-1]
        s[q] = BC.front.b/2.0 - BC.front.a/dz_end
        BCRHS[G[i,j,k]] = BC.front.c
    elseif BC.front.periodic || BC.back.periodic  # periodic
        # Back boundary
        k=1
        i = 2:Nx+1
        j=2:Ntheta+1
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,Nz+1]
        s[q] = -1.0
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,Nz+1]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0

        # Front boundary
        k=Nz+2
        i = 2:Nx+1
        j=2:Ntheta+1
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,k]
        s[q] = 1.0
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,2]
        s[q] = -1.0
        q = q[end]+[1:Nx*Ntheta;]
        ii[q] = G[i,j,k]
        jj[q] = G[i,j,2]
        s[q] = -1.0
        BCRHS[G[i,j,k]] = 0.0
    end

    # Build the sparse matrix of the boundary conditions
    BCMatrix = sparse(ii[1:q[end]], jj[1:q[end]], s[1:q[end]],
        (Nx+2)*(Ntheta+2)*(Nz+2), (Nx+2)*(Ntheta+2)*(Nz+2))

    (BCMatrix, BCRHS)
end

# %% ======================== Cell Boundary ==================================
"""
!!! INTERNAL !!!

Assigns values to ghost cells
"""
cellBoundary!{T<:Real}(phi :: CellValue{T}, bc :: BoundaryCondition{T}) =
    cellBoundary!(bc.domain.meshtype, phi, bc)

cellBoundary!(::Mesh1D, phi, bc) = cellBoundary1D!(phi, bc)
cellBoundary!(::Mesh1DPolar, phi, bc) = cellBoundary1D!(phi, bc)
cellBoundary!(::Mesh2D, phi, bc) = cellBoundary2D!(phi, bc)
cellBoundary!(::Mesh2DCylindrical, phi, bc) = cellBoundary2D!(phi, bc)
cellBoundary!(::Mesh2DPolar, phi, bc) = cellBoundaryRadial2D!(phi, bc)
cellBoundary!(::Mesh3D, phi, bc) = cellBoundary3D!(phi, bc)
cellBoundary!(::Mesh3DCylindrical, phi, bc) = cellBoundaryCylindrical3D!(phi, bc)

# =========================== 1D Cartesian ================================
function cellBoundary1D!(phi::CellValue, BC::BoundaryCondition)
    dx_1 = phi.domain.cellsize.x[1]
    dx_end = phi.domain.cellsize.x[end]
    # boundary condition (a d\phi/dx + b \phi = c, a column vector of [d a])
    # a (phi(i)-phi(i-1))/dx + b (phi(i)+phi(i-1))/2 = c
    # phi(i) (a/dx+b/2) + phi(i-1) (-a/dx+b/2) = c
    # Right boundary, i=m+2
    # phi(i) (a/dx+b/2) = c- phi(i-1) (-a/dx+b/2)
    # Left boundary, i=2
    #  phi(i-1) (-a/dx+b/2) = c - phi(i) (a/dx+b/2)
    # define the new phi
    if !BC.left.periodic && !BC.right.periodic
        phiBC = [
            (BC.left.c[1]-phi.value[2]*(BC.left.a[1]/dx_1+BC.left.b[1]/2.0))/
                (-BC.left.a[1]/dx_1+BC.left.b[1]/2.0);
            phi.value[2:end-1];
            (BC.right.c[1]-phi.value[end-1]*(-BC.right.a[1]/dx_end+BC.right.b[1]/2.0))/
                (BC.right.a[1]/dx_end+BC.right.b[1]/2.0)
        ]
    else
        phiBC = [phi.value[end]; phi.value[2:end-1]; phi.value[1]]
    end
    phi.value[:] = phiBC[:]
    phi
end

# =========================== 2D Cartesian ================================
function cellBoundary2D!{T<:Real}(phi::CellValue{T}, BC::BoundaryCondition{T})
    Nx, Ny = tuple(BC.domain.dims...)
    G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
    dx_1 = BC.domain.cellsize.x[1]
    dx_end = BC.domain.cellsize.x[end]
    dy_1 = BC.domain.cellsize.y[1]
    dy_end = BC.domain.cellsize.y[end]
    # define the output matrix
    phiBC = zeros(T,Nx+2, Ny+2)
    phiBC[:] = phi.value[:]
    # top and bottom
    if !BC.top.periodic && !BC.bottom.periodic
        # top boundary
        j=Ny+2
        for i = 2:Nx+1
            phiBC[i,j]=
                (BC.top.c[i-1]-phi.value[i,end-1].*(-BC.top.a[i-1]/dy_end+BC.top.b[i-1]/2.0))/
                    (BC.top.a[i-1]/dy_end+BC.top.b[i-1]/2.0)
        end
        # Bottom boundary
        j=1;
        for i = 2:Nx+1
            phiBC[i,j]=
                (BC.bottom.c[i-1]-phi.value[i,2].*(BC.bottom.a[i-1]/dy_1+BC.bottom.b[i-1]/2.0))/
                    (-BC.bottom.a[i-1]/dy_1+BC.bottom.b[i-1]/2.0)
        end
    else
        # top boundary
        j=Ny+2
        for i = 2:Nx+1
          phiBC[i,j]= phi.value[i,2]
        end
        # Bottom boundary
        j=1
        for i = 2:Nx+1
          phiBC[i,j]= phi.value[i,end-1]
        end
    end
    # left and right
    if !BC.left.periodic && !BC.right.periodic
        # Right boundary
        i = Nx+2
        for j = 2:Ny+1
            phiBC[i,j]=
                (BC.right.c[j-1]-phi.value[end-1,j]*(-BC.right.a[j-1]/dx_end+BC.right.b[j-1]/2.0))/
                    (BC.right.a[j-1]/dx_end+BC.right.b[j-1]/2.0)
        end

        # Left boundary
        i = 1
        for j = 2:Ny+1
            phiBC[i,j]=
                (BC.left.c[j-1]-phi.value[2,j]*(BC.left.a[j-1]/dx_1+BC.left.b[j-1]/2.0))/
                    (-BC.left.a[j-1]/dx_1+BC.left.b[j-1]/2)
        end
    else
        # Right boundary
        i = Nx+2
        for j = 2:Ny+1
            phiBC[i,j]= phi.value[2,j]
        end
        # Left boundary
        i = 1
        for j = 2:Ny+1
            phiBC[i,j]= phi.value[end-1,j]
        end
    end
    phi.value[:] = phiBC[:]
    phi
end

# ======================================== Cell Boundary 3D ===================================
function cellBoundary3D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
Nx, Ny, Nz = tuple(BC.domain.dims...)
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
dy_1 = BC.domain.cellsize.y[1]
dy_end = BC.domain.cellsize.y[end]
dz_1 = BC.domain.cellsize.z[1]
dz_end = BC.domain.cellsize.z[end]

# define the output matrix
phiBC = zeros(Nx+2, Ny+2, Nz+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ny+2
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.top.c-(phi.value[i,end-1:end-1,k]).*(-BC.top.a/dy_end+BC.top.b/2.0))./(BC.top.a/dy_end+BC.top.b/2.0)

    # Bottom boundary
    j=1
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.bottom.c-(phi.value[i,2:2,k]).*(BC.bottom.a/dy_1+BC.bottom.b/2.0))./(-BC.bottom.a/dy_1+BC.bottom.b/2.0)
else
    # top boundary
    j=Ny+2
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,2,k]

    # Bottom boundary
    j=1
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,end-1,k]
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nx+2
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.right.c-(phi.value[end-1:end-1,j,k]).*(-BC.right.a/dx_end+BC.right.b/2.0))./(BC.right.a/dx_end+BC.right.b/2.0)

    # Left boundary
    i = 1;
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.left.c-(phi.value[2:2,j,k]).*(BC.left.a/dx_1+BC.left.b/2.0))./(-BC.left.a/dx_1+BC.left.b/2.0)
else
    # Right boundary
    i = Nx+2
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[2,j,k]

    # Left boundary
    i = 1
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[end-1,j,k]
end

if !BC.bottom.periodic && !BC.top.periodic
    # front boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = Nz+2
    phiBC[i,j,k]= (BC.front.c-(phi.value[i,j,end-1:end-1]).*(-BC.front.a/dz_end+BC.front.b/2.0))./(BC.front.a/dz_end+BC.front.b/2.0)

    # back boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = 1
    phiBC[i,j,k]= (BC.back.c-(phi.value[i,j,2:2]).*(BC.back.a/dz_1+BC.back.b/2.0))./(-BC.back.a/dz_1+BC.back.b/2.0)
else
    # front boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = Nz+2
    phiBC[i,j,k]= phi.value[i,j,2]

    # back boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = 1
    phiBC[i,j,k]= phi.value[j,j,end-1]
end
phi.value[:] = phiBC[:]
phi
end


# ===================================== Cell Boundary Radial 2D ===========================
function cellBoundaryRadial2D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
Nr, Ntheta = tuple(BC.domain.dims...)
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
dr_1 = BC.domain.cellsize.x[1]
dr_end = BC.domain.cellsize.x[end]
dtheta_1 = BC.domain.cellsize.y[1]
dtheta_end = BC.domain.cellsize.y[end]
rp = BC.domain.cellcenters.x

# define the output matrix
phiBC = zeros(Nr+2, Ntheta+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    for i = 2:Nr+1
      phiBC[i,j]= (BC.top.c[i-1]-phi.value[i,end-1].*(-BC.top.a[i-1]/(dtheta_end*rp[i-1])+BC.top.b[i-1]/2.0))/(BC.top.a[i-1]/(dtheta_end*rp[i-1])+BC.top.b[i-1]/2.0)
    end
    # Bottom boundary
    j=1;
    for i = 2:Nr+1
      phiBC[i,j]= (BC.bottom.c[i-1]-phi.value[i,2].*(BC.bottom.a[i-1]/(dtheta_1*rp[i-1])+BC.bottom.b[i-1]/2.0))/(-BC.bottom.a[i-1]/(dtheta_1*rp[i-1])+BC.bottom.b[i-1]/2.0)
    end
else
    # top boundary
    j=Ntheta+2
    for i = 2:Nr+1
      phiBC[i,j]= phi.value[i,2]
    end
    # Bottom boundary
    j=1
    for i = 2:Nr+1
      phiBC[i,j]= phi.value[i,end-1]
    end
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nr+2
    for j = 2:Ntheta+1
      phiBC[i,j]= (BC.right.c[j-1]-phi.value[end-1,j]*(-BC.right.a[j-1]/dr_end+BC.right.b[j-1]/2.0))/(BC.right.a[j-1]/dr_end+BC.right.b[j-1]/2.0)
    end

    # Left boundary
    i = 1
    for j = 2:Ntheta+1
      phiBC[i,j]= (BC.left.c[j-1]-phi.value[2,j]*(BC.left.a[j-1]/dr_1+BC.left.b[j-1]/2.0))/(-BC.left.a[j-1]/dr_1+BC.left.b[j-1]/2.0)
    end
else
    # Right boundary
    i = Nr+2
    for j = 2:Ntheta+1
      phiBC[i,j]= phi.value[2,j]
    end
    # Left boundary
    i = 1
    for j = 2:Ntheta+1
      phiBC[i,j]= phi.value[end-1,j]
    end
end
phi.value[:] = phiBC[:]
phi
end


# ===================================== Cell Boundary Cylindrical 3D ===========================
function cellBoundaryCylindrical3D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
Nr, Ntheta, Nz = tuple(BC.domain.dims...)
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
dr_1 = BC.domain.cellsize.x[1]
dr_end = BC.domain.cellsize.x[end]
dtheta_1 = BC.domain.cellsize.y[1]
dtheta_end = BC.domain.cellsize.y[end]
dz_1 = BC.domain.cellsize.z[1]
dz_end = BC.domain.cellsize.z[end]
#rp = zeros(Nx,1,Nz)
#rp[:,1,:] = repmat(phi.domain.cellcenters.x, 1, Nz)
rp = phi.domain.cellcenters.x

# define the output matrix
phiBC = zeros(Nr+2, Ntheta+2, Nz+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.top.c-(phi.value[i,end-1:end-1,k]).*(-BC.top.a./(dtheta_end*rp)+BC.top.b/2.0))./(BC.top.a./(dtheta_end*rp)+BC.top.b/2.0)

    # Bottom boundary
    j=1
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.bottom.c-(phi.value[i,2:2,k]).*(BC.bottom.a./(dtheta_1*rp)+BC.bottom.b/2.0))./(-BC.bottom.a./(dtheta_1*rp)+BC.bottom.b/2.0)
else
    # top boundary
    j=Ntheta+2
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,2,k]

    # Bottom boundary
    j=1
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,end-1,k]
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nr+2
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.right.c-(phi.value[end-1:end-1,j,k]).*(-BC.right.a/dr_end+BC.right.b/2.0))./(BC.right.a/dr_end+BC.right.b/2.0)

    # Left boundary
    i = 1;
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.left.c-(phi.value[2:2,j,k]).*(BC.left.a/dr_1+BC.left.b/2.0))./(-BC.left.a/dr_1+BC.left.b/2.0)
else
    # Right boundary
    i = Nr+2
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[2,j,k]

    # Left boundary
    i = 1
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[end-1,j,k]
end

if !BC.bottom.periodic && !BC.top.periodic
    # front boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = Nz+2
    phiBC[i,j,k]= (BC.front.c-(phi.value[i,j,end-1:end-1]).*(-BC.front.a/dz_end+BC.front.b/2.0))./(BC.front.a/dz_end+BC.front.b/2.0)

    # back boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = 1
    phiBC[i,j,k]= (BC.back.c-(phi.value[i,j,2:2]).*(BC.back.a/dz_1+BC.back.b/2.0))./(-BC.back.a/dz_1+BC.back.b/2.0)
else
    # front boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = Nz+2
    phiBC[i,j,k]= phi.value[i,j,2]

    # back boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = 1
    phiBC[i,j,k]= phi.value[j,j,end-1]
end
phi.value[:] = phiBC[:]
phi
end

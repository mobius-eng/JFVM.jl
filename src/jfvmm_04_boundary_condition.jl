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
createBC(m::MeshStructure) = createBC(Val{length(dimensions(m))}, m)

createBC(::Type{Val{2}}, m) = create2DBC(m)
createBC(::Type{Val{3}}, m) = create3DBC(m)

function create2DBC{T<:Real}(m :: MeshStructure{T})
    Nx, Ny = dimensions(m)
    BoundaryCondition(domain = m,
        # Active: left, right, top, bottom
        left = BorderValue(ones(Bool,Ny), zeros(T,Ny)),
        right = BorderValue(ones(Bool,Ny), zeros(T,Ny)),
        bottom = BorderValue(ones(Bool,Nx), zeros(T,Nx)),
        top = BorderValue(ones(Bool,Nx), zeros(T,Nx)),
        # Inactive: front and back
        back = BorderValue([zero(Bool)], [zero(T)]),
        front = BorderValue([zero(Bool)], [zero(T)]))
end

function create3DBC{T<:Real}(m :: MeshStructure{T})
    Nx, Ny, Nz = dimensions(m)
    BoundaryCondition(domain = m,
        left = BorderValue(ones(Bool,Ny,Nz), zeros(T,Ny,Nz)),
        right = BorderValue(ones(Bool,Ny,Nz), zeros(T,Ny,Nz)),
        bottom = BorderValue(ones(Bool,Nx,Nz), zeros(T,Nx,Nz)),
        top = BorderValue(ones(Bool,Nx,Nz), zeros(T,Nx,Nz)),
        back = BorderValue(ones(Bool,Nx,Ny), zeros(T,Nx,Ny)),
        front = BorderValue(ones(Bool,Nx,Ny), zeros(T,Nx,Ny)))
end

"""
Creates the matrix of coefficients and RHS vector for a specified boundary
condition

```julia
(MBC, RHSBC) = boundaryConditionTerm(bc, v)
```
"""
function boundaryConditionTerm(bc :: BoundaryCondition, v)
    if length(dimensions(bc.domain)) == 2
        boundaryCondition2D(bc, v)
    else
        boundaryCondition3D(bc, v)
    end
end
# %% ========================= BOUNDARY CARTESIAN 2D =========================
# Unfortunately have to keep both even though the logic is the same
function assignBC2Di!(mesh, bv, i, v, fluxV)
    println("assignBC2Dx! with i = $i")
    Nj = dimensions(mesh)[2]
    for j=1:Nj
        println("j = $j")
        println("value = $(bv.value[j])")
        println("area = $(mesh.faces.ifaces[i,j].area)")
        println("n = $(globalindex(mesh, i, j))")
        fluxV.ivalue[i,j] += bv.isflux[j]? bv.value[j] * mesh.faces.ifaces[i,j].area : v.ivalue[i,j] * mesh.faces.ifaces[i,j].area * bv.value[j]
    end
end

function assignBC2Dj!(mesh, bv, j, v, fluxV)
    Ni = dimensions(mesh)[1]
    for i=1:Ni
        fluxV.jvalue[i,j] += bv.isflux[i]? bv.value[i]*mesh.faces.jfaces[i,j].area : v.jvalue[i,j]*mesh.faces.jfaces[i,j].area*bv.value[i]
    end
end

"""
!!! INTERNAL !!!

Creates the matrix of coefficients and RHS vector for a specified boundary
condition in 2D Cartesian and cylindrical domains
"""
function boundaryCondition2D!{T<:Real}(bc::BoundaryCondition{T}, convvel :: FaceValue{T}, fluxV)
    Ni, Nj = dimensions(bc.domain)
    mesh = bc.domain
    # FluxV: non-linear part of the flux resulted from BC on element
    assignBC2Di!(mesh, bc.left, 1, convvel, fluxV)
    assignBC2Di!(mesh, bc.right, Ni+1, convvel, fluxV)
    assignBC2Dj!(mesh, bc.bottom, 1, convvel, fluxV)
    assignBC2Dj!(mesh, bc.top, Nj+1, convvel, fluxV)
    return nothing
end

function boundaryCondition2D{T<:Real}(bc::BoundaryCondition{T}, convvel :: FaceValue{T})
    Ni, Nj = dimensions(bc.domain)
    fluxV = createFaceVariable(bc.domain, zero(T))
    boundaryCondition2D!(bc, convvel, fluxV)
    return fluxV
end

# %% ================================ BOUNDARY 3D =============================
function assignBC3Di!(mesh, bv, i, v, fluxV)
    # i refers to face index
    Nj, Nk = dimensions(mesh)[2:3]
    for k = 1:Nk
        for j=1:Nj
            fluxV.ivalue[i,j,k] += bv.isflux[j,k]? bv.value[j,k]*mesh.faces.ifaces[i,j,k].area : v.ivalue[i,j,k]*mesh.faces.ifaces[i,j,k].area*bv.value[j,k]
        end
    end
end

function assignBC3Dj!(mesh, bv, j, v, fluxV)
    Ni, Nk = dimensions(mesh)[(1,3)]
    for k = 1:Nk
        for i=1:Ni
            fluxV.jvalue[i,j,k] += bv.isflux[i,k]? bv.value[i,k]*mesh.faces.jfaces[i,j,k].area : v.jvalue[i,j,k]*mesh.faces.jfaces[i,j,k].area*bv.value[i,k]
        end
    end
end

function assignBC3Dk!(mesh, bv, k, v, fluxV)
    Ni, Nj = dimensions(mesh)[1:2]
    for j = 1:Nj
        for i=1:Ni
            fluxV.kvalue[i,j,k] += bv.isflux[i,j]? bv.value[i,j]*mesh.faces.kfaces[i,j,k].area : v.kvalue[i,j,k]*mesh.faces.kfaces[i,j,k].area*bv.value[i,j]
        end
    end
end

"""
!!! INTERNAL !!!

Creates the matrix of coefficients and RHS vector for a specified boundary
condition in 2D Cartesian and cylindrical domains
"""
function boundaryCondition3D!{T<:Real}(bc::BoundaryCondition{T}, convvel :: FaceValue{T}, fluxV)
    Ni, Nj, Nk = dimensions(bc.domain)
    mesh = bc.domain
    # FluxV: non-linear part of the flux resulted from BC on element
    assignBC3Di!(mesh, bc.left, 1, convvel, fluxV)
    assignBC3Di!(mesh, bc.right, Ni+1, convvel, fluxV)
    assignBC3Dj!(mesh, bc.bottom, 1, convvel, fluxV)
    assignBC3Dj!(mesh, bc.top, Nj+1, convvel, fluxV)
    assignBC3Dk!(mesh, bv.back, 1, convvel, fluxV)
    assignBC3Dk!(mesh, bv.front, Nk+1, convvel, fluxV)
    return nothing
end

function boundaryCondition3D{T<:Real}(bc::BoundaryCondition{T}, convvel :: FaceValue{T})
    fluxV = createFaceVariable(bc.domain, 0.0)
    boundaryCondition3D!(bc, convvel, fluxV)
    return fluxV
end

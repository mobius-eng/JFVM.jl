# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

"""
Produces matrix and RHS vector for transient term (∂φ/∂t)

```julia
Mt, RHSt = transientTerm(φ0, Δt, α)
```

Input:

- `φ0` (cell variable) : previous value of variable
- `Δt` (real) : time step
- `α` (optional) : coefficient in front of ∂φ/∂t term. `α` defaults to 1.
  Can also be an array over the domain or a cell variable.

Output:

- `Mt` : transient term matrix
- `RHSt` : transient term RHS vector

"""
function transientTerm{T<:Real}(phi_old::CellValue{T}, dt::T, alfa::T)
    transientTerm(phi_old, dt, alfa .* ones(T,tuple(phi_old.domain.dims...)))
end

transientTerm{T<:Real}(phi_old::CellValue{T}, dt::T) = transientTerm(phi_old, dt, one(T))

function transientTerm{T<:Real}(phi_old::CellValue{T}, dt::T, alfa::CellValue{T})
    d = length(phi_old.domain.dims)
    if d==1
      transientTerm1D(phi_old, dt, alfa.value[2:end-1])
    elseif d==2
      transientTerm2D(phi_old, dt, alfa.value[2:end-1,2:end-1])
    elseif d==3
      transientTerm3D(phi_old, dt, alfa[2:end-1,2:end-1,2:end-1])
    end
end

function transientTerm{T<:Real}(phi_old::CellValue{T}, dt::T, alfa::Array{T})
    d = length(phi_old.domain.dims)
    if d==1
      transientTerm1D(phi_old, dt, alfa)
    elseif d==2
      transientTerm2D(phi_old, dt, alfa)
    elseif d==3
      transientTerm3D(phi_old, dt, alfa)
    end
end

# %% Internal implementations

function transientTerm1D{T<:Real}(phi_old::CellValue, dt::T, alfa::Array{T})
    # returns the matrix and RHS for a d(phi)/dt term
    # extract data from the mesh structure
    Nx = phi_old.domain.dims[1]
    G = [1:Nx+2;]
    # rearrange the matrix of k and build the sparse matrix for internal cells
    row_index = reshape(G[2:Nx+1],Nx) # main diagonal (only internal cells)
    AP_diag = reshape(alfa/dt,Nx)
    M = sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
    # define the RHS Vector
    RHS = zeros(T,Nx+2)
    # assign the values of the RHS vector
    RHS[row_index] = reshape(alfa.*phi_old.value[2:Nx+1]/dt,Nx)
    (M, RHS)
end

function transientTerm2D{T<:Real}(phi_old::CellValue, dt::T, alfa::Array{T})
    # returns the matrix and RHS for a d(phi)/dt term
    # extract data from the mesh structure
    Nx = phi_old.domain.dims[1]
    Ny = phi_old.domain.dims[2]
    G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
    # rearrange the matrix of k and build the sparse matrix for internal cells
    row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny) # main diagonal (only internal cells)
    AP_diag = reshape(alfa/dt,Nx*Ny)
    M = sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
    # define the RHS Vector
    RHS = zeros(T,(Nx+2)*(Ny+2))
    # assign the values of the RHS vector
    RHS[row_index] = reshape(alfa.*phi_old.value[2:Nx+1,2:Ny+1]/dt,Nx*Ny)
    (M, RHS)
end

function transientTerm3D{T<:Real}(phi_old::CellValue, dt::T, alfa::Array{T})
    # returns the matrix and RHS for a d(phi)/dt term

    # extract data from the mesh structure
    Nx = phi_old.domain.dims[1]
    Ny = phi_old.domain.dims[2]
    Nz = phi_old.domain.dims[3]
    G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)

    # rearrange the matrix of k and build the sparse matrix for internal cells
    row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal (only internal cells)
    AP_diag = reshape(alfa/dt,Nx*Ny*Nz)
    M = sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))

    # define the RHS Vector
    RHS = zeros(T,(Nx+2)*(Ny+2)*(Nz+2))

    # assign the values of the RHS vector
    RHS[row_index] = reshape(alfa.*phi_old.value[2:Nx+1,2:Ny+1,2:Nz+1]/dt,Nx*Ny*Nz)

    (M, RHS)

end

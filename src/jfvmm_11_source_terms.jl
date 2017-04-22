# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================================================
# Changes:
#    2015-01-10 changed numberofcells to dims
# ================================================================

# %% ======================= Linear source term ========================
function linearSourceTerm{T<:Real}(betta0::CellValue{T})
    m = betta0.domain
    d = length(m.dims)
    if d ==1
      Nx = m.dims[1]
      G = [1:Nx+2;]
      b = betta0.value[2:end-1]
      row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
      AP_diag = reshape(b,Nx)
      sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
    elseif d == 2
      Nx, Ny = m.dims
      G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
      b = betta0.value[2:end-1,2:end-1]
      row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
      AP_diag = reshape(b,Nx*Ny)
      sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
    elseif d == 3
      Nx, Ny, Nz = m.dims
      G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
      b = betta0.value[2:end-1,2:end-1,2:end-1]
      row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
      AP_diag = reshape(b,Nx*Ny*Nz)
      sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
    end
end

linearSourceTerm{T<:Real}(m::MeshStructure{T}, betta0::T) =
    linearSourceTerm(createCellVariable(m, betta0))

linearSourceTerm{T<:Real}(m::MeshStructure{T}, betta0::Array{T}) =
    linearSourceTerm(CellValue(m, betta0))

# %% ================================== constant source term ================================
function constantSourceTerm{T<:Real}(phi0::CellValue{T})
    m = phi0.domain
    d = length(m.dims)
    if  d ==1
      Nx = m.dims[1]
      G = [1:Nx+2;]
      row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
      RHS = zeros(T, Nx+2)
      RHS[row_index] = reshape(phi0.value[2:end-1],Nx)
    elseif d == 2
      Nx, Ny = m.dims
      G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
      row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
      RHS = zeros(T, (Nx+2)*(Ny+2))
      RHS[row_index] = reshape(phi0.value[2:end-1,2:end-1],Nx*Ny)
    elseif d == 3
      Nx, Ny, Nz = m.dims
      G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
      row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
      RHS = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
      RHS[row_index] = reshape(phi0.value[2:end-1,2:end-1,2:end-1],Nx*Ny*Nz)
    end
    RHS
end

constantSourceTerm{T<:Real}(m::MeshStructure{T}, phi0::T) =
    constantSourceTerm(createCellVariable(m, phi0))

constantSourceTerm{T<:Real}(m::MeshStructure{T}, phi0::Array{T}) =
    constantSourceTerm(CellValue(m, phi0))

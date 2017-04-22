# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================================================
# Changes:
#    2014-12-30 added 2D radial and 3D cylindrical grids
#    2014-12-31 many many bugs fixed
#    2015-01-10 extended to accept nonuniform grids
# ================================================================

"""
Calculates the RHS vector of the (known) vector field divergence

```julia
RHSv = divergenceTerm(f)
```

Input:
- `f` : face vector field

Output:
- `RHSv` : RHS vector for global equation
"""
divergenceTerm{T<:Real}(F::FaceValue{T}) = divergenceTerm(F.domain.meshtype, F)

divergenceTerm{T<:Real}(::Mesh1D, F::FaceValue{T}) = divergenceTerm1D(F)

divergenceTerm{T<:Real}(::Mesh1DPolar, F::FaceValue{T}) = divergenceTermCylindrical1D(F)

divergenceTerm{T<:Real}(::Mesh2D, F::FaceValue{T}) = divergenceTerm2D(F)[1]

divergenceTerm{T<:Real}(::Mesh2DCylindrical, F::FaceValue{T}) = divergenceTermCylindrical2D(F)[1]

divergenceTerm{T<:Real}(::Mesh2DPolar, F::FaceValue{T}) = divergenceTermRadial2D(F)[1]

divergenceTerm{T<:Real}(::Mesh3D, F::FaceValue{T}) = divergenceTerm3D(F)[1]

divergenceTerm{T<:Real}(::Mesh3DCylindrical, F::FaceValue{T}) = divergenceTermCylindrical3D(F)[1]


# %% =============== Divergence 1D Term ============================
function divergenceTerm1D{T<:Real}(F::FaceValue{T})
    Nx = F.domain.dims[1]
    G = [1:Nx+2;]
    DX = F.domain.cellsize.x[2:end-1]
    # define the vector of cell index
    row_index = reshape(G[2:Nx+1],Nx) # main diagonal
    # compute the divergence
    div_x = (F.xvalue[2:Nx+1]-F.xvalue[1:Nx])./DX
    # define the RHS Vector
    RHSdiv = zeros(T, Nx+2)
    # assign the values of the RHS vector
    RHSdiv[row_index] = reshape(div_x,Nx)
    return RHSdiv
end

# %% ============== Divergence Cylindrical 1D Term ========================
function divergenceTermCylindrical1D{T<:Real}(F::FaceValue{T})
    Nx = F.domain.dims[1]
    G = [1:Nx+2;]
    DX = F.domain.cellsize.x[2:end-1]
    rp = F.domain.cellcenters.x
    rf = F.domain.facecenters.x
    # define the vector of cell index
    row_index = reshape(G[2:Nx+1],Nx) # main diagonal
    # reassign the east, west, north, and south flux vectors for the
    # code readability
    Fe = F.xvalue[2:Nx+1]
    Fw = F.xvalue[1:Nx]
    re = rf[2:Nx+1]
    rw = rf[1:Nx]
    # compute the divergence
    div_x = (re.*Fe-rw.*Fw)./(DX.*rp)
    # define the RHS Vector
    RHSdiv = zeros(T, Nx+2)
    # assign the values of the RHS vector
    RHSdiv[row_index] = reshape(div_x,Nx)
    return RHSdiv
end

# %% =============== Divergence 2D Term ============================
function divergenceTerm2D{T<:Real}(F::FaceValue{T})
    Nx, Ny = F.domain.dims
    G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
    DX = F.domain.cellsize.x[2:end-1]
    DY = Array(T, 1, Ny)
    DY[:] = F.domain.cellsize.y[2:end-1]
    # define the vector of cell index
    row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny) # main diagonal
    # reassign the east, west, north, and south flux vectors for the
    # code readability
    Fe = F.xvalue[2:Nx+1,:]
    Fw = F.xvalue[1:Nx,:]
    Fn = F.yvalue[:,2:Ny+1]
    Fs = F.yvalue[:,1:Ny]
    # compute the divergence
    div_x = (Fe - Fw)./DX
    div_y = (Fn - Fs)./DY
    # define the RHS Vector
    RHSdiv = zeros(T, (Nx+2)*(Ny+2))
    RHSdivx = zeros(T, (Nx+2)*(Ny+2))
    RHSdivy = zeros(T, (Nx+2)*(Ny+2))
    # assign the values of the RHS vector
    RHSdiv[row_index] = reshape(div_x+div_y,Nx*Ny)
    RHSdivx[row_index] = reshape(div_x,Nx*Ny)
    RHSdivy[row_index] = reshape(div_y,Nx*Ny)
    (RHSdiv, RHSdivx, RHSdivy)
end

# %% ============== Divergence 2D Cylindrical Term ======================
function divergenceTermCylindrical2D{T<:Real}(F::FaceValue{T})
    Nr, Nz = F.domain.dims
    G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
    dr = F.domain.cellsize.x[2:end-1]
    dz= Array(T, 1, Nz)
    dz[:] = F.domain.cellsize.y[2:end-1]
    rp = F.domain.cellcenters.x
    rf = F.domain.facecenters.x
    # define the vector of cell index
    row_index = reshape(G[2:Nr+1,2:Nz+1],Nr*Nz) # main diagonal
    # reassign the east, west, north, and south flux vectors for the
    # code readability
    Fe = F.xvalue[2:Nr+1,:]
    Fw = F.xvalue[1:Nr,:]
    Fn = F.yvalue[:,2:Nz+1]
    Fs = F.yvalue[:,1:Nz]
    re = rf[2:Nr+1]
    rw = rf[1:Nr]
    # compute the divergence
    div_x = (re.*Fe - rw.*Fw)./(dr.*rp)
    div_y = (Fn - Fs)./dz
    # define the RHS Vector
    RHSdiv = zeros(T, (Nr+2)*(Nz+2))
    RHSdivx = zeros(T, (Nr+2)*(Nz+2))
    RHSdivy = zeros(T, (Nr+2)*(Nz+2))
    # assign the values of the RHS vector
    RHSdiv[row_index] = reshape(div_x+div_y,Nr*Nz)
    RHSdivx[row_index] = reshape(div_x,Nr*Nz)
    RHSdivy[row_index] = reshape(div_y,Nr*Nz)
    (RHSdiv, RHSdivx, RHSdivy)
end

# %% ========== Divergence 2D Radial Term ========================
function divergenceTermRadial2D{T<:Real}(F::FaceValue{T})
    Nr, Ntheta = F.domain.dims
    G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
    dr = F.domain.cellsize.x[2:end-1]
    dtheta= Array(T, 1, Ntheta)
    dtheta[:]= F.domain.cellsize.y[2:end-1]
    rp = F.domain.cellcenters.x
    rf = F.domain.facecenters.x
    # define the vector of cell index
    row_index = reshape(G[2:Nr+1,2:Ntheta+1],Nr*Ntheta) # main diagonal
    # reassign the east, west, north, and south flux vectors for the
    # code readability
    Fe = F.xvalue[2:Nr+1,:]
    Fw = F.xvalue[1:Nr,:]
    Fn = F.yvalue[:,2:Ntheta+1]
    Fs = F.yvalue[:,1:Ntheta]
    re = rf[2:Nr+1]
    rw = rf[1:Nr]
    # compute the divergence
    div_x = (re.*Fe-rw.*Fw)./(dr.*rp)
    div_y = (Fn-Fs)./(dtheta.*rp)
    # define the RHS Vector
    RHSdiv = zeros(T, (Nr+2)*(Ntheta+2))
    RHSdivx = zeros(T, (Nr+2)*(Ntheta+2))
    RHSdivy = zeros(T, (Nr+2)*(Ntheta+2))
    # assign the values of the RHS vector
    RHSdiv[row_index] = reshape(div_x+div_y,Nr*Ntheta)
    RHSdivx[row_index] = reshape(div_x,Nr*Ntheta)
    RHSdivy[row_index] = reshape(div_y,Nr*Ntheta)
    (RHSdiv, RHSdivx, RHSdivy)
end

# %% =============== Divergence 3D Term ============================
function divergenceTerm3D{T<:Real}(F::FaceValue{T})
    Nx, Ny, Nz = F.domain.dims
    G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
    dx = F.domain.cellsize.x[2:end-1]
    dy = Array(T, 1, Ny)
    dy[:] = F.domain.cellsize.y[2:end-1]
    dz = Array(T, 1,1,Nz)
    dz[:] = F.domain.cellsize.z[2:end-1]
    # define the vector of cell index
    row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal
    # reassign the east, west, north, and south flux vectors for the
    # code readability
    Fe = F.xvalue[2:Nx+1,:,:]
    Fw = F.xvalue[1:Nx,:,:]
    Fn = F.yvalue[:,2:Ny+1,:]
    Fs = F.yvalue[:,1:Ny,:]
    Ff = F.zvalue[:,:,2:Nz+1]
    Fb = F.zvalue[:,:,1:Nz]
    # compute the divergence
    div_x = (Fe - Fw)./dx
    div_y = (Fn - Fs)./dy
    div_z = (Ff - Fb)./dz
    # define the RHS Vector
    RHSdiv = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    RHSdivx = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    RHSdivy = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    RHSdivz = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    # assign the values of the RHS vector
    RHSdiv[row_index] = reshape(div_x+div_y+div_z,Nx*Ny*Nz)
    RHSdivx[row_index] = reshape(div_x,Nx*Ny*Nz)
    RHSdivy[row_index] = reshape(div_y,Nx*Ny*Nz)
    RHSdivz[row_index] = reshape(div_z,Nx*Ny*Nz)
    (RHSdiv, RHSdivx, RHSdivy, RHSdivz)
end

# =============== Divergence 3D Cylindrical Term ============================
function divergenceTermCylindrical3D(F::FaceValue)
    Nx, Ny, Nz = F.domain.dims
    G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
    dx = F.domain.cellsize.x[2:end-1]
    dy = Array(T, 1, Ny)
    dy[:] = F.domain.cellsize.y[2:end-1]
    dz = Array(T, 1,1,Nz)
    dz[:] = F.domain.cellsize.z[2:end-1]
    rp = F.domain.cellcenters.x
    # define the vector of cell index
    row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal
    # reassign the east, west, north, and south flux vectors for the
    # code readability
    Fe = F.xvalue[2:Nx+1,:,:]
    Fw = F.xvalue[1:Nx,:,:]
    Fn = F.yvalue[:,2:Ny+1,:]
    Fs = F.yvalue[:,1:Ny,:]
    Ff = F.zvalue[:,:,2:Nz+1]
    Fb = F.zvalue[:,:,1:Nz]
    # compute the divergence
    div_x = (Fe - Fw)./dx
    div_y = (Fn - Fs)./(dy.*rp)
    div_z = (Ff - Fb)./dz
    # define the RHS Vector
    RHSdiv = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    RHSdivx = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    RHSdivy = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    RHSdivz = zeros(T, (Nx+2)*(Ny+2)*(Nz+2))
    # assign the values of the RHS vector
    RHSdiv[row_index] = reshape(div_x+div_y+div_z,Nx*Ny*Nz)
    RHSdivx[row_index] = reshape(div_x,Nx*Ny*Nz)
    RHSdivy[row_index] = reshape(div_y,Nx*Ny*Nz)
    RHSdivz[row_index] = reshape(div_z,Nx*Ny*Nz)
    (RHSdiv, RHSdivx, RHSdivy, RHSdivz)
end

# %% ===================== Gradient Term ============================
"""
Calculates the gradient of the known cell field

```julia
dx = gradientTerm(x)
```

Input:

- `x` : cell field variable

Output:

- `dx` : face field gradient of `x`

"""
gradientTerm{T<:Real}(phi :: CellValue{T}) = gradientTerm(phi.domain.meshtype, phi)

gradientTerm(::Mesh1D, phi :: CellValue) = gradientTerm1D(phi)
gradientTerm(::Mesh1DPolar, phi :: CellValue) = gradientTerm1D(phi)
gradientTerm(::Mesh2D, phi :: CellValue) = gradientTerm2D(phi)
gradientTerm(::Mesh2DCylindrical, phi :: CellValue) = gradientTerm2D(phi)
gradientTerm(::Mesh2DPolar, phi :: CellValue) = gradientTerm2DPolar(phi)
gradientTerm(::Mesh3D, phi :: CellValue) = gradientTerm3D(phi)
gradientTerm(::Mesh3DCylindrical, phi :: CellValue) = gradientTerm3DCylindrical(phi)

function gradientTerm1D{T<:Real}(phi :: CellValue{T})
    dx = (phi.domain.cellsize.x[1:end-1] .+ phi.domain.cellsize.x[2:end]) ./ 2
    FaceValue(phi.domain,
      (phi.value[2:end]-phi.value[1:end-1])./dx,
      [one(T)],
      [one(T)])
end

function gradientTerm2D{T<:Real}(phi :: CellValue{T})
    dx = (phi.domain.cellsize.x[1:end-1] .+ phi.domain.cellsize.x[2:end]) ./ 2
    Ny = phi.domain.dims[2]
    dy = Array(T, 1, Ny+1)
    dy[:] = (phi.domain.cellsize.y[1:end-1] .+ phi.domain.cellsize.y[2:end]) ./ 2
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1]-phi.value[1:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end]-phi.value[2:end-1,1:end-1])./dy,
      [one(T)])
end

function gradientTerm2DPolar{T<:Real}(phi :: CellValue{T})
    dx = (phi.domain.cellsize.x[1:end-1] .+ phi.domain.cellsize.x[2:end]) ./ 2
    Ntheta = phi.domain.dims[2]
    dtheta = Array(T, 1, Ntheta+1)
    dtheta[:] = (phi.domain.cellsize.y[1:end-1] .+ phi.domain.cellsize.y[2:end]) ./ 2
    rp = phi.domain.cellcenters.x
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1]-phi.value[1:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end]-phi.value[2:end-1,1:end-1])./(dtheta.*rp),
      [one(T)])
end

function gradientTerm3D{T<:Real}(phi :: CellValue{T})
    Ny = phi.domain.dims[2]
    Nz = phi.domain.dims[3]
    dx = (phi.domain.cellsize.x[1:end-1] .+ phi.domain.cellsize.x[2:end]) ./ 2
    dy= Array(T, 1, Ny+1)
    dy[:] = (phi.domain.cellsize.y[1:end-1] .+ phi.domain.cellsize.y[2:end]) ./ 2
    dz= Array(T, 1, 1, Nz+1)
    dz[:] = (phi.domain.cellsize.z[1:end-1]+phi.domain.cellsize.z[2:end]) ./ 2
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1,2:end-1].-phi.value[1:end-1,2:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end,2:end-1].-phi.value[2:end-1,1:end-1,2:end-1])./dy,
      (phi.value[2:end-1,2:end-1,2:end].-phi.value[2:end-1,2:end-1,1:end-1])./dz)
end

function gradientTerm3DCylindrical{T<:Real}(phi :: CellValue{T})
    Ntheta = phi.domain.dims[2]
    Nz = phi.domain.dims[3]
    dx = (phi.domain.cellsize.x[1:end-1] .+ phi.domain.cellsize.x[2:end]) ./ 2
    dy= Array(T, 1, Ntheta+1)
    dy[:] = (phi.domain.cellsize.y[1:end-1] .+ phi.domain.cellsize.y[2:end]) ./ 2
    dz= Array(T, 1, 1, Nz+1)
    dz[:] = (phi.domain.cellsize.z[1:end-1] .+ phi.domain.cellsize.z[2:end]) ./ 2
    rp = phi.domain.cellcenters.x
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1,2:end-1].-phi.value[1:end-1,2:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end,2:end-1].-phi.value[2:end-1,1:end-1,2:end-1])./(dy.*rp),
      (phi.value[2:end-1,2:end-1,2:end].-phi.value[2:end-1,2:end-1,1:end-1])./dz)
end

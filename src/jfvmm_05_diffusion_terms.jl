# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================================================
# Changes:
#    2014-12-30 added 2D radial and 3D cylindrical grids
#    2015-01-10 extended to accept nonuniform grids
# ================================================================

# %% ======================= DIFFUSION TERM ================================
"""
Constructs the matrix for the diffusion term.

```julia
MD = diffusionTerm(D)
```

Input:

- `D`: diffusion coefficent (FaceValue)

Output:

- `MD` : diffusion matrix
"""
diffusionTerm{T<:Real}(D::FaceValue{T}) = diffusionTerm(D.domain.meshtype, D)
diffusionTerm(:: Mesh1D, D) = diffusionTerm1D(D)
diffusionTerm(:: Mesh1DPolar, D) = diffusionTermCylindrical1D(D)
diffusionTerm(:: Mesh2D, D) = diffusionTerm2D(D)[1]
diffusionTerm(:: Mesh2DPolar, D) = diffusionTermRadial2D(D)[1]
diffusionTerm(:: Mesh2DCylindrical, D) = diffusionTermCylindrical2D(D)[1]
diffusionTerm(:: Mesh3D, D) = diffusionTerm3D(D)[1]
diffusionTerm(:: Mesh3DCylindrical, D) = diffusionTermCylindrical3D(D)[1]

# %% ======================== 2D CARTESIAN DIFFUSION =========================
function diffusionTerm2D{T<:Real}(D::FaceValue{T})
    # extract data from the mesh structure
    Nx, Ny = D.domain.dims
    G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
    # DX, DY: cell sizes
    # dx, dy: ???
    DX = D.domain.cellsize.x
    DY = D.domain.cellsize.y
    dx = (DX[1:end-1] .+ DX[2:end]) ./ 2
    dy = (DY[1:end-1] .+ DY[2:end]) ./ 2
    # define the vectors to store the sparse matrix data
    iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
    iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
    jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
    jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
    sx = zeros(T, 3*(Nx+2)*(Ny+2))
    sy = zeros(T, 3*(Nx+2)*(Ny+2))
    mnx = Nx*Ny
    mny = Nx*Ny
    # reassign the east, west for code readability (use broadcasting in Julia)
    # all these values are matrices of Nx × Ny size
    De = D.xvalue[2:Nx+1,:] ./ (dx[2:Nx+1] .* DX[2:Nx+1])
    Dw = D.xvalue[1:Nx,:] ./ (dx[1:Nx] .* DX[2:Nx+1])
    # need to transpose dy & DY to get right dimensions
    Dn = D.yvalue[:,2:Ny+1] ./ transpose(dy[2:Ny+1] .* DY[2:Ny+1])
    Ds = D.yvalue[:,1:Ny] ./ transpose(dy[1:Ny] .* DY[2:Ny+1])
    # calculate the coefficients for the internal cells
    # vectorise matrices
    AE = reshape(De,mnx)
    AW = reshape(Dw,mnx)
    AN = reshape(Dn,mny)
    AS = reshape(Ds,mny)
    APx = -(AE+AW)
    APy = -(AN+AS)
    # build the sparse matrix based on the numbering system
    rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
    iix[1:3*mnx] = repmat(rowx_index,3)
    rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
    iiy[1:3*mny] = repmat(rowy_index,3)
    jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
    jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
    sx[1:3*mnx] = [AW; APx; AE]
    sy[1:3*mny] = [AS; APy; AN]

    # build the sparse matrix
    kx = 3*mnx
    ky = 3*mny
    Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
    My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
    M = Mx + My
    (M, Mx, My)
end

# %% ======================== 3D CARTESIAN DIFFUSION =========================
function diffusionTerm3D{T<:Real}(D::FaceValue{T})
    Nx, Ny, Nz = D.domain.dims
    G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
    DX = D.domain.cellsize.x
    DY = Array(T, 1, Ny+2)
    DY[:] = D.domain.cellsize.y
    DZ = Array(T, 1,1,Nz+2)
    DZ[:] = D.domain.cellsize.z
    dx = 0.5*(DX[1:end-1]+DX[2:end])
    dy = 0.5*(DY[:,1:end-1]+DY[:,2:end])
    dz = 0.5*(DZ[:,:,1:end-1]+DZ[:,:,2:end])
    # define the vectors to store the sparse matrix data
    iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
    sz = zeros(T, 3*(Nx+2)*(Ny+2)*(Nz+2))
    mnx = Nx*Ny*Nz
    mny = Nx*Ny*Nz
    mnz = Nx*Ny*Nz
    # reassign the east, west, north, and south velocity vectors for the
    # code readability (use broadcasting)
    De = D.xvalue[2:Nx+1,:,:]./(dx[2:Nx+1].*DX[2:Nx+1])
    Dw = D.xvalue[1:Nx,:,:]./(dx[1:Nx].*DX[2:Nx+1])
    Dn = D.yvalue[:,2:Ny+1,:]./(dy[:,2:Ny+1].*DY[:,2:Ny+1])
    Ds = D.yvalue[:,1:Ny,:]./(dy[:,1:Ny].*DY[:,2:Ny+1])
    Df = D.zvalue[:,:,2:Nz+1]./(dz[:,:,2:Nz+1].*DZ[:,:,2:Nz+1])
    Db = D.zvalue[:,:,1:Nz]./(dz[:,:,1:Nz].*DZ[:,:,2:Nz+1])
    # calculate the coefficients for the internal cells
    AE = reshape(De,mnx)
    AW = reshape(Dw,mnx)
    AN = reshape(Dn,mny)
    AS = reshape(Ds,mny)
    AF = reshape(Df,mnz)
    AB = reshape(Db,mnz)
    APx = reshape(-(De+Dw),mnx)
    APy = reshape(-(Dn+Ds),mny)
    APz = reshape(-(Df+Db),mnz)

    # build the sparse matrix based on the numbering system
    rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
    iix[1:3*mnx] = repmat(rowx_index,3)
    rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
    iiy[1:3*mny] = repmat(rowy_index,3);
    rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
    iiz[1:3*mnz] = repmat(rowz_index,3)
    jjx[1:3*mnx] = [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx);
    		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx);
    		reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
    jjy[1:3*mny] = [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny);
    		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny);
    		reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
    jjz[1:3*mnz] = [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz);
    		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz);
    		reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)];
    sx[1:3*mnx] = [AW; APx; AE]
    sy[1:3*mny] = [AS; APy; AN]
    sz[1:3*mnz] = [AB; APz; AF]

    # build the sparse matrix
    kx = 3*mnx
    ky = 3*mny
    kz = 3*mnz
    Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
    My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
    Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
    M = Mx + My + Mz
    (M, Mx, My, Mz)
end

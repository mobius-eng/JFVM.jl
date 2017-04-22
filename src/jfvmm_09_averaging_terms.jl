# ===============================
# Written by AAE, TU Delft
# simulkade.com
# ===============================

# %%  ================== Linear averaging scheme ==================
"""
Calculates the face value approximation of cell field variable using
linear mean approximation

φ(face) ≈ 1/(V1 + V2) [φ(1) V2 + φ(2) V1]

```julia
φf = linearMean(φ)
```

Input:

- `φ` : cell field variable

Output:

- `φf` : face value field

"""
function linearMean{T<:Real}(phi::CellValue{T})
    d=length(phi.domain.dims)
    if d==1
      dx = phi.domain.cellsize.x
      FaceValue(phi.domain,
        (dx[2:end].*phi.value[1:end-1]+dx[1:end-1].*phi.value[2:end])./(dx[2:end]+dx[1:end-1]),
        [one(T)],
        [one(T)])
    elseif d==2
      dx = phi.domain.cellsize.x
      Ny = phi.domain.dims[2]
      dy = Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      FaceValue(phi.domain,
        (dx[2:end].*phi.value[1:end-1,2:end-1]+dx[1:end-1].*phi.value[2:end,2:end-1])./(dx[2:end]+dx[1:end-1]),
        (dy[:,2:end].*phi.value[2:end-1,1:end-1]+dy[:,1:end-1].*phi.value[2:end-1,2:end])./(dy[:,2:end]+dy[:,1:end-1]),
        [one(T)])
    elseif d==3
      Ny = phi.domain.dims[2]
      Nz = phi.domain.dims[3]
      dx = phi.domain.cellsize.x
      dy= Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      dz= Array(T, 1, 1, Nz+2)
      dz[:] = phi.domain.cellsize.z
      FaceValue(phi.domain,
        (dx[2:end].*phi.value[1:end-1,2:end-1,2:end-1]+dx[1:end-1].*phi.value[2:end,2:end-1,2:end-1])./(dx[2:end]+dx[1:end-1]),
        (dy[:,2:end].*phi.value[2:end-1,1:end-1,2:end-1]+dy[:,1:end-1].*phi.value[2:end-1,2:end,2:end-1])./(dy[:,1:end-1]+dy[:,2:end]),
        (dz[:,:,2:end].*phi.value[2:end-1,2:end-1,1:end-1]+dz[:,:,1:end-1].*phi.value[2:end-1,2:end-1,2:end])./(dz[:,:,1:end-1]+dz[:,:,2:end]))
    end
end

# %% ================== Aritmetic averaging scheme ==================
"""
Calculates the face value approximation of cell field variable using
arithmetic mean approximation

φ(face) ≈ 1/(V1 + V2) [φ(1) V1 + φ(2) V2]

```julia
φf = arithmeticMean(φ)
```

Input:

- `φ` : cell field variable

Output:

- `φf` : face value field

"""
function arithmeticMean{T<:Real}(phi::CellValue{T})
    d = length(phi.domain.dims)
    if d==1
      dx = phi.domain.cellsize.x
      FaceValue(phi.domain,
        (dx[1:end-1].*phi.value[1:end-1]+dx[2:end].*phi.value[2:end])./(dx[2:end]+dx[1:end-1]),
        [one(T)],
        [one(T)])
    elseif d==2
      dx = phi.domain.cellsize.x
      Ny = phi.domain.dims[2]
      dy = Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      FaceValue(phi.domain,
        (dx[1:end-1].*phi.value[1:end-1,2:end-1]+dx[2:end].*phi.value[2:end,2:end-1])./(dx[2:end]+dx[1:end-1]),
        (dy[:,1:end-1].*phi.value[2:end-1,1:end-1]+dy[:,2:end].*phi.value[2:end-1,2:end])./(dy[:,2:end]+dy[:,1:end-1]),
        [one(T)])
    elseif d==3
      Ny = phi.domain.dims[2]
      Nz = phi.domain.dims[3]
      dx = phi.domain.cellsize.x
      dy= Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      dz= Array(T, 1, 1, Nz+2)
      dz[:] = phi.domain.cellsize.z
      FaceValue(phi.domain,
        (dx[1:end-1].*phi.value[1:end-1,2:end-1,2:end-1]+dx[2:end].*phi.value[2:end,2:end-1,2:end-1])./(dx[2:end]+dx[1:end-1]),
        (dy[:,1:end-1].*phi.value[2:end-1,1:end-1,2:end-1]+dy[:,2:end].*phi.value[2:end-1,2:end,2:end-1])./(dy[:,1:end-1]+dy[:,2:end]),
        (dz[:,:,1:end-1].*phi.value[2:end-1,2:end-1,1:end-1]+dz[:,:,2:end].*phi.value[2:end-1,2:end-1,2:end])./(dz[:,:,1:end-1]+dz[:,:,2:end]))
    end
end

# ================== Geometric averaging scheme ==================
"""
Calculates the face value approximation of cell field variable using
geometric mean approximation

φ(face) ≈ φ(1)^(V1/(V1+V2)) ⋅ φ(2)^(V2/(V1+V2))

(φ(face) ≈ √[φ(1) ⋅ φ(2)] for V1 = V2)

```julia
φf = geometricMean(φ)
```

Input:

- `φ` : cell field variable (>= 0)

Output:

- `φf` : face value field

"""
function geometricMean{T<:Real}(phi::CellValue{T})
    d=length(phi.domain.dimens)
    if d==1
      dx = phi.domain.cellsize.x
      n= phi.domain.dims[1]
      phix=zeros(T, n+1)
      for i=1:n+1
          if !(phi.value[i] == 0 || phi.value[i+1] == 0)
              phix[i]=exp((dx[i]*log(phi.value[i])+dx[i+1]*log(phi.value[i+1]))/(dx[i+1]+dx[i]))
          end
      end
      FaceValue(phi.domain,
        phix,
        [one(T)],
        [one(T)])
    elseif d==2
      dx = phi.domain.cellsize.x
      Ny = phi.domain.dims[2]
      dy = Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      FaceValue(phi.domain,
        exp((dx[1:end-1].*log(phi.value[1:end-1,2:end-1])+dx[2:end].*log(phi.value[2:end,2:end-1]))./(dx[2:end]+dx[1:end-1])),
        exp((dy[:,1:end-1].*log(phi.value[2:end-1,1:end-1])+dy[:,2:end].*log(phi.value[2:end-1,2:end]))./(dy[:,2:end]+dy[:,1:end-1])),
        [one(T)])
    elseif d==3
      Ny = phi.domain.dims[2]
      Nz = phi.domain.dims[3]
      dx = phi.domain.cellsize.x
      dy= Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      dz= Array(T, 1, 1, Nz+2)
      dz[:] = phi.domain.cellsize.z
      FaceValue(phi.domain,
        exp((dx[1:end-1].*log(phi.value[1:end-1,2:end-1,2:end-1])+dx[2:end].*log(phi.value[2:end,2:end-1,2:end-1]))./(dx[2:end]+dx[1:end-1])),
        exp((dy[:,1:end-1].*log(phi.value[2:end-1,1:end-1,2:end-1])+dy[:,2:end].*log(phi.value[2:end-1,2:end,2:end-1]))./(dy[:,1:end-1]+dy[:,2:end])),
        exp((dz[:,:,1:end-1].*log(phi.value[2:end-1,2:end-1,1:end-1])+dz[:,:,2:end].*log(phi.value[2:end-1,2:end-1,2:end]))./(dz[:,:,1:end-1]+dz[:,:,2:end])))
    end
end

# %% ================== Harmonic averaging scheme ==================
"""
Calculates the face value approximation of cell field variable using
harmonic mean approximation

φ(face) ≈ (V₁+v₂) / [V₁/φ₁ + V₂/φ₂]

```julia
φf = harmonicMean(φ)
```

Input:

- `φ` : cell field variable (>= 0)

Output:

- `φf` : face value field

"""
function harmonicMean{T<:Real}(phi::CellValue{T})
    d=length(phi.domain.dims)
    if d==1
      dx = phi.domain.cellsize.x
      n=phi.domain.dims[1]
      phix=zeros(T,n+1)
      for i=1:n+1
          if !(phi.value[i] == 0 || phi.value[i+1] == 0)
              phix[i]=(dx[i+1]+dx[i])/(dx[i+1]/phi.value[i+1]+dx[i]/phi.value[i])
          end
      end
      FaceValue(phi.domain,
        phix,
        [one(T)],
        [one(T)])
    elseif d==2
      dx = phi.domain.cellsize.x
      Ny = phi.domain.dims[2]
      dy = Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      FaceValue(phi.domain,
        phi.value[2:end,2:end-1].*phi.value[1:end-1,2:end-1].*(dx[2:end]+dx[1:end-1])./(dx[2:end].*phi.value[1:end-1,2:end-1]+dx[1:end-1].*phi.value[2:end,2:end-1]),
        phi.value[2:end-1,2:end].*phi.value[2:end-1,1:end-1].*(dy[:,2:end]+dy[:,1:end-1])./(dy[:,2:end].*phi.value[2:end-1,1:end-1]+dy[:,1:end-1].*phi.value[2:end-1,2:end]),
        [one(T)])
    else # d==3
      Ny = phi.domain.dims[2]
      Nz = phi.domain.dims[3]
      dx = phi.domain.cellsize.x
      dy= Array(T, 1, Ny+2)
      dy[:] = phi.domain.cellsize.y
      dz= Array(T, 1, 1, Nz+2)
      dz[:] = phi.domain.cellsize.z
      FaceValue(phi.domain,
        phi.value[2:end,2:end-1,2:end-1].*phi.value[1:end-1,2:end-1,2:end-1].*(dx[2:end]+dx[1:end-1])./(dx[2:end].*phi.value[1:end-1,2:end-1,2:end-1]+dx[1:end-1].*phi.value[2:end,2:end-1,2:end-1]),
        phi.value[2:end-1,2:end,2:end-1].*phi.value[2:end-1,1:end-1,2:end-1].*(dy[:,1:end-1]+dy[:,2:end])./(dy[:,2:end].*phi.value[2:end-1,1:end-1,2:end-1]+dy[:,1:end-1].*phi.value[2:end-1,2:end,2:end-1]),
        phi.value[2:end-1,2:end-1,2:end].*phi.value[2:end-1,2:end-1,1:end-1].*(dz[:,:,1:end-1]+dz[:,:,2:end])./(dz[:,:,2:end].*phi.value[2:end-1,2:end-1,1:end-1]+dz[:,:,1:end-1].*phi.value[2:end-1,2:end-1,2:end]))
    end
end

# %% ================== Upwind averaging scheme ==================
"""
Calculates the face value approximation of cell field variable using
upwind mean as follows:
- Use cell value from the cell from which flow is occuring
- If the flow is 0, uses simple arithmetic average

```julia
φf = upwindMean(φ, u)
```

Input:

- `φ` : cell field variable (>= 0)
- `u` : face field of convective velocity

Output:

- `φf` : face value field

"""
function upwindMean{T<:Real}(phi::CellValue{T}, u::FaceValue{T})
    d=length(phi.domain.dims)
    phi_tmp = Base.copy(phi.value)
    if d==1
      ux=u.xvalue
      # assign the value of the left boundary to the left ghost cell
      phi_tmp[1] = (phi.value[1] .+ phi.value[2]) ./ 2
      # assign the value of the right boundary to the right ghost cell
      phi_tmp[end] = (phi.value[end] .+ phi.value[end-1]) ./ 2
      FaceValue(phi.domain,
        (ux.>0).*phi_tmp[1:end-1]+(ux.<0).*phi_tmp[2:end]+
        (ux.==0).*(phi.value[1:end-1]+phi.value[2:end]) ./ 2,
        [one(T)],
        [one(T)])
    elseif d==2
      ux = u.xvalue
      uy = u.yvalue
      # assign the value of the left boundary to the left ghost cells
      phi_tmp[1,:] = (phi.value[1,:] .+ phi.value[2,:]) ./ 2
      # assign the value of the right boundary to the right ghost cells
      phi_tmp[end,:] = (phi.value[end,:] .+ phi.value[end-1,:]) ./ 2
      # assign the value of the bottom boundary to the bottom ghost cells
      phi_tmp[:,1] = (phi.value[:,1] .+ phi.value[:,2]) ./ 2
      # assign the value of the top boundary to the top ghost cells
      phi_tmp[:,end] = (phi.value[:,end] .+ phi.value[:,end-1]) ./ 2
      FaceValue(phi.domain,
        (ux.>0).*phi_tmp[1:end-1,2:end-1]+
        (ux.<0).*phi_tmp[2:end,2:end-1]+
        (ux.==0).* (phi.value[1:end-1,2:end-1] .+ phi.value[2:end,2:end-1]) ./ 2,
        (uy.>0).*phi_tmp[2:end-1,1:end-1]+
        (uy.<0).*phi_tmp[2:end-1,2:end]+
        (uy.==0.0).*(phi.value[2:end-1,1:end-1]+phi.value[2:end-1,2:end]) ./ 2,
        [one(T)])
    else # d==3
      ux = u.xvalue
      uy = u.yvalue
      uz = u.zvalue
      # assign the value of the left boundary to the left ghost cells
      phi_tmp[1,:,:] = (phi.value[1,:,:] .+ phi.value[2,:,:]) ./ 2
      # assign the value of the right boundary to the right ghost cells
      phi_tmp[end,:,:] = (phi.value[end,:,:] .+ phi.value[end-1,:,:]) ./ 2
      # assign the value of the bottom boundary to the bottom ghost cells
      phi_tmp[:,1,:] = (phi.value[:,1,:] .+ phi.value[:,2,:]) ./ 2
      # assign the value of the top boundary to the top ghost cells
      phi_tmp[:,end,:] = (phi.value[:,end,:] .+ phi.value[:,end-1,:]) ./ 2
      # assign the value of the back boundary to the back ghost cells
      phi_tmp[:,:,1] = (phi.value[:,:,1] .+ phi.value[:,:,2]) ./ 2
      # assign the value of the front boundary to the front ghost cells
      phi_tmp[:,:,end] = (phi.value[:,:,end] .+ phi.value[:,:,end-1]) ./ 2
      FaceValue(phi.domain,
        (ux.>0).*phi_tmp[1:end-1,2:end-1,2:end-1]+
        (ux.<0).*phi_tmp[2:end,2:end-1,2:end-1]+
        (ux.==0).*(phi.value[1:end-1,2:end-1,2:end-1]+phi.value[2:end,2:end-1,2:end-1]) / 2,
        (uy.>0).*phi_tmp[2:end-1,1:end-1,2:end-1]+
        (uy.<0).*phi_tmp[2:end-1,2:end,2:end-1]+
        (uy.==0.0).*(phi.value[2:end-1,1:end-1,2:end-1]+phi.value[2:end-1,2:end,2:end-1]) / 2,
        (uz.>0).*phi_tmp[2:end-1,2:end-1,1:end-1]+
        (uz.<0).*phi_tmp[2:end-1,2:end-1,2:end]+
        (uz.==0).*(phi.value[2:end-1,2:end-1,1:end-1]+phi.value[2:end-1,2:end-1,2:end]) / 2)
    end
end

# %% ================== TVD averaging scheme ==================
function tvdMean(phi::CellValue, u::FaceValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs(phi_in).<eps1).*sign(phi_in)

d=length(phi.domain.dims)

if d == 1 || d == 1.5
  # extract data from the mesh structure
  Nx = u.domain.dims[1]
  dx = 0.5*(u.domain.cellsize.x[1:end-1]+u.domain.cellsize.x[2:end])
  phi_p = zeros(Float64, Nx+1)
  phi_m = zeros(Float64, Nx+1)

  # extract the velocity data
  ux = u.xvalue

  # calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
  dphi_p = (phi.value[2:Nx+2]-phi.value[1:Nx+1])./dx
  rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
  phi_p[2:Nx+1] = phi.value[2:Nx+1]+0.5*FL(rp).*(phi.value[3:Nx+2]-phi.value[2:Nx+1])
  phi_p[1] = (phi.value[1]+phi.value[2])/2.0 # left boundary

  # calculate the upstream to downstream gradient ratios for u<0 (- ratio)
  rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
  phi_m[1:Nx] = phi.value[2:Nx+1]+0.5*FL(rm).*(phi.value[1:Nx]-phi.value[2:Nx+1])
  phi_m[Nx+1] = (phi.value[end]+phi.value[end-1])/2.0 # right boundary

  FaceValue(phi.domain,
      (ux.>0.0).*phi_p+(ux.<0.0).*phi_m+
      0.5*(ux.==0).*(phi.value[1:end-1]+phi.value[2:end]),
      [1.0],
      [1.0])
elseif d==2 || d==2.5 || d==2.8
  # extract data from the mesh structure
  Nx = u.domain.dims[1]
  Ny = u.domain.dims[2]
  dx=0.5*(u.domain.cellsize.x[1:end-1]+u.domain.cellsize.x[2:end])
  dy=Array(Float64, 1, Ny+1)
  dy[:]=0.5*(u.domain.cellsize.y[1:end-1]+u.domain.cellsize.y[2:end])
  phi_p = zeros(Float64, Nx+1)
  phi_m = zeros(Float64, Nx+1)
  phiX_p = zeros(Float64, Nx+1, Ny)
  phiX_m = zeros(Float64, Nx+1,Ny)
  phiY_p = zeros(Float64, Nx,Ny+1)
  phiY_m = zeros(Float64, Nx,Ny+1)

  # extract the velocity data
  ux = u.xvalue
  uy = u.yvalue

  # calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
  # x direction
  dphiX_p = (phi.value[2:Nx+2, 2:Ny+1]-phi.value[1:Nx+1, 2:Ny+1])./dx
  rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
  phiX_p[2:Nx+1,:] = phi.value[2:Nx+1, 2:Ny+1]+0.5*FL(rX_p).*
  (phi.value[3:Nx+2,2:Ny+1]-phi.value[2:Nx+1, 2:Ny+1])
  phiX_p[1, :] = (phi.value[1, 2:Ny+1]+phi.value[2, 2:Ny+1])/2.0  # left boundary
  # y direction
  dphiY_p = (phi.value[2:Nx+1, 2:Ny+2]-phi.value[2:Nx+1, 1:Ny+1])./dy
  rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
  phiY_p[:,2:Ny+1] = phi.value[2:Nx+1, 2:Ny+1]+0.5*FL(rY_p).*
      (phi.value[2:Nx+1,3:Ny+2]-phi.value[2:Nx+1, 2:Ny+1])
  phiY_p[:,1] = (phi.value[2:Nx+1,1]+phi.value[2:Nx+1,2])/2.0  # Bottom boundary

  # calculate the upstream to downstream gradient ratios for u<0 (- ratio)
  # x direction
  rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
  phiX_m[1:Nx,:] = phi.value[2:Nx+1, 2:Ny+1]+0.5*FL(rX_m).*
      (phi.value[1:Nx, 2:Ny+1]-phi.value[2:Nx+1, 2:Ny+1])
  phiX_m[Nx+1,:] = (phi.value[end, 2:Ny+1]+phi.value[end-1, 2:Ny+1])/2.0  # right boundary
  # y direction
  rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
  phiY_m[:,1:Ny] = phi.value[2:Nx+1, 2:Ny+1]+0.5*FL(rY_m).*
      (phi.value[2:Nx+1, 1:Ny]-phi.value[2:Nx+1, 2:Ny+1])
  phiY_m[:, Ny+1] = (phi.value[2:Nx+1, end]+phi.value[2:Nx+1, end-1])/2.0  # top boundary

  FaceValue(phi.domain,
      (ux.>0.0).*phiX_p+(ux.<0.0).*phiX_m+
         0.5*(ux.==0.0).*(phi.value[1:Nx+1,2:Ny+1]+phi.value[2:Nx+2,2:Ny+1]),
      (uy.>0.0).*phiY_p+(uy.<0.0).*phiY_m+
         0.5*(uy.==0.0).*(phi.value[2:Nx+1,1:Ny+1]+phi.value[2:Nx+1,2:Ny+2]),
      [1.0])

elseif d==3 || d==3.2
  # extract data from the mesh structure
  Nx = u.domain.dims[1]
  Ny = u.domain.dims[2]
  Nz = u.domain.dims[3]
  dx=0.5*(u.domain.cellsize.x[1:end-1]+u.domain.cellsize.x[2:end])
  dy=Array(Float64, 1, Ny+1)
  dy[:]=0.5*(u.domain.cellsize.y[1:end-1]+u.domain.cellsize.y[2:end])
  dz=Array(Float64, 1, 1, Nz+1)
  dz[:]=0.5*(u.domain.cellsize.z[1:end-1]+u.domain.cellsize.z[2:end])
  # extract the velocity data
  ux = u.xvalue
  uy = u.yvalue
  uz = u.zvalue

  # define the tvd face vectors
  phiX_p = zeros(Float64, Nx+1,Ny,Nz)
  phiX_m = zeros(Float64, Nx+1,Ny,Nz)
  phiY_p = zeros(Float64, Nx,Ny+1,Nz)
  phiY_m = zeros(Float64, Nx,Ny+1,Nz)
  phiZ_p = zeros(Float64, Nx,Ny,Nz+1)
  phiZ_m = zeros(Float64, Nx,Ny,Nz+1)

  # calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
  # x direction
  dphiX_p = (phi.value[2:Nx+2, 2:Ny+1, 2:Nz+1]-phi.value[1:Nx+1, 2:Ny+1, 2:Nz+1])./dx
  rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
  phiX_p[2:Nx+1,:,:] = phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rX_p).*
      (phi.value[3:Nx+2,2:Ny+1,2:Nz+1]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
  phiX_p[1,:,:] = (phi.value[1,2:Ny+1,2:Nz+1]+phi.value[2,2:Ny+1,2:Nz+1])/2.0  # left boundary
  # y direction
  dphiY_p = (phi.value[2:Nx+1, 2:Ny+2, 2:Nz+1]-phi.value[2:Nx+1, 1:Ny+1, 2:Nz+1])./dy
  rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
  phiY_p[:,2:Ny+1,:] = phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rY_p).*
      (phi.value[2:Nx+1,3:Ny+2,2:Nz+1]-phi.value[2:Nx+1, 2:Ny+1,2:Nz+1])
  phiY_p[:,1,:] = (phi.value[2:Nx+1,1,2:Nz+1]+phi.value[2:Nx+1,2,2:Nz+1])/2.0  # Bottom boundary
  # z direction
  dphiZ_p = (phi.value[2:Nx+1, 2:Ny+1, 2:Nz+2]-phi.value[2:Nx+1, 2:Ny+1, 1:Nz+1])./dz
  rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
  phiZ_p[:,:,2:Nz+1] = phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rZ_p).*
      (phi.value[2:Nx+1,2:Ny+1,3:Nz+2]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
  phiZ_p[:,:,1] = (phi.value[2:Nx+1,2:Ny+1,1]+phi.value[2:Nx+1,2:Ny+1,2])/2.0  # Back boundary

  # calculate the upstream to downstream gradient ratios for u<0 (- ratio)
  # x direction
  rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
  phiX_m[1:Nx,:,:] = phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1]+0.5*FL(rX_m).*
      (phi.value[1:Nx, 2:Ny+1, 2:Nz+1]-phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1])
  phiX_m[Nx+1,:,:] = (phi.value[end,2:Ny+1,2:Nz+1]+phi.value[end-1,2:Ny+1,2:Nz+1])/2.0  # right boundary
  # y direction
  rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
  phiY_m[:,1:Ny,:] = phi.value[2:Nx+1,2:Ny+1,2:Nz+1]+0.5*FL(rY_m).*
      (phi.value[2:Nx+1,1:Ny,2:Nz+1]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
  phiY_m[:,Ny+1,:] = (phi.value[2:Nx+1, end,2:Nz+1]+phi.value[2:Nx+1, end-1,2:Nz+1])/2.0  # top boundary
  # z direction
  rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
  phiZ_m[:,:,1:Nz] = phi.value[2:Nx+1,2:Ny+1,2:Nz+1]+0.5*FL(rZ_m).*
      (phi.value[2:Nx+1,2:Ny+1,1:Nz]-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
  phiZ_m[:,:,Nz+1] = (phi.value[2:Nx+1,2:Ny+1,end]+phi.value[2:Nx+1,2:Ny+1,end-1])/2.0  # front boundary

  FaceValue(phi.domain,
      (ux.>0.0).*phiX_p+(ux.<0.0).*phiX_m+
         0.5*(ux.==0.0).*(phi.value[1:Nx+1,2:Ny+1,2:Nz+1]+phi.value[2:Nx+2,2:Ny+1,2:Nz+1]),
      (uy.>0.0).*phiY_p+(uy.<0).*phiY_m+
         0.5*(uy.==0.0).*(phi.value[2:Nx+1,1:Ny+1,2:Nz+1]+phi.value[2:Nx+1,2:Ny+2,2:Nz+1]),
      (uz.>0.0).*phiZ_p+(uz.<0).*phiZ_m+
         0.5*(uz.==0.0).*(phi.value[2:Nx+1,2:Ny+1,1:Nz+1]+phi.value[2:Nx+1,2:Ny+1,2:Nz+2]))
end

end

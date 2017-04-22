# ===============================
# Written by AAE
# TU Delft, Winter 2014
# simulkade.com
# ===============================

abstract FluxLimiter

type FlCHARM <: FluxLimiter end

const CHARM = FlCHARM()

type FlHCUS <: FluxLimiter end

const HCUS = FlHCUS()

type FlHQUICK <: FluxLimiter end

const HQUICK = FlHQUICK()

type FlOSPRE <: FluxLimiter end

const OSPRE = FlOSPRE()

type FlVanLeer <: FluxLimiter end

const VanLeer = FlVanLeer()

type FlVanAlbada1 <: FluxLimiter end

const VanAlbada1 = FlVanAlbada1()

type FlVanAlbada2 <: FluxLimiter end

const VanAlbada2 = FlVanAlbada2()

type FlMinMod <: FluxLimiter end

const MinMod = FlMinMod()

type FlSUPERBEE <: FluxLimiter end

const SUPERBEE = FlSUPERBEE()

type FlOsher <: FluxLimiter end

const Osher = FlOsher()

type FlSweby <: FluxLimiter end

const Sweby = FlSweby()

type FlSmart <: FluxLimiter end

const Smart = FlSmart()

type FlKoren <: FluxLimiter end

const Koren = FlKoren()

type FlMUSCL <: FluxLimiter end

const MUSCL = FlMUSCL()

type FlQUICK <: FluxLimiter end

const QUICK = FlQUICK()

type FlUMIST <: FluxLimiter end

const MIST = FlUMIST()

# This function returns a function handle to a flux limiter of user's
# choice.
# available flux limiters are: 'CHARM', 'HCUS', 'HQUICK', 'VanLeer',
# 'VanAlbada1', 'VanAlbada2', 'MinMod', 'SUPERBEE', 'Sweby', 'Osher',
# 'Koren', 'smart', 'MUSCL', 'QUICK', 'MC', and 'UMIST'.
# Default limiter is 'SUPERBEE'. See:
# <http://en.wikipedia.org/wiki/Flux_limiter>

# find the flux limiter function
fluxLimiter(::FlCHARM) = r->((r.>0.0).*r.*(3.0*r+1.0)./(((r+1.0).^2.0)+eps()*(r.==-1.0)))
fluxLimiter(::FlHCUS) = r->(1.5*(r+abs(r))./(r+2.0))

fluxLimiter(::FlHQUICK) = r->(2.0*(r+abs(r))./((r+3.0)+eps()*(r.==-3.0)))
fluxLimiter(::FlOSPRE) = r->((1.5*r.*(r+1.0))./(r.*(r+1.0)+1.0+eps()*((r.*(r+1.0)+1.0).==0.0)))
fluxLimiter(::FlVanLeer) = r->((r+abs(r))./(1.0+abs(r)))
fluxLimiter(::FlVanAlbada1) = r->((r+r.*r)./(1.0+r.*r))
fluxLimiter(::FlVanAlbada2) = r->(2.0*r./(1+r.*r))
fluxLimiter(::FlMinMod) = r->((r>0.0).*min(r,1.0))
fluxLimiter(::FlSUPERBEE) = r->(max(0.0, max(min(2.0*r,1.0), min(r,2.0))))
fluxLimiter(::FlOsher) = r->(max(0.0, min(r,1.5)))

function fluxLimiter(::FlSweby)
  b=1.5
  r->(max(0.0, max(min(b*r,1.0), min(r,b))))
end

fluxLimiter(::FlSmart) = r->(max(0.0, min(4.0,min(0.25+0.75*r, 2.0*r))))
fluxLimiter(::FlKoren) = r->(max(0.0, min(2.0*r, min((1.0+2.0*r)/3.0, 2.0))))
fluxLimiter(::FlMUSCL) = r->(max(0.0, min(2.0*r, min(0.5*(1+r), 2.0))))
fluxLimiter(::FlQUICK) = r->(max(0.0, min(2.0, min(2.0*r, (3.0+r)/4.0))))
fluxLimiter(::FlUMIST) = r->(max(0.0, min(2.0, min(2.0*r, min((1.0+3.0*r)/4.0, (3.0+r)/4.0)))))


# %% faceEval
"""
Applies function to a face value

```julia
xnew = faceEval(fn, x)
```

Input:

- `fn` : function to be applied, it will be applid to the whole array of values
- `x` : face field variable

Output:

- `xnew` : transformed face field variable

"""
function faceEval(f::Function, x::FaceValue)
    FaceValue(x.domain,
        f(x.xvalue),
        f(x.yvalue),
        f(x.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue)
    FaceValue(x1.domain,
        f(x1.xvalue, x2.xvalue),
        f(x1.yvalue, x2.yvalue),
        f(x1.zvalue, x2.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue)
    FaceValue(x1.domain,
        f(x1.xvalue, x2.xvalue, x3.xvalue),
        f(x1.yvalue, x2.yvalue, x3.yvalue),
        f(x1.zvalue, x2.zvalue, x3.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue, x4::FaceValue)
    FaceValue(x1.domain,
        f(x1.xvalue, x2.xvalue, x3.xvalue, x4.xvalue),
        f(x1.yvalue, x2.yvalue, x3.yvalue, x4.yvalue),
        f(x1.zvalue, x2.zvalue, x3.zvalue, x4.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue, x4::FaceValue, x5::FaceValue)
    FaceValue(x1.domain,
        f(x1.xvalue, x2.xvalue, x3.xvalue, x4.xvalue, x5.xvalue),
        f(x1.yvalue, x2.yvalue, x3.yvalue, x4.yvalue, x5.yvalue),
        f(x1.zvalue, x2.zvalue, x3.zvalue, x4.zvalue, x5.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue, x4::FaceValue, x5::FaceValue, x6::FaceValue)
    FaceValue(x1.domain,
        f(x1.xvalue, x2.xvalue, x3.xvalue, x4.xvalue, x5.xvalue, x6.xvalue),
        f(x1.yvalue, x2.yvalue, x3.yvalue, x4.yvalue, x5.yvalue, x6.yvalue),
        f(x1.zvalue, x2.zvalue, x3.zvalue, x4.zvalue, x5.zvalue, x6.zvalue))
end

"""
Applies function to cell field values

```julia
xnew = cellEval(fn, x)
```

Input:

- `fn` : function to be applied, it will be applid to the whole array of values
- `x` : cell field variable

Output:

- `xnew` : transformed cell field variable

"""
function cellEval(f::Function, x::CellValue)
    CellValue(x1.domain,
        f(x.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue)
    CellValue(x1.domain,
        f(x1.value, x2.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue)
    CellValue(x1.domain,
        f(x1.value, x2.value, x3.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue, x4::CellValue)
    CellValue(x1.domain,
        f(x1.value, x2.value, x3.value, x4.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue, x4::CellValue, x5::CellValue)
    CellValue(x1.domain,
        f(x1.value, x2.value, x3.value, x4.value, x5.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue, x4::CellValue, x5::CellValue, x6::CellValue)
    CellValue(x1.domain,
        f(x1.value, x2.value, x3.value, x4.value, x5.value, x6.value))
end

# %% ========================= Generate random perm field ======================
function permfieldlogrndg(Nx,k_avrg,V_dp,cl)
  # 1D random field generator:
  # hdf: Gaussian
  # acf: Gaussian
  Lx=1.0
  x = linspace(-Lx/2,Lx/2,Nx)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2.0
  Z = s*randn(Nx)

  # Gaussian filter
  F = exp(-x.^2/(cl^2/2.0))

  # correlation of surface using convolution (faltung), inverse
  # Fourier transform and normalizing prefactors
  f = sqrt(2.0/sqrt(pi))*sqrt(Lx/Nx/cl)*ifft(fft(Z).*fft(F))
  perm = exp(mu+real(f))
end

function permfieldlogrnde(Nx,k_avrg,V_dp,cl)
  # 1D random field generator:
  # hdf: Gaussian
  # acf: Exponential
  Lx=1.0
  x = linspace(-Lx/2,Lx/2,Nx)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2.0
  Z = s*randn(Nx)

  # Gaussian filter
  F = exp(-abs(x)/(cl/2.0))

  # correlation of surface using convolution (faltung), inverse
  # Fourier transform and normalizing prefactors
  f = sqrt(2.0)*sqrt(Lx/Nx/cl)*ifft(fft(Z).*fft(F))
  perm = exp(mu+real(f))
end


# ======================== 2D ==================================>
function permfieldlogrndg(Nx,Ny,k_avrg,V_dp,clx,cly)
# 2D random field generator:
# hdf: Gaussian
# acf: Gaussian
# The surface has a Gaussian height distribution function
# and Gaussian autocovariance functions
  Lx=1.0
  X = linspace(-Lx/2,Lx/2,Nx)
  Y = linspace(-Lx/2,Lx/2,Ny)'

  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  Z = s*randn(Nx,Ny)
  # Gaussian filter
  F = exp(-(X.^2/(clx^2/2.0).+Y.^2/(cly^2/2.0)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0/sqrt(pi)*Lx/sqrt(Nx*Ny)/sqrt(clx)/sqrt(cly)*ifft(fft(Z).*fft(F))
  perm = exp(mu+real(f))
end

function permfieldlogrnde(Nx,Ny,k_avrg,V_dp,clx,cly)
# 2D random field generator:
# hdf: Gaussian
# acf: Exponential
# The surface has a Gaussian height distribution function
# and Exponential autocovariance functions
  Lx=1.0
  X = linspace(-Lx/2,Lx/2,Nx)
  Y = linspace(-Lx/2,Lx/2,Ny)'

  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  Z = s*randn(Nx,Ny)
  # Gaussian filter
  F = exp(-(abs(X)/(clx/2).+abs(Y)/(cly/2)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0*Lx/sqrt(Nx*Ny)/sqrt(clx*cly)*ifft(fft(Z).*fft(F))
  perm = exp(mu+real(f))
end
# <========================= 2D ================================


# ======================== 3D ==================================>
function permfieldlogrndg(Nx,Ny,Nz,k_avrg,V_dp,clx,cly,clz)
# 2D random field generator:
# hdf: Gaussian
# acf: Gaussian
# The surface has a Gaussian height distribution function
# and Gaussian autocovariance functions
  Lx=1.0
  X = linspace(-Lx/2,Lx/2,Nx)
  Y = linspace(-Lx/2,Lx/2,Ny)'
  Z = zeros(1,1,Nz)
  Z[:]=linspace(-Lx/2,Lx/2,Nz)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  z = s*randn(Nx,Ny,Nz)
  # Gaussian filter
  F = exp(-(X.^2/(clx^2/2.0).+Y.^2/(cly^2/2.0).+Z.^2/(clz^2/2.0)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0/sqrt(pi)*Lx/(Nx*Ny*Nz)^(1/3)/(clx*cly*clz)^(1/3)*ifft(fft(z).*fft(F))
  perm = exp(mu+real(f))
end

function permfieldlogrnde(Nx,Ny,Nz,k_avrg,V_dp,clx,cly,clz)
# 2D random field generator:
# hdf: Gaussian
# acf: Exponential
# The surface has a Gaussian height distribution function
# and Exponential autocovariance functions
  Lx=1.0
  X = linspace(-Lx/2,Lx/2,Nx)
  Y = linspace(-Lx/2,Lx/2,Ny)'
  Z = zeros(1,1,Nz)
  Z[:]=linspace(-Lx/2,Lx/2,Nz)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  z = s*randn(Nx,Ny,Nz)
  # Gaussian filter
  F = exp(-(abs(X)/(clx/2.0).+abs(Y)/(cly/2.0).+abs(Z)/(clz/2.0)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0*Lx/(Nx*Ny*Nz)^(1/3)/(clx*cly*clz)^(1/3)*ifft(fft(z).*fft(F))
  perm = exp(mu+real(f))
end
# <============================== 3D ==============================

"""
Returns the volume of each cell in the form of a cell variable

```julia
cellvol = cellVolume(m::MeshStructure)
```
"""
cellVolume(m::MeshStructure) =
    createCellVariable(m, calcCellVolume(m.meshtype, m), createBC(m))

calcCellVolume(::Mesh1D, m) = m.cellsize.x[2:end-1]
calcCellVolume(::Mesh1DPolar, m) = 2π*m.cellsize.x[2:end-1].*m.cellcenters.x
calcCellVolume(::Mesh2D, m) = m.cellsize.x[2:end-1] * m.cellsize.y[2:end-1]'
calcCellVolume(::Mesh2DCylindrical, m) = 2π*m.cellcenters.x.*m.cellsize.x[2:end-1]*m.cellsize.y[2:end-1]'
calcCellVolume(::Mesh2DPolar, m) = cellcenters.x.*m.cellsize.x[2:end-1]*m.cellsize.y[2:end-1]'
calcCellVolume(::Mesh3D, m) = error("Not implemented yet")
calcCellVolume(::Mesh3DCylindrical, m) = error("Not implemented yet")

"""
this function reshapes a vetorized cell variable to its domain shape
matrix based on the mesh structure data; it is assumed that the phi
includes the ghost cell data as well.
"""
function reshapeCell(m, phi)
  reshape(full(x), tuple(m.dims+2...))
end

"""
this function reshapes a vetorized cell variable to its domain shape
matrix based on the mesh structure data; it is assumed that the phi
does NOT include the ghost cell data.
"""
function reshapeInternalCell(m, phi)
  reshape(full(x), tuple(m.dims...))
end

"""
returns the internal cells of a cell variable as an array of the same shape
"""
function internalCells{T<:Real}(phi::CellValue{T})
  d = length(phi.domain.dims)
  N = phi.domain.dims

  if d==1
  	cellvar= phi.value[2:N[1]+1]
  elseif d==2
  	cellvar= phi.value[2:N[1]+1, 2:N[2]+1]
  else #(d==3)
      cellvar= phi.value[2:N[1]+1, 2:N[2]+1, 2:N[3]+1]
  end
  return cellvar
end

"""
Integrate variable phi over the domain it is defined
"""
function domainInt(phi::CellValue)
  return sum(internalCells(phi).*internalCells(cellVolume(phi.domain)))
end

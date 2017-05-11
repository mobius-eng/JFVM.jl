__precompile__()
module BoundaryConditions

using JFVMM.Geometry
using JFVMM.MeshStructure

export AbstractBoundaryValue, DirichletBV,VonNeumannV, MixedBV
export BoundaryCondition, createBoundaryCondition
export boundaryvalue

"""
Abstract representation of the boundary value
"""
abstract AbstractBoundaryValue

"""
Dirichlet boundary value: the value of the quanitity
"""
immutable DirichletBV <: AbstractBoundaryValue
    value :: Float64
end

"""
Von Neumann (natural) boudnary value: flux perpendicular to boundary surface area.

Convention: positive flux is out (in the direction of the boundary surface vector).
"""
immutable VonNeumannBV <: AbstractBoundaryValue
    fluxout :: Float64
end

"""
Mixed boundary value: usually it is the flux due to some external convection of the form

```
h∞ (φ∞ - φ)
```
Transfer coefficient and the the external value must be supplied
"""
immutable MixedBV <: AbstractBoundaryValue
    coeff :: Float64
    value :: Float64
end

"""
Boundary condition over the domain
"""
type BoundaryCondition
    domain :: Mesh
    ifacemin :: Vector{AbstractBoundaryValue}
    ifacemax :: Vector{AbstractBoundaryValue}
    jfacemin :: Vector{AbstractBoundaryValue}
    jfacemax :: Vector{AbstractBoundaryValue}
end

"""
Creates default (natural with zero flux) boundary condition for the mesh
"""
function createBoundaryCondition(m :: Mesh)
    Ni, Nj = size(m)
    BoundaryCondition(m,
        [VonNeumannBV(0.0) for j=1:Nj],
        [VonNeumannBV(0.0) for j=1:Nj],
        [VonNeumannBV(0.0) for i=1:Ni],
        [VonNeumannBV(0.0) for i=1:Ni])
end


function boundaryvalue(bv :: DirichletBV, face :: Face, phi :: Float64, gamma :: Float64, gardphi :: Point2D = Point2D(0,0))
    bv.value
end

function boundaryvalue(bv :: VonNeumannBV, face :: Face, phi :: Float64, gamma :: Float64, gradphi :: Point2D = Point2D(0,0))
    coeff = gamma * face.gDiff
    (coeff * phi - bv.fluxout) / coeff
end

function boundaryvalue(bv :: MixedBV, face :: Face, phi :: Float64, gamma :: Float64, gradphi :: Point2D = Point2D(0,0))
    a = bv.coeff * face.area
    b = gamma * face.gDiff
    c = gamma * (gradphi ⋅ (face.Sf - face.Ef))
    (a*bv.value + b * phi - c) / (a+b)
end

end

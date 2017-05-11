__precompile__()
module DiffusionTerm

using JFVMM.Geometry
using JFVMM.MeshStructure
using JFVMM.BoundaryConditions
using JFVMM.Fields

"""
Assign matrix of coefficients and RHS vector values corresponding to diffusion term for one face
"""
function assignFaceDiffTerm!(A :: AbstractMatrix{Float64}, b :: AbstractVector{Float64}, m :: Mesh, Γ :: Array{Float64,2}, ∇φf :: Array{Point2D,2}, f :: Face)
    I = meshindex(f)
    owneri = faceowneri(f)
    neighbouri = faceneighbouri(f)
    aF = - Γ[I] * f.gDiff
    fluxVf = Γ[I] * (∇φf ⋅ (f.Sf - f.Ef))
    # owner:
    A[globalindex(m, owneri), globalindex(m, neighbouri)] += aF
    A[globalindex(m, owneri), globalindex(m, owneri)] += -aF
    b[globalindex(m,owneri)] += -fluxVf
    # neighbour
    A[globalindex(m,neighbouri), globalindex(m,owneri)] += aF
    A[globalindex(m,neighbouri), globalindex(m,neighbouri)] += -aF
    b[globalindex(m,neighbouri)] += fluxVf
end

"""
Assigns coefficient matrix and RHS vector values corresponding to diffusion term for internal
faces. BCs need to be incorporated separately
"""
function diffusionTerm!(A :: AbstractMatrix{Float64}, b :: AbstractVector{Float64}, m :: Mesh, Γ :: FaceField{Float64}, ∇φf :: FaceField{Point2D})
    ifcs = ifaces(m, true)
    jfcs = jfaces(m, true)
    for f in ifcs
        assignFaceDiffTerm!(A, b, m, Γ.ivals, ∇φf.ivals, f)
    end
    for f in jfcs
        assignFaceDiffTerm!(A, b, m, Γ.jvals, ∇φf.jvals, f)
    end
end

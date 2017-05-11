__precompile__()
module Fields

using JFVMM.Geometry
using JFVMM.MeshStructure
using JFVMM.BoundaryConditions

export CellField, createCellField
export linearMean!, cellGradient!
export linearMean

# %% Cell field
"""
Cell Field value

Cell field stores the values at cell centroids (= average over cell)
and values on the centroids of boundary faces. The latter are usually supplied via
boundary conditions. If BC are not supplied, the default assumption is
von Neumann BC (natural BC) with zero flux, in which case corresponding cell centroid
values are projected to boundary face centroids
"""
type CellField{V}
    domain :: Mesh
    cellvalues :: Array{V,2}
    iboundarymin :: Vector{V}
    iboundarymax :: Vector{V}
    jboundarymin :: Vector{V}
    jboundarymax :: Vector{V}
end

"""
```julia
cf[i,j]
```
applied to a cell field, gets the value of the cell field on the cell `(i,j)`
"""
Base.getindex{V}(cv :: CellField{V}, i :: Int64, j :: Int64) = cv.cellvalues[i,j]

Base.getindex{V}(cv :: CellField{V}, I :: CartesianIndex{2}) = cv.cellvalues[I]

Base.setindex!{V}(cv :: CellField{V}, x :: V, i :: Int64, j :: Int64) = cv.cellvalues[i,j] = x

Base.setindex!{V}(cv :: CellField{V}, x :: V, I :: CartesianIndex{2}) = cv.cellvalues[I] = x

"""
**INTERNAL!!!**

Creates arrays for boundary values of cell field. Is used internally to create cell field.
"""
function boundaryFaceValuesOfCellField(m :: Mesh, V :: DataType)
    Ni,Nj = dimensions(m)
    ibmin = Array{V}(Nj)
    ibmax = Array{V}(Nj)
    jbmin = Array{V}(Ni)
    jbmax = Array{V}(Ni)
    (ibmin, ibmax, jbmin, jbmax)
end

"""
**INTERNAL!!!**

Fills boundary values of the cell field based on a kind of value provided.
"""
function fillBoundaryValues!{V}(cf :: CellField{V}, x :: V)
    cf.iboundarymin[:] = x
    cf.iboundarymax[:] = x
    cf.jboundarymin[:] = x
    cf.jboundarymax[:] = x
end

function fillBoundaryValues!{V}(cf :: CellField{V}, T :: DataType)
    # no need to do anything: no value is provided
end

function fillBoundaryValues!{V}(cf :: CellField{V}, val :: Array{V,2})
    cf.iboundarymin[:] .= @view val[1,:]
    cf.iboundarymax[:] .= @view val[end,:]
    cf.jboudarymin[:] .= @view val[:,1]
    cf.jboundarymax[:] .= @view val[:,end]
end

function fillBoundaryValues!{V}(cf :: CellField{V}, f :: Function)
    m = cf.domain
    cf.iboundarymin[:] .= f.(centroid.(@view m.ifaces[1,:]))
    cf.iboundarymax[:] .= f.(centroid.(@view m.ifaces[end,:]))
    cf.jboundarybmin[:] .= f.(centroid.(@view m.jfaces[:,1]))
    cf.jboundarymax[:] .= f.(centroid.(@view m.jfaces[:,end]))
end

"""
**INTERNAL!!!**

Creates CellField with only cell values initialized: BC values are allocated but not filled.
"""
function createCellFieldInternal{V}(mesh :: Mesh, x :: V)
    val = fill(x, size(mesh))
    CellField(mesh, val, boundaryFaceValuesOfCellField(mesh,V)...)
end

function createCellFieldInternal(mesh :: Mesh, V :: DataType)
    val = Array{V}(size(mesh))
    CellField(mesh, val, boundaryFaceValuesOfCellField(mesh,V)...)
end


function createCellFieldInternal{V}(mesh :: Mesh, x :: Array{V,2})
    if size(mesh) != size(x)
        error("Incompatible size: the arrays size ($(size(x))) was expected to be equal to mesh size ($(size(mesh))")
    end
    val = copy(x)
    CellField(mesh, val, boundaryFaceValuesOfCellValue(m,V)...)
end

function createCellFieldInternal(mesh :: Mesh, f :: Function)
    val = [f(centroid(x)) for x in cells(mesh)]
    b = boundaryFaceValuesOfCellField(m, eltype(vals))
    CellField(mesh, val, b...)
end

"""
```julia
cv = createCellField(mesh, x)
```
Creates cell field for `mesh` and `x`, where `x` is
- a value, in which case, the uniform cell field is created
- a 2D array with cell values to be assigned
- a function of `Point2D` to calculate the value
"""
function createCellField(m :: Mesh, x)
    cf = createCellFieldInternal(m, x)
    fillBoundaryValues!(cf, x)
    cf
end

# %% Face Variable
"""
Face field

Represents the values defined on faces. Due to tetragon shape of the element, two disctinct
faces are present. Field `ivals` stores the values on the faces with normal in `i` direction.
Field `jvals` stores the values on the faces with normal in `j` direction.
"""
type FaceField{V}
    domain :: Mesh
    ivals :: Array{V,2}
    jvals :: Array{V,2}
end

"""
Produces in place linear mean of cell field for internal faces uncorrected for skewness
"""
function internalLinearMean!{U,V}(fv :: FaceField{U}, cv :: CellField{V})
    # loop over internal faces in i-direction
    for f in ifaces(cv.domain, true)
        phiC = cv[faceowneri(f)]
        phiF = cv[faceneighbouri(f)]
        gf = f.gf
        fv.ivals[meshindex(f)] = gf*phiF + (1-gf)*phiC
    end
    # loop over internal faces in j-direction
    for f in jfaces(cv.domain, true)
        phiC = cv[faceowneri(f)]
        phiF = cv[faceneighbouri(f)]
        gf = f.gf
        fv.jvals[meshindex(f)] = gf*phiF + (1-gf)*phiC
    end
    # use it as a rule: modification does not produce the value
    return nothing
end

"""
In-place linear mean projection of cell field values onto faces.
Boundary values of the cell field are used for boundary faces.
"""
function linearMean!{U,V}(fv :: FaceField{U}, cv :: CellField{V})
    internalLinearMean!(fv, cv)
    fv.ivals[1,:] .= cv.iboundarymin
    fv.ivals[end,:] .= cv.iboundarymax
    fv.jvals[:,1] .= cv.jboundarymin
    fv.jvals[:,end] .= cv.jboundarymax
    return nothing
end

"""
Creates a new face field as a result of the linear mean projection of the cell field
"""
function linearMean{V}(cv :: CellField{V})
    fv = FaceField(cv.domain,
        Array{V}(size(ifaces(cv.domain))),
        Array{V}(size(jfaces(cv.domain))))
    linearMean!(fv, cv)
    fv
end

function skewnessCorrection!(φf :: FaceField{Float64}, ∇φf :: FaceField{Point2D})
    m = φf.domain
    Ni,Nj = size(m)
    for i in 1:Ni
        φf.ivals[i,:] .= (@view φf.ivals[i,:]) .+ dot.((@view ∇φf.ivals[i,:]), centroid_correction.(@view m.ifaces[i,:]))
    end
    for j in 1:Nj
        φf.jvals[:,j] .= (@view φf.jvals[:,j]) .+ dot.((@view ∇φf.jvals[:,j]), centroid_correction.(@view m.jfaces[:,j]))
    end
    return nothing
end

"""
Correction to the gradient on the faces to base the gradient mostly on the cells straddling the face
"""
function gradCorrection(∇φf :: FaceFiled{Point2D}, φ :: CellField{Float64})
    m = φ.domain
    Ni,Nj = size(m)
    ifcs = ifaces(m, true)
    jfcs = jfaces(m, true)
    cls = cells(m)
    for face in ifcs
        I = meshindex(face)
        CF = face.CF
        φC = φ[faceowneri(face)]
        φF = φ[faceneighbouri(face)]
        ∇φf.ivals[I] = ∇φf.ivals[I] + (φF - φC - ∇φf.ivals[I] ⋅ CF)/(CF ⋅ CF) * CF
    end
    for face in jfcs
        I = meshindex(face)
        CF = face.CF
        φC = φ[faceowneri(face)]
        φF = φ[faceneighbouri(face)]
        ∇φf.jvals[I] = ∇φf.jvals[I] + (φF - φC - ∇φf.jvals[I] ⋅ CF)/(CF ⋅ CF) * CF
    end
end

# %% Cell gradient

function cellGradient!(gradphi :: CellField{Point2D}, phif :: FaceField{Float64})
    mesh = phif.domain
    ifcs = ifaces(mesh)
    jfcs = jfaces(mesh)
    for c in cells(mesh)
        VC = c.volume
        I1, I2 = cellifacesi(c)
        J1, J2 = celljfacesi(c)
        println("cellGradient!: cell $(c.meshindex)")
        # Gauss-Green theorem:
        gradphi[meshindex(c)] = ((c.ifaceSign[I1]*phif.ivals[I1]) * facevector(ifcs[I1]) + (c.ifaceSign[I2]*phif.ivals[I2]) * facevector(ifcs[I2]) + (c.jfaceSign[J1]*phif.jvals[J1]) * facevector(jfcs[J1]) + (c.jfaceSign[J2]*phif.jvals[J2]) * facevector(jfcs[J2])) ./ VC
    end
    return nothing
end

function cellGradient!(gradphi :: CellField{Point2D}, tmp :: FaceField{Float64}, phi :: CellField{Float64})
    linearMean!(tmp, phi)
    cellGradient!(gradphi, tmp)
end


function updateCellFieldBoundaryValues!{U}(cf :: CellField{Float64}, bc :: BoundaryCondition, gamma :: FaceField{U})
    mesh = cf.domain
    #imin
    cf.iboundarymin[:] .= boundaryvalue.(bc.imin, (@view mesh.ifaces[1,:]), (@view cf[1,:]), (@view gamma.ivals[1,:]))
    # imax
    cf.iboundarymax[:] .= boundaryvalue.(bc.imax, (@view mesh.ifaces[end,:]), (@view cf[end,:]), (@view gamma.ivals[end,:]))
    # jmin
    cf.jboundarymin[:] .= boundaryvalue.(bc.jmin, (@view mesh.jfaces[:,1]), (@view cf[:,1]), (@view gamma.jvals[:,1]))
    # jmax
    cf.jboundarymax[:] .= boundaryvalue.(bc.jmax, (@view mesh.jfaces[:,end]), (@view cf[:,end]), (@view gamma.jvals[:,end]))
end

function updateCellFieldBoundaryValues!{U}(cf :: CellField{Float64}, bc :: BoundaryCondition, gamma :: FaceField{U}, gradphi ::FaceField{Point2D})
    mesh = cf.domain
    # go over boundaries imin
    cf.iboundarymin[:] .= boundaryvalue.(bc.imin, (@view mesh.ifaces[1,:]), (@view cf[1,:]), (@view gamma.ivals[1,:]), (@view gradphi.ivals[1,:]))
    # imax
    cf.iboundarymax[:] .= boundaryvalue.(bc.imax, (@view mesh.ifaces[end,:]), (@view cf[end,:]), (@view gamma.ivals[end,:]), (@view gradphi.ivals[end,:]))
    # jmin
    cf.jboundarymin[:] .= boundaryvalue.(bc.jmin, (@view mesh.jfaces[:,1]), (@view cf[:,1]), (@view gamma.jvals[:,1]), (@view gradphi.jvals[:,1]))
    # jmax
    cf.jboundarymax[:] .= boundaryvalue.(bc.jmax, (@view mesh.jfaces[:,end]), (@view cf[:,end]), (@view gamma.jvals[:,end]), (@view gradphi.jvals[:,end]))
end

function createCellField{U,V}(x :: V, bc :: BoundaryCondition, gamma :: FaceField{U})
    mesh = bc.domain
    cf = createCellFieldInternal(mesh, x)
    updateCellFieldBoundaryValues!(cf, bc, gamma)
    cf
end

# Almost a repeat, but for now choose to express it explicitly
function createCellField{U}(x :: Float64, bc :: BoundaryCondition, gamma :: FaceField{U}, gradphi ::FaceField{Point2D})
    mesh = bc.domain
    cf = createCellFieldInternal(mesh, x)
    updateCellFieldBoundaryValues!(cf, bc, gamma, gradphi)
    cf
end


"""
```julia
updateCellFieldBoundaryValues!(φ, ∇tmp, tmpf, ∇tmpf, BC, Γ)
```

Updates cell field `φ` boundary values in accordance with provided BC and `Γ` (diffusivity) face
field. This function takes into account non-orthogonality of the mesh.

Implements a simple predictor-corrector: the first round is calculated assuming orthogonal mesh,
the second round performs the correction. Should be accurate for not too skewed meshes.

Temporary variables `tmpf` and `∇tmpf` are used to store temporary face fields for `φ` and `∇φ`.
"""
function updateCellFieldBoundaryValues!{U}(phi :: CellField{Float64}, dtmp :: CellField{Point2D}, tmpf :: FaceField{Float64}, dtmpf :: FaceField{Point2D}, bc :: BoundaryCondition, gamma :: FaceField{U})
    # Project phi into faces
    linearMean!(tmpf, phi, bc, gamma)
    # calculate cell gradients using this projection
    cellGradient!(dtmp, tmpf)
    # porject cell gradient to faces
    linearMean!(dtmpf, dtmp)
    # Use this to correct boundary values tmpf
    updateCellFieldBoundaryValues!(phi, bc, gamma, dtmpf)
end

end

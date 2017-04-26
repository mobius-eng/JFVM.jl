# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# Last edited: 30 December, 2014
# ===============================

# =====================================================================
# 2014-12-30: - cell variables can be created including bounddary cells
# =====================================================================

# %% ===================== Cell Variables ==================================
"""
Creates a cell variable on the mesh and assigns value `phi0` to it:

```julia
var = createCellVariable(mesh, phi0)
```

Inputs:

    + mesh: is a mesh structure created by one of `create*Mesh` functions
    + phi0 (real): is an initial value of the variable on each cell

Output:

    + var: cell variable

"""
function createCellVariable{T<:Real}(m::MeshStructure{T}, phi0::T)
    CellValue(m, phi0*ones(T,dimensions(m)...))
end

function createCellVariable{T<:Real, T2}(m::MeshStructure{T}, phi0::T2)
    dims = dimensions(m)
    if length(dims) == 2
        CellValue(m, [copy(phi0) for i in 1:dims[1], j in 1:dims[2]])
    else
        CellValue(m, [copy(phi0) for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]])
    end
end

"""
Creates a cell variable on the mesh from the existing values in the
array that follows the mesh structure shape:

```julia
var = createCellVariable(mesh, array)
```

Inputs:
    + mesh: is a mesh structure created by one of `create*Mesh` functions
    + array (array of reals): array of initial values

Output:

    + var: cell variable
"""
function createCellVariable{T<:Real, T2}(m::MeshStructure{T}, phi0::Array{T2})
    if prod(m.dims) == length(phi0)
        CellValue(m, phi0)
    else
        error("JFVMM: Matrix must be the same size ($(size(phi0))) as the domain.($(m.dims))")
    end
end

# %% ============== Face Variavle ======================================
"""
Creates a face variable (vector) on the mesh and assigns value `phi0` to it:

```julia
var = createFaceVariable(mesh, phi0)
```

Inputs:

    + mesh: is a mesh structure created by one of `create*Mesh` functions
    + phi0 (vector): is an initial value of the variable on each face;
        the dimension of phi0 must be the same as the domain (1, 2 or 3)

Output:

    + var: face variable

"""
function createFaceVariable{T<:Real}(m::MeshStructure{T}, phi0::Array{T,1})
    d = length(m.dims)
    if d == 2
        FaceValue(m,
	      ones(T, m.dims[1]+1, m.dims[2])*phi0[1],
	      ones(T, m.dims[1], m.dims[2]+1)*phi0[2],
	      ones(T,1))
    else # d == 3
        FaceValue(m,
	      ones(T,m.dims[1]+1, m.dims[2], m.dims[3])*phi0[1],
	      ones(T,m.dims[1], m.dims[2]+1, m.dims[3])*phi0[2],
	      ones(T,m.dims[1], m.dims[2], m.dims[3]+1)*phi0[3])
    end
end

"""
Creates a face variable (vector) on the mesh and assigns value `phi0` to it:

```julia
var = createFaceVariable(mesh, phi0)
```

Inputs:

- `mesh`: is a mesh structure created by one of `createMesh*` functions
- `phi0` (real): is an initial value of the variable on each face;
    each active direction in the face value will get this value

Output:

- `var`: face variable

"""
function createFaceVariable{T<:Real}(m::MeshStructure{T}, phi0::T)
    d=length(m.dims)
    if d==2
        FaceValue(m,
            ones(T,m.dims[1]+1, m.dims[2])*phi0,
            ones(T,m.dims[1], m.dims[2]+1)*phi0,
            ones(T,1))
    else # d==3
    FaceValue(m,
        ones(T,m.dims[1]+1, m.dims[2], m.dims[3])*phi0,
        ones(T,m.dims[1], m.dims[2]+1, m.dims[3])*phi0,
        ones(T,m.dims[1], m.dims[2], m.dims[3]+1)*phi0)
    end
end

# %% ============== Cell Vector =====================================
# TODO: docstring for createCellVector
function createCellVector{T<:Real}(m::MeshStructure{T}, phi0::Array{T,1})
    d=length(m.dims)
    if d==1
        CellVector(m,
    	      ones(T,m.dims[1])*phi0[1],
    	      ones(T,1),
    	      ones(T,1))
    elseif d==2
        CellVector(m,
    	      ones(T,m.dims[1], m.dims[2])*phi0[1],
    	      ones(T,m.dims[1], m.dims[2])*phi0[2],
    	      ones(T,1))
    elseif d==3
        CellVector(m,
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[1],
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[2],
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[3])
    end
end

function createCellVector{T<:Real}(m::MeshStructure{T}, phi0::T)
    d=length(m.dims)
    if d==1
        CellVector(m,
    	      ones(T,m.dims[1])*phi0,
    	      ones(T,1),
    	      ones(T,1))
    elseif d==2
        CellVector(m,
    	      ones(T,m.dims[1], m.dims[2])*phi0,
    	      ones(T,m.dims[1], m.dims[2])*phi0,
    	      ones(T,1))
    elseif d==3
        CellVector(m,
    	      ones(T,m.dims[1], m.dims[2], m.dims[3])*phi0,
    	      ones(T,m.dims[1], m.dims[2], m.dims[3])*phi0,
    	      ones(T,m.dims[1], m.dims[2], m.dims[3])*phi0)
    end
end

# %% Copy cell
"""
Copy cell value
"""
copyCell(phi::CellValue) = CellValue(phi.domain, Base.copy(phi.value))

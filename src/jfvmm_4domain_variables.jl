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
function createCellVariable(m::MeshStructure, phi0::Real)
    CellValue(m, phi0*ones(tuple(m.dims.+2...)))
end

"""
Creates a cell variable on the mesh and assigns value `phi0` to it.
Cell variable is forced to satisfy provided boundary conditions

```julia
var = createCellVariable(mesh, phi0, bc)
```

Inputs:

    + mesh: is a mesh structure created by one of `create*Mesh` functions
    + phi0 (real): is an initial value of the variable on each cell
    + bc: boundary conditions created on the mesh

Output:

    + var: cell variable

"""
function createCellVariable(m::MeshStructure, phi0::Real, BC::BoundaryCondition)
    phi = CellValue(m, phi0*ones(tuple(m.dims.+2...)))
    cellBoundary!(phi, BC)
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
function createCellVariable{T<:Real}(m::MeshStructure, phi0::Array{T})
    if prod(m.dims+2)==length(phi0)
      CellValue(m, phi0)
    elseif prod(m.dims)==length(phi0)
      d=m.dimension
      phival = zeros(tuple(m.dims.+2...))
      BC = createBC(m) # Neumann boundaries
      if d==1 || d==1.5
        phival[2:end-1] = phi0
      elseif d==2 d==2 || d==2.5 || d==2.8
        phival[2:end-1, 2:end-1] = phi0
      elseif d==3 || d==3.2
        phival[2:end-1,2:end-1,2:end-1] = phi0
      end
      phi = CellValue(m, phival)
      cellBoundary!(phi, BC)
    else
      error("JFVMM: Matrix must be the same size as the domain.")
    end
end


"""
Creates a cell variable on the mesh from the existing values in the
array that follows the mesh structure shape. The cell variable is forced to
satisfy provided boundary conditions:

```julia
var = createCellVariable(mesh, array, bc)
```

Inputs:
    + mesh: is a mesh structure created by one of `create*Mesh` functions
    + array (array of reals): array of initial values
    + bc: boundary conditions on the mesh

Output:

    + var: cell variable
"""
function createCellVariable{T<:Real}(m::MeshStructure, phi0::Array{T}, BC::BoundaryCondition)
    if prod(m.dims+2)==length(phi0)
      error("jFVT: Matrix must be the same size as the domain.")
    elseif prod(m.dims)==length(phi0)
      d=m.dimension
      phival = zeros(tuple(m.dims.+2...))
      if d==1 || d==1.5
        phival[2:end-1] = phi0
      elseif d==2 || d==2.5 || d==2.8
        phival[2:end-1, 2:end-1] = phi0
      elseif d==3 || d==3.2
        phival[2:end-1,2:end-1,2:end-1] = phi0
      end
      phi = CellValue(m, phival)
      cellBoundary!(phi, BC)
    else
      error("JFVMM: Matrix must be the same size as the domain.")
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
        the dimension of phi0 must be the same as the mesh (1,2 or 3)

Output:

    + var: face variable

"""
function createFaceVariable{T<:Real}(m::MeshStructure, phi0::Array{T,1})
d=m.dimension
  if d==1 || d==1.5
    FaceValue(m,
	      ones(m.dims[1]+1)*phi0[1],
	      [1.0],
	      [1.0])
  elseif d==2 || d==2.5 || d==2.8
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2])*phi0[1],
	      ones(m.dims[1], m.dims[2]+1)*phi0[2],
	      [1.0])
  elseif d==3 || d==3.2
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2], m.dims[3])*phi0[1],
	      ones(m.dims[1], m.dims[2]+1, m.dims[3])*phi0[2],
	      ones(m.dims[1], m.dims[2], m.dims[3]+1)*phi0[3])
  end
end

"""
Creates a face variable (vector) on the mesh and assigns value `phi0` to it:

```julia
var = createFaceVariable(mesh, phi0)
```

Inputs:

    + mesh: is a mesh structure created by one of `create*Mesh` functions
    + phi0 (real): is an initial value of the variable on each face;
        each active direction in the face value will get this value

Output:

    + var: face variable

"""
function createFaceVariable(m::MeshStructure, phi0::Real)
# creates a face variable based on the mesh structure
d=m.dimension
  if d==1 || d==1.5
    FaceValue(m,
	      ones(m.dims[1]+1)*phi0,
	      [1.0],
	      [1.0])
  elseif d==2 || d==2.5 || d==2.8
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2])*phi0,
	      ones(m.dims[1], m.dims[2]+1)*phi0,
	      [1.0])
  elseif d==3 || d==3.2
    FaceValue(m,
	      ones(m.dims[1]+1, m.dims[2], m.dims[3])*phi0,
	      ones(m.dims[1], m.dims[2]+1, m.dims[3])*phi0,
	      ones(m.dims[1], m.dims[2], m.dims[3]+1)*phi0)
  end
end


# %% ============== Cell Vector =====================================
# TODO: docstring for createCellVector
function createCellVector{T<:Real}(m::MeshStructure, phi0::Array{T,1})
    d=m.dimension
    if d==1 || d==1.5
        CellVector(m,
    	      ones(m.dims[1])*phi0[1],
    	      [1.0],
    	      [1.0])
    elseif d==2 || d==2.5 || d==2.8
        CellVector(m,
    	      ones(m.dims[1], m.dims[2])*phi0[1],
    	      ones(m.dims[1], m.dims[2])*phi0[2],
    	      [1.0])
    elseif d==3 || d==3.2
        CellVector(m,
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[1],
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[2],
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0[3])
    end
end

function createCellVector(m::MeshStructure, phi0::Real)
    d=m.dimension
    if d==1 || d==1.5
        CellVector(m,
    	      ones(m.dims[1])*phi0,
    	      [1.0],
    	      [1.0])
    elseif d==2 || d==2.5 || d==2.8
        CellVector(m,
    	      ones(m.dims[1], m.dims[2])*phi0,
    	      ones(m.dims[1], m.dims[2])*phi0,
    	      [1.0])
    elseif d==3 || d==3.2
        CellVector(m,
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0,
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0,
    	      ones(m.dims[1], m.dims[2], m.dims[3])*phi0)
    end
end

# %% Copy cell
"""
Copy cell value
"""
copyCell(phi::CellValue) = CellValue(phi.domain, Base.copy(phi.value))

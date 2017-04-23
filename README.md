# JFVMM(odified)

## New README under construction

### Finite volume method

JFVMM implements finite volume method for structured meshes. While structured
meshes are more limiting, the implementation of the method is (a) simpler and
(b) more efficient. Available geometries:
- 1D linear
- 1D polar with angular symmetry.
- 2D planar on the rectangle
- 2D polar on the disk segment.
- 2D cylindrical: cylindrical domain with angular symmetry
- 3D on rectangular block
- 3D cylindrical on the segment of the cylinder

### What's modified?

Original JFVM, which is awesome otherwise, is essentially a port from MATLAB
code, which lacks Julia's character. JFVMM tries to re-build it in a more
Julia-like manner.

### What is done and what still needs to be done
- MOSTLY DONE: document exported functions
- MOSTLY DONE: make the code generic over real type. This is mostly done, except
for TVD
- Factor the code: there are large bits of code that are repeated and need to
be factored out.
- The code allocates lots of arrays as it computes terms and solves equations:
this is the main area for optimization

## Bits of old README

### A simple finite volume tool written in Julia
This code is a Matlabesque implementation of my Matlab finite volume tool. The code is not in its most beautiful form, but it works if you believe my words. Please remember that the code is written by a chemical/petroleum engineer. Petroleum engineers are known for being simple-minded folks and chemical engineers have only one rule: "any answer is better than no answer". You can expect to easily discretize a linear transient convection-diffusion PDE into the matrix of coefficients and RHS vectors. Domain shape is limited to rectangles, circles (or a section of a circle), cylinders, and soon spheres. The mesh can be uniform or nonuniform:
  - Cartesian (1D, 2D, 3D)
  - Cylindrical (1D, 2D, 3D)
  - Radial (2D r and \theta)

You can have the following boundary conditions or a combination of them on each boundary:
  - Dirichlet (constant value)
  - Neumann (constant flux)
  - Robin (a linear combination of the above)
  - Periodic (so funny when visualize)

It is relatively easy to use the code to solve a system of coupled linear PDE's and not too difficult to solve nonlinear PDE's.

### Tutorial
I have written a short [tutorial](http://nbviewer.ipython.org/github/simulkade/JFVM.jl/blob/master/examples/jfvm_tutorial.ipynb), which will be extended gradually.

### In action
Copy and paste the following code to solve a transient diffusion equation:
```julia
using JFVM
Nx = 10
Lx = 1.0
m = createMesh1D(Nx, Lx)
BC = createBC(m)
BC.left.a[:]=BC.right.a[:]=0.0
BC.left.b[:]=BC.right.b[:]=1.0
BC.left.c[:]=1.0
BC.right.c[:]=0.0
c_init = 0.0 # initial value of the variable
c_old = createCellVariable(m, 0.0, BC)
D_val = 1.0 # value of the diffusion coefficient
D_cell = createCellVariable(m, D_val) # assigned to cells
# Harmonic average
D_face = harmonicMean(D_cell)
N_steps = 20 # number of time steps
dt= sqrt(Lx^2/D_val)/N_steps # time step
M_diff = diffusionTerm(D_face) # matrix of coefficient for diffusion term
(M_bc, RHS_bc)=boundaryConditionTerm(BC) # matrix of coefficient and RHS for the BC
for i =1:5
    (M_t, RHS_t)=transientTerm(c_old, dt, 1.0)
    M=M_t-M_diff+M_bc # add all the [sparse] matrices of coefficient
    RHS=RHS_bc+RHS_t # add all the RHS's together
    c_old = solveLinearPDE(m, M, RHS) # solve the PDE
    visualizeCells(c_old)
end
```

Now change the 4th line to `m=createMesh2D(Nx,2*Nx, Lx,2*Lx)` and see what happens.

### IJulia notebooks
  - [Introduction](http://nbviewer.ipython.org/github/simulkade/JFVM.jl/blob/master/examples/jfvm-a-finite-volume-tool-for-julia.ipynb)
  - [tutorial](http://nbviewer.ipython.org/github/simulkade/JFVM.jl/blob/master/examples/jfvm_tutorial.ipynb)
  - [compare analytical solution of a diffusion equation with uniform and nonuniform grids](http://nbviewer.ipython.org/github/simulkade/JFVM.jl/blob/master/examples/jfvm_diffusion_analytics.ipynb)
  - [Nonlinear PDE](http://nbviewer.ipython.org/github/simulkade/JFVM.jl/blob/master/examples/solving-nonlinear-pdes-with-fvm.ipynb)
  - [Foam flow in porous media](https://github.com/simulkade/FoamFingering/blob/master/foam_flow_stars_numeric.ipynb)
  - New notebooks soon...

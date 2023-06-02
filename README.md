# Multigrid
A multigrid solver for the 1D heat equation. It will try to estimate the condition number of the given Matrix by using the small grids. Depended on this it will either use a simple jacobi or the more complicated mutigrid solver. 

When running the multigrid solver it tries to find an optimal overdamping/overrelaxation to speed up convergence.

Neither the condition estimate nor the overrelaxation are based on any theroy and were just developed by me playing around.

# Building and Running
Eigen is required and was located in the directory above the main branch. 

Build and run by typing make

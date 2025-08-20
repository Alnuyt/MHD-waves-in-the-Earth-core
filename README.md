# Simplified model of magnetohydrodynamic waves in the Earth's core
In this thesis, the linear equations of hydrodynamics and magnetohydrodynamics are developed for a uniformly rotating conductive fluid impregnated with a constant magnetic field. A toroidal-poloidal decomposition is applied to the solenoidal fields, allowing us to obtain a numerically solvable system. The numerical solution of these equations using an inertial wave solver, developed as part of this thesis, allows us to calculate all the oscillation modes of these waves for the hydrodynamic case. The basis for such a solver was constructed during this study, and improvements will lead to a true magneto-Coriolis wave calculator applicable to geophysics and astrophysics.
## Author
Alexandre Nuyt - Catholic University of Louvain, Faculty of Science, School of Physics
## Objectives
As this thesis is exploratory, three major objectives are targeted:
- To study the fundamentals of rotating fluid dynamics and magnetohydrodynamics
in a geophysical context.
- To develop the differential equations of inertial waves for hydrodynamics
and MHD in toroidal-poloidal decomposition.
- Build the basis for a numerical tool: Hydro + MHD solver for inertial waves,
for all m and k modes, in cylindrical geometry.
## Description of the codes
During the study, three main codes were written:
- [Analytical_Inertial_waves.py](Analytical_Inertial_waves.py): Numerically solves the transcendental equation (1.56) which gives the analytical solutions of inertial waves propagating in a rotating non-conductive fluid for a cylindrical geometry. 
- [Helmoltz_cylindrical.py](Helmoltz_cylindrical.py): Solves the cylindrical Helmholtz equation using the finite difference method to ensure the functionality of the numerical method used to create the Hydro and MHD solvers for inertial waves.
- [Toroidal_Poloidal_Inertial_waves.py](Toroidal_Poloidal_Inertial_waves.py): Forms the basis of a Hydro + MHD solver for inertial waves, for all modes m and k, in cylindrical geometry.
## Conclusion
This work explores the resolution of magneto-Coriolis waves in cylindrical geometry based on theoretical foundations in hydrodynamics and magnetohydrodynamics. A numerical algorithm has been developed to address this problem, although the current results still require improvement. By optimising the discretisation, this tool could become a real asset for the study of geophysical and astrophysical phenomena.

Although my master's thesis is complete, this **project continues** and progress can be found in the repository ‘**Magneto-Coriolis Waves Solver**’.

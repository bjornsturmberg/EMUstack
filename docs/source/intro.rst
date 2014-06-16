Introduction
================


EMUstack calculates the scattering matrices of a multi-layered structure, where each layer can be homogeneous or structured (down to sub-wavelength dimensions) and the materials may have complex, dispersive refractive indices. The scattering matrices are powerful tools from which many physical quantities, such as the total transmission, absorption in each layer, and the resonances of the structure, can be derived.

-----
An advantage of EMUstack over other scattering matrix programs (for example `CAMFR <http://docutils.sf.net/rst.html>`_) is that the fields in each layer are considered in their natural basis with transmission scattering matrices converting fields between them. The fields in homogeneous layers are therefore expressed in terms of plane waves, while the natural basis in the periodically structured layers are Bloch modes. Expressing fields in their natural basis gives the terms of the scattering matrices intuitive meaning, providing access to greater physical insights. It is also advantages for the speed and accuracy of the numerical method.

Inherent to the scattering matrix approach is the requirement that the interfaces between layers be planar, *ie.* that each layer is uniform in one direction (here labelled {\it z}). In this nomenclature the incident field must have $k_z = k_{\perp} \ne 0$, but is unconstrained in $k_{\parallel} = k_{x,y}$.
In our implementation, the only constraint placed on each layer in the x-y plane, is that it must be periodic, at least at the supercell level. This is because the modes of structured media are calculated using a vectorial Finite Element Method (FEM) routine with periodic boundary conditions. 
The scattering matrices of homogeneous media are calculated analytically resulting in excellent accuracy and speed.

EMUstack has been designed to handle lossy media with dispersive refractive indices, with the complex refractive index at each frequency being taken directly from tabulated results of experimental measurements. This is an advantage of frequency domain methods over time domain methods such as the Finite Difference Time Domain (FDTD) where refractive indices are included by analytic approximations such as the Drude model (for example `MEEP <http://ab-initio.mit.edu/wiki/index.php/Meep>`_). It is also possible to include media with lossless and/or non-dispersive refractive indices in EMUstack.

Taking full advantage of the boundary-element nature of the scattering matrix method it is possible to vary the thickness of a layer by a single, numerically inexpensive, matrix multiplication. Furthermore, EMUstack recognises when interfaces are repeated so that their scattering matrices need not be recalculated but rather just retrieved from memory, which takes practically no computation time.

By integrating a 2D finite element method calculation into the scattering matrix method EMUstack provides a powerful, versatile tool for nanophotonic simulations that are both computationally efficient and physically insightful. 

EMUstack is a completely open source package, utilising free, open source compilers, meshing programs and libraries! The low-level numerical routines are written in Fortran for optimal performance, while higher-level processing is done in python. EMUstack currently comes with template FEM mesh for 1D and 2D gratings, such as lamellar gratings and Nanowire/Nanohole arrays. For these structures the EMUstack will automatically create FEM mesh with the specified parameters. For other structures, the open source program `gmsh <http://geuz.org/gmsh/>`_ may be used to create the FEM mesh. 

In summary, the advantages of EMUstack are;

* Calculates the scattering matrices between layers in their natural bases, for maximum physical insights.
* Designed to include lossy, dispersive materials, with frequency resolved (experimentally measured) refractive indices.
* FEM allows for arbitrary in-plane geometries, down to the periodicity of the supercell.
* Homogeneous layers are calculated analytically for optimal accuracy and speed.
* The scattering matrix method efficiently combines arbitrary number of layers into a stack.
* Synthesis of efficient Fortran routine with dynamic, high-level programming in Python.
* Completely open source package! Including FEM meshing program, Fortran FEM routine, Python multi-layered scattering matrix implementation. 
* Integrated with highly optimised libraries (but also functions without these at slower speeds), including; BLAS, LAPACK, ARPACK, UMFPACK
* Vectorial FEM advantages?
* Get band structure at same time as t,r,a.
* Both/all polarasations at once



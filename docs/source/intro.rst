.. role:: raw-math(raw)
    :format: latex html

Introduction
================

EMUstack is an open-source simulation package for calculating light propagation through multi-layered stacks of dispersive, lossy, nanostructured, optical media. It implements a generalised scattering matrix method, which extends the physical intuition of thin film optics to complex structures.

At the heart of the scattering matrix approach is the requirement that each layer is uniform in one direction, here labelled *z*. In this nomenclature the incident field is unconstrained in :raw-math:`$ k_{\parallel} = k_{x,y} $` but must have :raw-math:`$ k_{\perp} = k_z \ne 0 $`.

In-plane each layer can be homogeneous, periodic in x or y, or double periodic (periodic in x and y). The modes of periodic (structured layers) are calculated using the Finite Element Method in respectively 1 or 2 dimensions, while the modes of homogeneous media are calculated analytically. This approach maximises the speed and accuracy of the calculations.
These layers can be stacked in arbitrary order.

An advantage of EMUstack over other scattering matrix methods (for example `CAMFR <http://docutils.sf.net/rst.html>`_) is that the fields in each layer are considered in their natural basis with transmission scattering matrices converting fields between them. The fields in homogeneous layers are expressed in terms of plane waves, while the natural basis in the periodically structured layers are Bloch modes. Expressing fields in their natural basis gives the terms of the scattering matrices intuitive meaning, providing access to greater physical insights. It is also advantages for the speed and accuracy of the numerical method.

EMUstack has been designed to handle lossy media with dispersive refractive indices, with the complex refractive index at each frequency being taken directly from tabulated results of experimental measurements. This is an advantage of frequency domain methods over time domain methods such as the Finite Difference Time Domain (FDTD) where refractive indices are included by analytic approximations such as the Drude model. It is also possible to include media with lossless and/or non-dispersive refractive indices and EMUstack comes with a built in Drude model.

Taking full advantage of the boundary-element nature of the scattering matrix method it is possible to vary the thickness of a layer by a single, numerically inexpensive, matrix multiplication. Furthermore, EMUstack recognises when interfaces are repeated so that their scattering matrices need not be recalculated but rather just retrieved from memory, which takes practically no computation time.

EMUstack is a completely open source package, utilising free, open source compilers, meshing programs and libraries. 
All user interaction with EMUstack is done using the dynamic and easy to script language of python.
The low-level numerical routines are written in Fortran for optimal performance making use of the LAPACK, ARPACK, and UMFPACK libraries. The Fortran routines are compiled as python subroutines using f2py.
EMUstack currently comes with template FEM mesh for 1D and 2D gratings, Nanowire/Nanohole arrays, elliptical inclusions and split ring resonators. The mesh of other structures may be easily created using the open source program `gmsh <http://geuz.org/gmsh/>`_. 

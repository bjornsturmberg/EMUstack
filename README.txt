



Advantages;
- Fequency domain so can include real disperisive refractive index data.
- Designed to include lossy materials.
- FEM allows for arbitrary geometries in x-y plane.
- Vectorial FEM ... advantages.
- Analytic treatment of honogeneous films, including dispersive and lossy layers.
- Only need to calculate scattering matrices of each unique layer once per wavelength.
- Can then combine layers in arbitrary stack(s), including nanostructured layers directly onto one another.
- Repeated interfaces are recognised and not re-calcuated.
- Completely open source and compatibale with freely available compilers for fortran and python. 
- Integrated with highly optimised libries (but also functions without these at slower speeds), including;
    - BLAS
    - LAPACK
    - UMFPACK
    - ARPACK





Download PCPV
-------------------------------------------------
* github checkout ...




Hi!

BlochCode can calculate Bloch modes, Bloch factors and impedances for 2D
Photonic Crystals with up-down symmetry. To do this, it needs to be given
electric and magnetic field data sampled throughout the PC at the frequency,
incident angle, and polarization of interest.


The extensive human-readable documentation for PCPV lives in PCPV/Docs/

Example simulations are located in 000-simmo_templates/

License GNU GENERAL PUBLIC LICENSE found in LICENSE-GNU_GPL.txt


OFFICIALLY SUPPORTED PLATFORMS
------------------------------

BlochCode is most thoroughly tested on Mac OS X, running inside Sage <http://sagemath.org>.  On the Mac, it is also known to work with Enthought Python Distribution 7.1.

BlochCode has also been successfully run inside Sage on Linux, and should also work with Enthought Python Distribution (or other recent distributions of Numpy and Scipy).

It is currently unknown whether BlochCode works with Windows.  It can probably be made to work with Enthought on Windows with little effort -- but this is still uncharted territory.


TESTING
-------

BlochCode comes with a range of unit tests to ensure it is running correctly.  These are found in the tests directory.


HELP!
-----

Limited community support is available on the Launchpad site: <https://launchpad.net/blochcode>.  Feel free to ask questions.


CONTRIBUTING
------------

If you make an improvement to BlochCode, please share it with others by contributing back to the project, via the Launchpad site.

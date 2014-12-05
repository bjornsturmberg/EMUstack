#!/bin/bash

function makeEMUs() {
## compile Fortran routines
  cd backend/fortran/
  make
## run tests
  cd ../../tests/
  nosetests

##
  echo ''
  echo '#####################################################################'
  echo '        EMUstack and its dependencies have been installed        '
  echo ' Last tests may give errors if not using Ubuntu 12.04, gmsh 2.5  '
  echo '    ... as long as the tests run, things are okay.  '
  echo ''
  echo 'EMUstack is brought to you by Bjorn Sturmberg, Kokou Dossou,'
  echo 'Felix Lawrence and Lindsay Botton, with support from CUDOS and ARENA'
  echo '#####################################################################'
  echo ''

}

makeEMUs


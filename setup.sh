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
  echo ''
  echo 'EMUstack is brought to you by Bjorn Sturmberg, Kokou Dossou,'
  echo 'Felix Lawrence and Lindsay Botton, with support from CUDOS and ARENA'
  echo '#####################################################################'
  echo ''

}

makeEMUs

